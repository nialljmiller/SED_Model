# Extinction

`sed_model.sed_extinction` is a self-contained, pure-NumPy implementation of interstellar dust extinction — no `dust_extinction`, `extinction`, or `astropy` dependency.

## Supported laws

| Key | Reference | Notes |
|---|---|---|
| `ccm89` | Cardelli, Clayton & Mathis (1989), ApJ 345, 245 | classic \(R_V\)-parameterised Milky Way law |
| `odonnell94` | O'Donnell (1994), ApJ 422, 158 | CCM89 with revised optical coefficients |
| `fitzpatrick99` | Fitzpatrick (1999), PASP 111, 63 | \(R_V\)-dependent spline; the package default |
| `fm07` | Fitzpatrick & Massa (2007), ApJ 663, 320 | \(R_V = 3.1\) fixed by construction |
| `calzetti00` | Calzetti et al. (2000), ApJ 533, 682 | starburst attenuation; default \(R_V = 4.05\), here \(A_V = R_V \, E(B\!-\!V)_s\) |
| `gordon23` | Gordon et al. (2023), ApJ 950, 86 | piecewise UV–optical–IR averages; `gordon23_env` selects `"mw"`, `"lmc"`, or `"smc"` |

The tuple `sed_model.AVAILABLE_LAWS` lists these keys programmatically.

## Conventions

- Wavelengths in **Å**; extinction returned as \(A(\lambda)\) in **magnitudes**.
- Flux attenuation: \(F_{\rm obs}(\lambda) = F(\lambda) \, 10^{-0.4 A(\lambda)}\).
- `a_v = 0` means no extinction regardless of the law.

## Usage

```python
from sed_model import make_extinction_model, ExtinctionModel

ext = make_extinction_model(
    enabled=True,
    law="gordon23",
    a_v=0.8,
    r_v=3.1,            # ignored by fm07 (fixed) and gordon23 (per-environment)
    gordon23_env="mw",  # "mw", "lmc", or "smc" — gordon23 only
)

a_lambda  = ext.extinction_mag(wavelengths)   # A(λ) in magnitudes
reddened  = ext.apply(wavelengths, flux)      # F · 10^(−0.4 A)
deredden  = ext.remove(wavelengths, flux)     # F · 10^(+0.4 A)
```

A disabled model (`enabled=False`, or `make_extinction_model()` with no arguments) passes flux through unchanged — so you can keep a single code path and toggle extinction by configuration.

In the forward model, pass the model via `run_forward(..., extinction=ext)`; it is applied after distance dilution and before filter convolution.

## Extinction in the inverse fit

When the [`FitParams`](parameters.md) makes `Av` **free**, the likelihood rebuilds the extinction model at each MCMC step with the proposed `Av` value (preserving your chosen `law`, `r_v`, and `gordon23_env`). The `a_v` stored on the `ExtinctionModel` you pass in is only a placeholder in that case.

If `Av` is active and you don't supply an extinction model, `run_inverse` creates a Fitzpatrick (1999), \(R_V = 3.1\) model automatically.

The [Extinction Comparison demo](../demos/extinction-comparison.md) fits the same synthetic observations once per prescription (including the three Gordon+2023 environments and non-standard \(R_V\) values) to show how the choice of law propagates into the recovered parameters.
