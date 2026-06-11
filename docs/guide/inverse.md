# Inverse Model

The inverse model recovers a posterior over stellar parameters from observed broadband photometry, using [`emcee`](https://emcee.readthedocs.io/) ensemble MCMC.

## Statistical model

- **Prior** — flat within the bounds of each free parameter (from [`FitParams`](parameters.md)); \(-\infty\) outside.
- **Likelihood** — Gaussian per filter, summed:

\[
\ln \mathcal{L} = -\frac{1}{2} \sum_i \left[ \frac{(m_i^{\rm obs} - m_i^{\rm model})^2}{\sigma_i^2} + \ln\left(2\pi \sigma_i^2\right) \right]
\]

Model magnitudes come from `run_forward` in FitParams mode — the **same code path as the forward model**, with no separate likelihood implementation to drift out of sync. Any exception inside the forward evaluation, or a non-finite model magnitude, returns \(-\infty\) so the walker simply moves away.

## Basic usage

```python
from sed_model import run_inverse

posterior = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10,           # cm — always required, not sampled
    d=3.086e19,           # cm — convenience: builds default FitParams with d fixed
    grid=grid,
    filters=filters,
    mag_system="AB",
    n_walkers=32, n_steps=2000, n_burn=500, n_thin=1,
    seed=42,
    progress=True,
)
```

With no `fit_params` supplied, a default is built from the grid: `Teff`, `logg`, `[M/H]` free across the full grid, `Av` fixed at 0, distance fixed at `d` (default 1 pc). If `fit_params` **is** supplied, the `d` keyword is ignored — set distance through `fit_params`.

## Choosing what to sample

```python
from sed_model import fit_params_from_grid
from sed_model.params import PC_TO_CM

# Av free with flat prior over [0, 3] mag
params = fit_params_from_grid(grid, a_v=(0.0, 3.0), d_cm=500 * PC_TO_CM)

# Av and distance both free — expect a correlated Av–d posterior
params = fit_params_from_grid(
    grid,
    a_v=(0.0, 2.0),
    d_cm=(100 * PC_TO_CM, 2000 * PC_TO_CM),
)

posterior = run_inverse(..., fit_params=params)
```

When `Av` is active (free, or fixed > 0) and no `extinction` model is supplied, a default Fitzpatrick (1999) model is created automatically; pass `extinction_law=...` or a full `ExtinctionModel` to control the prescription. See [Extinction](extinction.md).

## Walker initialisation

Walkers start in a Gaussian ball of width `p0_scatter` (a fraction of each free parameter's range, default 0.02), centred on the midpoint of each parameter's bounds unless you supply a starting point:

```python
posterior = run_inverse(
    ...,
    p0_centre={"teff": 5800, "logg": 4.4, "meta": 0.0},
    # or the convenience aliases (p0_centre wins on conflict):
    p0_teff=5800, p0_logg=4.4, p0_meta=0.0,
    p0_scatter=0.02,
)
```

Initial positions are clipped to the prior bounds.

## Input validation

`run_inverse` raises `ValueError` before any sampling when:

- `obs_magnitudes` and `obs_uncertainties` differ in length, or `filter_names` does not match them;
- any uncertainty is ≤ 0 (the likelihood would be undefined);
- a filter name is not among the supplied filters;
- `n_walkers < 6`, `n_walkers < 2 × n_free`, or `n_burn >= n_steps`;
- every parameter is fixed (nothing to sample).

## Diagnostics

After sampling, a `UserWarning` is emitted if any walker's acceptance fraction falls below 0.1 or above 0.9. The integrated autocorrelation time is estimated with `sampler.get_autocorr_time(quiet=True)` and stored as `autocorr_time` (or `None` if the chain was too short for a reliable estimate).

## The result

`run_inverse` returns an [`InverseResult`](results-io.md#inverseresult) holding the flattened post-burn-in chain (`samples`, shape `(n_samples, n_free)` with columns in `param_names` order), per-sample `log_prob`, the observations, sampler settings, acceptance fractions, and the fixed-parameter values used. See [Results and I/O](results-io.md) for summaries, MAP estimates, and `.npz` persistence.
