# Forward Model

The forward model maps stellar parameters to an observer-frame SED and synthetic photometry.

## Pipeline

```text
(Teff, logg, [M/H])  →  SED interpolation   (Fortran, Hermite or linear)
                     →  distance dilution    (Fortran, (R/d)²)
                     →  extinction           (Python, optional)
                     →  bolometric           (Fortran)
                     →  filter convolution   (Fortran, per filter)
                     →  ForwardResult
```

Extinction is applied **after** dilution and **before** filter convolution — physically: star → distance → dust column → telescope → filter.

## Two calling conventions

### Classic

```python
from sed_model import run_forward

result = run_forward(
    teff=5778.0,
    logg=4.44,
    meta=0.0,
    R=6.957e10,        # cm
    d=3.086e19,        # cm
    grid=grid,
    filters=filters,
    mag_system="Vega",          # "AB" (default), "Vega", or "ST"
    interp_method="hermite",    # "hermite" (default) or "linear"
)
```

### FitParams (used internally by the inverse model)

```python
result = run_forward(
    fit_params=params,   # a FitParams
    theta=theta_vec,     # only the FREE parameters, canonical order
    R=R_sun,
    grid=grid,
    filters=filters,
)
```

`theta` contains only the free parameters in the canonical order `teff, logg, meta, a_v, d` (skipping fixed ones); `fit_params.unpack(theta)` fills in the fixed values. If `Av` or distance are free, they are taken from `theta`, not from keyword arguments. This is the exact code path the MCMC likelihood uses, so anything you verify in the forward direction holds inside the sampler too. See [Parameter Modes](parameters.md).

## Interpolation

Two schemes operate on the 3D `(Teff, logg, [M/H])` flux cube, vectorised over wavelength:

- **`hermite`** — cubic Hermite tensor interpolation with finite-difference tangents (central differences in the interior, one-sided at grid edges).
- **`linear`** — trilinear interpolation.

If the query point falls outside the grid, both schemes fall back to the **nearest grid node** and the result is flagged (`clamped=True`). Interpolated fluxes are floored at a tiny positive value, never negative.

The `interp_radius` diagnostic (distance in normalised parameter space from the nearest grid node) is computed for every call so you can gauge how far from tabulated models the interpolation reached.

## `ForwardResult`

| Field | Description |
|---|---|
| `wavelengths` | wavelength grid (Å) |
| `surface_flux` | interpolated surface flux (erg/s/cm²/Å) |
| `observed_flux` | after `(R/d)²` dilution and any extinction |
| `magnitudes` | `{filter_name: magnitude}` in the requested system |
| `band_fluxes` | `{filter_name: in-band flux}` (erg/s/cm²/Å) |
| `bol_flux`, `bol_mag` | bolometric flux and magnitude (`Mag_bol = -2.5 log10 F_bol`, MESA colors convention — no additional zero-point) |
| `interp_radius`, `clamped` | interpolation diagnostics |
| `teff`, `logg`, `meta`, `R`, `d`, `a_v` | the parameters that produced this result (self-describing) |
| `mag_system`, `extinction_applied` | bookkeeping |

## Synthetic photometry

For each filter, the transmission is linearly interpolated onto the SED wavelength grid, then the photon-counting in-band flux is computed:

\[
F_{\rm band} = \frac{\int F(\lambda)\, T(\lambda)\, \lambda\, d\lambda}{\int T(\lambda)\, \lambda\, d\lambda},
\qquad
m = -2.5 \log_{10}\!\left( F_{\rm band} / F_{zp} \right)
\]

with the zero-point \(F_{zp}\) precomputed per filter and system at load time ([Grids and Filters](grids-and-filters.md#zero-points)). Failed integrations or non-positive fluxes produce a sentinel magnitude of `-99.9` with a non-zero error flag rather than crashing.

## Batch evaluation

```python
from sed_model import run_forward_batch
import numpy as np

params = np.array([
    [5778.0, 4.44, 0.0],
    [6200.0, 4.30, -0.2],
])
results = run_forward_batch(params, R=R, d=d, grid=grid, filters=filters)
```

`run_forward_batch` accepts an `(N, 3)` array of `(teff, logg, meta)` rows (classic convention only) and returns a list of `ForwardResult`.
