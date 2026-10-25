# Forward Model

The forward model maps stellar parameters to an observer-frame SED and synthetic photometry.

## Pipeline

```text
(Teff, logg, [M/H])  →  SED interpolation   (Fortran, bounded Hermite, Hermite or linear)
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
    interp_method="hermite_bounded",  # "hermite_bounded" (default), "hermite" or "linear"
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

Three schemes operate on the 3D `(Teff, logg, [M/H])` flux cube, vectorised over wavelength. They follow the current MESA colors module:

- **`hermite_bounded`** (default) — Hermite interpolation with a safeguard. If the Hermite SED is finite and usable, negligible negative undershoot (below 0.1 % of the positive SED, by both integral and peak) is clipped to zero. Otherwise the call falls back to trilinear interpolation. If neither gives a valid SED, `run_forward` raises `ValueError` (the inverse model treats that as `-inf` log-likelihood).
- **`hermite`** — raw cubic Hermite tensor interpolation. Tangents on non-uniform grids are the spacing-weighted average of the one-sided slopes in the interior and one-sided at grid edges. Axes with a single grid point are treated as degenerate dimensions. No positivity floor is applied, so the result can undershoot below zero.
- **`linear`** — trilinear interpolation. A non-finite or non-positive result falls back to the nearest grid node; values are floored at `1e-300`.

If the query point falls outside the grid, it is clamped to the grid boundary and the result is flagged (`clamped=True`). `clamped` is also set when the bolometric flux is not positive, or when plain `hermite` reports an interpolation error. It is not set by `linear` filling a few wavelength bins from the nearest node, and `hermite_bounded` raises instead of flagging when it has no usable SED.

### Empty grid cells

A grid directory holds a full rectangular `(Teff, logg, [M/H])` cube, but not every node has a model. For `Kurucz2003all` about 40 % of the cube is empty, mostly the hot, low-gravity corners. Empty nodes are stored as all-zero spectra, and interpolating through them gives unreliable results.

`run_forward` emits a `UserWarning` ("lies in a grid cell with empty nodes") when:

- any corner of the cell containing the query is empty, for every interpolation method. A query exactly on a populated node never warns, because the empty neighbours carry no weight. The result can be wildly wrong here: with an empty corner the flux underflows and magnitudes come out near +770.
- for `hermite` and `hermite_bounded` only, any node the derivative estimates read is empty, even though all eight corners are populated. These are the axis-aligned neighbours of the corners: one extra node on each side of the cell along each axis. Hermite next to such a gap can be off by several percent in bolometric flux (up to about 0.08 mag in photometry in the comparisons described under [Testing](../testing.md)), which is MESA's behaviour and is reproduced exactly. `linear` does not use these nodes and is not affected.

The warning is advisory: the interpolated values are returned unchanged. `grid.cell_has_empty_nodes(teff, logg, meta, stencil=True)` runs the same check without a forward call. The inverse model avoids empty corners by treating them as zero posterior (see [Inverse Model](inverse.md)).

The `nearest_grid_distance` diagnostic (distance in normalised parameter space from the nearest tabulated model, MESA's `Interp_rad`) is computed for every call so you can gauge how far from tabulated models the interpolation reached.

## `ForwardResult`

| Field | Description |
|---|---|
| `wavelengths` | wavelength grid (Å) |
| `surface_flux` | interpolated surface flux (erg/s/cm²/Å) |
| `observed_flux` | after `(R/d)²` dilution and any extinction |
| `magnitudes` | `{filter_name: magnitude}` in the requested system |
| `band_fluxes` | `{filter_name: in-band flux}` (erg/s/cm²/Å) |
| `bol_flux`, `bol_mag` | bolometric flux and magnitude (`Mag_bol = -2.5 log10(F_bol / F_bol,0)`, where `F_bol,0` is the IAU 2015 B2 irradiance of an `M_bol = 0` source, so the Sun at 10 pc has `M_bol = 4.74`; samples with a non-finite or negative flux are zeroed before integrating) |
| `nearest_grid_distance`, `clamped` | interpolation diagnostics: distance to the nearest tabulated model, and whether the query was clamped or an interpolation failed |
| `teff`, `logg`, `meta`, `R`, `d`, `a_v` | the parameters that produced this result (self-describing) |
| `mag_system`, `extinction_applied` | bookkeeping |

## Synthetic photometry

For each filter, the transmission is linearly interpolated onto the SED wavelength grid, then the photon-counting in-band flux is computed:

\[
F_{\rm band} = \frac{\int F(\lambda)\, T(\lambda)\, \lambda\, d\lambda}{\int T(\lambda)\, \lambda\, d\lambda},
\qquad
m = -2.5 \log_{10}\!\left( F_{\rm band} / F_{zp} \right)
\]

with the zero-point \(F_{zp}\) precomputed per filter and system at load time ([Grids and Filters](grids-and-filters.md#zero-points)). Failed integrations or non-positive fluxes produce a non-zero error flag rather than crashing: band magnitudes become `NaN` in `ForwardResult.magnitudes`, and the bolometric magnitude is `99.0` (MESA's value) when the integrated flux is not positive.

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
