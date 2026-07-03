# Changelog

## [0.2.0] — 2026-07-03

### Breaking changes

- **`ForwardResult.interp_radius` renamed to `nearest_grid_distance`.**
  The old name implied a geometric radius; the new name describes what the
  value actually is: Euclidean distance in normalised parameter space from the
  nearest grid point.  Update any code that reads `result.interp_radius` or
  calls `grid.interp_radius(...)`.

- **`run_inverse` — `p0_teff`, `p0_logg`, `p0_meta` kwargs removed.**
  Use the unified `p0_centre` dict instead:
  ```python
  # before
  run_inverse(..., p0_teff=5800, p0_logg=4.4, p0_meta=0.0)
  # after
  run_inverse(..., p0_centre={'teff': 5800, 'logg': 4.4, 'meta': 0.0})
  ```

- **`fit_params_from_grid` — all overrides are now keyword-only.**
  `a_v` and `d_cm` were previously positional; they are now keyword-only
  alongside `teff`, `logg`, and `meta`:
  ```python
  # before (worked but fragile)
  fit_params_from_grid(grid, 0.3, 500*PC_TO_CM)
  # after (always required)
  fit_params_from_grid(grid, a_v=0.3, d_cm=500*PC_TO_CM)
  ```

- **`ExtinctionConfig` / `ExtinctionModel` / `make_extinction_model` —
  `distance_pc` and `scale_distance` removed.**
  Distance dilution belongs exclusively in `forward.py` via `FitParams`.
  Remove any `distance_pc=` or `scale_distance=` kwargs passed to these
  classes.  For standalone distance dilution, apply `(R/d)^2` to the flux
  directly.

- **NumPy runtime requirement raised to `>=2.0`.**
  The package now calls `np.trapezoid` (added in NumPy 2.0).

### New features

- **`ExtinctionModel.with_av(a_v)`** — returns a new `ExtinctionModel`
  identical to the caller except with a different A_V.  Intended for hot-path
  use inside the inverse model; avoids an import inside a tight loop.

- **`run_inverse` — `verbose` parameter (default `True`).**
  Set `verbose=False` to suppress the parameter-summary and extinction
  prints when running batch fits.

- **`run_inverse` — `extinction_law` now warns when `fit_params` is also
  supplied.**  Previously the argument was silently discarded.

- **`AtmosphereGrid.flux_fortran`** — Fortran-contiguous copy of the flux
  cube pre-computed at `load_grid` time.  Eliminates `np.asfortranarray`
  calls from the forward-model hot path.

- **`demos/plot_utils.py`** — shared posterior-visualisation utilities for
  all demo scripts.  Contains two functions:
  - `plot_posterior(result, corner_path, chains_path, truths=None)` — corner
    plot and chain traces for arbitrary free-parameter sets.
  - `plot_sed_posterior(result, fit_params, grid, filters, path, ...)` — SED
    credible-interval band with observed photometric data overlaid, plus a
    filter-transmission sub-panel.

- **`demos/simbad_query_inverse.py`** — new `--extinction-law` CLI flag
  (default `fitzpatrick99`); produces a `_sed.png` posterior SED plot
  alongside the corner and chain plots; `--mh` renamed to `--meta`.

### Performance

- **`_natural_cubic_spline` vectorised.** The scalar Python loop over
  evaluation points is replaced by a fully vectorised array expression.
  Pre-computed M-vectors (`_FM07_SPLINE_M`) eliminate the `np.linalg.solve`
  call from every FM07 extinction evaluation.

- **`interp_filter_onto_sed` (Fortran) two-pointer walk.** Replaces a binary
  search per SED wavelength point with a monotonic two-pointer scan —
  O(n_sed + n_filt) instead of O(n_sed × log n_filt).

- **`find_nearest_point` (Fortran) uses binary search.** Refactored to call
  `find_interval` (already binary) rather than three independent linear scans.

- **Filter and grid arrays no longer copied per call.** `np.asfortranarray`
  and `np.ascontiguousarray` calls were removed from the forward-model loop;
  the arrays are guaranteed to be in the correct layout by the loaders.

### Fixes

- **`gordon23` FUV bump was silently dropped.** Double-indexing `k[m][fuv]
  += ...` modified a temporary copy.  Fixed to accumulate into a local array
  before assigning back to `k[m]`.

- **Fortran `ierr` from bolometric and synthetic-magnitude calls is now
  checked.** Failures are stored as `nan` magnitudes so `np.isfinite` in the
  likelihood function correctly rejects them.

- **`FitParams.__post_init__` no longer mutates a frozen `ParamSpec`.**
  The name-stamping step now uses `dataclasses.replace` instead of
  `object.__setattr__` on a frozen dataclass.

- **`simpson_integration` (Fortran) — even-length arrays.** Instead of
  falling back to full trapezoidal, the routine now applies composite
  Simpson's rule to the first n−1 points (odd) and a single trapezoid panel
  for the last pair — O(h⁴) accurate everywhere except the final panel.

- **`demos/extinction_prescription_comparison.py`** — `posterior_stats`
  percentiles corrected from `15.87`/`84.13` to `15.865`/`84.135`.

### Cleanup

- **`FitParams._spec()` removed.** It was a one-line wrapper around
  `getattr`.  Call sites now use `getattr(self, name)` directly.

- **`FitParams.in_prior()` simplified** to use `unpack()` rather than a
  parallel free-iterator loop.

- **`bounded` is now a proper function** with its own docstring instead of a
  bare module-level alias.

- **`_get_cc_api()` unified** into `sed_model/_cc_ext.py` with a single
  `get_cc_api(required=False)` function used by both `filters.py` and
  `forward.py`.

- **`pandas` import is now lazy** in `grid.py` — imported inside
  `validate_lookup_table` only, not at module load time.

- **`BAD_MAG = -99.9_dp`** named constant added to `colors_def.f90`;
  `M_BOL_SUN`, `L_SUN_CGS`, `SIGMA_SB` (unreferenced) removed.

- **`lin_val` renamed to `interp_val`** in `linear_interp.f90`.

- **Fortran `private/cc_api.f90`** marked clearly as NOT COMPILED (it is a
  legacy reference copy; the compiled file is `fortran/cc_api.f90`).

- **Demo scripts** refactored to share `plot_utils.py`; duplicate
  `plot_posterior_diagnostics` implementations removed from
  `sed_playground_demo.py` and `simbad_query_inverse.py`; `demo_inverse.py`
  corner-plot code replaced by a two-line call.

---

## [0.1.7] — 2026-06-xx

Fortran restructure: mirrored MESA's `colors/private` and `colors/public`
module layout.  See `MIGRATION.md` for the full account of the restructure,
the two bugs it surfaced, and the step-by-step upgrade instructions.
