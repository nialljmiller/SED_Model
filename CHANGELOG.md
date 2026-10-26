# Changelog

## [Unreleased]

Documentation only: `docs/testing.md` now attributes the large grid-edge Hermite
differences from MSG to empty stencil nodes, which the 0.4.2 warning covers.

## [0.4.4] — 2026-10-25

Inverse-side guard rails; the forward model and the sampled posterior for a
given problem are unchanged except as noted in the first item.

- **The fitter rejects Hermite-stencil holes.** With `interp_method='hermite'`
  or `'hermite_bounded'`, `run_inverse` now treats a cell whose Hermite
  stencil touches an empty node (`AtmosphereGrid.cell_has_empty_nodes(...,
  stencil=True)`) like a cell with an empty corner: the posterior is -inf
  there, so neither the walker initialisation nor the chain visits it. The
  per-cell mask is built once, lazily, so the check is an array lookup.
  `interp_method='linear'` and grids without empty nodes are unaffected.
- **Non-convergence warning.** `run_inverse` emits a `UserWarning` when the
  chain ends with `converged False` (shorter than 50 autocorrelation times
  after burn-in), suggesting a larger `n_steps` or `auto_stop=True`. (With
  `auto_stop=True` the existing auto_stop warning is used instead.) Defaults
  and results are unchanged.
- **`fit_jitter` hint.** With `fit_jitter=False` and `chi2_best / n_data > 5`
  a `UserWarning` suggests `fit_jitter=True` (or `sys_floor`): the quoted
  uncertainties then probably under-estimate the model/photometry mismatch.
- Tests: `tests/test_inverse_guards.py`.

## [0.4.3] - 2026-10-25

Documentation only.

- README: the Quick Start now loads a Vega spectrum before asking for the Vega
  system; the inverse-model section, test instructions and repository layout
  describe the current code; added citation and credits sections.
- Documentation: warning behaviour for empty grid cells, the full set of
  `InverseResult` fields, calibration numbers from the 0.4.1 re-test (including
  the effect of a fixed radius), updated cross-code validation results, and a
  new Fortran Kernels page replacing `MIGRATION.md`.
- `grids-and-filters.md` listed a non-existent `interp_radius`; the method is
  `nearest_grid_distance`.
- `demos/api_usage.py` asked for the Vega system without loading a Vega
  spectrum and failed; it now passes `vega_sed_path`.
- `CONTRIBUTING.md`: current test and release workflow.
- `CITATION.cff` lists both licences. Package classifiers list Python
  3.10-3.14 and Linux/macOS instead of "OS Independent".

## [0.4.2] — 2026-10-25

- **The empty-node warning now covers the Hermite stencil.** With
  `interp_method='hermite'` or `'hermite_bounded'`, `run_forward` also warns
  when a node the derivative estimates read is empty: each corner's
  axis-aligned neighbours, one extra node either side of the cell along every
  axis. In the MSG/MIST re-test about one in six cells whose eight corners
  are all populated had such a gap, and Hermite there is off by several
  percent in bolometric flux (up to ~0.08 mag in photometry) with no warning.
  `AtmosphereGrid.cell_has_empty_nodes(..., stencil=True)` exposes the check.
  The warning is advisory: the interpolated values are returned unchanged
  (including the large unphysical magnitudes of a cell with an empty corner).
  `interp_method='linear'` is not affected.

## [0.4.1] — 2026-10-25

- **Filter detector type.** `load_filters` reads `DetectorType=photon|energy`
  from the column-name row that SED_Tools >= 0.4.0 writes
  (`Wavelength,Transmission,...;DetectorType=energy;...`). An energy-counting
  curve is converted to its photon-counting equivalent T/lambda, rescaled to
  its original peak, so that the lambda-weighted band integrals equal
  int(F T) / int(T), as in the MESA colors `load_filter`. A missing or
  unrecognised tag means photon counting (a warning for an unrecognised one).
  The new `Filter.detector_type` records the tag; for `'energy'` filters
  `Filter.transmission` holds the converted curve. Zero points and magnitudes
  both follow. Only the column row is read; leading `# DetectorType:` comment
  lines (SED_Tools 0.3.0) are ignored.
- `WavelengthUnit` on the column row (the unit of that file's wavelength
  column) is honoured: known units are converted to Angstrom with the same
  table as SED_Tools, an unknown unit raises `ValueError`. Files without the
  tag are read as Angstrom, as before.

## [0.4.0] — 2026-10-25

Validated against the MESA `colors` module, MSG and MIST (see
[Testing and Validation](docs/testing.md)); forward parity fixes and a rebuilt
inverse model.

### Forward model

Forward-model parity fixes found by sweeping against MESA `colors` (see
`tests/test_kernels.py`). SEDs and magnitudes otherwise agree with MESA to
rounding error.

- **`bol_mag` is `99.0` when the integrated flux is not positive** (MESA's
  value), instead of `-99.9`.
- **`nearest_grid_distance` now equals MESA's `Interp_rad`**: the distance to
  the nearest tabulated model (empty cube nodes excluded), with axes
  normalised by the range of the tabulated models. Previously it measured to
  every node of the rectangular cube.
- **`run_forward` warns** (`UserWarning`) when the query's grid cell touches a
  node with no model; the interpolated result there is unreliable.
  `AtmosphereGrid.cell_has_empty_nodes` exposes the check.
- **`clamped` is no longer set by `interp_method='linear'`** merely because some
  wavelength bins fell back to the nearest node at a populated, in-bounds point.

### Inverse model

The new fast evaluator calls the same Fortran kernels and agrees with
`run_forward` to 1e-9 mag.

### Changed (default behaviour)

- **`run_inverse` now initialises walkers with a global search + Levenberg-
  Marquardt refinement (`init='auto'`)** instead of a 2 % Gaussian ball at the
  middle of the prior box. For the full Kurucz grid that ball sat in an empty
  grid corner, and walkers that started there often never left it. In an
  inject-and-recover benchmark (16-24 stars, 0.005-0.05 mag noise, Johnson /
  optical+NIR / 11-band sets) the old default gave, with Av free, a Teff error
  MAD of 10.9 % (now 1.6 %) and 12 % of fits with logg off by more than 4 sigma
  (now 0.7 %); with Av and distance free the 95 % Teff coverage went from 0.72
  to 0.91. Giving `p0_centre` (or `init='ball'`) restores the old behaviour
  exactly, including emcee's stretch move.
- With `init='auto'` the default sampler move is a differential-evolution
  mixture (`moves='de'`), which halves the autocorrelation time and lets
  walkers rebalance between modes; `moves='stretch'` restores the old move.
- Grid nodes holding an all-zero SED (unphysical corners of Kurucz grids) now
  give a zero posterior instead of finite magnitudes of ~770.
- The likelihood no longer builds a full `ForwardResult` per call.

### Added

- `run_inverse` options: `init`, `sys_floor`, `fit_jitter` / `jitter_max`,
  `priors` (`normal`, `halfnormal`, `parallax`, `edsd`, or callables),
  `upper_limits`, `auto_stop`, `moves`, `pool`; `n_burn=None` for automatic
  burn-in. Non-finite magnitudes are dropped with a warning.
- `InverseResult`: `converged`, `n_eff`, `chi2_best`, `n_data`, `init_info`
  (saved and loaded with `.npz`), and `credible_interval(name, level, method)`
  with highest-posterior-density intervals.

## [0.3.2] — 2026-10-25

- **Filter files from SED_Tools 0.4.0 load again.** SED_Tools 0.4.0 puts its
  metadata on the column row (`Wavelength,Transmission,FilterID=...;DetectorType=...`);
  the loader now reads only the first two columns instead of failing with
  `expected two-column file`. Required for SED_Tools >= 0.4.0 filter files.
  The metadata itself is not read yet.
- Tests find their data in `~/.sed_tools/data` as well as the older locations
  and `$MESA_DIR/data/colors_data`.

## [0.3.1] — 2026-10-23

Packaging and CI only; no numerical changes.

- README logos use absolute URLs, so they render on the PyPI project page.
- Square `docs/favicon.png` for the documentation site.
- Release job: the tag guard no longer misreads the API error body as a commit
  when the tag does not yet exist.
- CI: actions moved to current majors (Node 24), and runners pinned to
  `ubuntu-24.04` instead of `ubuntu-latest`.
- Removed the legacy, uncompiled `fortran/private/cc_api.f90`.

## [0.3.0] — 2026-10-21

Fortran kernels synchronised with the current MESA `colors` module, and
prebuilt wheels for Linux and macOS published automatically from CI.

> Version 0.2.0 was tagged but never published to PyPI; 0.3.0 supersedes it and
> includes everything listed under 0.2.0 below.

### Changed (numerical behaviour)

- **Bolometric magnitude now uses the IAU 2015 B2 zero-point.**
  `Mag_bol = -2.5 log10(F_bol / F_bol,0)` with
  `F_bol,0 = Lsun 10^(0.4 * 4.74) / (4 pi (10 pc)^2)`, so a solar-luminosity
  source at 10 pc has `bol_mag = 4.74`. Previously `bol_mag = -2.5 log10(F_bol)`
  (about 11.5 mag fainter). `bol_flux` is unchanged.
- **Bolometric integration cleans its input.** Samples with a non-finite or
  non-positive wavelength, or a non-finite or negative flux, are zeroed before
  integrating.
- **Hermite tangents on non-uniform grids** use the spacing-weighted average of
  the one-sided slopes (exact for quadratics) instead of the unweighted central
  difference. Singleton axes are handled as degenerate dimensions. The
  `1e-10` positivity floor is removed from the raw `hermite` interpolant.
- **Linear interpolation** pins out-of-range queries to the boundary cell,
  falls back to the nearest grid node when the interpolant is non-finite or
  non-positive, and uses a `1e-300` floor (previously `1e-10`).
- **Filter transmission** is interpolated linearly onto the SED grid, held at
  the end values, then zeroed strictly outside the tabulated range with
  negative values clamped to zero (MESA `zero_filter_outside_support`). The
  filter's own edge wavelengths are no longer forced to zero. This also applies
  to the AB/ST/Vega zero-points computed in Fortran.
- **Simpson integration** is the non-uniform-spacing rule, with a trapezoid on
  the last interval when the number of points is even.
- **Grid location**: `find_interval` collapses dummy axes (all `0`, `999` or
  `-999`) and clamps out-of-range values to the boundary cell; the
  nearest-node search is a true minimum-distance search.

### Packaging and CI

- **Build backend is now meson-python** (`pyproject.toml`, `meson.build`,
  `meson_options.txt`); `setup.py` and `MANIFEST.in` are removed. The
  version lives in `meson.build` only; `sed_model.__version__` reads the
  installed package metadata.
- **Four deliverables from one source tree**, selected with
  `-Dcomponents={all,python,fortran}`: platform wheels (Fortran built, Python
  packaged), a pure `py3-none-any` wheel (Fortran not built), standalone
  Fortran static/shared libraries with `.mod` files and a Fortran smoke test
  (no Python needed), and the sdist.
- **`.github/workflows/build.yml`** builds and tests all four on every push and
  pull request (Linux x86_64/aarch64 via cibuildwheel and manylinux_2_28,
  macOS arm64/x86_64) and publishes on `v*` tags: platform wheels and sdist to
  PyPI, everything to the GitHub release. It replaces `publish.yml`, which only
  uploaded the sdist. The pure wheel is intentionally kept off PyPI.
- The missing-extension `ImportError` now explains the install options.
- **Fortran runtime linking.** macOS wheels link the Fortran runtime
  statically (no bundled dylibs, minimum macOS set by the build target, not
  by the Homebrew gfortran bottle); Linux wheels bundle `libgfortran` with
  auditwheel. The extension is linked by the Fortran driver
  (`link_language: 'fortran'`). The standalone Fortran libraries are
  self-contained on both: statically linked runtime on macOS, and on Linux
  `libgfortran`/`libquadmath` vendored beside the library with an `$ORIGIN`
  rpath (the distro `libgfortran.a` is not position-independent on x86_64).
- **Clearer source-build failures.** A missing C or Fortran compiler now stops
  meson with instructions (per platform) instead of a bare "unknown compiler",
  and building the pure-Python wheel needs no compiler at all.

### New features

- **`interp_method="hermite_bounded"`** (new default for `run_forward`,
  `run_forward_batch`, `run_inverse`), backed by the new
  `cc_api.interp_sed_hermite_bounded` and
  `fortran/private/hermite_interp_bounded.f90`. Negligible negative Hermite
  undershoot is clipped to zero; an unusable Hermite SED falls back to
  trilinear interpolation. If both fail, `run_forward` raises `ValueError`.
  Pass `interp_method="hermite"` for the previous behaviour.
- `tests/test_kernels.py`: data-free kernel tests reproducing the MESA colors
  unit-test cases (filter support, non-uniform Hermite, bounded repair and
  fallback, bolometric zero-point).

### Fixed

- **Source builds on musl/Alpine.** The kernels used the `ieee_arithmetic`
  intrinsic module, which gfortran on musl does not provide; they now use a
  portable elemental `is_finite` in `colors_utils`.
- Type-checking imports in `forward.py` / `inverse.py` pointed at the
  top-level `sed_extinction` instead of `.sed_extinction`.
- `ForwardResult.clamped` is no longer set by interpolation fallbacks inside
  the grid when `hermite_bounded` is used; it reflects out-of-grid queries.

### Docs and tests

- README, docs and build lists updated (version, Fortran source list,
  `nearest_grid_distance`, interpolation schemes).
- Group 1 forward tests now compare against the current MESA
  `test_output` to 1e-5, including the new `Mag_bol` values.

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
