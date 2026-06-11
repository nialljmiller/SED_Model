# Grids and Filters

## Atmosphere grids

```python
from sed_model import load_grid

grid = load_grid("/path/to/Kurucz2003all/")
```

`load_grid()` expects an SED_Tools-built directory containing `flux_cube.bin`. A co-located `lookup_table.csv` can be checked separately with `sed_model.grid.validate_lookup_table(grid)`; it is not read at load time.

### Binary format of `flux_cube.bin`

| Section | Contents |
|---|---|
| Header | 4 × int32 → `(n_teff, n_logg, n_meta, n_lambda)` |
| Axes | `n_teff + n_logg + n_meta + n_lambda` float64 values, written consecutively with no padding |
| Payload | `n_teff × n_logg × n_meta × n_lambda` float64 values stored in `(W, M, L, T)` order on disk, loaded into `(T, L, M, W)` in memory |

This matches the layout written by SED_Tools and read by the MESA colors module (`colors_utils.f90 :: load_flux_cube`). The loader validates the header dimensions and raises `ValueError` on truncated axis or payload data.

### The `AtmosphereGrid` object

`load_grid` returns an immutable `AtmosphereGrid` exposing:

| Attribute / method | Description |
|---|---|
| `teff_grid`, `logg_grid`, `meta_grid` | grid axes (K, dex, dex) |
| `wavelengths` | wavelength grid in Å |
| `flux` | surface flux cube, shape `(n_T, n_L, n_M, n_W)`, erg/s/cm²/Å |
| `teff_bounds`, `logg_bounds`, `meta_bounds` | `(min, max)` per axis |
| `in_bounds(teff, logg, meta)` | True if the point lies within all three axes |
| `clamp(teff, logg, meta)` | the point clipped to the grid boundary |
| `interp_radius(teff, logg, meta)` | Euclidean distance in normalised parameter space from the nearest grid node — mirrors the `Interp_rad` diagnostic of the MESA colors module |

## Filters

```python
from sed_model import load_filters, load_filters_from_instrument_dir

filters = load_filters([
    "/path/to/filters/Generic/Johnson/B.dat",
    "/path/to/filters/Generic/Johnson/V.dat",
])

filters = load_filters_from_instrument_dir(
    "/path/to/filters/Generic/Johnson",
    vega_sed_path="/path/to/stellar_models/vega_flam.csv",  # optional
)
```

Filter files are two-column (wavelength Å, transmission 0–1); `#` comments are skipped, and both whitespace-separated and CSV files (with or without a header row) are handled. Wavelengths are sorted ascending if needed.

`load_filters_from_instrument_dir` looks for an index file whose name matches the directory's last component (e.g. `Johnson/Johnson`) listing one `.dat` filename per line — the structure the MESA colors module expects. If no index file exists, every `.dat` in the directory is loaded.

### Zero-points

Each `Filter` carries precomputed photon-counting zero-points using the same integrals as MESA `colors/private/synthetic.f90`:

| System | Definition |
|---|---|
| Vega | \( F_{zp} = \int F_{\rm Vega}(\lambda)\, T(\lambda)\, \lambda\, d\lambda \;/\; \int T(\lambda)\, \lambda\, d\lambda \) |
| AB | same integral with \( F_{\rm AB}(\lambda) = 3.631 \times 10^{-20} \, c / \lambda^2 \) erg/s/cm²/Å |
| ST | \( F_{zp} = 3.63 \times 10^{-9} \) erg/s/cm²/Å (flat \( F_\lambda \); the constant cancels in the integral but it is computed explicitly for consistency with the Fortran implementation) |

The Vega zero-point requires a Vega reference SED (CSV: wavelength Å, flux erg/s/cm²/Å) passed via `vega_sed_path`; without it, `vega_zero_point` is `-1.0` and requesting the Vega system raises `ValueError`. Query a zero-point with `filt.zero_point("Vega")`, `"AB"`, or `"ST"`.

!!! note "Trapezoid vs Simpson"
    The Python zero-point integrals use the trapezoid rule (valid for non-uniform wavelength grids); the Fortran runtime uses adaptive Simpson for the in-band flux. The difference is negligible for the dense wavelength grids produced by SED_Tools, and is folded into the tolerances of the forward test suite.
