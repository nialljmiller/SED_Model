# Testing and Validation

```bash
pytest tests/test_forward.py
pytest tests/test_inverse.py
pytest tests/test_kernels.py
```

Run them from outside the source tree against an installed build; inside it `import sed_model` finds the uncompiled sources, has no extension, and every test skips.

All three suites skip, with an explanatory reason, when the Fortran extension is not built or the required data products are not on disk, so `pytest` is always safe to run. `test_kernels.py` needs no data files.

## Data locations

The tests look for SED_Tools data products in `~/.sed_tools/data`, then `~/SED_Tools/data`, then `$MESA_DIR/data/colors_data`, using the first location that exists. Each path can be overridden with an environment variable:

| Variable | Path under the data root |
|---|---|
| `CC_KURUCZ_DIR` | `stellar_models/Kurucz2003all` |
| `CC_JOHNSON_DIR` | `filters/Generic/Johnson` |
| `CC_VEGA_SED` | `stellar_models/vega_flam.csv` |

## Forward model: validated against MESA

`tests/test_forward.py` checks grid loading, bounds/clamping, both interpolation schemes, filter loading and zero-points, and then the headline validation: synthetic magnitudes and bolometric quantities are compared against **reference values produced by the MESA colors unit test** (`colors/test/test_output`, MESA Kurucz2003 grid) for three stars spanning the grid — solar (5778 K / 4.44 / 0.0), a hot main-sequence star (15000 K / 4.0, 5 R☉), and a cool giant (4000 K / 2.0, 20 R☉).

Tolerances are `TOL_MAG = 0.05` mag and `TOL_FLUX = 2%` relative — deliberately generous to absorb grid version differences and the trapezoid (Python zero-points) vs adaptive Simpson (Fortran in-band flux) integration difference, while still catching any real regression in the pipeline.

## Inverse model: self-consistent recovery

`tests/test_inverse.py` does not use external references. Instead it closes the loop on the package itself:

1. run the forward model at known parameters to produce synthetic magnitudes,
2. perturb them with a small known noise level (0.01 mag),
3. run the inverse model and require the posterior median to recover the truth within physically motivated tolerances (≲ 2 grid spacings: 500 K in Teff, 0.5 dex in logg and [M/H]).

Chains are kept deliberately short (hundreds of steps, 16 walkers) so the suite finishes quickly; the tolerances account for the resulting shot noise. The suite also covers the `InverseResult` API — summary structure and ordering, finite MAP estimates, `.npz` save/load round-trips, thinning behaviour — and the input-validation errors (mismatched array lengths, non-positive uncertainties, unknown filter names).

## Cross-code validation

The forward model was swept against the MESA `colors` module, MSG and MIST for 0.4.0 and swept again against the released 0.4.1 wheel; the inverse model was benchmarked by injection and recovery. What a user should expect:

| Against | Result |
|---|---|
| **MESA colors** (standalone driver built from the MESA sources; Kurucz2003all, its alpha-enhanced variant and blackbodies; Johnson, 2MASS, Gaia, TESS, JWST and Argus filters; AB, Vega and ST; all three interpolation methods; on-node, off-node, cell-edge, corner and out-of-range queries; about 21 000 queries per method) | SEDs, bolometric magnitudes and wavelengths bit-identical; band magnitudes within 3e-13 mag (double-precision rounding). Filters tagged `DetectorType=energy` agree with MESA's photon-counting-filters branch to about 1e-13 mag. |
| **MSG** (same atmospheres, built into an MSG grid from the SED_Model cube) | `linear` agrees to about 1e-7 relative in the spectra. With a complete Hermite stencil, Hermite photometry agrees to better than 6 mmag. Linear photometry agrees to 1-4 mmag; U and B reach 13 mmag from quadrature on the coarse Johnson files. Vega and AB zero points agree with an independent calculation to better than 1e-5 mag. |
| **MIST v2.5** bolometric corrections (C3K atmospheres, a different library) | Median offsets of 0.005-0.02 mag in V, I, J, Ks and Gaia, with a scatter of about 0.01 mag; about 0.06 mag in U. These are model differences, not code errors. The offset grows from -0.006 mag at 9 000 K to about -0.04 mag at 25 000 K in V and is nearly colour-neutral (B-V and V-I change by less than 0.003 mag). |

Things that matter when comparing with other work:

- **Weighting.** `synthetic_magnitude` counts photons, as MESA does. Energy weighting would change Gaia G by up to 0.25 mag for cool stars, for curves defined for energy-counting detectors (e.g. Bessell UBVRI). From 0.4.1 filters tagged `DetectorType=energy` are converted accordingly; see [Grids and Filters](guide/grids-and-filters.md#detector-type).
- **Filter definition.** The shipped `Generic/Johnson` curves differ from the Bessell (1990) curves MIST uses: mean wavelength 7025 vs 6650 Å in R and 8871 vs 8086 Å in I. Against MIST this gives residuals up to +0.15 mag in I for cool stars and -0.09 mag for hot ones. Use Bessell-type curves when comparing with MIST in these bands.
- **2MASS Vega zero points** follow Blanton and Roweis, about 0.025 mag from MIST's Cohen-based values.
- **Empty grid cells.** A rectangular grid cube can contain nodes with no model (for Kurucz2003all about 40 % of the cube). `run_forward` emits a `UserWarning` when the query's cell touches one (for the Hermite methods, also when a node in the derivative stencil around the cell is empty), and the result there is unreliable; the values are still returned unchanged. In a cell with an empty corner the magnitudes come out near +770. With a populated cell but an empty stencil node, Hermite is off by a median 5.6 % in bolometric flux (up to 14 %, up to about 80 mmag in photometry). On Kurucz2003all 455 of the 2 723 cells with populated corners have such a gap. This is MESA's behaviour and is reproduced exactly. Before 0.4.2 only the empty-corner case warned. MIST points that fall in a warned cell should be masked before comparing.
- **Grid edges.** Hermite photometry in grid-edge cells also differs from MSG by up to about 80 mmag. This has not been traced to empty nodes, and the warning does not claim to cover it.

Inverse recovery on the Kurucz grid (16 stars per cell, 0.005-0.02 mag noise): the old mid-box start left walkers stuck in an empty corner and chains unconverged. With the 0.4.x defaults, fitting Teff, logg, [M/H] and Av gives a Teff bias near 0.1 %, a scatter of 1.2 % and Teff 95 % intervals that cover the truth 97 % of the time; with distance also free the scatter is 4.9 % and coverage 90 %. Photometry from a different atmosphere library (the MIST test) is far outside the nominal intervals unless `fit_jitter` is used, which is recommended for real photometry. See [Inverse Model](guide/inverse.md#diagnostics) for the full table and for what does not work (three bands with Av and distance free).

The cross-code harnesses are not part of the repository; the in-repo regression tests (`tests/test_kernels.py`) encode the MESA cases that need no data, and `tests/test_reference.py` replays a frozen subset (`tests/data/reference_points.npz`: 33 MESA queries to 1e-10 mag, 24 MSG points at the tolerances above; skips if the extension or Kurucz data is absent).

## What this means for the documentation

Every numerical claim in these docs traces to one of three places: the source code (embedded directly in the [demo pages](demos/index.md)), the docstrings (rendered in the [API Reference](api.md)), or the test suite described here. Nothing is quoted from memory.
