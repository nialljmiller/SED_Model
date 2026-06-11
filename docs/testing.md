# Testing and Validation

```bash
pytest tests/test_forward.py
pytest tests/test_inverse.py
```

Both suites skip gracefully — with an explanatory reason — when the Fortran extension is not built or the required data products are not on disk, so `pytest` is always safe to run.

## Data locations

The tests look for SED_Tools data products at `~/SED_Tools/data/...` by default, overridable via environment variables:

| Variable | Default |
|---|---|
| `CC_KURUCZ_DIR` | `~/SED_Tools/data/stellar_models/Kurucz2003all` |
| `CC_JOHNSON_DIR` | `~/SED_Tools/data/filters/Generic/Johnson` |
| `CC_VEGA_SED` | `~/SED_Tools/data/stellar_models/vega_flam.csv` |

## Forward model: validated against MESA

`tests/test_forward.py` checks grid loading, bounds/clamping, both interpolation schemes, filter loading and zero-points, and then the headline validation: synthetic magnitudes and bolometric quantities are compared against **reference values produced by the MESA colors unit test** (`colors/test/test_output`, MESA Kurucz2003 grid) for three stars spanning the grid — solar (5778 K / 4.44 / 0.0), a hot main-sequence star (15000 K / 4.0, 5 R☉), and a cool giant (4000 K / 2.0, 20 R☉).

Tolerances are `TOL_MAG = 0.05` mag and `TOL_FLUX = 2%` relative — deliberately generous to absorb grid version differences and the trapezoid (Python zero-points) vs adaptive Simpson (Fortran in-band flux) integration difference, while still catching any real regression in the pipeline.

## Inverse model: self-consistent recovery

`tests/test_inverse.py` does not use external references. Instead it closes the loop on the package itself:

1. run the forward model at known parameters to produce synthetic magnitudes,
2. perturb them with a small known noise level (0.01 mag),
3. run the inverse model and require the posterior median to recover the truth within physically motivated tolerances (≲ 2 grid spacings: 500 K in Teff, 0.5 dex in logg and [M/H]).

Chains are kept deliberately short (hundreds of steps, 16 walkers) so the suite finishes quickly; the tolerances account for the resulting shot noise. The suite also covers the `InverseResult` API — summary structure and ordering, finite MAP estimates, `.npz` save/load round-trips, thinning behaviour — and the input-validation errors (mismatched array lengths, non-positive uncertainties, unknown filter names).

## What this means for the documentation

Every numerical claim in these docs traces to one of three places: the source code (embedded directly in the [demo pages](demos/index.md)), the docstrings (rendered in the [API Reference](api.md)), or the test suite described here. Nothing is quoted from memory.
