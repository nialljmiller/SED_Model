# Demonstrations

Every page in this section embeds the **actual source file** from the repository's `demos/` directory — the documentation build pulls the scripts in directly, so what you read here is byte-identical to what you run.

All demos follow the same pattern: paths to SED_Tools data products are configured at the top of the file, synthetic "observations" are generated from known true parameters, and the recovery is checked against that truth — so each script doubles as an end-to-end validation you can run on your own grid and filters.

| Demo | What it shows |
|---|---|
| [Forward Model Demo](forward.md) | parameters → SED + synthetic magnitudes, with an annotated SED plot |
| [Inverse Model Demo](inverse.md) | magnitudes → posterior on `(Teff, logg, [M/H])`, with corner and trace plots |
| [Parameter Modes Demo](param-modes.md) | the same data fitted under four fixed/bounded/free configurations |
| [Extinction Comparison](extinction-comparison.md) | one inverse fit per extinction prescription, compared against truth |
| [Playground](playground.md) | a single configurable script for arbitrary forward + inverse experiments |
| [API Usage Script](api-usage.md) | the minimal end-to-end workflow, suitable as a template |

## Requirements

The demos need the compiled Fortran extension (`make`), `matplotlib` (`pip install "sed-model[demo]"`), `emcee` for the inverse demos, and SED_Tools data products on disk. Each script's path block looks like:

```python
GRID_DIR   = "~/SED_Tools/data/stellar_models/Kurucz2003all"
FILTER_DIR = "~/SED_Tools/data/filters/Generic/Johnson"
VEGA_SED   = "~/SED_Tools/data/stellar_models/vega_flam.csv"
```

Edit these to match your machine and run with `python demos/<name>.py`.
