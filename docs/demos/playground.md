# Playground

`demos/sed_playground_demo.py` is a single configurable script for arbitrary forward + inverse experiments. Two dictionaries at the top control everything:

- **`TRUE_STAR`** — what the fake observation really is: atmosphere parameters, radius (R☉), distance (pc), and extinction (law, Av, R_V, Gordon+2023 environment).
- **`FIT_CONFIG`** — what the fitter is allowed to know: each of the five parameters is either `{"mode": "fixed", "value": ...}` or `{"mode": "fit", "bounds": (lo, hi), "start": ...}`.

Distances are given in parsecs in the config and converted to cm internally. The script prints the synthetic observations and a dynamic posterior summary (it adapts to whichever parameters you freed), and saves `sed_playground_sed.png` plus corner/trace diagnostics (`inverse_corner.png`, `inverse_chains.png`) sized to the number of free parameters.

Concepts demonstrated: translating a human-friendly config into [`FitParams`](../guide/parameters.md), keeping the fit's [extinction model](../guide/extinction.md) consistent with a fixed-or-free Av, dynamic posterior plotting.

## Source

```python
--8<-- "demos/sed_playground_demo.py"
```
