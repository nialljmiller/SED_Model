# Parameter Modes Demo

`demos/demo_param_modes.py` fits the **same** synthetic observations (a solar-type star at 500 pc with injected Av = 0.3) under four sampling configurations:

- **A — All open**: `Teff`/`logg`/`[M/H]` free across the full grid; `Av` and distance fixed.
- **B — Teff bounded, logg fixed**: `Teff` free within ±500 K of a prior guess; `logg` pinned to a catalogue value; `[M/H]` open.
- **C — All bounded + Av free**: physically motivated windows on the atmospheric parameters, `Av` free over [0, 2] mag — the Av–Teff correlation appears.
- **D — Everything free**: all five parameters sampled, including distance — the [Av–distance degeneracy](../guide/parameters.md#the-avdistance-degeneracy) is on full display.

The output figure `param_modes_demo.png` shows one row per scenario: 1D marginals for each free parameter, with the truth in red and the prior window shaded.

Concepts demonstrated: [`FitParams` and `ParamSpec`](../guide/parameters.md), `fit_params_from_grid`, mixing fixed/bounded/free modes, per-scenario walker counts.

## Source

```python
--8<-- "demos/demo_param_modes.py"
```
