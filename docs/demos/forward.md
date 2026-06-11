# Forward Model Demo

`demos/demo_forward.py` runs the forward model for solar parameters at 10 pc, prints the bolometric flux/magnitude and per-filter synthetic magnitudes, saves the SED and magnitudes to CSV (`forward_sed.csv`, `forward_magnitudes.csv`), and produces `forward_demo.png` — the observed SED on log–log axes with each filter's pivot wavelength marked.

Concepts demonstrated: [grid and filter loading](../guide/grids-and-filters.md), the classic `run_forward` convention, interpolation diagnostics (`interp_radius`, `clamped`), and [CSV output](../guide/results-io.md).

## Source

```python
--8<-- "demos/demo_forward.py"
```
