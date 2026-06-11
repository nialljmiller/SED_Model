# API Usage Script

`demos/api_usage.py` is the minimal end-to-end workflow — load grid and filters, forward model, save CSV outputs, build a small synthetic inverse problem, sample, print and save the posterior, then repeat the forward model with Av = 0.3 extinction. It is deliberately compact and a good template to copy as the starting point for your own pipeline.

Note one difference from the other demos: its data paths are relative (`data/...`) rather than `~/SED_Tools/...` — adjust to your layout.

## Source

```python
--8<-- "demos/api_usage.py"
```
