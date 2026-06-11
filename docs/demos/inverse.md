# Inverse Model Demo

`demos/demo_inverse.py` synthesises observations from known solar parameters (with 0.02 mag Gaussian noise per band), runs the MCMC inverse model on them, prints the posterior summary, saves the chain to `inverse_posterior.npz`, and produces two diagnostic figures: a corner-style plot (`inverse_corner.png`, with the truth overplotted in red) and per-parameter chain traces (`inverse_chains.png`).

To fit real data instead, set `SYNTHESISE = False` and fill in the `OBS_*` arrays — the rest of the script is unchanged.

Concepts demonstrated: [the inverse model](../guide/inverse.md), walker initialisation, posterior [persistence and summaries](../guide/results-io.md).

## Source

```python
--8<-- "demos/demo_inverse.py"
```
