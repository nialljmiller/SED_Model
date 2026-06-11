# Parameter Modes

`sed_model.params` is the single place that defines what a "stellar parameter" is, whether it is sampled or held fixed, and what its physical bounds are. Both the forward and inverse models speak this language, which is what makes the package genuinely bidirectional.

## The five parameters

Every fit involves the same five physical quantities, each described by a `ParamSpec`:

| Name | Quantity | Unit | Default in `fit_params_from_grid` |
|---|---|---|---|
| `teff` | effective temperature | K | free over the grid's Teff range |
| `logg` | log surface gravity | dex (log₁₀ g / cm s⁻²) | free over the grid's logg range |
| `meta` | metallicity [M/H] | dex | free over the grid's [M/H] range |
| `a_v` | V-band extinction | mag | fixed at 0.0 |
| `d` | distance | **cm** | fixed at 1 pc |

## Modes

```python
from sed_model import fixed, free, bounded

fixed(value)      # not sampled; passed unchanged to the forward model
free(lo, hi)      # sampled by MCMC, flat prior over [lo, hi]
bounded(lo, hi)   # alias for free — explicit about hard limits
```

Validation is enforced at construction: free specs require `lo < hi`; distance must be positive; `Av` must be non-negative.

## Building a `FitParams`

### From a grid, with overrides

```python
from sed_model import fit_params_from_grid
from sed_model.params import PC_TO_CM

# Defaults: Teff/logg/[M/H] free over the grid, Av fixed 0, d fixed 1 pc
params = fit_params_from_grid(grid)

# Pass a float to fix a parameter, a (lo, hi) tuple to sample it:
params = fit_params_from_grid(
    grid,
    teff=(5200, 6400),
    logg=4.4,
    meta=0.0,
    a_v=(0.0, 1.5),
    d_cm=500 * PC_TO_CM,
)
```

### Explicitly

```python
from sed_model import FitParams, fixed, free
from sed_model.params import PC_TO_CM

params = FitParams(
    teff = free(4000, 8000),
    logg = free(3.0, 5.0),
    meta = free(-1.0, 0.5),
    a_v  = fixed(0.3),
    d    = fixed(500 * PC_TO_CM),
)
```

## The theta vector

The MCMC samples only the **free** parameters, in the canonical order `teff, logg, meta, a_v, d` (skipping fixed ones). `FitParams` provides the conversion both ways:

```python
params.free_names          # e.g. ['teff', 'meta', 'a_v']
params.n_free              # 3
params.fixed_names         # ['logg', 'd']

theta = params.pack(teff=5800, logg=4.4, meta=0.0, a_v=0.5)
full  = params.unpack(theta)   # {'teff': 5800, 'logg': 4.4, ...} with fixed filled in
params.in_prior(theta)         # True/False
print(params.summary())        # human-readable fixed/free table
```

`initial_ball(n_walkers, centre=..., scatter=...)` generates the starting walker positions used by `run_inverse`.

## The Av–distance degeneracy

If `Av` and `d` are both free simultaneously the posterior **will** be correlated — the data genuinely cannot break the degeneracy without additional information, and the MCMC will show it honestly rather than hiding it. The [Parameter Modes demo](../demos/param-modes.md) makes this visible by running the same observations through four sampling configurations, from everything-fixed-but-atmosphere up to everything-free.
