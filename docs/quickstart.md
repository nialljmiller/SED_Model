# Quick Start

Everything below assumes you have an SED_Tools-built atmosphere grid and a set of filter files on disk — see [Installation](installation.md#data-products).

## Load a grid and filters

```python
from sed_model import load_grid, load_filters

grid = load_grid("/path/to/Kurucz2003all/")
filters = load_filters([
    "/path/to/filters/GAIA/G.dat",
    "/path/to/filters/2MASS/J.dat",
])
```

To use the **Vega** magnitude system, supply a Vega reference SED when loading filters (otherwise only AB and ST zero-points are computed):

```python
filters = load_filters(
    ["/path/to/filters/GAIA/G.dat"],
    vega_sed_path="/path/to/stellar_models/vega_flam.csv",
)
```

## Forward model: parameters → magnitudes

```python
from sed_model import run_forward

result = run_forward(
    teff=5777,
    logg=4.44,
    meta=0.0,
    R=6.957e10,      # stellar radius in cm (1 R_sun)
    d=3.086e19,      # distance in cm (10 pc)
    grid=grid,
    filters=filters,
    mag_system="AB",
)

print(result.magnitudes)   # {'G': ..., 'J': ...}
print(result.bol_mag)      # bolometric magnitude
```

All lengths are **cgs (cm)**; the convenience constants `RSUN_TO_CM` and `PC_TO_CM` are exported:

```python
from sed_model import RSUN_TO_CM, PC_TO_CM

R = 1.0 * RSUN_TO_CM
d = 10.0 * PC_TO_CM
```

## Inverse model: magnitudes → posterior

```python
from sed_model import run_inverse

posterior = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10,
    d=3.086e19,
    grid=grid,
    filters=filters,
    mag_system="AB",
    n_walkers=32,
    n_steps=1000,
    n_burn=300,
    seed=42,
)

posterior.print_summary()
posterior.save("posterior.npz")
```

By default this samples `Teff`, `logg`, and `[M/H]` across the full grid with flat priors, holding `Av = 0` and the supplied distance fixed. To change which parameters are sampled, pass a [`FitParams`](guide/parameters.md):

```python
from sed_model import fit_params_from_grid

params = fit_params_from_grid(
    grid,
    a_v=(0.0, 2.0),   # Av free, flat prior over [0, 2] mag
    d_cm=3.086e19,    # distance fixed
)

posterior = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10,
    grid=grid,
    filters=filters,
    fit_params=params,
)
```

## Add extinction to the forward model

```python
from sed_model import make_extinction_model

ext = make_extinction_model(
    enabled=True,
    law="fitzpatrick99",
    a_v=0.3,
    r_v=3.1,
)

reddened = run_forward(
    teff=5778, logg=4.44, meta=0.0,
    R=6.957e10, d=3.086e19,
    grid=grid, filters=filters,
    extinction=ext,
)
```

See the [Extinction guide](guide/extinction.md) for the six supported laws and the conventions used.

## Next steps

The [Demonstrations](demos/index.md) section walks through complete, runnable scripts covering each of these workflows end to end, including posterior corner plots and multi-scenario comparisons.
