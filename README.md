<p align="center">
  <img src="https://raw.githubusercontent.com/nialljmiller/SED_Model/main/docs/SED_Model_logo.jpeg" alt="SED_Model logo" width="600"/>
</p>

<h1 align="center">SED_Model</h1>

**Synthetic photometry and stellar parameter inference.**

SED_Model is a Python package for computing observer-ready synthetic spectral energy distributions (SEDs) and broadband magnitudes from stellar atmosphere grids, and for recovering stellar parameters from observed photometry via Bayesian MCMC inference. It is designed to work directly downstream of `SED_Tools` and as a standalone drop-in for workflows that use MESA-style atmosphere grids and filter transmission curves.

Version: `0.4.0`

**Documentation:** https://nialljmiller.github.io/SED_Model/

---

## Overview

SED_Model provides two paired models that share the same parameter language:

- **Forward model** — maps `(Teff, logg, [M/H], R, d)` to an interpolated SED, bolometric flux and magnitude, and synthetic magnitudes in any set of loaded filters.
- **Inverse model** — maps observed magnitudes with uncertainties to a posterior distribution over `(Teff, logg, [M/H])`, and optionally over `Av` and distance, via `emcee` MCMC sampling.

Both directions operate through a shared `FitParams` object. The same parameter specification flows from the forward model to the likelihood evaluator inside the sampler and back out again.

---

## Features

- Bounded Hermite (default), Hermite and linear interpolation over 3D atmosphere grids `(Teff, logg, [M/H])`.
- Photon-counting filter convolution in Vega, AB, and ST photometric systems.
- Precomputed zero-points matching MESA `colors/private/synthetic.f90` conventions.
- Pure-NumPy interstellar dust extinction with multiple supported laws.
- Bayesian stellar-parameter inference via `emcee` with fixed, bounded, and fully free parameter modes.
- Structured result containers with summaries, persistence to CSV and NPZ, and interpolation diagnostics.
- Fortran kernels for performance-critical paths, compiled via `f2py`/Meson and shipped as prebuilt wheels.
- Forward tests validated against MESA colors reference outputs; inverse tests verified by self-consistent synthetic recovery.

---

## Installation

```bash
pip install sed-model
```

Prebuilt wheels (Fortran extension included, runtime bundled) are published for Linux x86_64/aarch64 and macOS arm64/x86_64. Every release also provides a pure-Python wheel, standalone Fortran libraries (no Python needed), and the source distribution — see the [installation guide](https://nialljmiller.github.io/SED_Model/installation/) for choosing between them.

### From source

```bash
git clone https://github.com/nialljmiller/SED_Model.git
cd SED_Model
python -m pip install -e .      # needs gfortran; meson-python builds the extension
```

The extension compiles:

```text
fortran/public/colors_def.f90
fortran/private/{colors_utils,hermite_interp,linear_interp,hermite_interp_bounded,synthetic,bolometric}.f90
fortran/cc_api.f90
```

into `sed_model.cc_api`. `make` is a quick alternative that builds it in place.

### Compatibility

| Platform | `pip install sed-model` | Python | How it was verified |
|---|---|---|---|
| Linux x86_64 (glibc ≥ 2.28) | prebuilt wheel | 3.10 – 3.14 | CI on every push; built and tested in manylinux_2_28 (emulated) |
| Linux aarch64 (glibc ≥ 2.28) | prebuilt wheel | 3.10 – 3.14 | CI on every push; all five Pythons tested on Debian 11/12/13 |
| macOS arm64 (macOS 14+) | prebuilt wheel | 3.10 – 3.14 | CI on every push; build and repair path tested locally |
| macOS x86_64 (macOS 15+) | prebuilt wheel | 3.10 – 3.14 | CI on every push; not tested locally |
| Linux musl (Alpine) | builds from source | 3.10+ | tested; needs `gcc gfortran musl-dev` |
| macOS older than the above | builds from source | 3.10+ | needs `gfortran` |
| Windows | not supported natively | – | use WSL |

The Fortran runtime is bundled in the Linux wheels and linked statically into the macOS wheels, so nothing besides Python and NumPy ≥ 2.0 is needed to use a wheel. "Builds from source" means pip compiles the Fortran on your machine and requires `gfortran` and a C compiler; if either is missing the build stops with instructions.

Each release also provides, on the [GitHub releases page](https://github.com/nialljmiller/SED_Model/releases), a pure-Python wheel (no compiled extension) and standalone Fortran libraries (no Python) for the Linux and macOS platforms above. See the [installation guide](https://nialljmiller.github.io/SED_Model/installation/).

### Requirements

* Python >= 3.10
* `numpy` >= 2.0
* `pandas`
* `emcee`
* `gfortran` (building from source only)

`matplotlib` is only needed for plotting/demo workflows.

---

## Quick Start

```python
from sed_model import load_grid, load_filters, run_forward, run_inverse

grid = load_grid("/path/to/Kurucz2003all/")
filters = load_filters([
    "/path/to/filters/GAIA/G.dat",
    "/path/to/filters/2MASS/J.dat",
])

# Forward model
result = run_forward(
    teff=5777,
    logg=4.44,
    meta=0.0,
    R=6.957e10,
    d=3.086e19,
    grid=grid,
    filters=filters,
    mag_system="Vega",
)

print(result.magnitudes)

# Inverse model
posterior = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10,
    d=3.086e19,
    grid=grid,
    filters=filters,
    mag_system="Vega",
)

posterior.print_summary()
```

---

## Grid Loading

```python
from sed_model import load_grid

grid = load_grid("/path/to/Kurucz2003all/")
```

`load_grid()` expects a MESA/SED_Tools-style atmosphere grid containing a `flux_cube.bin`, with optional validation against a co-located `lookup_table.csv`.

The returned `AtmosphereGrid` exposes:

* `teff_grid`, `logg_grid`, `meta_grid`
* `teff_bounds`, `logg_bounds`, `meta_bounds`
* `wavelengths`
* `in_bounds(teff, logg, meta)`
* `clamp(teff, logg, meta)`
* `nearest_grid_distance(teff, logg, meta)`

---

## Filter Loading

```python
from sed_model import load_filters, load_filters_from_instrument_dir

filters = load_filters([
    "/path/to/filters/Generic/Johnson/B.dat",
    "/path/to/filters/Generic/Johnson/V.dat",
])

filters = load_filters_from_instrument_dir(
    "/path/to/filters/Generic/Johnson"
)
```

Each `Filter` stores:

* `name`
* `wavelengths`
* `transmission`
* Vega, AB, and ST zero-points
* `zero_point(system)`

Supported systems are:

```text
Vega
AB
ST
```

---

## Forward Model

```python
from sed_model import run_forward

result = run_forward(
    teff=5778.0,
    logg=4.44,
    meta=0.0,
    R=6.957e10,
    d=3.086e19,
    grid=grid,
    filters=filters,
    mag_system="Vega",
    interp_method="hermite_bounded",  # or "hermite", "linear"
)
```

`ForwardResult` contains:

* `wavelengths`
* `surface_flux`
* `observed_flux`
* `magnitudes`
* `band_fluxes`
* `bol_flux`
* `bol_mag`
* `nearest_grid_distance`
* `clamped`
* `teff`, `logg`, `meta`, `R`, `d`, `a_v`
* `mag_system`
* `extinction_applied`

---

## Inverse Model

```python
from sed_model import run_inverse

result = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10,
    d=3.086e19,
    grid=grid,
    filters=filters,
    n_walkers=32,
    n_steps=1000,
    n_burn=300,
    seed=42,
)
```

`InverseResult` contains:

* `samples`
* `log_prob`
* `filter_names`
* `obs_magnitudes`
* `obs_uncertainties`
* `R`, `d`, `mag_system`
* `n_walkers`, `n_steps`, `n_burn`, `n_thin`
* `acceptance_fraction`
* `autocorr_time`

It also provides:

```python
result.summary()
result.print_summary()
result.save("posterior.npz")
```

---

## Parameter Modes

`FitParams` and `ParamSpec` define the shared parameter language used by the forward and inverse models.

```python
from sed_model import fit_params_from_grid, fixed, free, bounded

params = fit_params_from_grid(
    grid,
    a_v=(0.0, 2.0),   # free Av
    d_cm=3.086e19,    # fixed distance
)
```

Available modes:

```python
fixed(value)
free(lo, hi)
bounded(lo, hi)
```

Example:

```python
from sed_model import FitParams, fixed, free

params = FitParams(
    teff=free(5000.0, 6500.0),
    logg=fixed(4.44),
    meta=free(-1.0, 0.5),
    a_v=fixed(0.0),
    d=fixed(3.086e19),
)
```

---

## Extinction

```python
from sed_model import make_extinction_model, run_forward

ext = make_extinction_model(
    enabled=True,
    law="fitzpatrick99",
    a_v=0.3,
    r_v=3.1,
)

result = run_forward(
    teff=5778,
    logg=4.44,
    meta=0.0,
    R=6.957e10,
    d=3.086e19,
    grid=grid,
    filters=filters,
    extinction=ext,
)
```

Supported extinction-law keys include:

```text
ccm89
odonnell94
fitzpatrick99
fm07
calzetti00
gordon23
```

---

## Output and I/O

```python
from sed_model import save_sed, save_magnitudes

save_sed(result, "forward_sed.csv")
save_magnitudes(result, "forward_magnitudes.csv")
```

Inverse results can be saved directly:

```python
posterior.save("posterior_samples.npz")
```

---

## Relationship to SED_Tools

SED_Model sits directly downstream of `SED_Tools`.

`SED_Tools` handles acquisition, standardization, filter downloads, grid construction, and SED completion/generation. SED_Model uses the resulting atmosphere grids and filter files for forward synthetic photometry and Bayesian stellar-parameter inference.

A typical workflow is:

```text
SED_Tools  ->  build/download SED grids and filters
SED_Model  ->  generate synthetic magnitudes and infer stellar parameters
```

---

## Documentation and Demos

Full documentation, including walkthroughs of every demo script in `demos/`, lives at:

https://nialljmiller.github.io/SED_Model/

The site is built with MkDocs from the `docs/` directory and embeds the demo sources directly, so the documented examples are always identical to the code in this repository.

---

## Testing

```bash
pytest tests/test_forward.py
pytest tests/test_inverse.py
```

The forward tests verify grid loading, bounds checking, interpolation, filter loading, zero-points, bolometric quantities, and synthetic magnitudes.

The inverse tests verify self-consistent posterior recovery and `InverseResult` persistence.

---

## Repository Structure

```text
sed_model/
  __init__.py
  grid.py
  filters.py
  forward.py
  inverse.py
  params.py
  io.py
  sed_extinction.py

fortran/
  public/colors_def.f90
  private/{colors_utils,hermite_interp,linear_interp,
           hermite_interp_bounded,synthetic,bolometric}.f90
  cc_api.f90

demos/
tests/
Makefile
meson.build
pyproject.toml
meson.build
meson_options.txt
tools/
```

---

## License

SED_Model uses a mixed-license structure. The Python source code, documentation, examples, tests, and packaging files are distributed under the MIT License. The Fortran numerical kernels in `fortran/` are distributed under LGPL-3.0-or-later. See `LICENSE` and source file headers for details.

---

<p align="center">
  <img src="https://raw.githubusercontent.com/nialljmiller/SED_Model/main/docs/SED_Model_Logo_full.png" alt="SED_Model: The Pythonic Inside Story of MCMC-Based Stellar Parameter Inference from Photometric Grids, by Niall Miller" width="360"/>
</p>
