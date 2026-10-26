<p align="center">
  <img src="https://raw.githubusercontent.com/nialljmiller/SED_Model/main/docs/SED_Model_logo.jpeg" alt="SED_Model logo" width="600"/>
</p>

<h1 align="center">SED_Model</h1>

**Synthetic photometry and stellar parameter inference.**

SED_Model is a Python package for computing observer-ready synthetic spectral energy distributions (SEDs) and broadband magnitudes from stellar atmosphere grids, and for recovering stellar parameters from observed photometry via Bayesian MCMC inference. It is designed to work directly downstream of `SED_Tools` and as a standalone drop-in for workflows that use MESA-style atmosphere grids and filter transmission curves.

Version: `0.4.4`

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
- Reads the `DetectorType` (photon or energy counting) and `WavelengthUnit` tags that SED_Tools writes on each filter file.
- Warns when a query falls in a part of the grid that has no model (see [Forward Model](#forward-model)).
- Forward model checked against the MESA `colors` module (magnitudes agree to double-precision rounding), MSG and MIST; inverse model benchmarked by injection and recovery. See [Testing](#testing).

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
filters = load_filters(
    [
        "/path/to/filters/GAIA/G.dat",
        "/path/to/filters/2MASS/J.dat",
    ],
    vega_sed_path="/path/to/stellar_models/vega_flam.csv",  # needed for mag_system="Vega"
)

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

`load_grid()` expects a MESA/SED_Tools-style atmosphere grid directory containing a `flux_cube.bin`. A co-located `lookup_table.csv` can be checked with `sed_model.grid.validate_lookup_table(grid)`; it is not read at load time.

The returned `AtmosphereGrid` exposes:

* `teff_grid`, `logg_grid`, `meta_grid`
* `teff_bounds`, `logg_bounds`, `meta_bounds`
* `wavelengths`
* `in_bounds(teff, logg, meta)`
* `clamp(teff, logg, meta)`
* `nearest_grid_distance(teff, logg, meta)` (same quantity as MESA's `Interp_rad`)
* `cell_has_empty_nodes(teff, logg, meta, stencil=False)`

---

## Filter Loading

```python
from sed_model import load_filters, load_filters_from_instrument_dir

filters = load_filters([
    "/path/to/filters/Generic/Johnson/B.dat",
    "/path/to/filters/Generic/Johnson/V.dat",
])

filters = load_filters_from_instrument_dir(
    "/path/to/filters/Generic/Johnson",
    vega_sed_path="/path/to/stellar_models/vega_flam.csv",  # optional
)
```

Each `Filter` stores:

* `name`
* `wavelengths` (Å)
* `transmission`
* `detector_type` (`"photon"` or `"energy"`, from the file's `DetectorType` tag)
* Vega, AB, and ST zero-points
* `zero_point(system)`

The Vega system needs a Vega reference SED (`vega_sed_path`); without it only AB and ST are available. An energy-counting curve is converted to its photon-counting equivalent when it is loaded; see [Grids and Filters](https://nialljmiller.github.io/SED_Model/guide/grids-and-filters/#detector-type).

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

If the query lies in a grid cell that touches a node with no model, `run_forward` emits a `UserWarning`. Kurucz-type grids are stored as a full rectangular cube in which a large fraction of the nodes are empty (about 40 % for `Kurucz2003all`). With the Hermite methods the warning also covers the extra nodes the derivative estimates read around the cell. The result is still returned, but it is unreliable there: in a cell with an empty corner the magnitudes can be wildly wrong, and next to an empty node Hermite can overshoot by several percent in flux. `linear` does not read the extra nodes. `grid.cell_has_empty_nodes(...)` runs the same check without a forward call.

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
    n_steps=2000,
    n_burn=500,
    seed=42,
)
```

By default the walkers are started from a global search plus Levenberg-Marquardt refinement (`init="auto"`) and moved with a differential-evolution mixture, rather than from a ball at the middle of the prior box. Optional extras: `sys_floor` and `fit_jitter` (systematic error terms), `priors` (e.g. a parallax prior on distance), and `upper_limits`. See the [Inverse Model guide](https://nialljmiller.github.io/SED_Model/guide/inverse/).

`InverseResult` contains:

* `samples`, `log_prob`, `param_names`, `fixed_params`
* `filter_names`, `obs_magnitudes`, `obs_uncertainties`
* `R`, `d`, `mag_system`
* `n_walkers`, `n_steps`, `n_burn`, `n_thin`
* `acceptance_fraction`
* `autocorr_time`, `converged`, `n_eff`
* `chi2_best`, `n_data`, `init_info`

It also provides:

```python
result.summary()
result.print_summary()
result.map_estimate()
result.credible_interval("teff", 0.68, method="hpd")
result.save("posterior.npz")
```

Intervals are only as good as the model. With real photometry the atmosphere library, filter curves and zero points add systematic errors that the likelihood does not know about; use `fit_jitter=True` or `sys_floor`, and check `converged`.

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
pytest tests/test_kernels.py     # needs no data
pytest tests/test_forward.py     # needs a Kurucz grid and filters
pytest tests/test_inverse.py     # needs a Kurucz grid and filters
```

Run these from outside the source tree against an installed build: inside it, `import sed_model` finds the uncompiled sources and every test is skipped. Tests that need data look in `~/.sed_tools/data`, `~/SED_Tools/data` and `$MESA_DIR/data/colors_data`.

* `test_kernels.py` covers the Fortran kernels and filter handling (MESA colors unit-test cases, detector type, wavelength units, empty-grid warnings) without any data files.
* `test_forward.py` compares magnitudes and bolometric quantities against reference values from the MESA colors unit test.
* `test_inverse.py` injects synthetic photometry and checks that the posterior recovers the truth, and covers `InverseResult` persistence and input validation.

For the comparison with MESA, MSG and MIST, and what to expect when comparing with other codes, see [Testing and Validation](https://nialljmiller.github.io/SED_Model/testing/).

---

## Repository Structure

```text
sed_model/
  __init__.py
  grid.py            atmosphere grids
  filters.py         filter files, detector type, zero-points
  forward.py         forward model
  inverse.py         inverse model (run_inverse)
  _inv_core.py       fast likelihood and initialisation for the inverse model
  params.py          FitParams / ParamSpec
  io.py              result containers and CSV / NPZ output
  sed_extinction.py  extinction laws
  _cc_ext.py         loader for the compiled extension

fortran/             kernels, built into sed_model.cc_api (see docs/fortran.md)
demos/               example scripts, embedded in the documentation
docs/                MkDocs sources
tests/
tools/               build helpers used by CI
Makefile             quick in-place extension build
meson.build, meson_options.txt, pyproject.toml
```

---

## Citing

If you use SED_Model in published work, please cite it using [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button reads it).

---

## Credits

* The Fortran kernels (interpolation, integration, zero-points, bolometric quantities) are ported from the `colors` module of [MESA](https://docs.mesastar.org/), and SED_Model is tested against it. They keep MESA's licence (LGPL-3.0-or-later).
* [SED_Tools](https://github.com/nialljmiller/SED_Tools) produces the atmosphere grids and filter files SED_Model reads.
* Sampling is done with [emcee](https://emcee.readthedocs.io/) (Foreman-Mackey et al. 2013, PASP 125, 306).
* The extinction laws are re-implemented from the original papers: Cardelli, Clayton & Mathis (1989), O'Donnell (1994), Fitzpatrick (1999), Fitzpatrick & Massa (2007), Calzetti et al. (2000) and Gordon et al. (2023). Full references are in the [Extinction guide](https://nialljmiller.github.io/SED_Model/guide/extinction/).
* The forward model was cross-checked against [MSG](https://github.com/rhdtownsend/msg) and the MIST bolometric-correction tables (Choi et al. 2016, ApJ 823, 102; Dotter 2016, ApJS 222, 8).

---

## License

SED_Model uses a mixed-license structure. The Python source code, documentation, examples, tests, and packaging files are distributed under the MIT License. The Fortran numerical kernels in `fortran/` are distributed under LGPL-3.0-or-later. See `LICENSE` and source file headers for details.

---

<p align="center">
  <img src="https://raw.githubusercontent.com/nialljmiller/SED_Model/main/docs/SED_Model_Logo_full.png" alt="SED_Model: The Pythonic Inside Story of MCMC-Based Stellar Parameter Inference from Photometric Grids, by Niall Miller" width="360"/>
</p>
