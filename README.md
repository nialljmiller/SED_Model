<p align="center">
  <img src="https://raw.githubusercontent.com/nialljmiller/SED_Model/main/docs/SED_Model_logo.jpeg" alt="SED_Model logo" width="600"/>
</p>

<h1 align="center">SED_Model</h1>

SED_Model computes synthetic spectral energy distributions and broadband magnitudes from stellar atmosphere grids, and infers stellar parameters from observed magnitudes with MCMC. It is a spin-off of the `colors` module of [MESA](https://docs.mesastar.org/) and reproduces its synthetic photometry. It reads the grids and filter files written by [SED_Tools](https://github.com/nialljmiller/SED_Tools).

Version `0.4.6` · [Documentation](https://nialljmiller.github.io/SED_Model/) · [Changelog](CHANGELOG.md)

---

## What it does

- **Forward model:** `(Teff, logg, [M/H], R, d)` to an interpolated SED, bolometric flux and magnitude, and synthetic magnitudes in the loaded filters (Vega, AB or ST system), with optional dust extinction.
- **Inverse model:** observed magnitudes with uncertainties to a posterior over `(Teff, logg, [M/H])`, and optionally `Av` and distance, sampled with [emcee](https://emcee.readthedocs.io/).

Both use the same parameter object (`FitParams`). The interpolation, zero points and band integrals are MESA's; the forward model agrees with MESA to double-precision rounding (see [Testing](#testing)).

---

## Installation

```bash
pip install sed-model
```

Requires Python 3.10 or newer and NumPy 2.0 or newer. `pandas` and `emcee` are installed with it; `matplotlib` is only needed for the demos.

| Platform | Install | Python |
|---|---|---|
| Linux x86_64 and aarch64 (glibc 2.28 or newer) | prebuilt wheel | 3.10 – 3.14 |
| macOS arm64 (14 or newer) and x86_64 (15 or newer) | prebuilt wheel | 3.10 – 3.14 |
| Linux musl (Alpine), older macOS | builds from source; needs `gfortran` and a C compiler | 3.10+ |
| Windows | not supported natively; use WSL | – |

The wheels include the Fortran extension and its runtime. Each release also has a pure-Python wheel (no extension) and standalone Fortran libraries on the [releases page](https://github.com/nialljmiller/SED_Model/releases); see the [installation guide](https://nialljmiller.github.io/SED_Model/installation/).

### From source

```bash
git clone https://github.com/nialljmiller/SED_Model.git
cd SED_Model
python -m pip install -e .      # needs gfortran; meson-python builds the extension
```

`make` builds the extension in place as a quicker alternative.

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

`R` and `d` are in cm. The grid and filter files come from SED_Tools.

---

## Inputs

### Grid

`load_grid(path)` reads a SED_Tools or MESA-style atmosphere directory containing `flux_cube.bin`. The returned `AtmosphereGrid` has the axes (`teff_grid`, `logg_grid`, `meta_grid`), their bounds, `wavelengths`, and these methods:

- `in_bounds(teff, logg, meta)` and `clamp(teff, logg, meta)`
- `nearest_grid_distance(teff, logg, meta)`: the same quantity as MESA's `Interp_rad`
- `cell_has_empty_nodes(teff, logg, meta, stencil=False)`: see [Empty grid cells](#empty-grid-cells)

`sed_model.grid.validate_lookup_table(grid)` checks a co-located `lookup_table.csv`; it is not read at load time.

### Filters

`load_filters([paths], vega_sed_path=...)` loads individual files; `load_filters_from_instrument_dir(directory, vega_sed_path=...)` loads a whole instrument directory. Each `Filter` holds `name`, `wavelengths` (Å), `transmission`, `detector_type`, and the Vega, AB and ST zero points.

- The Vega system needs `vega_sed_path`. Without it only AB and ST are available, and asking for Vega raises `ValueError`.
- The `DetectorType` (`photon` or `energy`) and `WavelengthUnit` tags are read from the column-name row of the file, as SED_Tools writes them. Energy-counting curves are converted to photon-counting when loaded. Details are in [Grids and Filters](https://nialljmiller.github.io/SED_Model/guide/grids-and-filters/#detector-type).

---

## Forward model

```python
result = run_forward(
    teff=5778.0, logg=4.44, meta=0.0, R=6.957e10, d=3.086e19,
    grid=grid, filters=filters,
    mag_system="Vega",
    interp_method="hermite_bounded",  # default; or "hermite", "linear"
)
```

`ForwardResult` has `wavelengths`, `surface_flux`, `observed_flux`, `magnitudes`, `band_fluxes`, `bol_flux`, `bol_mag`, `nearest_grid_distance`, `clamped`, the input parameters, `mag_system` and `extinction_applied`. `save_sed(result, path)` and `save_magnitudes(result, path)` write CSV.

### Empty grid cells

Kurucz-type grids are stored as a full rectangular cube in which many nodes have no model (about 40 % of `Kurucz2003all`). If the query lies in a cell that touches such a node, `run_forward` emits a `UserWarning` and still returns the result. For the Hermite methods the warning also covers the extra nodes their derivative estimates read around the cell; `linear` does not read them.

The values in a warned cell are unreliable. With an empty corner the magnitudes can be wildly wrong; with only an empty neighbour node, Hermite can be off by several percent in flux. `grid.cell_has_empty_nodes(...)` runs the same check without a forward call.

### Extinction

```python
from sed_model import make_extinction_model

ext = make_extinction_model(enabled=True, law="fitzpatrick99", a_v=0.3, r_v=3.1)
result = run_forward(..., extinction=ext)
```

Laws: `ccm89`, `odonnell94`, `fitzpatrick99` (R_V 2 to 6), `gordon23` (Milky Way, R_V 2.3 to 5.6), `fm07` (R_V = 3.1), `g03_smcbar` (R_V = 2.74), `g03_lmcavg` (R_V = 3.41) and `calzetti00`. They follow the MESA `colors` extinction module. An R_V outside a law's calibrated range raises `ValueError`; outside the fitted wavelength range k is held constant shortward and continues its functional form longward. Extinction reddens the filter magnitudes only: `bol_flux` and `bol_mag` are those of the intrinsic SED and the zero points are not reddened.

---

## Inverse model

```python
result = run_inverse(
    obs_magnitudes=[5.03, 4.17],
    obs_uncertainties=[0.01, 0.02],
    filter_names=["G", "J"],
    R=6.957e10, d=3.086e19,
    grid=grid, filters=filters,
    n_walkers=32, n_steps=2000, n_burn=500, seed=42,
)
```

Walkers start from a global search plus Levenberg-Marquardt refinement (`init="auto"`) and move with a differential-evolution mixture. Options: `sys_floor` and `fit_jitter` (systematic error terms), `priors` (for example a parallax prior on distance), `upper_limits`, `auto_stop`. See the [Inverse Model guide](https://nialljmiller.github.io/SED_Model/guide/inverse/).

`InverseResult` has the chain (`samples`, `log_prob`, `param_names`, `fixed_params`), the inputs, the run settings, and the diagnostics `acceptance_fraction`, `autocorr_time`, `converged`, `n_eff`, `chi2_best`, `n_data` and `init_info`. Methods:

```python
result.summary()
result.print_summary()
result.map_estimate()
result.credible_interval("teff", 0.68, method="hpd")
result.save("posterior.npz")
```

Things to check before trusting a fit:

- `converged`. `run_inverse` warns when the chain is not converged; the default 2000 steps are too few when `Av` and distance are both free.
- With real photometry the atmosphere library, filter curves and zero points add systematic errors that the likelihood does not include. Use `fit_jitter=True` or `sys_floor`. `run_inverse` warns when `chi2_best / n_data` exceeds 5 with `fit_jitter=False`.
- For the Hermite methods, cells that touch an empty node, including the extra derivative nodes, are rejected during sampling.

### Parameter modes

```python
from sed_model import FitParams, fixed, free, bounded

params = FitParams(
    teff=free(5000.0, 6500.0),
    logg=fixed(4.44),
    meta=free(-1.0, 0.5),
    a_v=fixed(0.0),
    d=fixed(3.086e19),
)
```

`fit_params_from_grid(grid, a_v=(0.0, 2.0), d_cm=3.086e19)` builds the same object from the grid bounds.

---

## Testing

```bash
pytest tests/test_kernels.py     # Fortran kernels and filter handling; no data needed
pytest tests/test_forward.py     # forward model against MESA reference values
pytest tests/test_inverse.py     # injection and recovery
pytest tests/test_reference.py   # frozen MESA and MSG reference points
```

Run them from outside the source tree against an installed build. Inside it, `import sed_model` finds the uncompiled sources and every test skips. Tests that need data look in `~/.sed_tools/data` and `~/SED_Tools/data`; the `CC_*` variables in [Testing and Validation](docs/testing.md) override these. Without the extension or the data, tests skip with a reason.

Cross-code results (0.4.1 wheel):

- **MESA colors:** SEDs and bolometric magnitudes bit-identical; band magnitudes within 3e-13 mag.
- **MSG:** linear interpolation matches to about 1e-7 in the spectra; Hermite photometry to better than 6 mmag where the derivative stencil is complete.
- **MIST:** median offsets of 0.005 – 0.02 mag in V, I, J, Ks and Gaia. These come from the different atmosphere library, not from the code.

What to expect when comparing with other codes, and the inverse-model calibration, are in [Testing and Validation](https://nialljmiller.github.io/SED_Model/testing/).

---

## Relationship to MESA and SED_Tools

SED_Model is a spin-off of the `colors` module of [MESA](https://docs.mesastar.org/), developed by the same author. That module computes bolometric corrections and synthetic photometry inside a stellar-evolution run; SED_Model moves it into a standalone package for use on observed photometry. The Fortran kernels follow the module file by file. The inverse model, the Python interface, the filter and extinction handling and the empty-cell warnings are specific to SED_Model. SED_Model does not need MESA and does not read a MESA installation.

[SED_Tools](https://github.com/nialljmiller/SED_Tools) downloads and standardises atmosphere grids and filters, builds the grid cubes and completes SEDs. SED_Model reads its output.

---

## Repository layout

```text
sed_model/    Python package (grid, filters, forward, inverse, params, io, extinction)
fortran/      kernels, built into sed_model.cc_api (see docs/fortran.md)
demos/        example scripts, embedded in the documentation
docs/         MkDocs sources
tests/
tools/        build helpers used by CI
```

---

## Citing and credits

Cite SED_Model with [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button reads it).

- **MESA:** the Fortran kernels are ported from its `colors` module and keep its licence (LGPL-3.0-or-later). Please also cite MESA: Paxton et al. 2011, ApJS 192, 3; 2013, ApJS 208, 4; 2015, ApJS 220, 15; 2018, ApJS 234, 34; 2019, ApJS 243, 10; Jermyn et al. 2023, ApJS 265, 15.
- **SED_Tools** produces the grids and filter files.
- **emcee:** Foreman-Mackey et al. 2013, PASP 125, 306.
- **Extinction laws**, ported from the MESA `colors` extinction module (coefficients as in `dust_extinction`): Cardelli, Clayton & Mathis (1989), O'Donnell (1994), Fitzpatrick (1999), Fitzpatrick & Massa (2007), Gordon et al. (2003), Calzetti et al. (2000), Gordon et al. (2023). References are in the [Extinction guide](https://nialljmiller.github.io/SED_Model/guide/extinction/).
- **Cross-checks:** [MSG](https://github.com/rhdtownsend/msg) and the MIST bolometric-correction tables (Choi et al. 2016, ApJ 823, 102; Dotter 2016, ApJS 222, 8).

## License

The Python source, documentation, examples, tests and packaging files are MIT. The Fortran kernels in `fortran/` are LGPL-3.0-or-later. See `LICENSE` and the source file headers.

---

<p align="center">
  <img src="https://raw.githubusercontent.com/nialljmiller/SED_Model/main/docs/SED_Model_Logo_full.png" alt="SED_Model: The Pythonic Inside Story of MCMC-Based Stellar Parameter Inference from Photometric Grids, by Niall Miller" width="360"/>
</p>
