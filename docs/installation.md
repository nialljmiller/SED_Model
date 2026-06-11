# Installation

## From PyPI

```bash
pip install sed-model
```

The PyPI release is a source distribution: pip builds the Fortran extension on your machine during installation, so you need a Fortran compiler available (see [Requirements](#requirements)).

## From source

```bash
git clone https://github.com/nialljmiller/SED_Model.git
cd SED_Model
python -m pip install -e .
```

### Building the Fortran extension

Performance-critical paths (SED interpolation, flux dilution, integration, filter convolution, zero-points) live in two Fortran files:

```text
fortran/cc_kernels.f90
fortran/cc_api.f90
```

These are compiled into the `sed_model.cc_api` extension module. For an editable install, build in place with either:

```bash
python setup.py build_ext --inplace
```

or:

```bash
make
```

Both routes drive `f2py`; with NumPy ≥ 1.26 the `f2py` backend uses Meson and Ninja under the hood.

!!! note "The Fortran extension is required at run time, not import time"
    `import sed_model` succeeds without the compiled extension — `cc_api` is loaded lazily on the first `run_forward` call. If the extension is missing you will get an `ImportError` at that point with instructions to run `make`.

## Requirements

| Requirement | Purpose |
|---|---|
| Python ≥ 3.10 | |
| `numpy` ≥ 1.22 | arrays, `f2py` |
| `pandas` ≥ 1.5 | lookup-table validation |
| `emcee` ≥ 3.1 | inverse-model MCMC sampling |
| `gfortran` | compiling the Fortran kernels |
| `meson`, `ninja` | `f2py` build backend |
| `matplotlib` | optional — plotting in the demo scripts only |

The optional extras defined in `pyproject.toml`:

```bash
pip install "sed-model[dev]"    # adds pytest
pip install "sed-model[demo]"   # adds matplotlib
```

## Data products

SED_Model consumes data prepared by [SED_Tools](https://github.com/nialljmiller/SED_Tools):

- an **atmosphere grid directory** containing `flux_cube.bin` (and a `lookup_table.csv` for optional validation), and
- **filter transmission files** (two-column `.dat`: wavelength in Å, transmission 0–1), optionally organised into instrument directories with an index file.

A Vega reference SED (CSV: wavelength Å, flux erg/s/cm²/Å) is required only if you want Vega-system zero-points. See [Grids and Filters](guide/grids-and-filters.md) for the formats in detail.
