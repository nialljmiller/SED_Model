# Installation

SED_Model needs a compiled Fortran extension. You can get it in whichever form suits you; every release publishes all four combinations of Fortran/Python, built and unbuilt:

| You want | Get | Fortran | Python | Needs gfortran? |
|---|---|---|---|---|
| Just use it | **Platform wheel** — `pip install sed-model` | built, bundled | packaged | no |
| Python only (bring your own extension) | **Pure-Python wheel** `sed_model-*-py3-none-any.whl` (GitHub release) | not built | packaged | no |
| Fortran without Python | **Fortran libraries** `sed_model-fortran-*-<platform>.tar.gz` (GitHub release) | built | not needed | no |
| Build everything yourself | **Source distribution** or `git clone` | built locally | from source | yes |

## Compatibility

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

Each release also provides, on the [GitHub releases page](https://github.com/nialljmiller/SED_Model/releases), a pure-Python wheel (no compiled extension) and standalone Fortran libraries (no Python) for the Linux and macOS platforms above.

## From PyPI (platform wheels)

```bash
pip install sed-model
```

Prebuilt wheels are published for Linux (x86_64, aarch64) and macOS (arm64 on macOS 14+, x86_64 on macOS 15+), Python 3.10–3.14. The Fortran runtime is bundled inside the wheel, so nothing else needs installing. On any other platform (including Windows), pip falls back to the source distribution and builds locally — see [From source](#from-source).

## Pure-Python wheel

`sed_model-<version>-py3-none-any.whl` on the [GitHub releases page](https://github.com/nialljmiller/SED_Model/releases) contains the Python package without the compiled extension. It is deliberately **not** on PyPI, where it would be picked on platforms that have no platform wheel and then fail at run time. Use it when you supply the extension yourself — for instance a locally built `cc_api*.so` dropped into the installed `sed_model/` directory.

## Standalone Fortran libraries

`sed_model-fortran-<version>-<platform>.tar.gz` (GitHub releases) holds the Fortran kernels with no Python dependency:

```text
lib/libsed_model_kernels.a       static library
lib/libsed_model_kernels.so      shared library (.dylib on macOS)
lib/libgfortran.so.5, ...        bundled Fortran runtime (Linux only)
include/sed_model/*.mod          Fortran module files (use colors_lib)
```

No system `libgfortran` is needed to use the shared library: on macOS the runtime is linked in statically, and on Linux it is shipped next to the library and found through an `$ORIGIN` rpath, so keep the files in `lib/` together. From Fortran:

```bash
gfortran -I<prefix>/include/sed_model my_program.f90 <prefix>/lib/libsed_model_kernels.a
```

with `use colors_lib` in your program. The `.mod` files are gfortran module-format files (version 15 for the current builds) and can only be read by compilers that understand that format; they have been verified with gfortran 12. If your compiler rejects them, build from source instead (`-Dcomponents=fortran`). The module is the same aggregation point MESA's `colors` module provides (`hermite_interp_vector`, `hermite_bounded_interp_vector`, `trilinear_interp_vector`, `calculate_synthetic_flux`, `calculate_bolometric_phot`, the zero-point functions, …). Linux builds target glibc 2.28 or newer.

## From source

```bash
git clone https://github.com/nialljmiller/SED_Model.git
cd SED_Model
python -m pip install -e .
```

The build backend is [meson-python](https://mesonbuild.com/meson-python/): it runs `f2py` over the Fortran sources, compiles the extension, and installs it with the Python package. You need `gfortran`; `meson`, `ninja` and NumPy are fetched automatically into the isolated build environment. An editable install rebuilds the extension when the Fortran changes.

The same source tree builds each product directly with Meson:

| Product | Command |
|---|---|
| Platform wheel | `python -m build --wheel` |
| Pure-Python wheel | `python -m build --wheel -Csetup-args=-Dcomponents=python` |
| Fortran libraries | `meson setup build -Dcomponents=fortran --prefix=$PWD/stage && meson compile -C build && meson test -C build && meson install -C build` |
| Source distribution | `python -m build --sdist` |

### The Fortran sources

Performance-critical paths (SED interpolation, flux dilution, integration, filter convolution, zero-points) live in:

```text
fortran/public/colors_def.f90
fortran/private/{colors_utils,hermite_interp,linear_interp,hermite_interp_bounded,synthetic,bolometric}.f90
fortran/cc_api.f90
```

These are compiled into the `sed_model.cc_api` extension module. `make` is kept as a quick alternative that drives `f2py` directly and writes `sed_model/cc_api*.so` in place.

!!! note "The Fortran extension is required at run time, not import time"
    `import sed_model` succeeds without the compiled extension — `cc_api` is loaded lazily on the first `run_forward` call. If the extension is missing you will get an `ImportError` at that point explaining the options above.

## Continuous integration and releases

`.github/workflows/build.yml` builds and tests all four products on every push and pull request (uploaded as workflow artifacts). Releasing needs no manual step:

1. Bump `version` in `meson.build` and merge to `main`.
2. Once every platform build and test is green, CI publishes the platform wheels and the source distribution to PyPI, but only if that version is not already there. Pushes that do not change the version publish nothing.
3. CI then installs the new release from PyPI, binary-only, on Linux and macOS (Python 3.10 and 3.14) and runs the kernel tests against it.
4. Finally it creates the `v<version>` tag and GitHub release with every deliverable attached, including the pure-Python wheel and the Fortran libraries.

PyPI credentials: the workflow uses the `PYPI_API_TOKEN` repository secret if it exists, and otherwise [trusted publishing](https://docs.pypi.org/trusted-publishers/) (OIDC, no stored secret). To use trusted publishing, add a publisher for this repository on PyPI once, with workflow `build.yml` and environment `pypi`.

## Requirements

| Requirement | Purpose |
|---|---|
| Python ≥ 3.10 | |
| `numpy` ≥ 2.0 | arrays, `f2py` |
| `pandas` ≥ 1.5 | lookup-table validation |
| `emcee` ≥ 3.1 | inverse-model MCMC sampling |
| `gfortran` | building from source only — not needed for wheels |
| `meson-python`, `ninja` | build backend, fetched automatically |
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
