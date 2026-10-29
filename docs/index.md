<p align="center">
  <img src="SED_Model_logo.jpeg" alt="SED_Model logo" width="520"/>
</p>

# SED_Model

**Synthetic photometry and stellar parameter inference.**

SED_Model is a Python package for computing observer-ready synthetic spectral energy distributions (SEDs) and broadband magnitudes from stellar atmosphere grids, and for recovering stellar parameters from observed photometry via Bayesian MCMC inference. It works directly downstream of [SED_Tools](https://github.com/nialljmiller/SED_Tools) and as a standalone drop-in for workflows that use MESA-style atmosphere grids and filter transmission curves. It is a spin-off of the `colors` module of [MESA](https://docs.mesastar.org/), developed by the same author, and reproduces its synthetic photometry; it does not require MESA.

Current version: **0.4.7** ([PyPI](https://pypi.org/project/sed-model/))

<video controls loop muted playsinline width="100%" poster="assets/gallery/sed_fit_still.png">
  <source src="assets/gallery/sed_fit.mp4" type="video/mp4">
</video>

*A run_inverse fit of UBVRI + JHKs photometry. Each curve is one walker's model SED, coloured by its own Teff; the right panel is the same walkers in the Teff–log g plane.*

## Two paired models, one parameter language

- **Forward model** — maps `(Teff, logg, [M/H], R, d)` to an interpolated SED, bolometric flux and magnitude, and synthetic magnitudes in any set of loaded filters.
- **Inverse model** — maps observed magnitudes with uncertainties to a posterior distribution over `(Teff, logg, [M/H])`, and optionally over `Av` and distance, via [`emcee`](https://emcee.readthedocs.io/) MCMC sampling.

Both directions operate through a shared [`FitParams`](guide/parameters.md) object. The same parameter specification flows from the forward model to the likelihood evaluator inside the sampler and back out again.

```text
SED_Tools  ->  build/download SED grids and filters
SED_Model  ->  generate synthetic magnitudes and infer stellar parameters
```

## Features

- Bounded Hermite (default), Hermite and linear interpolation over 3D atmosphere grids `(Teff, logg, [M/H])`.
- Photon-counting filter convolution in Vega, AB, and ST photometric systems.
- Precomputed zero-points matching MESA `colors/private/synthetic.f90` conventions.
- Pure-NumPy interstellar dust extinction with eight supported laws, including Gordon et al. (2003, 2023).
- Bayesian stellar-parameter inference via `emcee` with fixed, bounded, and fully free parameter modes.
- Structured result containers with summaries, persistence to CSV and NPZ, and interpolation diagnostics.
- Fortran kernels for performance-critical paths, compiled via `f2py`/Meson and shipped as prebuilt wheels.
- Reads the `DetectorType` and `WavelengthUnit` tags that SED_Tools writes on each filter file.
- Warns when a query falls in a part of the grid that has no model.
- Forward model checked against the MESA `colors` module (magnitudes agree to double-precision rounding), MSG and MIST; inverse model benchmarked by injection and recovery. See [Testing and Validation](testing.md).

## Where to go next

- [Installation](installation.md) — install from PyPI or source, and build the Fortran extension.
- [Quick Start](quickstart.md) — a forward and an inverse run in a few lines.
- [User Guide](guide/grids-and-filters.md) — grids, filters, both models, parameter modes, extinction, and I/O.
- [Demonstrations](demos/index.md) — full walkthroughs of every script in `demos/`, with the actual source embedded.
- [API Reference](api.md) — generated directly from the docstrings in the source.
- [Fortran kernels](fortran.md) — how the Fortran sources map onto MESA's `colors` module.
- [Testing and Validation](testing.md) — what is tested, and how the results compare with MESA, MSG and MIST.

## Credits

The Fortran kernels are ported from the `colors` module of [MESA](https://docs.mesastar.org/). If you use SED_Model, please also cite MESA: Paxton et al. 2011, ApJS 192, 3; 2013, ApJS 208, 4; 2015, ApJS 220, 15; 2018, ApJS 234, 34; 2019, ApJS 243, 10; Jermyn et al. 2023, ApJS 265, 15. The README lists the other works SED_Model builds on.

## License

SED_Model uses a mixed-license structure. The Python source code, documentation, examples, tests, and packaging files are distributed under the MIT License. The Fortran numerical kernels in `fortran/` are distributed under LGPL-3.0-or-later. See [`LICENSE`](https://github.com/nialljmiller/SED_Model/blob/main/LICENSE) and source file headers for details.
