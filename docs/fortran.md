# Fortran Kernels

The numerically heavy parts of SED_Model are Fortran, laid out to mirror the `colors` module of [MESA](https://docs.mesastar.org/) so that results can be compared with MESA directly and changes in MESA can be followed. The kernels are LGPL-3.0-or-later (see the licence section of the README).

## Layout

```text
fortran/
  public/
    colors_def.f90              kind parameter and shared constants
    colors_lib.f90              aggregates the private kernels under one module
                                name, as MESA's colors_lib does
  private/
    colors_utils.f90            interval search, nearest-model search, flux dilution,
                                trapezoid and Simpson integration, filter resampling
    hermite_interp.f90          cubic Hermite interpolation over (Teff, logg, [M/H])
    hermite_interp_bounded.f90  Hermite with undershoot clipping and a linear fallback
    linear_interp.f90           trilinear interpolation
    synthetic.f90               band fluxes, magnitudes, Vega/AB/ST zero points
    bolometric.f90              bolometric flux and magnitude
  cc_api.f90                    thin wrapper exposed to Python through f2py
```

`cc_api.f90` is compiled into the `sed_model.cc_api` extension. It uses the private kernels directly, not `colors_lib`, because f2py cannot wrap a module that only re-exports procedures from other modules. `colors_lib.f90` is built into the standalone Fortran libraries (`use colors_lib` from your own program) but not into the Python extension.

## Correspondence with MESA

The interpolation and integration kernels use the same module and subroutine names, argument order and assumed-shape array interfaces as the files of the same name in MESA's `colors/private/`. The one addition is a trailing optional `ierr`, which MESA call sites can omit. `cc_api` also takes an optional `photon_counting` argument on the band kernels, used for filters tagged as energy-counting.

Parts of MESA's module with no SED_Model counterpart:

- the handle lifecycle (`colors_init`, `alloc_colors_handle`, ...) and the stencil cache in `construct_sed_hermite` / `construct_sed_linear`, because SED_Model always works on a preloaded cube from SED_Tools;
- `colors_history.f90` and `colors_ctrls_io.f90` (MESA namelist and history I/O; SED_Model's own I/O is in `sed_model/io.py`);
- `knn_interp.f90`; its `interpolate_array` is reproduced in `colors_utils`.

How closely the two agree, and where, is described in [Testing and Validation](testing.md).

## Building

See [Installation](installation.md#from-source) for building the extension, the standalone libraries and the wheels. `make` in the repository root is a quick alternative that drives `f2py` directly and writes `sed_model/cc_api*.so` in place.

## A gotcha when editing the interfaces

f2py resolves `real(kind=dp)` through a kind map. `dp` is defined two `use` levels away from `cc_api.f90`, and without help f2py maps it to C `float`, a 4-byte/8-byte mismatch that corrupts memory on the first call. The repository root therefore carries `.f2py_f2cmap`, which f2py reads from the directory it runs in. Keep it when moving or rebuilding things.
