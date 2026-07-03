# Fortran restructure: mirroring MESA's colors module

## Integration steps (do these, in order)

1. Replace your `SED_Model/fortran/` directory entirely with the
   `fortran/` in this bundle (deletes the old `cc_kernels.f90` +
   `cc_api.f90`, adds the `public/` and `private/` split below). If
   you're doing this by hand: **delete `fortran/` first, then copy
   this bundle's `fortran/` in** — don't merge directory-by-directory,
   that's almost certainly how `fortran/cc_api.f90` went missing last
   time (git showed it as deleted, and it wasn't only in
   `fortran/private/` or `fortran/public/` afterward, so the file was
   the exact one your `make`/`pip install` errors were pointing at
   before it got dropped).
2. Replace `Makefile`, `meson.build`, and `setup.py` with the versions
   in this bundle.
3. Copy `.f2py_f2cmap` to the repo root. Required — see "Bug #1" below.
4. In `pyproject.toml`, change the build-system numpy pin:
   ```toml
   [build-system]
   requires = ["setuptools<70", "wheel", "numpy>=2.0", "meson", "ninja"]
   ```
   (was `"numpy<2.0"`). See "Bug #2" below for why.
5. `pip install -e .` — that's it, no `make` step. Confirmed below.

## New layout

```
fortran/
  public/
    colors_def.f90     <- dp kind + shared constants
    colors_lib.f90      <- aggregates the private kernels (source-level
                            parity with MESA only -- NOT compiled into
                            the extension, see its file header)
  private/
    colors_utils.f90    <- find_interval, find_containing_cell, find_nearest_point,
                            dilute_flux, trapezoidal_integration, simpson_integration,
                            interp_filter_onto_sed
    hermite_interp.f90  <- hermite_interp_vector
    linear_interp.f90   <- trilinear_interp_vector
    synthetic.f90        <- calculate_synthetic_flux, magnitude,
                            compute_{vega,ab,st}_zero_point
    bolometric.f90        <- bolometric_flux, bolometric_magnitude,
                            calculate_bolometric_phot
  cc_api.f90             <- f2py shim (unchanged public surface), `use`s
                            the private kernels directly
```

`hermite_interp.f90` and `linear_interp.f90` use the same module name,
subroutine name, argument order, and assumed-shape array interface as
MESA's `colors/private/hermite_interp.f90` / `linear_interp.f90`. The
one addition is a trailing `optional` `ierr` -- MESA call sites that
don't pass it still work; `cc_api.f90` passes it to get the existing
`clamped`-outside-grid diagnostic.

Not ported, per your call to skip it: `hermite_interp_bounded.f90`,
`knn_interp.f90`. Also not ported, no SED_Model equivalent:
`colors_history.f90`, `colors_ctrls_io.f90` (MESA namelist/history I/O
-- SED_Model's own I/O lives in `sed_model/io.py`), and the
handle-based stencil cache in MESA's `construct_sed_hermite` /
`construct_sed_linear` (SED_Model always gets a preloaded cube from
SED_Tools, nothing to cache).

## Correctness, checked two ways

1. **Fortran-level regression test.** A standalone driver calls both
   the old `cc_kernels` module and the new `colors_utils` /
   `hermite_interp` / `linear_interp` / `synthetic` / `bolometric`
   chain on identical synthetic grid data (in-bounds and
   out-of-bounds query points, filters, zero-points, integration).
   Every output matches to `0.000E+00`, `ierr` flags included.
2. **Built and ran the actual extension.** Not just `gfortran -c` --
   `f2py -c --backend meson`, imported from Python, called every
   public `cc_api` function, checked dtypes and values.

## Two real bugs this surfaced

**Bug #1 -- f2py silently mismapped `dp` to `float`.**
With `dp` defined two `use`-hops from `cc_api.f90`, f2py's kind
resolution fell back to mapping `real(kind=dp)` to C `float` instead
of `double` -- a 4-byte/8-byte mismatch that corrupted memory on the
first call (crashed with `malloc(): invalid size`). Fixed with the
`.f2py_f2cmap` file in this bundle. `f2py` picks it up automatically
from the directory it runs in, so no Makefile/meson changes needed
beyond having the file present at the repo root.

**Bug #2 -- `pip install -e .` failed under the project's own pinned
build requirements.** This one I only caught because I tested the
actual `pip install` path (not just a manual `f2py -c` in a scratch
directory), and it's two separate issues stacked on top of each other:

- `public/colors_lib.f90` -- a module that only re-exports procedures
  from other modules, no `contains` block of its own -- crashes f2py's
  module wrapper with `KeyError: 'void'` under the `numpy<2.0` this
  project's `pyproject.toml` pins for the build. f2py tries to build a
  Python getter for every public name in a wrapped module, can't find
  a variable type for a re-exported subroutine, and falls over. Fix:
  `cc_api.f90` now `use`s the private kernels (`hermite_interp`,
  `linear_interp`, `colors_utils`, `synthetic`, `bolometric`,
  `colors_def`) directly instead of routing through `colors_lib`, and
  `colors_lib.f90` is excluded from the Makefile/meson.build/setup.py
  source lists that build the extension. It's still in the repo,
  still compiles standalone, still matches MESA's file -- it's just
  not part of what gets built into `cc_api`.
- Separately: `pyproject.toml`'s `[build-system] requires` pins
  `numpy<2.0` for the *build*, while `[project] dependencies` allows
  `numpy>=1.22` at *runtime* -- so a fresh `pip install` builds the
  extension against NumPy 1.x's ABI and then installs NumPy 2.x to
  satisfy the runtime dependency, and the import fails with `"A module
  that was compiled using NumPy 1.x cannot be run in NumPy 2.x"`. This
  is pre-existing and independent of the Fortran restructure -- it
  would have broken the original `cc_kernels.f90`/`cc_api.f90` under
  the same `pyproject.toml` just as much. Fix: bump the build-time pin
  to `numpy>=2.0` so build-time and run-time NumPy are ABI-consistent.

Verified after both fixes, in a throwaway venv, starting from nothing:
```
pip install -e .     # editable install
pip install .         # plain install, separate fresh venv
```
Both exit 0, no `make` invoked at any point, and
`from sed_model import cc_api` immediately works with correct,
bit-for-bit-matching `float64` results.
