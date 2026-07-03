project(
  'sed_model',
  'c', 'fortran',
  version: '0.2.0',
  meson_version: '>=1.1.0',
)

py = import('python').find_installation(pure: false)

np_inc = run_command(
  py, '-c', 'import numpy; print(numpy.get_include())',
  check: true,
).stdout().strip()

f2py_inc = run_command(
  py, '-c', 'import numpy.f2py; import pathlib; print(pathlib.Path(numpy.f2py.__file__).parent / "src")',
  check: true,
).stdout().strip()

# public/colors_lib.f90 is deliberately NOT in this list -- see its
# file header and fortran/cc_api.f90's header for why (f2py can't wrap
# a re-export-only module under this project's pinned numpy<2.0).
py.extension_module(
  'cc_api',
  [
    'fortran/public/colors_def.f90',
    'fortran/private/colors_utils.f90',
    'fortran/private/hermite_interp.f90',
    'fortran/private/linear_interp.f90',
    'fortran/private/synthetic.f90',
    'fortran/private/bolometric.f90',
    'fortran/cc_api.f90',
    f2py_inc / 'fortranobject.c',
  ],
  include_directories: [np_inc, f2py_inc],
  dependencies: py.dependency(),
  subdir: 'sed_model',
  install: true,
  fortran_args: ['-O2', '-fPIC'],
  c_args: ['-O2', '-fPIC'],
)
