# No languages are declared here: the python-only product needs no compiler,
# and for the others the checks below give a readable error if one is missing.
project(
  'sed_model',
  version: '0.4.5',
  meson_version: '>=1.1.0',
)

# The same source tree produces three products, selected with
# -Dcomponents=...  (see meson_options.txt):
#   all     Python package + compiled f2py extension  -> platform wheel
#   python  Python sources only                       -> pure py3-none-any wheel
#   fortran standalone Fortran libraries + .mod files -> no Python required
components = get_option('components')

# public/colors_lib.f90 is NOT part of the f2py extension -- see its file
# header and fortran/cc_api.f90's header for why (f2py cannot wrap a module
# that only re-exports procedures from other modules).  It is part of the
# standalone Fortran library, which is its intended use.
kernel_sources = files(
  'fortran/public/colors_def.f90',
  'fortran/private/colors_utils.f90',
  'fortran/private/hermite_interp.f90',
  'fortran/private/linear_interp.f90',
  'fortran/private/hermite_interp_bounded.f90',
  'fortran/private/synthetic.f90',
  'fortran/private/bolometric.f90',
)

if components in ['all', 'python']
  py = import('python').find_installation(pure: components == 'python')

  py.install_sources(
    files(
      'sed_model/__init__.py',
      'sed_model/_cc_ext.py',
      'sed_model/_inv_core.py',
      'sed_model/filters.py',
      'sed_model/forward.py',
      'sed_model/grid.py',
      'sed_model/inverse.py',
      'sed_model/io.py',
      'sed_model/params.py',
      'sed_model/sed_extinction.py',
    ),
    subdir: 'sed_model',
  )
endif

# pip falls back to building the sdist on platforms with no prebuilt wheel
# (e.g. Windows, musl Linux, older macOS).  Fail with instructions rather than
# a bare "compiler not found".  The extension needs C as well (f2py's
# fortranobject.c); the standalone Fortran libraries need only Fortran.
if components in ['all', 'fortran']
  needed_languages = components == 'all' ? ['c', 'fortran'] : ['fortran']
  if not add_languages(needed_languages, native: false, required: false)
    error('''
SED_Model needs C and Fortran compilers to build from source, and at least one
was not found.  Prebuilt wheels exist for Linux (x86_64, aarch64; glibc) and
macOS (arm64 on 14+, x86_64 on 15+) -- on those platforms, upgrade pip
(python -m pip install -U pip) and use a supported Python (3.10-3.14).
Otherwise install a toolchain and retry:
  macOS:         xcode-select --install; brew install gfortran
  Debian/Ubuntu: sudo apt install gcc gfortran
  Alpine:        apk add gcc gfortran musl-dev
  conda:         conda install gfortran
  Windows:       not supported natively -- use WSL
''')
  endif
endif

if components == 'all'

  incdir_numpy = run_command(
    py, '-c', 'import numpy; print(numpy.get_include())',
    check: true,
  ).stdout().strip()
  incdir_f2py = run_command(
    py, '-c',
    'import numpy.f2py; import pathlib; print(pathlib.Path(numpy.f2py.__file__).parent / "src")',
    check: true,
  ).stdout().strip()

  extension_sources = [kernel_sources, files('fortran/cc_api.f90')]

  # On macOS, link the Fortran runtime into the extension instead of bundling
  # Homebrew's dylibs: delocate would otherwise tie the wheel's minimum macOS
  # to whatever version the installed gfortran bottle targets.  (On Linux,
  # auditwheel bundles libgfortran, which is the established practice there.)
  # The flags are passed unconditionally: gfortran always accepts them, and
  # meson's usual "is this argument supported" probe links with
  # -Wl,-fatal_warnings, so an harmless "library built for a newer macOS"
  # warning would make it answer NO and silently fall back to dynamic linking.
  extension_link_args = []
  if host_machine.system() == 'darwin'
    extension_link_args = ['-static-libgfortran', '-static-libquadmath', '-static-libgcc']
  endif

  # f2py generates the C wrapper and the Fortran-module glue.  The kind map
  # makes `real(dp)` a C double.  Meson resolves Fortran module order itself.
  cc_api_wrappers = custom_target(
    'cc_api_f2py',
    input: extension_sources,
    output: ['cc_apimodule.c', 'cc_api-f2pywrappers2.f90'],
    command: [
      py, '-m', 'numpy.f2py', '@INPUT@',
      '-m', 'cc_api', '--lower',
      '--f2cmap', meson.current_source_dir() / '.f2py_f2cmap',
    ],
  )

  py.extension_module(
    'cc_api',
    [extension_sources, cc_api_wrappers, incdir_f2py / 'fortranobject.c'],
    include_directories: include_directories(incdir_numpy, incdir_f2py),
    dependencies: py.dependency(),
    subdir: 'sed_model',
    install: true,
    link_args: extension_link_args,
    # Link with the Fortran driver (as f2py does) so it supplies the Fortran
    # runtime and honours the -static-lib* flags above.
    link_language: 'fortran',
    fortran_args: ['-O2', '-fPIC'],
    c_args: ['-O2', '-fPIC', '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'],
  )
endif

if components == 'fortran'
  fc = meson.get_compiler('fortran')

  # Ship the libraries self-contained.  On macOS the Fortran runtime is linked
  # in statically.  On Linux that is not possible for a shared library on
  # every architecture (the distro's libgfortran.a is not position-independent
  # on x86_64), so tools/ci/build_fortran_lib.sh copies libgfortran and
  # libquadmath next to the library and points its rpath at them instead.
  link_args = []
  if host_machine.system() == 'darwin'
    link_args = ['-static-libgfortran', '-static-libquadmath', '-static-libgcc']
  endif

  kernels = both_libraries(
    'sed_model_kernels',
    [kernel_sources, files('fortran/public/colors_lib.f90')],
    fortran_args: ['-O2', '-fPIC'],
    link_args: link_args,
    install: true,
  )

  install_data('LICENSE', 'README.md', install_dir: get_option('datadir') / 'sed_model')

  meson.add_install_script('tools/post_install_fortran.py')

  smoke = executable(
    'fortran_smoke',
    'tests/fortran/smoke.f90',
    link_with: kernels.get_static_lib(),
    link_args: link_args,
  )
  test('fortran_kernels', smoke)
endif
