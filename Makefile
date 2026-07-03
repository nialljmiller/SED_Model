FC      ?= gfortran
F2PY    ?= f2py
FFLAGS  := -O2 -fPIC -Wall

OUTDIR  := sed_model

# Compile order matters: each file only appears after the modules it
# `use`s. This mirrors MESA's colors/public + colors/private split.
#
# public/colors_lib.f90 is deliberately NOT in this list: f2py (under
# the numpy<2.0 pinned below) cannot wrap a module that only
# re-exports procedures from other modules. cc_api.f90 uses the
# private kernels directly instead -- see fortran/cc_api.f90 and
# fortran/public/colors_lib.f90 headers, and MIGRATION.md.
SOURCES := \
	fortran/public/colors_def.f90 \
	fortran/private/colors_utils.f90 \
	fortran/private/hermite_interp.f90 \
	fortran/private/linear_interp.f90 \
	fortran/private/synthetic.f90 \
	fortran/private/bolometric.f90 \
	fortran/cc_api.f90

.PHONY: all clean

all: $(OUTDIR)/cc_api$(shell python3-config --extension-suffix 2>/dev/null || echo .so)

$(OUTDIR)/cc_api$(shell python3-config --extension-suffix 2>/dev/null || echo .so): $(SOURCES)
	$(F2PY) -c \
		--fcompiler=$(FC) \
		--f90flags="$(FFLAGS)" \
		-m cc_api \
		$(SOURCES) \
		--build-dir /tmp/cc_api_build
	mv cc_api*.so $(OUTDIR)/ 2>/dev/null || true
	mv cc_api*.pyd $(OUTDIR)/ 2>/dev/null || true

clean:
	rm -f $(OUTDIR)/cc_api*.so $(OUTDIR)/cc_api*.pyd
	rm -rf /tmp/cc_api_build
