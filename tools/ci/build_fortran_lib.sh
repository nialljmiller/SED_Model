#!/usr/bin/env bash
# Build, test and install the standalone Fortran kernels (no Python product).
# Usage: tools/ci/build_fortran_lib.sh <install-prefix>
# Runs unchanged on a dev machine, a macOS runner, or inside a manylinux
# container (where it installs gfortran itself).
set -euo pipefail

OUT=${1:?usage: build_fortran_lib.sh <install-prefix>}
OUT=$(mkdir -p "$OUT" && cd "$OUT" && pwd)

if [ -d /opt/python/cp312-cp312 ]; then          # manylinux container
  # EPEL is not needed for gfortran, and its mirrors are occasionally flaky.
  yum install -y --disablerepo=epel gcc-gfortran >/dev/null \
    || { sleep 30; yum install -y --disablerepo=epel gcc-gfortran >/dev/null; }
  export PATH=/opt/python/cp312-cp312/bin:$PATH
fi

# Install the build tools into a throwaway venv: system Pythons (Homebrew,
# Debian) are "externally managed" and refuse pip installs (PEP 668).
rm -rf .venv-fortran
python3 -m venv .venv-fortran
# shellcheck disable=SC1091
. .venv-fortran/bin/activate
python -m pip install --quiet meson ninja

rm -rf build-fortran
meson setup build-fortran -Dcomponents=fortran --buildtype=release --prefix="$OUT" --libdir=lib
meson compile -C build-fortran
meson test -C build-fortran --print-errorlogs
meson install -C build-fortran

# Guard against silently checking nothing (e.g. a different libdir).
libs=()
for f in "$OUT"/lib/libsed_model_kernels.so "$OUT"/lib/libsed_model_kernels.dylib; do
  [ -e "$f" ] && libs+=("$f")
done
if [ ${#libs[@]} -eq 0 ] || [ ! -e "$OUT/lib/libsed_model_kernels.a" ]; then
  echo "expected libraries not found under $OUT/lib:" >&2
  ls -R "$OUT" >&2
  exit 1
fi

# Make the shared library self-contained.
#  macOS: the runtime is already linked in statically.
#  Linux: vendor libgfortran/libquadmath next to the library with an $ORIGIN
#         rpath (what auditwheel does for wheels).  glibc, libm and libgcc_s
#         are provided by every manylinux-compatible system.
if [ "$(uname)" = "Linux" ]; then
  python -m pip install --quiet patchelf
  so="$OUT/lib/libsed_model_kernels.so"
  for dep in $(ldd "$so" | awk '/libgfortran|libquadmath/ {print $3}'); do
    cp -L "$dep" "$OUT/lib/$(basename "$dep")"
  done
  # libgfortran itself needs libquadmath: bundle that too if it was not direct.
  for dep in $(ldd "$OUT"/lib/libgfortran.so* | awk '/libquadmath/ {print $3}' | sort -u); do
    [ -e "$OUT/lib/$(basename "$dep")" ] || cp -L "$dep" "$OUT/lib/$(basename "$dep")"
  done
  for f in "$OUT"/lib/*.so*; do
    [ -L "$f" ] || patchelf --set-rpath '$ORIGIN' "$f"
  done
fi

# Verify: nothing Fortran-runtime related may resolve from outside the install.
case "$(uname)" in
  Darwin)
    bad=$(otool -L "${libs[@]}" | grep -v ':$' | awk '{print $1}' \
          | grep -v -E '^(/usr/lib/|/System/|@rpath/libsed_model)' || true)
    ;;
  *)
    bad=$(ldd "$OUT/lib/libsed_model_kernels.so" | awk '/libgfortran|libquadmath/ {print $3}' \
          | grep -v -F "$OUT/lib/" || true)
    unresolved=$(ldd "$OUT/lib/libsed_model_kernels.so" | grep -i "not found" || true)
    bad="$bad$unresolved"
    ;;
esac
if [ -n "$bad" ]; then
  echo "Fortran library is not self-contained; unexpected runtime dependencies:" >&2
  echo "$bad" >&2
  exit 1
fi
echo "Fortran library OK: $(ls "$OUT"/lib | tr '\n' ' ')"
