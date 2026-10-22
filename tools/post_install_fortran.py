#!/usr/bin/env python3
"""Meson post-install step for the standalone Fortran libraries.

1. Copy the compiled .mod files next to the installed libraries so downstream
   Fortran code can `use colors_lib`.
2. On macOS, reset the dylib's install name to @rpath/... (meson rewrites it
   to the absolute build-time prefix, which breaks relocating the tarball).
"""

import os
import pathlib
import shutil
import subprocess
import sys

build_root = pathlib.Path(os.environ["MESON_BUILD_ROOT"])
dest = pathlib.Path(os.environ["MESON_INSTALL_DESTDIR_PREFIX"]) / "include" / "sed_model"
dest.mkdir(parents=True, exist_ok=True)

seen = set()
for mod in sorted(build_root.rglob("*.mod")):
    if mod.name in seen:
        continue
    seen.add(mod.name)
    shutil.copy2(mod, dest / mod.name)
print(f"installed {len(seen)} Fortran module files to {dest}")

if sys.platform == "darwin":
    for dylib in (dest.parent.parent / "lib").glob("libsed_model_kernels*.dylib"):
        if not dylib.is_symlink():
            subprocess.run(
                ["install_name_tool", "-id", f"@rpath/{dylib.name}", str(dylib)],
                check=True,
            )
            print(f"set install name of {dylib.name} to @rpath/{dylib.name}")
