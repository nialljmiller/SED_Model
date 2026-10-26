"""The README Quick Start and demos/api_usage.py actually run."""

from __future__ import annotations

import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from _data_paths import KURUCZ_DIR, JOHNSON_DIR, VEGA_SED

pytest.importorskip("sed_model.cc_api", reason="cc_api Fortran extension not built")

REPO = Path(__file__).resolve().parents[1]

for _p, _what in ((KURUCZ_DIR, "Kurucz2003all grid"), (JOHNSON_DIR, "Johnson filters"),
                  (VEGA_SED, "Vega SED")):
    if not os.path.exists(_p):
        pytest.skip(f"{_what} not found at {_p}", allow_module_level=True)

DATA_ROOT = Path(KURUCZ_DIR).resolve().parents[1]      # .../data


def _find_filter(name: str) -> str:
    hits = glob.glob(str(DATA_ROOT / "filters" / "*" / "*" / name))
    if not hits:
        pytest.skip(f"filter {name} not found under {DATA_ROOT / 'filters'}")
    return hits[0]


def _quick_start() -> str:
    text = (REPO / "README.md").read_text()
    section = text.split("## Quick Start", 1)[1]
    return re.search(r"```python\n(.*?)```", section, re.S).group(1)


def test_readme_quick_start(monkeypatch, tmp_path):
    code = _quick_start()
    code = code.replace("/path/to/Kurucz2003all/", str(KURUCZ_DIR))
    code = code.replace("/path/to/filters/GAIA/G.dat", _find_filter("G.dat"))
    code = code.replace("/path/to/filters/2MASS/J.dat", _find_filter("J.dat"))
    code = code.replace("/path/to/stellar_models/vega_flam.csv", str(VEGA_SED))
    assert "/path/to" not in code
    monkeypatch.chdir(tmp_path)
    exec(compile(code, "README.md:QuickStart", "exec"), {"__name__": "__quickstart__"})


def test_demo_api_usage(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "stellar_models").symlink_to(DATA_ROOT / "stellar_models")
    (tmp_path / "data" / "filters").symlink_to(DATA_ROOT / "filters")
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.run([sys.executable, str(REPO / "demos" / "api_usage.py")],
                          cwd=tmp_path, env=env, capture_output=True, text=True,
                          timeout=900)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-3000:]
    assert (tmp_path / "forward_magnitudes.csv").exists()
    assert (tmp_path / "posterior_samples.npz").exists()
