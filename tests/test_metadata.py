"""Release metadata stays consistent: one version everywhere, valid Zenodo file.

These read repository files only (no extension or data needed).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if not (ROOT / "meson.build").exists():            # tests run from an unpacked sdist elsewhere
    pytest.skip("repository files not available", allow_module_level=True)


def _meson_version() -> str:
    return re.search(r"version:\s*'([^']+)'", (ROOT / "meson.build").read_text()).group(1)


def test_version_is_the_same_everywhere():
    v = _meson_version()
    cff = (ROOT / "CITATION.cff").read_text()
    assert re.search(r'^version:\s*"%s"' % re.escape(v), cff, re.M), "CITATION.cff version"
    assert f"Version `{v}`" in (ROOT / "README.md").read_text(), "README.md version"
    assert f"**{v}**" in (ROOT / "docs" / "index.md").read_text(), "docs/index.md version"
    first = re.search(r"^## \[([^\]]+)\]", (ROOT / "CHANGELOG.md").read_text(), re.M).group(1)
    assert first == v, f"top CHANGELOG entry is [{first}], expected [{v}]"


def test_citation_date_matches_changelog():
    changelog = (ROOT / "CHANGELOG.md").read_text()
    v = _meson_version()
    day = re.search(r"^## \[%s\]\s*[—-]\s*(\d{4}-\d{2}-\d{2})" % re.escape(v), changelog, re.M).group(1)
    cff = (ROOT / "CITATION.cff").read_text()
    assert f"date-released: {day}" in cff


def test_zenodo_file_is_valid():
    z = json.loads((ROOT / ".zenodo.json").read_text())
    assert z["upload_type"] == "software"
    assert z["creators"] and all("name" in c for c in z["creators"])
    assert z["license"]
