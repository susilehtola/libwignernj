#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""Verify that every version string in the repository agrees.

libwignernj's version lives in four places that the build system,
the wheel pipeline, and the runtime extension each read independently:

  1. ``CMakeLists.txt``       — ``project(wignernj VERSION ...)`` →
                                ``libwignernj.so.<MAJOR>.<MINOR>.<PATCH>``
                                + CMake package config files.
  2. ``pyproject.toml``       — PEP 621 ``[project] version = "..."`` →
                                what ``pip install`` records in the
                                wheel METADATA.
  3. ``setup.py``             — fallback for setuptools<77, kept in
                                sync as a defensive measure.
  4. ``wignernj/__init__.py`` — ``__version__`` exposed to Python
                                consumers via the package object.

When these disagree, downstream consumers see inconsistent results
(`pip list` says one version, `import wignernj; wignernj.__version__`
says another, the installed `.so` SONAME says a third).  This script
exits 0 only when all four agree; otherwise it prints a per-source
table and exits 1.  CI runs it in the same fast job as
``check_source_lists.py``; run locally before pushing to catch
drift.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib  # type: ignore[import-not-found]

REPO = Path(__file__).resolve().parent.parent


def read_cmake_version() -> str:
    text = (REPO / "CMakeLists.txt").read_text(encoding="utf-8")
    m = re.search(r"\bproject\(\s*wignernj\s+VERSION\s+([0-9][0-9.]*)", text)
    if not m:
        sys.exit("error: could not find project(wignernj VERSION ...) in CMakeLists.txt")
    return m.group(1)


def read_pyproject_version() -> str:
    with (REPO / "pyproject.toml").open("rb") as f:
        data = tomllib.load(f)
    try:
        return data["project"]["version"]
    except KeyError:
        sys.exit("error: pyproject.toml has no [project] version field")


def read_setup_py_version() -> str:
    text = (REPO / "setup.py").read_text(encoding="utf-8")
    m = re.search(r"""version\s*=\s*["']([0-9][0-9.]*)["']""", text)
    if not m:
        sys.exit("error: could not find version=\"...\" in setup.py")
    return m.group(1)


def read_init_version() -> str:
    text = (REPO / "wignernj" / "__init__.py").read_text(encoding="utf-8")
    m = re.search(r"""__version__\s*=\s*["']([0-9][0-9.]*)["']""", text)
    if not m:
        sys.exit("error: could not find __version__ in wignernj/__init__.py")
    return m.group(1)


def main() -> int:
    versions = {
        "CMakeLists.txt": read_cmake_version(),
        "pyproject.toml": read_pyproject_version(),
        "setup.py": read_setup_py_version(),
        "wignernj/__init__.py": read_init_version(),
    }

    distinct = set(versions.values())
    if len(distinct) == 1:
        (version,) = distinct
        print(f"OK: all 4 version strings agree at {version}.")
        return 0

    print("FAIL: version strings disagree across the build:", file=sys.stderr)
    width = max(len(name) for name in versions)
    for name, ver in versions.items():
        print(f"  {name:<{width}}  {ver}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
