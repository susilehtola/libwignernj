#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""
Verify that the C source lists in CMakeLists.txt and setup.py stay in sync.

The CMake build (LIB_SOURCES in CMakeLists.txt) and the Python extension
build (sources in setup.py) compile the same library twice, with the
Python build adding one extra file: src/python/wignermodule.c, the CPython
entry point.  When a new file is added to LIB_SOURCES but not to setup.py
the Python extension links with an undefined-symbol error at import time,
which manifests in CI as "0 tests collected" rather than a clear link
failure -- a footgun this script exists to prevent.

Run from the repository root:
    python3 tools/check_source_lists.py

Exits 0 on agreement, 1 with a diff on divergence.
"""

import ast
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PYTHON_ENTRY_POINT = "src/python/wignermodule.c"


def parse_cmake_lib_sources(text: str) -> list[str]:
    """Extract the LIB_SOURCES list from CMakeLists.txt."""
    m = re.search(r"set\(LIB_SOURCES\b(.*?)\)", text, re.DOTALL)
    if not m:
        sys.exit("error: could not find set(LIB_SOURCES ...) in CMakeLists.txt")
    body = m.group(1)
    # One token per line, ignore comments and blank lines.
    return [
        tok
        for line in body.splitlines()
        for tok in [line.split("#", 1)[0].strip()]
        if tok
    ]


def parse_setup_py_sources(text: str) -> list[str]:
    """Extract the `sources = [...]` list from setup.py via the AST."""
    tree = ast.parse(text)
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "sources"
            and isinstance(node.value, ast.List)
        ):
            return [ast.literal_eval(elt) for elt in node.value.elts]
    sys.exit("error: could not find `sources = [...]` in setup.py")


def main() -> int:
    cmake = parse_cmake_lib_sources((REPO / "CMakeLists.txt").read_text())
    setup = parse_setup_py_sources((REPO / "setup.py").read_text())

    expected = [PYTHON_ENTRY_POINT, *cmake]
    if setup == expected:
        print(f"OK: {len(cmake)} library source(s) + 1 Python entry point in sync.")
        return 0

    print("ERROR: setup.py `sources` does not match CMake LIB_SOURCES + Python entry point")
    print()
    print("CMake LIB_SOURCES:")
    for s in cmake:
        print(f"    {s}")
    print()
    print("setup.py sources:")
    for s in setup:
        print(f"    {s}")
    print()
    only_in_cmake = [s for s in cmake if s not in setup]
    only_in_setup = [s for s in setup if s not in cmake and s != PYTHON_ENTRY_POINT]
    if only_in_cmake:
        print("Missing from setup.py:")
        for s in only_in_cmake:
            print(f"    {s}")
    if only_in_setup:
        print("Extra in setup.py:")
        for s in only_in_setup:
            print(f"    {s}")
    if not only_in_cmake and not only_in_setup:
        print("Same files but different order — setup.py must list "
              f"{PYTHON_ENTRY_POINT} first, then LIB_SOURCES in CMake order.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
