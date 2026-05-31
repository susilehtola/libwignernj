#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""Extract a single version's section from CHANGELOG.md.

Used by `.github/workflows/publish-pypi.yml` to feed the release-notes
body to `gh release create` on every v*-tag push, and by the
release-backfill flow that re-creates GitHub releases from existing
tags.  Prints the section body (everything between `## [X.Y.Z]` and
the next `## ` header), stripped of the version header line itself.
Exits 1 if the version is not found.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def extract(text: str, version: str) -> str:
    header_pat = re.compile(
        r"^## \[" + re.escape(version) + r"\][^\n]*\n", re.MULTILINE
    )
    m = header_pat.search(text)
    if not m:
        raise SystemExit(f"version '{version}' not found in CHANGELOG")
    start = m.end()
    nxt = re.search(r"^## ", text[start:], re.MULTILINE)
    end = start + nxt.start() if nxt else len(text)
    return text[start:end].strip("\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("version", help="Version string, e.g. 0.6.2 (no leading v)")
    ap.add_argument("--changelog", type=Path, default=Path("CHANGELOG.md"))
    args = ap.parse_args()
    print(extract(args.changelog.read_text(), args.version))
    return 0


if __name__ == "__main__":
    sys.exit(main())
