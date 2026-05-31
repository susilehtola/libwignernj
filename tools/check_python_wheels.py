#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""Detect currently-supported CPython releases missing from the wheel matrix.

Reads the endoflife.date Python feed (`releaseDate <= today < eol`) and
checks each active cycle against
  - `CIBW_BUILD` in `.github/workflows/publish-pypi.yml`
  - `Programming Language :: Python :: 3.X` classifiers in
    `pyproject.toml`

Emits `missing_wheel=...` and `missing_classifier=...` lines to the
file pointed to by `--github-output`, suitable for $GITHUB_OUTPUT in
GitHub Actions, and prints a one-line summary to stdout.  Exits 0
unless an argument is malformed; the calling workflow decides what to
do with a non-empty gap.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import date
from pathlib import Path


def active_cycles(eol_json: Path, today: date) -> list[str]:
    cycles = json.loads(eol_json.read_text())
    out: list[str] = []
    for c in cycles:
        cyc = c.get("cycle", "")
        if not cyc.startswith("3."):
            continue
        try:
            rel = date.fromisoformat(c["releaseDate"])
            eol = date.fromisoformat(c["eol"])
        except (KeyError, TypeError, ValueError):
            continue
        if rel <= today < eol:
            out.append(cyc)
    # Sort numerically (3.9 < 3.10 < 3.14, not lexicographic).
    out.sort(key=lambda v: tuple(int(p) for p in v.split(".")))
    return out


def cibw_tags(workflow: Path) -> set[str]:
    """Extract every `cp3XX-*` tag mentioned on a CIBW_BUILD line."""
    tags: set[str] = set()
    pat = re.compile(r"cp(3\d+)-")
    in_build = False
    for line in workflow.read_text().splitlines():
        if "CIBW_BUILD" in line:
            in_build = True
        if in_build:
            tags.update(pat.findall(line))
            # CIBW_BUILD may be a single line or YAML-folded across
            # several; stop accumulating once we see a key that is
            # clearly a sibling (a new uppercase env var or a
            # dedented YAML key).
            stripped = line.strip()
            if stripped and "CIBW_BUILD" not in line and ":" in stripped \
                    and not stripped.startswith(("cp", "\"cp", "'cp", "-")):
                in_build = False
    return tags


def classifier_versions(pyproject: Path) -> set[str]:
    pat = re.compile(r'Programming Language :: Python :: (3\.\d+)\b')
    return set(pat.findall(pyproject.read_text()))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eol-json", required=True, type=Path)
    ap.add_argument("--workflow", required=True, type=Path)
    ap.add_argument("--pyproject", required=True, type=Path)
    ap.add_argument("--github-output", type=Path,
                    help="If given, append `missing_wheel=` and "
                         "`missing_classifier=` lines for $GITHUB_OUTPUT.")
    ap.add_argument("--today", type=date.fromisoformat,
                    default=date.today(),
                    help="Override 'today' for testing (YYYY-MM-DD).")
    args = ap.parse_args()

    active = active_cycles(args.eol_json, args.today)
    tags = cibw_tags(args.workflow)
    classified = classifier_versions(args.pyproject)

    missing_wheel = [v for v in active if v.replace(".", "") not in tags]
    missing_class = [v for v in active if v not in classified]

    print(f"Active CPython cycles: {' '.join(active) or '(none)'}")
    print(f"Missing wheel tags:    {' '.join(missing_wheel) or '(none)'}")
    print(f"Missing classifiers:   {' '.join(missing_class) or '(none)'}")

    if args.github_output:
        with args.github_output.open("a") as f:
            f.write(f"missing_wheel={' '.join(missing_wheel)}\n")
            f.write(f"missing_classifier={' '.join(missing_class)}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
