#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""Smoke-test that the CMake-installed wignernj .dist-info is well-formed.

Run from CI right after `cmake --install` against a staging prefix:

    python3 tools/verify_python_metadata.py STAGED_SITEARCH EXPECTED_VERSION

Verifies that

  1. importlib.metadata.distribution("wignernj") resolves under
     STAGED_SITEARCH and reports EXPECTED_VERSION;
  2. every file listed in RECORD exists on disk;
  3. the recorded sha256 hash matches the file's actual contents
     (with the convention that the RECORD row for RECORD itself
     carries an empty hash);
  4. METADATA carries the required PEP 566 fields
     (Name, Version, Summary, License);
  5. WHEEL, INSTALLER, top_level.txt are present and non-empty.

Exits 0 on success, 1 with a diagnostic on any mismatch.
"""
from __future__ import annotations

import base64
import hashlib
import sys
from pathlib import Path


def _record_hash(data: bytes) -> str:
    digest = hashlib.sha256(data).digest()
    return "sha256=" + base64.urlsafe_b64encode(digest).rstrip(b"=").decode("ascii")


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(
            f"usage: {argv[0]} STAGED_SITEARCH EXPECTED_VERSION", file=sys.stderr
        )
        return 2
    sitearch = Path(argv[1]).resolve()
    expected_version = argv[2]

    if not sitearch.is_dir():
        print(f"FAIL: staging sitearch {sitearch} does not exist", file=sys.stderr)
        return 1

    # Use importlib.metadata via PathDistribution so the test is scoped
    # to the staging tree alone, without depending on the ambient sys.path.
    sys.path.insert(0, str(sitearch))
    import importlib.metadata as md
    importlib_md_path_cls = md.PathDistribution

    candidates = sorted(sitearch.glob("wignernj-*.dist-info"))
    if not candidates:
        print(
            f"FAIL: no wignernj-*.dist-info directory under {sitearch}",
            file=sys.stderr,
        )
        return 1
    if len(candidates) > 1:
        print(
            f"FAIL: multiple dist-info directories under {sitearch}: {candidates}",
            file=sys.stderr,
        )
        return 1
    dist_info = candidates[0]

    dist = importlib_md_path_cls(dist_info)
    name = dist.metadata["Name"]
    version = dist.version
    if name != "wignernj":
        print(f"FAIL: Name = {name!r}, expected 'wignernj'", file=sys.stderr)
        return 1
    if version != expected_version:
        print(
            f"FAIL: Version = {version!r}, expected {expected_version!r}",
            file=sys.stderr,
        )
        return 1

    required_metadata = ("Name", "Version", "Summary", "License")
    for field in required_metadata:
        if not dist.metadata[field]:
            print(f"FAIL: METADATA missing required field {field}", file=sys.stderr)
            return 1

    for fname in ("WHEEL", "INSTALLER", "top_level.txt", "RECORD"):
        f = dist_info / fname
        if not f.is_file():
            print(f"FAIL: missing {fname}", file=sys.stderr)
            return 1
        if f.stat().st_size == 0:
            print(f"FAIL: {fname} is empty", file=sys.stderr)
            return 1

    record_path = dist_info / "RECORD"
    record_lines = record_path.read_text(encoding="utf-8").splitlines()
    record_row_for_record = (
        f"{dist_info.name}/RECORD"
    )
    saw_record_row = False
    for line in record_lines:
        if not line.strip():
            continue
        parts = line.split(",")
        if len(parts) != 3:
            print(f"FAIL: malformed RECORD row {line!r}", file=sys.stderr)
            return 1
        rel, recorded_hash, recorded_size = parts
        if rel == record_row_for_record:
            if recorded_hash != "" or recorded_size != "":
                print(
                    f"FAIL: RECORD row for RECORD has non-empty hash/size: {line!r}",
                    file=sys.stderr,
                )
                return 1
            saw_record_row = True
            continue
        path = sitearch / rel
        if not path.is_file():
            print(f"FAIL: RECORD references missing file {rel}", file=sys.stderr)
            return 1
        data = path.read_bytes()
        actual_hash = _record_hash(data)
        if actual_hash != recorded_hash:
            print(
                f"FAIL: hash mismatch for {rel}: "
                f"RECORD has {recorded_hash}, actual is {actual_hash}",
                file=sys.stderr,
            )
            return 1
        if int(recorded_size) != len(data):
            print(
                f"FAIL: size mismatch for {rel}: "
                f"RECORD has {recorded_size}, actual is {len(data)}",
                file=sys.stderr,
            )
            return 1
    if not saw_record_row:
        print(
            "FAIL: RECORD is missing its self-row (PEP 376 requires "
            "'<dist-info>/RECORD,,')",
            file=sys.stderr,
        )
        return 1

    print(
        f"OK: wignernj-{version} dist-info under {sitearch} is well-formed; "
        f"{len(record_lines) - 1} files verified against RECORD."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
