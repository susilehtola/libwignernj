#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""Write PEP 376 / PEP 427 .dist-info metadata for a CMake-installed wignernj.

CMake's BUILD_PYTHON path drops `_wignernj.<EXT_SUFFIX>.so` and
`__init__.py` under ${Python3_SITEARCH}/wignernj/ but does not produce
the `.dist-info` directory that pip and importlib.metadata require.
This script writes that directory at `cmake --install` time so that

  * `pip list`                    sees the installed package
  * `pip uninstall wignernj`      knows what to remove
  * `importlib.metadata.version`  returns the installed version
  * `pip install wignernj` later  detects the existing install

The script is invoked via CMake `install(CODE ...)` with arguments

    install_python_metadata.py SITEARCH VERSION PYPROJECT INSTALLER

where SITEARCH is the install destination with DESTDIR pre-applied,
VERSION is the project version (from CMake PROJECT_VERSION), PYPROJECT
is the absolute path to pyproject.toml (single source of truth for
metadata fields), and INSTALLER is the string to write to INSTALLER
(typically `cmake`).

Metadata follows PEP 566 (METADATA), PEP 376 (RECORD, INSTALLER,
top_level.txt) and PEP 427 (WHEEL).  The wheel tag is computed from
the running Python (sysconfig + sys.implementation) so the .dist-info
matches whatever .so suffix CMake's Python3 module actually produced.
"""
from __future__ import annotations

import base64
import hashlib
import os
import sys
import sysconfig
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # type: ignore[import-not-found]


def _wheel_tag() -> str:
    """Compute the wheel filename triple `<python>-<abi>-<platform>`.

    Mirrors `packaging.tags.sys_tags()[0]` without taking a dependency
    on the `packaging` PyPI distribution.  The platform component is
    normalised the same way the wheel spec mandates (dashes/dots to
    underscores).
    """
    impl = sys.implementation.name
    impl_short = {"cpython": "cp", "pypy": "pp"}.get(impl, impl[:2])
    py = f"{impl_short}{sys.version_info.major}{sys.version_info.minor}"
    soabi = sysconfig.get_config_var("SOABI") or ""
    if soabi.startswith(("cpython-", "pypy")):
        # 'cpython-313-x86_64-linux-gnu' -> 'cp313'
        # 'pypy39-pp73-darwin' -> 'pp39_pp73'
        head = soabi.split("-", 2)[:2]
        abi = "_".join(head).replace("cpython_", "cp")
    else:
        abi = py
    plat = (sysconfig.get_platform() or "any").replace("-", "_").replace(".", "_")
    return f"{py}-{abi}-{plat}"


def _record_hash(data: bytes) -> str:
    """sha256 hash in the urlsafe-base64 form RECORD requires."""
    digest = hashlib.sha256(data).digest()
    return "sha256=" + base64.urlsafe_b64encode(digest).rstrip(b"=").decode("ascii")


def _flatten_authors(entries: list[dict]) -> tuple[str, str]:
    """Return (Author, Author-email) lines for PEP 566 METADATA."""
    names = ", ".join(e["name"] for e in entries if "name" in e)
    emails = ", ".join(
        f'{e["name"]} <{e["email"]}>' if "name" in e else e["email"]
        for e in entries
        if "email" in e
    )
    return names, emails


def _build_metadata(version: str, pyproject: dict) -> str:
    """Render METADATA as a PEP 566 message body."""
    proj = pyproject.get("project", {})
    lines: list[str] = []
    lines.append("Metadata-Version: 2.1")
    lines.append(f"Name: {proj.get('name', 'wignernj')}")
    lines.append(f"Version: {version}")
    if "description" in proj:
        lines.append(f"Summary: {proj['description']}")
    if "license" in proj:
        lic = proj["license"]
        if isinstance(lic, str):
            lines.append(f"License: {lic}")
        elif isinstance(lic, dict) and "text" in lic:
            lines.append(f"License: {lic['text']}")
    authors = proj.get("authors") or []
    if authors:
        names, emails = _flatten_authors(authors)
        if names:
            lines.append(f"Author: {names}")
        if emails:
            lines.append(f"Author-email: {emails}")
    if "requires-python" in proj:
        lines.append(f"Requires-Python: {proj['requires-python']}")
    for kw in proj.get("keywords", []):
        lines.append(f"Keywords: {kw}")
    for cls in proj.get("classifiers", []):
        lines.append(f"Classifier: {cls}")
    for label, url in (proj.get("urls") or {}).items():
        lines.append(f"Project-URL: {label}, {url}")
    readme = proj.get("readme")
    if isinstance(readme, str):
        readme_path = Path(readme)
        if not readme_path.is_absolute():
            # readme is relative to pyproject.toml's directory
            readme_path = Path(os.environ.get("WIGNERNJ_PYPROJECT_DIR", ".")) / readme_path
        if readme_path.is_file():
            lines.append("Description-Content-Type: text/markdown")
            body = readme_path.read_text(encoding="utf-8")
            return "\n".join(lines) + "\n\n" + body
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    if len(argv) != 5:
        print(
            f"usage: {argv[0]} SITEARCH VERSION PYPROJECT INSTALLER",
            file=sys.stderr,
        )
        return 2
    sitearch = Path(argv[1])
    version = argv[2]
    pyproject_path = Path(argv[3])
    installer = argv[4]

    # Single source of truth for human-readable metadata.
    with pyproject_path.open("rb") as f:
        pyproject = tomllib.load(f)
    os.environ["WIGNERNJ_PYPROJECT_DIR"] = str(pyproject_path.parent)

    pkg_dir = sitearch / "wignernj"
    if not pkg_dir.is_dir():
        print(
            f"error: package directory {pkg_dir} not found; "
            "did the wignernj install step run first?",
            file=sys.stderr,
        )
        return 1

    dist_info = sitearch / f"wignernj-{version}.dist-info"
    dist_info.mkdir(parents=True, exist_ok=True)

    metadata_text = _build_metadata(version, pyproject)
    (dist_info / "METADATA").write_text(metadata_text, encoding="utf-8")

    wheel_text = (
        "Wheel-Version: 1.0\n"
        "Generator: cmake (libwignernj install_python_metadata.py)\n"
        "Root-Is-Purelib: false\n"
        f"Tag: {_wheel_tag()}\n"
    )
    (dist_info / "WHEEL").write_text(wheel_text, encoding="utf-8")

    (dist_info / "INSTALLER").write_text(installer + "\n", encoding="utf-8")
    (dist_info / "top_level.txt").write_text("wignernj\n", encoding="utf-8")

    # RECORD lists every file in wignernj/ and every dist-info file
    # except RECORD itself, with sha256 hash and byte size.  The RECORD
    # row for RECORD has empty hash and size, per PEP 376.
    record_lines: list[str] = []
    for path in sorted(pkg_dir.rglob("*")):
        if path.is_file():
            rel = path.relative_to(sitearch)
            data = path.read_bytes()
            record_lines.append(f"{rel.as_posix()},{_record_hash(data)},{len(data)}")
    for name in ("METADATA", "WHEEL", "INSTALLER", "top_level.txt"):
        f = dist_info / name
        rel = f.relative_to(sitearch)
        data = f.read_bytes()
        record_lines.append(f"{rel.as_posix()},{_record_hash(data)},{len(data)}")
    record_lines.append(f"{(dist_info / 'RECORD').relative_to(sitearch).as_posix()},,")
    (dist_info / "RECORD").write_text("\n".join(record_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
