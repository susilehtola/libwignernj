#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# End-to-end test of the libwignernj CMake package:
#   1. Configure, build, and install libwignernj into a temporary prefix.
#   2. Configure the downstream consumer in tests/cmake_downstream
#      (using -Dwignernj_DIR or CMAKE_PREFIX_PATH) and verify that
#      find_package(wignernj REQUIRED [COMPONENTS Fortran]) succeeds.
#   3. Build the C, C++, and (optionally) Fortran consumer executables.
#   4. Run them through ctest.
#
# Usage:
#     tests/cmake_downstream/run.sh [--with-fortran]
#
# Returns non-zero on any failure.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

WITH_FORTRAN=OFF
SHARED=ON
for arg in "$@"; do
    case "$arg" in
        --with-fortran) WITH_FORTRAN=ON ;;
        --static)       SHARED=OFF ;;
        *) echo "Unknown option: $arg" >&2; exit 2 ;;
    esac
done

WORK="$(mktemp -d -t wignernj-downstream.XXXXXX)"
trap 'rm -rf "$WORK"' EXIT
BUILD="$WORK/build"
PREFIX="$WORK/install"
DOWN="$WORK/downstream"

echo "=== 1. Build & install libwignernj (BUILD_SHARED_LIBS=$SHARED, BUILD_FORTRAN=$WITH_FORTRAN) ==="
cmake -S "$SOURCE_DIR" -B "$BUILD" \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DBUILD_SHARED_LIBS=$SHARED \
      -DBUILD_FORTRAN=$WITH_FORTRAN \
      -DBUILD_TESTS=OFF -DBUILD_CXX_TESTS=OFF
cmake --build "$BUILD" --parallel --verbose
cmake --install "$BUILD"

echo "=== 2. Configure the downstream consumer ==="
cmake -S "$SCRIPT_DIR" -B "$DOWN" \
      -DCMAKE_PREFIX_PATH="$PREFIX" \
      -DWITH_FORTRAN=$WITH_FORTRAN

echo "=== 3. Build the downstream consumer ==="
cmake --build "$DOWN" --parallel --verbose

echo "=== 4. Run the downstream tests ==="
ctest --test-dir "$DOWN" --output-on-failure
