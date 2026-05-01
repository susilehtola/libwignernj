# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
from setuptools import setup, Extension

sources = [
    "src/python/wignermodule.c",
    "src/xalloc.c",
    "src/primes.c",
    "src/bigint.c",
    "src/pfrac.c",
    "src/wigner_exact.c",
    "src/wigner3j.c",
    "src/wigner6j.c",
    "src/wigner9j.c",
    "src/clebsch.c",
    "src/racah.c",
    "src/gaunt.c",
]

ext = Extension(
    name="wigner._wigner",
    sources=sources,
    include_dirs=["src", "include"],
    extra_compile_args=["-O2", "-std=c99"],
    libraries=["m"],
)

setup(
    name="wigner",
    version="0.1.0",
    description="Exact Wigner 3j/6j/9j symbols and related coefficients",
    packages=["wigner"],
    ext_modules=[ext],
    python_requires=">=3.7",
)
