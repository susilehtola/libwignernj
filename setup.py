# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
import sys

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

# MSVC has no -std=c99 / -O2 spelling and bundles libm into the
# default C runtime; the GCC/Clang flags would either be passed
# silently as cl warnings or, in the libm case, produce a hard
# linker error (LNK1181: cannot open input file 'm.lib').
if sys.platform == "win32":
    extra_compile_args = ["/O2"]
    libraries = []
else:
    extra_compile_args = ["-O2", "-std=c99"]
    libraries = ["m"]

ext = Extension(
    name="wigner._wigner",
    sources=sources,
    include_dirs=["src", "include"],
    extra_compile_args=extra_compile_args,
    libraries=libraries,
)

setup(
    name="wigner",
    version="0.2.0",
    description="Exact Wigner 3j/6j/9j symbols and related coefficients via prime factorization",
    packages=["wigner"],
    ext_modules=[ext],
    python_requires=">=3.7",
)
