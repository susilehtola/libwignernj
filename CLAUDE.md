# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this library is

**libwignernj** — exact evaluation of Wigner 3j, 6j, 9j symbols, Clebsch-Gordan coefficients, Racah W coefficients, and Gaunt coefficients in C99, following the prime-factorization algorithm of Johansson & Forssén (SIAM J. Sci. Comput. 38(1), A376–A384, 2016; doi:10.1137/15M1021908). The key property: all intermediate arithmetic is exact integer arithmetic; floating-point conversion happens only at the final step. Results are accurate to the last bit of the chosen output precision.

Language interfaces: C (primary), C++ (header-only wrapper that links against `libwignernj`), Python (CPython extension), Fortran 90 (iso_c_binding).

All public API arguments use `2*j` integers (so `j=3/2` is passed as `tj=3`). This avoids floating-point half-integer representation throughout.

## Build commands

```bash
# Standard build (out-of-tree, recommended)
cmake -B build && cmake --build build

# With all options explicit
cmake -B build -DBUILD_FORTRAN=ON -DBUILD_TESTS=ON -DBUILD_CXX_TESTS=ON
cmake --build build

# Run all tests
ctest --test-dir build

# Run a single C test binary
./build/tests/test_3j

# Run a single named test via ctest
ctest --test-dir build -R test_3j

# Python extension (pip, self-contained — does not need a prior CMake build)
pip install -e .
python -m pytest tests/python/

# Python extension via CMake
cmake -B build -DBUILD_PYTHON=ON && cmake --build build

# MPFR arbitrary-precision interface
cmake -B build -DBUILD_MPFR=ON && cmake --build build
```

## Architecture: computation pipeline

Every public symbol function follows the same pipeline:

1. **Selection rules** — return `is_zero=1` immediately if violated.
2. **`wigner*_exact()`** (internal, in `src/wigner3j.c` etc.) — produces a `wigner_exact_t`:
   - Builds the outer sqrt-factor as a `pfrac_t` (prime-factored rational representing the argument under the outer square root — triangle Δ coefficients and m-dependent factorials for 3j; four triangle Δ's for 6j).
   - Two-pass Racah sum: pass 1 finds `lcm_exp[i]` (max Legendre valuation per prime across all sum terms); pass 2 converts each term to a `bigint_t` scaled integer and accumulates the signed sum.
   - Calls `pfrac_to_sqrt_rational()` to split the outer pfrac into `int_num`, `int_den`, `sqrt_num`, `sqrt_den` bigints.
3. **`wigner_exact_to_{float,double,long_double}()`** (in `src/wigner_exact.c`) — the only floating-point step: `sign * sum * int_num / int_den * sqrt(sqrt_num / sqrt_den)`.

The 9j symbol is implemented as a sum over the intermediate quantum number `k` of products of three 6j exact-path evaluations. At each `k` the three k-dependent Δ factors each appear twice (making them Δ² = rational), and the six k-independent Δ factors form the outer sqrt(C). This enables the same `sqrt(C) × exact_integer` structure as 3j and 6j.

Clebsch-Gordan and Racah W are thin wrappers over the Wigner symbols. Gaunt has its own exact arithmetic pipeline: it builds a single combined `pfrac_t` for `[Δ(l1,l2,l3)]² × (li!)² × m-factorials × (2li+1)/4`, runs two independent Racah sums (for m=(0,0,0) and m=(m1,m2,m3)), multiplies the sums as bigints, and applies `1/sqrt(π)` only at the final float conversion.

## Module responsibilities

| File | Role |
|------|------|
| `src/primes.c` | Sieve of Eratosthenes; `legendre_valuation(n, pi)` = v_p(n!) via Legendre's formula |
| `src/bigint.c` | Unsigned multiword integer (little-endian `uint64_t` words) + sign; `bigint_to_{float,double,long_double}` with correct IEEE 754 round-to-nearest |
| `src/bigint_arith.h` | 64-bit arithmetic primitives (mul/add/sub/div with carry); native `__uint128_t` path on GCC/Clang/ICC, pure-C99 fallback (32×32 partial products, 64/32 long division) on other compilers including MSVC. `-DBIGINT_FORCE_PORTABLE` forces the fallback. |
| `src/pfrac.c` | Prime-factored rational: signed `int exp[]` indexed by prime table index; `pfrac_mul_factorial` / `pfrac_div_factorial`; `pfrac_to_sqrt_rational` |
| `src/wigner_exact.c` | `wigner_exact_t` struct and `wigner_exact_to_*` conversion |
| `src/wigner3j.c` | `wigner3j_exact()` + public `wigner3j_f/wigner3j/wigner3j_l` |
| `src/wigner6j.c` | `wigner6j_exact()` + public variants |
| `src/wigner9j.c` | `wigner9j_exact()` (sum over k of three 6j exact-path products) + public variants |
| `src/clebsch.c` | `clebsch_gordan*`: CG via 3j, formula `(-1)^(j1-j2+M) sqrt(2J+1) * 3j(j1,j2,J;m1,m2,-M)` |
| `src/racah.c` | `racah_w*`: Racah W via 6j, formula `(-1)^(j1+j2+J+j3) * 6j{j1,j2,j12;j3,J,j23}` |
| `src/gaunt.c` | `gaunt*`: own exact pipeline — combined pfrac for [Δ]²·factorials·normalization, two Racah sums multiplied as bigints, then `÷sqrt(π)` at float step |
| `include/wigner.h` | Public C API — all functions, all precisions |
| `include/wigner_mpfr.h` | MPFR API — requires `BUILD_MPFR=ON`; set precision on `rop` before calling |
| `include/wigner.hpp` | C++11 header-only wrapper (links `wignernj`): `wigner::symbol3j<T>()`, real-valued overloads, `std::invalid_argument` for non-half-integer inputs |
| `src/fortran/wigner_f90.F90` | Fortran module `wigner`: raw `bind(c)` interfaces + `w3j/w6j/w9j/wcg/wracah/wgaunt` real-valued wrappers |
| `src/python/wignermodule.c` | CPython extension `_wigner`: parses int/float/Fraction, `precision=` kwarg |
| `wigner/__init__.py` | Re-exports from `_wigner` |

## Key conventions

- **`exp[]` storage in `pfrac_t`**: exponent of prime `p_i` in the *argument* of the outer sqrt is stored directly as `exp[i] = v_p(numerator) - v_p(denominator)`. For the rational sum terms (no sqrt), exponents are exact integers stored as-is. The split is: even `exp[i]` → `p_i^(exp[i]/2)` goes into `int_num`/`int_den`; odd `exp[i]` → remaining `sqrt(p_i)` goes into `sqrt_num`/`sqrt_den`.
- **`wigner_exact_t` value**: `sign * sum * int_num / int_den * sqrt(sqrt_num / sqrt_den)` where `int_den` absorbs both the outer-sqrt integer denominator and the Racah sum LCM.
- **`bigint_to_double` rounding**: extracts `DBL_MANT_DIG` bits with round-to-nearest-even using explicit round/sticky bits; `bigint_to_long_double` uses `LDBL_MANT_DIG` (80-bit extended: 64 on x86-64; double-equiv: 53 on MSVC/ARM; 128-bit quad: 113 on aarch64/POWER).
- **Factorial arguments as integers**: because the triangle condition forces `tj1+tj2+tj3` even, all `(tja ± tjb ± tmc)/2` expressions in the Racah sum are guaranteed integers. No runtime checks needed.
- **SPDX headers**: every source file begins with `SPDX-License-Identifier: BSD-3-Clause` and `Copyright (c) 2026 Susi Lehtola` in the appropriate comment syntax.
