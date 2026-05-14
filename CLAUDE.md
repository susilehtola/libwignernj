# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this library is

**libwignernj** — exact evaluation of Wigner 3j, 6j, 9j symbols, Clebsch-Gordan coefficients, Racah W coefficients, Fano X-coefficients, and Gaunt coefficients in C99, following the prime-factorization technique introduced for the angular-momentum coefficients by Dodds & Wiechers (Comput. Phys. Commun. 4, 268, 1972; doi:10.1016/0010-4655(72)90019-7) and refined in subsequent work, combined with the multiword-integer Racah sum of Johansson & Forssén (SIAM J. Sci. Comput. 38(1), A376–A384, 2016; doi:10.1137/15M1021908). The key property: all intermediate arithmetic is exact integer arithmetic; floating-point conversion happens only at the final step. Results are accurate to the last bit of the chosen output precision.

Language interfaces: C (primary), C++ (header-only wrapper `wignernj.hpp` that links against `libwignernj`), Python (CPython extension `wignernj` — `pip install wignernj`, `import wignernj`), Fortran 90 (`module wignernj` from `libwignernj_f03`, via `iso_c_binding`). The `wignernj` name is used uniformly across every binding; the only `wigner`-prefixed identifiers that survive are the C math-symbol functions (`wigner3j`, `wigner6j`, `wigner9j`), which name the mathematical Wigner symbols rather than the library namespace.

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

# libquadmath / __float128 binary128 interface
cmake -B build -DBUILD_QUADMATH=ON && cmake --build build

# FLINT bigint backend (closes large-j speed gap; replaces in-tree schoolbook+Karatsuba)
cmake -B build -DBUILD_FLINT=ON && cmake --build build
```

## Architecture: computation pipeline

Every public symbol function follows the same pipeline:

1. **Selection rules** — return `is_zero=1` immediately if violated.
2. **`wigner*_exact()`** (internal, in `src/wigner3j.c` etc.) — produces a `wignernj_exact_t`:
   - Builds the outer sqrt-factor as a `pfrac_t` (prime-factored rational representing the argument under the outer square root — triangle Δ coefficients and m-dependent factorials for 3j; four triangle Δ's for 6j).
   - Two-pass Racah sum: pass 1 finds `lcm_exp[i]` (max Legendre valuation per prime across all sum terms); pass 2 converts each term to a `bigint_t` scaled integer and accumulates the signed sum.
   - Calls `pfrac_to_sqrt_rational()` to split the outer pfrac into `int_num`, `int_den`, `sqrt_num`, `sqrt_den` bigints.
3. **`wignernj_exact_to_{float,double,long_double,float128}()`** (in `src/wigner_exact.c`) — the only floating-point step: `sign * sum * int_num / int_den * sqrt(sqrt_num / sqrt_den)`. The `_float128` variant is gated on `WIGNERNJ_HAVE_QUADMATH` (set by `-DBUILD_QUADMATH=ON`).

The 9j symbol is implemented as a sum over the intermediate quantum number `k` of products of three 6j exact-path evaluations. At each `k` the three k-dependent Δ factors each appear twice (making them Δ² = rational), and the six k-independent Δ factors form the outer sqrt(C). This enables the same `sqrt(C) × exact_integer` structure as 3j and 6j.

Clebsch-Gordan, Racah W, and Fano X are thin wrappers over the Wigner symbols (CG and Racah W over 3j and 6j respectively; Fano X over 9j with four `sqrt(2j+1)` factors folded into the existing prime-decomposed outer-sqrt tuple). Gaunt has its own exact arithmetic pipeline: it builds a single combined `pfrac_t` for `[Δ(l1,l2,l3)]² × (li!)² × m-factorials × (2li+1)/4`, runs two independent Racah sums (for m=(0,0,0) and m=(m1,m2,m3)), multiplies the sums as bigints, and applies `1/sqrt(π)` only at the final float conversion.

## Module responsibilities

| File | Role |
|------|------|
| `src/primes.c` | `legendre_valuation(n, pi)` = v_p(n!) via Legendre's formula. The prime list is hard-coded into the compiled library by `#include "prime_table.inc"`, generated at build time by `tools/gen_prime_table.py` (the runtime sieve was removed in 0.2.0) |
| `src/bigint.c` | Default in-tree backend: unsigned multiword integer (little-endian `uint64_t` words) + sign; schoolbook multiplication below `KARATSUBA_THRESHOLD` (default 32 limbs), Karatsuba above; `bigint_to_{float,double,long_double}` with correct IEEE 754 round-to-nearest, plus `bigint_to_float128` / `bigint_frexp_q` (gated on `WIGNERNJ_HAVE_QUADMATH`) using top-3-words Horner-style binary128 evaluation |
| `src/bigint_flint.c` | Optional FLINT backend (gated on `BUILD_FLINT`): the same bigint API delegated to FLINT's `fmpz_t`. Floating-point conversions go through MPFR for correct rounding; binary128 uses `mpfr_get_float128`. Replaces (does not augment) `bigint.c` when enabled |
| `src/bigint_arith.h` | 64-bit arithmetic primitives (mul/add/sub/div with carry); native `__uint128_t` path on GCC/Clang/ICC, pure-C99 fallback (32×32 partial products, 64/32 long division) on other compilers including MSVC. `-DBIGINT_FORCE_PORTABLE` forces the fallback. |
| `src/pfrac.c` | Prime-factored rational: signed `int exp[]` indexed by prime table index; `pfrac_mul_factorial` / `pfrac_div_factorial`; `pfrac_to_sqrt_rational` |
| `src/wigner_exact.c` | `wignernj_exact_t` struct and `wignernj_exact_to_*` conversion |
| `src/wigner3j.c` | `wigner3j_exact()` + public `wigner3j_f/wigner3j/wigner3j_l` (and `wigner3j_q` when built with `BUILD_QUADMATH=ON`) |
| `src/wigner6j.c` | `wigner6j_exact()` + public variants |
| `src/wigner9j.c` | `wigner9j_exact()` (sum over k of three 6j exact-path products) + public variants |
| `src/clebsch.c` | `clebsch_gordan*`: CG via 3j, formula `(-1)^(j1-j2+M) sqrt(2J+1) * 3j(j1,j2,J;m1,m2,-M)` |
| `src/racah.c` | `racah_w*`: Racah W via 6j, formula `(-1)^(j1+j2+J+j3) * 6j{j1,j2,j12;j3,J,j23}` |
| `src/fano_x.c` | `fano_x*`: Fano X via 9j, formula `sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * 9j{j1,j2,j12;j3,j4,j34;j13,j24,J}` |
| `src/gaunt.c` | `gaunt*`: own exact pipeline — combined pfrac for [Δ]²·factorials·normalization, two Racah sums multiplied as bigints, then `÷sqrt(π)` at float step |
| `src/real_ylm_in_complex_ylm.c` | `wignernj_real_ylm_in_complex_ylm*`: fills the column-major `(2l+1)×(2l+1)` unitary `C` for `S = C Y` under the same real-Y convention as `gaunt_real`; direct entry-by-entry fill of `±1`, `±1/√2`, `±i/√2` — no Racah sum, no factorials |
| `include/wignernj.h` | Public C API — all functions, `_f`/`/`_l` precisions |
| `include/wignernj_quadmath.h` | libquadmath API — requires `BUILD_QUADMATH=ON`; declares `_q` (`__float128`) variant of every public symbol |
| `include/wignernj_mpfr.h` | MPFR API — requires `BUILD_MPFR=ON`; set precision on `rop` before calling |
| `include/wignernj.hpp` | C++11 header-only wrapper (links `wignernj`): `wignernj::symbol3j<T>()`, real-valued overloads, `std::invalid_argument` for non-half-integer inputs |
| `src/fortran/wignernj_f90.F90` | Fortran module `wignernj`: raw `bind(c)` interfaces + `w3j/w6j/w9j/wcg/wracahw/wfanox/wgaunt/wgaunt_real` real-valued scalar wrappers and `wreal_ylm_in_complex_ylm(l, c_out)` matrix-filling subroutine; `_q` interfaces and `w3jq`/etc.\ wrappers plus `wreal_ylm_in_complex_ylmq` gated on `WIGNERNJ_HAVE_QUADMATH` |
|  `src/python/wignernjmodule.c` | CPython extension `_wignernj`: parses int/float/Fraction, `precision=` kwarg |
| `wignernj/__init__.py` | Re-exports from `_wignernj` |

## Key conventions

- **`exp[]` storage in `pfrac_t`**: exponent of prime `p_i` in the *argument* of the outer sqrt is stored directly as `exp[i] = v_p(numerator) - v_p(denominator)`. For the rational sum terms (no sqrt), exponents are exact integers stored as-is. The split is: even `exp[i]` → `p_i^(exp[i]/2)` goes into `int_num`/`int_den`; odd `exp[i]` → remaining `sqrt(p_i)` goes into `sqrt_num`/`sqrt_den`.
- **`wignernj_exact_t` value**: `sign * sum * int_num / int_den * sqrt(sqrt_num / sqrt_den)` where `int_den` absorbs both the outer-sqrt integer denominator and the Racah sum LCM.
- **`bigint_to_double` rounding**: extracts `DBL_MANT_DIG` bits with round-to-nearest-even using explicit round/sticky bits; `bigint_to_long_double` uses `LDBL_MANT_DIG` (80-bit extended: 64 on x86-64; double-equiv: 53 on MSVC/ARM; 128-bit quad: 113 on aarch64/POWER); `bigint_to_float128` Horner-evaluates the top three 64-bit words in `__float128` (192 input bits feeding a 113-bit mantissa), good to within 2 ulp at quad precision.
- **Factorial arguments as integers**: because the triangle condition forces `tj1+tj2+tj3` even, all `(tja ± tjb ± tmc)/2` expressions in the Racah sum are guaranteed integers. No runtime checks needed.
- **SPDX headers**: every source file begins with `SPDX-License-Identifier: BSD-3-Clause` and `Copyright (c) 2026 Susi Lehtola` in the appropriate comment syntax.

## Project policy

These constraints apply to every change, not just to typical edits.

- **Semantic versioning.** The library follows [Semantic Versioning 2.0.0](https://semver.org). The project version (set by the top-level `project(wignernj VERSION ...)` declaration in `CMakeLists.txt`) is `MAJOR.MINOR.PATCH`:
  - **PATCH** (x.y.Z): backwards-compatible bug fixes — corrected results, performance improvements that don't change observable behaviour, internal refactoring with no public-API impact.
  - **MINOR** (x.Y.z): backwards-compatible additions — new symbol families, new precision variants, new build options, new public functions; no removals or signature changes.
  - **MAJOR** (X.y.z): incompatible API changes — removed functions, changed signatures, behavioural changes that break existing callers.

  The library is currently 0.x.y, i.e. pre-1.0; per the semver convention, public-API stability is not yet guaranteed and any minor bump may break compatibility. Callers should pin to a specific minor version until 1.0.0 ships. Bumping the version is part of the change that introduces the new behaviour, not a separate later commit.

- **Clean-room implementation.** libwignernj is BSD-3-Clause and depends on no GPL/LGPL competitor source. Do not read the source of WIGXJPF, FASTWIGXJ, or any other GPL/LGPL coupling-coefficient library — the library is to remain a derivation from the published Dodds–Wiechers prime-factorization scheme, the Johansson–Forssén multiword-integer Racah sum, and other openly described methods, never from competitor code.

- **Out of scope: intra-call threading.** OpenMP, pthreads, threaded BLAS, and any other multi-threading inside one library call are excluded by the embeddability mandate (single-threaded callers, host-managed thread pools, and language bindings whose own threading model is unknown all need to be able to use the library without surprise multi-threading). Per-thread caches and `wignernj_warmup_to` make concurrent calls from a host thread pool safe and scale well — that is the right threading-related work. Don't propose intra-call parallelism.

- **MAX_FACTORIAL_ARG extensibility.** The default-build prime-table ceiling (MAX_FACTORIAL_ARG = 20020 derived from the sieve limit) is one configuration; consumers regenerate `src/prime_table.inc` with a larger sieve when they need higher j. Optimisations must scale to that — don't bake "fits in uint64", "fits in 32 bits", or "valuation ≤ 127" assumptions tied only to the default ceiling unless they're guarded by a `#error` or `_Static_assert` that fires at build time when the bound is exceeded.

- **CI path coverage on every architecture.** The `-DBIGINT_NO_DIVQ` and `-DBIGINT_FORCE_PORTABLE` cells on x86-64 are *smoke tests* for the alternative dispatch arms — they exercise the same code paths the library selects in production on aarch64, ppc64le, MSVC, etc. They are not a substitute for testing on the native architecture where each path is selected. Any conditional code path needs CI coverage on every architecture where it can be selected.

- **Bench-driven optimisation.** Every performance-relevant change must be benched against the parent commit using paired alternating runs (alternate the variant under test with main on each iteration to control for thermal/frequency drift), and the headline numbers go in the commit message. Single-shot timings at small j are noise-dominated; use either many-trial microbenchmarks (`bench_term_cache.c`) or longer-budget paired sweeps (`bench_sweep.c`).

- **Library-versioning benchmarks belong in the repo; comparative benches against external libraries don't.** The `benchmarks/` directory hosts library-versioning microbenches and profile drivers — `bench_term_cache.c`, `bench_sweep.c`, `bench_mul.c`, `bench_div128.c`, `profile_*.c`, etc. — used to validate libwignernj changes against prior versions of libwignernj itself.  They have no external dependencies and are part of the project's bench-driven optimisation methodology.  Head-to-head benches against WIGXJPF, GSL, or FLINT belong with the paper they support (as supplementary material), not in this repo: they pull in dependency surface libwignernj itself doesn't have, and they exist to publish a comparison rather than to develop the library.

- **Known dead-end optimisations.** The following have all been implemented end-to-end, paired-benched, and reverted as net-negative on libwignernj's typical operating range (j ≤ 500); see `docs/optimization_notes.md` for the bench data and the underlying mechanism in each case. Do not re-propose them without new methodology that addresses the failure modes recorded there.
  - **Toom-3 multiplication** for the in-tree bigint backend. Asymptotic crossover with Karatsuba is past where libwignernj operates; the additional code size pollutes icache for the rest of the binary and slows the typical-j hot path. FLINT covers the very-large-j case.
  - **Compressed `pfrac_t.exp` / factorial-cache rows** (split int8/int16/int32 tier on the prime index). The hot-path dispatch cost (extra struct dereferences, two loops instead of one) outweighs the cache-pressure savings at j ≤ 500.
  - **Small-bigint specialisation** for Pass 2 (stack-allocated `uint64_t scaled[N]` + static-inline `bigint_mul_u64_small` etc., dispatching on the seed size). LTO is already inlining the existing primitives; the duplicate code path bloats icache and regresses the typical-j workload.

  The shape of optimisations that *did* work (Möller-Granlund division, Gaunt Pass-2 ratio recurrence, Gaunt Pass-1 cache skip) was different: each one *removed work* (replaced an O(n²) per-term sweep with an O(n) recurrence, gated a slow path behind a faster alternative on backends without hardware divq). They didn't add code complexity to the per-call hot path. Optimisation proposals should match this shape.

- **Code and documentation must stay consistent.** Every concrete fact in the docs (function name, build option, file path, symbol family, byte size, j-limit, performance number, attribution) must match the current source. When a change touches a public surface — a renamed symbol, a new build option, a new public family — sweep `README.md`, `CHANGELOG.md`, `CLAUDE.md`, `docs/reference.md`, `examples/README.md`, the public-header docstrings, and the wrapper docs in the same PR; don't defer the doc fix to a later commit. Treat every doc claim as a hypothesis and verify against the live source before pushing.

- **Run a full repository consistency audit before tagging a release.** Sweep every tracked file — sources, comments, public-header docstrings, every doc, every example, every test, every build-system / CI-config file. Audit by claim: every concrete fact is a hypothesis; verify against live source / build / CI / test output. The minimum coverage list: versions across `CMakeLists.txt` / `pyproject.toml` / `CHANGELOG.md`; public-symbol enumerations; build-option tables; factorial / sieve ceilings (`MAX_FACTORIAL_ARG`, `PRIME_SIEVE_LIMIT`); removed-symbol references; the CLAUDE.md "Module responsibilities" table; latent bugs in the diff against the previous tag. Past misses caught only by focused audits include `g_prime_index` listed in a CHANGELOG export-symbols entry after removal, `bench_compare.c` "siblings" claimed moved out of the repo when never tracked, `BUILD_FLINT=ON Python build path` (should have been `BUILD_PYTHON=ON`), the `wigner_warmup` entry inheriting a stale 19999 ceiling from a now-stale public-header docstring, and `src/primes.c` described as containing a Sieve of Eratosthenes when the runtime sieve was removed in 0.2.0.

- **Amend `CHANGELOG.md` as part of every pull request.** The matching `### Added` / `### Changed` / `### Fixed` / `### Deprecated` / `### Removed` / `### Security` entry under `## [Unreleased]` lands in the same commit (or at least the same PR) as the code change it describes. Don't defer the CHANGELOG to a release-time housekeeping commit. Pure-internal refactors with no observable behaviour change can skip the CHANGELOG, but err on the side of including an entry when in doubt; anything that touches public headers, build options, supported toolchains, build instructions, or test infrastructure that contributors interact with should get one. Release tagging is then mechanical: rename `[Unreleased]` → `[X.Y.Z] – DATE`, add a fresh empty `[Unreleased]`, update the compare-link footnote.

- **Every commit must build without warnings.** No `-Wall -Wextra` warnings on any supported toolchain (gcc, clang, Apple Clang, MSVC, gfortran, Intel ifx). New warnings are real bugs (or near-bugs) until proven otherwise: real bugs get fixed; standard-defying intentional code (e.g. `__int128` firing `-Wpedantic` in ISO C99 mode) gets the canonical narrowly-scoped workaround with a one-line comment so a reader doesn't think the cast is a mistake; compiler-version regressions get the new canonical idiom (e.g. `c_float128` for the gfortran-16 `-Wc-binding-type` on `real(real128)` `bind(c)` returns). Don't downgrade the compiler. The pre-tag and per-PR audits should rebuild with `-Wall -Wextra` and confirm the diagnostic output is empty before the tag/PR lands.

- **Workflow.** Three smaller conventions worth flagging:
  - **Amend, don't stack, fixups for the same topic.** A CI failure or a typo discovered in the just-pushed commit gets amended into that commit (followed by a `--force-with-lease` push), not stacked as a separate "fix typo" follow-up.
  - **CircleCI for OSS arm64 + macOS off-load,** not Cirrus. Cirrus shuts down June 2026; CircleCI is the durable alternative for native arm64 Linux + musl + i686 cells when GHA quotas don't suffice. (macOS stayed on GHA in the 0.4.0 cycle because CircleCI's macOS executor is too slow on the OSS plan.)
  - **`export VAR=val` for shell env vars,** not inline `VAR=val command`. Preferred in spec `%check`/`%build` blocks, README snippets, and CI recipes; more reliable across macro expansions and subprocesses.
