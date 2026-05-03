# Changelog

All notable changes to **libwignernj** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Optional FLINT bigint backend (`-DBUILD_FLINT=ON`) that replaces
  the in-tree schoolbook multiword integer with a thin wrapper
  around FLINT's `fmpz_t`. The asymptotic motivation is
  sub-quadratic multiplication (Karatsuba / Toom-Cook /
  Schönhage--Strassen via FLINT/GMP), which closes most of the
  large-`j` performance gap to `WIGXJPF`. Floating-point
  conversions are routed through MPFR (correct round-to-nearest-
  even at every IEEE 754 binary precision); the binary128
  conversion uses `mpfr_get_float128` (MPFR ≥ 4.1.0 with
  `--enable-float128`). Bit-identical output against the
  schoolbook is verified in CI via a dedicated `Ubuntu / GCC /
  FLINT backend` matrix cell. The default build is unchanged and
  remains dependency-free.
- Fano X-coefficient (`fano_x`, `fano_x_f`, `fano_x_l`, `fano_x_q`,
  `fano_x_mpfr`), implemented as a thin wrapper over the 9j exact
  pipeline that folds the four `sqrt(2j+1)` factors into the existing
  prime-decomposed outer-sqrt tuple via the same scheme used by the
  Clebsch--Gordan `sqrt(2J+1)` factor. Fortran (`wfanox`, `wfanoxq`),
  C++ (`wigner::fanox<T>`), and Python (`wigner.fano_x`) bindings
  follow the same pattern as the other derived symbols.
- IEEE 754 binary128 (`__float128`) back-end via libquadmath, enabled with
  `-DBUILD_QUADMATH=ON`. Adds `wigner3j_q`, `wigner6j_q`, `wigner9j_q`,
  `clebsch_gordan_q`, `racah_w_q`, `gaunt_q`, and `gaunt_real_q` in a new
  public header `wigner_quadmath.h`. Conversion path
  `wigner_exact_to_float128` mirrors the long-double variant; the underlying
  `bigint_to_float128` Horner-evaluates the top three 64-bit words in
  `__float128`, feeding 192 input bits into a 113-bit mantissa.
- Fortran 90 quadruple-precision interface gated on `WIGNERNJ_HAVE_QUADMATH`:
  raw `bind(c)` interfaces (`wigner3j_q`, …, `gaunt_real_q`) and real-valued
  convenience wrappers `w3jq`, `w6jq`, `w9jq`, `wcgq`, `wracahwq`, `wgauntq`,
  `wgaunt_realq`, all returning `real(real128)` and exposed through the same
  `use wigner` that imports the double-precision routines.
- CMake auto-discovery of `<quadmath.h>` on toolchains where the gcc-shipped
  header is not on the C compiler's default search path (Clang/Linux): a
  fallback probe queries `gcc -print-file-name=include` and propagates the
  directory to consumers.
- Per-term accuracy benchmark drivers (`benchmarks/maxerr_3j.c`,
  `benchmarks/sweep_3j_j100.c`, `benchmarks/sweep_3j_j100_reference.py`)
  that compare libwignernj/WIGXJPF/GSL term by term across the full
  `m`-range, replacing the cancellation-prone loop sum used by
  `bench_compare`.
- Microsoft Windows build via the MSVC toolchain (C, C++, and Python
  interfaces), exercised on every push by the `windows-latest` GitHub
  Actions runner. The Fortran 90 wrapper, the libquadmath back-end, and
  the GNU MPFR back-end are not built on Windows because the default
  MSVC toolchain ships with neither `gfortran`, `libquadmath`, nor
  `libmpfr`.
- CI matrix coverage of static and shared library builds (Linux GCC,
  Linux Clang, macOS Clang) alongside downstream `find_package` smoke
  tests for both library configurations and for the C-only and Fortran
  consumers.

### Fixed
- Windows MSVC link errors: `WINDOWS_EXPORT_ALL_SYMBOLS=ON` for
  auto-generated import libraries; explicit `__declspec(dllexport)`/
  `dllimport` decoration for the prime-table data symbols
  (`g_nprimes`, `g_primes`, `g_prime_index`); co-location of test
  executables and `wignernj.dll` in a single output directory so the
  Windows dynamic linker resolves the DLL at test time
  (`STATUS_DLL_NOT_FOUND`).
- `pyproject.toml` license metadata modernised to PEP 639 (SPDX-string
  `license` field, `license-files`, OSI classifier dropped) so
  `setuptools >= 77` no longer warns.
- Fortran wrapper now ships its own `libwignernj_f03.pc` (declares
  `Requires: libwignernj`, `Cflags: -I${moddir}` for the `.mod` file
  directory) so downstream Fortran consumers can locate both the
  `wigner` module and the C transitive dependency through pkg-config.
- Stop linking `libm` on Windows (where it does not exist); guard
  `target_link_libraries(... PRIVATE m)` and the corresponding
  pkg-config `Libs:` line behind `if(NOT WIN32)`.

## [0.2.0] – 2026-05-02

### Added
- Exact-arithmetic rewrites of the 9j and Gaunt coefficient pipelines, so
  every public symbol now goes through the same prime-factorization +
  multiword-integer Racah sum + final-cast pipeline as the 3j and 6j.
- `gaunt_real()` for the Gaunt coefficient over real spherical harmonics
  (Wikipedia/Condon–Shortley convention), implemented via a single complex
  Gaunt evaluation per call.
- Optional GNU MPFR back-end (`-DBUILD_MPFR=ON`) exposing `wigner3j_mpfr`,
  `wigner6j_mpfr`, `wigner9j_mpfr`, `clebsch_gordan_mpfr`, `racah_w_mpfr`,
  `gaunt_mpfr`, and `gaunt_real_mpfr`. Reuses the same exact tuple as the
  C-precision routines; precision is set on the output `mpfr_t` via
  `mpfr_init2`, rounding mode is the last argument.
- Pure-C99 fallback for the multiword-integer arithmetic (selected
  automatically on compilers without `__uint128_t`, forced via
  `-DBIGINT_FORCE_PORTABLE`); CI cell verifies bit-identical output on
  both code paths.
- GitHub Actions CI pipeline: Ubuntu/GCC (shared and static), Ubuntu/Clang,
  macOS/Clang, plus a no-optional-features cell. Downstream
  `find_package(wignernj)` smoke test included.
- Symmetry-oracle tests for 3j/6j/9j; allocation-failure injection harness
  (`test_oom`) using a fork-based `xmalloc` interceptor.
- Full API reference manual at `docs/reference.md`.
- Comparative benchmark harness in `benchmarks/` (libwignernj vs. WIGXJPF
  vs. GSL), with a Makefile that uses Fedora's `rpm -E %optflags` so
  all three libraries are built with identical compile flags.
- Public-API documentation of the phase conventions (Condon–Shortley for
  Clebsch–Gordan and for `Y_l^m`).

### Changed
- Allocations routed through `xmalloc`/`xrealloc` helpers that abort with
  a clear diagnostic instead of silently returning `NULL`.
- Prime sieve replaced by precomputed `static const` tables; eliminates
  the runtime sieve and the global-initialization race for concurrent
  callers.
- `wignerXj_exact()` callers now reuse a workspace (`bigint_ws_t`)
  pre-allocated once per call, removing the per-multiplication
  `malloc`/`realloc` traffic from the inner Racah-sum loop.
- `pfrac_t` tracks the highest active prime index so zero-tail iterations
  in inner loops are skipped, on top of the existing prime-by-prime
  representation.
- Fortran library renamed to `libwignernj_f03`; ELF `RUNPATH=$ORIGIN`
  baked in so co-located shared libraries resolve transitive dependencies
  on Linux distributions whose default linker uses
  `--no-copy-dt-needed-entries` (Ubuntu).
- C++ wrapper clarified as a header-only template that links against
  `libwignernj` rather than a stand-alone library.
- C standard moved from a global setting to per-target so consumers
  embedding the library are not forced into C99 themselves.

### Fixed
- `install(... DESTINATION)` error when configuring with
  `-DBUILD_MPFR=ON` against a fresh build tree.
- Diagnostic abort (rather than silent corruption) when the requested
  angular momenta exceed the prime-table ceiling.
- `test_cpp` long-double tolerance now scales with `LDBL_EPSILON`,
  fixing spurious failures on platforms whose `long double` has the
  same precision as `double`.
- macOS / Homebrew MPFR build: GMP is now pulled in alongside MPFR so
  the multi-prefix Homebrew layout finds `gmp.h`.

## [0.1.0] – 2026-05-01

### Added
- Initial release of **libwignernj** — exact evaluation of Wigner 3j,
  6j, and 9j symbols and of Clebsch–Gordan and Racah W coefficients in
  C99, following the prime-factorization scheme of Johansson and
  Forssén [SIAM J. Sci. Comput. 38, A376 (2016)]. All intermediate
  arithmetic is exact integer arithmetic; floating-point conversion
  happens only at the final step. Results are accurate to the last
  bit of the chosen output precision.
- Three IEEE 754 binary precisions (`float`, `double`, `long double`)
  for every public symbol.
- C++11 header-only wrapper, Fortran 90 module via `iso_c_binding`,
  and CPython extension module exposing the same routines.
- BSD 3-Clause licence.

[Unreleased]: https://github.com/susilehtola/libwignernj/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/susilehtola/libwignernj/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/susilehtola/libwignernj/releases/tag/v0.1.0
