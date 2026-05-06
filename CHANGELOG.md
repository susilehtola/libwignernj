# Changelog

All notable changes to **libwignernj** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- gfortran 16's new `-Wc-binding-type` diagnostic, which fired once per
  `_q` `bind(c)` interface (`wigner3j_q`, `wigner6j_q`, `wigner9j_q`,
  `clebsch_gordan_q`, `racah_w_q`, `fano_x_q`, `gaunt_q`, `gaunt_real_q`)
  because `real128` from `iso_fortran_env` is not a formally
  C-interoperable kind. The Fortran module now imports `c_float128` from
  `iso_c_binding` (a gfortran/Intel ifx extension that maps directly to
  C's `__float128`) and uses it for every quadmath `bind(c)` declaration.
  No public-API change — `c_float128` and `real128` are the same physical
  binary128 kind on every supported toolchain, so callers that pass or
  receive `real(real128)` continue to work.

## [0.3.0] – 2026-05-06

### Added
- Optional FLINT bigint backend (`-DBUILD_FLINT=ON`) that replaces
  the in-tree multiword integer (schoolbook with Karatsuba above a
  measured 32-limb crossover) with a thin wrapper around FLINT's
  `fmpz_t`. The asymptotic motivation is sub-quadratic multiplication
  beyond the Karatsuba range (Toom-Cook / Schönhage--Strassen via
  FLINT/GMP), which closes the remaining large-`j` performance gap to
  `WIGXJPF`. Floating-point conversions are routed through MPFR
  (correct round-to-nearest-even at every IEEE 754 binary precision);
  the binary128 conversion uses `mpfr_get_float128` (MPFR ≥ 4.1.0
  with `--enable-float128`). Bit-identical output against the
  in-tree backend is verified in CI via a dedicated
  `Ubuntu / GCC / FLINT backend` matrix cell. The default build is
  unchanged and remains dependency-free.
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
- Library-versioning microbench harnesses and profile drivers in
  `benchmarks/` (`bench_term_cache.c`, `bench_sweep.c`, `bench_mul.c`,
  `bench_div128.c`, `bench_div128_gm.c`, `profile_3j_4000.c`,
  `profile_6j_9j.c`) used to validate libwignernj changes against
  prior versions of itself. Comparative benchmarks against external
  libraries (`bench_compare.c` plus its Makefile) were removed from
  the repo and now ship with the paper as supplementary material,
  since they link against WIGXJPF and GSL and exist to publish a
  comparison rather than to develop the library.
- `examples/` directory with a single-file demonstration of every
  public symbol in C, C++, Fortran, and Python; built and run as
  ctest tests when `BUILD_EXAMPLES=ON` (the default).
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
- `wigner_warmup()` public function that pre-allocates the calling
  thread's cached scratch up to the default-build absolute maximum
  (`j1+j2+j3 ≤ 20019` for 3j/6j/CG/Racah W/Gaunt; equal-`j` ≤ 5004 for
  9j and Fano X), so latency-sensitive callers can pay the lazy-init
  cost up front. Companion `wigner_thread_local_scratch_available()`
  returns 1 when the per-thread cache is active and 0 when each call
  allocates fresh (the no-TLS fallback).
- Per-thread caching of bigint scratch and of the prime-factorial
  decomposition table, gated on a runtime probe of TLS keyword
  availability (`__thread`, `_Thread_local`, MSVC `__declspec(thread)`)
  so the same source compiles where TLS is missing; a no-TLS fallback
  cell in CI verifies the slow path.
- Karatsuba multiplication in the in-tree schoolbook bigint above a
  measured 32-limb crossover; the schoolbook + Karatsuba combination
  is what `BUILD_FLINT=OFF` ships and is what the in-tree benchmark
  drivers measure against.
- Hardware `divq` 128/64-bit division on x86-64 via inline assembly,
  with Knuth's Algorithm D as the portable fallback for arbitrary
  64-bit divisors on aarch64, ppc64le, MSVC, and i686. A
  `-DBIGINT_NO_DIVQ` smoke-test cell on x86-64 exercises the
  Algorithm-D path in addition to its native execution on every other
  architecture.
- Hensel exact division for the 6j and 9j Pass-2 ratio recurrence,
  replacing per-limb 128/64 division with a multiplication by the
  modular inverse of the divisor mod `2^64`.
- Möller--Granlund "improved division by invariant integers" (IEEE TC
  60(2):165–175, 2011) on backends that route through the Algorithm-D
  fallback (no hardware `divq`); replaces the per-limb long division
  with a precomputed reciprocal and two multiplications.
- Inner Racah-sum optimisations: small-integer ratio recurrence in
  3j/6j/9j and Gaunt Pass 2 (replaces per-term `O(π(j_max))` prime-
  power expansion with one batched multiply + one batched divide per
  term), `uint64`-batched prime-power accumulator in
  `pfrac_lcm_scaled_product`, π(N) lookup + `p=2` shift fast path,
  trial division of `pfrac_mul_int` up to `√k`, `restrict` qualifier
  on pfrac vector adds (unblocks SSE2 `vpaddd` auto-vectorisation),
  and elimination of the per-term pfrac cache in Gaunt Pass 1.
- `BUILD_LTO=ON` is now the default. Probes the toolchain via
  `CheckIPOSupported` and silently falls back if LTO is not available;
  also disabled on MSVC, where `/GL` emits IL `.obj` files that the
  `WINDOWS_EXPORT_ALL_SYMBOLS` auto-export step cannot parse.
- `BUILD_EXAMPLES=ON` (the default) builds and runs `examples/c`,
  `examples/cpp`, `examples/fortran`, and `examples/python` as ctest
  tests on every push, so a binding-side regression that the rest of
  the test suite happens not to cover gets caught.
- `BUILD_COVERAGE=ON` builds with `--coverage -O0`, runs the full
  ctest + pytest suite, and uploads the merged `coverage.info` to
  Codecov.io. A new `Ubuntu / GCC / coverage` GHA cell drives the
  build. Forked `_exit` paths in `tests/test_oom.c` flush the gcov
  buffers via weakly-referenced `__gcov_dump` /
  `__llvm_profile_write_file` so children that complete the workload
  contribute coverage data.
- `BUILD_PYTHON=ON` Python build path: `Python3_add_library` produces
  `_wigner.so` and links it against the shared `libwignernj` rather
  than re-compiling library sources. The Fedora packaging policy
  (no vendored copies inside the Python package) is satisfied. The
  `pip install -e .` self-contained path is unchanged.
- CircleCI configuration with three architecture / libc cells absent
  from GitHub Actions or metered against its small free quota:
  aarch64 Linux GCC (`arm.medium` machine executor), musl libc x86-64
  (Alpine `docker:` image), and i686 (`i386/debian:stable docker:`
  image). The aarch64 cell is the production native run for the
  Möller--Granlund + Algorithm-D path; the musl cell exercises a non-
  glibc TLS implementation; the i686 cell exercises `bigint_arith.h`'s
  pure-C99 fallback natively (32-bit toolchains have no `__uint128_t`)
  rather than via `BIGINT_FORCE_PORTABLE` on x86-64.
- arm64 Windows MSVC CI cell on a weekly schedule + tag pushes only,
  to catch arm64-MSVC codegen regressions before a release without
  burning the smallest of GHA's free arm64 budgets.

### Changed
- Default CMake build type is now `Release` when `CMAKE_BUILD_TYPE` is
  not specified. Previously the code dropped to its `add_compile_options`
  default (no optimisation flags); users who didn't pass
  `-DCMAKE_BUILD_TYPE` got an unintended Debug-shaped build that hid
  the library's performance.
- `MAX_FACTORIAL_ARG` is no longer a hand-maintained `#define` in
  `src/primes.h`; the regenerator (`tools/gen_prime_table.py`) now
  computes it from the prime-table sieve `LIMIT` and emits a generated
  `src/prime_table_macros.h` consumed by the source. Contributors who
  need a higher `j` ceiling regenerate the prime table with a single
  larger `LIMIT` and the symbol-arg ceiling tracks automatically.
- New project policy section in `CLAUDE.md` codifies the SemVer
  convention, the clean-room implementation requirement, the no-intra-
  call-threading constraint, the `MAX_FACTORIAL_ARG` extensibility
  requirement, and the bench-driven optimisation methodology;
  `docs/optimization_notes.md` records the optimisations that landed
  and the three (Toom-3, compressed `pfrac_t.exp` rows, small-bigint
  Pass-2 specialisation) that were paired-benched and reverted as
  net-negative.

### Fixed
- Windows MSVC link errors: `WINDOWS_EXPORT_ALL_SYMBOLS=ON` for
  auto-generated import libraries; explicit `__declspec(dllexport)`/
  `dllimport` decoration for the prime-table data symbols
  (`g_nprimes`, `g_primes`, `g_pi_table`); co-location of test
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
- macOS Clang link error in non-coverage builds: the weak external
  references to `__gcov_dump` and `__llvm_profile_write_file` in
  `tests/test_oom.c` are now gated on `WIGNERNJ_COVERAGE`, which
  CMake's `BUILD_COVERAGE=ON` defines. ELF leaves an unresolved weak
  reference null at link time, but Mach-O demands the symbol be
  defined, so non-coverage builds on macOS would fail the link step.

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

[Unreleased]: https://github.com/susilehtola/libwignernj/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/susilehtola/libwignernj/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/susilehtola/libwignernj/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/susilehtola/libwignernj/releases/tag/v0.1.0
