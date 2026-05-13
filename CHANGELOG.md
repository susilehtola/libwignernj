# Changelog

All notable changes to **libwignernj** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **PEP 376 / PEP 427 `.dist-info` for the CMake-installed Python
  extension.**  When the library is built with `-DBUILD_PYTHON=ON`,
  the install step (`cmake --install`) now writes a
  `wignernj-<version>.dist-info/` directory next to the package
  containing METADATA, WHEEL, INSTALLER, top_level.txt, and RECORD.
  Result: the CMake-installed package is fully recognised by
  `pip list`, `pip uninstall wignernj`, and
  `importlib.metadata.version("wignernj")`, the same way a
  `pip install wignernj` wheel would be.  Metadata fields
  (description, license, authors, classifiers, URLs) come from
  `pyproject.toml` so the wheel and CMake install paths share a
  single source of truth.  Implemented in
  `tools/install_python_metadata.py`, invoked via `install(CODE ...)`
  so RECORD's sha256/size pairs are computed against the
  actually-installed files (including the Python-ABI-specific
  `_wignernj.<EXT_SUFFIX>.so` and any DESTDIR staging).
- **CI guards against build-system drift.**  The existing
  `check-source-lists` job is extended (and renamed to
  `check-build-consistency`) with a new
  `tools/check_versions.py` step that verifies the four
  independently-maintained version strings (`CMakeLists.txt`,
  `pyproject.toml`, `setup.py`, `wignernj/__init__.py`) all agree.
  A new step in the coverage job stages a `cmake --install` into a
  temporary prefix and runs `tools/verify_python_metadata.py` to
  confirm the produced `.dist-info` is well-formed (RECORD entries
  hash-match the installed files; required PEP 566 fields present).
- **Real ↔ complex spherical-harmonic basis-overlap matrix as a
  public API.**  New `wignernj_real_ylm_in_complex_ylm` family fills
  the `(2l+1) × (2l+1)` unitary matrix `C` whose entries are
  `C[m_r, m_c] = <Y_l^{m_c} | S_{l, m_r}>`, equivalently the
  basis-vector relation `S_{l,m_r} = sum_{m_c} C[m_r, m_c] Y_l^{m_c}`
  under the same real-Y convention used internally by `gaunt_real`.
  Layout is column-major (Fortran / BLAS / LAPACK) with leading
  dimension `2l+1`.  A typedef shim in `wignernj.h` exposes
  `wignernj_c{float,double,ldouble}_t` (and `wignernj_cfloat128_t`
  in `wignernj_quadmath.h`) that maps to `T _Complex` on
  gcc/clang/Apple-Clang/Intel-icx, `_Tcomplex` on MSVC C, and a
  layout-compatible `struct {T _pair[2];}` in C++; all three
  representations have identical (re, im)-interleaved memory layout
  per C99 §6.2.5/13 and C++11 [complex.numbers]/4, so callers
  holding `T _Complex *` or `std::complex<T> *` (via the C++
  wrapper) pass them without a cast.  Variants:
  `wignernj_real_ylm_in_complex_ylm_f` / `_` / `_l` / `_q` (gated on
  `BUILD_QUADMATH`) in `wignernj.h` and `wignernj_quadmath.h`;
  `wignernj_real_ylm_in_complex_ylm_mpfr` (two parallel `mpfr_t`
  arrays for real and imaginary parts) in `wignernj_mpfr.h`.
  Header-only C++ overloads
  `wignernj::real_ylm_in_complex_ylm<T>()` (in-place fill and
  convenience `std::vector<std::complex<T>>` return) in
  `wignernj.hpp`.  Fortran `bind(c)` interfaces plus a typed
  `wreal_ylm_in_complex_ylm(l, c_out)` wrapper (and
  `wreal_ylm_in_complex_ylmq` under `WIGNERNJ_HAVE_QUADMATH`) in
  `wignernj_f90.F90`.  Python binding
  `wignernj.real_ylm_in_complex_ylm(l, precision='double')`
  returning a list-of-lists of `complex`.
- Unit test `tests/test_real_ylm_in_complex_ylm.c` verifying explicit l=0/1/2
  entries, unitarity for l = 0..10, and bit-equivalent reconstruction
  of `gaunt_real` by sandwiching three rows of `C` against complex
  Gaunts.
- Parallel `real_basis_lz` examples in all four language bindings
  (`examples/c/real_basis_lz.c`, `examples/cpp/real_basis_lz.cpp`,
  `examples/fortran/real_basis_lz.f90`,
  `examples/python/real_basis_lz.py`) building the orbital
  angular-momentum operator `l_z` in the real-Y basis from its
  diagonal complex-basis form via the similarity transform
  `O_real = conj(C) @ O_complex @ transpose(C)` (consequence of the
  basis-vector form `S = C Y` returned by libwignernj).  All four
  examples are registered as ctest targets and reproduce the
  textbook l = 1 result `l_z = ((0,0,+i),(0,0,0),(−i,0,0))` with
  Hermiticity residual 0.
- `wignernj_real_ylm_in_complex_ylm` is also exercised in each binding's
  `all_symbols.*` example so a downstream user of any one language
  sees the function alongside every other public symbol.

### Fixed
- **Python-side version strings synced with the project version.**
  `wignernj/__init__.py`'s `__version__` and `setup.py`'s `version=`
  argument had been lagging at `0.4.2` since the 0.5.0 tag, even
  though `CMakeLists.txt` and `pyproject.toml` advanced to `0.5.0`.
  Both are now at `0.5.0`; the new `check_versions.py` CI step
  prevents a recurrence.

## [0.5.0] – 2026-05-07

### Changed
- **Public warmup API consolidated.**  The two warmup entry points
  `wignernj_warmup(void)` (sized to the absolute prime-table ceiling)
  and `wignernj_warmup_factorial_cache(int N_max)` (factorial cache
  only) are replaced by a single `wignernj_warmup_to(int N_max)`
  that grows both per-thread caches---the Racah-pipeline scratch and
  the factorial-decomposition cache---to fit factorial arguments up
  to `N_max`.  Pass `0` to size to the absolute prime-table ceiling
  (equivalent to the old `wignernj_warmup()`); pass the result of a
  per-symbol `wigner*_max_factorial(...)` helper to size precisely
  for a workload.  Combined with the existing helpers, callers no
  longer have to issue both a `wignernj_warmup()` (for scratch) and
  a `wignernj_warmup_factorial_cache(N)` (for the factorial cache)
  to pre-allocate ahead of a hot loop.  This is a breaking
  source-level change; rename `wignernj_warmup()` to
  `wignernj_warmup_to(0)` and `wignernj_warmup_factorial_cache(N)`
  to `wignernj_warmup_to(N)` in caller code.  No wrappers are
  kept---the library is still pre-1.0, so a small public-API rename
  is preferable to deprecated shims.

### Removed
- **`wignernj_thread_local_scratch_available()` demoted to internal.**
  The runtime introspection helper had a single in-tree caller
  (`tests/test_warmup.c`) and no production use case; the same
  answer is available at compile time as the `WIGNERNJ_HAVE_TLS`
  macro from the internal `src/wignernj_tls.h`.  The prototype
  moved to `src/scratch.h`; the test now includes the internal
  header directly.  Not exposed in the C++ or Fortran bindings.

### Added
- **C++ and Fortran cache-control bindings.**  The C++ wrapper
  (`include/wignernj.hpp`) gains
  `wignernj::warmup_to`, `wignernj::thread_cleanup`,
  `wignernj::max_factorial_arg`, and the eight per-symbol
  `wignernj::max_factorial_*` helpers.  The Fortran 90 module
  (`src/fortran/wignernj_f90.F90`) gains the same set as `bind(c)`
  interfaces.  Both bindings now have full access to the per-thread
  cache-control surface that was previously C-only.  The Python
  binding deliberately does not expose these (Python callers do not
  typically pre-warm hot loops at the C extension level).
- **Closed-form fast path for the all-$m$-zero 3j symbol.**
  `wigner3j_exact` now detects the $(m_1,m_2,m_3) = (0,0,0)$ case
  and short-circuits the Racah single sum with the closed form
  $\binom{j_1\,j_2\,j_3}{0\,0\,0} = (-1)^g\,\sqrt{\Delta^2}\,
  g!/[(g-j_1)!(g-j_2)!(g-j_3)!]$, $g = (j_1+j_2+j_3)/2$ (vanishing
  by parity for odd $g$).  Bit-identical to the general Racah path
  on every input and verified against an mpmath reference at
  $j = 2,4,10,50,100,200$.  Paired bench against the parent commit
  on a 12th-gen Intel Core i5-1235U: speedup of $\sim 22\times$ at
  $j = 5$, $\sim 45\times$ at $j = 200$, and $\sim 190\times$ at
  $j = 1000$ for the $(j\,j\,j; 0\,0\,0)$ family.  `tests/gen_refs.py`
  gains a deterministic $m=(0,0,0)$ block at $2j \in \{12, 20, 24,
  40, 60, 100, 200, 300, 400\}$ so the regression suite exercises
  the new path at moderate-to-large $j$.
- **Closed-form fast path for the all-$m$-zero sub-3j inside Gaunt.**
  `gaunt_3j_racah_sum` is called twice from `gaunt_exact` (and from
  `gaunt_real_exact` via the same back-end): once with
  $(m_1,m_2)=(0,0)$ for the $\Delta$-only 3j, and once with the
  actual $(m_1,m_2,m_3)$ for the $m$-dependent 3j.  The first call
  is the all-$m$-zero case, where the Racah single sum collapses to
  $S_{3j}(0,0,0) = (-1)^{g-j_1+j_2}\, g!/[j_1! j_2! j_3!
  (g-j_1)! (g-j_2)! (g-j_3)!]$ with $g=(j_1+j_2+j_3)/2$.  Gaunt does
  not call `wigner3j_exact`, so the closed-form path added above is
  not inherited automatically; this commit ports the same
  short-circuit to `gaunt_3j_racah_sum`.  Paired bench
  for `gaunt(l, +2, l, -2, l, 0)`: $1.2\times$ at $\ell=2$ growing
  to $\sim 1.6\times$ at $\ell=80$--$200$.  The speedup is bounded
  by the cost of the remaining $(m_1,m_2,m_3)$ sub-sum which
  Gaunt still evaluates via the full Racah loop.

## [0.4.2] – 2026-05-07

### Added
- PyPI project metadata in `pyproject.toml`: `readme = "README.md"`
  for the long-form description shown on the project page;
  `authors` for the sidebar; a `[project.urls]` block with
  Homepage, Repository, Documentation, Changelog, and Issues
  links; an expanded keyword list (wigner, clebsch-gordan, racah,
  fano, gaunt, spherical-harmonics, angular-momentum,
  exact-arithmetic, prime-factorization, …); and a fuller
  classifier set including `Development Status :: 5 -
  Production/Stable`, `Intended Audience :: Science/Research`, the
  per-version Python classifiers `3.9`–`3.13`, the per-OS
  classifiers, and additional `Topic :: Scientific/Engineering`
  sub-categories.  After this lands, the PyPI sidebar at
  https://pypi.org/p/wignernj will display the complete metadata
  starting with the next published release.

### Changed
- `requires-python` raised from `>=3.7` to `>=3.9` in both
  `pyproject.toml` and `setup.py`, matching the `CIBW_BUILD` set
  in the publish workflow (CPython 3.9..3.13) and the per-push
  Python CI matrix.  Python 3.7 reached end-of-life in 2023 and
  3.8 in 2024.
- Expanded the Python C extension's per-function docstrings.  Each
  of the eight public Python functions (`wigner3j`, `wigner6j`,
  `wigner9j`, `clebsch_gordan`, `racah_w`, `fano_x`, `gaunt`,
  `gaunt_real`) now has a NumPy/SciPy-style docstring covering the
  signature, a one-line summary, the underlying mathematical
  definition (where helpful), a `Parameters` block explaining
  argument types and unit conventions, a `Returns` block, a
  `Raises` block where applicable, brief notes on phase
  conventions, and runnable `Examples` that compute textbook-known
  values.  The `Examples` blocks double as doctests verified
  against analytic references at 1e-14 tolerance during the audit
  for this change.  Visible to users via `help(wignernj.foo)` or
  `wignernj.foo.__doc__`; the `repr(wignernj.foo)` output remains
  CPython's default `<built-in function foo>` (the docstring does
  not change the repr).

## [0.4.1] – 2026-05-07

### Added
- PyPI publishing workflow `.github/workflows/publish-pypi.yml`.
  On every `v*`-tag push, the workflow uses
  [cibuildwheel](https://github.com/pypa/cibuildwheel) to build
  CPython 3.9..3.13 wheels for Linux x86_64 (manylinux + musllinux),
  Linux aarch64 (native `ubuntu-24.04-arm` runner), macOS Apple
  Silicon, and Windows x86_64, plus an sdist for any platform /
  interpreter not covered by a wheel.  Each wheel is smoke-tested
  before upload.  Authentication uses PyPI Trusted Publishing (OIDC)
  via the `id-token: write` permission and the `pypi` environment;
  no API token in repository secrets.  Manual `workflow_dispatch`
  runs build the wheels but skip the publish step, so a maintainer
  can dry-run the wheel matrix without releasing.  Wheels installed
  via `pip install wignernj` are the self-contained
  `setup.py`-driven build (every `src/*.c` recompiled into
  `_wignernj.so`); the optional libquadmath / MPFR / FLINT back-ends
  remain CMake-only.
- Five new bullets in the `## Project policy` section of `CLAUDE.md`,
  codifying conventions that had been operating only as session-level
  rules until now: **code and documentation must stay consistent**,
  **run a full repository consistency audit before tagging a release**,
  **amend `CHANGELOG.md` as part of every pull request**, **every
  commit must build without warnings**, and a `Workflow` umbrella
  bullet covering amend-don't-stack fixups, CircleCI off-load over
  Cirrus, and `export VAR=val` for shell env vars in CI recipes.
  The "What this library is" opening section also gains a one-sentence
  note that `wignernj` is the uniform name across every binding (the
  surviving `wigner`-prefixed C symbols are the math-symbol functions,
  not library-namespace identifiers).

### Fixed
- Documentation consistency audit covering attribution, public-symbol
  enumeration, default-build numeric ceilings, the build-options table,
  and the warmup/cache API:
  - **Attribution.** The prime-factorization technique was attributed
    solely to Johansson & Forssén in `README.md`, `docs/reference.md`,
    `CLAUDE.md`, and the source comments of `src/wigner3j.c` and
    `src/wigner6j.c`.  Prime factorization for the angular-momentum
    coefficients was introduced earlier by Dodds & Wiechers (*Comput.
    Phys. Commun.* **4**, 268 (1972)) and refined in subsequent work;
    Johansson & Forssén's specific contribution is the multiword-integer
    Racah sum.  Every site now credits both lineages; the Dodds &
    Wiechers DOI is hyperlinked in `README.md` and `docs/reference.md`.
    The `CHANGELOG.md` `[0.1.0]` historical entry is intentionally not
    rewritten.
  - **Fano X-coefficient.** Added to the public-symbol family list in
    six places that had been missed when Fano X landed in 0.3.0: the
    "wrappers over the Wigner symbols" paragraph in `CLAUDE.md`, the
    phase-conventions paragraph in `include/wignernj.h`,
    `include/wignernj.hpp`, `src/fortran/wignernj_f90.F90`, and
    `wignernj/__init__.py`, the j-limit and performance-scaling
    comments in `src/primes.h`, the Python-section example block and
    performance scaling row in `docs/reference.md`, and the Fortran
    wrapper-list comment in `src/fortran/wignernj_f90.F90` (which had
    been missing both `wfanox` and `wgaunt_real`).
  - **Default-build ceiling.**  Stale `MAX_FACTORIAL_ARG = 20000`
    references in `src/wigner3j.c`, `src/pfrac.c`, and
    `include/wignernj.h` updated to `20020`, matching the value
    derived from `PRIME_SIEVE_LIMIT` in `src/prime_table_macros.h`.
  - **Build-options table.**  `docs/reference.md` was missing
    `BUILD_LTO`, `BUILD_EXAMPLES`, and `BUILD_COVERAGE`, all of
    which landed in 0.3.0/0.4.0.  Added them so the table now matches
    `CMakeLists.txt` and the `README.md` table.
  - **Warmup / cache API.**  `wignernj_warmup`,
    `wignernj_thread_local_scratch_available`,
    `wignernj_warmup_factorial_cache`, `wignernj_max_factorial_arg`,
    `wignernj_thread_cleanup`, and the eight per-symbol
    `wigner*_max_factorial` companions are public since 0.3.0 but
    were undocumented in `docs/reference.md`.  Added a
    "Per-thread caches and warmup" subsection to the C API chapter
    that covers each function with its semantics, memory cost, and
    no-TLS fallback behaviour.

  No code or behaviour change; this PR amends documentation only.

## [0.4.0] – 2026-05-07

### Changed
- **Breaking (all bindings): public-facing names renamed `wigner` → `wignernj`
  for naming consistency with the C library (`libwignernj`) and the Fortran
  shared library (`libwignernj_f03`).** The C symbol names (`wigner3j`,
  `wigner6j`, `wigner9j`, …, `gaunt_real`) keep their `wigner` prefix
  because that prefix denotes the mathematical object, not the library;
  every other library-namespace identifier moves:
  - **C/C++ public headers**: `<wigner.h>` → `<wignernj.h>`,
    `<wigner.hpp>` → `<wignernj.hpp>`, `<wigner_quadmath.h>` →
    `<wignernj_quadmath.h>`, `<wigner_mpfr.h>` → `<wignernj_mpfr.h>`.
    Include guards renamed in step (`WIGNER_H` → `WIGNERNJ_H`, …).
  - **C++ namespace**: `wigner::symbol3j<>`, `wigner::cg`, …, `wigner::fanox<>`
    → `wignernj::symbol3j<>`, `wignernj::cg`, …, `wignernj::fanox<>`.
  - **Fortran module**: `module wigner` → `module wignernj`. Callers must
    change `use wigner` to `use wignernj`. The shared-library name
    (`libwignernj_f03`) is unchanged.
  - **Python package**: `pip install wigner` / `import wigner` →
    `pip install wignernj` / `import wignernj`. The CPython extension
    binary moves from `_wigner.so` to `_wignernj.so`; `PyInit__wigner` →
    `PyInit__wignernj`.
  - **C library-namespace identifiers**: every `wigner_*` /
    `WIGNER_*` identifier outside the mathematical-symbol functions
    moves to `wignernj_*` / `WIGNERNJ_*`. This includes the public
    helpers `wignernj_warmup`, `wignernj_thread_local_scratch_available`,
    `wignernj_max_factorial_arg`; the public exact-arithmetic tuple
    `wignernj_exact_t` plus `wignernj_exact_init`/`reset`/`free`; the
    conversion routines `wignernj_exact_to_{float,double,long_double,
    float128,mpfr,T}`; the internal scratch allocator
    (`wignernj_scratch_t`, `wignernj_scratch_acquire`/`release`/…);
    the feature macros `WIGNERNJ_HAVE_QUADMATH` / `WIGNERNJ_HAVE_MPFR`
    (the C side now matches the Fortran side, which already used the
    `WIGNERNJ_` spelling); and the include guards on every renamed
    header. The internal header `src/wigner_exact.h` was renamed to
    `src/wignernj_exact.h` in step.
  - **Source files** (cosmetic, only the moved ones surface in build
    rules): `src/python/wignermodule.c` → `wignernjmodule.c`,
    `src/fortran/wigner_f90.F90` → `wignernj_f90.F90`,
    `tests/python/test_wigner_python.py` → `test_wignernj_python.py`,
    `tests/fortran/test_wigner_fortran.f90` → `test_wignernj_fortran.f90`,
    `tests/fortran/test_wigner_quadmath.F90` → `test_wignernj_quadmath.F90`.
    The `src/wigner3j.c`, `src/wigner6j.c`, `src/wigner9j.c`, and
    `src/wigner_exact.c` source files keep their names because they reflect
    the public C function names rather than the library namespace.

  Migration is a search-and-replace at every call site:

  ```sed
  s/\bwigner\.h\b/wignernj.h/g
  s/\bwigner\.hpp\b/wignernj.hpp/g
  s/\bwigner_quadmath\.h\b/wignernj_quadmath.h/g
  s/\bwigner_mpfr\.h\b/wignernj_mpfr.h/g
  s/\bwigner::/wignernj::/g
  s/\buse wigner\b/use wignernj/g
  s/\bimport wigner\b/import wignernj/g
  s/\bwigner\.\(wigner3j\|wigner6j\|wigner9j\|clebsch_gordan\|racah_w\|fano_x\|gaunt\|gaunt_real\)/wignernj.\1/g
  s/\bwigner_\(warmup\|thread_local_scratch_available\|max_factorial_arg\|exact_\)/wignernj_\1/g
  s/\bWIGNER_\(HAVE_\|EXACT_H\|SCRATCH_\)/WIGNERNJ_\1/g
  ```

  No code logic changes; only identifiers move.

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

[Unreleased]: https://github.com/susilehtola/libwignernj/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/susilehtola/libwignernj/compare/v0.4.2...v0.5.0
[0.4.2]: https://github.com/susilehtola/libwignernj/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/susilehtola/libwignernj/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/susilehtola/libwignernj/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/susilehtola/libwignernj/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/susilehtola/libwignernj/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/susilehtola/libwignernj/releases/tag/v0.1.0
