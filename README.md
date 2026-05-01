# libwignernj

Exact evaluation of Wigner 3j, 6j, and 9j symbols, Clebsch-Gordan coefficients,
Racah W-coefficients, and Gaunt coefficients in ANSI C.

All intermediate arithmetic is exact integer arithmetic using a prime-factorization
representation; floating-point conversion happens only at the final step.  Results
are accurate to the last bit of the chosen output precision.

Algorithm: Johansson & Forssén, SIAM J. Sci. Comput. 38(1), A376–A384, 2016.
doi:[10.1137/15M1021908](https://doi.org/10.1137/15M1021908)

Language interfaces: C (primary), C++ (header-only wrapper, links against `libwignernj`), Python (CPython extension),
Fortran 90 (iso_c_binding).

## Argument convention

All angular momentum arguments are passed as **twice their value** so that
half-integers are represented exactly as odd integers:

```
j = 3/2  →  tj = 3
m = -1/2 →  tm = -1
```

This applies to the C, C++, and Fortran interfaces.  The Python interface also
accepts plain floats (e.g. `0.5`) and `fractions.Fraction` objects.

## Building

```sh
cmake -B build && cmake --build build
ctest --test-dir build
```

CMake options (all `ON` by default except `BUILD_PYTHON`):

| Option | Default | Description |
|---|---|---|
| `BUILD_SHARED_LIBS` | `ON` | Shared library |
| `BUILD_FORTRAN` | `ON` | Fortran interface |
| `BUILD_TESTS` | `ON` | C/Fortran test suite |
| `BUILD_CXX_TESTS` | `ON` | C++ header tests |
| `BUILD_PYTHON` | `OFF` | Python extension |
| `BUILD_MPFR` | `OFF` | MPFR arbitrary-precision interface |

## C API

```c
#include "wigner.h"

/* Wigner 3j:  ( j1  j2  j3 ) */
/*             ( m1  m2  m3 ) */
double      wigner3j  (int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
float       wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);

/* Wigner 6j:  { j1 j2 j3 } */
/*             { j4 j5 j6 } */
double      wigner6j  (int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
float       wigner6j_f(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
long double wigner6j_l(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);

/* Wigner 9j (row-major order) */
double      wigner9j  (int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);

/* Clebsch-Gordan:  <j1 m1; j2 m2 | J M> */
double      clebsch_gordan  (int tj1, int tm1, int tj2, int tm2, int tJ, int tM);

/* Racah W-coefficient:  W(j1 j2 J j3; j12 j23) */
double      racah_w  (int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);

/* Gaunt coefficient:  integral Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ */
double      gaunt  (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
```

Each function is available in three precisions: `double` (no suffix), `float`
(`_f`), and `long double` (`_l`).  Functions return 0 for symbols that vanish
by selection rules; selection-rule violations are not errors.

Linking: `pkg-config --libs libwignernj` or `-lwignernj -lm`.

## MPFR API

Build with `-DBUILD_MPFR=ON` (requires libmpfr).  Include `wigner_mpfr.h` in
addition to `wigner.h`.  Set the output precision on `rop` before calling;
the rounding mode is the last argument.

```c
#include "wigner.h"
#include "wigner_mpfr.h"

mpfr_t v;
mpfr_init2(v, 256);   /* 256-bit precision */

wigner3j_mpfr(v, 4, 4, 0,  2, -2, 0,  MPFR_RNDN);
wigner6j_mpfr(v, 2, 2, 0,  2,  2, 0,  MPFR_RNDN);
wigner9j_mpfr(v, 2, 2, 2,  2, 2, 2,  2, 2, 2,  MPFR_RNDN);
clebsch_gordan_mpfr(v, 2, 2,  2, -2,  4, 0,  MPFR_RNDN);
racah_w_mpfr(v, 2, 2, 4,  2,  4, 4,  MPFR_RNDN);
gaunt_mpfr(v, 4, 2,  4, -2,  4, 0,  MPFR_RNDN);

mpfr_clear(v);
```

Link with `-lwignernj -lmpfr -lm`.

## C++ API

The header-only wrapper `wigner.hpp` provides a template interface that accepts
either `2*j` integers or real-valued doubles (half-integers).  The wrapper has
no separate translation unit, but every function forwards to a C symbol in
`libwignernj`, so you still need to link the C library (`-lwignernj -lm`):

```cpp
#include "wigner.hpp"

// Integer 2*j form
double v = wigner::symbol3j<double>(2, 2, 0,  0, 0, 0);
float  f = wigner::symbol6j<float> (2, 2, 2,  2, 2, 2);

// Real-valued form (throws std::invalid_argument if not a half-integer)
double v = wigner::symbol3j(1.0, 1.0, 0.0,  0.0, 0.0, 0.0);
double c = wigner::cg(0.5, 0.5, 0.5, -0.5, 1.0, 0.0);

// Available functions: symbol3j, symbol6j, symbol9j, cg, racahw, gaunt
```

Link with `-lwignernj -lm` (and `-lmpfr` if `BUILD_MPFR=ON`).

## Python API

```sh
pip install -e .                    # build and install in-place
```

```python
import wigner

wigner.wigner3j(1, 1, 0, 0, 0, 0)             # integer 2*j form
wigner.wigner3j(0.5, 0.5, 1, 0.5, -0.5, 0)   # real half-integer form
wigner.wigner6j(1, 1, 2, 1, 1, 2)
wigner.wigner9j(1, 1, 2, 1, 1, 2, 2, 2, 4)
wigner.clebsch_gordan(1, 1, 1, -1, 2, 0)
wigner.racah_w(1, 1, 2, 1, 2, 2)
wigner.gaunt(2, 1, 2, -1, 2, 0, precision='longdouble')
```

The optional `precision=` keyword selects `'float'`, `'double'` (default), or
`'longdouble'`.  Arguments may be integers, floats, or `fractions.Fraction`.

## Fortran API

The `wigner` module provides real-valued wrappers `w3j`, `w6j`, `w9j`, `wcg`,
`wracahw`, `wgaunt` that accept double-precision real arguments:

```fortran
use wigner
real(8) :: v
v = w3j(1.0d0, 1.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0)
v = w6j(1.0d0, 1.0d0, 2.0d0,  1.0d0, 1.0d0, 2.0d0)
v = wcg(0.5d0, 0.5d0, 0.5d0, -0.5d0, 1.0d0, 0.0d0)
```

Raw `bind(c)` interfaces using `2*j` integers are also available for all
three precisions (the `long double` variants require gfortran or Cray Fortran).

Link with `-lwignernj_fortran -lwignernj -lm`.

## Limits

The prime table covers factorials up to 20000!, which translates to:

- 3j / 6j / CG / Racah W / Gaunt: j1+j2+j3 ≤ 19999 (equal-j: **j ≤ 6666**)
- 9j: equal-j **j ≤ 4999** (k-dependent triangle denominators reach (4j+1)!)

Exceeding these limits prints a diagnostic to stderr and aborts.  The 9j is also
O(j²) in computation time; evaluations with j > a few hundred can be slow.
See [docs/reference.md](docs/reference.md#limitations) for details.

## Documentation

Full API reference with mathematical definitions, selection rules, and
per-language examples: [docs/reference.md](docs/reference.md).

## License

BSD 3-Clause — see [LICENSE](LICENSE).
