# libwignernj — Reference Manual

## Overview

**libwignernj** computes Wigner 3j, 6j, and 9j symbols, Clebsch-Gordan
coefficients, Racah W-coefficients, and Gaunt coefficients.  All intermediate
arithmetic is exact integer arithmetic using a prime-factorization
representation; floating-point conversion happens only at the final step.
Results are accurate to the last bit of the chosen output precision.

Algorithm: Johansson & Forssén, *SIAM J. Sci. Comput.* **38**(1), A376–A384,
2016. [doi:10.1137/15M1021908](https://doi.org/10.1137/15M1021908)

---

## Building

```sh
# Minimal build (shared library + tests)
cmake -B build && cmake --build build
ctest --test-dir build

# All features
cmake -B build \
    -DBUILD_FORTRAN=ON \
    -DBUILD_CXX_TESTS=ON \
    -DBUILD_QUADMATH=ON \
    -DBUILD_MPFR=ON \
    -DBUILD_PYTHON=ON
cmake --build build

# Python package (self-contained, no prior CMake build needed)
pip install -e .
python -m pytest tests/python/
```

### CMake options

| Option | Default | Description |
|---|---|---|
| `BUILD_SHARED_LIBS` | `ON` | Shared library (`.so` / `.dylib` / `.dll`) |
| `BUILD_FORTRAN` | `ON` | Fortran 90 interface (`libwignernj_f03`) |
| `BUILD_TESTS` | `ON` | C and Fortran test suite |
| `BUILD_CXX_TESTS` | `ON` | C++ header tests |
| `BUILD_QUADMATH` | `OFF` | libquadmath / IEEE 754 binary128 (`__float128`) interface (requires GCC, Clang, or Intel ICC/ICX on Linux/macOS) |
| `BUILD_MPFR` | `OFF` | MPFR arbitrary-precision interface (requires libmpfr) |
| `BUILD_FLINT` | `OFF` | FLINT/GMP/MPFR bigint backend (sub-quadratic multiplication; replaces the in-tree schoolbook bigint, requires libflint, libgmp, and libmpfr) |
| `BUILD_PYTHON` | `OFF` | CPython extension module |

### Linking

```sh
# pkg-config
pkg-config --cflags --libs libwignernj

# Manual
-lwignernj -lm

# With libquadmath interface
-lwignernj -lquadmath -lm

# With MPFR interface
-lwignernj -lmpfr -lm

# With Fortran interface
-lwignernj_f03 -lwignernj -lm
```

---

## Argument convention

All angular momentum arguments are passed as **twice their value** so that
half-integers are represented as odd integers and all arguments are plain C
`int`:

| Physical value | Argument |
|---|---|
| j = 0 | tj = 0 |
| j = 1/2 | tj = 1 |
| j = 1 | tj = 2 |
| j = 3/2 | tj = 3 |
| m = −1/2 | tm = −1 |

This convention applies to the C, C++, and Fortran `bind(c)` interfaces.
The Fortran convenience wrappers and the Python interface also accept
real-valued half-integers directly.

---

## Coefficients

### Wigner 3j symbol

```
( j1  j2  j3 )
( m1  m2  m3 )
```

**Selection rules** (returns 0 if any is violated):
- m1 + m2 + m3 = 0
- |mi| ≤ ji for each i
- Triangle inequality: |j1 − j2| ≤ j3 ≤ j1 + j2
- j1 + j2 + j3 is a non-negative integer (sum is integer)

**Special values:**

```
( j  j  0 )  =  (−1)^(j−m) / sqrt(2j+1)
( m −m  0 )
```

---

### Wigner 6j symbol

```
{ j1  j2  j3 }
{ j4  j5  j6 }
```

**Selection rules** (returns 0 if any is violated):

Each of the four triples must satisfy the triangle inequality and have
integer sum:

- (j1, j2, j3)
- (j1, j5, j6)
- (j4, j2, j6)
- (j4, j5, j3)

**Special value:**

```
{ j1  j2   0  }  =  (−1)^(j1+j2) / sqrt[(2j1+1)(2j2+1)]   (if j1, j2 satisfy triangle)
{ j2  j1  j12 }     (for all allowed j12)
```

---

### Wigner 9j symbol

```
{ j11  j12  j13 }
{ j21  j22  j23 }   (arguments in row-major order)
{ j31  j32  j33 }
```

**Selection rules** (returns 0 if any is violated):

Each row and each column must satisfy the triangle inequality:

- Rows: (j11, j12, j13), (j21, j22, j23), (j31, j32, j33)
- Columns: (j11, j21, j31), (j12, j22, j32), (j13, j23, j33)

---

### Clebsch-Gordan coefficient

```
<j1 m1; j2 m2 | J M>
```

Related to the 3j symbol by:

```
<j1 m1; j2 m2 | J M> = (−1)^(j1−j2+M) sqrt(2J+1) * ( j1  j2   J )
                                                       ( m1  m2  −M )
```

**Selection rules:**
- m1 + m2 = M
- |mi| ≤ ji, |M| ≤ J
- Triangle inequality: |j1 − j2| ≤ J ≤ j1 + j2
- j1 + j2 + J is a non-negative integer

---

### Racah W-coefficient

```
W(j1 j2 J j3; j12 j23)
```

Related to the 6j symbol by:

```
W(j1 j2 J j3; j12 j23) = (−1)^(j1+j2+J+j3) * { j1   j2   j12 }
                                                  { j3   J    j23 }
```

The same selection rules as the 6j symbol apply.

---

### Fano X-coefficient

```
X(j1 j2 j12; j3 j4 j34; j13 j24 J)
```

A normalisation variant of the 9j symbol used in the analysis of
polarisation correlations (Fano 1953; Edmonds 1957 §6.4):

```
X(j1 j2 j12; j3 j4 j34; j13 j24 J)
  = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)]
    * { j1   j2   j12 }
      { j3   j4   j34 }
      { j13  j24  J   }
```

The same selection rules as the 9j symbol apply.

---

### Gaunt coefficient

```
G(l1,m1,l2,m2,l3,m3) = ∫ Y_{l1}^{m1}(Ω) Y_{l2}^{m2}(Ω) Y_{l3}^{m3}(Ω) dΩ
```

In terms of 3j symbols:

```
G = sqrt[(2l1+1)(2l2+1)(2l3+1) / (4π)]
    * ( l1  l2  l3 ) * ( l1  l2  l3 )
      (  0   0   0 )   ( m1  m2  m3 )
```

**Selection rules:**
- m1 + m2 + m3 = 0
- |mi| ≤ li for each i
- l1 + l2 + l3 is even (parity)
- Triangle inequality: |l1 − l2| ≤ l3 ≤ l1 + l2
- All li must be non-negative integers (tli = 2*li is always even)

---

## C API

Include `wigner.h`.  Link with `-lwignernj -lm`.

All functions return 0 (exactly) for selection-rule violations.

### Wigner 3j

```c
float       wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
double      wigner3j  (int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
```

### Wigner 6j

```c
float       wigner6j_f(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
double      wigner6j  (int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
long double wigner6j_l(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
```

### Wigner 9j

Arguments are in row-major order: first row, then second row, then third row.

```c
float       wigner9j_f(int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);
double      wigner9j  (int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);
long double wigner9j_l(int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);
```

### Clebsch-Gordan

```c
float       clebsch_gordan_f(int tj1, int tm1, int tj2, int tm2, int tJ, int tM);
double      clebsch_gordan  (int tj1, int tm1, int tj2, int tm2, int tJ, int tM);
long double clebsch_gordan_l(int tj1, int tm1, int tj2, int tm2, int tJ, int tM);
```

### Racah W

```c
float       racah_w_f(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
double      racah_w  (int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
long double racah_w_l(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
```

### Fano X

```c
float       fano_x_f(int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);
double      fano_x  (int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);
long double fano_x_l(int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);
```

### Gaunt

```c
float       gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
double      gaunt  (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
```

### Example (C)

```c
#include "wigner.h"
#include <stdio.h>

int main(void)
{
    /* (1  1  0)  = -1/sqrt(3) */
    /* (0  0  0)               */
    double v = wigner3j(2, 2, 0,  0, 0, 0);
    printf("%.15g\n", v);   /* -0.577350269189626 */

    /* Clebsch-Gordan <1/2, 1/2; 1/2, -1/2 | 1, 0> = 1/sqrt(2) */
    double c = clebsch_gordan(1, 1,  1, -1,  2, 0);
    printf("%.15g\n", c);   /* 0.707106781186548 */

    /* Gaunt G(1,0,1,0,2,0) */
    double g = gaunt(2, 0,  2, 0,  4, 0);
    printf("%.15g\n", g);

    return 0;
}
```

---

## C++ API

Header-only; include `wigner.hpp`.  No additional linking required.
Requires C++11.

### Integer-argument form

```cpp
#include "wigner.hpp"

// Default precision is double
double v3 = wigner::symbol3j<double>(2, 2, 0,  0, 0, 0);
float  v6 = wigner::symbol6j<float> (2, 2, 2,  2, 2, 2);
long double v9 = wigner::symbol9j<long double>(
                     2,2,2, 2,2,2, 2,2,2);

double cg = wigner::cg<double>    (1, 1,  1, -1,  2, 0);
double rw = wigner::racahw<double>(2, 2,  4,  2,  4, 4);
double ga = wigner::gaunt<double> (2, 0,  2,  0,  4, 0);
```

### Real-valued form

Accepts integer or half-integer `double` arguments. Throws
`std::invalid_argument` if an argument is not a half-integer.

```cpp
// j = 1/2 written as 0.5
double v = wigner::symbol3j(0.5, 0.5, 1.0,  0.5, -0.5, 0.0);

// Explicit precision template
float  f = wigner::symbol3j<float>(1.0, 1.0, 0.0,  0.0, 0.0, 0.0);
```

### Available functions

| Function | Arguments |
|---|---|
| `wigner::symbol3j<T>(...)` | `(tj1,tj2,tj3, tm1,tm2,tm3)` |
| `wigner::symbol6j<T>(...)` | `(tj1,tj2,tj3, tj4,tj5,tj6)` |
| `wigner::symbol9j<T>(...)` | `(tj11..tj33)` row-major |
| `wigner::cg<T>(...)` | `(tj1,tm1, tj2,tm2, tJ,tM)` |
| `wigner::racahw<T>(...)` | `(tj1,tj2,tJ, tj3,tj12,tj23)` |
| `wigner::gaunt<T>(...)` | `(tl1,tm1, tl2,tm2, tl3,tm3)` |

`T` is `float`, `double`, or `long double`.

---

## libquadmath API

Build with `-DBUILD_QUADMATH=ON` (requires GCC, Clang, or Intel
ICC/ICX with `__float128` support; not available on Apple Clang or
MSVC).  Include `wigner_quadmath.h` in addition to `wigner.h`.  Link
with `-lwignernj -lquadmath -lm`.

```c
#include "wigner.h"
#include "wigner_quadmath.h"

__float128 v = wigner6j_q(4, 4, 4, 4, 4, 4);
```

### Function signatures

```c
__float128 wigner3j_q(int tj1, int tj2, int tj3,
                      int tm1, int tm2, int tm3);
__float128 wigner6j_q(int tj1, int tj2, int tj3,
                      int tj4, int tj5, int tj6);
__float128 wigner9j_q(int tj11, int tj12, int tj13,
                      int tj21, int tj22, int tj23,
                      int tj31, int tj32, int tj33);
__float128 clebsch_gordan_q(int tj1, int tm1, int tj2, int tm2,
                            int tJ,  int tM);
__float128 racah_w_q(int tj1, int tj2, int tJ,
                     int tj3, int tj12, int tj23);
__float128 fano_x_q(int tj1, int tj2, int tj12,
                    int tj3, int tj4, int tj34,
                    int tj13, int tj24, int tJ);
__float128 gaunt_q     (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
__float128 gaunt_real_q(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
```

The Fortran `wigner` module exposes the same routines (and real-valued
convenience wrappers `w3jq`, `w6jq`, `w9jq`, `wcgq`, `wracahwq`,
`wgauntq`, `wgaunt_realq`) returning `real(real128)`.

Results are accurate to within 2 ulp at quad precision; the conversion
path Horner-evaluates the top three 64-bit words of the bigint in
`__float128` arithmetic, feeding 192 input bits into a 113-bit
mantissa.

---

## MPFR API

Build with `-DBUILD_MPFR=ON` (requires libmpfr ≥ 4).
Include `wigner_mpfr.h` in addition to `wigner.h`.
Link with `-lwignernj -lmpfr -lm`.

The precision of the result is determined by the precision set on `rop`
before calling.  The rounding mode `rnd` is applied to each elementary
floating-point operation.

```c
#include "wigner.h"
#include "wigner_mpfr.h"

mpfr_t v;
mpfr_init2(v, 512);          /* 512-bit precision */

wigner3j_mpfr(v,
    4, 4, 0,                 /* tj1, tj2, tj3 */
    2, -2, 0,                /* tm1, tm2, tm3 */
    MPFR_RNDN);
/* v now holds (-1)^(j-m)/sqrt(2j+1) at 512-bit precision */

mpfr_clear(v);
```

### Function signatures

```c
void wigner3j_mpfr(mpfr_t rop,
                   int tj1, int tj2, int tj3,
                   int tm1, int tm2, int tm3,
                   mpfr_rnd_t rnd);

void wigner6j_mpfr(mpfr_t rop,
                   int tj1, int tj2, int tj3,
                   int tj4, int tj5, int tj6,
                   mpfr_rnd_t rnd);

void wigner9j_mpfr(mpfr_t rop,
                   int tj11, int tj12, int tj13,
                   int tj21, int tj22, int tj23,
                   int tj31, int tj32, int tj33,
                   mpfr_rnd_t rnd);

void clebsch_gordan_mpfr(mpfr_t rop,
                         int tj1, int tm1,
                         int tj2, int tm2,
                         int tJ,  int tM,
                         mpfr_rnd_t rnd);

void racah_w_mpfr(mpfr_t rop,
                  int tj1, int tj2, int tJ, int tj3,
                  int tj12, int tj23,
                  mpfr_rnd_t rnd);

void fano_x_mpfr(mpfr_t rop,
                 int tj1, int tj2, int tj12,
                 int tj3, int tj4, int tj34,
                 int tj13, int tj24, int tJ,
                 mpfr_rnd_t rnd);

void gaunt_mpfr(mpfr_t rop,
                int tl1, int tm1,
                int tl2, int tm2,
                int tl3, int tm3,
                mpfr_rnd_t rnd);
```

All functions set `rop` to zero for selection-rule violations.

---

## Python API

```sh
pip install -e .           # build and install in-place
python -m pytest tests/python/
```

```python
import wigner

# Integer 2*j form
wigner.wigner3j(2, 2, 0,  0, 0, 0)          # -1/sqrt(3) ≈ -0.5774
wigner.wigner6j(2, 2, 2,  2, 2, 2)
wigner.wigner9j(2, 2, 2,  2, 2, 2,  2, 2, 2)

# Real half-integer form
wigner.wigner3j(0.5, 0.5, 1.0,  0.5, -0.5, 0.0)

# fractions.Fraction form
from fractions import Fraction
wigner.wigner3j(Fraction(1,2), Fraction(1,2), 1,
                Fraction(1,2), Fraction(-1,2), 0)

# Clebsch-Gordan, Racah W, Gaunt
wigner.clebsch_gordan(1, 1,  1, -1,  2, 0)
wigner.racah_w(2, 2,  4,  2,  4, 4)
wigner.gaunt(4, 2,  4, -2,  4, 0)
```

The optional `precision=` keyword selects the output type:

| Value | Type | Digits |
|---|---|---|
| `'float'` | 32-bit float | ~7 |
| `'double'` | 64-bit double (default) | ~15 |
| `'longdouble'` | 80/128-bit extended | ~18–33 |

```python
wigner.gaunt(2, 1,  2, -1,  2, 0, precision='longdouble')
```

Returns `0.0` for selection-rule violations (not an error).
Raises `ValueError` for non-half-integer arguments.

---

## Fortran API

Build with `-DBUILD_FORTRAN=ON`.  Link with `-lwignernj_f03 -lwignernj -lm`.

The `wigner` module provides two layers:

### Convenience wrappers (recommended)

Accept `real(8)` j/m arguments directly:

```fortran
use wigner
real(8) :: v

! Wigner 3j
v = w3j(1.0d0, 1.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0)

! Wigner 6j
v = w6j(1.0d0, 1.0d0, 2.0d0,  1.0d0, 1.0d0, 2.0d0)

! Wigner 9j
v = w9j(1.0d0, 1.0d0, 2.0d0, &
        1.0d0, 1.0d0, 2.0d0, &
        2.0d0, 2.0d0, 4.0d0)

! Clebsch-Gordan <1/2, 1/2; 1/2, -1/2 | 1, 0>
v = wcg(0.5d0, 0.5d0,  0.5d0, -0.5d0,  1.0d0, 0.0d0)

! Racah W
v = wracahw(1.0d0, 1.0d0, 2.0d0, 1.0d0,  2.0d0, 2.0d0)

! Gaunt
v = wgaunt(1.0d0, 0.0d0,  1.0d0, 0.0d0,  2.0d0, 0.0d0)
```

On invalid input (non-half-integer argument), the wrappers write a
diagnostic to stderr and return 0.

### Raw `bind(c)` interfaces

Use the `2*j` integer convention directly.  Available in `float`,
`double`, and (where supported by the compiler) `long double`:

```fortran
use wigner
real(8) :: v
v = wigner3j(2, 2, 0,  0, 0, 0)   ! (1 1 0; 0 0 0)
v = wigner6j(2, 2, 2,  2, 2, 2)
```

The `long double` interfaces (`wigner3j_l`, `wigner6j_l`, etc.) are
conditionally compiled and only available with gfortran or Cray Fortran
(`c_long_double > 0`).

---

## Accuracy and precision

### Fixed-width outputs

The `double` functions are accurate to within 1 ULP of the true value for
all inputs tested.  The exact-integer intermediate result is converted to
floating-point with round-to-nearest-even using the correct number of
mantissa bits (`DBL_MANT_DIG = 53`).

The `long double` functions use `LDBL_MANT_DIG` bits:
- 64 bits on x86-64 with glibc (80-bit extended precision)
- 53 bits on MSVC and most ARM targets (same as `double`)
- 113 bits on AArch64/POWER with glibc (128-bit quad precision)

The `float` functions use 24 bits (`FLT_MANT_DIG = 24`).

### MPFR output

The result is computed exactly in integer arithmetic up to the bigint
multiply and LCM scaling; the final `sqrt` and division are performed in
MPFR at the caller's requested precision.  The total rounding error is at
most a few ULP at the requested precision.

### Selection-rule zeros

Functions return exactly 0 (not a small float) when selection rules are
violated.  No computation is performed for zero results.

---

## Limitations

### Angular momentum range (correctness)

All intermediate factorials are factored by trial division against a prime
table sieved up to `PRIME_SIEVE_LIMIT = 20011`.  A factorial argument `n`
that has a prime factor larger than 20011 would be silently mis-factored.
The first such `n` is 20021 (which is prime).  `MAX_FACTORIAL_ARG = 20000`
is set conservatively below this boundary.

The prime list (`g_primes`) and its inverse-lookup index (`g_prime_index`)
are hard-coded into the compiled library (~9 kB and ~40 kB respectively,
totalling ~49 kB) rather than built at run time.  This is a deliberate
design choice: it means the library has **no caller-side initialization
step**, is safely usable from concurrent threads with no initialization
race, and lets the tables live in read-only program data.  The default
sieve limit was chosen on a rule of thumb that ~50 kB is a reasonable
upper bound for compile-time constant tables, while still admitting
angular momenta well beyond what most application domains require.

The binding factorial argument in each coefficient is the triangle-coefficient
denominator `(j1+j2+j3+1)!`.  The resulting default-build limits on the sum of
angular momenta are:

| Symbol | Binding factorial | Limit on sum | Equal-j limit |
|---|---|---|---|
| 3j, 6j, CG, Racah W, complex Gaunt, real Gaunt | `(j1+j2+j3+1)!` | j1+j2+j3 ≤ 19999 | **j ≤ 6666** |
| 9j | `(4j+1)!` (k-dependent Δ at k = k_max) | 4j ≤ 19999 | **j ≤ 4999** |

These limits apply to the sum of the three largest angular momenta appearing
in any single triangle.  Symbols with very unequal arguments can involve
larger individual `j` values as long as no triangle sum exceeds 19999.

Exceeding these limits causes the library to print a diagnostic to stderr
and call `abort()`.  The check is unconditional (not gated on `NDEBUG`).
The ceiling is **not architectural**: it is determined by the size of the
compile-time prime table (`PRIME_SIEVE_LIMIT` in `src/primes.h`) and can
be raised by regenerating the table with `tools/gen_prime_table.py` for a
larger sieve limit and rebuilding the library.

### Performance scaling

The Racah summation has O(j) terms for 3j and 6j; the 9j outer loop over the
intermediate quantum number k adds another factor of O(j), so the *number* of
sum terms is O(j) for 3j/6j and O(j²) for 9j.  Intermediate bigints grow
proportional to j/ln 2 bits (the LCM denominator grows as the primorial), so
each elementary bigint operation costs O(j).  The resulting end-to-end cost
per symbol is

- **3j, 6j, CG, Racah W, complex / real Gaunt:** O(j²)
- **9j:** O(j⁴) (each k step contains three multiplications of size-O(j²) bigints)

Approximate wall-clock times on a modern single core:

| Symbol | j ~ 100 | j ~ 1000 | j ~ 5000 |
|---|---|---|---|
| 3j / 6j | < 1 ms | < 1 s | minutes |
| 9j | ~ 1 ms | ~ hours | impractical |

The 9j is the most expensive because it multiplies three large bigints at every
k step.  For the highest 9j angular momenta supported by the prime table
(j ~ 5000), a single evaluation may take hours or more.

### Argument type

All C API arguments are plain `int`.  Angular momenta up to j = 6666
correspond to `tj = 13332`, well within `INT_MAX`.  No overflow can occur in
intermediate index calculations for inputs within the correctness range.

### Thread safety

`primes_init()` (called automatically on the first symbol evaluation) is not
thread-safe: it writes global arrays on first call.  Call it explicitly once
from the main thread before spawning workers, or accept the benign data race
(writes are deterministic and idempotent).

---

## License

BSD 3-Clause — see [LICENSE](../LICENSE).

Copyright (c) 2026 Susi Lehtola.
