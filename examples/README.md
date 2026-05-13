# libwignernj examples

Single-file demonstrations of the libwignernj public API in each of
the four supported language bindings.  Every `all_symbols.*` calls
the nine public families exposed by the library — `wigner3j`,
`wigner6j`, `wigner9j`, `clebsch_gordan`, `racah_w`, `fano_x`,
`gaunt`, `gaunt_real`, and `real_ylm_in_complex_ylm` — at small
textbook-scale arguments and verifies the results against analytic
references (factors of 1/√3, 1/6, √(2/3), 1/√(5π), 1/√2, …).

| Language | File                                 | API surface                                |
|----------|--------------------------------------|--------------------------------------------|
| C        | `c/all_symbols.c`                    | `#include "wignernj.h"`                      |
| C++      | `cpp/all_symbols.cpp`                | `#include "wignernj.hpp"`                    |
| Fortran  | `fortran/all_symbols.f90`            | `use wignernj` from `libwignernj_f03`        |
| Python   | `python/all_symbols.py`              | `import wignernj`                            |

Parallel focused demos of the real ↔ complex Y_lm change of basis
(building the orbital `l_z` operator in the real-Y basis from its
diagonal complex-basis form via the similarity transform
`O_real = conj(C) @ O_complex @ transpose(C)`) ship alongside the
`all_symbols.*` files:

| Language | File                                   |
|----------|----------------------------------------|
| C        | `c/real_basis_lz.c`                    |
| C++      | `cpp/real_basis_lz.cpp`                |
| Fortran  | `fortran/real_basis_lz.f90`            |
| Python   | `python/real_basis_lz.py`              |

All four reproduce the textbook `l_z = ((0,0,+i),(0,0,0),(−i,0,0))`
at l = 1 with Hermiticity residual 0.

## Building and running from this source tree

The examples are built as part of the standard CMake build whenever
`BUILD_EXAMPLES=ON` (the default), and are registered as ctest tests:

```sh
cmake -B build
cmake --build build
ctest --test-dir build -R '^example_'      # run only the example tests
```

Each example exits with status 0 on success and prints a per-symbol
table of computed-vs-expected values.  Non-zero exit indicates a
tolerance violation (1e-14 against the analytic reference).

The Python example is registered only when the Python extension is
built in-tree (`-DBUILD_PYTHON=ON`); the source file ships in the
distribution either way and can be run directly against an installed
package (`pip install wignernj`) by

```sh
python examples/python/all_symbols.py
```

## Building against an installed libwignernj

Each example compiles standalone against the installed library:

```sh
# C
cc -o all_symbols_c c/all_symbols.c -lwignernj -lm

# C++
c++ -std=c++11 -o all_symbols_cpp cpp/all_symbols.cpp -lwignernj -lm

# Fortran (note the f03 wrapper)
gfortran -o all_symbols_f fortran/all_symbols.f90 \
        -lwignernj_f03 -lwignernj -lm

# Python (requires pip-installed wignernj package)
python python/all_symbols.py
```

## Argument-passing conventions

The C, C++, and Python APIs each have their own convention for
half-integer arguments:

* **C / Fortran C-binding interfaces** for the coupling-coefficient
  routines (`wigner3j`, `wigner6j`, `wigner9j`, `clebsch_gordan`,
  `racah_w`, `fano_x`, `gaunt`, `gaunt_real`) take `2j` as an
  integer, so `j = 1/2` is `tj = 1` and `j = 1` is `tj = 2`.  This
  avoids floating-point representation of half-integer quantum
  numbers entirely.  `wignernj_real_ylm_in_complex_ylm` is the lone
  exception — its argument is always an integer orbital angular
  momentum and is therefore taken as plain `l`.
* **C++ wrapper** (`wignernj::symbol3j(double, …)`) and **Fortran
  module wrappers** (`w3j`, `w6j`, …) accept `j` as a real number;
  the implementation doubles it internally and forwards to the C
  API.  Non-half-integer real arguments raise `std::invalid_argument`
  in C++ and are reduced to the nearest half-integer in Fortran.
* **Python wrapper** accepts integers, floats, or `fractions.Fraction`
  objects directly.  `wignernj.wigner3j(0.5, 0.5, 1, 0.5, -0.5, 0)`
  returns `1/√6`.

The examples in this directory all use the most idiomatic surface for
their language — `tj` integers in C, real-valued overloads in C++ and
Fortran, native ints in Python.
