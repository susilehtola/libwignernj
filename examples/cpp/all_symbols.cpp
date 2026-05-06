// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// libwignernj C++ API demonstration.
//
// Calls every public symbol family exposed by wignernj.hpp once with small,
// textbook-scale arguments and prints the result alongside an analytic
// reference value.  The program exits 0 on success, non-zero if any
// computed value disagrees with its reference by more than 1e-14.
//
// The C++ wrapper is header-only and links against the C library; it
// adds two conveniences over the C API:
//   1. a real-valued (double j) overload that accepts j directly rather
//      than 2*j, and throws std::invalid_argument on non-half-integer
//      arguments;
//   2. compile-time selection of float/double/long double through a
//      template parameter.
//
// Build (out-of-tree, against an installed libwignernj):
//     c++ -std=c++11 -o all_symbols_cpp all_symbols.cpp -lwignernj -lm
#include "wignernj.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

static int check(const char *label, double computed, double expected,
                 double tol = 1e-14)
{
    std::printf("  %-38s = %+.15f   (expected %+.15f)\n",
                label, computed, expected);
    if (std::fabs(computed - expected) > tol) {
        std::fprintf(stderr,
                     "  FAIL: |%g - %g| = %g exceeds tolerance %g\n",
                     computed, expected, std::fabs(computed - expected), tol);
        return 1;
    }
    return 0;
}

int main()
{
    int failed = 0;

    std::printf("libwignernj C++ API demonstration -- one call per symbol family\n");
    std::printf("---------------------------------------------------------------\n");

    // 1. Wigner 3j symbol.   Real-valued overload: pass j directly.
    failed |= check("wignernj::symbol3j(1,1,0; 0,0,0)",
                    wignernj::symbol3j<double>(1.0, 1.0, 0.0,
                                              0.0, 0.0, 0.0),
                    -1.0 / std::sqrt(3.0));

    // 2. Wigner 6j symbol.
    failed |= check("wignernj::symbol6j{1,1,1; 1,1,1}",
                    wignernj::symbol6j<double>(1.0, 1.0, 1.0,
                                              1.0, 1.0, 1.0),
                    1.0 / 6.0);

    // 3. Wigner 9j symbol.
    failed |= check("wignernj::symbol9j{1,1,0; 1,1,0; 0,0,0}",
                    wignernj::symbol9j<double>(1.0, 1.0, 0.0,
                                              1.0, 1.0, 0.0,
                                              0.0, 0.0, 0.0),
                    1.0 / 3.0);

    // 4. Clebsch-Gordan coefficient.
    failed |= check("wignernj::cg(1,0; 1,0 | 2,0)",
                    wignernj::cg<double>(1.0, 0.0,  1.0, 0.0,  2.0, 0.0),
                    std::sqrt(2.0 / 3.0));

    // 5. Racah W coefficient.
    failed |= check("wignernj::racahw(1,1,1,1; 1,1)",
                    wignernj::racahw<double>(1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0),
                    1.0 / 6.0);

    // 6. Fano X coefficient.
    failed |= check("wignernj::fanox(1,1,1; 1,1,1; 1,1,2)",
                    wignernj::fanox<double>(1.0, 1.0, 1.0,
                                           1.0, 1.0, 1.0,
                                           1.0, 1.0, 2.0),
                    0.5);

    // 7. Gaunt coefficient (complex Y_l^m).
    failed |= check("wignernj::gaunt(1,0; 1,0; 2,0)",
                    wignernj::gaunt<double>(1.0, 0.0,  1.0, 0.0,  2.0, 0.0),
                    1.0 / std::sqrt(5.0 * M_PI));

    // 8. Gaunt over real spherical harmonics.
    failed |= check("wignernj::gauntreal(1,0; 1,0; 2,0)",
                    wignernj::gauntreal<double>(1.0, 0.0,  1.0, 0.0,  2.0, 0.0),
                    1.0 / std::sqrt(5.0 * M_PI));

    if (!failed)
        std::printf("\nAll symbols agree with their analytic references.\n");
    return failed;
}
