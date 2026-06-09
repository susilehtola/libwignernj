// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// libwignernj C++ API quadmath (__float128) demonstration.
//
// Exercises the wignernj::symbol3j<__float128>(), wignernj::cg<__float128>(),
// ... template specialisations added by wignernj_quadmath.hpp.  Calls one
// representative entry point from each public symbol family at binary128
// precision and prints the result alongside an analytic reference value.
//
// Requires libwignernj built with -DBUILD_QUADMATH=ON.  __float128 is a
// GCC extension, and <quadmath.h> macros (FLT128_EPSILON, ...) expand to
// `Q`-suffixed literals; this file therefore relies on GNU C++ mode
// (-std=gnu++11).  Build (out-of-tree, against an installed libwignernj):
//     c++ -std=gnu++11 -o quadmath_cpp quadmath.cpp -lwignernj -lquadmath -lm
#include "wignernj_quadmath.hpp"

#include <quadmath.h>
#include <cstdio>

static int g_fail = 0;

static void check(const char *label, __float128 computed, __float128 expected,
                  __float128 tol)
{
    char cb[64], eb[64];
    quadmath_snprintf(cb, sizeof cb, "%+.34Qg", computed);
    quadmath_snprintf(eb, sizeof eb, "%+.34Qg", expected);
    std::printf("  %-44s = %s\n", label, cb);
    std::printf("  %-44s   (expected %s)\n", "", eb);

    __float128 diff = fabsq(computed - expected);
    __float128 ref  = fabsq(expected) > 1e-300Q ? fabsq(expected) : 1.0Q;
    if (diff > tol * ref) {
        char db[64];
        quadmath_snprintf(db, sizeof db, "%.3Qg", diff);
        std::fprintf(stderr, "  FAIL: %s exceeds tolerance (diff = %s)\n",
                     label, db);
        ++g_fail;
    }
}

int main()
{
    const __float128 tol = 4.0Q * FLT128_EPSILON;

    std::printf("libwignernj C++ quadmath API demonstration\n");
    std::printf("------------------------------------------\n");

    // 3j(1,1,0; 0,0,0) = -1/sqrt(3)
    check("symbol3j<__float128>(1,1,0; 0,0,0)",
          wignernj::symbol3j<__float128>(2, 2, 0, 0, 0, 0),
          -1.0Q / sqrtq(3.0Q), tol);

    // Real-valued (double-arg) overload at __float128: agrees with the
    // integer-argument call to the last quad ulp.
    check("symbol3j<__float128>(0.5,0.5,0.0; 0.5,-0.5,0.0)",
          wignernj::symbol3j<__float128>(0.5, 0.5, 0.0,
                                         0.5,-0.5, 0.0),
          wignernj::symbol3j<__float128>(1, 1, 0, 1,-1, 0),
          tol);

    // 6j{1,1,1; 1,1,1} = 1/6
    check("symbol6j<__float128>{1,1,1; 1,1,1}",
          wignernj::symbol6j<__float128>(2, 2, 2, 2, 2, 2),
          1.0Q / 6.0Q, tol);

    // CG(1/2 1/2; 1/2 -1/2 | 1 0) = 1/sqrt(2)
    check("cg<__float128>(1/2,1/2; 1/2,-1/2 | 1,0)",
          wignernj::cg<__float128>(1, 1, 1, -1, 2, 0),
          1.0Q / sqrtq(2.0Q), tol);

    // Real Gaunt at all-m=0 equals complex Gaunt.
    {
        __float128 gr = wignernj::gauntreal<__float128>(2, 0, 2, 0, 4, 0);
        __float128 gc = wignernj::gaunt<__float128>    (2, 0, 2, 0, 4, 0);
        check("gauntreal<__float128>(.,0; .,0; .,0) - gaunt<__float128>",
              gr - gc, 0.0Q, tol);
    }

    // real_ylm_in_complex_ylm non-template overload at l=1: C[+1,-1]
    // real part is 1/sqrt(2).  Column-major (m_r=+1, m_c=-1) lives at
    // complex slot 0*3+2 = 2, whose real part is float index 4.
    {
        wignernj_cfloat128_t C[9];
        wignernj::real_ylm_in_complex_ylm(1, C);
        const __float128 *F = reinterpret_cast<const __float128 *>(C);
        check("real_ylm_in_complex_ylm[l=1] re(+1,-1)",
              F[4], 1.0Q / sqrtq(2.0Q), tol);
    }

    // Vector-returning convenience form.
    {
        std::vector<wignernj_cfloat128_t> C =
            wignernj::real_ylm_in_complex_ylm_q(1);
        std::printf("  vector form size = %zu (expected 9)\n", C.size());
        if (C.size() != 9) ++g_fail;
    }

    if (g_fail) {
        std::fprintf(stderr, "%d check(s) failed\n", g_fail);
        return 1;
    }
    std::printf("\nAll quadmath checks passed.\n");
    return 0;
}
