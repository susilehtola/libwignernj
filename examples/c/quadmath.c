/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * libwignernj C API quadmath (__float128) demonstration.
 *
 * Calls one representative entry point from each public symbol family at
 * binary128 precision and prints the result alongside an analytic
 * reference value formatted with the libquadmath %Qg conversion.
 * Exits 0 on success, non-zero if any computed value disagrees with its
 * reference by more than 4 * FLT128_EPSILON of the reference.
 *
 * Requires libwignernj built with -DBUILD_QUADMATH=ON.
 *
 * Build (out-of-tree, against an installed libwignernj):
 *     cc -o quadmath_c quadmath.c -lwignernj -lquadmath -lm
 */
#include "wignernj_quadmath.h"   /* pulls in wignernj.h transitively */

#include <quadmath.h>
#include <float.h>
#include <stdio.h>

static int g_fail = 0;

static void check(const char *label, __float128 computed, __float128 expected,
                  __float128 tol)
{
    char cb[64], eb[64];
    quadmath_snprintf(cb, sizeof cb, "%+.34Qg", computed);
    quadmath_snprintf(eb, sizeof eb, "%+.34Qg", expected);
    printf("  %-44s = %s\n", label, cb);
    printf("  %-44s   (expected %s)\n", "", eb);

    __float128 diff = fabsq(computed - expected);
    __float128 ref  = fabsq(expected) > 1e-300Q ? fabsq(expected) : 1.0Q;
    if (diff > tol * ref) {
        char db[64];
        quadmath_snprintf(db, sizeof db, "%.3Qg", diff);
        fprintf(stderr, "  FAIL: %s exceeds tolerance (diff = %s)\n", label, db);
        ++g_fail;
    }
}

int main(void)
{
    const __float128 tol = 4.0Q * FLT128_EPSILON;

    printf("libwignernj C quadmath API demonstration -- one call per symbol family\n");
    printf("----------------------------------------------------------------------\n");

    /* 3j(1,1,0; 0,0,0) = -1/sqrt(3) */
    check("wigner3j_q(1,1,0; 0,0,0)",
          wigner3j_q(2, 2, 0,  0, 0, 0),
          -1.0Q / sqrtq(3.0Q), tol);

    /* 6j{1,1,1; 1,1,1} = 1/6 */
    check("wigner6j_q{1,1,1; 1,1,1}",
          wigner6j_q(2, 2, 2,  2, 2, 2),
          1.0Q / 6.0Q, tol);

    /* 9j{1/2,1/2,1; 1/2,1/2,0; 1,1,1} = sqrt(6)/18 */
    check("wigner9j_q{1/2,1/2,1; 1/2,1/2,0; 1,1,1}",
          wigner9j_q(1, 1, 2,  1, 1, 0,  2, 2, 2),
          sqrtq(6.0Q) / 18.0Q, tol);

    /* CG(1/2 1/2; 1/2 -1/2 | 1 0) = 1/sqrt(2) */
    check("clebsch_gordan_q(1/2,1/2; 1/2,-1/2 | 1,0)",
          clebsch_gordan_q(1, 1, 1, -1, 2, 0),
          1.0Q / sqrtq(2.0Q), tol);

    /* Racah W(1,1,1,1; 1,1) = (-1)^4 * 6j{1,1,1; 1,1,1} = 1/6 */
    check("racah_w_q(1,1,1,1; 1,1)",
          racah_w_q(2, 2, 2, 2, 2, 2),
          1.0Q / 6.0Q, tol);

    /* Fano X at all-tj=4 = 25 * 9j{4,4,4; 4,4,4; 4,4,4} */
    check("fano_x_q(2,2,2; 2,2,2; 2,2,2)",
          fano_x_q(4, 4, 4,  4, 4, 4,  4, 4, 4),
          25.0Q * wigner9j_q(4,4,4, 4,4,4, 4,4,4), tol);

    /* Real Gaunt at all-m=0 equals complex Gaunt. */
    check("gaunt_real_q(1,0; 1,0; 2,0) - gaunt_q(.)",
          gaunt_real_q(2, 0, 2, 0, 4, 0) - gaunt_q(2, 0, 2, 0, 4, 0),
          0.0Q, tol);

    /* Real <-> complex Y_lm basis at l=1: C[+1,-1] real part is 1/sqrt(2). */
    {
        wignernj_cfloat128_t C[9];
        wignernj_real_ylm_in_complex_ylm_q(1, C);
        /* Column-major (m_r=+1, m_c=-1) at idx 0*3+2 = 2. */
        __float128 re = ((__float128 *)C)[4];
        check("real_ylm_in_complex_ylm_q[l=1] re(+1,-1)",
              re, 1.0Q / sqrtq(2.0Q), tol);
    }

    if (g_fail) {
        fprintf(stderr, "%d check(s) failed\n", g_fail);
        return 1;
    }
    printf("\nAll quadmath checks passed.\n");
    return 0;
}
