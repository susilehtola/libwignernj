/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * libquadmath / __float128 interface verification.
 *
 * Compares the _q variants against:
 *   - closed-form analytic values where available, at 4-ulp quad
 *     tolerance (~ 4 * FLT128_EPSILON);
 *   - the long-double variants, at 2-ulp long-double tolerance, since
 *     long double has fewer bits than float128 and the two cannot be
 *     expected to agree past long-double precision.
 *
 * The expected behaviour is that whenever long double and quad share an
 * analytic value, the quad result has at least as many correct bits as
 * the long-double result.
 */
#include "../include/wignernj.h"
#include "../include/wignernj_quadmath.h"
#include <quadmath.h>
#include <float.h>
#include <stdio.h>

static int g_run = 0, g_pass = 0, g_fail = 0;

static int near_q(__float128 a, __float128 b, __float128 tol)
{
    __float128 d  = fabsq(a - b);
    __float128 sc = fabsq(b) > 1e-300Q ? fabsq(b) : 1.0Q;
    return d <= tol * sc;
}

#define EXPECT_NEAR_Q(expr, ref, tol) do {                                 \
    __float128 _v = (expr);                                                \
    __float128 _r = (ref);                                                 \
    g_run++;                                                               \
    if (near_q(_v, _r, (tol))) { g_pass++; }                               \
    else {                                                                 \
        char _bv[64], _br[64];                                             \
        quadmath_snprintf(_bv, sizeof _bv, "%.34Qg", _v);                  \
        quadmath_snprintf(_br, sizeof _br, "%.34Qg", _r);                  \
        fprintf(stderr,                                                    \
                "FAIL  %s:%d  %s vs %s\n      got=%s\n expected=%s\n",     \
                __FILE__, __LINE__, #expr, #ref, _bv, _br);                \
        g_fail++;                                                          \
    }                                                                      \
} while(0)

int main(void)
{
    /* Quad ulp tolerance: FLT128_EPSILON = 2^-112 ≈ 1.93e-34 */
    const __float128 q_tol = 4.0Q * FLT128_EPSILON;
    /* Long-double-precision floor (used when comparing against _l results) */
    const __float128 ld_tol = 2.0Q * (__float128)LDBL_EPSILON;

    /* ── Closed-form analytic values, at quad tolerance ───────────────── */

    /* 3j(1,1,0; 0,0,0) = -1/sqrt(3) */
    EXPECT_NEAR_Q(wigner3j_q(2,2,0, 0,0,0), -1.0Q / sqrtq(3.0Q), q_tol);

    /* 3j(50 50 0; 0 0 0) = 1/sqrt(101) */
    EXPECT_NEAR_Q(wigner3j_q(100,100,0, 0,0,0),
                  1.0Q/sqrtq(101.0Q), q_tol);

    /* Very large j: exercises multi-word bigint conversion */
    EXPECT_NEAR_Q(wigner3j_q(1000,1000,0, 0,0,0),
                  1.0Q/sqrtq(1001.0Q), q_tol);

    /* 6j{1,1,1;1,1,1} = 1/6 */
    EXPECT_NEAR_Q(wigner6j_q(2,2,2, 2,2,2), 1.0Q/6.0Q, q_tol);

    /* 6j{j j 0; j j 0} = (-1)^(2j+0)/(2j+1) = +1/3 for j=1 */
    EXPECT_NEAR_Q(wigner6j_q(2,2,0, 2,2,0), 1.0Q/3.0Q, q_tol);

    /* CG(1/2 1/2; 1/2 -1/2 | 1 0) = 1/sqrt(2) */
    EXPECT_NEAR_Q(clebsch_gordan_q(1, 1, 1,-1, 2, 0),
                  1.0Q/sqrtq(2.0Q), q_tol);

    /* ── Selection-rule zeros ─────────────────────────────────────────── */

    EXPECT_NEAR_Q(wigner3j_q(2,2,2, 4,0,-4), 0.0Q, q_tol);
    EXPECT_NEAR_Q(wigner6j_q(2,2,2, 2,2,5),  0.0Q, q_tol);
    EXPECT_NEAR_Q(gaunt_q   (2,1, 2,1, 2,1), 0.0Q, q_tol);  /* m sum != 0 */

    /* ── Quad agrees with long double to long-double precision ────────── */
    /* (Quad will typically be more accurate; tolerance is the LD floor.) */

    EXPECT_NEAR_Q(wigner9j_q(2,2,4, 2,2,0, 4,4,4),
                  (__float128)wigner9j_l(2,2,4, 2,2,0, 4,4,4), ld_tol);

    EXPECT_NEAR_Q(wigner3j_q(20, 20, 20, 4, -8, 4),
                  (__float128)wigner3j_l(20, 20, 20, 4, -8, 4), ld_tol);

    EXPECT_NEAR_Q(gaunt_real_q(4, 0, 4, 0, 4, 0),
                  (__float128)gaunt_real_l(4, 0, 4, 0, 4, 0), ld_tol);

    EXPECT_NEAR_Q(racah_w_q(8, 6, 4, 6, 8, 10),
                  (__float128)racah_w_l(8, 6, 4, 6, 8, 10), ld_tol);

    /* Fano X = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * 9j.
     * Closed-form check: at all-equal-tj=4 the explicit normalisation is
     * exactly 5^2 = 25 times the underlying 9j. */
    EXPECT_NEAR_Q(fano_x_q(4,4,4, 4,4,4, 4,4,4),
                  25.0Q * wigner9j_q(4,4,4, 4,4,4, 4,4,4), q_tol);

    /* Cross-precision: quad agrees with long-double on a non-trivial value. */
    EXPECT_NEAR_Q(fano_x_q(8,6,2, 6,8,2, 2,2,4),
                  (__float128)fano_x_l(8,6,2, 6,8,2, 2,2,4), ld_tol);

    /* Large-j 9j stresses three independent bigints in the float128
     * conversion path. */
    EXPECT_NEAR_Q(wigner9j_q(20,20,20, 20,20,20, 20,20,20),
                  (__float128)wigner9j_l(20,20,20, 20,20,20, 20,20,20),
                  ld_tol);

    /* Real <-> complex Y_lm basis-overlap matrix at quad precision.  For
     * l = 1, the (m_r=+1, m_c=-1) entry is 1/sqrt(2) (real) and the
     * (m_r=-1, m_c=-1) entry is 1/sqrt(2) (imaginary). */
    {
        wignernj_cfloat128_t C[9];
        __float128 *F = (__float128 *)C;
        wignernj_real_ylm_in_complex_ylm_q(1, C);
        const __float128 inv_sqrt2 = 1.0Q / sqrtq(2.0Q);
        /* Column-major (m_r=+1, m_c=-1) -> idx 2*(0*3 + 2) = 4. */
        EXPECT_NEAR_Q(F[4], inv_sqrt2, q_tol);
        /* Column-major (m_r=-1, m_c=-1) -> idx 2*(0*3 + 0) + 1 = 1. */
        EXPECT_NEAR_Q(F[1], inv_sqrt2, q_tol);
        /* C[0,0] real part at idx 2*(1*3 + 1) = 8 == 1. */
        EXPECT_NEAR_Q(F[8], 1.0Q, q_tol);
    }

    /* ── Corner case: all j=0 → exact 1.0 ─────────────────────────────── */

    EXPECT_NEAR_Q(wigner3j_q(0,0,0, 0,0,0), 1.0Q, 0.0Q);
    EXPECT_NEAR_Q(wigner6j_q(0,0,0, 0,0,0), 1.0Q, 0.0Q);
    EXPECT_NEAR_Q(clebsch_gordan_q(0,0,0,0,0,0), 1.0Q, 0.0Q);

    printf("%s: %d/%d passed%s\n", __FILE__, g_pass, g_run,
           g_fail ? "  (FAILED)" : "");
    return g_fail ? 1 : 0;
}
