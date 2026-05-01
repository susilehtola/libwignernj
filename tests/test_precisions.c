/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Verify float, double, and long double variants agree with expected values.
 */
#include "run_tests.h"
#include "../include/wigner.h"
#include <float.h>
#include <math.h>

int main(void)
{
    /* Expected: wigner3j(1,1,0; 0,0,0) = -1/sqrt(3) */
    double expected = -1.0 / sqrt(3.0);

    /* float: tolerance ~ 2 * FLT_EPSILON */
    {
        float vf = wigner3j_f(2,2,0, 0,0,0);
        float ef = (float)expected;
        float tol = 2.0f * FLT_EPSILON * fabsf(ef);
        TEST_ASSERT(fabsf(vf - ef) <= tol || fabsf(vf - ef) < 1e-6f);
    }

    /* double: tolerance ~ 2 * DBL_EPSILON */
    {
        double vd = wigner3j(2,2,0, 0,0,0);
        TEST_NEAR(vd, expected, 2.0 * DBL_EPSILON);
    }

    /* long double: tolerance ~ 2 * LDBL_EPSILON */
    {
        long double vl = wigner3j_l(2,2,0, 0,0,0);
        long double el = -1.0L / sqrtl(3.0L);
        long double tol = 2.0L * LDBL_EPSILON;
        TEST_ASSERT(fabsl(vl - el) <= tol * fabsl(el) + 1e-18L);
    }

    /* All three agree */
    {
        float  vf = wigner3j_f(4,4,0, 0,0,0);
        double vd = wigner3j  (4,4,0, 0,0,0);
        long double vl = wigner3j_l(4,4,0, 0,0,0);
        TEST_NEAR((double)vf, vd, 1e-6);
        TEST_NEAR((double)vl, vd, 2.0 * DBL_EPSILON);
    }

    /* 6j symbol: all precisions */
    {
        float  vf = wigner6j_f(2,2,2, 2,2,2);
        double vd = wigner6j  (2,2,2, 2,2,2);
        long double vl = wigner6j_l(2,2,2, 2,2,2);
        TEST_DBL(vd, 1.0/6.0);
        TEST_NEAR((double)vf,  1.0f/6.0f, 1e-6);
        TEST_NEAR((double)vl, 1.0/6.0, 2.0 * DBL_EPSILON);
    }

    /* 9j symbol: all precisions consistent */
    {
        float  vf = wigner9j_f(1,1,2, 1,1,0, 2,2,2);
        double vd = wigner9j  (1,1,2, 1,1,0, 2,2,2);
        long double vl = wigner9j_l(1,1,2, 1,1,0, 2,2,2);
        TEST_DBL(vd, sqrt(6.0)/18.0);
        TEST_NEAR((double)vf, (float)(sqrt(6.0)/18.0), 5e-6);
        TEST_NEAR((double)vl, sqrt(6.0)/18.0, 2.0*DBL_EPSILON);
    }

    /* Large j: (50 50 0; 0 0 0) = 1/sqrt(101) */
    {
        double vd = wigner3j  (100,100,0, 0,0,0);
        long double vl = wigner3j_l(100,100,0, 0,0,0);
        TEST_DBL(vd,  1.0/sqrt(101.0));
        TEST_NEAR((double)vl, 1.0/sqrt(101.0), 2.0*DBL_EPSILON);
    }

    /* Very large j: (500 500 0; 0 0 0) = 1/sqrt(1001)
     * exercises 1000! prime factorization and multi-word bigint arithmetic */
    {
        double vd = wigner3j(1000,1000,0, 0,0,0);
        TEST_DBL(vd, 1.0/sqrt(1001.0));
    }

    /* Large 6j: {50 50 0; 50 50 50} = (-1)^(3*50)/(2*50+1) = (-1)^150/101 = +1/101 */
    {
        double vd = wigner6j(100,100,0, 100,100,100);
        TEST_DBL(vd, 1.0/101.0);
    }

    /* Corner case: all j=0 → exact 1 */
    TEST_ASSERT(wigner3j(0,0,0, 0,0,0) == 1.0);
    TEST_ASSERT(wigner6j(0,0,0, 0,0,0) == 1.0);
    TEST_ASSERT(clebsch_gordan(0,0,0,0,0,0) == 1.0);

    /* Selection-rule zero: all precisions return exactly 0 */
    TEST_ABS(wigner3j_f(2,2,2, 4,0,-4),   0.0, 1e-30);
    TEST_ABS(wigner3j  (2,2,2, 4,0,-4),   0.0, 1e-30);
    TEST_ABS(wigner3j_l(2,2,2, 4,0,-4),   0.0, 1e-30);
    TEST_ABS(wigner6j_f(2,2,2, 2,2,5),    0.0, 1e-30); /* triangle violated */
    TEST_ABS(wigner6j  (2,2,2, 2,2,5),    0.0, 1e-30);

    /* Sum-cancellation zero: passes selection rules but Racah sum cancels.
     * 3j(1,1,1;0,0,0): j1+j2+j3=3 odd, all-zero-m parity rule → 0 */
    TEST_ABS(wigner3j_f(2,2,2, 0,0,0),    0.0, 1e-30);
    TEST_ABS(wigner3j  (2,2,2, 0,0,0),    0.0, 1e-30);
    TEST_ABS(wigner3j_l(2,2,2, 0,0,0),    0.0, 1e-30);

    /* ── Symmetry relations ─────────────────────────────────────────── */

    /* 3j column permutation: even permutation (123→231) leaves value unchanged */
    {
        double v123 = wigner3j(4,6,2, 2,-4,2);
        double v231 = wigner3j(6,2,4, -4,2,2);
        double v312 = wigner3j(2,4,6, 2,2,-4);
        TEST_NEAR(v123, v231, 2e-14);
        TEST_NEAR(v123, v312, 2e-14);
    }

    /* 3j odd permutation: (j1+j2+j3) = (4+6+2)/2 = 6 (even) → same sign */
    {
        double v123 = wigner3j(4,6,2, 2,-4,2);
        double v132 = wigner3j(4,2,6, 2,2,-4);  /* swap cols 2,3 → (-1)^6 = +1 */
        TEST_NEAR(v123, v132, 2e-14);
    }

    /* 3j odd permutation: (j1+j2+j3)=(3+1+2)/2=3 (odd) → swap cols 2,3 flips sign */
    {
        /* (3/2, 1/2, 1; 1/2, -1/2, 0): non-zero, triangle is satisfied */
        double v123 = wigner3j(3,1,2,  1,-1, 0);
        double v132 = wigner3j(3,2,1,  1, 0,-1);  /* (-1)^3 = -1 */
        TEST_NEAR(v123, -v132, 2e-14);
    }

    /* 3j m-reversal: (j1 j2 j3;m1 m2 m3) = (-1)^(j1+j2+j3) * (...;-m1 -m2 -m3) */

    /* j1+j2+j3 = 2+3+1 = 6 (even): same value */
    {
        double vp = wigner3j(4,6,2,  2,-4,2);
        double vm = wigner3j(4,6,2, -2, 4,-2);
        TEST_NEAR(vp, vm, 2e-14);
    }

    /* j1+j2+j3 = 1+1+1 = 3 (odd): negated */
    {
        double vp = wigner3j(2,2,2,  2,0,-2);
        double vm = wigner3j(2,2,2, -2,0, 2);
        TEST_NEAR(vp, -vm, 2e-14);
    }

    /* 6j permutation: swap first two columns → same value */
    {
        double v1 = wigner6j(4,6,2, 2,4,6);
        double v2 = wigner6j(6,4,2, 4,2,6);
        TEST_NEAR(v1, v2, 2e-14);
    }

    /* 6j cyclic column permutation: {a b c; d e f} = {b c a; e f d} */
    {
        double v1 = wigner6j(2,4,6, 2,4,6);  /* {1 2 3; 1 2 3} */
        double v2 = wigner6j(4,6,2, 4,6,2);  /* {2 3 1; 2 3 1} — cyclic shift */
        TEST_NEAR(v1, v2, 2e-14);
    }

    /* 9j transposition symmetry: 9j{j11..j33} = 9j{transpose} (non-symmetric matrix) */
    {
        /* {1 1 2; 1 1 0; 0 2 2} in 2j notation: {2 2 4; 2 2 0; 0 4 4} */
        double v  = wigner9j(2,2,4, 2,2,0, 0,4,4);
        /* Transpose of j-matrix: {1 1 0; 1 1 2; 2 0 2} → {2 2 0; 2 2 4; 4 0 4} */
        double vt = wigner9j(2,2,0, 2,2,4, 4,0,4);
        TEST_NEAR(v, vt, 2e-14);
    }

    /* 9j row/column reversal: reverse all rows AND columns → (-1)^(sum of all j) */
    {
        /* {1 1 2; 1 1 0; 2 0 2}: sum = (2+2+4+2+2+0+4+0+4)/2 = 10 (even) */
        double v  = wigner9j(2,2,4, 2,2,0, 4,0,4);
        double vr = wigner9j(4,0,4, 2,2,0, 2,2,4);  /* rows reversed */
        TEST_NEAR(v, vr, 2e-14);
    }

    /* CG orthogonality: sum_m1,m2 |<j1m1;j2m2|JM>|^2 = 1 for fixed J,M */
    {
        /* j1=j2=1/2, J=1, M=0: two terms: <1/2,1/2;1/2,-1/2|1,0> and <1/2,-1/2;1/2,1/2|1,0> */
        double c1 = clebsch_gordan(1, 1, 1,-1, 2, 0);  /* <1/2 1/2; 1/2 -1/2 | 1 0> */
        double c2 = clebsch_gordan(1,-1, 1, 1, 2, 0);  /* <1/2 -1/2; 1/2 1/2 | 1 0> */
        TEST_NEAR(c1*c1 + c2*c2, 1.0, 2e-14);
    }

    /* CG vs 3j: <j1 m1; j2 m2 | J M> = (-1)^(j1-j2+M) sqrt(2J+1) * 3j(j1,j2,J;m1,m2,-M) */
    {
        int tj1=2, tm1=1, tj2=4, tm2=-1, tJ=4, tM=0;
        double cg  = clebsch_gordan(tj1, tm1, tj2, tm2, tJ, tM);
        double w3  = wigner3j(tj1, tj2, tJ, tm1, tm2, -tM);
        double phase = (((tj1 - tj2 + tM) / 2) & 1) ? -1.0 : 1.0;
        double cg_expected = phase * sqrt((double)(tJ + 1)) * w3;
        TEST_NEAR(cg, cg_expected, 2e-14);
    }

    /* ── Racah W zeros and precision ────────────────────────────────── */

    /* Triangle violated (tj12 > tj1+tj2) → 0 */
    TEST_ABS(racah_w_f(2,2,2, 2,5,2), 0.0, 1e-30);
    TEST_ABS(racah_w  (2,2,2, 2,5,2), 0.0, 1e-30);
    TEST_ABS(racah_w_l(2,2,2, 2,5,2), 0.0, 1e-30);

    /* All-precision consistency: W(4,3,2,3;4,5) — all triangles satisfied */
    {
        float      vf = racah_w_f(8,6,4,6, 8,10);
        double     vd = racah_w  (8,6,4,6, 8,10);
        long double vl = racah_w_l(8,6,4,6, 8,10);
        TEST_NEAR((double)vf, vd, 1e-6);
        TEST_NEAR((double)vl, vd, 2.0 * DBL_EPSILON);
    }

    /* Racah W vs 6j: W(j1,j2,J,j3;j12,j23) = (-1)^(j1+j2+J+j3) * {j1 j2 j12; j3 J j23} */
    {
        int tj1=6, tj2=4, tJ=4, tj3=6, tj12=6, tj23=8;
        double w   = racah_w(tj1, tj2, tJ, tj3, tj12, tj23);
        double w6  = wigner6j(tj1, tj2, tj12, tj3, tJ, tj23);
        int ph     = ((((tj1 + tj2 + tJ + tj3) / 2) & 1)) ? -1 : 1;
        TEST_NEAR(w, ph * w6, 2e-14);
    }

    /* ── Gaunt zeros ─────────────────────────────────────────────────── */

    /* m-conservation violated → 0 */
    TEST_ABS(gaunt_f(2,1, 2,1, 2,1),  0.0, 1e-30);
    TEST_ABS(gaunt  (2,1, 2,1, 2,1),  0.0, 1e-30);
    TEST_ABS(gaunt_l(2,1, 2,1, 2,1),  0.0, 1e-30);

    /* l1+l2+l3 odd → 3j(l1,l2,l3;0,0,0) = 0 → gaunt = 0 */
    TEST_ABS(gaunt_f(2,1, 2,-1, 2,0), 0.0, 1e-30);
    TEST_ABS(gaunt  (2,1, 2,-1, 2,0), 0.0, 1e-30);
    TEST_ABS(gaunt_l(2,1, 2,-1, 2,0), 0.0, 1e-30);

    /* ── Gaunt all-precision consistency ─────────────────────────────── */

    /* G(4,2, 2,0, 4,-2): l1+l2+l3=5 (even), m-sum=0, non-zero */
    {
        float      vf = gaunt_f(8, 4, 4, 0, 8, -4);
        double     vd = gaunt  (8, 4, 4, 0, 8, -4);
        long double vl = gaunt_l(8, 4, 4, 0, 8, -4);
        TEST_NEAR((double)vf, vd, 1e-6);
        TEST_NEAR((double)vl, vd, 2.0 * DBL_EPSILON);
    }

    SUMMARY();
}
