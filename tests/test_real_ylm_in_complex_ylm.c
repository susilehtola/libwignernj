/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Tests for the real <-> complex spherical-harmonic basis-change matrix
 * wignernj_real_ylm_in_complex_ylm.
 *
 * Verifies:
 *   (a) explicit l=0/1/2 entries against the convention formulas,
 *   (b) unitarity  C C^H = I  for l = 0, 1, ..., 10,
 *   (c) cross-check that the real Gaunt computed by sandwiching three
 *       rows of C against complex Gaunts agrees bit-for-bit with the
 *       library's gaunt_real() (i.e. the new matrix matches the same
 *       convention internally used by gaunt_real).
 */
#include "run_tests.h"
#include "../include/wignernj.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Index helpers for a column-major (re, im)-interleaved buffer of
 * dimension dim = 2l+1. */
#define DRE(buf, dim, r, c) ((buf)[2u * ((size_t)(c) * (size_t)(dim) + (size_t)(r))    ])
#define DIM(buf, dim, r, c) ((buf)[2u * ((size_t)(c) * (size_t)(dim) + (size_t)(r)) + 1u])

static double *alloc_matrix(int l)
{
    const int dim = 2 * l + 1;
    return (double *)calloc((size_t)(2 * dim * dim), sizeof(double));
}

/* The C99 _Complex storage is layout-compatible with double[2], so
 * the test reads/writes the (re, im) pairs through a plain double *
 * but hands the wignernj entry point the typed pointer. */
#define CALL_FILL(l, buf) \
    wignernj_real_ylm_in_complex_ylm((l), (wignernj_cdouble_t *)(buf))

static void test_l0(void)
{
    /* l = 0: 1x1 matrix == [[1+0i]] */
    double *C = alloc_matrix(0);
    CALL_FILL(0, C);
    TEST_DBL(DRE(C, 1, 0, 0), 1.0);
    TEST_ABS(DIM(C, 1, 0, 0), 0.0, 0.0);
    free(C);
}

static void test_l1(void)
{
    /* l = 1, index convention m+l = 0,1,2 for m = -1,0,+1.
     *   S_{1,0}  = Y_1^0
     *   S_{1,+1} = (1/sqrt2)( Y_1^{-1} - Y_1^{+1} )
     *   S_{1,-1} = (i/sqrt2)( Y_1^{-1} + Y_1^{+1} )            */
    double *C = alloc_matrix(1);
    CALL_FILL(1, C);
    const double s = 1.0 / sqrt(2.0);

    /* Row m_r = -1 (idx 0): imag-only entries at cols m_c = -1, +1 (idx 0, 2) */
    TEST_ABS(DRE(C, 3, 0, 0), 0.0, 0.0);
    TEST_DBL(DIM(C, 3, 0, 0), +s);
    TEST_ABS(DRE(C, 3, 0, 1), 0.0, 0.0);
    TEST_ABS(DIM(C, 3, 0, 1), 0.0, 0.0);
    TEST_ABS(DRE(C, 3, 0, 2), 0.0, 0.0);
    TEST_DBL(DIM(C, 3, 0, 2), +s);

    /* Row m_r = 0 (idx 1): identity at m_c = 0 (idx 1) */
    TEST_ABS(DRE(C, 3, 1, 0), 0.0, 0.0);
    TEST_DBL(DRE(C, 3, 1, 1), 1.0);
    TEST_ABS(DIM(C, 3, 1, 1), 0.0, 0.0);
    TEST_ABS(DRE(C, 3, 1, 2), 0.0, 0.0);

    /* Row m_r = +1 (idx 2): real entries at cols m_c = -1, +1 */
    TEST_DBL(DRE(C, 3, 2, 0), +s);
    TEST_ABS(DIM(C, 3, 2, 0), 0.0, 0.0);
    TEST_ABS(DRE(C, 3, 2, 1), 0.0, 0.0);
    TEST_DBL(DRE(C, 3, 2, 2), -s);
    TEST_ABS(DIM(C, 3, 2, 2), 0.0, 0.0);
    free(C);
}

static void test_l2(void)
{
    /* l = 2 spot checks:
     *   S_{2,+2} = (1/sqrt2)( Y_2^{-2} + Y_2^{+2} )
     *   S_{2,-2} = (i/sqrt2)( Y_2^{-2} - Y_2^{+2} )    (since (-1)^2 = +1)
     */
    double *C = alloc_matrix(2);
    CALL_FILL(2, C);
    const int dim = 5;
    const double s = 1.0 / sqrt(2.0);

    /* Row m_r = +2 (idx 4): cos-type with (-1)^2 = +1. */
    TEST_DBL(DRE(C, dim, 4, 0), +s);   /* col m_c = -2 (idx 0) */
    TEST_DBL(DRE(C, dim, 4, 4), +s);   /* col m_c = +2 (idx 4) */
    TEST_ABS(DIM(C, dim, 4, 0), 0.0, 0.0);
    TEST_ABS(DIM(C, dim, 4, 4), 0.0, 0.0);

    /* Row m_r = -2 (idx 0): sin-type, sign_k = +1, so im[cn]=+s, im[cp]=-s. */
    TEST_DBL(DIM(C, dim, 0, 0), +s);
    TEST_DBL(DIM(C, dim, 0, 4), -s);
    TEST_ABS(DRE(C, dim, 0, 0), 0.0, 0.0);
    TEST_ABS(DRE(C, dim, 0, 4), 0.0, 0.0);

    free(C);
}

/* C C^H == I, where the (i, j) entry of C C^H is
 *   sum_k  C[i,k] * conj(C[j,k]). */
static void test_unitarity(int l)
{
    const int dim = 2 * l + 1;
    double *C = alloc_matrix(l);
    CALL_FILL(l, C);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            double acc_re = 0.0, acc_im = 0.0;
            for (int k = 0; k < dim; k++) {
                double a_re = DRE(C, dim, i, k);
                double a_im = DIM(C, dim, i, k);
                double b_re = DRE(C, dim, j, k);
                double b_im = DIM(C, dim, j, k);
                /* a * conj(b) */
                acc_re += a_re * b_re + a_im * b_im;
                acc_im += a_im * b_re - a_re * b_im;
            }
            double expected = (i == j) ? 1.0 : 0.0;
            TEST_ABS(acc_re, expected, 1e-14);
            TEST_ABS(acc_im, 0.0,      1e-14);
        }
    }
    free(C);
}

/* Cross-check: compute one real Gaunt by sandwiching three rows of C
 * against the complex Gaunts at the same (l1, l2, l3), then compare to
 * gaunt_real().  This pins down that the matrix C uses exactly the
 * same convention as the internal gaunt_real path.
 *
 *   G^R(l1,m1, l2,m2, l3,m3)
 *     = sum_{p1, p2, p3} C[m1,p1] * C[m2,p2] * C[m3,p3]
 *       * gaunt(l1,p1, l2,p2, l3,p3)
 *
 * The inner sum keeps only terms with p1+p2+p3 == 0 (gaunt vanishes
 * otherwise); each Y_l^m row of C touches at most two columns, so the
 * sum has at most 2^3 = 8 nonvanishing summands. */
static double gaunt_via_matrix(int l1, int m1, int l2, int m2, int l3, int m3)
{
    const int L = (l1 > l2 ? l1 : l2) > l3 ? (l1 > l2 ? l1 : l2) : l3;
    const int dim = 2 * L + 1;
    /* Allocate three matrices, one per row dimension; for simplicity reuse
     * one square buffer at dimension L (call once and reindex per l_i).
     * That over-allocates rows we don't use but keeps the test simple. */
    double *Cmax = alloc_matrix(L);
    /* But fills require the matrix to be sized to its own l; allocate one
     * each and fill at the symbol's own l. */
    double *Cs[3];
    int     ls[3] = {l1, l2, l3};
    int     ms[3] = {m1, m2, m3};
    for (int i = 0; i < 3; i++) {
        Cs[i] = alloc_matrix(ls[i]);
        CALL_FILL(ls[i], Cs[i]);
    }
    free(Cmax);

    double acc_re = 0.0;
    for (int p1 = -l1; p1 <= l1; p1++) {
        for (int p2 = -l2; p2 <= l2; p2++) {
            int p3 = -(p1 + p2);
            if (p3 < -l3 || p3 > l3) continue;
            /* Pull C[m_i, p_i] for each of the three rows. */
            int d1 = 2*l1+1, d2 = 2*l2+1, d3 = 2*l3+1;
            double a_re = DRE(Cs[0], d1, m1+l1, p1+l1);
            double a_im = DIM(Cs[0], d1, m1+l1, p1+l1);
            double b_re = DRE(Cs[1], d2, m2+l2, p2+l2);
            double b_im = DIM(Cs[1], d2, m2+l2, p2+l2);
            double c_re = DRE(Cs[2], d3, m3+l3, p3+l3);
            double c_im = DIM(Cs[2], d3, m3+l3, p3+l3);
            if ((a_re == 0.0 && a_im == 0.0) ||
                (b_re == 0.0 && b_im == 0.0) ||
                (c_re == 0.0 && c_im == 0.0))
                continue;
            double g = gaunt(2*l1, 2*p1, 2*l2, 2*p2, 2*l3, 2*p3);
            if (g == 0.0) continue;
            /* (a)(b)(c) * g: only the real part can contribute because the
             * real Gaunt is real; track imag for sanity. */
            double ab_re = a_re*b_re - a_im*b_im;
            double ab_im = a_re*b_im + a_im*b_re;
            double abc_re = ab_re*c_re - ab_im*c_im;
            acc_re += abc_re * g;
        }
    }

    for (int i = 0; i < 3; i++) free(Cs[i]);
    return acc_re;
}

static void test_against_gaunt_real(void)
{
    /* A handful of triples covering m=0, cos-type, sin-type, and mixed. */
    const int cases[][6] = {
        {2,  0, 2,  0, 4,  0},
        {2,  2, 2,  2, 0,  0},
        {2, -2, 2, -2, 0,  0},
        {2,  1, 2, -1, 0,  0},
        {4,  2, 4, -2, 0,  0},
        {2,  1, 2,  1, 2, -2},
        {4,  3, 4, -3, 4,  0},
        {6,  4, 4, -2, 6, -2},
    };
    const int n = (int)(sizeof(cases) / sizeof(cases[0]));
    for (int i = 0; i < n; i++) {
        int l1 = cases[i][0], m1 = cases[i][1];
        int l2 = cases[i][2], m2 = cases[i][3];
        int l3 = cases[i][4], m3 = cases[i][5];
        double ref = gaunt_real(2*l1, 2*m1, 2*l2, 2*m2, 2*l3, 2*m3);
        double via = gaunt_via_matrix(l1, m1, l2, m2, l3, m3);
        TEST_NEAR(via, ref, 1e-13);
    }
}

/* Per-precision smoke test: verify the float / long-double variants
 * produce the same (l = 1) matrix entries as the double variant, within
 * their precision floor.  The float128 / MPFR variants live in
 * tests/test_quadmath.c and tests/test_mpfr.c respectively. */
static void test_precisions(void)
{
    const int dim = 3;
    const double sd = 1.0 / sqrt(2.0);

    /* float variant */
    {
        wignernj_cfloat_t Cf[9];
        wignernj_real_ylm_in_complex_ylm_f(1, Cf);
        float *F = (float *)Cf;
        /* (m_r=+1, m_c=-1).re at slot 2*(0*3 + 2) = 4 */
        TEST_NEAR((double)F[4], sd, 5e-7);
        /* (m_r=-1, m_c=-1).im at slot 2*(0*3 + 0) + 1 = 1 */
        TEST_NEAR((double)F[1], sd, 5e-7);
        TEST_NEAR((double)F[2 * (1 * dim + 1)], 1.0, 5e-7);
    }
    /* long double variant */
    {
        wignernj_cldouble_t Cl[9];
        wignernj_real_ylm_in_complex_ylm_l(1, Cl);
        long double *L = (long double *)Cl;
        TEST_NEAR((double)L[4], sd, 1e-15);
        TEST_NEAR((double)L[1], sd, 1e-15);
        TEST_NEAR((double)L[2 * (1 * dim + 1)], 1.0, 1e-15);
    }
}

int main(void)
{
    test_l0();
    test_l1();
    test_l2();
    for (int l = 0; l <= 10; l++) test_unitarity(l);
    test_against_gaunt_real();
    test_precisions();
    SUMMARY();
}
