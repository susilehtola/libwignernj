/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "run_tests.h"
#include "../include/wignernj_mpfr.h"
#include <mpfr.h>

/*
 * Three test categories:
 *
 * 1. CORRECTNESS AT 64 BITS  -- mpfr result converted to double should
 *    match wigner*() double to within 1e-13 relative error.
 *
 * 2. HIGH-PRECISION ACCURACY -- for a symbol with a simple algebraic value,
 *    verify that at 256 bits the result agrees with the closed form to within
 *    2 ULP at 256-bit precision.  Example: (j j 0; m -m 0) = (-1)^(j-m)/sqrt(2j+1)
 *
 * 3. ZEROS                   -- selection-rule violations must produce mpfr_zero_p.
 */

/* ── helpers ─────────────────────────────────────────────────────────────── */

static double mpfr_to_d(mpfr_t v) { return mpfr_get_d(v, MPFR_RNDN); }

/* Return |a - b| / 2^(exp-1) where exp is the binary exponent of rop. */
static int mpfr_agrees_to_bits(mpfr_t a, mpfr_t b, int bits)
{
    mpfr_t diff;
    int ok;
    mpfr_init2(diff, mpfr_get_prec(a));
    mpfr_sub(diff, a, b, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);
    /* Allow 2 ULP at the given precision */
    mpfr_mul_2si(b, b, -(bits - 1), MPFR_RNDN);
    mpfr_abs(b, b, MPFR_RNDN);
    mpfr_mul_ui(b, b, 2, MPFR_RNDN);
    ok = (mpfr_cmp(diff, b) <= 0);
    mpfr_clear(diff);
    return ok;
}

/* ── 1. Correctness at 64 bits ───────────────────────────────────────────── */

static void test_3j_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    wigner3j_mpfr(v, 2,2,0,  0, 0,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner3j(2,2,0, 0,0,0), 1e-13);

    wigner3j_mpfr(v, 6,4,4,  2,-2,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner3j(6,4,4, 2,-2,0), 1e-13);

    wigner3j_mpfr(v, 10,6,8, -8,4,4, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner3j(10,6,8,-8,4,4), 1e-13);

    mpfr_clear(v);
}

static void test_6j_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    wigner6j_mpfr(v, 2,2,2,  2,2,2, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner6j(2,2,2, 2,2,2), 1e-13);

    wigner6j_mpfr(v, 4,6,4,  6,4,6, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner6j(4,6,4, 6,4,6), 1e-13);

    wigner6j_mpfr(v, 10,8,6, 4,6,8, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner6j(10,8,6, 4,6,8), 1e-13);

    mpfr_clear(v);
}

static void test_9j_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    wigner9j_mpfr(v, 2,2,2, 2,2,2, 2,2,2, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner9j(2,2,2, 2,2,2, 2,2,2), 1e-13);

    wigner9j_mpfr(v, 2,4,4, 4,2,4, 4,4,2, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), wigner9j(2,4,4, 4,2,4, 4,4,2), 1e-13);

    mpfr_clear(v);
}

static void test_cg_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    clebsch_gordan_mpfr(v, 2,2, 2,-2, 4,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), clebsch_gordan(2,2, 2,-2, 4,0), 1e-13);

    clebsch_gordan_mpfr(v, 1,1, 1,-1, 2,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), clebsch_gordan(1,1, 1,-1, 2,0), 1e-13);

    mpfr_clear(v);
}

static void test_racah_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    racah_w_mpfr(v, 2,2,4, 2, 4,4, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), racah_w(2,2,4, 2, 4,4), 1e-13);

    racah_w_mpfr(v, 4,4,4, 4, 4,4, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), racah_w(4,4,4, 4, 4,4), 1e-13);

    mpfr_clear(v);
}

static void test_gaunt_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    gaunt_mpfr(v, 2,0, 2,0, 2,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), gaunt(2,0, 2,0, 2,0), 1e-13);

    gaunt_mpfr(v, 4,2, 4,-2, 4,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), gaunt(4,2, 4,-2, 4,0), 1e-13);

    mpfr_clear(v);
}

static void test_gaunt_real_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    /* All-m=0: real Gaunt equals complex Gaunt. */
    gaunt_real_mpfr(v, 2,0, 2,0, 4,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), gaunt_real(2,0, 2,0, 4,0), 1e-13);

    /* Generic non-zero real-Gaunt case. */
    gaunt_real_mpfr(v, 4,2, 4,-2, 4,0, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), gaunt_real(4,2, 4,-2, 4,0), 1e-13);

    mpfr_clear(v);
}

static void test_fano_x_64bit(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    /* X(j1,j2,j12;j3,j4,j34;j13,j24,J)
     *   = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * 9j */
    fano_x_mpfr(v, 2,2,2, 2,2,2, 2,2,2, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), fano_x(2,2,2, 2,2,2, 2,2,2), 1e-13);

    fano_x_mpfr(v, 4,4,4, 4,4,4, 4,4,4, MPFR_RNDN);
    TEST_NEAR(mpfr_to_d(v), fano_x(4,4,4, 4,4,4, 4,4,4), 1e-13);

    /* Selection-rule zero (one triangle violated) */
    fano_x_mpfr(v, 2,2,2, 2,2,2, 2,2,5, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    mpfr_clear(v);
}

/* ── 2. High-precision accuracy ─────────────────────────────────────────── */

/*
 * (j j 0; m -m 0) = (-1)^(j-m) / sqrt(2j+1).
 * Use j=2 (tj=4), m=1 (tm=2): value = (-1)^1 / sqrt(5) = -1/sqrt(5).
 * Verify at 256 bits that the result agrees with -1/sqrt(5) to within 2 ULP.
 */
static void test_high_precision_3j(void)
{
    mpfr_t got, ref;
    mpfr_init2(got, 256);
    mpfr_init2(ref, 256);

    wigner3j_mpfr(got, 4,4,0, 2,-2,0, MPFR_RNDN);

    /* ref = -1/sqrt(5) */
    mpfr_set_ui(ref, 5, MPFR_RNDN);
    mpfr_sqrt(ref, ref, MPFR_RNDN);
    mpfr_ui_div(ref, 1, ref, MPFR_RNDN);
    mpfr_neg(ref, ref, MPFR_RNDN);

    TEST_ASSERT(mpfr_agrees_to_bits(got, ref, 256));

    mpfr_clear(got);
    mpfr_clear(ref);
}

/*
 * (j j 0; m -m 0) = (-1)^(j-m) / sqrt(2j+1).
 * Use j=5 (tj=10), m=2 (tm=4): value = (-1)^3 / sqrt(11) = -1/sqrt(11).
 */
static void test_high_precision_3j_2(void)
{
    mpfr_t got, ref;
    mpfr_init2(got, 256);
    mpfr_init2(ref, 256);

    wigner3j_mpfr(got, 10,10,0, 4,-4,0, MPFR_RNDN);

    mpfr_set_ui(ref, 11, MPFR_RNDN);
    mpfr_sqrt(ref, ref, MPFR_RNDN);
    mpfr_ui_div(ref, 1, ref, MPFR_RNDN);
    mpfr_neg(ref, ref, MPFR_RNDN);

    TEST_ASSERT(mpfr_agrees_to_bits(got, ref, 256));

    mpfr_clear(got);
    mpfr_clear(ref);
}

/*
 * 6j orthogonality: {1 1 j12; 1 1 j12} summed over j12=0,1,2 with weight
 * (2*j12+1) equals 1 (if the formula is right this is just checking the sum).
 * Instead, use a simple exact value: {1 1 0; 1 1 0} = 1/3 exactly.
 * (This is a well-known entry in 6j tables.)
 *
 * Actually {j1 j2 0; j2 j1 0} = (-1)^(j1+j2) / sqrt[(2j1+1)(2j2+1)]
 * For j1=j2=1: (-1)^2 / sqrt(9) = 1/3.
 */
static void test_high_precision_6j(void)
{
    mpfr_t got, ref;
    mpfr_init2(got, 256);
    mpfr_init2(ref, 256);

    wigner6j_mpfr(got, 2,2,0, 2,2,0, MPFR_RNDN);

    mpfr_set_ui(ref, 3, MPFR_RNDN);
    mpfr_ui_div(ref, 1, ref, MPFR_RNDN);

    TEST_ASSERT(mpfr_agrees_to_bits(got, ref, 256));

    mpfr_clear(got);
    mpfr_clear(ref);
}

/* ── 3. Selection-rule zeros ─────────────────────────────────────────────── */

static void test_zeros(void)
{
    mpfr_t v;
    mpfr_init2(v, 64);

    /* 3j: m1+m2+m3 != 0 */
    wigner3j_mpfr(v, 2,2,2, 2,2,0, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    /* 3j: triangle rule violated */
    wigner3j_mpfr(v, 2,2,8, 0,0,0, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    /* 6j: triangle rule violated */
    wigner6j_mpfr(v, 2,2,8, 2,2,2, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    /* 9j: triangle rule violated */
    wigner9j_mpfr(v, 2,2,8, 2,2,2, 2,2,2, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    /* CG: m conservation violated */
    clebsch_gordan_mpfr(v, 2,2, 2,2, 4,0, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    /* Gaunt: m1+m2+m3 != 0 */
    gaunt_mpfr(v, 2,2, 2,2, 2,0, MPFR_RNDN);
    TEST_ASSERT(mpfr_zero_p(v));

    mpfr_clear(v);
}

/* ── main ────────────────────────────────────────────────────────────────── */

int main(void)
{
    test_3j_64bit();
    test_6j_64bit();
    test_9j_64bit();
    test_cg_64bit();
    test_racah_64bit();
    test_gaunt_64bit();
    test_gaunt_real_64bit();
    test_fano_x_64bit();
    test_high_precision_3j();
    test_high_precision_3j_2();
    test_high_precision_6j();
    test_zeros();
    SUMMARY();
}
