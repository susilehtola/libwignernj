/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Gaunt coefficients: integral of three complex spherical harmonics.
 *
 * Definition:
 *   G(l1,m1,l2,m2,l3,m3)
 *       = integral Y_{l1}^{m1}(Ω) Y_{l2}^{m2}(Ω) Y_{l3}^{m3}(Ω) dΩ
 *       = sqrt[(2l1+1)(2l2+1)(2l3+1) / (4π)]
 *         * (l1 l2 l3; 0 0 0) * (l1 l2 l3; m1 m2 m3)
 *
 * All l arguments must be non-negative integers; m must satisfy
 * m1+m2+m3=0 and |mi|<=li.  Arguments passed as tl=2*l, tm=2*m.
 *
 * Exact arithmetic path:
 *   The product of the two 3j symbols has the form
 *     (-1)^m3 * sqrt[outer_combined] * R0 * Rm
 *   where outer_combined = [Δ(l1,l2,l3)]² * (l1!)²(l2!)²(l3!)²
 *                          * (l1+m1)!(l1-m1)!(l2+m2)!(l2-m2)!(l3+m3)!(l3-m3)!
 *                          * (2l1+1)(2l2+1)(2l3+1) / 4
 *   and R0, Rm are the Racah integer sums for m=(0,0,0) and m=(m1,m2,m3).
 *   All factorial/integer factors are prime-factored exactly; only 1/sqrt(π)
 *   remains as a floating-point factor at the final conversion step.
 *
 * Phase derivation:
 *   phase_0 = (-1)^(l1-l2),   phase_m = (-1)^(l1-l2-m3)
 *   combined = (-1)^(2(l1-l2)-m3) = (-1)^(-m3) = (-1)^m3.
 */
#include "wigner_exact.h"
#include "pfrac.h"
#include "primes.h"
#include "wigner.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/* ── helpers ─────────────────────────────────────────────────────────────── */

static int gaunt_selection_rules(int tl1, int tm1, int tl2, int tm2,
                                  int tl3, int tm3)
{
    if (tm1 + tm2 + tm3 != 0)                          return 0;
    if (abs(tm1) > tl1 || abs(tm2) > tl2 || abs(tm3) > tl3) return 0;
    if ((tl1 + tl2 + tl3) & 1)                         return 0;
    if (tl3 < abs(tl1 - tl2) || tl3 > tl1 + tl2)      return 0;
    return 1;
}

/*
 * Summation bounds for the 3j Racah sum (same logic as bounds_3j in wigner3j.c).
 */
static void gaunt_racah_bounds(int tj1, int tj2, int tj3, int tm1, int tm2,
                                int *s_min, int *s_max)
{
    int a = (tj2 - tj3 - tm1) / 2;
    int b = (tj1 + tm2 - tj3) / 2;
    int c;
    *s_min = (a > b) ? a : b;
    if (*s_min < 0) *s_min = 0;
    a = (tj1 + tj2 - tj3) / 2;
    b = (tj1 - tm1) / 2;
    c = (tj2 + tm2) / 2;
    *s_max = a;
    if (b < *s_max) *s_max = b;
    if (c < *s_max) *s_max = c;
}

/*
 * Racah integer sum for 3j(j1,j2,j3; m1,m2,m3), no outer sqrt factor.
 * Fills sum_out (LCM-scaled magnitude), sum_sign_out, and lcm_exp[g_nprimes].
 * lcm_exp must be zeroed by the caller.
 */
static void gaunt_3j_racah_sum(int tj1, int tj2, int tj3,
                                int tm1, int tm2,
                                bigint_t *sum_out, int *sum_sign_out, int *lcm_exp)
{
    int s, s_min, s_max, pi;
    pfrac_t term;
    bigint_t sum_pos, sum_neg, scaled;

    gaunt_racah_bounds(tj1, tj2, tj3, tm1, tm2, &s_min, &s_max);

    if (s_min > s_max) {
        bigint_set_zero(sum_out);
        *sum_sign_out = 1;
        return;
    }

    pfrac_init(&term);
    bigint_init(&sum_pos);
    bigint_init(&sum_neg);
    bigint_init(&scaled);

    /* Pass 1: LCM exponents */
    for (s = s_min; s <= s_max; s++) {
        pfrac_zero(&term);
        pfrac_mul_factorial(&term, s);
        pfrac_mul_factorial(&term, (tj1 + tj2 - tj3) / 2 - s);
        pfrac_mul_factorial(&term, (tj1 - tm1) / 2 - s);
        pfrac_mul_factorial(&term, (tj2 + tm2) / 2 - s);
        pfrac_mul_factorial(&term, (tj3 - tj2 + tm1) / 2 + s);
        pfrac_mul_factorial(&term, (tj3 - tj1 - tm2) / 2 + s);
        for (pi = 0; pi < g_nprimes; pi++) {
            if (term.exp[pi] > lcm_exp[pi]) lcm_exp[pi] = term.exp[pi];
        }
    }

    /* Pass 2: accumulate scaled terms */
    for (s = s_min; s <= s_max; s++) {
        pfrac_zero(&term);
        pfrac_mul_factorial(&term, s);
        pfrac_mul_factorial(&term, (tj1 + tj2 - tj3) / 2 - s);
        pfrac_mul_factorial(&term, (tj1 - tm1) / 2 - s);
        pfrac_mul_factorial(&term, (tj2 + tm2) / 2 - s);
        pfrac_mul_factorial(&term, (tj3 - tj2 + tm1) / 2 + s);
        pfrac_mul_factorial(&term, (tj3 - tj1 - tm2) / 2 + s);
        bigint_set_u64(&scaled, 1);
        for (pi = 0; pi < g_nprimes; pi++) {
            int diff = lcm_exp[pi] - term.exp[pi];
            if (diff > 0)
                bigint_mul_prime_pow(&scaled, (uint64_t)g_primes[pi], diff);
        }
        if ((s & 1) == 0) bigint_add(&sum_pos, &sum_pos, &scaled);
        else              bigint_add(&sum_neg, &sum_neg, &scaled);
    }

    *sum_sign_out = bigint_sub_signed(sum_out, &sum_pos, &sum_neg);

    pfrac_free(&term);
    bigint_free(&sum_pos);
    bigint_free(&sum_neg);
    bigint_free(&scaled);
}

/* ── exact Gaunt (without the 1/sqrt(π) factor) ─────────────────────────── */

/*
 * Fills *out such that G(l1,m1,l2,m2,l3,m3)
 *   = wigner_exact_to_T(out) / sqrt(π).
 *
 * When is_zero is set, the Gaunt coefficient is zero and no further work
 * need be done.
 */
static void gaunt_exact(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3,
                         wigner_exact_t *out)
{
    int pi;
    pfrac_t outer;
    bigint_t sum0, summ;
    int ss0, ssm;
    int *lcm0, *lcmm;

    wigner_exact_init(out);

    if (!gaunt_selection_rules(tl1, tm1, tl2, tm2, tl3, tm3)) {
        out->is_zero = 1;
        return;
    }

    /* Phase = (-1)^m3 */
    out->sign = (((tm3 / 2) & 1) == 0) ? 1 : -1;

    /*
     * Combined outer pfrac (argument under the outer sqrt, without π):
     *
     *   [Δ(l1,l2,l3)]²:
     *     each factor contributes (l1+l2-l3)!(l1-l2+l3)!(-l1+l2+l3)! / (l1+l2+l3+1)!
     *     called twice → all exponents even → folds entirely into int_num/int_den
     *
     *   m=0 factor:  (l1!)²(l2!)²(l3!)²
     *
     *   general-m factor: (l1+m1)!(l1-m1)!(l2+m2)!(l2-m2)!(l3+m3)!(l3-m3)!
     *
     *   normalization: (2l1+1)(2l2+1)(2l3+1)/4
     */
    pfrac_init(&outer);

    /* Δ(l1,l2,l3) — first copy */
    pfrac_mul_factorial(&outer, (tl1 + tl2 - tl3) / 2);
    pfrac_mul_factorial(&outer, (tl1 - tl2 + tl3) / 2);
    pfrac_mul_factorial(&outer, (-tl1 + tl2 + tl3) / 2);
    pfrac_div_factorial(&outer, (tl1 + tl2 + tl3) / 2 + 1);
    /* Δ(l1,l2,l3) — second copy (makes [Δ]² rational) */
    pfrac_mul_factorial(&outer, (tl1 + tl2 - tl3) / 2);
    pfrac_mul_factorial(&outer, (tl1 - tl2 + tl3) / 2);
    pfrac_mul_factorial(&outer, (-tl1 + tl2 + tl3) / 2);
    pfrac_div_factorial(&outer, (tl1 + tl2 + tl3) / 2 + 1);

    /* (l1!)²(l2!)²(l3!)² from the m=0 3j */
    pfrac_mul_factorial(&outer, tl1 / 2); pfrac_mul_factorial(&outer, tl1 / 2);
    pfrac_mul_factorial(&outer, tl2 / 2); pfrac_mul_factorial(&outer, tl2 / 2);
    pfrac_mul_factorial(&outer, tl3 / 2); pfrac_mul_factorial(&outer, tl3 / 2);

    /* (l1+m1)!(l1-m1)!(l2+m2)!(l2-m2)!(l3+m3)!(l3-m3)! from the general-m 3j */
    pfrac_mul_factorial(&outer, (tl1 + tm1) / 2);
    pfrac_mul_factorial(&outer, (tl1 - tm1) / 2);
    pfrac_mul_factorial(&outer, (tl2 + tm2) / 2);
    pfrac_mul_factorial(&outer, (tl2 - tm2) / 2);
    pfrac_mul_factorial(&outer, (tl3 + tm3) / 2);
    pfrac_mul_factorial(&outer, (tl3 - tm3) / 2);

    /* (2l1+1)(2l2+1)(2l3+1) / 4 : integer factors folded in via prime factorization.
     * g_primes[0] = 2, so subtracting 2 from exp[0] divides by 4 = 2². */
    pfrac_mul_int(&outer, tl1 + 1);
    pfrac_mul_int(&outer, tl2 + 1);
    pfrac_mul_int(&outer, tl3 + 1);
    outer.exp[0] -= 2;

    /* Racah integer sums for both 3j symbols */
    lcm0 = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
    lcmm = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
    bigint_init(&sum0);
    bigint_init(&summ);

    gaunt_3j_racah_sum(tl1, tl2, tl3, 0, 0, &sum0, &ss0, lcm0);
    if (bigint_is_zero(&sum0)) { out->is_zero = 1; goto cleanup; }

    gaunt_3j_racah_sum(tl1, tl2, tl3, tm1, tm2, &summ, &ssm, lcmm);
    if (bigint_is_zero(&summ)) { out->is_zero = 1; goto cleanup; }

    /* Product of the two Racah sums */
    bigint_mul(&out->sum, &sum0, &summ);
    out->sum_sign = ss0 * ssm;

    /* Convert outer pfrac */
    bigint_set_u64(&out->int_num,  1);
    bigint_set_u64(&out->int_den,  1);
    bigint_set_u64(&out->sqrt_num, 1);
    bigint_set_u64(&out->sqrt_den, 1);
    pfrac_to_sqrt_rational(&outer, &out->int_num, &out->int_den,
                           &out->sqrt_num, &out->sqrt_den);

    /* Combined LCM denominator → int_den */
    for (pi = 0; pi < g_nprimes; pi++) {
        int e = lcm0[pi] + lcmm[pi];
        if (e > 0)
            bigint_mul_prime_pow(&out->int_den, (uint64_t)g_primes[pi], e);
    }

cleanup:
    free(lcm0);
    free(lcmm);
    pfrac_free(&outer);
    bigint_free(&sum0);
    bigint_free(&summ);
}

/* ── public API ──────────────────────────────────────────────────────────── */

float gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3)
{
    wigner_exact_t e;
    float result;
    gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &e);
    result = wigner_exact_to_float(&e) / sqrtf((float)M_PI);
    wigner_exact_free(&e);
    return result;
}

double gaunt(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3)
{
    wigner_exact_t e;
    double result;
    gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &e);
    result = wigner_exact_to_double(&e) / sqrt(M_PI);
    wigner_exact_free(&e);
    return result;
}

long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3)
{
    wigner_exact_t e;
    long double result;
    gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &e);
    result = wigner_exact_to_long_double(&e) / sqrtl(acosl(-1.0L));
    wigner_exact_free(&e);
    return result;
}

#ifdef WIGNER_HAVE_MPFR
#include "wigner_mpfr.h"
void gaunt_mpfr(mpfr_t rop, int tl1, int tm1, int tl2, int tm2,
                             int tl3, int tm3, mpfr_rnd_t rnd)
{
    wigner_exact_t e;
    mpfr_t pi;

    gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &e);
    wigner_exact_to_mpfr(rop, &e, rnd);
    wigner_exact_free(&e);

    if (!mpfr_zero_p(rop)) {
        mpfr_init2(pi, mpfr_get_prec(rop));
        mpfr_const_pi(pi, rnd);
        mpfr_sqrt(pi, pi, rnd);
        mpfr_div(rop, rop, pi, rnd);
        mpfr_clear(pi);
    }
}
#endif
