/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Wigner 6j symbol via the Racah formula with prime-factorization exact
 * arithmetic following Johansson & Forssén (doi:10.1137/15M1021908).
 *
 * Racah formula:
 *   {j1 j2 j3}  = Delta(j1,j2,j3)*Delta(j1,j5,j6)*Delta(j4,j2,j6)*Delta(j4,j5,j3)
 *   {j4 j5 j6}
 *     * sum_s (-1)^s (s+1)! / [(s-j1-j2-j3)!(s-j1-j5-j6)!(s-j4-j2-j6)!
 *                               (s-j4-j5-j3)!(j1+j2+j4+j5-s)!
 *                               (j2+j3+j5+j6-s)!(j1+j3+j4+j6-s)!]
 */
#include "wigner_exact.h"
#include "pfrac.h"
#include "primes.h"
#include "wigner.h"
#include "xalloc.h"
#include <stdlib.h>

static int triangle_ok(int ta, int tb, int tc)
{
    if ((ta + tb + tc) & 1) return 0;
    if (tc < abs(ta - tb)) return 0;
    if (tc > ta + tb)      return 0;
    return 1;
}

static int selection_rules_6j(int tj1, int tj2, int tj3,
                               int tj4, int tj5, int tj6)
{
    return triangle_ok(tj1, tj2, tj3) &&
           triangle_ok(tj1, tj5, tj6) &&
           triangle_ok(tj4, tj2, tj6) &&
           triangle_ok(tj4, tj5, tj3);
}

/* Summation bounds for the 6j Racah sum.
 * s_min = max(j1+j2+j3, j1+j5+j6, j4+j2+j6, j4+j5+j3)
 * s_max = min(j1+j2+j4+j5, j2+j3+j5+j6, j1+j3+j4+j6)
 * All in half-integer 2x representation.
 */
static void bounds_6j(int tj1, int tj2, int tj3,
                      int tj4, int tj5, int tj6,
                      int *s_min, int *s_max)
{
    int a, b, c, d;
    a = (tj1 + tj2 + tj3) / 2;
    b = (tj1 + tj5 + tj6) / 2;
    c = (tj4 + tj2 + tj6) / 2;
    d = (tj4 + tj5 + tj3) / 2;
    *s_min = a;
    if (b > *s_min) *s_min = b;
    if (c > *s_min) *s_min = c;
    if (d > *s_min) *s_min = d;

    a = (tj1 + tj2 + tj4 + tj5) / 2;
    b = (tj2 + tj3 + tj5 + tj6) / 2;
    c = (tj1 + tj3 + tj4 + tj6) / 2;
    *s_max = a;
    if (b < *s_max) *s_max = b;
    if (c < *s_max) *s_max = c;
}

/*
 * Add all four Delta triangle-coefficient sqrt factors to outer pfrac.
 * Each Delta(a,b,c) contributes:
 *   numerator:   (a+b-c)! (a-b+c)! (-a+b+c)!
 *   denominator: (a+b+c+1)!
 * (all divided by 2 since ta,tb,tc are 2x values)
 */
static void add_delta_sqrt(pfrac_t *outer, int ta, int tb, int tc)
{
    pfrac_mul_factorial(outer, (ta + tb - tc) / 2);
    pfrac_mul_factorial(outer, (ta - tb + tc) / 2);
    pfrac_mul_factorial(outer, (-ta + tb + tc) / 2);
    pfrac_div_factorial(outer, (ta + tb + tc) / 2 + 1);
}

void wigner6j_exact(int tj1, int tj2, int tj3,
                    int tj4, int tj5, int tj6,
                    wigner_exact_t *out)
{
    int s, s_min, s_max, pi;
    int *lcm_exp;
    pfrac_t outer, term;
    bigint_t sum_pos, sum_neg, scaled;
    bigint_ws_t ws;
    size_t mw;

    wigner_exact_init(out);

    if (!selection_rules_6j(tj1, tj2, tj3, tj4, tj5, tj6)) {
        out->is_zero = 1;
        return;
    }

    bounds_6j(tj1, tj2, tj3, tj4, tj5, tj6, &s_min, &s_max);
    if (s_min > s_max) { out->is_zero = 1; return; }

    out->sign = 1; /* phase is carried by (-1)^s in the sum */

    /* Largest factorial argument: max β_u + 1 (the (s+1)! numerator).  Use a
     * safe upper bound on N for sizing the workspace and bigints. */
    mw = bigint_words_for_factorial(
            (tj1 + tj2 + tj3 + tj4 + tj5 + tj6) / 2 + 5);

    bigint_ws_init(&ws);
    bigint_ws_reserve(&ws, mw);

    /* Outer sqrt = product of 4 Delta triangle coefficients */
    pfrac_init(&outer);
    add_delta_sqrt(&outer, tj1, tj2, tj3);
    add_delta_sqrt(&outer, tj1, tj5, tj6);
    add_delta_sqrt(&outer, tj4, tj2, tj6);
    add_delta_sqrt(&outer, tj4, tj5, tj3);

    /* ── Pass 1: find LCM exponents ── */
    lcm_exp = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
    pfrac_init(&term);

    for (s = s_min; s <= s_max; s++) {
        /* Denominator factorials of term s:
         * (s-j1-j2-j3)! (s-j1-j5-j6)! (s-j4-j2-j6)! (s-j4-j5-j3)!
         * (j1+j2+j4+j5-s)! (j2+j3+j5+j6-s)! (j1+j3+j4+j6-s)!
         * Numerator factorial: (s+1)!
         */
        pfrac_zero(&term);
        pfrac_div_factorial(&term, s - (tj1 + tj2 + tj3) / 2);
        pfrac_div_factorial(&term, s - (tj1 + tj5 + tj6) / 2);
        pfrac_div_factorial(&term, s - (tj4 + tj2 + tj6) / 2);
        pfrac_div_factorial(&term, s - (tj4 + tj5 + tj3) / 2);
        pfrac_div_factorial(&term, (tj1 + tj2 + tj4 + tj5) / 2 - s);
        pfrac_div_factorial(&term, (tj2 + tj3 + tj5 + tj6) / 2 - s);
        pfrac_div_factorial(&term, (tj1 + tj3 + tj4 + tj6) / 2 - s);
        pfrac_mul_factorial(&term, s + 1);

        /* term.exp[pi] is the net exponent; for LCM we want -min(term.exp[pi]) */
        for (pi = 0; pi < g_nprimes; pi++) {
            int neg = -term.exp[pi]; /* positive = denominator contribution */
            if (neg > lcm_exp[pi]) lcm_exp[pi] = neg;
        }
    }

    /* ── Pass 2: accumulate sum (all bigints pre-sized to mw) ── */
    bigint_init(&sum_pos); bigint_reserve(&sum_pos, mw);
    bigint_init(&sum_neg); bigint_reserve(&sum_neg, mw);
    bigint_init(&scaled);  bigint_reserve(&scaled,  mw);

    for (s = s_min; s <= s_max; s++) {
        pfrac_zero(&term);
        pfrac_div_factorial(&term, s - (tj1 + tj2 + tj3) / 2);
        pfrac_div_factorial(&term, s - (tj1 + tj5 + tj6) / 2);
        pfrac_div_factorial(&term, s - (tj4 + tj2 + tj6) / 2);
        pfrac_div_factorial(&term, s - (tj4 + tj5 + tj3) / 2);
        pfrac_div_factorial(&term, (tj1 + tj2 + tj4 + tj5) / 2 - s);
        pfrac_div_factorial(&term, (tj2 + tj3 + tj5 + tj6) / 2 - s);
        pfrac_div_factorial(&term, (tj1 + tj3 + tj4 + tj6) / 2 - s);
        pfrac_mul_factorial(&term, s + 1);

        /* scaled = prod p_i^(lcm_exp[i] + term.exp[i])  [all non-negative] */
        bigint_set_u64(&scaled, 1);
        for (pi = 0; pi < g_nprimes; pi++) {
            int e = lcm_exp[pi] + term.exp[pi];
            if (e > 0)
                bigint_mul_prime_pow_ws(&scaled, (uint64_t)g_primes[pi], e, &ws);
        }

        if ((s & 1) == 0)
            bigint_add(&sum_pos, &sum_pos, &scaled);
        else
            bigint_add(&sum_neg, &sum_neg, &scaled);
    }

    {
        int sign_diff = bigint_sub_signed(&out->sum, &sum_pos, &sum_neg);
        out->sum_sign = sign_diff;
    }

    /* Build output bigints (also pre-sized) */
    bigint_set_u64(&out->int_num,  1);
    bigint_set_u64(&out->int_den,  1);
    bigint_set_u64(&out->sqrt_num, 1);
    bigint_set_u64(&out->sqrt_den, 1);
    bigint_reserve(&out->int_num,  mw);
    bigint_reserve(&out->int_den,  mw);
    bigint_reserve(&out->sqrt_num, mw);
    bigint_reserve(&out->sqrt_den, mw);

    pfrac_to_sqrt_rational_ws(&outer,
                               &out->int_num, &out->int_den,
                               &out->sqrt_num, &out->sqrt_den, &ws);

    /* LCM denominator → int_den */
    for (pi = 0; pi < g_nprimes; pi++) {
        if (lcm_exp[pi] > 0)
            bigint_mul_prime_pow_ws(&out->int_den,
                                    (uint64_t)g_primes[pi], lcm_exp[pi], &ws);
    }

    free(lcm_exp);
    pfrac_free(&outer);
    pfrac_free(&term);
    bigint_free(&sum_pos);
    bigint_free(&sum_neg);
    bigint_free(&scaled);
    bigint_ws_free(&ws);
}

/* ── public API ──────────────────────────────────────────────────────────── */

float wigner6j_f(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6)
{
    wigner_exact_t e;
    float result;
    wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &e);
    result = wigner_exact_to_float(&e);
    wigner_exact_free(&e);
    return result;
}

double wigner6j(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6)
{
    wigner_exact_t e;
    double result;
    wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &e);
    result = wigner_exact_to_double(&e);
    wigner_exact_free(&e);
    return result;
}

long double wigner6j_l(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6)
{
    wigner_exact_t e;
    long double result;
    wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &e);
    result = wigner_exact_to_long_double(&e);
    wigner_exact_free(&e);
    return result;
}

#ifdef WIGNER_HAVE_MPFR
void wigner6j_mpfr(mpfr_t rop, int tj1, int tj2, int tj3,
                               int tj4, int tj5, int tj6, mpfr_rnd_t rnd)
{
    wigner_exact_t e;
    wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &e);
    wigner_exact_to_mpfr(rop, &e, rnd);
    wigner_exact_free(&e);
}
#endif
