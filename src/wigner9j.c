/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Wigner 9j symbol via exact prime-factorization arithmetic.
 *
 * Formula (Racah W-coefficient representation):
 *   {j11 j12 j13}
 *   {j21 j22 j23} = (-1)^(tj13+tj22+tj31) * sum_k (2k+1)
 *   {j31 j32 j33}
 *       * {j11 j21 j31} * {j11 j12 j13} * {j22 j21 j23}
 *         {j32 j33  k }   {j23 j33  k }   { k  j12 j32}
 *
 * (Phase exponent is in 2j units; tj values are always integers so the phase
 * is always ±1.)
 *
 * Exactness argument: the product of the three 6j symbols at each k is
 *   sqrt(C) * rational(k) * R1(k) * R2(k) * R3(k)
 * where sqrt(C) is formed by the six k-independent Delta factors (two per 6j):
 *   Δ(j11,j21,j31), Δ(j32,j33,j31), Δ(j11,j12,j13), Δ(j23,j33,j13),
 *   Δ(j22,j21,j23), Δ(j22,j12,j32)
 * and rational(k) = Δ²(j11,j33,k)·Δ²(j21,j32,k)·Δ²(j12,j23,k) is purely
 * rational (each k-dependent Δ appears exactly twice across the three 6j
 * symbols, making it a perfect square and thus rational).  Ri(k) is the
 * Racah integer sum for the i-th 6j.  The 9j therefore has the form
 * sqrt(C) × exact_integer, enabling the same exact pipeline as 3j and 6j.
 *
 * k range: max(|j11-j33|,|j21-j32|,|j12-j23|) <= k <= min(j11+j33,j21+j32,j12+j23)
 */
#include "wigner_exact.h"
#include "pfrac.h"
#include "primes.h"
#include "wigner.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>

/* ── helpers ─────────────────────────────────────────────────────────────── */

static int triangle_ok(int ta, int tb, int tc)
{
    if ((ta + tb + tc) & 1) return 0;
    if (tc < abs(ta - tb))  return 0;
    if (tc > ta + tb)       return 0;
    return 1;
}

static int max3(int a, int b, int c)
{ int m = a; if (b > m) m = b; if (c > m) m = c; return m; }

static int min3(int a, int b, int c)
{ int m = a; if (b < m) m = b; if (c < m) m = c; return m; }

static int selection_rules_9j(int tj11, int tj12, int tj13,
                               int tj21, int tj22, int tj23,
                               int tj31, int tj32, int tj33)
{
    return triangle_ok(tj11, tj12, tj13) &&
           triangle_ok(tj21, tj22, tj23) &&
           triangle_ok(tj31, tj32, tj33) &&
           triangle_ok(tj11, tj21, tj31) &&
           triangle_ok(tj12, tj22, tj32) &&
           triangle_ok(tj13, tj23, tj33);
}

/* Δ(a,b,c) = sqrt[(a+b-c)!(a-b+c)!(-a+b+c)!/(a+b+c+1)!]; add to pfrac. */
static void add_delta_sqrt(pfrac_t *f, int ta, int tb, int tc)
{
    pfrac_mul_factorial(f, (ta + tb - tc) / 2);
    pfrac_mul_factorial(f, (ta - tb + tc) / 2);
    pfrac_mul_factorial(f, (-ta + tb + tc) / 2);
    pfrac_div_factorial(f, (ta + tb + tc) / 2 + 1);
}

/*
 * Compute the Racah integer sum for 6j {a b c; d e f} (2j-integer args).
 * Fills sum_out with the LCM-scaled magnitude, sum_sign_out with its sign,
 * and lcm_exp[0..g_nprimes-1] with the LCM denominator exponents.
 * lcm_exp must be a caller-allocated array of g_nprimes ints (need not be
 * zeroed; this function initialises it).
 */
/* Inner Racah sum for one of the three 6j symbols.  Caller-owned scratch
 * bigints (sum_pos, sum_neg, scaled) are reused across all calls of one 9j
 * evaluation, eliminating per-call bigint allocations. */
static void racah_6j_sum(int tj1, int tj2, int tj3,
                          int tj4, int tj5, int tj6,
                          bigint_t *sum_out, int *sum_sign_out, int *lcm_exp,
                          bigint_t *sum_pos, bigint_t *sum_neg,
                          bigint_t *scaled, bigint_ws_t *ws)
{
    int s, s_min, s_max, pi;
    pfrac_t term;

    {
        int a = (tj1+tj2+tj3)/2, b = (tj1+tj5+tj6)/2;
        int c = (tj4+tj2+tj6)/2, d = (tj4+tj5+tj3)/2;
        s_min = a;
        if (b > s_min) s_min = b;
        if (c > s_min) s_min = c;
        if (d > s_min) s_min = d;
    }
    {
        int a = (tj1+tj2+tj4+tj5)/2;
        int b = (tj2+tj3+tj5+tj6)/2;
        int c = (tj1+tj3+tj4+tj6)/2;
        s_max = a;
        if (b < s_max) s_max = b;
        if (c < s_max) s_max = c;
    }

    memset(lcm_exp, 0, (size_t)g_nprimes * sizeof(int));

    if (s_min > s_max) {
        bigint_set_zero(sum_out);
        *sum_sign_out = 1;
        return;
    }

    pfrac_init(&term);
    bigint_set_zero(sum_pos);
    bigint_set_zero(sum_neg);

    /* Pass 1: find LCM exponents */
    for (s = s_min; s <= s_max; s++) {
        pfrac_zero(&term);
        pfrac_div_factorial(&term, s - (tj1+tj2+tj3)/2);
        pfrac_div_factorial(&term, s - (tj1+tj5+tj6)/2);
        pfrac_div_factorial(&term, s - (tj4+tj2+tj6)/2);
        pfrac_div_factorial(&term, s - (tj4+tj5+tj3)/2);
        pfrac_div_factorial(&term, (tj1+tj2+tj4+tj5)/2 - s);
        pfrac_div_factorial(&term, (tj2+tj3+tj5+tj6)/2 - s);
        pfrac_div_factorial(&term, (tj1+tj3+tj4+tj6)/2 - s);
        pfrac_mul_factorial(&term, s + 1);
        for (pi = 0; pi < g_nprimes; pi++) {
            int neg = -term.exp[pi];
            if (neg > lcm_exp[pi]) lcm_exp[pi] = neg;
        }
    }

    /* Pass 2: accumulate */
    for (s = s_min; s <= s_max; s++) {
        pfrac_zero(&term);
        pfrac_div_factorial(&term, s - (tj1+tj2+tj3)/2);
        pfrac_div_factorial(&term, s - (tj1+tj5+tj6)/2);
        pfrac_div_factorial(&term, s - (tj4+tj2+tj6)/2);
        pfrac_div_factorial(&term, s - (tj4+tj5+tj3)/2);
        pfrac_div_factorial(&term, (tj1+tj2+tj4+tj5)/2 - s);
        pfrac_div_factorial(&term, (tj2+tj3+tj5+tj6)/2 - s);
        pfrac_div_factorial(&term, (tj1+tj3+tj4+tj6)/2 - s);
        pfrac_mul_factorial(&term, s + 1);
        bigint_set_u64(scaled, 1);
        for (pi = 0; pi < g_nprimes; pi++) {
            int e = lcm_exp[pi] + term.exp[pi];
            if (e > 0) bigint_mul_prime_pow_ws(scaled, (uint64_t)g_primes[pi], e, ws);
        }
        if ((s & 1) == 0) bigint_add(sum_pos, sum_pos, scaled);
        else              bigint_add(sum_neg, sum_neg, scaled);
    }

    *sum_sign_out = bigint_sub_signed(sum_out, sum_pos, sum_neg);

    pfrac_free(&term);
}

/* ── exact 9j ────────────────────────────────────────────────────────────── */

void wigner9j_exact(int tj11, int tj12, int tj13,
                    int tj21, int tj22, int tj23,
                    int tj31, int tj32, int tj33,
                    wigner_exact_t *out)
{
    int tk, tk_min, tk_max, pi;
    pfrac_t outer, kdep;
    bigint_t sum_pos, sum_neg, scaled, s1, s2, s3, tmp;
    int ss1, ss2, ss3;
    int *global_lcm, *lcm1, *lcm2, *lcm3;
    bigint_ws_t ws;
    size_t mw, mw_prod;

    wigner_exact_init(out);

    if (!selection_rules_9j(tj11,tj12,tj13,tj21,tj22,tj23,tj31,tj32,tj33)) {
        out->is_zero = 1;
        return;
    }

    tk_min = max3(abs(tj11 - tj33), abs(tj12 - tj23), abs(tj21 - tj32));
    tk_max = min3(tj11 + tj33,       tj12 + tj23,       tj21 + tj32);
    if (((tk_min + tj11 + tj33) & 1) != 0) tk_min++;
    if (((tk_max + tj11 + tj33) & 1) != 0) tk_max--;
    if (tk_min > tk_max) { out->is_zero = 1; return; }

    out->sign = (((tj13 + tj22 + tj31) & 1) == 0) ? 1 : -1;

    /* Largest factorial argument is bounded by max(β_u + 1) over all three
     * inner 6j calls, which in turn is at most (sum of all 9 tj plus 2*tk_max)/2
     * + 1.  We use a uniform safe upper bound for sizing. */
    mw = bigint_words_for_factorial(
            (tj11 + tj12 + tj13 + tj21 + tj22 + tj23
             + tj31 + tj32 + tj33 + 2 * tk_max) / 2 + 5);
    /* The 9j accumulator also holds products of three 6j-sized bigints,
     * so allocate up to ~5x that for the long-lived bigints. */
    mw_prod = 5 * mw;

    bigint_ws_init(&ws);
    bigint_ws_reserve(&ws, mw_prod);

    /*
     * Outer sqrt: product of the six k-independent Δ factors, two per 6j:
     *   6j_1 = {j11 j21 j31; j32 j33 k}: Δ(j11,j21,j31)  Δ(j32,j33,j31)
     *   6j_2 = {j11 j12 j13; j23 j33 k}: Δ(j11,j12,j13)  Δ(j23,j33,j13)
     *   6j_3 = {j22 j21 j23; k  j12 j32}: Δ(j22,j21,j23)  Δ(j22,j12,j32)
     */
    pfrac_init(&outer);
    add_delta_sqrt(&outer, tj11, tj21, tj31);
    add_delta_sqrt(&outer, tj32, tj33, tj31);
    add_delta_sqrt(&outer, tj11, tj12, tj13);
    add_delta_sqrt(&outer, tj23, tj33, tj13);
    add_delta_sqrt(&outer, tj22, tj21, tj23);
    add_delta_sqrt(&outer, tj22, tj12, tj32);

    global_lcm = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
    lcm1 = (int *)xmalloc((size_t)g_nprimes * sizeof(int));
    lcm2 = (int *)xmalloc((size_t)g_nprimes * sizeof(int));
    lcm3 = (int *)xmalloc((size_t)g_nprimes * sizeof(int));

    pfrac_init(&kdep);
    bigint_init(&s1); bigint_reserve(&s1, mw);
    bigint_init(&s2); bigint_reserve(&s2, mw);
    bigint_init(&s3); bigint_reserve(&s3, mw);
    bigint_init(&tmp);    bigint_reserve(&tmp,    mw_prod);
    bigint_init(&scaled); bigint_reserve(&scaled, mw_prod);
    bigint_init(&sum_pos); bigint_reserve(&sum_pos, mw_prod);
    bigint_init(&sum_neg); bigint_reserve(&sum_neg, mw_prod);
    /* Per-call scratch for the inner racah_6j_sum, hoisted out so it is
     * allocated once per 9j evaluation rather than once per (tk, sub-6j). */
    bigint_t inner_pos, inner_neg, inner_scaled;
    bigint_init(&inner_pos);    bigint_reserve(&inner_pos,    mw);
    bigint_init(&inner_neg);    bigint_reserve(&inner_neg,    mw);
    bigint_init(&inner_scaled); bigint_reserve(&inner_scaled, mw);

    /* ── Pass 1: global LCM ───────────────────────────────────────────────── */
    for (tk = tk_min; tk <= tk_max; tk += 2) {
        /*
         * k-dependent rational: Δ²(j11,j33,k)·Δ²(j21,j32,k)·Δ²(j12,j23,k).
         * Two add_delta_sqrt calls per pair → all pfrac exponents are even →
         * purely rational, no sqrt residual.
         */
        pfrac_zero(&kdep);
        add_delta_sqrt(&kdep, tj11, tj33, tk); add_delta_sqrt(&kdep, tj11, tj33, tk);
        add_delta_sqrt(&kdep, tj21, tj32, tk); add_delta_sqrt(&kdep, tj21, tj32, tk);
        add_delta_sqrt(&kdep, tj12, tj23, tk); add_delta_sqrt(&kdep, tj12, tj23, tk);

        racah_6j_sum(tj11,tj21,tj31, tj32,tj33,tk,  &tmp,&ss1,lcm1,
                     &inner_pos, &inner_neg, &inner_scaled, &ws);
        racah_6j_sum(tj11,tj12,tj13, tj23,tj33,tk,  &tmp,&ss2,lcm2,
                     &inner_pos, &inner_neg, &inner_scaled, &ws);
        racah_6j_sum(tj22,tj21,tj23, tk,tj12,tj32,  &tmp,&ss3,lcm3,
                     &inner_pos, &inner_neg, &inner_scaled, &ws);

        for (pi = 0; pi < g_nprimes; pi++) {
            /* kdep.exp[pi] is the sqrt-argument exponent; divide by 2 for value */
            int net = kdep.exp[pi] / 2 - lcm1[pi] - lcm2[pi] - lcm3[pi];
            if (-net > global_lcm[pi]) global_lcm[pi] = -net;
        }
    }

    /* ── Pass 2: accumulate ───────────────────────────────────────────────── */
    for (tk = tk_min; tk <= tk_max; tk += 2) {
        pfrac_zero(&kdep);
        add_delta_sqrt(&kdep, tj11, tj33, tk); add_delta_sqrt(&kdep, tj11, tj33, tk);
        add_delta_sqrt(&kdep, tj21, tj32, tk); add_delta_sqrt(&kdep, tj21, tj32, tk);
        add_delta_sqrt(&kdep, tj12, tj23, tk); add_delta_sqrt(&kdep, tj12, tj23, tk);

        racah_6j_sum(tj11,tj21,tj31, tj32,tj33,tk,  &s1,&ss1,lcm1,
                     &inner_pos, &inner_neg, &inner_scaled, &ws);
        racah_6j_sum(tj11,tj12,tj13, tj23,tj33,tk,  &s2,&ss2,lcm2,
                     &inner_pos, &inner_neg, &inner_scaled, &ws);
        racah_6j_sum(tj22,tj21,tj23, tk,tj12,tj32,  &s3,&ss3,lcm3,
                     &inner_pos, &inner_neg, &inner_scaled, &ws);

        if (bigint_is_zero(&s1) || bigint_is_zero(&s2) || bigint_is_zero(&s3))
            continue;

        /* scale[pi] = global_lcm[pi] + net[pi] >= 0 by construction */
        bigint_set_u64(&scaled, (uint64_t)(tk + 1)); /* factor (2k+1) */
        for (pi = 0; pi < g_nprimes; pi++) {
            int e = global_lcm[pi] + kdep.exp[pi] / 2 - lcm1[pi] - lcm2[pi] - lcm3[pi];
            if (e > 0) bigint_mul_prime_pow_ws(&scaled, (uint64_t)g_primes[pi], e, &ws);
        }

        bigint_mul_ws(&tmp,    &scaled, &s1, &ws);
        bigint_mul_ws(&scaled, &tmp,    &s2, &ws);
        bigint_mul_ws(&tmp,    &scaled, &s3, &ws);

        /* outer_phase absorbed into out->sign; accumulate by ss1*ss2*ss3 */
        if (ss1 * ss2 * ss3 > 0) bigint_add(&sum_pos, &sum_pos, &tmp);
        else                      bigint_add(&sum_neg, &sum_neg, &tmp);
    }

    out->sum_sign = bigint_sub_signed(&out->sum, &sum_pos, &sum_neg);

    bigint_set_u64(&out->int_num,  1); bigint_reserve(&out->int_num,  mw_prod);
    bigint_set_u64(&out->int_den,  1); bigint_reserve(&out->int_den,  mw_prod);
    bigint_set_u64(&out->sqrt_num, 1); bigint_reserve(&out->sqrt_num, mw_prod);
    bigint_set_u64(&out->sqrt_den, 1); bigint_reserve(&out->sqrt_den, mw_prod);
    pfrac_to_sqrt_rational_ws(&outer,
                              &out->int_num, &out->int_den,
                              &out->sqrt_num, &out->sqrt_den, &ws);
    for (pi = 0; pi < g_nprimes; pi++) {
        if (global_lcm[pi] > 0)
            bigint_mul_prime_pow_ws(&out->int_den,
                                    (uint64_t)g_primes[pi], global_lcm[pi], &ws);
    }

    free(global_lcm);
    free(lcm1); free(lcm2); free(lcm3);
    pfrac_free(&outer);
    pfrac_free(&kdep);
    bigint_free(&s1); bigint_free(&s2); bigint_free(&s3); bigint_free(&tmp);
    bigint_free(&scaled);
    bigint_free(&sum_pos);
    bigint_free(&sum_neg);
    bigint_free(&inner_pos);
    bigint_free(&inner_neg);
    bigint_free(&inner_scaled);
    bigint_ws_free(&ws);
}

/* ── public API ──────────────────────────────────────────────────────────── */

float wigner9j_f(int tj11, int tj12, int tj13,
                 int tj21, int tj22, int tj23,
                 int tj31, int tj32, int tj33)
{
    wigner_exact_t e;
    float result;
    wigner9j_exact(tj11,tj12,tj13, tj21,tj22,tj23, tj31,tj32,tj33, &e);
    result = wigner_exact_to_float(&e);
    wigner_exact_free(&e);
    return result;
}

double wigner9j(int tj11, int tj12, int tj13,
                int tj21, int tj22, int tj23,
                int tj31, int tj32, int tj33)
{
    wigner_exact_t e;
    double result;
    wigner9j_exact(tj11,tj12,tj13, tj21,tj22,tj23, tj31,tj32,tj33, &e);
    result = wigner_exact_to_double(&e);
    wigner_exact_free(&e);
    return result;
}

long double wigner9j_l(int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33)
{
    wigner_exact_t e;
    long double result;
    wigner9j_exact(tj11,tj12,tj13, tj21,tj22,tj23, tj31,tj32,tj33, &e);
    result = wigner_exact_to_long_double(&e);
    wigner_exact_free(&e);
    return result;
}

#ifdef WIGNER_HAVE_MPFR
void wigner9j_mpfr(mpfr_t rop,
                   int tj11, int tj12, int tj13,
                   int tj21, int tj22, int tj23,
                   int tj31, int tj32, int tj33,
                   mpfr_rnd_t rnd)
{
    wigner_exact_t e;
    wigner9j_exact(tj11,tj12,tj13, tj21,tj22,tj23, tj31,tj32,tj33, &e);
    wigner_exact_to_mpfr(rop, &e, rnd);
    wigner_exact_free(&e);
}
#endif
