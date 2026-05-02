/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Wigner 3j symbol via the Racah formula with prime-factorization exact
 * arithmetic following Johansson & Forssén (doi:10.1137/15M1021908).
 *
 * All arguments are passed as 2*j integers (e.g. j=3/2 → tj=3).
 *
 * Racah formula:
 *   (j1 j2 j3)  = (-1)^(j1-j2-m3) * Delta(j1,j2,j3)
 *   (m1 m2 m3)    * sqrt[(j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j3+m3)!(j3-m3)!]
 *                 * sum_s (-1)^s / [s!(j1+j2-j3-s)!(j1-m1-s)!
 *                                    (j2+m2-s)!(j3-j2+m1+s)!(j3-j1-m2+s)!]
 *
 * Delta(a,b,c) = sqrt[(a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!]
 */
#include "wigner_exact.h"
#include "pfrac.h"
#include "primes.h"
#include "wigner.h"
#include "xalloc.h"
#include <stdlib.h>

/* ── selection rules ─────────────────────────────────────────────────────── */

static int triangle_ok(int ta, int tb, int tc)
{
    /* ta+tb+tc must be even (so that j-sums are integers) */
    if ((ta + tb + tc) & 1) return 0;
    if (tc < abs(ta - tb)) return 0;
    if (tc > ta + tb)      return 0;
    return 1;
}

static int selection_rules_3j(int tj1, int tj2, int tj3,
                               int tm1, int tm2, int tm3)
{
    /* m conservation */
    if (tm1 + tm2 + tm3 != 0) return 0;
    /* m within j */
    if (abs(tm1) > tj1 || abs(tm2) > tj2 || abs(tm3) > tj3) return 0;
    /* same half-integer type: tmi must have same parity as tji */
    if (((tj1 - tm1) & 1) || ((tj2 - tm2) & 1) || ((tj3 - tm3) & 1)) return 0;
    /* triangle */
    if (!triangle_ok(tj1, tj2, tj3)) return 0;
    return 1;
}

/* ── Racah summation bounds ──────────────────────────────────────────────── */

static void bounds_3j(int tj1, int tj2, int tj3,
                      int tm1, int tm2,
                      int *s_min, int *s_max)
{
    /* All factorial arguments in the Racah sum must be non-negative.
     * Denominator factorials: s!, (j1+j2-j3-s)!, (j1-m1-s)!, (j2+m2-s)!,
     *                         (j3-j2+m1+s)!, (j3-j1-m2+s)!
     * Lower bounds (from last two): s >= j2-j3-m1 and s >= j1+m2-j3
     * Upper bounds (from first four): s <= j1+j2-j3, s <= j1-m1, s <= j2+m2
     * All in 2x representation, divided by 2.
     */
    int a = (tj2 - tj3 - tm1) / 2;   /* lower from (j3-j2+m1+s)! */
    int b = (tj1 + tm2 - tj3) / 2;   /* lower from (j3-j1-m2+s)! */
    *s_min = (a > b) ? a : b;
    if (*s_min < 0) *s_min = 0;

    a = (tj1 + tj2 - tj3) / 2;
    b = (tj1 - tm1) / 2;
    int c = (tj2 + tm2) / 2;
    *s_max = a;
    if (b < *s_max) *s_max = b;
    if (c < *s_max) *s_max = c;
}

/* ── core exact computation ──────────────────────────────────────────────── */

/*
 * Build the outer sqrt factor pfrac_t for the 3j symbol.
 * This represents the argument under the outer sqrt sign:
 *   (j1+j2-j3)!(j1-j2+j3)!(-j1+j2+j3)! * (j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j3+m3)!(j3-m3)!
 *   / (j1+j2+j3+1)!
 */
static void build_outer_sqrt_3j(pfrac_t *outer,
                                 int tj1, int tj2, int tj3,
                                 int tm1, int tm2, int tm3)
{
    pfrac_zero(outer);
    /* Delta numerator factorials */
    pfrac_mul_factorial(outer, (tj1 + tj2 - tj3) / 2);
    pfrac_mul_factorial(outer, (tj1 - tj2 + tj3) / 2);
    pfrac_mul_factorial(outer, (-tj1 + tj2 + tj3) / 2);
    /* m-dependent sqrt numerator factorials */
    pfrac_mul_factorial(outer, (tj1 + tm1) / 2);
    pfrac_mul_factorial(outer, (tj1 - tm1) / 2);
    pfrac_mul_factorial(outer, (tj2 + tm2) / 2);
    pfrac_mul_factorial(outer, (tj2 - tm2) / 2);
    pfrac_mul_factorial(outer, (tj3 + tm3) / 2);
    pfrac_mul_factorial(outer, (tj3 - tm3) / 2);
    /* Delta denominator factorial */
    pfrac_div_factorial(outer, (tj1 + tj2 + tj3) / 2 + 1);
}

void wigner3j_exact(int tj1, int tj2, int tj3,
                    int tm1, int tm2, int tm3,
                    wigner_exact_t *out)
{
    int s, s_min, s_max, pi;
    int *lcm_exp;        /* max denominator exponent per prime across all s */
    pfrac_t outer, term;
    bigint_t sum_pos, sum_neg, scaled;
    bigint_ws_t ws;
    size_t mw;

    wigner_exact_init(out);

    if (!selection_rules_3j(tj1, tj2, tj3, tm1, tm2, tm3)) {
        out->is_zero = 1;
        return;
    }

    bounds_3j(tj1, tj2, tj3, tm1, tm2, &s_min, &s_max);
    if (s_min > s_max) { out->is_zero = 1; return; }

    /* Overall phase: (-1)^(j1-j2-m3) */
    out->sign = ((((tj1 - tj2 - tm3) / 2) & 1) == 0) ? 1 : -1;

    /* Largest factorial argument is bounded by (tj1+tj2+tj3)/2+1 plus
     * (tji+|tmi|)/2 ≤ tji.  Use a safe upper bound for sizing. */
    mw = bigint_words_for_factorial((tj1 + tj2 + tj3) / 2 + 5);

    bigint_ws_init(&ws);
    bigint_ws_reserve(&ws, mw);

    /* Compute outer sqrt factor */
    pfrac_init(&outer);
    build_outer_sqrt_3j(&outer, tj1, tj2, tj3, tm1, tm2, tm3);

    /* ── Pass 1: find LCM exponents (max denom exponent per prime) ── */
    int lcm_max_idx = 0;
    lcm_exp = (int *)xcalloc((size_t)g_nprimes, sizeof(int));

    pfrac_init(&term);
    for (s = s_min; s <= s_max; s++) {
        /* Build denominator pfrac for term s:
         * s! * (j1+j2-j3-s)! * (j1-m1-s)! * (j2+m2-s)! *
         * (j3-j2+m1+s)! * (j3-j1-m2+s)!
         */
        pfrac_zero(&term);
        pfrac_mul_factorial(&term, s);
        pfrac_mul_factorial(&term, (tj1 + tj2 - tj3) / 2 - s);
        pfrac_mul_factorial(&term, (tj1 - tm1) / 2 - s);
        pfrac_mul_factorial(&term, (tj2 + tm2) / 2 - s);
        pfrac_mul_factorial(&term, (tj3 - tj2 + tm1) / 2 + s);
        pfrac_mul_factorial(&term, (tj3 - tj1 - tm2) / 2 + s);

        for (pi = 0; pi < term.max_idx; pi++) {
            if (term.exp[pi] > lcm_exp[pi])
                lcm_exp[pi] = term.exp[pi];
        }
        if (term.max_idx > lcm_max_idx) lcm_max_idx = term.max_idx;
    }

    /* ── Pass 2: accumulate scaled integer Racah sum (pre-sized bigints) ── */
    bigint_init(&sum_pos); bigint_reserve(&sum_pos, mw);
    bigint_init(&sum_neg); bigint_reserve(&sum_neg, mw);
    bigint_init(&scaled);  bigint_reserve(&scaled,  mw);

    for (s = s_min; s <= s_max; s++) {
        /* Rebuild term s denominator */
        pfrac_zero(&term);
        pfrac_mul_factorial(&term, s);
        pfrac_mul_factorial(&term, (tj1 + tj2 - tj3) / 2 - s);
        pfrac_mul_factorial(&term, (tj1 - tm1) / 2 - s);
        pfrac_mul_factorial(&term, (tj2 + tm2) / 2 - s);
        pfrac_mul_factorial(&term, (tj3 - tj2 + tm1) / 2 + s);
        pfrac_mul_factorial(&term, (tj3 - tj1 - tm2) / 2 + s);

        /* scaled_term = LCM / term_denom = prod p_i^(lcm_exp[i] - term.exp[i]).
         * lcm_exp[i] is zero for i >= lcm_max_idx; term.exp[i] is zero for
         * i >= term.max_idx; so iterating to lcm_max_idx covers every nonzero
         * difference (term.max_idx <= lcm_max_idx by construction). */
        bigint_set_u64(&scaled, 1);
        for (pi = 0; pi < lcm_max_idx; pi++) {
            int diff = lcm_exp[pi] - term.exp[pi];
            if (diff > 0)
                bigint_mul_prime_pow_ws(&scaled, (uint64_t)g_primes[pi], diff, &ws);
        }

        /* Accumulate with sign (-1)^s */
        if ((s & 1) == 0)
            bigint_add(&sum_pos, &sum_pos, &scaled);
        else
            bigint_add(&sum_neg, &sum_neg, &scaled);
    }

    /* sum = sum_pos - sum_neg */
    {
        int sign_diff = bigint_sub_signed(&out->sum, &sum_pos, &sum_neg);
        out->sum_sign = sign_diff;
    }

    /* ── Build output bigints from outer sqrt factor and LCM ── */
    bigint_set_u64(&out->int_num,  1);
    bigint_set_u64(&out->int_den,  1);
    bigint_set_u64(&out->sqrt_num, 1);
    bigint_set_u64(&out->sqrt_den, 1);
    bigint_reserve(&out->int_num,  mw);
    bigint_reserve(&out->int_den,  mw);
    bigint_reserve(&out->sqrt_num, mw);
    bigint_reserve(&out->sqrt_den, mw);

    /* Outer sqrt factor → int_num, int_den, sqrt_num, sqrt_den */
    pfrac_to_sqrt_rational_ws(&outer,
                               &out->int_num, &out->int_den,
                               &out->sqrt_num, &out->sqrt_den, &ws);

    /* LCM denominator: prod p_i^lcm_exp[i] → absorbed into int_den */
    for (pi = 0; pi < lcm_max_idx; pi++) {
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

float wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_exact_t e;
    float result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &e);
    result = wigner_exact_to_float(&e);
    wigner_exact_free(&e);
    return result;
}

double wigner3j(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_exact_t e;
    double result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &e);
    result = wigner_exact_to_double(&e);
    wigner_exact_free(&e);
    return result;
}

long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_exact_t e;
    long double result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &e);
    result = wigner_exact_to_long_double(&e);
    wigner_exact_free(&e);
    return result;
}

#ifdef WIGNER_HAVE_QUADMATH
__float128 wigner3j_q(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_exact_t e;
    __float128 result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &e);
    result = wigner_exact_to_float128(&e);
    wigner_exact_free(&e);
    return result;
}
#endif

#ifdef WIGNER_HAVE_MPFR
void wigner3j_mpfr(mpfr_t rop, int tj1, int tj2, int tj3,
                               int tm1, int tm2, int tm3, mpfr_rnd_t rnd)
{
    wigner_exact_t e;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &e);
    wigner_exact_to_mpfr(rop, &e, rnd);
    wigner_exact_free(&e);
}
#endif
