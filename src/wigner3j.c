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
#include "scratch.h"
#include "wigner.h"
#include <stdlib.h>      /* abs */

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
    wigner_scratch_t *scratch;
    bigint_ws_t *ws;
    pfrac_t  *outer, *term;
    int      *lcm_exp;
    bigint_t *sum_pos, *sum_neg, *scaled;
    size_t mw;
    int lcm_max_idx;
    int n_terms;

    /* Reset (rather than init) `out`: callers always pass a previously-
     * initialised wigner_exact_t (either freshly allocated or held in
     * the cached scratch).  Reset preserves the bigint capacities so
     * subsequent reserves are no-ops once they have been sized once. */
    wigner_exact_reset(out);

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

    /* Acquire the thread's cached scratch and bring every buffer up to
     * the size required by the current call (no-op if it is already
     * big enough from a prior, larger call). */
    scratch  = wigner_scratch_acquire();
    ws       = &scratch->ws;
    outer    = &scratch->pfracs[0];
    lcm_exp  =  scratch->lcm_exp[0];
    sum_pos  = &scratch->bigints[0];
    sum_neg  = &scratch->bigints[1];
    scaled   = &scratch->bigints[2];

    bigint_ws_reserve(ws, mw);
    pfrac_zero(outer);
    wigner_scratch_lcm_clear(scratch, 0);
    bigint_set_zero(sum_pos);
    bigint_set_zero(sum_neg);
    bigint_set_zero(scaled);
    bigint_reserve(sum_pos, mw);
    bigint_reserve(sum_neg, mw);
    bigint_reserve(scaled,  mw);

    /* Compute outer sqrt factor */
    build_outer_sqrt_3j(outer, tj1, tj2, tj3, tm1, tm2, tm3);

    /* ── Pass 1: build each term pfrac once into the cache, find LCM ── */
    n_terms = s_max - s_min + 1;
    wigner_scratch_terms_reserve(scratch, n_terms);
    lcm_max_idx = 0;
    for (s = s_min; s <= s_max; s++) {
        /* Build denominator pfrac for term s:
         * s! * (j1+j2-j3-s)! * (j1-m1-s)! * (j2+m2-s)! *
         * (j3-j2+m1+s)! * (j3-j1-m2+s)!
         */
        term = &scratch->terms[s - s_min];
        pfrac_zero(term);
        pfrac_mul_factorial(term, s);
        pfrac_mul_factorial(term, (tj1 + tj2 - tj3) / 2 - s);
        pfrac_mul_factorial(term, (tj1 - tm1) / 2 - s);
        pfrac_mul_factorial(term, (tj2 + tm2) / 2 - s);
        pfrac_mul_factorial(term, (tj3 - tj2 + tm1) / 2 + s);
        pfrac_mul_factorial(term, (tj3 - tj1 - tm2) / 2 + s);

        for (pi = 0; pi < term->max_idx; pi++) {
            if (term->exp[pi] > lcm_exp[pi])
                lcm_exp[pi] = term->exp[pi];
        }
        if (term->max_idx > lcm_max_idx) lcm_max_idx = term->max_idx;
    }
    wigner_scratch_lcm_dirty(scratch, 0, lcm_max_idx);

    /* ── Pass 2: walk cached pfracs, accumulate scaled integer sum ── */
    for (s = s_min; s <= s_max; s++) {
        term = &scratch->terms[s - s_min];

        /* scaled_term = LCM / term_denom = prod p_i^(lcm_exp[i] - term.exp[i]).
         * term.max_idx <= lcm_max_idx by construction. */
        bigint_set_u64(scaled, 1);
        for (pi = 0; pi < lcm_max_idx; pi++) {
            int diff = lcm_exp[pi] - term->exp[pi];
            if (diff > 0)
                bigint_mul_prime_pow_ws(scaled, (uint64_t)g_primes[pi], diff, ws);
        }

        /* Accumulate with sign (-1)^s */
        if ((s & 1) == 0) bigint_add(sum_pos, sum_pos, scaled);
        else              bigint_add(sum_neg, sum_neg, scaled);
    }

    out->sum_sign = bigint_sub_signed(&out->sum, sum_pos, sum_neg);

    /* ── Build output bigints from outer sqrt factor and LCM ── */
    /* Reserve before set_u64: bigint_set_u64 calls bigint_ensure(.., 1)
     * which would otherwise grow cap from 0 to 4 first (one realloc),
     * and the subsequent bigint_reserve(.., mw) would then grow
     * 4→mw (a second realloc).  Reserving first collapses the two
     * reallocs into one. */
    bigint_reserve(&out->int_num,  mw);
    bigint_reserve(&out->int_den,  mw);
    bigint_reserve(&out->sqrt_num, mw);
    bigint_reserve(&out->sqrt_den, mw);
    bigint_set_u64(&out->int_num,  1);
    bigint_set_u64(&out->int_den,  1);
    bigint_set_u64(&out->sqrt_num, 1);
    bigint_set_u64(&out->sqrt_den, 1);

    pfrac_to_sqrt_rational_ws(outer,
                               &out->int_num, &out->int_den,
                               &out->sqrt_num, &out->sqrt_den, ws);

    for (pi = 0; pi < lcm_max_idx; pi++) {
        if (lcm_exp[pi] > 0)
            bigint_mul_prime_pow_ws(&out->int_den,
                                    (uint64_t)g_primes[pi], lcm_exp[pi], ws);
    }

    wigner_scratch_relinquish(scratch);
}

/* ── public API ──────────────────────────────────────────────────────────── */

int wigner3j_max_factorial(int tj1, int tj2, int tj3,
                           int tm1, int tm2, int tm3)
{
    /* The largest factorial referenced in wigner3j_exact is the Δ
     * denominator (j1+j2+j3+1)!.  All m-dependent and Racah-sum
     * factorials have arguments <= (j1+j2+j3)/2. */
    (void)tm1; (void)tm2; (void)tm3;
    return (tj1 + tj2 + tj3) / 2 + 1;
}

/* The public wrappers compute through the cached wigner_exact_t held in
 * the thread's scratch.  wigner3j_exact resets it (preserving bigint
 * capacities) before each call, so on every call after the first the
 * output bigints reuse their previously-allocated buffers and no
 * allocation happens.  wigner_exact_free is therefore not called -- the
 * scratch's destructor takes care of it on thread exit. */

float wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    float result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
    result = wigner_exact_to_float(&s->exact);
    wigner_scratch_relinquish(s);
    return result;
}

double wigner3j(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    double result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
    result = wigner_exact_to_double(&s->exact);
    wigner_scratch_relinquish(s);
    return result;
}

long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    long double result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
    result = wigner_exact_to_long_double(&s->exact);
    wigner_scratch_relinquish(s);
    return result;
}

#ifdef WIGNER_HAVE_QUADMATH
__float128 wigner3j_q(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    __float128 result;
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
    result = wigner_exact_to_float128(&s->exact);
    wigner_scratch_relinquish(s);
    return result;
}
#endif

#ifdef WIGNER_HAVE_MPFR
void wigner3j_mpfr(mpfr_t rop, int tj1, int tj2, int tj3,
                               int tm1, int tm2, int tm3, mpfr_rnd_t rnd)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
    wigner_exact_to_mpfr(rop, &s->exact, rnd);
    wigner_scratch_relinquish(s);
}
#endif
