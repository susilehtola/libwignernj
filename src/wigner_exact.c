/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "wigner_exact.h"
#include <math.h>

void wigner_exact_init(wigner_exact_t *e)
{
    bigint_init(&e->sum);
    bigint_init(&e->int_num);
    bigint_init(&e->int_den);
    bigint_init(&e->sqrt_num);
    bigint_init(&e->sqrt_den);
    e->sum_sign = 1;
    e->sign     = 1;
    e->is_zero  = 0;
}

void wigner_exact_free(wigner_exact_t *e)
{
    if (e->is_zero) return;
    bigint_free(&e->sum);
    bigint_free(&e->int_num);
    bigint_free(&e->int_den);
    bigint_free(&e->sqrt_num);
    bigint_free(&e->sqrt_den);
}

/*
 * Floating-point conversion helpers.
 *
 * The exact value is:  sign * sum_sign * sum * int_num / int_den
 *                      * sqrt(sqrt_num / sqrt_den)
 *
 * For large angular momenta the individual bigints can exceed the range of
 * double (or float), so we use bigint_frexp* to extract a normalised mantissa
 * in [0.5, 1) and a binary exponent, combine the exponents in integer
 * arithmetic, and call ldexp only at the final step.
 */

/* Apply the sqrt exponent adjustment.  Returns the adjusted mantissa and
 * adds the integer part of sqrt_diff/2 to *iexp. */
static double apply_sqrt_exp_d(double m, int sqrt_diff, int *iexp)
{
    if (sqrt_diff & 1) {
        if (sqrt_diff > 0) {
            m      *= sqrt(2.0);
            *iexp  += (sqrt_diff - 1) / 2;
        } else {
            m      /= sqrt(2.0);
            *iexp  += (sqrt_diff + 1) / 2;
        }
    } else {
        *iexp += sqrt_diff / 2;
    }
    return m;
}

static float apply_sqrt_exp_f(float m, int sqrt_diff, int *iexp)
{
    if (sqrt_diff & 1) {
        if (sqrt_diff > 0) {
            m      *= sqrtf(2.0f);
            *iexp  += (sqrt_diff - 1) / 2;
        } else {
            m      /= sqrtf(2.0f);
            *iexp  += (sqrt_diff + 1) / 2;
        }
    } else {
        *iexp += sqrt_diff / 2;
    }
    return m;
}

static long double apply_sqrt_exp_l(long double m, int sqrt_diff, int *iexp)
{
    if (sqrt_diff & 1) {
        if (sqrt_diff > 0) {
            m      *= sqrtl(2.0L);
            *iexp  += (sqrt_diff - 1) / 2;
        } else {
            m      /= sqrtl(2.0L);
            *iexp  += (sqrt_diff + 1) / 2;
        }
    } else {
        *iexp += sqrt_diff / 2;
    }
    return m;
}

float wigner_exact_to_float(const wigner_exact_t *e)
{
    int es, en, ed, esn, esd, iexp;
    float ms, mn, md, msn, msd, m;

    if (e->is_zero) return 0.0f;
    if (bigint_is_zero(&e->sum)) return 0.0f;

    ms  = bigint_frexpf(&e->sum,      &es);
    mn  = bigint_frexpf(&e->int_num,  &en);
    md  = bigint_frexpf(&e->int_den,  &ed);
    msn = bigint_frexpf(&e->sqrt_num, &esn);
    msd = bigint_frexpf(&e->sqrt_den, &esd);

    m    = ms * (mn / md) * sqrtf(msn / msd);
    iexp = es + en - ed;
    m    = apply_sqrt_exp_f(m, esn - esd, &iexp);

    return (float)(e->sign * e->sum_sign) * ldexpf(m, iexp);
}

double wigner_exact_to_double(const wigner_exact_t *e)
{
    int es, en, ed, esn, esd, iexp;
    double ms, mn, md, msn, msd, m;

    if (e->is_zero) return 0.0;
    if (bigint_is_zero(&e->sum)) return 0.0;

    ms  = bigint_frexp(&e->sum,      &es);
    mn  = bigint_frexp(&e->int_num,  &en);
    md  = bigint_frexp(&e->int_den,  &ed);
    msn = bigint_frexp(&e->sqrt_num, &esn);
    msd = bigint_frexp(&e->sqrt_den, &esd);

    m    = ms * (mn / md) * sqrt(msn / msd);
    iexp = es + en - ed;
    m    = apply_sqrt_exp_d(m, esn - esd, &iexp);

    return (double)(e->sign * e->sum_sign) * ldexp(m, iexp);
}

#ifdef WIGNER_HAVE_QUADMATH
static __float128 apply_sqrt_exp_q(__float128 m, int sqrt_diff, int *iexp)
{
    if (sqrt_diff & 1) {
        if (sqrt_diff > 0) {
            m      *= sqrtq(2.0Q);
            *iexp  += (sqrt_diff - 1) / 2;
        } else {
            m      /= sqrtq(2.0Q);
            *iexp  += (sqrt_diff + 1) / 2;
        }
    } else {
        *iexp += sqrt_diff / 2;
    }
    return m;
}

__float128 wigner_exact_to_float128(const wigner_exact_t *e)
{
    int es, en, ed, esn, esd, iexp;
    __float128 ms, mn, md, msn, msd, m;

    if (e->is_zero) return 0.0Q;
    if (bigint_is_zero(&e->sum)) return 0.0Q;

    ms  = bigint_frexp_q(&e->sum,      &es);
    mn  = bigint_frexp_q(&e->int_num,  &en);
    md  = bigint_frexp_q(&e->int_den,  &ed);
    msn = bigint_frexp_q(&e->sqrt_num, &esn);
    msd = bigint_frexp_q(&e->sqrt_den, &esd);

    m    = ms * (mn / md) * sqrtq(msn / msd);
    iexp = es + en - ed;
    m    = apply_sqrt_exp_q(m, esn - esd, &iexp);

    return (__float128)(e->sign * e->sum_sign) * ldexpq(m, iexp);
}
#endif

long double wigner_exact_to_long_double(const wigner_exact_t *e)
{
    int es, en, ed, esn, esd, iexp;
    long double ms, mn, md, msn, msd, m;

    if (e->is_zero) return 0.0L;
    if (bigint_is_zero(&e->sum)) return 0.0L;

    ms  = bigint_frexpl(&e->sum,      &es);
    mn  = bigint_frexpl(&e->int_num,  &en);
    md  = bigint_frexpl(&e->int_den,  &ed);
    msn = bigint_frexpl(&e->sqrt_num, &esn);
    msd = bigint_frexpl(&e->sqrt_den, &esd);

    m    = ms * (mn / md) * sqrtl(msn / msd);
    iexp = es + en - ed;
    m    = apply_sqrt_exp_l(m, esn - esd, &iexp);

    return (long double)(e->sign * e->sum_sign) * ldexpl(m, iexp);
}

/* ── MPFR conversion ─────────────────────────────────────────────────────── */

#ifdef WIGNER_HAVE_MPFR
void wigner_exact_to_mpfr(mpfr_t rop, const wigner_exact_t *e, mpfr_rnd_t rnd)
{
    mpfr_prec_t prec = mpfr_get_prec(rop);
    mpfr_t tmp, aux;

    if (e->is_zero || bigint_is_zero(&e->sum)) {
        mpfr_set_zero(rop, +1);
        return;
    }

    mpfr_init2(tmp, prec);
    mpfr_init2(aux, prec);

    /* rop = sum * int_num / int_den */
    bigint_to_mpfr(rop, &e->sum,     rnd);
    bigint_to_mpfr(tmp, &e->int_num, rnd);
    mpfr_mul(rop, rop, tmp, rnd);
    bigint_to_mpfr(tmp, &e->int_den, rnd);
    mpfr_div(rop, rop, tmp, rnd);

    /* rop *= sqrt(sqrt_num / sqrt_den) */
    bigint_to_mpfr(tmp, &e->sqrt_num, rnd);
    bigint_to_mpfr(aux, &e->sqrt_den, rnd);
    mpfr_div(tmp, tmp, aux, rnd);
    mpfr_sqrt(tmp, tmp, rnd);
    mpfr_mul(rop, rop, tmp, rnd);

    if (e->sign * e->sum_sign < 0)
        mpfr_neg(rop, rop, MPFR_RNDN);

    mpfr_clear(tmp);
    mpfr_clear(aux);
}
#endif
