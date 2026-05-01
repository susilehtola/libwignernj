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
 * The two 3j symbols are each computed via the exact integer arithmetic path.
 * The irrational factor sqrt[(2l1+1)(2l2+1)(2l3+1)/(4π)] is applied as a
 * floating-point multiplication after the exact intermediate step.
 */
#include "wigner.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

float gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3)
{
    float w0, norm;
    if (tm1 + tm2 + tm3 != 0) return 0.0f;
    w0 = wigner3j_f(tl1, tl2, tl3, 0, 0, 0);
    if (w0 == 0.0f) return 0.0f;
    norm = sqrtf((float)(tl1 + 1) * (float)(tl2 + 1) * (float)(tl3 + 1)
                 / (4.0f * (float)M_PI));
    return norm * w0 * wigner3j_f(tl1, tl2, tl3, tm1, tm2, tm3);
}

double gaunt(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3)
{
    double w0, norm;
    if (tm1 + tm2 + tm3 != 0) return 0.0;
    w0 = wigner3j(tl1, tl2, tl3, 0, 0, 0);
    if (w0 == 0.0) return 0.0;
    norm = sqrt((double)(tl1 + 1) * (double)(tl2 + 1) * (double)(tl3 + 1)
                / (4.0 * M_PI));
    return norm * w0 * wigner3j(tl1, tl2, tl3, tm1, tm2, tm3);
}

long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3)
{
    long double w0, norm;
    if (tm1 + tm2 + tm3 != 0) return 0.0L;
    w0 = wigner3j_l(tl1, tl2, tl3, 0, 0, 0);
    if (w0 == 0.0L) return 0.0L;
    norm = sqrtl((long double)(tl1 + 1) * (long double)(tl2 + 1)
                 * (long double)(tl3 + 1) / (4.0L * acosl(-1.0L)));
    return norm * w0 * wigner3j_l(tl1, tl2, tl3, tm1, tm2, tm3);
}

#ifdef WIGNER_HAVE_MPFR
#include "wigner_mpfr.h"
void gaunt_mpfr(mpfr_t rop, int tl1, int tm1, int tl2, int tm2,
                             int tl3, int tm3, mpfr_rnd_t rnd)
{
    mpfr_prec_t prec;
    mpfr_t w0, norm, pi;

    if (tm1 + tm2 + tm3 != 0) { mpfr_set_zero(rop, +1); return; }

    prec = mpfr_get_prec(rop);
    wigner3j_mpfr(rop, tl1, tl2, tl3, 0, 0, 0, rnd);
    if (mpfr_zero_p(rop)) return;

    mpfr_init2(w0,   prec);
    mpfr_init2(norm, prec);
    mpfr_init2(pi,   prec);

    wigner3j_mpfr(w0, tl1, tl2, tl3, tm1, tm2, tm3, rnd);

    /* norm = sqrt[(tl1+1)(tl2+1)(tl3+1) / (4*pi)] */
    mpfr_const_pi(pi, rnd);
    mpfr_mul_ui(pi, pi, 4, rnd);
    mpfr_set_ui(norm, (unsigned long)(tl1 + 1), rnd);
    mpfr_mul_ui(norm, norm, (unsigned long)(tl2 + 1), rnd);
    mpfr_mul_ui(norm, norm, (unsigned long)(tl3 + 1), rnd);
    mpfr_div(norm, norm, pi, rnd);
    mpfr_sqrt(norm, norm, rnd);

    mpfr_mul(rop, rop, w0,   rnd);
    mpfr_mul(rop, rop, norm, rnd);

    mpfr_clear(w0);
    mpfr_clear(norm);
    mpfr_clear(pi);
}
#endif
