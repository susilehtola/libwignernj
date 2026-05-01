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
