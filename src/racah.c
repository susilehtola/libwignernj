/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Racah W-coefficient via the Wigner 6j symbol.
 *
 * Definition:
 *   W(j1 j2 J j3; j12 j23) = (-1)^(j1+j2+J+j3) * {j1  j2  j12}
 *                                                    {j3  J   j23}
 *
 * Arguments are passed as twice their value (tj = 2*j).
 */
#include "wigner.h"

float racah_w_f(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23)
{
    int phase = ((((tj1 + tj2 + tJ + tj3) / 2) & 1) == 0) ? 1 : -1;
    return (float)phase * wigner6j_f(tj1, tj2, tj12, tj3, tJ, tj23);
}

double racah_w(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23)
{
    int phase = ((((tj1 + tj2 + tJ + tj3) / 2) & 1) == 0) ? 1 : -1;
    return (double)phase * wigner6j(tj1, tj2, tj12, tj3, tJ, tj23);
}

long double racah_w_l(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23)
{
    int phase = ((((tj1 + tj2 + tJ + tj3) / 2) & 1) == 0) ? 1 : -1;
    return (long double)phase * wigner6j_l(tj1, tj2, tj12, tj3, tJ, tj23);
}
