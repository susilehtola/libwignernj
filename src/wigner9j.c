/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Wigner 9j symbol as a sum over an intermediate quantum number k of
 * products of three Wigner 6j symbols (each evaluated in long double via the
 * exact path).  All three precisions share this implementation.
 *
 * Formula (derived from the Racah W-coefficient representation):
 *   {j11 j12 j13}
 *   {j21 j22 j23} = (-1)^(j13+j22+j31) * sum_k (2k+1)
 *   {j31 j32 j33}
 *       * {j11 j21 j31} * {j11 j12 j13} * {j22 j21 j23}
 *         {j32 j33  k }   {j23 j33  k }   { k  j12 j32}
 *
 * The overall phase (-1)^(j13+j22+j31) [in 2j units: (-1)^(tj13+tj22+tj31)]
 * comes from the Racah W to 6j conversion: W(a,b,c,d;e,f)=(-1)^(a+b+c+d)*6j.
 *
 * k range: max(|j11-j33|,|j21-j32|,|j12-j23|) <= k <= min(j11+j33,j21+j32,j12+j23)
 */
#include "wigner.h"
#include "primes.h"
#include <stdlib.h>

static int triangle_ok(int ta, int tb, int tc)
{
    if ((ta + tb + tc) & 1) return 0;
    if (tc < abs(ta - tb)) return 0;
    if (tc > ta + tb)      return 0;
    return 1;
}

static int max3(int a, int b, int c)
{
    int m = a; if (b > m) m = b; if (c > m) m = c; return m;
}
static int min3(int a, int b, int c)
{
    int m = a; if (b < m) m = b; if (c < m) m = c; return m;
}

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

static long double wigner9j_long_double(
    int tj11, int tj12, int tj13,
    int tj21, int tj22, int tj23,
    int tj31, int tj32, int tj33)
{
    int tk, tk_min, tk_max;
    long double result = 0.0L;

    primes_init();

    if (!selection_rules_9j(tj11,tj12,tj13,tj21,tj22,tj23,tj31,tj32,tj33))
        return 0.0L;

    /* k range: k must satisfy all three triangle conditions in the 6j factors */
    tk_min = max3(abs(tj11 - tj33), abs(tj12 - tj23), abs(tj21 - tj32));
    tk_max = min3(tj11 + tj33,       tj12 + tj23,       tj21 + tj32);

    /* k must have same parity constraints: tk_min parity check */
    /* Ensure tk_min has the right parity (step by 2) */
    if (((tk_min + tj11 + tj33) & 1) != 0) tk_min++;
    if (((tk_max + tj11 + tj33) & 1) != 0) tk_max--;

    if (tk_min > tk_max) return 0.0L;

    /* Overall phase: (-1)^(tj13+tj22+tj31) from Racah W to 6j conversion */
    {
        int outer_phase = (((tj13 + tj22 + tj31) & 1) == 0) ? 1 : -1;

        for (tk = tk_min; tk <= tk_max; tk += 2) {
            long double w1, w2, w3;

            /* {j11 j21 j31; j32 j33 k} */
            w1 = wigner6j_l(tj11, tj21, tj31, tj32, tj33, tk);
            /* {j11 j12 j13; j23 j33 k} */
            w2 = wigner6j_l(tj11, tj12, tj13, tj23, tj33, tk);
            /* {j22 j21 j23; k j12 j32} */
            w3 = wigner6j_l(tj22, tj21, tj23, tk,   tj12, tj32);

            result += (long double)(tk + 1) * w1 * w2 * w3;
        }
        result *= (long double)outer_phase;
    }
    return result;
}

float wigner9j_f(int tj11, int tj12, int tj13,
                 int tj21, int tj22, int tj23,
                 int tj31, int tj32, int tj33)
{
    return (float)wigner9j_long_double(tj11,tj12,tj13,
                                       tj21,tj22,tj23,
                                       tj31,tj32,tj33);
}

double wigner9j(int tj11, int tj12, int tj13,
                int tj21, int tj22, int tj23,
                int tj31, int tj32, int tj33)
{
    return (double)wigner9j_long_double(tj11,tj12,tj13,
                                        tj21,tj22,tj23,
                                        tj31,tj32,tj33);
}

long double wigner9j_l(int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33)
{
    return wigner9j_long_double(tj11,tj12,tj13,
                                tj21,tj22,tj23,
                                tj31,tj32,tj33);
}
