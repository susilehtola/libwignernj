/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef WIGNER_EXACT_H
#define WIGNER_EXACT_H

#include "bigint.h"

/*
 * Exact symbolic result of a Wigner symbol computation.
 *
 * The floating-point value is:
 *   sign * sum_sign * bigint_to_T(sum) * bigint_to_T(int_num)
 *                   / bigint_to_T(int_den)
 *                   * sqrt(bigint_to_T(sqrt_num) / bigint_to_T(sqrt_den))
 *
 * where T is float, double, or long double.
 *
 * int_den absorbs both the denominator from the outer sqrt factor and
 * the LCM of all Racah sum term denominators.
 *
 * When is_zero is set the symbol vanishes by selection rules; all bigints
 * hold their init-time zero state and carry no useful data.
 */
typedef struct {
    bigint_t sum;        /* magnitude of Racah integer sum                */
    int      sum_sign;   /* sign of Racah sum: +1 or -1 (1 when sum == 0) */
    bigint_t int_num;    /* integer-part numerator of outer sqrt factor   */
    bigint_t int_den;    /* integer-part denominator (includes LCM)       */
    bigint_t sqrt_num;   /* squarefree sqrt numerator                     */
    bigint_t sqrt_den;   /* squarefree sqrt denominator                   */
    int      sign;       /* overall phase: +1 or -1                       */
    int      is_zero;    /* 1 if symbol is identically zero               */
} wigner_exact_t;

void wigner_exact_init(wigner_exact_t *e);
void wigner3j_exact(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3,
                    wigner_exact_t *out);
void wigner6j_exact(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6,
                    wigner_exact_t *out);
void wigner9j_exact(int tj11, int tj12, int tj13,
                    int tj21, int tj22, int tj23,
                    int tj31, int tj32, int tj33,
                    wigner_exact_t *out);
void wigner_exact_free(wigner_exact_t *e);

float       wigner_exact_to_float      (const wigner_exact_t *e);
double      wigner_exact_to_double     (const wigner_exact_t *e);
long double wigner_exact_to_long_double(const wigner_exact_t *e);

#ifdef WIGNER_HAVE_MPFR
#include <mpfr.h>
void wigner_exact_to_mpfr(mpfr_t rop, const wigner_exact_t *e, mpfr_rnd_t rnd);
#endif

#endif /* WIGNER_EXACT_H */
