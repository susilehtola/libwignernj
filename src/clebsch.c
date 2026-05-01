/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Clebsch-Gordan coefficients via the Wigner 3j symbol.
 *
 * Definition:
 *   <j1 m1; j2 m2 | J M> = (-1)^(j1-j2+M) * sqrt(2J+1)
 *                           * (j1 j2  J )
 *                             (m1 m2 -M)
 *
 * Arguments are passed as twice their value (tj = 2*j, tm = 2*m).
 * Requires m1+m2 = M (else zero).
 */
#include "wigner_exact.h"
#include "primes.h"
#include "wigner.h"

/*
 * Internal: compute the exact CG coefficient.
 * The 3j symbol is computed exactly; sqrt(2J+1) is folded in by factoring
 * (tJ+1) and distributing prime powers into int_num and sqrt_num.
 */
static void clebsch_gordan_exact(int tj1, int tm1, int tj2, int tm2,
                                  int tJ, int tM,
                                  wigner_exact_t *out)
{
    int phase;
    wigner_exact_init(out);

    /* m conservation */
    if (tm1 + tm2 != tM) { out->is_zero = 1; return; }

    /* Compute the underlying 3j symbol exactly */
    wigner3j_exact(tj1, tj2, tJ, tm1, tm2, -tM, out);
    if (out->is_zero) return;

    /* Phase: (-1)^(j1-j2+M) = (-1)^((tj1-tj2+tM)/2) */
    phase = (((tj1 - tj2 + tM) / 2) & 1) ? -1 : 1;
    out->sign *= phase;

    /* Multiply by sqrt(2J+1): factor (tJ+1) directly.
     * Even prime powers p^(2k) → p^k into int_num (outside sqrt).
     * Odd prime powers p^(2k+1) → p^k into int_num, p into sqrt_num. */
    {
        int k = tJ + 1, pi;
        for (pi = 0; pi < g_nprimes && g_primes[pi] <= k; pi++) {
            int cnt = 0;
            while (k % g_primes[pi] == 0) { cnt++; k /= g_primes[pi]; }
            if (cnt == 0) continue;
            if (cnt / 2 > 0)
                bigint_mul_prime_pow(&out->int_num,  (uint64_t)g_primes[pi], cnt / 2);
            if (cnt & 1)
                bigint_mul_prime_pow(&out->sqrt_num, (uint64_t)g_primes[pi], 1);
            if (k == 1) break;
        }
    }
}

float clebsch_gordan_f(int tj1, int tm1, int tj2, int tm2, int tJ, int tM)
{
    wigner_exact_t e; float r;
    clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &e);
    r = wigner_exact_to_float(&e);
    wigner_exact_free(&e);
    return r;
}

double clebsch_gordan(int tj1, int tm1, int tj2, int tm2, int tJ, int tM)
{
    wigner_exact_t e; double r;
    clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &e);
    r = wigner_exact_to_double(&e);
    wigner_exact_free(&e);
    return r;
}

long double clebsch_gordan_l(int tj1, int tm1, int tj2, int tm2, int tJ, int tM)
{
    wigner_exact_t e; long double r;
    clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &e);
    r = wigner_exact_to_long_double(&e);
    wigner_exact_free(&e);
    return r;
}

#ifdef WIGNER_HAVE_MPFR
#include <mpfr.h>
void clebsch_gordan_mpfr(mpfr_t rop, int tj1, int tm1, int tj2, int tm2,
                                     int tJ, int tM, mpfr_rnd_t rnd)
{
    wigner_exact_t e;
    clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &e);
    wigner_exact_to_mpfr(rop, &e, rnd);
    wigner_exact_free(&e);
}
#endif
