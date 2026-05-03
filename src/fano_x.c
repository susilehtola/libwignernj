/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Fano X-coefficient via the Wigner 9j symbol.
 *
 * Definition (Fano 1953; Edmonds 1957 §6.4; Tamura 1970):
 *   X(j1 j2 j12; j3 j4 j34; j13 j24 J)
 *     = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)]
 *       * { j1   j2   j12 }
 *         { j3   j4   j34 }
 *         { j13  j24  J   }
 *
 * The X-coefficient is a normalisation variant of the 9j symbol used in
 * the analysis of polarisation correlations and in the Russell--Saunders
 * coupling literature; it carries no information beyond the underlying
 * 9j symbol but is produced directly because some downstream
 * applications expect that normalisation.
 *
 * Arguments are passed as twice their value (tj = 2*j); the four
 * (2j+1) factors above appear in the exact pipeline as (tj12+1),
 * (tj34+1), (tj13+1), and (tj24+1) and are folded into the existing
 * wigner_exact_t in the same way the (2J+1) factor is folded into the
 * Clebsch--Gordan coefficient (see src/clebsch.c).
 */
#include "wigner_exact.h"
#include "primes.h"
#include "scratch.h"
#include "wigner.h"

/*
 * Internal: compute the exact Fano X coefficient.  The 9j symbol is
 * computed exactly by the existing wigner9j_exact pipeline; the four
 * (2j+1) factors above are then folded in by factoring each of them
 * into prime powers and distributing those prime powers over the
 * int_num and sqrt_num bigints of the wigner_exact_t:
 *
 *   even prime power p^(2k) -> p^k into int_num
 *   odd  prime power p^(2k+1) -> p^k into int_num, p into sqrt_num
 *
 * This is the same scheme that src/clebsch.c uses for its single
 * sqrt(2J+1) factor; here it is iterated over the four factors of
 * the Fano X normalisation.
 */
static void fano_x_exact(int tj1, int tj2, int tj12,
                         int tj3, int tj4, int tj34,
                         int tj13, int tj24, int tJ,
                         wigner_exact_t *out)
{
    int factors[4];
    int f;

    /* The underlying 9j */
    wigner9j_exact(tj1, tj2, tj12,
                   tj3, tj4, tj34,
                   tj13, tj24, tJ, out);
    if (out->is_zero) return;

    /* Fold sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] into the exact tuple. */
    factors[0] = tj12 + 1;
    factors[1] = tj34 + 1;
    factors[2] = tj13 + 1;
    factors[3] = tj24 + 1;
    for (f = 0; f < 4; ++f) {
        int k = factors[f];
        int pi;
        for (pi = 0; pi < g_nprimes && g_primes[pi] <= k; ++pi) {
            int cnt = 0;
            while (k % g_primes[pi] == 0) { cnt++; k /= g_primes[pi]; }
            if (cnt == 0) continue;
            if (cnt / 2 > 0)
                bigint_mul_prime_pow(&out->int_num,
                                     (uint64_t)g_primes[pi], cnt / 2);
            if (cnt & 1)
                bigint_mul_prime_pow(&out->sqrt_num,
                                     (uint64_t)g_primes[pi], 1);
            if (k == 1) break;
        }
    }
}

/* ── public API ──────────────────────────────────────────────────────────── */

int fano_x_max_factorial(int tj1, int tj2, int tj12,
                         int tj3, int tj4, int tj34,
                         int tj13, int tj24, int tJ)
{
    /* X = sqrt[...] * 9j; the sqrt prefactor is folded in via
     * bigint_mul_prime_pow without touching the factorial cache. */
    return wigner9j_max_factorial(tj1,  tj2,  tj12,
                                  tj3,  tj4,  tj34,
                                  tj13, tj24, tJ);
}

float fano_x_f(int tj1, int tj2, int tj12,
               int tj3, int tj4, int tj34,
               int tj13, int tj24, int tJ)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    float r;
    fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
    r = wigner_exact_to_float(&s->exact);
    wigner_scratch_relinquish(s);
    return r;
}

double fano_x(int tj1, int tj2, int tj12,
              int tj3, int tj4, int tj34,
              int tj13, int tj24, int tJ)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    double r;
    fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
    r = wigner_exact_to_double(&s->exact);
    wigner_scratch_relinquish(s);
    return r;
}

long double fano_x_l(int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    long double r;
    fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
    r = wigner_exact_to_long_double(&s->exact);
    wigner_scratch_relinquish(s);
    return r;
}

#ifdef WIGNER_HAVE_QUADMATH
__float128 fano_x_q(int tj1, int tj2, int tj12,
                    int tj3, int tj4, int tj34,
                    int tj13, int tj24, int tJ)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    __float128 r;
    fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
    r = wigner_exact_to_float128(&s->exact);
    wigner_scratch_relinquish(s);
    return r;
}
#endif

#ifdef WIGNER_HAVE_MPFR
#include "wigner_mpfr.h"
void fano_x_mpfr(mpfr_t rop,
                 int tj1, int tj2, int tj12,
                 int tj3, int tj4, int tj34,
                 int tj13, int tj24, int tJ,
                 mpfr_rnd_t rnd)
{
    wigner_scratch_t *s = wigner_scratch_acquire();
    fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
    wigner_exact_to_mpfr(rop, &s->exact, rnd);
    wigner_scratch_relinquish(s);
}
#endif
