/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "pfrac.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void pfrac_init(pfrac_t *f)
{
    f->exp = (int *)calloc((size_t)g_nprimes, sizeof(int));
}

void pfrac_free(pfrac_t *f)
{
    free(f->exp);
    f->exp = NULL;
}

void pfrac_zero(pfrac_t *f)
{
    memset(f->exp, 0, (size_t)g_nprimes * sizeof(int));
}

void pfrac_copy(pfrac_t *dst, const pfrac_t *src)
{
    memcpy(dst->exp, src->exp, (size_t)g_nprimes * sizeof(int));
}

void pfrac_mul_factorial(pfrac_t *f, int n)
{
    int i;
    if (n <= 1) return;
    for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
        f->exp[i] += legendre_valuation(n, i);
}

void pfrac_div_factorial(pfrac_t *f, int n)
{
    int i;
    if (n <= 1) return;
    for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
        f->exp[i] -= legendre_valuation(n, i);
}

void pfrac_mul_int(pfrac_t *f, int k)
{
    int pk, i;
    if (k <= 1) return;
    /* Factor k into primes and add to exp[] */
    for (i = 0; i < g_nprimes && g_primes[i] <= k; i++) {
        pk = g_primes[i];
        while (k % pk == 0) {
            f->exp[i]++;
            k /= pk;
        }
        if (k == 1) break;
    }
    assert(k == 1); /* k must have no prime factor > PRIME_SIEVE_LIMIT */
}

void pfrac_to_sqrt_rational(const pfrac_t *f,
                             bigint_t *int_num,  bigint_t *int_den,
                             bigint_t *sqrt_num, bigint_t *sqrt_den)
{
    int i, e, ae, half;
    for (i = 0; i < g_nprimes; i++) {
        e = f->exp[i];
        if (e == 0) continue;
        ae   = (e > 0) ? e : -e;
        half = ae / 2;
        if (e > 0) {
            if (half > 0) bigint_mul_prime_pow(int_num,  (uint64_t)g_primes[i], half);
            if (ae & 1)   bigint_mul_prime_pow(sqrt_num, (uint64_t)g_primes[i], 1);
        } else {
            if (half > 0) bigint_mul_prime_pow(int_den,  (uint64_t)g_primes[i], half);
            if (ae & 1)   bigint_mul_prime_pow(sqrt_den, (uint64_t)g_primes[i], 1);
        }
    }
}
