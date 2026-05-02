/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "pfrac.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void pfrac_fatal(const char *msg)
{
    fprintf(stderr, "libwignernj: %s\n", msg);
    abort();
}

void pfrac_init(pfrac_t *f)
{
    f->exp = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
    f->max_idx = 0;
}

void pfrac_free(pfrac_t *f)
{
    free(f->exp);
    f->exp = NULL;
    f->max_idx = 0;
}

void pfrac_zero(pfrac_t *f)
{
    /* Zero only [0, max_idx); the invariant guarantees the rest is zero. */
    if (f->max_idx > 0)
        memset(f->exp, 0, (size_t)f->max_idx * sizeof(int));
    f->max_idx = 0;
}

void pfrac_copy(pfrac_t *dst, const pfrac_t *src)
{
    /* If dst was wider than src, zero the excess so the invariant holds. */
    if (dst->max_idx > src->max_idx) {
        memset(dst->exp + src->max_idx, 0,
               (size_t)(dst->max_idx - src->max_idx) * sizeof(int));
    }
    if (src->max_idx > 0)
        memcpy(dst->exp, src->exp, (size_t)src->max_idx * sizeof(int));
    dst->max_idx = src->max_idx;
}

void pfrac_mul_factorial(pfrac_t *f, int n)
{
    int i;
    if (n <= 1) return;
    if (n > MAX_FACTORIAL_ARG) {
        fprintf(stderr,
            "libwignernj: factorial argument %d exceeds MAX_FACTORIAL_ARG=%d "
            "(angular momenta too large for prime table)\n",
            n, MAX_FACTORIAL_ARG);
        pfrac_fatal("aborting");
    }
    for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
        f->exp[i] += legendre_valuation(n, i);
    if (i > f->max_idx) f->max_idx = i;
}

void pfrac_div_factorial(pfrac_t *f, int n)
{
    int i;
    if (n <= 1) return;
    if (n > MAX_FACTORIAL_ARG) {
        fprintf(stderr,
            "libwignernj: factorial argument %d exceeds MAX_FACTORIAL_ARG=%d "
            "(angular momenta too large for prime table)\n",
            n, MAX_FACTORIAL_ARG);
        pfrac_fatal("aborting");
    }
    for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
        f->exp[i] -= legendre_valuation(n, i);
    if (i > f->max_idx) f->max_idx = i;
}

void pfrac_mul_int(pfrac_t *f, int k)
{
    int pk, i;
    int last = -1;
    if (k <= 1) return;
    if (k > PRIME_SIEVE_LIMIT) {
        fprintf(stderr,
            "libwignernj: pfrac_mul_int argument %d exceeds PRIME_SIEVE_LIMIT=%d "
            "(angular momenta too large for prime table)\n",
            k, PRIME_SIEVE_LIMIT);
        pfrac_fatal("aborting");
    }
    for (i = 0; i < g_nprimes && g_primes[i] <= k; i++) {
        pk = g_primes[i];
        while (k % pk == 0) {
            f->exp[i]++;
            last = i;
            k /= pk;
        }
        if (k == 1) break;
    }
    /* Unreachable when k <= PRIME_SIEVE_LIMIT, but kept as a defensive check. */
    if (k != 1)
        pfrac_fatal("pfrac_mul_int: residual factor after full factorization (internal error)");
    if (last + 1 > f->max_idx) f->max_idx = last + 1;
}

void pfrac_to_sqrt_rational(const pfrac_t *f,
                             bigint_t *int_num,  bigint_t *int_den,
                             bigint_t *sqrt_num, bigint_t *sqrt_den)
{
    int i, e, ae, half;
    for (i = 0; i < f->max_idx; i++) {
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

void pfrac_to_sqrt_rational_ws(const pfrac_t *f,
                                bigint_t *int_num,  bigint_t *int_den,
                                bigint_t *sqrt_num, bigint_t *sqrt_den,
                                bigint_ws_t *ws)
{
    int i, e, ae, half;
    for (i = 0; i < f->max_idx; i++) {
        e = f->exp[i];
        if (e == 0) continue;
        ae   = (e > 0) ? e : -e;
        half = ae / 2;
        if (e > 0) {
            if (half > 0) bigint_mul_prime_pow_ws(int_num,  (uint64_t)g_primes[i], half, ws);
            if (ae & 1)   bigint_mul_prime_pow_ws(sqrt_num, (uint64_t)g_primes[i], 1,    ws);
        } else {
            if (half > 0) bigint_mul_prime_pow_ws(int_den,  (uint64_t)g_primes[i], half, ws);
            if (ae & 1)   bigint_mul_prime_pow_ws(sqrt_den, (uint64_t)g_primes[i], 1,    ws);
        }
    }
}
