/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "primes.h"
#include <string.h>

int   g_nprimes = 0;
int   g_primes[MAX_PRIME_COUNT];
short g_prime_index[PRIME_SIEVE_LIMIT + 1];

static int g_initialized = 0;

void primes_init(void)
{
    static unsigned char sieve[PRIME_SIEVE_LIMIT + 1];
    int i, j;

    if (g_initialized) return;
    g_initialized = 1;

    memset(sieve, 0, sizeof(sieve));
    memset(g_prime_index, -1, sizeof(g_prime_index));  /* -1 = not prime */

    sieve[0] = sieve[1] = 1; /* mark composites */
    for (i = 2; (long)i * i <= PRIME_SIEVE_LIMIT; i++) {
        if (!sieve[i]) {
            for (j = i * i; j <= PRIME_SIEVE_LIMIT; j += i)
                sieve[j] = 1;
        }
    }

    g_nprimes = 0;
    for (i = 2; i <= PRIME_SIEVE_LIMIT; i++) {
        if (!sieve[i]) {
            g_primes[g_nprimes] = i;
            g_prime_index[i]    = (short)g_nprimes;
            g_nprimes++;
        }
    }
}

int legendre_valuation(int n, int pi)
{
    /* v_p(n!) = sum_{k>=1} floor(n / p^k) */
    int p   = g_primes[pi];
    int val = 0;
    int pk  = p;
    while (pk <= n) {
        val += n / pk;
        if (pk > n / p) break;  /* prevent overflow before multiplication */
        pk *= p;
    }
    return val;
}
