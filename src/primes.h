/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef PRIMES_H
#define PRIMES_H

/*
 * Prime table and Legendre factorial valuation.
 *
 * Covers all primes up to PRIME_SIEVE_LIMIT, which is sufficient for
 * factorials up to MAX_FACTORIAL_ARG! as needed by Wigner symbol computations
 * with angular momenta up to ~j_max = 1000.
 */

#define PRIME_SIEVE_LIMIT   20011   /* sieve all primes up to this value */
#define MAX_PRIME_COUNT     2263    /* upper bound: pi(20011) = 2262        */
#define MAX_FACTORIAL_ARG   20000   /* maximum n for which n! may be needed */

/* Number of primes found (set by primes_init). */
extern int g_nprimes;

/* g_primes[i] = i-th prime (0-indexed: g_primes[0]=2, g_primes[1]=3, ...) */
extern int g_primes[MAX_PRIME_COUNT];

/*
 * g_prime_index[p] = index i such that g_primes[i] == p, or -1.
 * Valid for 0 <= p <= PRIME_SIEVE_LIMIT.
 */
extern short g_prime_index[PRIME_SIEVE_LIMIT + 1];

/*
 * Initialize prime table.  Safe to call multiple times (idempotent).
 * Called automatically by the first Wigner symbol computation.
 * Not thread-safe: call once from a single thread before spawning workers,
 * or accept the benign data race (writes are deterministic and idempotent).
 */
void primes_init(void);

/*
 * Return v_p(n!) = sum_{k>=1} floor(n / p^k), the p-adic valuation of n!.
 * pi is the index of p in g_primes.  n must satisfy 0 <= n <= MAX_FACTORIAL_ARG.
 */
int legendre_valuation(int n, int pi);

#endif /* PRIMES_H */
