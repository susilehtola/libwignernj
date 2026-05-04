/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "pfrac.h"
#include "wigner.h"
#include "wignernj_tls.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Public accessor: absolute compile-time ceiling on factorial arguments.
 * Determined by the prime sieve in src/prime_table.inc. */
int wigner_max_factorial_arg(void) { return MAX_FACTORIAL_ARG; }

/* ── factorial-decomposition cache ────────────────────────────────────────────
 *
 * pfrac_mul_factorial / pfrac_div_factorial used to call
 * legendre_valuation(N, p_i) for every prime p_i <= N on every entry,
 * paying a per-call cost proportional to the number of primes below N
 * times the per-Legendre-valuation divide loop.  Because the prime
 * decomposition of N! depends only on N, this is wasteful: at the
 * scale of a single 3j evaluation the same N is decomposed dozens of
 * times across the Racah sum.
 *
 * The cache is a per-thread table indexed by N: g_fact_cache[N], when
 * non-NULL, points to a malloc'd array of fact_width(N) ints holding
 * the exponent of each g_primes[i] (i = 0 .. fact_width(N)-1) in N!.
 * pfrac_mul_factorial / pfrac_div_factorial are now a vector add /
 * subtract over that row instead of a per-prime divide loop.
 *
 * Lifecycle:  lazy-allocated on first use, persists for the lifetime
 * of the calling thread, never freed (the OS reclaims at thread exit).
 * Memory grows with the largest N that thread sees; at the absolute
 * default-build ceiling (MAX_FACTORIAL_ARG = 20000) the per-thread
 * cost is bounded by sum_{N=1}^{20000} fact_width(N) ints ~ 80 MB.
 * For typical workloads (j <~ a few thousand) it is a small fraction
 * of that.
 *
 * Toolchains without thread-local storage skip the cache entirely and
 * fall back to the per-call legendre_valuation loop -- still
 * thread-safe (no shared state), just slower.
 */

#if WIGNERNJ_HAVE_TLS

/* Number of primes <= n in g_primes.  O(1) lookup into the
 * compile-time table generated alongside g_primes.  Only used by the
 * TLS-cached factorial-decomposition cache; the no-TLS fallback walks
 * the prime list directly via g_primes[i] <= n. */
static int fact_width(int n)
{
    return (int)g_pi_table[n];
}

/* Pointer table: indexed by N, lazy-allocated rows.  The table itself
 * is allocated lazily on first cache access. */
static WIGNERNJ_TLS int **g_fact_cache = NULL;

static void fact_cache_ensure_table(void)
{
    if (!g_fact_cache)
        g_fact_cache =
            (int **)xcalloc((size_t)(MAX_FACTORIAL_ARG + 1), sizeof(int *));
}

static const int *fact_cache_get(int n)
{
    fact_cache_ensure_table();
    if (g_fact_cache[n]) return g_fact_cache[n];

    int width = fact_width(n);
    int *row  = (int *)xmalloc((size_t)width * sizeof(int));
    int i;
    for (i = 0; i < width; i++)
        row[i] = legendre_valuation(n, i);
    g_fact_cache[n] = row;
    return row;
}

void wigner_warmup_factorial_cache(int N_max)
{
    int n;
    if (N_max < 1) return;
    if (N_max > MAX_FACTORIAL_ARG) N_max = MAX_FACTORIAL_ARG;
    fact_cache_ensure_table();
    for (n = 2; n <= N_max; n++) {
        if (!g_fact_cache[n])
            (void)fact_cache_get(n);
    }
}

/* Release the calling thread's factorial-decomposition cache.  Frees
 * every cached row plus the pointer table itself.  After return, the
 * next pfrac_mul_factorial / pfrac_div_factorial call in this thread
 * starts fresh (re-allocates the table on first use, rebuilds rows
 * on demand). */
void wigner_factorial_cache_release(void)
{
    int n;
    if (!g_fact_cache) return;
    for (n = 0; n <= MAX_FACTORIAL_ARG; n++) {
        free(g_fact_cache[n]);
    }
    free(g_fact_cache);
    g_fact_cache = NULL;
}

#else  /* !WIGNERNJ_HAVE_TLS */

void wigner_warmup_factorial_cache(int N_max)
{
    /* No persistent cache to populate -- the no-TLS path already
     * computes legendre_valuation per call.  No-op for source
     * compatibility with the cached path. */
    (void)N_max;
}

void wigner_factorial_cache_release(void)
{
    /* No persistent cache; nothing to free. */
}

#endif

/* ─────────────────────────────────────────────────────────────────────────── */

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
    int i, width;
    if (n <= 1) return;
    if (n > MAX_FACTORIAL_ARG) {
        fprintf(stderr,
            "libwignernj: factorial argument %d exceeds MAX_FACTORIAL_ARG=%d "
            "(angular momenta too large for prime table)\n",
            n, MAX_FACTORIAL_ARG);
        pfrac_fatal("aborting");
    }
#if WIGNERNJ_HAVE_TLS
    {
        const int *row = fact_cache_get(n);
        width = fact_width(n);
        for (i = 0; i < width; i++) f->exp[i] += row[i];
    }
#else
    for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
        f->exp[i] += legendre_valuation(n, i);
    width = i;
#endif
    if (width > f->max_idx) f->max_idx = width;
}

void pfrac_div_factorial(pfrac_t *f, int n)
{
    int i, width;
    if (n <= 1) return;
    if (n > MAX_FACTORIAL_ARG) {
        fprintf(stderr,
            "libwignernj: factorial argument %d exceeds MAX_FACTORIAL_ARG=%d "
            "(angular momenta too large for prime table)\n",
            n, MAX_FACTORIAL_ARG);
        pfrac_fatal("aborting");
    }
#if WIGNERNJ_HAVE_TLS
    {
        const int *row = fact_cache_get(n);
        width = fact_width(n);
        for (i = 0; i < width; i++) f->exp[i] -= row[i];
    }
#else
    for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
        f->exp[i] -= legendre_valuation(n, i);
    width = i;
#endif
    if (width > f->max_idx) f->max_idx = width;
}

void pfrac_mul_int(pfrac_t *f, int k)
{
    int i;
    int last = -1;
    if (k <= 1) return;
    if (k > PRIME_SIEVE_LIMIT) {
        fprintf(stderr,
            "libwignernj: pfrac_mul_int argument %d exceeds PRIME_SIEVE_LIMIT=%d "
            "(angular momenta too large for prime table)\n",
            k, PRIME_SIEVE_LIMIT);
        pfrac_fatal("aborting");
    }
    /*
     * Trial-divide k by primes up to floor(sqrt(k)).  When the loop
     * exits with k > 1, the residual factor must itself be prime (it
     * has no divisor below sqrt(k_original)), so its index can be
     * looked up in O(1) via prime_index_of() rather than continuing
     * the trial-division loop.  For prime or near-prime k this turns
     * an O(pi(k)) sweep into an O(pi(sqrt(k))) one -- e.g. for
     * k = 1031 (prime, just above 32^2), 173 prime tests become 11.
     *
     * The bound check uses the *current* k (which decreases as we
     * divide it), so smooth k -- which collapse to 1 quickly -- exit
     * the loop as soon as no remaining factor is small enough to
     * matter, matching the speed of the original tight loop.
     */
    for (i = 0; i < g_nprimes; i++) {
        long long p = g_primes[i];
        if (p * p > (long long)k) break;
        while (k % g_primes[i] == 0) {
            f->exp[i]++;
            last = i;
            k /= g_primes[i];
        }
        if (k == 1) break;
    }
    if (k > 1) {
        int idx = prime_index_of(k);
        if (idx < 0)
            pfrac_fatal("pfrac_mul_int: residual factor not in prime table "
                        "(internal error)");
        f->exp[idx]++;
        if (idx > last) last = idx;
    }
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
