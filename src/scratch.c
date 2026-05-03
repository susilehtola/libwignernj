/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "scratch.h"
#include "primes.h"
#include "wigner.h"
#include "wignernj_tls.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>

/* ── shared helpers ──────────────────────────────────────────────────────── */

static void scratch_init(wigner_scratch_t *s)
{
    int i;
    bigint_ws_init(&s->ws);
    for (i = 0; i < WIGNER_SCRATCH_BIGINTS; i++)
        bigint_init(&s->bigints[i]);
    for (i = 0; i < WIGNER_SCRATCH_PFRACS; i++)
        pfrac_init(&s->pfracs[i]);
    for (i = 0; i < WIGNER_SCRATCH_LCMEXP; i++) {
        s->lcm_exp[i] = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
        s->lcm_max_dirty[i] = 0;
    }
    wigner_exact_init(&s->exact);
}

static void scratch_destroy(wigner_scratch_t *s)
{
    int i;
    bigint_ws_free(&s->ws);
    for (i = 0; i < WIGNER_SCRATCH_BIGINTS; i++)
        bigint_free(&s->bigints[i]);
    for (i = 0; i < WIGNER_SCRATCH_PFRACS; i++)
        pfrac_free(&s->pfracs[i]);
    for (i = 0; i < WIGNER_SCRATCH_LCMEXP; i++)
        free(s->lcm_exp[i]);
    /* Mark the cached exact as never-zero so wigner_exact_free actually
     * frees its bigints (the is_zero short-circuit would otherwise leak). */
    s->exact.is_zero = 0;
    wigner_exact_free(&s->exact);
}

void wigner_scratch_lcm_clear(wigner_scratch_t *s, int idx)
{
    if (s->lcm_max_dirty[idx] > 0)
        memset(s->lcm_exp[idx], 0,
               (size_t)s->lcm_max_dirty[idx] * sizeof(int));
    s->lcm_max_dirty[idx] = 0;
}

void wigner_scratch_lcm_dirty(wigner_scratch_t *s, int idx, int new_max)
{
    if (new_max > s->lcm_max_dirty[idx])
        s->lcm_max_dirty[idx] = new_max;
}

/* ── TLS-cached path (the common case) ────────────────────────────────────── */

#if WIGNERNJ_HAVE_TLS

static WIGNERNJ_TLS wigner_scratch_t *g_scratch = NULL;

wigner_scratch_t *wigner_scratch_acquire(void)
{
    wigner_scratch_t *s = g_scratch;
    if (s) return s;
    s = (wigner_scratch_t *)xcalloc(1, sizeof(*s));
    scratch_init(s);
    g_scratch = s;
    return s;
}

void wigner_scratch_relinquish(wigner_scratch_t *s)
{
    /* Cached path: nothing to do, the scratch persists for the
     * lifetime of the calling thread.  Argument unused. */
    (void)s;
}

void wigner_scratch_release(void)
{
    wigner_scratch_t *s = g_scratch;
    if (!s) return;
    scratch_destroy(s);
    free(s);
    g_scratch = NULL;
}

/* ── per-call fallback (no TLS available) ─────────────────────────────────── */

#else /* !WIGNERNJ_HAVE_TLS */

wigner_scratch_t *wigner_scratch_acquire(void)
{
    wigner_scratch_t *s = (wigner_scratch_t *)xcalloc(1, sizeof(*s));
    scratch_init(s);
    return s;
}

void wigner_scratch_relinquish(wigner_scratch_t *s)
{
    if (!s) return;
    scratch_destroy(s);
    free(s);
}

void wigner_scratch_release(void)
{
    /* No persistent state to drop.  Each acquire/relinquish pair is
     * already independent and thread-safe by construction. */
}

#endif /* WIGNERNJ_HAVE_TLS */

/* ── public warmup ──────────────────────────────────────────────────────── */

void wigner_warmup(void)
{
#if WIGNERNJ_HAVE_TLS
    /* Pre-grow every cached buffer to the absolute default-build maximum.
     * The dominant factor is mw = bigint_words_for_factorial(MAX_FACTORIAL_ARG),
     * which is the size of the largest factorial the prime table can index.
     * 9j and Gaunt size their long-lived bigints to mw_prod = 5*mw to cover
     * triple-product cross terms; we use that as the universal upper bound.
     *
     * After this call, every subsequent symbol evaluation in the calling
     * thread is allocation-free, regardless of the symbol family or input
     * angular momenta (within the prime-table ceiling). */
    wigner_scratch_t *s = wigner_scratch_acquire();
    size_t mw      = bigint_words_for_factorial(MAX_FACTORIAL_ARG);
    size_t mw_prod = 5 * mw;
    int i;

    bigint_ws_reserve(&s->ws, mw_prod);
    for (i = 0; i < WIGNER_SCRATCH_BIGINTS; i++)
        bigint_reserve(&s->bigints[i], mw_prod);
    bigint_reserve(&s->exact.sum,      mw_prod);
    bigint_reserve(&s->exact.int_num,  mw_prod);
    bigint_reserve(&s->exact.int_den,  mw_prod);
    bigint_reserve(&s->exact.sqrt_num, mw_prod);
    bigint_reserve(&s->exact.sqrt_den, mw_prod);
    /* lcm_exp arrays are already sized to g_nprimes ints in scratch_init,
     * which is the worst-case for the prime table.  Pfracs are similarly
     * sized once on the first pfrac_mul_factorial call. */

    wigner_scratch_relinquish(s);
#else
    /* No-TLS fallback: every public-API call already allocates a fresh
     * scratch and frees it on return, so there is no persistent state
     * for the warmup to populate.  The function is therefore a no-op on
     * toolchains without thread-local storage; correctness is unaffected,
     * but the caller will continue to pay the per-call allocation cost. */
#endif
}

int wigner_thread_local_scratch_available(void)
{
#if WIGNERNJ_HAVE_TLS
    return 1;
#else
    return 0;
#endif
}

/* Forward declaration -- implemented in src/pfrac.c.  Not exposed in
 * any internal header because it is a pure cleanup hook. */
extern void wigner_factorial_cache_release(void);

void wigner_thread_cleanup(void)
{
    /* Release every per-thread cache held by the calling thread:
     * the per-symbol scratch (bigint workspace, exact-output bigints,
     * pfrac scratch, lcm-exponent arrays) and the factorial-decomposition
     * cache (per-N prime-exponent rows + the pointer table).  No-op
     * on no-TLS toolchains. */
    wigner_scratch_release();
    wigner_factorial_cache_release();
}
