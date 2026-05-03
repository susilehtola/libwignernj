/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Verify that wigner_warmup() pre-grows the calling thread's cached
 * scratch so that every subsequent symbol evaluation is allocation-
 * free, regardless of which symbol family is called or what angular
 * momenta are passed (within the prime-table ceiling).
 *
 * This test wraps glibc's malloc/calloc/realloc family via a counter-
 * pair (start/stop_counting), counts the libwignernj-internal
 * allocations across a representative call mix, and asserts the count
 * is exactly zero.
 */
#include "run_tests.h"
#include "../include/wigner.h"

#if defined(__APPLE__) || defined(_WIN32)
/* dlsym/RTLD_NEXT-based interception is fragile on macOS and absent on
 * Windows; the warmup behaviour is exercised on Linux only, which is
 * sufficient for CI coverage. */
int main(void)
{
    printf("test_warmup: SKIP (allocation interception requires GNU libc)\n");
    return 0;
}
#else

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

static int  g_counting = 0;
static long g_n_alloc  = 0;

void *malloc(size_t n) {
    static void *(*real)(size_t);
    if (!real) real = dlsym(RTLD_NEXT, "malloc");
    if (g_counting) __atomic_add_fetch(&g_n_alloc, 1, __ATOMIC_RELAXED);
    return real(n);
}
void *calloc(size_t n, size_t s) {
    static void *(*real)(size_t, size_t);
    if (!real) real = dlsym(RTLD_NEXT, "calloc");
    if (g_counting) __atomic_add_fetch(&g_n_alloc, 1, __ATOMIC_RELAXED);
    return real(n, s);
}
void *realloc(void *p, size_t n) {
    static void *(*real)(void *, size_t);
    if (!real) real = dlsym(RTLD_NEXT, "realloc");
    if (g_counting) __atomic_add_fetch(&g_n_alloc, 1, __ATOMIC_RELAXED);
    return real(p, n);
}

int main(void)
{
    /* The warmup is a no-op on a toolchain without thread-local storage
     * (every call there already allocates fresh).  In that case the
     * zero-allocations assertion below would always fail, so skip it
     * with a clear note. */
    if (!wigner_thread_local_scratch_available()) {
        printf("test_warmup: SKIP (no thread-local storage on this build)\n");
        return 0;
    }

    /* Warm up both the cached scratch and the factorial-decomposition
     * cache.  We size the latter precisely from the call mix below
     * via the per-symbol max_factorial helpers, taking the maximum so
     * that every factorial reached by the test is pre-populated. */
    wigner_warmup();
    {
        int N = 0, n;
        n = wigner3j_max_factorial(2,    2,    2,    0, 0, 0);    if (n > N) N = n;
        n = wigner3j_max_factorial(200,  200,  200,  2, 2,-4);    if (n > N) N = n;
        n = wigner3j_max_factorial(13332,13332,13332,2, 2,-4);    if (n > N) N = n;
        n = wigner6j_max_factorial(2,   2,   2,   2,   2,   2);   if (n > N) N = n;
        n = wigner6j_max_factorial(200, 200, 200, 200, 200, 200); if (n > N) N = n;
        n = wigner9j_max_factorial(2, 2, 2, 2, 2, 2, 2, 2, 2);    if (n > N) N = n;
        n = wigner9j_max_factorial(20,20,20,20,20,20,20,20,20);   if (n > N) N = n;
        n = clebsch_gordan_max_factorial(2, 1, 2, -1, 0, 0);      if (n > N) N = n;
        n = racah_w_max_factorial       (2, 2, 2, 2, 2, 2);       if (n > N) N = n;
        n = fano_x_max_factorial        (2, 2, 2, 2, 2, 2,
                                         2, 2, 2);                if (n > N) N = n;
        n = gaunt_max_factorial         (4, 0, 4, 0, 8, 0);       if (n > N) N = n;
        n = gaunt_real_max_factorial    (4, 0, 4, 0, 8, 0);       if (n > N) N = n;
        wigner_warmup_factorial_cache(N);
    }

    /* Now count allocations across a representative call mix. */
    g_n_alloc = 0;
    g_counting = 1;

    /* 3j: small and at the upper-end of the supported range. */
    (void)wigner3j(2,    2,    2,     0, 0, 0);
    (void)wigner3j(200,  200,  200,   2, 2, -4);
    (void)wigner3j(13332,13332,13332, 2, 2, -4);

    /* 6j: small and at the practical upper end (avoid the absolute
     * ceiling at 2j=9998 which would take many minutes). */
    (void)wigner6j(2,   2,   2,   2,   2,   2  );
    (void)wigner6j(200, 200, 200, 200, 200, 200);

    /* 9j: small and moderate. */
    (void)wigner9j(2, 2, 2,  2, 2, 2,  2, 2, 2);
    (void)wigner9j(20,20,20, 20,20,20, 20,20,20);

    /* Derived symbols. */
    (void)clebsch_gordan(2, 1, 2, -1, 0, 0);
    (void)racah_w       (2, 2, 2, 2,  2,  2);
    (void)fano_x        (2, 2, 2, 2, 2, 2, 2, 2, 2);
    (void)gaunt         (4, 0, 4, 0, 8, 0);
    (void)gaunt_real    (4, 0, 4, 0, 8, 0);

    g_counting = 0;

    printf("test_warmup: %ld allocations across the representative mix\n",
           g_n_alloc);
    TEST_ASSERT(g_n_alloc == 0);

    /* wigner_thread_cleanup() must drop both the scratch and the
     * factorial cache, so a subsequent symbol evaluation re-enters
     * the lazy-init path and allocates again. */
    wigner_thread_cleanup();
    g_n_alloc = 0;
    g_counting = 1;
    (void)wigner3j(2, 2, 2, 0, 0, 0);
    g_counting = 0;
    printf("test_warmup: %ld allocations on the call after cleanup "
           "(expected > 0)\n", g_n_alloc);
    TEST_ASSERT(g_n_alloc > 0);

    return g_tests_failed ? 1 : 0;
}

#endif /* !__APPLE__ && !_WIN32 */
