/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Microbenchmark for bigint_div128: cycles per call across all three
 * dispatch arms.  Pulls the static-inline header in directly so it can
 * be timed without the overhead of bigint_div_u64 wrapper iteration.
 *
 * Build:  gcc -O3 -I src benchmarks/bench_div128.c -o bench_div128 \
 *             [extra: -DBIGINT_NO_DIVQ | -DBIGINT_FORCE_PORTABLE]
 */
#include "bigint_arith.h"
#include <stdio.h>
#include <stdint.h>
#include <time.h>

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static double bench_div(uint64_t hi_seed, uint64_t lo_seed, uint64_t d, int N)
{
    int i;
    volatile uint64_t acc_q = 0, acc_r = 0;
    /* Mix the seed each iteration to keep hi < d (the bigint precondition)
     * and prevent the compiler from hoisting the call out. */
    double t0 = now();
    for (i = 0; i < N; i++) {
        uint64_t r;
        uint64_t lo = lo_seed ^ ((uint64_t)i * 0x9E3779B97F4A7C15ULL);
        uint64_t hi = hi_seed % d;            /* ensure hi < d */
        uint64_t q  = bigint_div128(hi, lo, d, &r);
        acc_q ^= q;
        acc_r ^= r;
    }
    (void)acc_q; (void)acc_r;
    return (now() - t0) / N;
}

int main(void)
{
    const int N = 5000000;

    struct { uint64_t d; const char *label; } cases[] = {
        { 7,                       "d = 7              (32-bit fast path)" },
        { 0x12345,                 "d = 75077          (32-bit fast path)" },
        { 0xFFFFFFFFULL,           "d = 2^32 - 1       (32-bit fast path)" },
        { 0x100000007ULL,          "d = 2^32 + 7       (Algorithm D)     " },
        { 0xDEADBEEFCAFEBABEULL,   "d = 0xDEADBEEFCAFEBABE (Algorithm D)" },
        { 0xFFFFFFFFFFFFFFFFULL,   "d = 2^64 - 1       (Algorithm D + sat)" },
    };
    size_t k;

    printf("%-50s   ns/call\n", "case");
    printf("%-50s   -------\n", "");
    for (k = 0; k < sizeof(cases)/sizeof(cases[0]); k++) {
        double tmin = 1e18;
        int trial;
        for (trial = 0; trial < 5; trial++) {
            double t = bench_div(0xDEADBEEF12345678ULL, 0xCAFEBABE87654321ULL,
                                  cases[k].d, N);
            if (t < tmin) tmin = t;
        }
        printf("%-50s   %7.2f\n", cases[k].label, tmin * 1e9);
    }
    return 0;
}
