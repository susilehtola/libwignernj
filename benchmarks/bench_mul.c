/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Microbench for the multiword-bigint multiplication primitive.
 * Sweeps m (operand size in 64-bit limbs) and reports ns/call for
 * symmetric m × m multiplication.  Used to (i) characterise the
 * Karatsuba constant factor and asymptotic exponent on this machine
 * and (ii) decide whether Toom-3 would close the gap to FLINT/GMP at
 * m ≳ 100.
 *
 * Build:
 *   gcc -O3 -DNDEBUG -Iinclude -Isrc benchmarks/bench_mul.c \
 *       build/libwignernj.a -o bench_mul
 */
#include "bigint.h"
#include "scratch.h"
#include "wignernj.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Fill a bigint with m limbs of pseudo-random data. */
static void fill_random(bigint_t *r, size_t m, uint64_t seed)
{
    bigint_reserve(r, m);
    uint64_t x = seed;
    for (size_t i = 0; i < m; i++) {
        /* SplitMix64 step */
        x += 0x9E3779B97F4A7C15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z =  z ^ (z >> 31);
        r->words[i] = z;
    }
    /* ensure top limb non-zero so .size = m */
    if (r->words[m - 1] == 0) r->words[m - 1] = 1;
    r->size = m;
}

static double bench_mul(size_t m, int N, bigint_ws_t *ws)
{
    bigint_t a, b, r;
    bigint_init(&a); bigint_init(&b); bigint_init(&r);
    fill_random(&a, m, 0xDEADBEEF12345678ULL);
    fill_random(&b, m, 0xCAFEBABE87654321ULL);
    bigint_reserve(&r, 2 * m);

    /* Warm-up to size r and ensure ws scratch is sized. */
    bigint_mul_ws(&r, &a, &b, ws);

    double t0 = now();
    for (int i = 0; i < N; i++) {
        /* Vary one bit per iteration so the compiler can't hoist the
         * call (the operand changes between iterations) and so the
         * computation is genuinely repeated. */
        a.words[0] ^= (uint64_t)i;
        bigint_mul_ws(&r, &a, &b, ws);
        a.words[0] ^= (uint64_t)i;
    }
    double t = (now() - t0) / N;
    bigint_free(&a); bigint_free(&b); bigint_free(&r);
    return t;
}

int main(int argc, char **argv)
{
    int trials = (argc > 1) ? atoi(argv[1]) : 5;
    bigint_ws_t ws;
    bigint_ws_init(&ws);
    bigint_ws_reserve(&ws, 8192);

    /* Sweep m over realistic libwignernj sizes:
     *   m=8    typical at j ~ 50
     *   m=16   typical at j ~ 200
     *   m=32   typical at j ~ 500     (Karatsuba threshold is 32)
     *   m=64   typical at j ~ 1000
     *   m=128  typical at j ~ 2000
     *   m=256  typical at j ~ 4000
     *   m=512  beyond library bound, asymptotic check
     *   m=1024 ditto */
    size_t ms[] = { 4, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024 };
    int Ns[]   = { 2000000, 2000000, 1000000, 500000, 500000,
                   200000,  200000,   100000, 100000, 50000,
                    50000,   20000,    10000,   5000,   2000 };

    printf("%-6s  %12s  %10s  %10s\n", "m", "ns/call", "ns/limb²", "regime");
    for (size_t i = 0; i < sizeof(ms)/sizeof(ms[0]); i++) {
        double tmin = 1e18;
        for (int t = 0; t < trials; t++) {
            double tt = bench_mul(ms[i], Ns[i], &ws);
            if (tt < tmin) tmin = tt;
        }
        double per_limb_sq = tmin * 1e9 / (double)(ms[i] * ms[i]);
        const char *regime = (ms[i] < 32) ? "schoolbook" : "Karatsuba";
        printf("%-6zu  %12.1f  %10.4f  %10s\n",
               ms[i], tmin * 1e9, per_limb_sq, regime);
    }
    bigint_ws_free(&ws);
    return 0;
}
