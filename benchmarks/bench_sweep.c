/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Logarithmic-spaced bench sweep for one symbol family.  Prints
 *   j  ns/call
 * one line per j, flushing after each row so a sweep can be
 * interrupted with Ctrl-C and the partial data is still on disk.
 *
 * Usage:
 *   ./bench_sweep 3j 4000 0.1   # symbol family, j_max, budget seconds per j
 *
 * j sweep: stride 1 up to j=30, then geometric factor 1.15.
 */
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

static double bench3j(int j, int N)
{
    int i;
    int tj = 2 * j;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++) {
        int m = (i % (j + 1)) - j / 2;
        acc += wigner3j(tj, tj, tj, 2 * m, -2 * m, 0);
    }
    (void)acc;
    return (now() - t0) / N;
}

static double bench6j(int j, int N)
{
    int i;
    int tj = 2 * j;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++)
        acc += wigner6j(tj, tj, tj, tj, tj, tj);
    (void)acc;
    return (now() - t0) / N;
}

static double bench9j(int j, int N)
{
    int i;
    int tj = 2 * j;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++)
        acc += wigner9j(tj, tj, tj, tj, tj, tj, tj, tj, tj);
    (void)acc;
    return (now() - t0) / N;
}

static double benchgaunt(int j, int N)
{
    /* Diagonal Gaunt sweep at l1=l2=l3=l with m=(0,0,0).  Gaunt
     * vanishes when l1+l2+l3 is odd, so map the sweep parameter
     * j to l = 2j so the selection rule is always satisfied; the
     * inner Racah sum range is then n_terms = l+1 = 2j+1, which
     * exercises the Pass-2 ratio recurrence comparably to the
     * 3j(j,j,j;0,0,0) sweep. */
    int i;
    int tl = 4 * j;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++)
        acc += gaunt(tl, 0, tl, 0, tl, 0);
    (void)acc;
    return (now() - t0) / N;
}

int main(int argc, char **argv)
{
    const char *sym = (argc > 1) ? argv[1] : "3j";
    int j_max      = (argc > 2) ? atoi(argv[2]) : 4000;
    double budget  = (argc > 3) ? atof(argv[3]) : 0.1;
    double single_call_ceiling = (argc > 4) ? atof(argv[4]) : 30.0;

    double (*bench)(int, int);
    if (strcmp(sym, "3j") == 0)         bench = bench3j;
    else if (strcmp(sym, "6j") == 0)    bench = bench6j;
    else if (strcmp(sym, "9j") == 0)    bench = bench9j;
    else if (strcmp(sym, "gaunt") == 0) bench = benchgaunt;
    else { fprintf(stderr, "symbol must be 3j, 6j, 9j, or gaunt\n"); return 1; }

    wignernj_warmup_to(0);

    /* Build j list. */
    int j_list[5000];
    int n = 0;
    int j;
    for (j = 1; j <= 30 && j <= j_max && n < 5000; j++)
        j_list[n++] = j;
    double next = (n > 0 ? j_list[n-1] : 1) * 1.15;
    while (n < 5000) {
        int jn = (int)next;
        if (jn > j_max) break;
        if (jn > j_list[n-1]) j_list[n++] = jn;
        next *= 1.15;
    }

    fprintf(stderr, "# %s sweep: %d j values, j_max=%d, budget=%g s, ceiling=%g s\n",
            sym, n, j_max, budget, single_call_ceiling);
    printf("# %s   j   ns/call\n", sym);
    fflush(stdout);

    int i;
    for (i = 0; i < n; i++) {
        j = j_list[i];

        /* Calibrate N from a single-call probe. */
        double t1 = bench(j, 1);
        if (t1 > single_call_ceiling) {
            fprintf(stderr, "# stop at j=%d (single call %.3g s exceeds ceiling)\n",
                    j, t1);
            break;
        }
        int N = (int)(budget / t1);
        if (N < 3) N = 3;
        if (N > 1000000) N = 1000000;

        /* Min over 3 trials. */
        double tmin = 1e18;
        int trial;
        for (trial = 0; trial < 3; trial++) {
            double t = bench(j, N);
            if (t < tmin) tmin = t;
        }
        printf("%d %.1f\n", j, tmin * 1e9);
        fflush(stdout);
    }
    return 0;
}
