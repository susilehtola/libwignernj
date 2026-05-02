/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Comparative benchmark of libwignernj against WIGXJPF and the GNU
 * Scientific Library on a sweep of "all-equal-j" inputs:
 *
 *   3j(j,j,j; m,-m,0)         with m sweeping over -j..+j
 *   6j{j j j; j j j}          (single fixed input per j)
 *   9j{j j j; j j j; j j j}   (single fixed input per j)
 *
 * Each measurement is repeated BENCH_REPS times consecutively; the
 * minimum and median per-evaluation wall times are reported, together
 * with the relative spread (max-min)/min as a stability indicator.  The
 * minimum is the conventional figure of merit for microbenchmarks
 * because it represents the run least disturbed by scheduler/cache/IRQ
 * noise.  The floating-point sum of the returned values is also reported
 * as a sanity check that the libraries agree on the numerical result.
 * GSL's gamma-function implementation overflows at moderate j; the
 * resulting GSL error is caught by a custom error handler and reported
 * as "GSL OVERFLOW".
 *
 * Build:  see the Makefile in this directory.
 *
 * Run on a quiet machine with the CPU frequency governor pinned to
 * "performance" for best repeatability.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "wigner.h"
#include "wigxjpf.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_errno.h>

#ifndef BENCH_REPS
#define BENCH_REPS 50   /* number of consecutive repetitions per measurement */
#endif

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* GSL's default error handler aborts the program; catch overflows here. */
static int g_gsl_overflow = 0;
static void gsl_err_handler(const char *r, const char *f, int l, int e)
{
    (void)r; (void)f; (void)l; (void)e;
    g_gsl_overflow = 1;
}

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/* Internal: invoke the BENCH macro body REPS times, collect per-eval ns,
 * and print {min, median, spread%}. */
#define BENCH(label, expr, N) do {                                       \
    double ns[BENCH_REPS];                                               \
    double s = 0.0;                                                      \
    for (int rep = 0; rep < BENCH_REPS; rep++) {                         \
        double t0 = now();                                               \
        for (long long i = 0; i < (long long)(N); i++) s += (expr);      \
        double t1 = now();                                               \
        ns[rep] = 1e9 * (t1 - t0) / (double)(N);                         \
    }                                                                    \
    qsort(ns, BENCH_REPS, sizeof(double), cmp_double);                   \
    double tmin = ns[0];                                                 \
    double tmed = ns[BENCH_REPS / 2];                                    \
    double tmax = ns[BENCH_REPS - 1];                                    \
    double spread = (tmax - tmin) / tmin * 100.0;                        \
    printf("  %-26s  min %10.0f  med %10.0f ns/eval"                     \
           "  spread %5.1f%%  (sum=%.6g)\n",                              \
           label, tmin, tmed, spread, s);                                \
} while (0)

#define BENCH_GSL(label, expr, N) do {                                   \
    g_gsl_overflow = 0;                                                  \
    double ns[BENCH_REPS];                                               \
    int reps_done = 0;                                                   \
    double s = 0.0;                                                      \
    for (int rep = 0; rep < BENCH_REPS && !g_gsl_overflow; rep++) {      \
        double t0 = now();                                               \
        for (long long i = 0; i < (long long)(N) && !g_gsl_overflow; i++)\
            s += (expr);                                                 \
        double t1 = now();                                               \
        if (g_gsl_overflow) break;                                       \
        ns[reps_done++] = 1e9 * (t1 - t0) / (double)(N);                 \
    }                                                                    \
    if (g_gsl_overflow) {                                                \
        printf("  %-26s  GSL OVERFLOW\n", label);                        \
    } else {                                                             \
        qsort(ns, reps_done, sizeof(double), cmp_double);                \
        double tmin = ns[0];                                             \
        double tmed = ns[reps_done / 2];                                 \
        double tmax = ns[reps_done - 1];                                 \
        double spread = (tmax - tmin) / tmin * 100.0;                    \
        printf("  %-26s  min %10.0f  med %10.0f ns/eval"                 \
               "  spread %5.1f%%  (sum=%.6g)\n",                          \
               label, tmin, tmed, spread, s);                            \
    }                                                                    \
} while (0)

int main(void)
{
    gsl_set_error_handler(gsl_err_handler);

    /* WIGXJPF requires the maximum 2j across all calls to be set up
     * before any evaluation.  Allocate generously. */
    const int max_two_j = 2 * 8000;
    wig_table_init(max_two_j, 9);
    wig_temp_init (max_two_j);

    /* j values and iteration counts (chosen so that each timing runs
     * for at least a few milliseconds, which is enough to outrun
     * clock-resolution noise). */
    int Js[]    = {10, 30, 60, 200, 1000, 2000, 4000, 8000};
    int N_3j[]  = {200000, 50000, 10000,  500,   50,   10,    5,    2};
    int N_6j[]  = {100000, 20000,  5000,  500,   50,   10,    1,    1};

    /* ── 3j sweep ─────────────────────────────────────────────────── */
    printf("=== 3j(j,j,j; m,-m,0)  --  ns/eval ===\n");
    for (size_t k = 0; k < sizeof(Js) / sizeof(Js[0]); k++) {
        int TJ = Js[k];
        int N  = N_3j[k];
        char Llib[40], Lwig[40], Lgsl[40];
        snprintf(Llib, sizeof Llib, "libwignernj j=%d", TJ / 2);
        snprintf(Lwig, sizeof Lwig, "WIGXJPF     j=%d", TJ / 2);
        snprintf(Lgsl, sizeof Lgsl, "GSL         j=%d", TJ / 2);
        BENCH    (Llib, wigner3j           (TJ, TJ, TJ,
                  2 * ((int)(i % (TJ + 1))) - TJ,
                  TJ - 2 * ((int)(i % (TJ + 1))), 0), N);
        BENCH    (Lwig, wig3jj             (TJ, TJ, TJ,
                  2 * ((int)(i % (TJ + 1))) - TJ,
                  TJ - 2 * ((int)(i % (TJ + 1))), 0), N);
        BENCH_GSL(Lgsl, gsl_sf_coupling_3j (TJ, TJ, TJ,
                  2 * ((int)(i % (TJ + 1))) - TJ,
                  TJ - 2 * ((int)(i % (TJ + 1))), 0), N);
    }

    /* ── 6j sweep ─────────────────────────────────────────────────── */
    printf("\n=== 6j{j,j,j; j,j,j}  --  ns/eval ===\n");
    for (size_t k = 0; k < sizeof(Js) / sizeof(Js[0]); k++) {
        int TJ = Js[k];
        int N  = N_6j[k];
        char Llib[40], Lwig[40], Lgsl[40];
        snprintf(Llib, sizeof Llib, "libwignernj j=%d", TJ / 2);
        snprintf(Lwig, sizeof Lwig, "WIGXJPF     j=%d", TJ / 2);
        snprintf(Lgsl, sizeof Lgsl, "GSL         j=%d", TJ / 2);
        BENCH    (Llib, wigner6j           (TJ,TJ,TJ, TJ,TJ,TJ), N);
        BENCH    (Lwig, wig6jj             (TJ,TJ,TJ, TJ,TJ,TJ), N);
        BENCH_GSL(Lgsl, gsl_sf_coupling_6j (TJ,TJ,TJ, TJ,TJ,TJ), N);
    }

    /* ── 9j sweep ─────────────────────────────────────────────────── */
    printf("\n=== 9j{j,j,j; j,j,j; j,j,j}  --  ns/eval ===\n");
    int Js9[] = {6, 10, 20, 80, 160};
    int N_9j[] = {1000, 100, 20, 1, 1};
    for (size_t k = 0; k < sizeof(Js9) / sizeof(Js9[0]); k++) {
        int TJ = Js9[k];
        int N  = N_9j[k];
        char Llib[40], Lwig[40], Lgsl[40];
        snprintf(Llib, sizeof Llib, "libwignernj j=%d", TJ / 2);
        snprintf(Lwig, sizeof Lwig, "WIGXJPF     j=%d", TJ / 2);
        snprintf(Lgsl, sizeof Lgsl, "GSL         j=%d", TJ / 2);
        BENCH    (Llib, wigner9j           (TJ,TJ,TJ, TJ,TJ,TJ, TJ,TJ,TJ), N);
        BENCH    (Lwig, wig9jj             (TJ,TJ,TJ, TJ,TJ,TJ, TJ,TJ,TJ), N);
        BENCH_GSL(Lgsl, gsl_sf_coupling_9j (TJ,TJ,TJ, TJ,TJ,TJ, TJ,TJ,TJ), N);
    }

    wig_temp_free();
    wig_table_free();
    return 0;
}
