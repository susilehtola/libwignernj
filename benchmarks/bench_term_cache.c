/* Quick microbench for the cached-term Racah-sum change.
 * Build:  gcc -O3 -I include -I src benchmarks/bench_term_cache.c \
 *             build/libwignernj.a -o bench_term_cache
 * Run:    ./bench_term_cache
 */
#include "wignernj.h"
#include <stdio.h>
#include <time.h>

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static double bench3j(int tj, int N)
{
    int i;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++) {
        int m = (i % (tj + 1)) - tj/2;
        acc += wigner3j(tj, tj, tj, 2*m, -2*m, 0);
    }
    (void)acc;
    return (now() - t0) / N;
}

static double bench6j(int tj, int N)
{
    int i;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++)
        acc += wigner6j(tj, tj, tj, tj, tj, tj);
    (void)acc;
    return (now() - t0) / N;
}

static double bench9j(int tj, int N)
{
    int i;
    volatile double acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++)
        acc += wigner9j(tj, tj, tj, tj, tj, tj, tj, tj, tj);
    (void)acc;
    return (now() - t0) / N;
}

int main(void)
{
    struct { int tj; int N3j; int N6j; int N9j; } cases[] = {
        {  2, 2000000, 1000000,  100000},
        { 10,  500000,  300000,    5000},
        { 30,  100000,   30000,     200},
        { 60,   30000,    5000,      20},
        {200,    3000,     200,       0},
    };
    size_t k;
    wignernj_warmup();
    printf("%-6s %12s %12s %12s\n", "j", "3j ns", "6j ns", "9j ns");
    for (k = 0; k < sizeof(cases)/sizeof(cases[0]); k++) {
        int tj = cases[k].tj;
        double t3 = 1e18, t6 = 1e18, t9 = 1e18;
        int trial;
        for (trial = 0; trial < 3; trial++) {
            if (cases[k].N3j) { double t = bench3j(tj, cases[k].N3j); if (t < t3) t3 = t; }
            if (cases[k].N6j) { double t = bench6j(tj, cases[k].N6j); if (t < t6) t6 = t; }
            if (cases[k].N9j) { double t = bench9j(tj, cases[k].N9j); if (t < t9) t9 = t; }
        }
        printf("%-6d", tj/2);
        if (cases[k].N3j) printf(" %12.1f", t3*1e9); else printf(" %12s", "-");
        if (cases[k].N6j) printf(" %12.1f", t6*1e9); else printf(" %12s", "-");
        if (cases[k].N9j) printf(" %12.1f", t9*1e9); else printf(" %12s", "-");
        printf("\n");
    }
    return 0;
}
