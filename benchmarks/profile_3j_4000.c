/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Single-purpose profile driver: many wigner3j(j=4000, j=4000, j=4000, m, -m, 0)
 * evaluations.  Used with perf record / callgrind to identify the
 * dominant cost path of libwignernj at large j, where it lags WIGXJPF
 * by about two orders of magnitude.
 */
#include "wigner.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    int N = (argc > 1) ? atoi(argv[1]) : 30;
    int j = (argc > 2) ? atoi(argv[2]) : 4000;
    int tj = 2 * j;
    int i;
    volatile double acc = 0.0;
    wigner_warmup();
    for (i = 0; i < N; i++) {
        int m = (i % (j + 1)) - j / 2;
        acc += wigner3j(tj, tj, tj, 2 * m, -2 * m, 0);
    }
    printf("acc=%.20e\n", (double)acc);
    return 0;
}
