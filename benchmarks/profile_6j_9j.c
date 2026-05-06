/* SPDX-License-Identifier: BSD-3-Clause */
#include "wigner.h"
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
    int N = (argc > 1) ? atoi(argv[1]) : 30;
    int j = (argc > 2) ? atoi(argv[2]) : 4000;
    char k = (argc > 3) ? argv[3][0] : '6';
    int tj = 2*j, i;
    volatile double acc = 0;
    wigner_warmup();
    for (i = 0; i < N; i++) {
        if (k == '6') acc += wigner6j(tj, tj, tj, tj, tj, tj);
        else          acc += wigner9j(tj, tj, tj, tj, tj, tj, tj, tj, tj);
    }
    printf("%.20e\n", (double)acc);
    return 0;
}
