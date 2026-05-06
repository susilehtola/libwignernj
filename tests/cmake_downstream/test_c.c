/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Smoke test for the C interface of an installed libwignernj. */
#include <math.h>
#include <stdio.h>
#include "wignernj.h"

int main(void)
{
    /* wigner3j(1,1,0; 0,0,0) = -1/sqrt(3) */
    double v   = wigner3j(2, 2, 0,  0, 0, 0);
    double ref = -1.0 / sqrt(3.0);
    if (fabs(v - ref) > 1e-14) {
        fprintf(stderr, "C: wigner3j returned %.17g, expected %.17g\n", v, ref);
        return 1;
    }
    printf("C: wigner3j(1,1,0;0,0,0) = %.15g  (expected %.15g)\n", v, ref);
    return 0;
}
