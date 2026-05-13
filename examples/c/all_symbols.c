/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * libwignernj C API demonstration.
 *
 * Calls every public symbol family exposed by wignernj.h once with small,
 * textbook-scale arguments and prints the result alongside an analytic
 * reference value.  The program exits 0 on success, non-zero if any
 * computed value disagrees with its reference by more than 1e-14.
 *
 * All angular-momentum arguments are passed as 2*j integers (so
 * j = 1/2 is encoded as tj = 1, j = 1 as tj = 2, etc.); this avoids
 * floating-point representation of half-integer quantum numbers.
 *
 * Build (out-of-tree, against an installed libwignernj):
 *     cc -o all_symbols_c all_symbols.c -lwignernj -lm
 */
#include "wignernj.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define CHECK(label, computed, expected, tol) do {                       \
    double _c = (computed), _e = (expected);                             \
    printf("  %-38s = %+.15f   (expected %+.15f)\n", label, _c, _e);    \
    if (fabs(_c - _e) > (tol)) {                                         \
        fprintf(stderr,                                                   \
                "  FAIL: |%g - %g| = %g exceeds tolerance %g\n",         \
                _c, _e, fabs(_c - _e), (tol));                           \
        return 1;                                                        \
    }                                                                    \
} while (0)

int main(void)
{
    const double tol = 1e-14;

    printf("libwignernj C API demonstration -- one call per symbol family\n");
    printf("--------------------------------------------------------------\n");

    /* 1. Wigner 3j symbol.   ( 1 1 0 )
     *                        ( 0 0 0 )  =  -1/sqrt(3)                  */
    CHECK("wigner3j(1,1,0; 0,0,0)",
          wigner3j(2, 2, 0,  0, 0, 0),
          -1.0 / sqrt(3.0), tol);

    /* 2. Wigner 6j symbol.   { 1 1 1 }
     *                        { 1 1 1 }  =  1/6                          */
    CHECK("wigner6j{1,1,1; 1,1,1}",
          wigner6j(2, 2, 2,  2, 2, 2),
          1.0 / 6.0, tol);

    /* 3. Wigner 9j symbol.   { 1 1 0 }
     *                        { 1 1 0 }  =  1/3
     *                        { 0 0 0 }                                  */
    CHECK("wigner9j{1,1,0; 1,1,0; 0,0,0}",
          wigner9j(2, 2, 0,  2, 2, 0,  0, 0, 0),
          1.0 / 3.0, tol);

    /* 4. Clebsch-Gordan coefficient.   <1,0; 1,0 | 2,0> = sqrt(2/3)     */
    CHECK("clebsch_gordan(1,0; 1,0 | 2,0)",
          clebsch_gordan(2, 0,  2, 0,  4, 0),
          sqrt(2.0 / 3.0), tol);

    /* 5. Racah W coefficient.   W(1,1,1,1; 1,1) = 1/6                   */
    CHECK("racah_w(1,1,1,1; 1,1)",
          racah_w(2, 2, 2,  2, 2, 2),
          1.0 / 6.0, tol);

    /* 6. Fano X coefficient.   X(1,1,1; 1,1,1; 1,1,2) = 1/2             */
    CHECK("fano_x(1,1,1; 1,1,1; 1,1,2)",
          fano_x(2, 2, 2,  2, 2, 2,  2, 2, 4),
          0.5, tol);

    /* 7. Gaunt coefficient over complex Y_l^m.
     *    integral(Y_1^0 Y_1^0 Y_2^0) dΩ = 1/sqrt(5*pi)                  */
    CHECK("gaunt(1,0; 1,0; 2,0)",
          gaunt(2, 0,  2, 0,  4, 0),
          1.0 / sqrt(5.0 * M_PI), tol);

    /* 8. Gaunt coefficient over real spherical harmonics.
     *    For m1 = m2 = m3 = 0 the real and complex Gaunt coincide.      */
    CHECK("gaunt_real(1,0; 1,0; 2,0)",
          gaunt_real(2, 0,  2, 0,  4, 0),
          1.0 / sqrt(5.0 * M_PI), tol);

    /* 9. Real <-> complex Y_lm basis transformation.  At l = 1 the
     *    (m_r=+1, m_c=-1) entry is +1/sqrt(2) (real part), the
     *    (m_r=-1, m_c=-1) entry is +1/sqrt(2) (imag part), and the
     *    matrix is unitary.  Storage is column-major with leading
     *    dimension 2l+1; element C[r,c] sits at index c*(2l+1)+r. */
    {
        wignernj_cdouble_t C[9];
        wignernj_real_ylm_in_complex_ylm(1, C);
        double *F = (double *)C;            /* re/im pairs interleaved */
        /* row m_r = +1, col m_c = -1: real part = 1/sqrt(2). */
        CHECK("real_ylm_in_complex_ylm[+1,-1].re (l=1)",
              F[2 * (0 * 3 + 2)], 1.0 / sqrt(2.0), tol);
        /* row m_r = -1, col m_c = -1: imag part = 1/sqrt(2). */
        CHECK("real_ylm_in_complex_ylm[-1,-1].im (l=1)",
              F[2 * (0 * 3 + 0) + 1], 1.0 / sqrt(2.0), tol);
        /* row m_r = 0, col m_c = 0: real part = 1, imag part = 0. */
        CHECK("real_ylm_in_complex_ylm[ 0, 0].re (l=1)",
              F[2 * (1 * 3 + 1)], 1.0, tol);
    }

    printf("\nAll symbols agree with their analytic references.\n");
    return 0;
}
