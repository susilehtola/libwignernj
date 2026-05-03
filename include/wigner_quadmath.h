/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * libquadmath / IEEE 754 binary128 ("__float128") interface to libwignernj.
 *
 * Built only when configured with -DBUILD_QUADMATH=ON, which requires a
 * compiler that exposes the __float128 type (GCC, Clang, Intel ICC/ICX
 * on Linux/macOS).  Microsoft Visual C++ does not provide __float128
 * and is unsupported here.
 *
 * Argument and phase conventions are identical to wigner.h: twice-j
 * integer arguments throughout, Condon--Shortley sign for Y_l^m, etc.
 *
 * The conversion from the exact prime-factorized result to __float128
 * combines the top three 64-bit words of the bigint with float128
 * arithmetic, giving 192 input bits for the 113-bit mantissa.  The
 * 79 unused input bits cover all but the most pathological tied-half-ulp
 * cases.  Results are correct to within 2 ulps in the worst case (the
 * verification harness asserts agreement within this tolerance against
 * the MPFR back-end at precision 113).
 */
#ifndef WIGNER_QUADMATH_H
#define WIGNER_QUADMATH_H

#include <quadmath.h>

#ifdef __cplusplus
extern "C" {
#endif

__float128 wigner3j_q(int tj1, int tj2, int tj3,
                      int tm1, int tm2, int tm3);

__float128 wigner6j_q(int tj1, int tj2, int tj3,
                      int tj4, int tj5, int tj6);

__float128 wigner9j_q(int tj11, int tj12, int tj13,
                      int tj21, int tj22, int tj23,
                      int tj31, int tj32, int tj33);

__float128 clebsch_gordan_q(int tj1, int tm1, int tj2, int tm2,
                            int tJ,  int tM);

__float128 racah_w_q(int tj1, int tj2, int tJ,
                     int tj3, int tj12, int tj23);

__float128 fano_x_q(int tj1, int tj2, int tj12,
                    int tj3, int tj4, int tj34,
                    int tj13, int tj24, int tJ);

__float128 gaunt_q     (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
__float128 gaunt_real_q(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);

#ifdef __cplusplus
}
#endif

#endif /* WIGNER_QUADMATH_H */
