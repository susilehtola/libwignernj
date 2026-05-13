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
 * Argument and phase conventions are identical to wignernj.h: twice-j
 * integer arguments for every coupling-coefficient routine,
 * Condon--Shortley sign for Y_l^m, plain integer `l` for the
 * wignernj_real_ylm_in_complex_ylm family.
 *
 * The conversion from the exact prime-factorized result to __float128
 * combines the top three 64-bit words of the bigint with float128
 * arithmetic, giving 192 input bits for the 113-bit mantissa.  The
 * 79 unused input bits cover all but the most pathological tied-half-ulp
 * cases.  Results are correct to within 2 ulps in the worst case (the
 * verification harness asserts agreement within this tolerance against
 * the MPFR back-end at precision 113).
 */
#ifndef WIGNERNJ_QUADMATH_H
#define WIGNERNJ_QUADMATH_H

#include <quadmath.h>
#include "wignernj.h"

/* libquadmath provides `__complex128` as `__float128 _Complex` on
 * gcc/clang/Intel (where the `_Complex` keyword is recognised); no
 * MSVC complication here because MSVC does not implement __float128
 * at all (and BUILD_QUADMATH is gated off on MSVC).  The C++ branch
 * exposes a layout-compatible struct so the C-side declaration
 * parses in C++ TUs; C++ consumers should go through the
 * `wignernj::real_ylm_in_complex_ylm<__float128>(...)` template overload
 * in wignernj.hpp. */
#ifndef __cplusplus
typedef __complex128 wignernj_cfloat128_t;
#else
typedef struct { __float128 _pair[2]; } wignernj_cfloat128_t;
#endif

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

/* Real <-> complex spherical-harmonic basis transformation at __float128
 * precision.  Layout and semantics identical to the
 * wignernj_real_ylm_in_complex_ylm family in wignernj.h: column-major, leading
 * dimension (2l+1), interleaved (re, im) inside each complex slot. */
void wignernj_real_ylm_in_complex_ylm_q(int l, wignernj_cfloat128_t *C_out);

#ifdef __cplusplus
}
#endif

#endif /* WIGNERNJ_QUADMATH_H */
