/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * MPFR arbitrary-precision interface for libwignernj.
 *
 * Include this header (in addition to wignernj.h) when MPFR output is needed.
 * The library must have been built with -DBUILD_MPFR=ON.
 *
 * Convention: rop is the first argument; set its precision before calling.
 * rnd is the MPFR rounding mode applied to each elementary operation.
 *
 * Phase conventions: identical to the floating-point API; see the header
 * comment of wignernj.h.  In short, 3j/6j/9j/CG/Racah W are convention-free
 * SU(2) algebraic objects (CG signs use Condon–Shortley); Gaunt, real-Gaunt
 * and the real ↔ complex Y_lm basis-overlap matrix assume the
 * Condon–Shortley phase for Y_l^m.
 *
 * Example:
 *   mpfr_t v;
 *   mpfr_init2(v, 256);
 *   wigner3j_mpfr(v, 2, 2, 0,  0, 0, 0,  MPFR_RNDN);
 *   mpfr_clear(v);
 */
#ifndef WIGNERNJ_MPFR_H
#define WIGNERNJ_MPFR_H

#include <mpfr.h>
#include "wignernj.h"

/* Wigner 3j */
void wigner3j_mpfr(mpfr_t rop,
                   int tj1, int tj2, int tj3,
                   int tm1, int tm2, int tm3,
                   mpfr_rnd_t rnd);

/* Wigner 6j */
void wigner6j_mpfr(mpfr_t rop,
                   int tj1, int tj2, int tj3,
                   int tj4, int tj5, int tj6,
                   mpfr_rnd_t rnd);

/* Wigner 9j (row-major order) */
void wigner9j_mpfr(mpfr_t rop,
                   int tj11, int tj12, int tj13,
                   int tj21, int tj22, int tj23,
                   int tj31, int tj32, int tj33,
                   mpfr_rnd_t rnd);

/* Clebsch-Gordan  <j1 m1; j2 m2 | J M> */
void clebsch_gordan_mpfr(mpfr_t rop,
                         int tj1, int tm1,
                         int tj2, int tm2,
                         int tJ,  int tM,
                         mpfr_rnd_t rnd);

/* Racah W-coefficient  W(j1 j2 J j3; j12 j23) */
void racah_w_mpfr(mpfr_t rop,
                  int tj1, int tj2, int tJ, int tj3,
                  int tj12, int tj23,
                  mpfr_rnd_t rnd);

/* Fano X-coefficient  X(j1 j2 j12; j3 j4 j34; j13 j24 J) */
void fano_x_mpfr(mpfr_t rop,
                 int tj1, int tj2, int tj12,
                 int tj3, int tj4, int tj34,
                 int tj13, int tj24, int tJ,
                 mpfr_rnd_t rnd);

/* Gaunt coefficient  ∫ Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ */
void gaunt_mpfr(mpfr_t rop,
                int tl1, int tm1,
                int tl2, int tm2,
                int tl3, int tm3,
                mpfr_rnd_t rnd);

/* Real-spherical-harmonic Gaunt  ∫ S_{l1,m1} S_{l2,m2} S_{l3,m3} dΩ */
void gaunt_real_mpfr(mpfr_t rop,
                     int tl1, int tm1,
                     int tl2, int tm2,
                     int tl3, int tm3,
                     mpfr_rnd_t rnd);

/* Real <-> complex spherical-harmonic basis transformation at MPFR
 * precision.  Fills two parallel column-major mpfr_t arrays of length
 * (2l+1)^2 each (real and imaginary parts).  Element C[m_r, m_c]
 * lives at flat index (m_c+l)*(2l+1) + (m_r+l) in both C_re and C_im.
 * Each mpfr_t in both arrays must be mpfr_init2'd to the desired
 * precision before the call; the function only writes values.
 * Semantics otherwise identical to wignernj_real_ylm_in_complex_ylm in
 * wignernj.h (same convention, same column-major layout). */
void wignernj_real_ylm_in_complex_ylm_mpfr(int l,
                                     mpfr_t *C_re,
                                     mpfr_t *C_im,
                                     mpfr_rnd_t rnd);

#endif /* WIGNERNJ_MPFR_H */
