/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Public C API for libwignernj.
 *
 * All angular momentum arguments are passed as TWICE their value so that
 * half-integer values are represented as odd integers:
 *   j = 3/2  →  tj = 3
 *   m = -1/2 →  tm = -1
 *
 * Functions returning 0 / 0.0f / 0.0L for symbols that vanish by selection
 * rules are not errors; only programmer-visible violations (wrong parity,
 * |m| > j) are treated as returning zero without any diagnostics.
 */
#ifndef WIGNER_H
#define WIGNER_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── Wigner 3j symbol ────────────────────────────────────────────────────── */
/* ( j1  j2  j3 )                                                            */
/* ( m1  m2  m3 )                                                            */

float       wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
double      wigner3j  (int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);

/* ── Wigner 6j symbol ────────────────────────────────────────────────────── */
/* { j1 j2 j3 }                                                              */
/* { j4 j5 j6 }                                                              */

float       wigner6j_f(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
double      wigner6j  (int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
long double wigner6j_l(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);

/* ── Wigner 9j symbol ────────────────────────────────────────────────────── */
/* { j11 j12 j13 }  (row-major order)                                        */
/* { j21 j22 j23 }                                                           */
/* { j31 j32 j33 }                                                           */

float       wigner9j_f(int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);
double      wigner9j  (int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);
long double wigner9j_l(int tj11, int tj12, int tj13,
                       int tj21, int tj22, int tj23,
                       int tj31, int tj32, int tj33);

/* ── Clebsch-Gordan coefficient ──────────────────────────────────────────── */
/* <j1 m1; j2 m2 | J M>  =  (-1)^(j1-j2+M) sqrt(2J+1) * (j1 j2  J )       */
/*                                                          (m1 m2 -M)       */

float       clebsch_gordan_f(int tj1, int tm1, int tj2, int tm2,
                              int tJ, int tM);
double      clebsch_gordan  (int tj1, int tm1, int tj2, int tm2,
                              int tJ, int tM);
long double clebsch_gordan_l(int tj1, int tm1, int tj2, int tm2,
                              int tJ, int tM);

/* ── Racah W-coefficient ─────────────────────────────────────────────────── */
/* W(j1 j2 J j3; j12 j23)  =  (-1)^(j1+j2+J+j3) * {j1  j2  j12}           */
/*                                                    {j3  J   j23}           */

float       racah_w_f(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
double      racah_w  (int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
long double racah_w_l(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);

/* ── Gaunt coefficient ───────────────────────────────────────────────────── */
/* G(l1,m1,l2,m2,l3,m3)                                                      */
/*   = integral Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ                     */
/*   = sqrt[(2l1+1)(2l2+1)(2l3+1)/(4π)]                                     */
/*     * (l1 l2 l3; 0 0 0) * (l1 l2 l3; m1 m2 m3)                          */
/* l arguments must be non-negative integers; tl = 2*l (always even).       */

float       gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
double      gaunt  (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);

#ifdef __cplusplus
}
#endif

#endif /* WIGNER_H */
