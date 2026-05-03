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
 *
 * Phase conventions
 * -----------------
 * The Wigner 3j, 6j, and 9j symbols, the Clebsch–Gordan coefficient, and
 * the Racah W coefficient are pure algebraic SU(2) objects.  Their values
 * are fixed entirely by the Racah/Wigner combinatorial formulas; no
 * spherical-harmonic phase convention enters anywhere.  The Clebsch–Gordan
 * coefficient is defined with the Condon–Shortley sign convention used by
 * Edmonds (1957) and Varshalovich et al. (1988):
 *   <j1 m1; j2 m2 | J M> = (-1)^(j1-j2+M) sqrt(2J+1) * 3j(j1,j2,J;m1,m2,-M),
 * which makes every Clebsch–Gordan coefficient real.
 *
 * The phase convention for Y_l^m enters only the Gaunt routines below
 * (gaunt and gaunt_real), since both are defined as integrals of three
 * spherical harmonics.  We adopt the Condon–Shortley phase for Y_l^m,
 *   Y_l^m(θ,φ) = (-1)^m sqrt[(2l+1)/(4π) · (l-m)!/(l+m)!] P_l^m(cos θ) e^{i m φ}
 *                  (m ≥ 0),
 *   Y_l^{-m}   = (-1)^m (Y_l^m)*,
 * which is the convention of Edmonds, Varshalovich, and the bulk of the
 * atomic-, molecular-, and nuclear-physics literature.  The real spherical
 * harmonics used by gaunt_real follow the Wikipedia/Condon–Shortley
 * construction; see the comment block above gaunt_real() below.
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

/* ── Fano X-coefficient ──────────────────────────────────────────────────── */
/* X(j1 j2 j12; j3 j4 j34; j13 j24 J)                                        */
/*   = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)]                                */
/*     * { j1   j2   j12 }                                                   */
/*       { j3   j4   j34 }                                                   */
/*       { j13  j24  J   }                                                   */
/* This is a normalisation variant of the 9j symbol (Fano 1953;             */
/* Edmonds 1957 §6.4) used in the analysis of polarisation correlations.   */

float       fano_x_f(int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);
double      fano_x  (int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);
long double fano_x_l(int tj1, int tj2, int tj12,
                     int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);

/* ── Gaunt coefficient ───────────────────────────────────────────────────── */
/* G(l1,m1,l2,m2,l3,m3)                                                      */
/*   = integral Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ                     */
/*   = sqrt[(2l1+1)(2l2+1)(2l3+1)/(4π)]                                     */
/*     * (l1 l2 l3; 0 0 0) * (l1 l2 l3; m1 m2 m3)                          */
/* The closed-form expression in terms of two 3j symbols assumes the         */
/* Condon–Shortley phase for Y_l^m (see file header above).                  */
/* l arguments must be non-negative integers; tl = 2*l (always even).       */

float       gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
double      gaunt  (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);

/* ── Real-spherical-harmonic Gaunt coefficient ──────────────────────────── */
/* G^R(l1,m1,l2,m2,l3,m3) = integral S_{l1,m1} S_{l2,m2} S_{l3,m3} dΩ        */
/*                                                                           */
/* with the Condon–Shortley / Wikipedia convention for the real spherical    */
/* harmonics S_{l,m}:                                                        */
/*   S_{l, 0}    = Y_l^0                                                     */
/*   S_{l,  m>0} = (1/√2)(Y_l^{-m} + (-1)^m Y_l^m)                           */
/*   S_{l,  m<0} = (i/√2)(Y_l^{ m} - (-1)^|m| Y_l^{-m})                      */
/*                                                                           */
/* ℓ and |m| arguments must be non-negative integers; tl = 2*l, tm = 2*m,    */
/* both always even.  The implementation runs at the cost of a single        */
/* complex-Gaunt evaluation per call: see src/gaunt.c for the algorithm.     */

float       gaunt_real_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
double      gaunt_real  (int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
long double gaunt_real_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);

/* ── Optional warmup and scratch introspection ──────────────────────────────
 *
 * libwignernj caches a per-thread scratch buffer that is lazy-allocated on
 * the first call from a given thread and reused for every subsequent call.
 * After the first call, no further allocations happen unless the angular
 * momenta exceed what the cached buffer was last sized for; in that case
 * the buffer grows in place.  The lazy growth is invisible to the caller.
 *
 * The cache requires thread-local storage and is therefore only active on
 * toolchains that provide one of `__thread` (GCC/Clang/Intel),
 * `__declspec(thread)` (MSVC), or C11 `_Thread_local`.  On a toolchain
 * that exposes none of these, libwignernj falls back to allocating a
 * fresh scratch on every public-API call and freeing it on return.  That
 * fallback is thread-safe by construction (each call owns its own
 * scratch) but pays the historical per-call allocation cost.
 *
 * `wigner_warmup` is an optional convenience that pre-allocates the
 * calling thread's cached scratch for the absolute default-build maximum
 * (j1+j2+j3 ≤ 19999 for 3j/6j/CG/Racah W/Gaunt; equal-j ≤ 4999 for 9j and
 * Fano X).  After it returns, every subsequent symbol evaluation in this
 * thread is guaranteed allocation-free, provided the cache is available.
 * Useful in benchmarks where the first-call lazy-init cost would
 * otherwise be observed, in hot loops where allocation overhead matters,
 * or in applications that prefer predictable up-front memory usage to
 * lazy growth.
 *
 * Memory cost: ~6 MB per thread at the absolute ceiling.  The function is
 * thread-safe and only touches the calling thread's scratch.  Calling it
 * is purely an optimisation; correctness is unchanged whether or not it
 * is called, and it is harmless to call it more than once.
 *
 * On a toolchain without thread-local storage, `wigner_warmup` is a
 * no-op: there is no persistent scratch to pre-grow, so the call returns
 * immediately and subsequent symbol evaluations continue to allocate
 * per-call.  Use `wigner_thread_local_scratch_available()` to detect at
 * runtime whether the cache is in effect; it returns 1 when the per-
 * thread cache is active and 0 when each call allocates fresh. */
void wigner_warmup(void);
int  wigner_thread_local_scratch_available(void);

#ifdef __cplusplus
}
#endif

#endif /* WIGNER_H */
