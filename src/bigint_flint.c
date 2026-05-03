/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * FLINT-backed implementation of the bigint API.  Selected by
 * -DBUILD_FLINT=ON; otherwise the schoolbook implementation in
 * src/bigint.c is used.  Both backends share the public bigint.h API
 * and produce bit-identical floating-point output, verified in CI.
 *
 * The asymptotic motivation is sub-quadratic multiplication:  FLINT
 * dispatches large multiword integer products through the same
 * Karatsuba / Toom-Cook / Schönhage--Strassen algorithms that GMP
 * implements (FLINT's `fmpz_t' wraps an `mpz_t' for values that do
 * not fit in a single inline 62-bit slot), so this backend closes
 * most of the speed gap that the schoolbook path opens against
 * \texttt{WIGXJPF} at very large angular momenta.  At small or
 * moderate angular momenta the inline-int representation also avoids
 * the heap allocation that an `mpz_t' would incur.
 *
 * Floating-point conversions are routed through MPFR (always
 * available because FLINT depends on it), which guarantees correct
 * round-to-nearest-even rounding at every IEEE 754 binary precision.
 */
#include "bigint.h"

#ifndef WIGNERNJ_USE_FLINT
#error "bigint_flint.c compiled without WIGNERNJ_USE_FLINT defined"
#endif

#include <mpfr.h>
#include <math.h>
#include <float.h>

/* ── lifecycle ───────────────────────────────────────────────────────────── */

void bigint_init(bigint_t *a)            { fmpz_init(a->v); }
void bigint_free(bigint_t *a)            { fmpz_clear(a->v); }
void bigint_set_zero(bigint_t *a)        { fmpz_zero(a->v); }
int  bigint_is_zero(const bigint_t *a)   { return fmpz_is_zero(a->v); }

void bigint_copy(bigint_t *dst, const bigint_t *src)
{
    fmpz_set(dst->v, src->v);
}

void bigint_set_u64(bigint_t *a, uint64_t v)
{
    /* Build the value via mpz to avoid assuming 64-bit `ulong'. */
    mpz_t m;
    mpz_init(m);
    mpz_set_ui(m, (unsigned long)(v >> 32));
    mpz_mul_2exp(m, m, 32);
    mpz_add_ui(m, m, (unsigned long)(v & 0xFFFFFFFFUL));
    fmpz_set_mpz(a->v, m);
    mpz_clear(m);
}

void bigint_reserve(bigint_t *a, size_t cap)
{
    /* No-op:  fmpz_t auto-grows.  The schoolbook path uses this hint to
     * avoid realloc traffic in inner loops; FLINT manages its own
     * memory so the hint is irrelevant here. */
    (void)a;
    (void)cap;
}

/* ── comparison ──────────────────────────────────────────────────────────── */

int bigint_cmp(const bigint_t *a, const bigint_t *b)
{
    /* The schoolbook variant compares unsigned magnitudes.  In
     * libwignernj all bigints stored in wigner_exact_t are
     * non-negative, so a signed compare is equivalent.  Normalise
     * to {-1, 0, +1}. */
    int c = fmpz_cmp(a->v, b->v);
    return (c < 0) ? -1 : (c > 0) ? 1 : 0;
}

/* ── arithmetic ──────────────────────────────────────────────────────────── */

void bigint_add(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    fmpz_add(r->v, a->v, b->v);
}

int bigint_sub_signed(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    /* Schoolbook contract: r = ||a| - |b||, return sign of (|a| - |b|),
     * with +1 for the equal-magnitude case.  In libwignernj a and b
     * are always non-negative so this reduces to abs(a - b). */
    fmpz_sub(r->v, a->v, b->v);
    if (fmpz_sgn(r->v) < 0) {
        fmpz_neg(r->v, r->v);
        return -1;
    }
    return +1;
}

void bigint_mul(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    fmpz_mul(r->v, a->v, b->v);
}

void bigint_mul_u64(bigint_t *r, const bigint_t *a, uint64_t b)
{
    fmpz_mul_ui(r->v, a->v, (unsigned long)b);
}

void bigint_div_u64(bigint_t *r, const bigint_t *a, uint64_t b)
{
    /* Caller guarantees exactness:  a is a multiple of b. */
    fmpz_divexact_ui(r->v, a->v, (unsigned long)b);
}

void bigint_mul_prime_pow(bigint_t *a, uint64_t p, int k)
{
    if (k <= 0) return;
    if (p == 2) {
        fmpz_mul_2exp(a->v, a->v, (unsigned long)k);
        return;
    }
    /* Compute p^k once and multiply in.  fmpz_pow_ui is FLINT's
     * binary-exponentiation routine, so we don't need to roll our
     * own. */
    fmpz_t pk;
    fmpz_init(pk);
    fmpz_set_ui(pk, (unsigned long)p);
    fmpz_pow_ui(pk, pk, (unsigned long)k);
    fmpz_mul(a->v, a->v, pk);
    fmpz_clear(pk);
}

void bigint_div_prime_pow(bigint_t *a, uint64_t p, int k)
{
    if (k <= 0) return;
    if (p == 2) {
        /* a is non-negative and exactly divisible, so a right shift
         * is exact-floor-division. */
        fmpz_tdiv_q_2exp(a->v, a->v, (unsigned long)k);
        return;
    }
    fmpz_t pk;
    fmpz_init(pk);
    fmpz_set_ui(pk, (unsigned long)p);
    fmpz_pow_ui(pk, pk, (unsigned long)k);
    fmpz_divexact(a->v, a->v, pk);
    fmpz_clear(pk);
}

/* ── workspaces (no-op for FLINT) ────────────────────────────────────────── */

void bigint_ws_init(bigint_ws_t *ws)
{
    fmpz_init(ws->mul_temp.v);
    fmpz_init(ws->pp_pw.v);
    fmpz_init(ws->pp_base.v);
    fmpz_init(ws->kar_scratch.v);   /* unused on this backend */
}

void bigint_ws_free(bigint_ws_t *ws)
{
    fmpz_clear(ws->mul_temp.v);
    fmpz_clear(ws->pp_pw.v);
    fmpz_clear(ws->pp_base.v);
    fmpz_clear(ws->kar_scratch.v);
}

void bigint_ws_reserve(bigint_ws_t *ws, size_t max_words)
{
    (void)ws;
    (void)max_words;
}

void bigint_mul_ws(bigint_t *r, const bigint_t *a, const bigint_t *b,
                    bigint_ws_t *ws)
{
    (void)ws;
    fmpz_mul(r->v, a->v, b->v);
}

void bigint_mul_prime_pow_ws(bigint_t *a, uint64_t p, int k, bigint_ws_t *ws)
{
    (void)ws;
    bigint_mul_prime_pow(a, p, k);
}

/* ── factorial-size hint ─────────────────────────────────────────────────── */

size_t bigint_words_for_factorial(int N)
{
    /* Used only as a pre-allocation hint by the schoolbook path.
     * fmpz_t auto-grows, so the value here is irrelevant; return a
     * small positive constant. */
    (void)N;
    return 4;
}

/* ── bit length ──────────────────────────────────────────────────────────── */

int bigint_bit_length(const bigint_t *a)
{
    if (fmpz_is_zero(a->v)) return 0;
    return (int)fmpz_bits(a->v);
}

/* ── floating-point conversions, routed through MPFR for correct rounding ── */

static void bigint_to_mpfr_internal(mpfr_t v, const bigint_t *a, mpfr_rnd_t rnd)
{
    mpz_t m;
    mpz_init(m);
    fmpz_get_mpz(m, a->v);
    mpfr_set_z(v, m, rnd);
    mpz_clear(m);
}

float bigint_to_float(const bigint_t *a)
{
    if (fmpz_is_zero(a->v)) return 0.0f;
    mpfr_t v;
    mpfr_init2(v, FLT_MANT_DIG);
    bigint_to_mpfr_internal(v, a, MPFR_RNDN);
    float r = mpfr_get_flt(v, MPFR_RNDN);
    mpfr_clear(v);
    return r;
}

double bigint_to_double(const bigint_t *a)
{
    if (fmpz_is_zero(a->v)) return 0.0;
    mpfr_t v;
    mpfr_init2(v, DBL_MANT_DIG);
    bigint_to_mpfr_internal(v, a, MPFR_RNDN);
    double r = mpfr_get_d(v, MPFR_RNDN);
    mpfr_clear(v);
    return r;
}

long double bigint_to_long_double(const bigint_t *a)
{
    if (fmpz_is_zero(a->v)) return 0.0L;
    mpfr_t v;
    mpfr_init2(v, LDBL_MANT_DIG);
    bigint_to_mpfr_internal(v, a, MPFR_RNDN);
    long double r = mpfr_get_ld(v, MPFR_RNDN);
    mpfr_clear(v);
    return r;
}

/* The frexp helpers must work for bigints far beyond the dynamic range
 * of the target floating-point type (factorials at j ~ 100 already
 * overflow IEEE 754 binary64).  Going through bigint_to_<T> first would
 * collapse these to infinity; instead we extract a normalised
 * mantissa-exponent pair via MPFR (which has unbounded exponent
 * range) and let wigner_exact_to_<T> recombine the exponents in
 * integer arithmetic. */
static double frexp_via_mpfr(const bigint_t *a, mpfr_prec_t prec,
                             int *out_exp)
{
    mpz_t m;
    mpz_init(m);
    fmpz_get_mpz(m, a->v);
    mpfr_t v;
    mpfr_init2(v, prec);
    mpfr_set_z(v, m, MPFR_RNDN);
    long exp = 0;
    double mant = mpfr_get_d_2exp(&exp, v, MPFR_RNDN);
    *out_exp = (int)exp;
    mpfr_clear(v);
    mpz_clear(m);
    return mant;
}

float bigint_frexpf(const bigint_t *a, int *out_exp)
{
    if (fmpz_is_zero(a->v)) { *out_exp = 0; return 0.0f; }
    return (float)frexp_via_mpfr(a, FLT_MANT_DIG, out_exp);
}

double bigint_frexp(const bigint_t *a, int *out_exp)
{
    if (fmpz_is_zero(a->v)) { *out_exp = 0; return 0.0; }
    return frexp_via_mpfr(a, DBL_MANT_DIG, out_exp);
}

long double bigint_frexpl(const bigint_t *a, int *out_exp)
{
    if (fmpz_is_zero(a->v)) { *out_exp = 0; return 0.0L; }
    /* mpfr_get_ld_2exp expects MPFR with at least LDBL_MANT_DIG bits;
     * fall back to the double path's mantissa, scaled to long double,
     * since wigner_exact_to_long_double recombines the exponents
     * separately and only needs a normalised mantissa. */
    mpz_t m;
    mpz_init(m);
    fmpz_get_mpz(m, a->v);
    mpfr_t v;
    mpfr_init2(v, LDBL_MANT_DIG);
    mpfr_set_z(v, m, MPFR_RNDN);
    long exp = 0;
    long double mant = mpfr_get_ld_2exp(&exp, v, MPFR_RNDN);
    *out_exp = (int)exp;
    mpfr_clear(v);
    mpz_clear(m);
    return mant;
}

/* ── binary128 conversion (when the compiler provides __float128) ────────── */

#ifdef WIGNER_HAVE_QUADMATH
/* mpfr_get_float128 is available in MPFR ≥ 4.1.0, but only when MPFR
 * was configured with --enable-float128 (Ubuntu's libmpfr-dev is not).
 * Avoid the dependency entirely by extracting the rounded 113-bit
 * mantissa with mpfr_get_z_2exp and Horner-assembling __float128
 * from its GMP limbs:  the mantissa fits in 113 bits, well within
 * the 113-bit __float128 significand, so the assembly is exact and
 * the only rounding step is the mpfr_set_z above. */
static __float128 mpfr_to_float128(mpfr_t v, int *out_exp2)
{
    mpz_t mant;
    mpz_init(mant);
    mpfr_exp_t e2 = mpfr_get_z_2exp(mant, v);
    int sign = mpz_sgn(mant);
    if (sign < 0) mpz_neg(mant, mant);
    const __float128 limb_radix = ldexpq(1.0Q, GMP_LIMB_BITS);
    __float128 r = 0.0Q;
    for (size_t i = mpz_size(mant); i-- > 0;) {
        r = r * limb_radix + (__float128)mpz_getlimbn(mant, i);
    }
    mpz_clear(mant);
    if (sign < 0) r = -r;
    *out_exp2 = (int)e2;
    return r;
}

__float128 bigint_to_float128(const bigint_t *a)
{
    if (fmpz_is_zero(a->v)) return 0.0Q;
    mpfr_t v;
    mpfr_init2(v, 113);            /* binary128 mantissa bits */
    bigint_to_mpfr_internal(v, a, MPFR_RNDN);
    int e2;
    __float128 r = mpfr_to_float128(v, &e2);
    mpfr_clear(v);
    return ldexpq(r, e2);
}

__float128 bigint_frexp_q(const bigint_t *a, int *out_exp)
{
    if (fmpz_is_zero(a->v)) { *out_exp = 0; return 0.0Q; }
    mpfr_t v;
    mpfr_init2(v, 113);
    bigint_to_mpfr_internal(v, a, MPFR_RNDN);
    int e2;
    __float128 r = mpfr_to_float128(v, &e2);
    mpfr_clear(v);
    /* r holds an integer in [2^112, 2^113) (or a single limb when the
     * value is small).  Renormalise to [0.5, 1). */
    int e3;
    __float128 mant = frexpq(r, &e3);
    *out_exp = e2 + e3;
    return mant;
}
#endif

/* ── MPFR conversion (when the user requests the MPFR back-end) ──────────── */

#ifdef WIGNER_HAVE_MPFR
void bigint_to_mpfr(mpfr_t rop, const bigint_t *a, mpfr_rnd_t rnd)
{
    bigint_to_mpfr_internal(rop, a, rnd);
}
#endif
