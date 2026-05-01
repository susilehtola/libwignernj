/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */

/* __uint128_t is required; available on GCC, Clang, and ICC.
 * MSVC does not support this type — see issue tracker for status. */
#if defined(_MSC_VER)
#  error "bigint.c requires __uint128_t (GCC/Clang/ICC). MSVC is not yet supported."
#endif

#include "bigint.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ── internal helpers ────────────────────────────────────────────────────── */

static void bigint_ensure(bigint_t *a, size_t need)
{
    if (need > a->cap) {
        size_t nc = a->cap ? a->cap * 2 : 4;
        while (nc < need) nc *= 2;
        a->words = (uint64_t *)xrealloc(a->words, nc * sizeof(uint64_t));
        memset(a->words + a->cap, 0, (nc - a->cap) * sizeof(uint64_t));
        a->cap = nc;
    }
}

static void bigint_trim(bigint_t *a)
{
    while (a->size > 0 && a->words[a->size - 1] == 0)
        a->size--;
}

/* ── lifecycle ───────────────────────────────────────────────────────────── */

void bigint_init(bigint_t *a)     { a->words = NULL; a->size = a->cap = 0; }
void bigint_free(bigint_t *a)     { free(a->words); a->words = NULL; a->size = a->cap = 0; }
void bigint_set_zero(bigint_t *a) { a->size = 0; }
int  bigint_is_zero(const bigint_t *a) { return a->size == 0; }

void bigint_copy(bigint_t *dst, const bigint_t *src)
{
    if (src->size == 0) { dst->size = 0; return; }
    bigint_ensure(dst, src->size);
    memcpy(dst->words, src->words, src->size * sizeof(uint64_t));
    dst->size = src->size;
}

void bigint_set_u64(bigint_t *a, uint64_t v)
{
    bigint_ensure(a, 1);
    a->words[0] = v;
    a->size = (v != 0) ? 1 : 0;
}

/* ── comparison ──────────────────────────────────────────────────────────── */

int bigint_cmp(const bigint_t *a, const bigint_t *b)
{
    size_t i;
    if (a->size != b->size) return (a->size > b->size) ? 1 : -1;
    for (i = a->size; i-- > 0;)
        if (a->words[i] != b->words[i])
            return (a->words[i] > b->words[i]) ? 1 : -1;
    return 0;
}

/* ── addition ────────────────────────────────────────────────────────────── */

void bigint_add(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    const bigint_t *big = (a->size >= b->size) ? a : b;
    const bigint_t *sml = (a->size >= b->size) ? b : a;
    size_t n = big->size, i;
    uint64_t carry = 0;
    bigint_t tmp;
    bigint_init(&tmp);
    bigint_ensure(&tmp, n + 1);
    tmp.size = n + 1;
    for (i = 0; i < sml->size; i++) {
        __uint128_t s = (__uint128_t)big->words[i] + sml->words[i] + carry;
        tmp.words[i] = (uint64_t)s;
        carry        = (uint64_t)(s >> 64);
    }
    for (; i < n; i++) {
        __uint128_t s = (__uint128_t)big->words[i] + carry;
        tmp.words[i] = (uint64_t)s;
        carry        = (uint64_t)(s >> 64);
    }
    tmp.words[n] = carry;
    bigint_trim(&tmp);
    bigint_copy(r, &tmp);
    bigint_free(&tmp);
}

/* ── subtraction (unsigned, a >= b assumed) ─────────────────────────────── */

static void bigint_sub_uu(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    uint64_t borrow = 0;
    size_t i;
    bigint_t tmp;
    bigint_init(&tmp);
    bigint_ensure(&tmp, a->size);
    tmp.size = a->size;
    for (i = 0; i < b->size; i++) {
        __uint128_t d = (__uint128_t)a->words[i] - b->words[i] - borrow;
        tmp.words[i] = (uint64_t)d;
        borrow       = (uint64_t)(-(int64_t)(d >> 64)); /* 1 if underflow */
    }
    for (; i < a->size; i++) {
        __uint128_t d = (__uint128_t)a->words[i] - borrow;
        tmp.words[i] = (uint64_t)d;
        borrow       = (uint64_t)(-(int64_t)(d >> 64));
    }
    bigint_trim(&tmp);
    bigint_copy(r, &tmp);
    bigint_free(&tmp);
}

int bigint_sub_signed(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    int c = bigint_cmp(a, b);
    if (c == 0)  { bigint_set_zero(r); return 1; }
    if (c > 0)   { bigint_sub_uu(r, a, b); return  1; }
    bigint_sub_uu(r, b, a); return -1;
}

/* ── multiplication ──────────────────────────────────────────────────────── */

void bigint_mul_u64(bigint_t *r, const bigint_t *a, uint64_t b)
{
    size_t i;
    uint64_t carry = 0;
    bigint_t tmp;
    if (b == 0 || a->size == 0) { bigint_set_zero(r); return; }
    bigint_init(&tmp);
    bigint_ensure(&tmp, a->size + 1);
    tmp.size = a->size + 1;
    for (i = 0; i < a->size; i++) {
        __uint128_t p = (__uint128_t)a->words[i] * b + carry;
        tmp.words[i] = (uint64_t)p;
        carry        = (uint64_t)(p >> 64);
    }
    tmp.words[a->size] = carry;
    bigint_trim(&tmp);
    bigint_copy(r, &tmp);
    bigint_free(&tmp);
}

void bigint_mul(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    size_t i, j;
    bigint_t res;
    if (a->size == 0 || b->size == 0) { bigint_set_zero(r); return; }
    bigint_init(&res);
    bigint_ensure(&res, a->size + b->size);
    memset(res.words, 0, (a->size + b->size) * sizeof(uint64_t));
    res.size = a->size + b->size;
    for (i = 0; i < a->size; i++) {
        uint64_t carry = 0;
        for (j = 0; j < b->size; j++) {
            __uint128_t p = (__uint128_t)a->words[i] * b->words[j]
                          + res.words[i+j] + carry;
            res.words[i+j] = (uint64_t)p;
            carry          = (uint64_t)(p >> 64);
        }
        res.words[i + b->size] += carry;
    }
    bigint_trim(&res);
    bigint_copy(r, &res);
    bigint_free(&res);
}

/* ── division by scalar ──────────────────────────────────────────────────── */

void bigint_div_u64(bigint_t *r, const bigint_t *a, uint64_t b)
{
    size_t i;
    __uint128_t rem = 0;
    bigint_t tmp;
    bigint_init(&tmp);
    bigint_ensure(&tmp, a->size);
    tmp.size = a->size;
    for (i = a->size; i-- > 0;) {
        __uint128_t cur = (rem << 64) | a->words[i];
        tmp.words[i] = (uint64_t)(cur / b);
        rem          = cur % b;
    }
    bigint_trim(&tmp);
    bigint_copy(r, &tmp);
    bigint_free(&tmp);
}

void bigint_mul_prime_pow(bigint_t *a, uint64_t p, int k)
{
    bigint_t pw, base;
    int i;
    if (k <= 0) return;
    if (k <= 8) {
        for (i = 0; i < k; i++) bigint_mul_u64(a, a, p);
        return;
    }
    /* Binary exponentiation: compute p^k, then a *= p^k. */
    bigint_init(&pw);
    bigint_init(&base);
    bigint_set_u64(&pw, 1);
    bigint_set_u64(&base, p);
    while (k > 0) {
        if (k & 1) bigint_mul(&pw, &pw, &base);
        k >>= 1;
        if (k > 0) bigint_mul(&base, &base, &base);
    }
    bigint_mul(a, a, &pw);
    bigint_free(&pw);
    bigint_free(&base);
}

void bigint_div_prime_pow(bigint_t *a, uint64_t p, int k)
{
    int i;
    for (i = 0; i < k; i++) bigint_div_u64(a, a, p);
}

/* ── bit length ──────────────────────────────────────────────────────────── */

int bigint_bit_length(const bigint_t *a)
{
    uint64_t top;
    int bits;
    if (a->size == 0) return 0;
    top  = a->words[a->size - 1];
    bits = (int)(a->size - 1) * 64;
    while (top) { bits++; top >>= 1; }
    return bits;
}

/* ── IEEE 754 float conversion ───────────────────────────────────────────── */

/*
 * Extract the top `mant` bits (mant <= 64) of `a` as a uint64_t, applying
 * round-to-nearest-even.  Sets *out_exp so that result == mantissa * 2^out_exp.
 */
static uint64_t extract_mantissa(const bigint_t *a, int mant, int *out_exp)
{
    int n, shift, ws, bs, tsticky_word, tsticky_off, i;
    uint64_t m, lo, hi, mask;
    int round_bit, sticky;

    n = bigint_bit_length(a);
    if (n == 0) { *out_exp = 0; return 0; }

    shift    = n - mant;   /* bits to discard; negative means number is short */
    *out_exp = shift;

    if (shift <= 0) {
        /* Number has fewer bits than mantissa; left-shift to fill mant slots.
         * -shift < 64 is guaranteed since mant <= 64 and n >= 1. */
        m = a->words[0];
        if (shift < 0) m <<= -shift;
        return m;
    }

    /* Extract bits [shift .. shift+mant-1] from a */
    ws = shift / 64;
    bs = shift % 64;
    if (bs == 0) {
        lo = (ws < (int)a->size) ? a->words[ws] : 0;
        m  = lo;
    } else {
        lo = (ws     < (int)a->size) ? a->words[ws]     : 0;
        hi = (ws + 1 < (int)a->size) ? a->words[ws + 1] : 0;
        m  = (lo >> bs) | (hi << (64 - bs));
    }
    if (mant < 64)
        m &= ((uint64_t)1 << mant) - 1;

    /* Round bit: bit at position shift-1 */
    {
        int rw = (shift - 1) / 64, rb = (shift - 1) % 64;
        round_bit = (rw < (int)a->size) ? (int)((a->words[rw] >> rb) & 1) : 0;
    }

    /* Sticky bit: OR of bits 0 .. shift-2 */
    sticky = 0;
    if (shift >= 2) {
        tsticky_word = (shift - 2) / 64;
        tsticky_off  = (shift - 2) % 64;
        /* Bits 0..tsticky_off in word tsticky_word */
        mask = (tsticky_off == 63) ? UINT64_MAX
                                   : (((uint64_t)1 << (tsticky_off + 1)) - 1);
        if ((int)a->size > tsticky_word)
            sticky = (a->words[tsticky_word] & mask) != 0;
        /* All lower words */
        for (i = 0; i < tsticky_word && !sticky; i++)
            sticky = (a->words[i] != 0);
    }

    /* Round-to-nearest-even */
    if (round_bit && (sticky || (m & 1))) {
        m++;
        if (mant < 64 && (m >> mant)) { m >>= 1; (*out_exp)++; }
        else if (mant == 64 && m == 0) { m = (uint64_t)1 << 63; (*out_exp)++; }
    }
    return m;
}

#if LDBL_MANT_DIG > 64
/*
 * Extract the 64 bits of `a` immediately below the top 64 bits, normalized
 * so that bit 63 of the result corresponds to bit (n-65) of `a`.
 * The positional value of the result is:  return_value * 2^(n-128).
 * Requires n > 64 and no rounding overflow in the caller's extract_mantissa.
 */
static uint64_t lo64_from_words(const bigint_t *a, int n)
{
    int lo_top = n - 65;   /* 0-indexed bit position of the top of the second group */
    if (lo_top < 64) {
        uint64_t w0 = (a->size > 0) ? a->words[0] : 0;
        return (lo_top < 63) ? (w0 << (63 - lo_top)) : w0;
    } else {
        int word_hi = lo_top / 64;
        int off_hi  = lo_top % 64;
        uint64_t wh = (word_hi   < (int)a->size) ? a->words[word_hi]     : 0;
        uint64_t wl = (word_hi-1 < (int)a->size) ? a->words[word_hi - 1] : 0;
        return (off_hi == 63) ? wh
                              : ((wh << (63 - off_hi)) | (wl >> (off_hi + 1)));
    }
}
#endif

float bigint_to_float(const bigint_t *a)
{
    int exp2; uint64_t m;
    if (a->size == 0) return 0.0f;
    m = extract_mantissa(a, FLT_MANT_DIG, &exp2);
    return ldexpf((float)m, exp2);
}

float bigint_frexpf(const bigint_t *a, int *out_exp)
{
    int exp2; uint64_t m;
    if (a->size == 0) { *out_exp = 0; return 0.0f; }
    m = extract_mantissa(a, FLT_MANT_DIG, &exp2);
    *out_exp = exp2 + FLT_MANT_DIG;
    return ldexpf((float)m, -FLT_MANT_DIG);
}

double bigint_frexp(const bigint_t *a, int *out_exp)
{
    int exp2; uint64_t m;
    if (a->size == 0) { *out_exp = 0; return 0.0; }
    m = extract_mantissa(a, DBL_MANT_DIG, &exp2);
    *out_exp = exp2 + DBL_MANT_DIG;
    return ldexp((double)m, -DBL_MANT_DIG);
}

long double bigint_frexpl(const bigint_t *a, int *out_exp)
{
    int mant = (LDBL_MANT_DIG <= 64) ? LDBL_MANT_DIG : 64;
    int exp2; uint64_t m;
    if (a->size == 0) { *out_exp = 0; return 0.0L; }
    m = extract_mantissa(a, mant, &exp2);
    *out_exp = exp2 + mant;
#if LDBL_MANT_DIG > 64
    {
        /* For quad-precision long double: add a second 64-bit chunk below the top 64.
         * Contribution to the normalized mantissa is always m_lo * 2^-128
         * (since *out_exp = exp2+64 = n, so dividing by 2^n shifts by -n,
         *  and the lower bits sit at 2^(n-128), giving 2^(n-128)/2^n = 2^-128). */
        int n = bigint_bit_length(a);
        long double result = ldexpl((long double)m, -mant);
        if (n > 64 && exp2 == n - 64) {
            result += ldexpl((long double)lo64_from_words(a, n), -128);
        }
        return result;
    }
#else
    return ldexpl((long double)m, -mant);
#endif
}

double bigint_to_double(const bigint_t *a)
{
    int exp2; uint64_t m;
    if (a->size == 0) return 0.0;
    m = extract_mantissa(a, DBL_MANT_DIG, &exp2);
    return ldexp((double)m, exp2);
}

long double bigint_to_long_double(const bigint_t *a)
{
    int mant = (LDBL_MANT_DIG <= 64) ? LDBL_MANT_DIG : 64;
    int exp2; uint64_t m;
    if (a->size == 0) return 0.0L;
    m = extract_mantissa(a, mant, &exp2);
#if LDBL_MANT_DIG > 64
    {
        /* For quad-precision long double: add a second 64-bit chunk.
         * m's positional value is m * 2^exp2; the lower chunk adds m_lo * 2^(n-128). */
        int n = bigint_bit_length(a);
        long double hi = ldexpl((long double)m, exp2);
        if (n > 64 && exp2 == n - 64) {
            hi += ldexpl((long double)lo64_from_words(a, n), n - 128);
        }
        return hi;
    }
#else
    return ldexpl((long double)m, exp2);
#endif
}

/* ── MPFR conversion ─────────────────────────────────────────────────────── */

#ifdef WIGNER_HAVE_MPFR
void bigint_to_mpfr(mpfr_t rop, const bigint_t *a, mpfr_rnd_t rnd)
{
    size_t i;
    if (a->size == 0) { mpfr_set_zero(rop, +1); return; }
    /* Accumulate most-significant word first: rop = rop*2^32 + half */
    mpfr_set_zero(rop, +1);
    for (i = a->size; i-- > 0;) {
        mpfr_mul_2ui(rop, rop, 32, rnd);
        mpfr_add_ui(rop, rop, (unsigned long)(a->words[i] >> 32), rnd);
        mpfr_mul_2ui(rop, rop, 32, rnd);
        mpfr_add_ui(rop, rop, (unsigned long)(a->words[i] & 0xFFFFFFFFUL), rnd);
    }
}
#endif
