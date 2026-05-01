/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * 64-bit and 128-bit arithmetic primitives are abstracted into
 * bigint_arith.h, which provides a native __uint128_t-based path on
 * GCC/Clang/ICC and a pure-C99 fallback on other compilers (e.g. MSVC).
 */
#include "bigint.h"
#include "bigint_arith.h"
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

void bigint_reserve(bigint_t *a, size_t cap)
{
    bigint_ensure(a, cap);
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

/*
 * The four "linear" routines below (add, sub_uu, mul_u64, div_u64) all use
 * a read-then-write pattern at index i: each iteration first reads the
 * operand word(s) at i and only then writes r->words[i].  This makes them
 * correct even when r aliases an operand, so we can write the result
 * directly into r and avoid the per-call temp + copy that the previous
 * implementation always performed.
 *
 * bigint_mul cannot be aliased away as easily: its inner loop reads
 * operand words ahead of where it writes, so when r aliases a or b a
 * separate destination is still required.  We use one in that case only.
 */

/* ── addition ────────────────────────────────────────────────────────────── */

void bigint_add(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    const bigint_t *big = (a->size >= b->size) ? a : b;
    const bigint_t *sml = (a->size >= b->size) ? b : a;
    size_t n = big->size, i;
    uint64_t carry = 0;
    if (n == 0) { bigint_set_zero(r); return; }
    bigint_ensure(r, n + 1);
    for (i = 0; i < sml->size; i++)
        r->words[i] = bigint_addc(big->words[i], sml->words[i], carry, &carry);
    for (; i < n; i++)
        r->words[i] = bigint_addc(big->words[i], 0, carry, &carry);
    r->words[n] = carry;
    r->size = n + 1;
    bigint_trim(r);
}

/* ── subtraction (unsigned, a >= b assumed) ─────────────────────────────── */

static void bigint_sub_uu(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    uint64_t borrow = 0;
    size_t i;
    if (a->size == 0) { bigint_set_zero(r); return; }
    bigint_ensure(r, a->size);
    for (i = 0; i < b->size; i++)
        r->words[i] = bigint_subb(a->words[i], b->words[i], borrow, &borrow);
    for (; i < a->size; i++)
        r->words[i] = bigint_subb(a->words[i], 0, borrow, &borrow);
    r->size = a->size;
    bigint_trim(r);
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
    if (b == 0 || a->size == 0) { bigint_set_zero(r); return; }
    bigint_ensure(r, a->size + 1);
    for (i = 0; i < a->size; i++) {
        uint64_t lo, hi;
        bigint_mul_add(a->words[i], b, carry, &lo, &hi);
        r->words[i] = lo;
        carry       = hi;
    }
    r->words[a->size] = carry;
    r->size = a->size + 1;
    bigint_trim(r);
}

void bigint_mul(bigint_t *r, const bigint_t *a, const bigint_t *b)
{
    size_t i, j;
    bigint_t local;
    bigint_t *dest;
    int aliased;
    if (a->size == 0 || b->size == 0) { bigint_set_zero(r); return; }
    /* The inner loop reads res.words[i+j] *after* writes to lower indices
     * within the same outer-i sweep have moved on, so destination aliasing
     * with an operand corrupts that operand mid-loop.  Only use a temp
     * when the destination is actually aliased. */
    aliased = (r == a || r == b);
    if (aliased) {
        bigint_init(&local);
        bigint_ensure(&local, a->size + b->size);
        memset(local.words, 0, (a->size + b->size) * sizeof(uint64_t));
        local.size = a->size + b->size;
        dest = &local;
    } else {
        bigint_ensure(r, a->size + b->size);
        memset(r->words, 0, (a->size + b->size) * sizeof(uint64_t));
        r->size = a->size + b->size;
        dest = r;
    }
    for (i = 0; i < a->size; i++) {
        uint64_t carry = 0;
        for (j = 0; j < b->size; j++) {
            uint64_t lo, hi;
            bigint_mul_add2(a->words[i], b->words[j],
                            dest->words[i+j], carry, &lo, &hi);
            dest->words[i+j] = lo;
            carry            = hi;
        }
        dest->words[i + b->size] += carry;
    }
    bigint_trim(dest);
    if (aliased) {
        bigint_copy(r, &local);
        bigint_free(&local);
    }
}

/* ── division by scalar ──────────────────────────────────────────────────── */

void bigint_div_u64(bigint_t *r, const bigint_t *a, uint64_t b)
{
    size_t i;
    uint64_t rem = 0;
    if (a->size == 0) { bigint_set_zero(r); return; }
    bigint_ensure(r, a->size);
    for (i = a->size; i-- > 0;)
        r->words[i] = bigint_div128(rem, a->words[i], b, &rem);
    r->size = a->size;
    bigint_trim(r);
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

/* ── workspace-aware variants ────────────────────────────────────────────── */

void bigint_ws_init(bigint_ws_t *ws)
{
    bigint_init(&ws->mul_temp);
    bigint_init(&ws->pp_pw);
    bigint_init(&ws->pp_base);
}

void bigint_ws_free(bigint_ws_t *ws)
{
    bigint_free(&ws->mul_temp);
    bigint_free(&ws->pp_pw);
    bigint_free(&ws->pp_base);
}

void bigint_ws_reserve(bigint_ws_t *ws, size_t max_words)
{
    bigint_ensure(&ws->mul_temp, max_words);
    bigint_ensure(&ws->pp_pw,    max_words);
    bigint_ensure(&ws->pp_base,  max_words);
}

void bigint_mul_ws(bigint_t *r, const bigint_t *a, const bigint_t *b,
                   bigint_ws_t *ws)
{
    size_t i, j;
    bigint_t *dest;
    int aliased;
    if (a->size == 0 || b->size == 0) { bigint_set_zero(r); return; }
    aliased = (r == a || r == b);
    if (aliased) {
        bigint_ensure(&ws->mul_temp, a->size + b->size);
        memset(ws->mul_temp.words, 0, (a->size + b->size) * sizeof(uint64_t));
        ws->mul_temp.size = a->size + b->size;
        dest = &ws->mul_temp;
    } else {
        bigint_ensure(r, a->size + b->size);
        memset(r->words, 0, (a->size + b->size) * sizeof(uint64_t));
        r->size = a->size + b->size;
        dest = r;
    }
    for (i = 0; i < a->size; i++) {
        uint64_t carry = 0;
        for (j = 0; j < b->size; j++) {
            uint64_t lo, hi;
            bigint_mul_add2(a->words[i], b->words[j],
                            dest->words[i+j], carry, &lo, &hi);
            dest->words[i+j] = lo;
            carry            = hi;
        }
        dest->words[i + b->size] += carry;
    }
    bigint_trim(dest);
    if (aliased) bigint_copy(r, dest);
}

void bigint_mul_prime_pow_ws(bigint_t *a, uint64_t p, int k, bigint_ws_t *ws)
{
    int i;
    if (k <= 0) return;
    if (k <= 8) {
        for (i = 0; i < k; i++) bigint_mul_u64(a, a, p);
        return;
    }
    /* Binary exponentiation using ws->pp_pw and ws->pp_base as scratch.
     * No allocations occur if ws was reserved to the worst-case size. */
    bigint_set_u64(&ws->pp_pw,   1);
    bigint_set_u64(&ws->pp_base, p);
    while (k > 0) {
        if (k & 1) bigint_mul_ws(&ws->pp_pw,   &ws->pp_pw,   &ws->pp_base, ws);
        k >>= 1;
        if (k > 0)  bigint_mul_ws(&ws->pp_base, &ws->pp_base, &ws->pp_base, ws);
    }
    bigint_mul_ws(a, a, &ws->pp_pw, ws);
}

size_t bigint_words_for_factorial(int N)
{
    /* Stirling: log2(N!) ≈ N*log2(N) - N*log2(e) + 0.5*log2(N) + log2(sqrt(2π))
     *
     * The returned size is generous: we double the Stirling bound so that
     * any LCM-sized bigint, plus the cross-product temporaries it produces
     * in bigint_mul_ws (whose size is operand1.size + operand2.size words),
     * fits without forcing the workspace to realloc inside the hot loop.
     */
    double bits, log2N;
    size_t base;
    if (N <= 1) return 8;
    log2N = log((double)N) * 1.4426950408889634;
    bits  = (double)N * log2N
          - (double)N * 1.4426950408889634
          + 0.5 * log2N
          + 1.33
          + 64.0;
    if (bits < 64.0) bits = 64.0;
    base = (size_t)(bits / 64.0) + 4;
    return 2 * base;  /* headroom for mid-multiplication temporaries */
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
