/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * 64-bit arithmetic primitives used by bigint.c.
 *
 * Two implementations are selected at compile time:
 *
 *   - When __SIZEOF_INT128__ is defined (GCC, Clang, ICC, and most other
 *     Unix compilers), the primitives use the __uint128_t extension for
 *     64x64 -> 128 products and 128/64 divisions.
 *
 *   - Otherwise, a pure-C99 fallback is used: 64x64 products are formed
 *     from four 32x32 -> 64 partial products combined with explicit carry
 *     tracking, and 128/64 divisions (whose divisor in this code is always
 *     a prime <= PRIME_SIEVE_LIMIT and therefore fits in 32 bits) are
 *     implemented as two 64/32 long-division steps.
 *
 * Both code paths are bit-exact and produce identical results.
 */
#ifndef BIGINT_ARITH_H
#define BIGINT_ARITH_H

#include <stdint.h>

/* Defining BIGINT_FORCE_PORTABLE at compile time forces use of the pure-C99
 * fallback path even on compilers that support __uint128_t.  Useful for
 * exercising the fallback in CI. */
#if defined(__SIZEOF_INT128__) && !defined(BIGINT_FORCE_PORTABLE)
#  define BIGINT_HAVE_UINT128 1
typedef unsigned __int128 bigint_u128_t;
#else
#  define BIGINT_HAVE_UINT128 0
#endif

/* {*plo, *phi} = a * b */
static inline void bigint_mul64x64(uint64_t a, uint64_t b,
                                    uint64_t *plo, uint64_t *phi)
{
#if BIGINT_HAVE_UINT128
    bigint_u128_t p = (bigint_u128_t)a * b;
    *plo = (uint64_t)p;
    *phi = (uint64_t)(p >> 64);
#else
    uint32_t a0 = (uint32_t)a, a1 = (uint32_t)(a >> 32);
    uint32_t b0 = (uint32_t)b, b1 = (uint32_t)(b >> 32);
    uint64_t p00 = (uint64_t)a0 * b0;
    uint64_t p01 = (uint64_t)a0 * b1;
    uint64_t p10 = (uint64_t)a1 * b0;
    uint64_t p11 = (uint64_t)a1 * b1;
    /* Combine partial products: a*b = p11 * 2^64 + (p10+p01) * 2^32 + p00 */
    uint64_t mid = (p00 >> 32) + (uint32_t)p01 + (uint32_t)p10;
    *plo = ((uint32_t)p00) | (mid << 32);
    *phi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
#endif
}

/* {*plo, *phi} = a * b + c */
static inline void bigint_mul_add(uint64_t a, uint64_t b, uint64_t c,
                                   uint64_t *plo, uint64_t *phi)
{
#if BIGINT_HAVE_UINT128
    bigint_u128_t p = (bigint_u128_t)a * b + c;
    *plo = (uint64_t)p;
    *phi = (uint64_t)(p >> 64);
#else
    uint64_t lo, hi;
    bigint_mul64x64(a, b, &lo, &hi);
    *plo = lo + c;
    *phi = hi + (*plo < lo);
#endif
}

/* {*plo, *phi} = a * b + c + d */
static inline void bigint_mul_add2(uint64_t a, uint64_t b,
                                    uint64_t c, uint64_t d,
                                    uint64_t *plo, uint64_t *phi)
{
#if BIGINT_HAVE_UINT128
    bigint_u128_t p = (bigint_u128_t)a * b + c + d;
    *plo = (uint64_t)p;
    *phi = (uint64_t)(p >> 64);
#else
    uint64_t lo, hi, t;
    bigint_mul64x64(a, b, &lo, &hi);
    t = lo + c;
    hi += (t < lo);
    *plo = t + d;
    *phi = hi + (*plo < t);
#endif
}

/* a + b + cin.  Returns the sum and writes the carry-out (0 or 1) to *cout. */
static inline uint64_t bigint_addc(uint64_t a, uint64_t b,
                                    uint64_t cin, uint64_t *cout)
{
    uint64_t s = a + b;
    uint64_t c1 = (s < a);
    uint64_t r = s + cin;
    uint64_t c2 = (r < s);
    *cout = c1 | c2;
    return r;
}

/* a - b - bin.  Returns the difference and writes the borrow-out (0 or 1)
 * to *bout. */
static inline uint64_t bigint_subb(uint64_t a, uint64_t b,
                                    uint64_t bin, uint64_t *bout)
{
    uint64_t d = a - b;
    uint64_t b1 = (a < b);
    uint64_t r = d - bin;
    uint64_t b2 = (d < bin);
    *bout = b1 | b2;
    return r;
}

/*
 * Divide [hi:lo] (128 bits, with hi < d) by d.  Returns the 64-bit quotient
 * and writes the remainder (< d) to *rem.
 *
 * On the C99 fallback path d must satisfy d < 2^32, otherwise the
 * intermediate (r1 << 32) below could overflow.  In libwignernj the divisor
 * is always a prime <= PRIME_SIEVE_LIMIT (20011), so this is always
 * satisfied.
 */
static inline uint64_t bigint_div128(uint64_t hi, uint64_t lo,
                                      uint64_t d, uint64_t *rem)
{
#if BIGINT_HAVE_UINT128
    bigint_u128_t numer = ((bigint_u128_t)hi << 64) | lo;
    *rem = (uint64_t)(numer % d);
    return (uint64_t)(numer / d);
#else
    uint64_t numer, q1, r1, q2;
    /* Two-step 64/32 long division. */
    numer = (hi << 32) | (lo >> 32);
    q1 = numer / d;
    r1 = numer % d;
    numer = (r1 << 32) | (lo & 0xFFFFFFFFu);
    q2 = numer / d;
    *rem = numer % d;
    return (q1 << 32) | (q2 & 0xFFFFFFFFu);
#endif
}

#endif /* BIGINT_ARITH_H */
