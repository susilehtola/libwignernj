/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef BIGINT_H
#define BIGINT_H

#include <stdint.h>
#include <stddef.h>

/*
 * Arbitrary-precision non-negative integer stored as a little-endian array
 * of uint64_t words.  Sign is tracked externally by callers (wigner_exact_t).
 *
 * Invariant: words[size-1] != 0 unless size == 0 (representing zero).
 */
typedef struct {
    uint64_t *words;   /* little-endian base-2^64 digits */
    size_t    size;    /* number of significant words    */
    size_t    cap;     /* allocated capacity             */
} bigint_t;

/* Lifecycle */
void bigint_init(bigint_t *a);
void bigint_free(bigint_t *a);
void bigint_copy(bigint_t *dst, const bigint_t *src);
void bigint_set_u64(bigint_t *a, uint64_t v);
void bigint_set_zero(bigint_t *a);
int  bigint_is_zero(const bigint_t *a);

/* Comparison (unsigned magnitudes) */
int  bigint_cmp(const bigint_t *a, const bigint_t *b); /* -1, 0, +1 */

/* Arithmetic (all operands are non-negative; sign is caller's problem) */
void bigint_add(bigint_t *r, const bigint_t *a, const bigint_t *b);
/* r = |a| - |b|, returns sign (+1 if a>=b, -1 if a<b). r may alias a or b. */
int  bigint_sub_signed(bigint_t *r, const bigint_t *a, const bigint_t *b);
void bigint_mul(bigint_t *r, const bigint_t *a, const bigint_t *b);
void bigint_mul_u64(bigint_t *r, const bigint_t *a, uint64_t b);
/* Exact (no remainder) division by a small divisor. */
void bigint_div_u64(bigint_t *r, const bigint_t *a, uint64_t b);
/* Multiply a by p^k in-place (p is a prime that fits in uint64_t). */
void bigint_mul_prime_pow(bigint_t *a, uint64_t p, int k);
/* Exact division by p^k in-place. */
void bigint_div_prime_pow(bigint_t *a, uint64_t p, int k);

/* Bit operations */
int  bigint_bit_length(const bigint_t *a); /* floor(log2(a))+1, 0 for zero */

/* Conversion to floating-point.  IEEE 754 round-to-nearest-even. */
float       bigint_to_float(const bigint_t *a);
double      bigint_to_double(const bigint_t *a);
long double bigint_to_long_double(const bigint_t *a);

/*
 * Normalised frexp-style decomposition: returns m ∈ [0.5, 1) and sets
 * *out_exp such that  a ≈ m * 2^(*out_exp).
 * Returns 0.0 and *out_exp=0 for zero.  Rounds the mantissa to the
 * precision of the respective type (FLT_MANT_DIG / DBL_MANT_DIG / LDBL_MANT_DIG).
 * Used by wigner_exact_to_* to avoid intermediate overflow/underflow.
 */
float       bigint_frexpf(const bigint_t *a, int *out_exp);
double      bigint_frexp (const bigint_t *a, int *out_exp);
long double bigint_frexpl(const bigint_t *a, int *out_exp);

#endif /* BIGINT_H */
