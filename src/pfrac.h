/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef PFRAC_H
#define PFRAC_H

#include "bigint.h"
#include "primes.h"

/*
 * Prime-factored exact representation.
 *
 * A pfrac_t stores, for each prime p_i in g_primes[], the signed exponent
 * of p_i in the represented value.  Positive exponents are in the numerator;
 * negative exponents are in the denominator.
 *
 * Two distinct uses:
 *   - "rational" pfrac: represents an exact rational number q = prod p_i^exp[i].
 *   - "sqrt" pfrac:     represents sqrt(q) = prod p_i^(exp[i]/2).
 *     In both cases the same struct is used; the interpretation is the caller's.
 *
 * exp[] is sized to g_nprimes and zero-initialised at allocation.
 */
typedef struct {
    int *exp;   /* exp[i] = exponent of g_primes[i]; len = g_nprimes */
} pfrac_t;

void pfrac_init(pfrac_t *f);
void pfrac_free(pfrac_t *f);
void pfrac_zero(pfrac_t *f);          /* set all exponents to 0 */
void pfrac_copy(pfrac_t *dst, const pfrac_t *src);

/* Multiply the rational by n! (adds Legendre valuations to exp[]).  */
void pfrac_mul_factorial(pfrac_t *f, int n);

/* Divide the rational by n! (subtracts Legendre valuations from exp[]). */
void pfrac_div_factorial(pfrac_t *f, int n);

/* Multiply the rational by the positive integer k.
 * k must have no prime factor exceeding PRIME_SIEVE_LIMIT (20011). */
void pfrac_mul_int(pfrac_t *f, int k);

/*
 * Convert a pfrac_t representing sqrt(Q) into four bigint components such
 * that   sqrt(Q) = (int_num / int_den) * sqrt(sqrt_num / sqrt_den).
 *
 * For each prime p_i with exponent e = f->exp[i]:
 *   If e > 0:  p_i^floor(e/2) → int_num;  if e is odd, p_i → sqrt_num.
 *   If e < 0:  p_i^floor(|e|/2) → int_den; if |e| is odd, p_i → sqrt_den.
 *   If e == 0: skipped.
 *
 * All four bigints are initialised by the caller; this function multiplies
 * into them (does not reset to 1).  The caller must initialise them to 1
 * (bigint_set_u64(x, 1)) before calling.
 */
void pfrac_to_sqrt_rational(const pfrac_t *f,
                             bigint_t *int_num,  bigint_t *int_den,
                             bigint_t *sqrt_num, bigint_t *sqrt_den);

#endif /* PFRAC_H */
