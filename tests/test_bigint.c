/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "run_tests.h"
#include "../src/bigint.h"
#include <float.h>

int main(void)
{
    bigint_t a, b, r;
    bigint_init(&a); bigint_init(&b); bigint_init(&r);

    /* zero */
    bigint_set_zero(&a);
    TEST_ASSERT(bigint_is_zero(&a));
    TEST_ASSERT(bigint_bit_length(&a) == 0);

    /* set and compare */
    bigint_set_u64(&a, 12);
    bigint_set_u64(&b, 8);
    TEST_ASSERT(!bigint_is_zero(&a));
    TEST_ASSERT(bigint_cmp(&a, &b) > 0);
    TEST_ASSERT(bigint_cmp(&b, &a) < 0);

    /* addition */
    bigint_add(&r, &a, &b);
    TEST_ASSERT(bigint_to_double(&r) == 20.0);

    /* subtraction */
    bigint_sub_signed(&r, &a, &b);
    TEST_ASSERT(bigint_to_double(&r) == 4.0);

    /* mul_u64 */
    bigint_set_u64(&a, 100);
    bigint_mul_u64(&r, &a, 200);
    TEST_ASSERT(bigint_to_double(&r) == 20000.0);

    /* mul (schoolbook) */
    bigint_set_u64(&a, 123456789ULL);
    bigint_set_u64(&b, 987654321ULL);
    bigint_mul(&r, &a, &b);
    TEST_DBL(bigint_to_double(&r), 123456789.0 * 987654321.0);

    /* mul: multi-word × multi-word → multi-word (2^65 * 2^65 = 2^130) */
    bigint_set_u64(&a, 1);
    bigint_mul_prime_pow(&a, 2, 65);        /* a = 2^65 */
    bigint_set_u64(&b, 1);
    bigint_mul_prime_pow(&b, 2, 65);        /* b = 2^65 */
    bigint_mul(&r, &a, &b);                 /* r = 2^130 */
    TEST_NEAR(bigint_to_double(&r), pow(2.0, 130.0), 1e-14);
    TEST_ASSERT(bigint_bit_length(&r) == 131);

    /* mul: zero × anything = zero */
    bigint_set_zero(&a);
    bigint_set_u64(&b, 42);
    bigint_mul(&r, &a, &b);
    TEST_ASSERT(bigint_is_zero(&r));

    /* mul_prime_pow / div_prime_pow round-trip */
    bigint_set_u64(&a, 1);
    bigint_mul_prime_pow(&a, 2, 10);   /* a = 2^10 = 1024 */
    TEST_ASSERT(bigint_to_double(&a) == 1024.0);
    bigint_div_prime_pow(&a, 2, 10);
    TEST_ASSERT(bigint_to_double(&a) == 1.0);

    /* large prime power */
    bigint_set_u64(&a, 1);
    bigint_mul_prime_pow(&a, 3, 40);   /* a = 3^40 */
    {
        double expected = pow(3.0, 40.0);
        TEST_NEAR(bigint_to_double(&a), expected, 1e-12);
    }

    /* div_u64 */
    bigint_set_u64(&a, 100);
    bigint_div_u64(&r, &a, 4);
    TEST_ASSERT(bigint_to_double(&r) == 25.0);

    bigint_set_u64(&a, 1000000ULL);
    bigint_div_u64(&r, &a, 1000);
    TEST_ASSERT(bigint_to_double(&r) == 1000.0);

    /* div_u64: multi-word dividend (2^64 / 2 = 2^63) */
    bigint_set_u64(&a, 1);
    bigint_mul_prime_pow(&a, 2, 64);  /* a = 2^64 */
    bigint_div_u64(&r, &a, 2);
    TEST_NEAR(bigint_to_double(&r), pow(2.0, 63.0), 1e-14);

    /* bigint_div_u64_exact: Hensel-style exact division.  Compare
     * against the generic bigint_div_u64 across a range of dividend
     * sizes and 64-bit divisor shapes.  Both routines must produce
     * identical output for any dividend that is a multiple of d. */
    {
        const uint64_t divisors[] = {
            3ULL,                             /* small odd */
            6ULL,                             /* small even (2*3) */
            123456789ULL,                     /* 32-bit ish */
            0x100000007ULL,                   /* > 2^32, prime-like */
            0xDEADBEEFCAFEBABEULL,            /* > 2^32, mixed */
            0xFFFFFFFFFFFFFFFFULL,            /* 2^64-1 (odd) */
            0x8000000000000000ULL,            /* 2^63 (high power of 2) */
        };
        size_t n_div = sizeof(divisors) / sizeof(divisors[0]);
        size_t k;
        bigint_t base, prod, q1, q2;
        bigint_init(&base); bigint_init(&prod);
        bigint_init(&q1);   bigint_init(&q2);
        /* Build a deliberately multi-word base, then for each divisor
         * compute (base * d) and verify (base * d) / d == base via
         * both code paths. */
        bigint_set_u64(&base, 0xCAFEBABE12345678ULL);
        bigint_mul_prime_pow(&base, 2, 50);
        bigint_mul_u64(&base, &base, 999983);   /* mid-32-bit prime */
        for (k = 0; k < n_div; k++) {
            bigint_mul_u64(&prod, &base, divisors[k]);
            bigint_div_u64      (&q1, &prod, divisors[k]);
            bigint_div_u64_exact(&q2, &prod, divisors[k]);
            TEST_ASSERT(bigint_cmp(&q1, &base) == 0);
            TEST_ASSERT(bigint_cmp(&q2, &base) == 0);
        }
        /* In-place aliasing: r == a must work. */
        bigint_mul_u64(&prod, &base, 0x100000007ULL);
        bigint_div_u64_exact(&prod, &prod, 0x100000007ULL);
        TEST_ASSERT(bigint_cmp(&prod, &base) == 0);
        bigint_free(&base); bigint_free(&prod);
        bigint_free(&q1);   bigint_free(&q2);
    }

    /* div_u64 with a 64-bit divisor (exercises the Algorithm-D path
     * of bigint_div128 on every backend that doesn't have a hardware
     * 128/64 instruction; the x86-64 hardware divq backend handles
     * arbitrary 64-bit divisors directly).  Compute 7 * D / D = 7
     * for D = 0xDEADBEEFCAFEBABE (a representative > 2^32 value). */
    {
        const uint64_t big_d = 0xDEADBEEFCAFEBABEULL;
        bigint_set_u64(&a, big_d);
        bigint_mul_u64(&a, &a, 7);
        bigint_div_u64(&r, &a, big_d);
        TEST_ASSERT(bigint_bit_length(&r) > 0);
        TEST_ASSERT(bigint_to_double(&r) == 7.0);
    }

    /* div_u64 with the largest 64-bit divisor (2^64-1).  Stress-tests
     * the trial-quotient saturation branch of Algorithm D: the high
     * 32 bits of the normalised dividend equal d1, forcing q_hat to
     * the cap 2^32-1, then refinement and add-back as needed. */
    {
        const uint64_t big_d = 0xFFFFFFFFFFFFFFFFULL;
        bigint_set_u64(&a, big_d);
        bigint_mul_u64(&a, &a, 13);
        bigint_div_u64(&r, &a, big_d);
        TEST_ASSERT(bigint_to_double(&r) == 13.0);
    }

    /* div_u64 with a 64-bit divisor whose product with a multi-word
     * dividend yields a non-trivial 96-bit running numerator at the
     * Algorithm-D long-division step boundary.  Computes
     * (2^96 + 1) / 0x100000001 (= 2^32 + 1), which has a known closed
     * form: floor((2^96 + 1) / (2^32 + 1)) = 2^64 - 2^32 + 1, with
     * remainder 0. */
    {
        const uint64_t big_d = 0x100000001ULL;
        bigint_t two_to_96, one;
        bigint_init(&two_to_96); bigint_init(&one);
        bigint_set_u64(&two_to_96, 1);
        bigint_mul_prime_pow(&two_to_96, 2, 96);
        bigint_set_u64(&one, 1);
        bigint_add(&a, &two_to_96, &one);                /* a = 2^96 + 1 */
        bigint_div_u64(&r, &a, big_d);
        /* Expected: 2^64 - 2^32 + 1.  Check via a round-trip. */
        bigint_mul_u64(&b, &r, big_d);
        TEST_ASSERT(bigint_cmp(&a, &b) == 0);
        bigint_free(&two_to_96); bigint_free(&one);
    }

    /* bigint_copy */
    bigint_set_u64(&a, 9999999999ULL);
    bigint_copy(&r, &a);
    TEST_ASSERT(bigint_to_double(&r) == 9999999999.0);
    bigint_set_u64(&a, 1);           /* modify a, r must be independent */
    TEST_ASSERT(bigint_to_double(&r) == 9999999999.0);

    /* bigint_to_long_double */
    bigint_set_u64(&a, 1);
    bigint_mul_prime_pow(&a, 2, 80);   /* a = 2^80 */
    TEST_NEAR((double)bigint_to_long_double(&a), pow(2.0, 80.0), 1e-14);

    /* bit_length */
    bigint_set_u64(&a, 1);
    TEST_ASSERT(bigint_bit_length(&a) == 1);
    bigint_set_u64(&a, 8);
    TEST_ASSERT(bigint_bit_length(&a) == 4);
    bigint_set_u64(&a, UINT64_MAX);
    TEST_ASSERT(bigint_bit_length(&a) == 64);

    /* to_double: exact for small values */
    bigint_set_u64(&a, 1ULL << 52);
    TEST_ASSERT(bigint_to_double(&a) == (double)(1ULL << 52));

    /* to_float */
    bigint_set_u64(&a, 1234567);
    TEST_NEAR(bigint_to_float(&a), 1234567.0f, 1e-6);

    /* large value: 2^100 */
    bigint_set_u64(&a, 1);
    bigint_mul_prime_pow(&a, 2, 100);
    TEST_NEAR(bigint_to_double(&a), pow(2.0, 100.0), 1e-12);

    /* factorial 10! = 3628800 */
    {
        int i;
        bigint_set_u64(&a, 1);
        for (i = 2; i <= 10; i++) bigint_mul_u64(&a, &a, (uint64_t)i);
        TEST_ASSERT(bigint_to_double(&a) == 3628800.0);
    }

    /* round-to-nearest-even: 2^54 - 1 ties between 2^54-2 and 2^54;
     * 2^54-2 has odd stored mantissa (all 52 bits set), 2^54 is even → rounds up */
    {
        bigint_set_u64(&a, 1);
        bigint_mul_prime_pow(&a, 2, 54);  /* a = 2^54 */
        bigint_set_u64(&b, 1);
        bigint_sub_signed(&r, &a, &b);    /* r = 2^54 - 1 */
        TEST_ASSERT(bigint_to_double(&r) == pow(2.0, 54.0));
    }

    /* round-to-nearest-even: 2^53 + 1 ties between 2^53 and 2^53+2;
     * 2^53 has even stored mantissa (all zeros) → rounds down */
    {
        bigint_set_u64(&a, (uint64_t)1 << 53);
        bigint_set_u64(&b, 1);
        bigint_add(&r, &a, &b);           /* r = 2^53 + 1 */
        TEST_ASSERT(bigint_to_double(&r) == pow(2.0, 53.0));
    }

    /* multi-word: 2^64 * 7 + 3  (two words) */
    {
        bigint_set_u64(&a, 1);
        bigint_mul_prime_pow(&a, 2, 64);  /* a = 2^64 */
        bigint_mul_u64(&a, &a, 7);        /* a = 7 * 2^64 */
        bigint_set_u64(&b, 3);
        bigint_add(&r, &a, &b);           /* r = 7 * 2^64 + 3 */
        TEST_NEAR(bigint_to_double(&r), 7.0 * pow(2.0, 64.0) + 3.0, 1e-9);
        TEST_ASSERT(r.size == 2);
    }

    /* bigint_frexp: 2^63 → mantissa 0.5, exp 64 */
    {
        int exp2;
        double m;
        bigint_set_u64(&a, (uint64_t)1 << 63);
        m = bigint_frexp(&a, &exp2);
        TEST_NEAR(m, 0.5, 1e-15);
        TEST_ASSERT(exp2 == 64);
    }

    /* bigint_frexpf: 2^32 → mantissa 0.5f, exp 33 */
    {
        int exp2;
        float m;
        bigint_set_u64(&a, (uint64_t)1 << 32);
        m = bigint_frexpf(&a, &exp2);
        TEST_NEAR((double)m, 0.5, 1e-6);
        TEST_ASSERT(exp2 == 33);
    }

    /* bigint_frexpl: 2^80 → mantissa 0.5L, exp 81 */
    {
        int exp2;
        long double m;
        bigint_set_u64(&a, 1);
        bigint_mul_prime_pow(&a, 2, 80);
        m = bigint_frexpl(&a, &exp2);
        TEST_NEAR((double)m, 0.5, 1e-14);
        TEST_ASSERT(exp2 == 81);
    }

    /* frexp invariant: a == m * 2^exp */
    {
        int exp2;
        double m;
        bigint_set_u64(&a, 123456789012ULL);
        m = bigint_frexp(&a, &exp2);
        TEST_NEAR(m * pow(2.0, (double)exp2), 123456789012.0, 1e-10);
    }

    /* extract_mantissa mant=64 rounding overflow: 2^65-1 rounds to 2^65.
     * Top 64 bits = 0xFFFF...FFFF, round bit = 1 → m++ wraps to 0.
     * Regression: frexpl must return 0.5 (not 0) with exp=66. */
    {
        int exp2;
        long double m;
#if LDBL_MANT_DIG >= 64
        bigint_set_u64(&a, 1);
        bigint_mul_prime_pow(&a, 2, 65);   /* a = 2^65 */
        bigint_set_u64(&b, 1);
        bigint_sub_signed(&r, &a, &b);     /* r = 2^65 - 1 */
        m = bigint_frexpl(&r, &exp2);
        TEST_NEAR((double)m, 0.5, 1e-14); /* must not be 0 */
        TEST_ASSERT(exp2 == 66);
#endif
    }

    /* bigint_add aliasing: r == a (Racah sum accumulation pattern) */
    bigint_set_u64(&a, 100);
    bigint_set_u64(&b, 200);
    bigint_add(&a, &a, &b);  /* a = 100 + 200 = 300, result aliases first input */
    TEST_ASSERT(bigint_to_double(&a) == 300.0);

    /* bigint_add aliasing: r == b */
    bigint_set_u64(&a, 7);
    bigint_set_u64(&b, 13);
    bigint_add(&b, &a, &b);  /* b = 7 + 13 = 20 */
    TEST_ASSERT(bigint_to_double(&b) == 20.0);

    /* bigint_sub_signed aliasing: r == a */
    bigint_set_u64(&a, 500);
    bigint_set_u64(&b, 123);
    bigint_sub_signed(&a, &a, &b);  /* a = 500 - 123 = 377 */
    TEST_ASSERT(bigint_to_double(&a) == 377.0);

    bigint_free(&a); bigint_free(&b); bigint_free(&r);
    SUMMARY();
}
