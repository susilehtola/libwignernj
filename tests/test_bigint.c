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
