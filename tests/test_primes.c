/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "run_tests.h"
#include "../src/primes.h"

int main(void)
{
    primes_init();

    /* First few primes */
    TEST_ASSERT(g_primes[0] == 2);
    TEST_ASSERT(g_primes[1] == 3);
    TEST_ASSERT(g_primes[2] == 5);
    TEST_ASSERT(g_primes[3] == 7);
    TEST_ASSERT(g_primes[4] == 11);

    /* Prime count: pi(20011) >= 2262 */
    TEST_ASSERT(g_nprimes >= 2262);

    /* Index lookup */
    TEST_ASSERT(g_prime_index[2]  == 0);
    TEST_ASSERT(g_prime_index[3]  == 1);
    TEST_ASSERT(g_prime_index[4]  == -1);
    TEST_ASSERT(g_prime_index[97] >= 0);
    TEST_ASSERT(g_primes[g_prime_index[97]] == 97);

    /* Legendre valuations: known values */
    /* v_2(8!) = 7 */
    TEST_ASSERT(legendre_valuation(8, 0) == 7);
    /* v_2(10!) = 8 */
    TEST_ASSERT(legendre_valuation(10, 0) == 8);
    /* v_3(9!) = 4 */
    TEST_ASSERT(legendre_valuation(9, 1) == 4);
    /* v_5(25!) = 6 */
    TEST_ASSERT(legendre_valuation(25, 2) == 6);
    /* v_7(49!) = 8 (floor(49/7)+floor(49/49) = 7+1 = 8) */
    TEST_ASSERT(legendre_valuation(49, 3) == 8);

    /* v_p(0!) = 0, v_p(1!) = 0 */
    TEST_ASSERT(legendre_valuation(0, 0) == 0);
    TEST_ASSERT(legendre_valuation(1, 0) == 0);

    /* v_2(1000!) = 994 */
    TEST_ASSERT(legendre_valuation(1000, 0) == 994);

    /* Factorisation check: 10! = 2^8 * 3^4 * 5^2 * 7^1 */
    TEST_ASSERT(legendre_valuation(10, 0) == 8);  /* 2 */
    TEST_ASSERT(legendre_valuation(10, 1) == 4);  /* 3 */
    TEST_ASSERT(legendre_valuation(10, 2) == 2);  /* 5 */
    TEST_ASSERT(legendre_valuation(10, 3) == 1);  /* 7 */
    TEST_ASSERT(legendre_valuation(10, 4) == 0);  /* 11 > 10 */

    SUMMARY();
}
