/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "run_tests.h"
#include "../src/pfrac.h"
#include <math.h>

int main(void)
{
    pfrac_t f, g;
    bigint_t in, id, sn, sd;
    int i2, i3, i5, i;

    pfrac_init(&f); pfrac_init(&g);
    bigint_init(&in); bigint_init(&id);
    bigint_init(&sn); bigint_init(&sd);

    i2 = prime_index_of(2);
    i3 = prime_index_of(3);
    i5 = prime_index_of(5);
    TEST_ASSERT(i2 == 0);
    TEST_ASSERT(i3 == 1);
    TEST_ASSERT(i5 == 2);

    /* pfrac_mul_factorial: 6! = 720 = 2^4 * 3^2 * 5 */
    pfrac_zero(&f);
    pfrac_mul_factorial(&f, 6);
    TEST_ASSERT(f.exp[i2] == 4);
    TEST_ASSERT(f.exp[i3] == 2);
    TEST_ASSERT(f.exp[i5] == 1);

    /* round-trip: mul then div */
    pfrac_div_factorial(&f, 6);
    for (i = 0; i < g_nprimes; i++)
        TEST_ASSERT(f.exp[i] == 0);

    /* 10! = 2^8 * 3^4 * 5^2 * 7 */
    pfrac_zero(&f);
    pfrac_mul_factorial(&f, 10);
    TEST_ASSERT(f.exp[i2] == 8);
    TEST_ASSERT(f.exp[i3] == 4);
    TEST_ASSERT(f.exp[i5] == 2);
    TEST_ASSERT(f.exp[prime_index_of(7)] == 1);

    /* pfrac_mul_int: 12 = 2^2 * 3 */
    pfrac_zero(&f);
    pfrac_mul_int(&f, 12);
    TEST_ASSERT(f.exp[i2] == 2);
    TEST_ASSERT(f.exp[i3] == 1);
    TEST_ASSERT(f.exp[i5] == 0);

    /* pfrac_mul_int: 1 (no change) */
    pfrac_zero(&f);
    pfrac_mul_int(&f, 1);
    for (i = 0; i < g_nprimes; i++)
        TEST_ASSERT(f.exp[i] == 0);

    /* pfrac_mul_int: prime 97 */
    pfrac_zero(&f);
    pfrac_mul_int(&f, 97);
    TEST_ASSERT(f.exp[prime_index_of(97)] == 1);

    /* pfrac_copy */
    pfrac_zero(&f);
    f.exp[i2] = 5; f.exp[i3] = -3; f.max_idx = i3 + 1;
    pfrac_copy(&g, &f);
    TEST_ASSERT(g.exp[i2] == 5);
    TEST_ASSERT(g.exp[i3] == -3);
    pfrac_zero(&g);
    TEST_ASSERT(g.exp[i2] == 0);
    TEST_ASSERT(f.exp[i2] == 5);  /* f unchanged after g zeroed */

    /* ── pfrac_to_sqrt_rational ─────────────────────────────────────── */

    /* sqrt(2^4) = 4 : even positive exponent → integer part only */
    pfrac_zero(&f);
    f.exp[i2] = 4; f.max_idx = i2 + 1;
    bigint_set_u64(&in, 1); bigint_set_u64(&id, 1);
    bigint_set_u64(&sn, 1); bigint_set_u64(&sd, 1);
    pfrac_to_sqrt_rational(&f, &in, &id, &sn, &sd);
    TEST_ASSERT(bigint_to_double(&in) == 4.0);
    TEST_ASSERT(bigint_to_double(&id) == 1.0);
    TEST_ASSERT(bigint_to_double(&sn) == 1.0);
    TEST_ASSERT(bigint_to_double(&sd) == 1.0);

    /* sqrt(2^3) = 2*sqrt(2) : odd positive exponent */
    pfrac_zero(&f);
    f.exp[i2] = 3; f.max_idx = i2 + 1;
    bigint_set_u64(&in, 1); bigint_set_u64(&id, 1);
    bigint_set_u64(&sn, 1); bigint_set_u64(&sd, 1);
    pfrac_to_sqrt_rational(&f, &in, &id, &sn, &sd);
    TEST_ASSERT(bigint_to_double(&in) == 2.0);
    TEST_ASSERT(bigint_to_double(&id) == 1.0);
    TEST_ASSERT(bigint_to_double(&sn) == 2.0);
    TEST_ASSERT(bigint_to_double(&sd) == 1.0);

    /* sqrt(1/4) = 1/2 : even negative exponent */
    pfrac_zero(&f);
    f.exp[i2] = -2; f.max_idx = i2 + 1;
    bigint_set_u64(&in, 1); bigint_set_u64(&id, 1);
    bigint_set_u64(&sn, 1); bigint_set_u64(&sd, 1);
    pfrac_to_sqrt_rational(&f, &in, &id, &sn, &sd);
    TEST_ASSERT(bigint_to_double(&in) == 1.0);
    TEST_ASSERT(bigint_to_double(&id) == 2.0);
    TEST_ASSERT(bigint_to_double(&sn) == 1.0);
    TEST_ASSERT(bigint_to_double(&sd) == 1.0);

    /* sqrt(1/3) = 1/sqrt(3) : odd negative exponent → squarefree denominator */
    pfrac_zero(&f);
    f.exp[i3] = -1; f.max_idx = i3 + 1;
    bigint_set_u64(&in, 1); bigint_set_u64(&id, 1);
    bigint_set_u64(&sn, 1); bigint_set_u64(&sd, 1);
    pfrac_to_sqrt_rational(&f, &in, &id, &sn, &sd);
    TEST_ASSERT(bigint_to_double(&in) == 1.0);
    TEST_ASSERT(bigint_to_double(&id) == 1.0);
    TEST_ASSERT(bigint_to_double(&sn) == 1.0);
    TEST_ASSERT(bigint_to_double(&sd) == 3.0);

    /* sqrt(12) = 2*sqrt(3) : two primes, one even one odd */
    pfrac_zero(&f);
    f.exp[i2] = 2;
    f.exp[i3] = 1; f.max_idx = i3 + 1;
    bigint_set_u64(&in, 1); bigint_set_u64(&id, 1);
    bigint_set_u64(&sn, 1); bigint_set_u64(&sd, 1);
    pfrac_to_sqrt_rational(&f, &in, &id, &sn, &sd);
    TEST_ASSERT(bigint_to_double(&in) == 2.0);
    TEST_ASSERT(bigint_to_double(&id) == 1.0);
    TEST_ASSERT(bigint_to_double(&sn) == 3.0);
    TEST_ASSERT(bigint_to_double(&sd) == 1.0);
    {
        double val = bigint_to_double(&in) / bigint_to_double(&id)
                   * sqrt(bigint_to_double(&sn) / bigint_to_double(&sd));
        TEST_NEAR(val, sqrt(12.0), 1e-14);
    }

    /* sqrt(3/(2^5)) = sqrt(3/32) = sqrt(3)/(4*sqrt(2)) :
     * exp[i2]=-5 (odd neg), exp[i3]=1 (odd pos)
     * → int_den gets 2^2=4, sqrt_den gets 2; sqrt_num gets 3 */
    pfrac_zero(&f);
    f.exp[i2] = -5;
    f.exp[i3] =  1; f.max_idx = i3 + 1;
    bigint_set_u64(&in, 1); bigint_set_u64(&id, 1);
    bigint_set_u64(&sn, 1); bigint_set_u64(&sd, 1);
    pfrac_to_sqrt_rational(&f, &in, &id, &sn, &sd);
    TEST_ASSERT(bigint_to_double(&in) == 1.0);
    TEST_ASSERT(bigint_to_double(&id) == 4.0);
    TEST_ASSERT(bigint_to_double(&sn) == 3.0);
    TEST_ASSERT(bigint_to_double(&sd) == 2.0);
    {
        double val = bigint_to_double(&in) / bigint_to_double(&id)
                   * sqrt(bigint_to_double(&sn) / bigint_to_double(&sd));
        TEST_NEAR(val, sqrt(3.0/32.0), 1e-14);
    }

    pfrac_free(&f); pfrac_free(&g);
    bigint_free(&in); bigint_free(&id);
    bigint_free(&sn); bigint_free(&sd);
    SUMMARY();
}
