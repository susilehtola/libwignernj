/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Microbench: Möller-Granlund "improved division by invariant integers"
 * (IEEE TC 60(2), 2011) vs the existing bigint_div128 (inline-asm divq
 * on x86-64; Algorithm D otherwise).
 *
 * The bench mirrors the inner loop of bigint_div_u64: sweep m limbs,
 * each step doing one 128/64 division with the previous limb's
 * remainder carried as `hi`.  The reciprocal is computed ONCE before
 * the sweep, matching the real call site.
 *
 * Build:
 *   gcc -O3 -I src benchmarks/bench_div128_gm.c -o bench_div128_gm
 *   gcc -O3 -DBIGINT_NO_DIVQ      -I src benchmarks/bench_div128_gm.c -o bench_div128_gm_nodivq
 *   gcc -O3 -DBIGINT_FORCE_PORTABLE -I src benchmarks/bench_div128_gm.c -o bench_div128_gm_portable
 *
 * Output: ns per limb-division for each backend, swept across a few
 * representative divisors and a few limb counts m matching realistic
 * 3j Pass-2 scratch sizes.
 */
#include "bigint_arith.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ── Möller-Granlund primitives ─────────────────────────────────────────── */

/* Precomputed reciprocal of a 64-bit divisor.
 *
 *   d_norm = d << s   where s = leading-zero count of d  (so d_norm's top bit is set)
 *   v      = floor((B^2 - 1) / d_norm) - B               where B = 2^64
 *
 * The shift s is needed to renormalise the dividend [hi:lo] before the
 * recip-based division and to shift the remainder back at the end.
 */
typedef struct {
    uint64_t v;        /* (B^2 - 1) / d_norm - B           */
    uint64_t d_norm;   /* d << s                            */
    int      s;        /* leading-zero count of d (0..63)   */
} gm_recip_t;

/* Count leading zeros for 64-bit; assumes d != 0. */
static inline int gm_clz64(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(x);
#else
    int n = 0;
    if ((x >> 32) == 0) { n += 32; x <<= 32; }
    if ((x >> 48) == 0) { n += 16; x <<= 16; }
    if ((x >> 56) == 0) { n +=  8; x <<=  8; }
    if ((x >> 60) == 0) { n +=  4; x <<=  4; }
    if ((x >> 62) == 0) { n +=  2; x <<=  2; }
    if ((x >> 63) == 0) { n +=  1; }
    return n;
#endif
}

/* Initialise reciprocal.  d > 0. */
static inline gm_recip_t gm_recip_init(uint64_t d)
{
    gm_recip_t r;
    uint64_t d_norm, v, rem;
    r.s      = gm_clz64(d);
    d_norm   = d << r.s;
    r.d_norm = d_norm;

    /* v = floor((B^2 - 1) / d_norm) - B
     *   = floor((~0 - d_norm) / d_norm) + something tricky? No --
     * exact form per Möller-Granlund: build it from a 128/64 division
     * of (B^2 - 1) = {2^64 - 1, 2^64 - 1} by d_norm.  Since d_norm has
     * its MSB set, the high half (2^64 - 1) >= d_norm, so we'd violate
     * the divq precondition hi < d.  Use the rearrangement:
     *
     *   v = floor( (~0ULL - d_norm) * 2^64 + ~0ULL ) / d_norm
     *
     * which fits hi < d_norm.  The closed form:
     *   v = ((B - d_norm) * B + (B - 1)) / d_norm
     * is what bigint_div128 computes when called with
     *   hi = ~0ULL - d_norm, lo = ~0ULL.
     */
    v = bigint_div128(~0ULL - d_norm, ~0ULL, d_norm, &rem);
    (void)rem;
    r.v = v;
    return r;
}

/* Divide [hi:lo] by d (the original, unnormalised divisor) using the
 * precomputed reciprocal.  hi < d required, just like bigint_div128.
 * Returns the 64-bit quotient and writes the remainder to *rem.
 *
 * Algorithm (Möller-Granlund Fig. 4 / Algorithm 4, adapted to operate
 * on the unnormalised [hi:lo] by shifting both the dividend and the
 * remainder around the normalised core).
 */
static inline uint64_t gm_div128(uint64_t hi, uint64_t lo,
                                  uint64_t d, gm_recip_t r,
                                  uint64_t *rem_out)
{
    uint64_t u1, u0, q1, q0, r_n;
    (void)d;

    /* Renormalise dividend by the same shift applied to the divisor.
     * After shifting, the top bit of the (now 128-bit) dividend may be
     * spread between u1 and u0; u1 < d_norm holds because hi < d. */
    if (r.s == 0) {
        u1 = hi;
        u0 = lo;
    } else {
        u1 = (hi << r.s) | (lo >> (64 - r.s));
        u0 = lo << r.s;
    }

    /* Core Möller-Granlund step on normalised inputs:
     *   {q1, q0} = v * u1
     *   {q1, q0} += {u1, u0}
     *   q1 += 1
     *   r_n = u0 - q1 * d_norm
     *   if r_n > q0: q1 -= 1; r_n += d_norm
     *   if r_n >= d_norm: q1 += 1; r_n -= d_norm   (rare)
     */
    bigint_mul64x64(r.v, u1, &q0, &q1);
    {
        uint64_t carry;
        q0 = bigint_addc(q0, u0, 0, &carry);
        q1 = bigint_addc(q1, u1, carry, &carry);
        (void)carry;
    }
    q1 += 1;
    r_n = u0 - q1 * r.d_norm;
    if (r_n > q0) {
        q1 -= 1;
        r_n += r.d_norm;
    }
    if (r_n >= r.d_norm) {
        q1 += 1;
        r_n -= r.d_norm;
    }

    /* Remainder must be denormalised by shifting back. */
    *rem_out = r_n >> r.s;
    return q1;
}

/* ── benches ────────────────────────────────────────────────────────────── */

/* Sweep m limbs, divq-style (current bigint_div128 chain). */
static double bench_divq_chain(uint64_t d, int m, int N)
{
    int i, k;
    volatile uint64_t acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++) {
        uint64_t rem = 0;
        for (k = 0; k < m; k++) {
            uint64_t lo = (uint64_t)k * 0x9E3779B97F4A7C15ULL ^ (uint64_t)i;
            uint64_t q  = bigint_div128(rem, lo, d, &rem);
            acc ^= q;
        }
        acc ^= rem;
    }
    (void)acc;
    return (now() - t0) / ((double)N * (double)m);
}

/* Sweep m limbs, Möller-Granlund chain (one recip init, m gm_div128 calls).
 * The divisor is varied per outer iteration so the compiler can't hoist the
 * recip_init out of the loop, accurately modelling the real bigint_div_u64
 * call site where each call computes its own reciprocal. */
static double bench_gm_chain(uint64_t d_seed, int m, int N)
{
    int i, k;
    volatile uint64_t acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++) {
        /* Vary low bits of d so it's not loop-invariant; keep top bits
         * stable so the divisor magnitude class doesn't drift. */
        uint64_t d = d_seed ^ ((uint64_t)i & 0x3F);
        if (d == 0) d = d_seed;
        gm_recip_t r = gm_recip_init(d);
        uint64_t rem = 0;
        for (k = 0; k < m; k++) {
            uint64_t lo = (uint64_t)k * 0x9E3779B97F4A7C15ULL ^ (uint64_t)i;
            uint64_t q  = gm_div128(rem, lo, d, r, &rem);
            acc ^= q;
        }
        acc ^= rem;
    }
    (void)acc;
    return (now() - t0) / ((double)N * (double)m);
}

/* Same shape, but for the divq-only chain (also varying d to keep the
 * comparison fair under the per-iter divisor variation). */
static double bench_divq_chain_var(uint64_t d_seed, int m, int N)
{
    int i, k;
    volatile uint64_t acc = 0;
    double t0 = now();
    for (i = 0; i < N; i++) {
        uint64_t d = d_seed ^ ((uint64_t)i & 0x3F);
        if (d == 0) d = d_seed;
        uint64_t rem = 0;
        for (k = 0; k < m; k++) {
            uint64_t lo = (uint64_t)k * 0x9E3779B97F4A7C15ULL ^ (uint64_t)i;
            uint64_t q  = bigint_div128(rem, lo, d, &rem);
            acc ^= q;
        }
        acc ^= rem;
    }
    (void)acc;
    return (now() - t0) / ((double)N * (double)m);
}

/* Verify gm matches divq on every limb. */
static int verify(uint64_t d, int trials)
{
    int i, mismatches = 0;
    gm_recip_t r = gm_recip_init(d);
    for (i = 0; i < trials; i++) {
        uint64_t hi = ((uint64_t)i * 0x9E3779B97F4A7C15ULL) % d;
        uint64_t lo = (uint64_t)i * 0xBF58476D1CE4E5B9ULL;
        uint64_t r1 = 0, r2 = 0;
        uint64_t q1 = bigint_div128(hi, lo, d, &r1);
        uint64_t q2 = gm_div128(hi, lo, d, r, &r2);
        if (q1 != q2 || r1 != r2) mismatches++;
    }
    return mismatches;
}

int main(void)
{
    /* Realistic 3j Pass-2 divisors: small (s+1)^3 or up to (s+1)*(d+s+1)*(e+s+1)
     * for various j's.  Sweep tiny -> medium -> large.  m is the number of
     * limbs in the running scaled bigint; 1-3 covers small j (1-100), 4-30
     * covers j 100-1000, 250 covers j=4000. */
    struct { uint64_t d; const char *label; } cases[] = {
        { 7,                           "d=7 (tiny)             " },
        { 0x12345,                     "d=75077 (32-bit)       " },
        { 0xFFFFFFFFULL,               "d=2^32-1 (32-bit edge) " },
        { 0x100000007ULL,              "d=2^32+7 (33-bit)      " },
        { 0xDEADBEEFCAFEBABEULL,       "d=64-bit dense         " },
        { 0xFFFFFFFFFFFFFFFFULL,       "d=2^64-1 (saturated)   " },
    };
    int ms[] = { 1, 2, 3, 4, 8, 30, 250 };
    const int N_outer_small = 1000000;   /* for small m */
    const int N_outer_large = 30000;     /* for large m */
    size_t k;
    int mi;

    /* Quick correctness check first. */
    for (k = 0; k < sizeof(cases)/sizeof(cases[0]); k++) {
        int mm = verify(cases[k].d, 10000);
        if (mm) {
            fprintf(stderr, "VERIFY FAIL: d=%016lx mismatches=%d\n",
                    (unsigned long)cases[k].d, mm);
            return 1;
        }
    }
    printf("# verify OK on all divisors x 10000 trials\n");

#if defined(__x86_64__) && !defined(BIGINT_FORCE_PORTABLE) && !defined(BIGINT_NO_DIVQ)
    printf("# backend: x86-64 hardware divq\n");
#elif defined(BIGINT_FORCE_PORTABLE) || !defined(__SIZEOF_INT128__)
    printf("# backend: pure-C99 Algorithm D\n");
#else
    printf("# backend: __uint128_t Algorithm D\n");
#endif

    printf("\n%-26s %5s %12s %12s %8s\n",
           "case", "m", "divq ns/limb", "gm ns/limb", "speedup");
    printf("%-26s %5s %12s %12s %8s\n",
           "----", "-", "------------", "----------", "-------");
    for (k = 0; k < sizeof(cases)/sizeof(cases[0]); k++) {
        for (mi = 0; mi < (int)(sizeof(ms)/sizeof(ms[0])); mi++) {
            int m = ms[mi];
            int N = (m >= 30) ? N_outer_large
                  : (m >= 8)  ? (N_outer_small/4)
                              : N_outer_small;
            double td = 1e18, tg = 1e18;
            int trial;
            for (trial = 0; trial < 5; trial++) {
                double t1 = bench_divq_chain_var(cases[k].d, m, N);
                double t2 = bench_gm_chain      (cases[k].d, m, N);
                if (t1 < td) td = t1;
                if (t2 < tg) tg = t2;
            }
            printf("%-26s %5d %12.2f %12.2f %7.2fx\n",
                   cases[k].label, m, td*1e9, tg*1e9, td/tg);
        }
        printf("\n");
    }
    return 0;
}
