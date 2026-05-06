/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Symmetry-oracle tests for 3j, 6j, 9j symbols.
 *
 * The reference-table tests in test_3j/6j/9j.c pin specific values against
 * SymPy.  These tests are complementary: they sweep many randomized
 * (deterministic-seed) configurations and check that the well-known
 * permutation/phase symmetries of each symbol hold.  Sign-bit errors,
 * factor-of-(-1)^J phase mistakes, and Racah-sum-bound asymmetries that
 * would slip past a hand-picked table are caught here cheaply.
 *
 * All arguments are 2j integers.  Triangle and m-parity constraints are
 * built into the generators so every accepted draw is a valid input.
 */
#include "run_tests.h"
#include "../include/wignernj.h"
#include <stdint.h>
#include <stdlib.h>

/* ── deterministic PRNG ─────────────────────────────────────────────────── */
static uint64_t g_state = 0x9E3779B97F4A7C15ULL;
static int rand_int(int lo, int hi)
{
    g_state = g_state * 6364136223846793005ULL + 1442695040888963407ULL;
    uint32_t r = (uint32_t)(g_state >> 33);
    return lo + (int)(r % (uint32_t)(hi - lo + 1));
}

/* ── 3j symmetries ──────────────────────────────────────────────────────── */
/*
 *  Cyclic columns:        (j1 j2 j3; m1 m2 m3) =  (j2 j3 j1; m2 m3 m1)
 *  Odd column swap:       (j1 j2 j3; m1 m2 m3) =  (-1)^J (j2 j1 j3; m2 m1 m3)
 *  Sign-flip of all m:    (j1 j2 j3; m1 m2 m3) =  (-1)^J (j1 j2 j3;-m1 -m2 -m3)
 *  where J = j1+j2+j3 (always integer because triangle ⇒ tj1+tj2+tj3 even).
 */
static void test_3j_symmetries(void)
{
    const int TRIALS = 400;
    const int MAX_TJ = 16;
    int valid = 0;
    for (int t = 0; t < TRIALS; t++) {
        int tj1 = rand_int(0, MAX_TJ);
        int tj2 = rand_int(0, MAX_TJ);
        int tj3 = rand_int(0, MAX_TJ);
        if ((tj1 + tj2 + tj3) & 1) continue;
        if (tj3 < abs(tj1 - tj2) || tj3 > tj1 + tj2) continue;
        /* tm has the same parity as tj iff tm = -tj + 2k for k ∈ [0, tj]. */
        int tm1 = -tj1 + 2 * rand_int(0, tj1);
        int tm2 = -tj2 + 2 * rand_int(0, tj2);
        int tm3 = -tm1 - tm2;
        if (abs(tm3) > tj3) continue;
        if ((tj3 - tm3) & 1) continue; /* should follow from sums; cheap check */

        valid++;
        int J = (tj1 + tj2 + tj3) / 2;
        double phase = (J & 1) ? -1.0 : 1.0;
        double v = wigner3j(tj1, tj2, tj3, tm1, tm2, tm3);

        TEST_DBL(v, wigner3j(tj2, tj3, tj1, tm2, tm3, tm1));
        TEST_DBL(v, wigner3j(tj3, tj1, tj2, tm3, tm1, tm2));
        TEST_DBL(v, phase * wigner3j(tj2, tj1, tj3, tm2, tm1, tm3));
        TEST_DBL(v, phase * wigner3j(tj1, tj2, tj3, -tm1, -tm2, -tm3));
    }
    printf("  3j: %d valid configurations exercised\n", valid);
    TEST_ASSERT(valid > 50);
}

/* ── 6j symmetries ──────────────────────────────────────────────────────── */
/*
 *  Column permutations (any S_3): {j1 j2 j3; j4 j5 j6} invariant.
 *  Upper-lower swap of any two columns: e.g.
 *       {j1 j2 j3; j4 j5 j6} = {j4 j5 j3; j1 j2 j6}.
 */
static int triangle(int ta, int tb, int tc)
{
    if ((ta + tb + tc) & 1) return 0;
    return tc >= abs(ta - tb) && tc <= ta + tb;
}

static void test_6j_symmetries(void)
{
    /* Independent random draws hit all four triangles rarely; bump trials. */
    const int TRIALS = 20000;
    const int MAX_TJ = 8;
    int valid = 0;
    for (int t = 0; t < TRIALS; t++) {
        int tj1 = rand_int(0, MAX_TJ), tj2 = rand_int(0, MAX_TJ);
        int tj3 = rand_int(0, MAX_TJ), tj4 = rand_int(0, MAX_TJ);
        int tj5 = rand_int(0, MAX_TJ), tj6 = rand_int(0, MAX_TJ);
        if (!triangle(tj1, tj2, tj3)) continue;
        if (!triangle(tj1, tj5, tj6)) continue;
        if (!triangle(tj4, tj2, tj6)) continue;
        if (!triangle(tj4, tj5, tj3)) continue;

        valid++;
        double v = wigner6j(tj1, tj2, tj3, tj4, tj5, tj6);

        /* Column swaps (1↔2, 1↔3, 2↔3) — the 6j is invariant. */
        TEST_DBL(v, wigner6j(tj2, tj1, tj3, tj5, tj4, tj6));
        TEST_DBL(v, wigner6j(tj3, tj2, tj1, tj6, tj5, tj4));
        TEST_DBL(v, wigner6j(tj1, tj3, tj2, tj4, tj6, tj5));
        /* Cyclic columns. */
        TEST_DBL(v, wigner6j(tj2, tj3, tj1, tj5, tj6, tj4));
        /* Upper-lower swap of two columns: cols (1,2). */
        TEST_DBL(v, wigner6j(tj4, tj5, tj3, tj1, tj2, tj6));
        /* Upper-lower swap of cols (2,3). */
        TEST_DBL(v, wigner6j(tj1, tj5, tj6, tj4, tj2, tj3));
        /* Upper-lower swap of cols (1,3). */
        TEST_DBL(v, wigner6j(tj4, tj2, tj6, tj1, tj5, tj3));
    }
    printf("  6j: %d valid configurations exercised\n", valid);
    TEST_ASSERT(valid > 50);
}

/* ── 9j symmetries ──────────────────────────────────────────────────────── */
/*
 *  Transpose:                  9j[ij] = 9j[ji]
 *  Odd row or column swap:     factor (-1)^S, S = sum of all nine j-values.
 *  (S = sum(tj)/2; integer because every row and every column is a triangle.)
 */
static void test_9j_symmetries(void)
{
    /* Six triangles to satisfy; brute-force needs many trials. */
    const int TRIALS = 200000;
    const int MAX_TJ = 4;
    int valid = 0;
    for (int t = 0; t < TRIALS; t++) {
        int tj[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                tj[i][j] = rand_int(0, MAX_TJ);
        /* Six triangles: three rows, three columns. */
        int ok = 1;
        for (int i = 0; i < 3 && ok; i++)
            ok = triangle(tj[i][0], tj[i][1], tj[i][2]);
        for (int j = 0; j < 3 && ok; j++)
            ok = triangle(tj[0][j], tj[1][j], tj[2][j]);
        if (!ok) continue;

        valid++;
        int S = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) S += tj[i][j];
        S /= 2;
        double phase = (S & 1) ? -1.0 : 1.0;

        double v = wigner9j(tj[0][0], tj[0][1], tj[0][2],
                             tj[1][0], tj[1][1], tj[1][2],
                             tj[2][0], tj[2][1], tj[2][2]);

        /* Transpose. */
        TEST_DBL(v, wigner9j(tj[0][0], tj[1][0], tj[2][0],
                              tj[0][1], tj[1][1], tj[2][1],
                              tj[0][2], tj[1][2], tj[2][2]));
        /* Swap rows 1↔2: factor (-1)^S. */
        TEST_DBL(v, phase * wigner9j(tj[1][0], tj[1][1], tj[1][2],
                                       tj[0][0], tj[0][1], tj[0][2],
                                       tj[2][0], tj[2][1], tj[2][2]));
        /* Swap cols 1↔2: factor (-1)^S. */
        TEST_DBL(v, phase * wigner9j(tj[0][1], tj[0][0], tj[0][2],
                                       tj[1][1], tj[1][0], tj[1][2],
                                       tj[2][1], tj[2][0], tj[2][2]));
        /* Swap rows 1↔3: factor (-1)^S. */
        TEST_DBL(v, phase * wigner9j(tj[2][0], tj[2][1], tj[2][2],
                                       tj[1][0], tj[1][1], tj[1][2],
                                       tj[0][0], tj[0][1], tj[0][2]));
        /* Cyclic row permutation (even, no phase). */
        TEST_DBL(v, wigner9j(tj[1][0], tj[1][1], tj[1][2],
                              tj[2][0], tj[2][1], tj[2][2],
                              tj[0][0], tj[0][1], tj[0][2]));
    }
    printf("  9j: %d valid configurations exercised\n", valid);
    TEST_ASSERT(valid > 50);
}

int main(void)
{
    test_3j_symmetries();
    test_6j_symmetries();
    test_9j_symmetries();
    SUMMARY();
}
