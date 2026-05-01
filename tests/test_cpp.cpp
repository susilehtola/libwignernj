// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// Tests for the C++11 header-only wrapper (wigner.hpp).

#include "../include/wigner.hpp"
#include <cmath>
#include <cstdio>
#include <stdexcept>

static int g_pass = 0;
static int g_fail = 0;

#define CHECK_NEAR(got, exp, tol) do {                                  \
    double _g = (double)(got), _e = (double)(exp), _t = (double)(tol); \
    double _d = std::abs(_g - _e);                                      \
    double _ref = std::abs(_e) > 1e-300 ? std::abs(_e) : 1.0;          \
    if (_d <= _t * _ref + 1e-300) {                                     \
        ++g_pass;                                                       \
    } else {                                                            \
        ++g_fail;                                                       \
        std::printf("FAIL %s:%d  got=%.15g  exp=%.15g  diff=%.3g\n",   \
                    __FILE__, __LINE__, _g, _e, _d);                    \
    }                                                                   \
} while(0)

#define CHECK_ABS(got, exp, tol) do {                                   \
    double _g = (double)(got), _e = (double)(exp), _t = (double)(tol); \
    double _d = std::abs(_g - _e);                                      \
    if (_d <= _t) {                                                     \
        ++g_pass;                                                       \
    } else {                                                            \
        ++g_fail;                                                       \
        std::printf("FAIL %s:%d  got=%.15g  exp=%.15g  diff=%.3g\n",   \
                    __FILE__, __LINE__, _g, _e, _d);                    \
    }                                                                   \
} while(0)

#define CHECK_THROWS(expr) do {                                         \
    bool _threw = false;                                                \
    try { (void)(expr); }                                               \
    catch (const std::invalid_argument &) { _threw = true; }           \
    if (_threw) { ++g_pass; }                                           \
    else { ++g_fail;                                                    \
        std::printf("FAIL %s:%d  expected std::invalid_argument\n",    \
                    __FILE__, __LINE__); }                              \
} while(0)

int main(void)
{
    using std::sqrt;

    /* ── symbol3j integer API ─────────────────────────────────────────── */
    CHECK_NEAR(wigner::symbol3j<float> (2,2,0, 0,0,0), -1.0f/sqrtf(3.0f), 2e-6f);
    CHECK_NEAR(wigner::symbol3j<double>(2,2,0, 0,0,0), -1.0/sqrt(3.0),    2e-15);
    CHECK_NEAR(wigner::symbol3j<long double>(2,2,0, 0,0,0), -1.0L/sqrtl(3.0L), 2e-18L);

    /* ── symbol3j real-valued overload ───────────────────────────────── */
    CHECK_NEAR(wigner::symbol3j(1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
               -1.0/sqrt(3.0), 2e-15);
    CHECK_NEAR(wigner::symbol3j<float>(0.5, 0.5, 0.0, 0.5, -0.5, 0.0),
               (float)(1.0/sqrt(2.0)), 2e-6f);

    /* ── symbol6j ─────────────────────────────────────────────────────── */
    CHECK_NEAR(wigner::symbol6j<float> (2,2,2, 2,2,2), 1.0f/6.0f, 2e-6f);
    CHECK_NEAR(wigner::symbol6j<double>(2,2,2, 2,2,2), 1.0/6.0,   2e-15);
    CHECK_NEAR(wigner::symbol6j<long double>(2,2,2, 2,2,2), 1.0L/6.0L, 2e-18L);

    CHECK_NEAR(wigner::symbol6j(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
               1.0/6.0, 2e-15);

    /* ── symbol9j ─────────────────────────────────────────────────────── */
    CHECK_NEAR(wigner::symbol9j<double>(1,1,2, 1,1,0, 2,2,2), sqrt(6.0)/18.0, 2e-15);
    CHECK_NEAR(wigner::symbol9j(0.5, 0.5, 1.0,
                                0.5, 0.5, 0.0,
                                1.0, 1.0, 1.0), sqrt(6.0)/18.0, 2e-15);

    /* ── cg ───────────────────────────────────────────────────────────── */
    CHECK_NEAR(wigner::cg<double>(1, 1, 1,-1, 2, 0), 1.0/sqrt(2.0), 2e-15);
    CHECK_NEAR(wigner::cg<float> (1, 1, 1,-1, 2, 0), (float)(1.0/sqrt(2.0)), 2e-6f);
    CHECK_NEAR(wigner::cg(0.5, 0.5, 0.5, -0.5, 1.0, 0.0), 1.0/sqrt(2.0), 2e-15);

    /* ── racahw ───────────────────────────────────────────────────────── */
    {
        /* racah_w(j1=1,j2=1,J=0,j3=1; j12=0,j23=1)
         * phase = (-1)^(j1+j2+J+j3) = (-1)^((2+2+0+2)/2) = (-1)^3 = -1 */
        double w  = wigner::racahw<double>(2,2,0,2,0,2);
        double w6 = wigner::symbol6j<double>(2,2,0, 2,0,2);
        CHECK_NEAR(w, -w6, 2e-15);
    }

    /* ── gaunt ────────────────────────────────────────────────────────── */
    {
        double g  = wigner::gaunt<double>(2,0,2,0,0,0);
        double w0 = wigner::symbol3j<double>(2,2,0, 0,0,0);
        double norm = sqrt(3.0 * 3.0 * 1.0 / (4.0 * 3.14159265358979323846));
        CHECK_NEAR(g, norm * w0 * w0, 2e-15);
    }

    /* ── selection-rule zeros ─────────────────────────────────────────── */
    CHECK_ABS(wigner::symbol3j<double>(2,2,2, 2,0,0), 0.0, 1e-15);
    CHECK_ABS(wigner::symbol6j<double>(0,0,2, 0,0,2), 0.0, 1e-15);
    CHECK_ABS(wigner::symbol9j<double>(0,0,2, 0,0,2, 0,0,2), 0.0, 1e-14);

    /* ── exception for non-half-integer real arg ──────────────────────── */
    CHECK_THROWS(wigner::symbol3j(0.3, 1.0, 1.0, 0.0, 0.0, 0.0));
    CHECK_THROWS(wigner::cg(0.5, 0.3, 0.5, 0.5, 1.0, 1.0));

    /* ── long double symbol9j ─────────────────────────────────────────── */
    CHECK_NEAR(wigner::symbol9j<long double>(1,1,2, 1,1,0, 2,2,2),
               sqrtl(6.0L)/18.0L, 2e-18L);

    /* ── large j ──────────────────────────────────────────────────────── */
    CHECK_NEAR(wigner::symbol3j<double>(100,100,0, 0,0,0),
               1.0/sqrt(101.0), 2e-15);

    int total = g_pass + g_fail;
    std::printf("%s: %d/%d tests passed\n",
                g_fail == 0 ? "PASS" : "FAIL", g_pass, total);
    return g_fail ? 1 : 0;
}
