// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// Tests for the C++11 header-only wrapper (wignernj.hpp).

#include "../include/wignernj.hpp"
#ifdef WIGNERNJ_HAVE_QUADMATH
#  include "../include/wignernj_quadmath.hpp"
#  include <quadmath.h>
#endif
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <vector>

// long double has 80-bit extended precision on x86-64 ELF, but is just an
// alias for double on Apple's arm64 / x86-64 ABIs (LDBL_MANT_DIG == 53) and
// on MSVC.  Scale the long-double tolerance accordingly so the test is not
// over-tight on platforms where long double is double in disguise.
static constexpr long double LD_TOL = 4 * LDBL_EPSILON;

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
    CHECK_NEAR(wignernj::symbol3j<float> (2,2,0, 0,0,0), -1.0f/sqrtf(3.0f), 2e-6f);
    CHECK_NEAR(wignernj::symbol3j<double>(2,2,0, 0,0,0), -1.0/sqrt(3.0),    2e-15);
    CHECK_NEAR(wignernj::symbol3j<long double>(2,2,0, 0,0,0), -1.0L/sqrtl(3.0L), LD_TOL);

    /* ── symbol3j real-valued overload ───────────────────────────────── */
    CHECK_NEAR(wignernj::symbol3j(1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
               -1.0/sqrt(3.0), 2e-15);
    CHECK_NEAR(wignernj::symbol3j<float>(0.5, 0.5, 0.0, 0.5, -0.5, 0.0),
               (float)(1.0/sqrt(2.0)), 2e-6f);

    /* ── symbol6j ─────────────────────────────────────────────────────── */
    CHECK_NEAR(wignernj::symbol6j<float> (2,2,2, 2,2,2), 1.0f/6.0f, 2e-6f);
    CHECK_NEAR(wignernj::symbol6j<double>(2,2,2, 2,2,2), 1.0/6.0,   2e-15);
    CHECK_NEAR(wignernj::symbol6j<long double>(2,2,2, 2,2,2), 1.0L/6.0L, LD_TOL);

    CHECK_NEAR(wignernj::symbol6j(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
               1.0/6.0, 2e-15);

    /* ── symbol9j ─────────────────────────────────────────────────────── */
    CHECK_NEAR(wignernj::symbol9j<double>(1,1,2, 1,1,0, 2,2,2), sqrt(6.0)/18.0, 2e-15);
    CHECK_NEAR(wignernj::symbol9j(0.5, 0.5, 1.0,
                                0.5, 0.5, 0.0,
                                1.0, 1.0, 1.0), sqrt(6.0)/18.0, 2e-15);

    /* ── cg ───────────────────────────────────────────────────────────── */
    CHECK_NEAR(wignernj::cg<double>(1, 1, 1,-1, 2, 0), 1.0/sqrt(2.0), 2e-15);
    CHECK_NEAR(wignernj::cg<float> (1, 1, 1,-1, 2, 0), (float)(1.0/sqrt(2.0)), 2e-6f);
    CHECK_NEAR(wignernj::cg(0.5, 0.5, 0.5, -0.5, 1.0, 0.0), 1.0/sqrt(2.0), 2e-15);

    /* ── racahw ───────────────────────────────────────────────────────── */
    {
        /* racah_w(j1=1,j2=1,J=0,j3=1; j12=0,j23=1)
         * phase = (-1)^(j1+j2+J+j3) = (-1)^((2+2+0+2)/2) = (-1)^3 = -1 */
        double w  = wignernj::racahw<double>(2,2,0,2,0,2);
        double w6 = wignernj::symbol6j<double>(2,2,0, 2,0,2);
        CHECK_NEAR(w, -w6, 2e-15);
    }

    /* ── gaunt ────────────────────────────────────────────────────────── */
    {
        double g  = wignernj::gaunt<double>(2,0,2,0,0,0);
        double w0 = wignernj::symbol3j<double>(2,2,0, 0,0,0);
        double norm = sqrt(3.0 * 3.0 * 1.0 / (4.0 * 3.14159265358979323846));
        CHECK_NEAR(g, norm * w0 * w0, 2e-15);
    }

    /* ── gauntreal ────────────────────────────────────────────────────── */
    {
        /* Real Gaunt at all-m=0 equals complex Gaunt. */
        double gr = wignernj::gauntreal<double>(2,0, 2,0, 4,0);
        double gc = wignernj::gaunt<double>     (2,0, 2,0, 4,0);
        CHECK_NEAR(gr, gc, 2e-15);
        CHECK_NEAR(wignernj::gauntreal<float>(2,0,2,0,4,0),
                   (float)gc, 2e-6f);
        CHECK_NEAR(wignernj::gauntreal<long double>(2,0,2,0,4,0),
                   (long double)gc, LD_TOL);
        /* Real-valued overload */
        CHECK_NEAR(wignernj::gauntreal(1.0, 0.0, 1.0, 0.0, 2.0, 0.0),
                   gc, 2e-15);
    }

    /* ── fanox ────────────────────────────────────────────────────────── */
    {
        /* X = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * 9j.
         * At all-equal-tj=4 the four factors are 5 each, so norm=25. */
        double w9 = wignernj::symbol9j<double>(4,4,4, 4,4,4, 4,4,4);
        CHECK_NEAR(wignernj::fanox<double>(4,4,4, 4,4,4, 4,4,4),
                   25.0 * w9, 2e-13);
        CHECK_NEAR(wignernj::fanox<float>(4,4,4, 4,4,4, 4,4,4),
                   (float)(25.0 * w9), 2e-5f);
        long double w9_ld = wignernj::symbol9j<long double>(4,4,4, 4,4,4, 4,4,4);
        CHECK_NEAR(wignernj::fanox<long double>(4,4,4, 4,4,4, 4,4,4),
                   25.0L * w9_ld, LD_TOL * 50);
        /* Real-valued overload */
        CHECK_NEAR(wignernj::fanox(2.0, 2.0, 2.0,
                                 2.0, 2.0, 2.0,
                                 2.0, 2.0, 2.0),
                   25.0 * w9, 2e-13);
        /* Selection-rule zero */
        CHECK_ABS(wignernj::fanox<double>(2,2,2, 2,2,2, 2,2,5), 0.0, 1e-15);
    }

    /* ── selection-rule zeros ─────────────────────────────────────────── */
    CHECK_ABS(wignernj::symbol3j<double>(2,2,2, 2,0,0), 0.0, 1e-15);
    CHECK_ABS(wignernj::symbol6j<double>(0,0,2, 0,0,2), 0.0, 1e-15);
    CHECK_ABS(wignernj::symbol9j<double>(0,0,2, 0,0,2, 0,0,2), 0.0, 1e-14);

    /* ── exception for non-half-integer real arg ──────────────────────── */
    CHECK_THROWS(wignernj::symbol3j(0.3, 1.0, 1.0, 0.0, 0.0, 0.0));
    CHECK_THROWS(wignernj::cg(0.5, 0.3, 0.5, 0.5, 1.0, 1.0));

    /* ── long double symbol9j ─────────────────────────────────────────── */
    CHECK_NEAR(wignernj::symbol9j<long double>(1,1,2, 1,1,0, 2,2,2),
               sqrtl(6.0L)/18.0L, LD_TOL);

    /* ── large j ──────────────────────────────────────────────────────── */
    CHECK_NEAR(wignernj::symbol3j<double>(100,100,0, 0,0,0),
               1.0/sqrt(101.0), 2e-15);

    /* ── real_ylm_in_complex_ylm ──────────────────────────────────────── */
    {
        // std::vector<std::complex<T>> return form, l = 1.
        auto C = wignernj::real_ylm_in_complex_ylm<double>(1);
        const double s = 1.0 / std::sqrt(2.0);
        // Column-major (m_r=+1, m_c=-1) at flat index 0*3 + 2 = 2.
        CHECK_NEAR(C[2].real(), s,   2e-15);
        CHECK_ABS (C[2].imag(), 0.0, 1e-300);
        // Column-major (m_r=-1, m_c=-1) at flat index 0*3 + 0 = 0.
        CHECK_ABS (C[0].real(), 0.0, 1e-300);
        CHECK_NEAR(C[0].imag(), s,   2e-15);
        // C[0, 0] at flat index 1*3 + 1 = 4.
        CHECK_NEAR(C[4].real(), 1.0, 1e-15);
        // float overload, in-place fill.
        std::complex<float> Cf[9];
        wignernj::real_ylm_in_complex_ylm<float>(1, Cf);
        CHECK_NEAR(Cf[2].real(), (float)s, 5e-7f);
        // long double overload.
        std::complex<long double> Cl[9];
        wignernj::real_ylm_in_complex_ylm<long double>(1, Cl);
        CHECK_NEAR((double)Cl[2].real(), s, 2.0 * DBL_EPSILON);
        // l = 0 edge case: 1x1 matrix == [[1+0i]].
        auto C0 = wignernj::real_ylm_in_complex_ylm<double>(0);
        CHECK_NEAR(C0[0].real(), 1.0, 1e-15);
        CHECK_ABS (C0[0].imag(), 0.0, 1e-300);
    }

#ifdef WIGNERNJ_HAVE_QUADMATH
    /* ── __float128 specialisations from wignernj_quadmath.hpp ──────────── */
    {
        /* Quad ulp tolerance: FLT128_EPSILON = 2^-112 ≈ 1.93e-34 */
        const __float128 q_tol = 4.0Q * FLT128_EPSILON;

        /* symbol3j<__float128>: closed form -1/sqrt(3). */
        {
            __float128 v   = wignernj::symbol3j<__float128>(2,2,0, 0,0,0);
            __float128 ref = -1.0Q / sqrtq(3.0Q);
            __float128 d   = fabsq(v - ref);
            if (d <= q_tol * fabsq(ref)) ++g_pass;
            else { ++g_fail; std::printf("FAIL %s:%d  quad symbol3j\n",
                                        __FILE__, __LINE__); }
        }

        /* Real-valued (double-arg) overload at __float128: agrees with
         * the integer-API specialisation to the last quad ulp. */
        {
            __float128 v1 = wignernj::symbol3j<__float128>(0.5, 0.5, 0.0,
                                                          0.5,-0.5, 0.0);
            __float128 v2 = wignernj::symbol3j<__float128>(1, 1, 0, 1,-1, 0);
            if (v1 == v2) ++g_pass;
            else { ++g_fail; std::printf("FAIL %s:%d  quad real-overload\n",
                                        __FILE__, __LINE__); }
        }

        /* symbol6j, symbol9j, cg, racahw, fanox, gaunt, gauntreal:
         * agree with long double specialisations within long-double ulps,
         * since quad has strictly more bits than long double. */
        #define CHECK_QUAD_AGAINST_LD(call_q, call_ld) do {                \
            __float128 q   = (call_q);                                    \
            __float128 ld  = (__float128)(call_ld);                       \
            __float128 d   = fabsq(q - ld);                               \
            __float128 sc  = fabsq(ld) > 1e-300Q ? fabsq(ld) : 1.0Q;      \
            __float128 tol = 4.0Q * (__float128)LDBL_EPSILON;             \
            if (d <= tol * sc) ++g_pass;                                  \
            else { ++g_fail;                                              \
                std::printf("FAIL %s:%d  %s\n",                           \
                            __FILE__, __LINE__, #call_q); }               \
        } while(0)

        CHECK_QUAD_AGAINST_LD(
            wignernj::symbol6j<__float128>(2,2,2, 2,2,2),
            wignernj::symbol6j<long double>(2,2,2, 2,2,2));
        CHECK_QUAD_AGAINST_LD(
            wignernj::symbol9j<__float128>(1,1,2, 1,1,0, 2,2,2),
            wignernj::symbol9j<long double>(1,1,2, 1,1,0, 2,2,2));
        CHECK_QUAD_AGAINST_LD(
            wignernj::cg<__float128>(1, 1, 1,-1, 2, 0),
            wignernj::cg<long double>(1, 1, 1,-1, 2, 0));
        CHECK_QUAD_AGAINST_LD(
            wignernj::racahw<__float128>(2,2,0, 2,0,2),
            wignernj::racahw<long double>(2,2,0, 2,0,2));
        CHECK_QUAD_AGAINST_LD(
            wignernj::fanox<__float128>(4,4,4, 4,4,4, 4,4,4),
            wignernj::fanox<long double>(4,4,4, 4,4,4, 4,4,4));
        CHECK_QUAD_AGAINST_LD(
            wignernj::gaunt<__float128>(2,0, 2,0, 0,0),
            wignernj::gaunt<long double>(2,0, 2,0, 0,0));
        CHECK_QUAD_AGAINST_LD(
            wignernj::gauntreal<__float128>(2,0, 2,0, 4,0),
            wignernj::gauntreal<long double>(2,0, 2,0, 4,0));
        #undef CHECK_QUAD_AGAINST_LD

        /* real_ylm_in_complex_ylm non-template overload at l=1: same
         * inv_sqrt2 entries as the C and Fortran tests. */
        {
            wignernj_cfloat128_t C[9];
            wignernj::real_ylm_in_complex_ylm(1, C);
            const __float128 *F = reinterpret_cast<const __float128 *>(C);
            const __float128 inv_sqrt2 = 1.0Q / sqrtq(2.0Q);
            __float128 d1 = fabsq(F[4] - inv_sqrt2);  /* re part, idx 4 */
            __float128 d2 = fabsq(F[1] - inv_sqrt2);  /* im part, idx 1 */
            __float128 d3 = fabsq(F[8] - 1.0Q);
            if (d1 <= q_tol && d2 <= q_tol && d3 <= q_tol) ++g_pass;
            else { ++g_fail; std::printf("FAIL %s:%d  quad real_ylm fill\n",
                                        __FILE__, __LINE__); }
        }

        /* Vector-returning real_ylm_in_complex_ylm_q form. */
        {
            std::vector<wignernj_cfloat128_t> C =
                wignernj::real_ylm_in_complex_ylm_q(1);
            if (C.size() == 9) ++g_pass;
            else { ++g_fail; std::printf("FAIL %s:%d  quad vector size\n",
                                        __FILE__, __LINE__); }
        }
    }
#endif

    int total = g_pass + g_fail;
    std::printf("%s: %d/%d tests passed\n",
                g_fail == 0 ? "PASS" : "FAIL", g_pass, total);
    return g_fail ? 1 : 0;
}
