# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# pytest suite for the Python wignernj extension.
# Run from the build directory after: pip install -e . --no-build-isolation
# or after cmake --build build with PYTHONPATH set.

import math
import pytest

try:
    from wignernj import (wigner3j, wigner6j, wigner9j,
                        clebsch_gordan, racah_w, fano_x,
                        gaunt, gaunt_real, real_ylm_in_complex_ylm)
except ImportError:
    pytest.skip("wignernj extension not installed", allow_module_level=True)

# ── helpers ──────────────────────────────────────────────────────────────────

def near(a, b, rtol=2e-13):
    ref = abs(b) if abs(b) > 1e-300 else 1.0
    return abs(a - b) <= rtol * ref

def near_abs(a, b, atol=1e-15):
    return abs(a - b) <= atol

# ── wigner3j ─────────────────────────────────────────────────────────────────

class TestWigner3j:
    def test_m_conservation_zero(self):
        assert near_abs(wigner3j(1, 1, 1, 1, 0, 0), 0.0)

    def test_triangle_violated(self):
        assert near_abs(wigner3j(0, 0, 1, 0, 0, 0), 0.0)

    def test_parity_zero(self):
        # j1+j2+j3=3 odd, m=0 → 0
        assert near_abs(wigner3j(1, 1, 1, 0, 0, 0), 0.0)

    def test_abs_m_exceeds_j(self):
        assert near_abs(wigner3j(1, 1, 1, 2, 0, -2), 0.0)

    def test_1_1_0(self):
        assert near(wigner3j(1, 1, 0, 0, 0, 0), -1.0/math.sqrt(3.0))

    def test_2_2_0(self):
        assert near(wigner3j(2, 2, 0, 0, 0, 0), 1.0/math.sqrt(5.0))

    def test_half_integer_inputs(self):
        # (1/2 1/2 1; 1/2 -1/2 0) = 1/sqrt(6)
        assert near(wigner3j(0.5, 0.5, 1, 0.5, -0.5, 0), 1.0/math.sqrt(6.0))

    def test_half_integer_float(self):
        assert near(wigner3j(1.0, 1.0, 0.0, 0.0, 0.0, 0.0), -1.0/math.sqrt(3.0))

    def test_fraction_input(self):
        from fractions import Fraction
        assert near(wigner3j(Fraction(1,2), Fraction(1,2), 1,
                             Fraction(1,2), Fraction(-1,2), 0),
                    1.0/math.sqrt(6.0))

    def test_large_j(self):
        # (50 50 0; 0 0 0) = 1/sqrt(101)
        assert near(wigner3j(50, 50, 0, 0, 0, 0), 1.0/math.sqrt(101.0))

    def test_integer_inputs(self):
        # (2 2 4; 0 0 0) = sqrt(2/35)
        assert near(wigner3j(2, 2, 4, 0, 0, 0), math.sqrt(2.0/35.0))

    def test_precision_float(self):
        v = wigner3j(1, 1, 0, 0, 0, 0, precision='float')
        assert near(v, -1.0/math.sqrt(3.0), rtol=1e-6)

    def test_precision_long_double(self):
        v = wigner3j(1, 1, 0, 0, 0, 0, precision='longdouble')
        assert near(v, -1.0/math.sqrt(3.0))

    def test_invalid_half_integer(self):
        with pytest.raises((ValueError, TypeError)):
            wigner3j(0.3, 1, 1, 0, 0, 0)

    def test_invalid_precision(self):
        with pytest.raises((ValueError, TypeError)):
            wigner3j(1, 1, 0, 0, 0, 0, precision='quad')

    def test_precision_non_string(self):
        with pytest.raises((ValueError, TypeError)):
            wigner3j(1, 1, 0, 0, 0, 0, precision=42)

    def test_regge_symmetry(self):
        v1 = wigner3j(2, 1, 2, 1, -1, 0)
        v2 = wigner3j(1, 2, 2, -1, 1, 0)
        # (-1)^(2+1+2) = -1
        assert near(v1, -v2)

# ── wigner6j ─────────────────────────────────────────────────────────────────

class TestWigner6j:
    def test_triangle_zero(self):
        assert near_abs(wigner6j(0, 0, 1, 0, 0, 1), 0.0)

    def test_1_1_1(self):
        assert near(wigner6j(1, 1, 1, 1, 1, 1), 1.0/6.0)

    def test_half_j(self):
        # {1/2 1/2 1; 1/2 1/2 0} = 1/2
        assert near(wigner6j(0.5, 0.5, 1.0, 0.5, 0.5, 0.0), 0.5)

    def test_closed_form(self):
        # {j j 0; j j j} = (-1)^(3j)/(2j+1), integer j
        for j in (1, 2, 5):
            val = wigner6j(j, j, 0, j, j, j)
            expected = ((-1)**(3*j)) / (2*j + 1)
            assert near(val, expected), f"j={j}: {val} != {expected}"

    def test_known_value(self):
        # {1 2 3; 3 2 1} = sqrt(2/175)
        assert near(wigner6j(1, 2, 3, 3, 2, 1), math.sqrt(2.0/175.0))

    def test_precision_float(self):
        v = wigner6j(1, 1, 1, 1, 1, 1, precision='float')
        assert near(v, 1.0/6.0, rtol=1e-6)

# ── wigner9j ─────────────────────────────────────────────────────────────────

class TestWigner9j:
    def test_triangle_zero(self):
        assert near_abs(wigner9j(0, 0, 1, 0, 0, 1, 0, 0, 1), 0.0, atol=1e-14)

    def test_known_value(self):
        # {1/2 1/2 1; 1/2 1/2 0; 1 1 1} = sqrt(6)/18
        assert near(wigner9j(0.5, 0.5, 1, 0.5, 0.5, 0, 1, 1, 1), math.sqrt(6.0)/18.0)

    def test_reduction_to_6j(self):
        # {1 1 0; 1 1 0; 0 0 0} = 6j{1,1,0;1,1,0}/sqrt(1)
        w9 = wigner9j(1, 1, 0, 1, 1, 0, 0, 0, 0)
        w6 = wigner6j(1, 1, 0, 1, 1, 0)
        assert near(w9, w6)

# ── Clebsch-Gordan ────────────────────────────────────────────────────────────

class TestClebschGordan:
    def test_1_2(self):
        # <1/2 1/2; 1/2 -1/2 | 1 0> = 1/sqrt(2)
        assert near(clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0), 1.0/math.sqrt(2.0))

    def test_singlet(self):
        # <1/2 -1/2; 1/2 1/2 | 0 0> = -1/sqrt(2)
        assert near(clebsch_gordan(0.5, -0.5, 0.5, 0.5, 0, 0), -1.0/math.sqrt(2.0))

    def test_triplet_max(self):
        # <1/2 1/2; 1/2 1/2 | 1 1> = 1
        assert near(clebsch_gordan(0.5, 0.5, 0.5, 0.5, 1, 1), 1.0)

    def test_m_conservation(self):
        assert near_abs(clebsch_gordan(1, 1, 1, 0, 2, 0), 0.0)

    def test_completeness(self):
        c1 = clebsch_gordan(0.5,  0.5, 0.5, -0.5, 1, 0)
        c2 = clebsch_gordan(0.5, -0.5, 0.5,  0.5, 1, 0)
        assert near(c1**2 + c2**2, 1.0)

    def test_integer_j(self):
        # <1 0; 1 0 | 2 0> = sqrt(2/3)
        assert near(clebsch_gordan(1, 0, 1, 0, 2, 0), math.sqrt(2.0/3.0))

    def test_precision_float(self):
        v = clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0, precision='float')
        assert near(v, 1.0/math.sqrt(2.0), rtol=1e-6)

# ── Racah W ───────────────────────────────────────────────────────────────────

class TestRacahW:
    def test_vs_6j(self):
        # W(j1,j2,J,j3;j12,j23) = (-1)^(j1+j2+J+j3) * 6j{j1,j2,j12;j3,J,j23}
        j1, j2, J, j3, j12, j23 = 1, 1, 0, 1, 0, 1
        w  = racah_w(j1, j2, J, j3, j12, j23)
        w6 = wigner6j(j1, j2, j12, j3, J, j23)
        phase = (-1)**int(j1 + j2 + J + j3)
        assert near(w, phase * w6)

    def test_vs_6j_2(self):
        j1, j2, J, j3, j12, j23 = 1, 1, 1, 1, 1, 1
        w  = racah_w(j1, j2, J, j3, j12, j23)
        w6 = wigner6j(j1, j2, j12, j3, J, j23)
        phase = (-1)**int(j1 + j2 + J + j3)
        assert near(w, phase * w6)

# ── Fano X ────────────────────────────────────────────────────────────────────

class TestFanoX:
    def test_vs_9j(self):
        # X(j1,j2,j12;j3,j4,j34;j13,j24,J)
        #   = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * {9j}.
        # all-equal-j=2 → norm = 5^2 = 25
        j = 2
        x = fano_x(j, j, j, j, j, j, j, j, j)
        w = wigner9j(j, j, j, j, j, j, j, j, j)
        assert near(x, 25.0 * w)

    def test_vs_9j_half_integer(self):
        # all-equal-j=1/2 → norm = 2^2 = 4
        x = fano_x(0.5, 0.5, 1.0, 0.5, 0.5, 0.0, 1.0, 1.0, 1.0)
        w = wigner9j(0.5, 0.5, 1.0, 0.5, 0.5, 0.0, 1.0, 1.0, 1.0)
        norm = math.sqrt(3.0 * 1.0 * 3.0 * 3.0)  # (2*1+1)(2*0+1)(2*1+1)(2*1+1)
        assert near(x, norm * w)

    def test_triangle_zero(self):
        # one triangle violated → X = 0
        assert near_abs(fano_x(1, 1, 1, 1, 1, 1, 1, 1, 3), 0.0)

    def test_precision_keyword(self):
        v = fano_x(2, 2, 2, 2, 2, 2, 2, 2, 2, precision='float')
        assert isinstance(v, float)


# ── Gaunt ─────────────────────────────────────────────────────────────────────

class TestGaunt:
    def test_basic(self):
        # G(1,0,1,0,0,0) = norm * 3j^2
        g   = gaunt(1, 0, 1, 0, 0, 0)
        w0  = wigner3j(1, 1, 0, 0, 0, 0)
        norm = math.sqrt(3.0 * 3.0 * 1.0 / (4.0 * math.pi))
        assert near(g, norm * w0 * w0)

    def test_m_sum_zero(self):
        assert near_abs(gaunt(1, 1, 1, 0, 1, 0), 0.0)

    def test_odd_parity_zero(self):
        # l1+l2+l3=3 odd: G=0
        assert near_abs(gaunt(1, 0, 1, 0, 1, 0), 0.0)


# ── Real-spherical-harmonic Gaunt ─────────────────────────────────────────────

class TestGauntReal:
    def test_S00_cubed(self):
        # ∫ (1/(2 sqrt π))^3 dΩ = 1 / (2 sqrt π)
        assert near(gaunt_real(0, 0, 0, 0, 0, 0), 1.0 / (2 * math.sqrt(math.pi)))

    def test_l_parity_zero(self):
        # l1+l2+l3 odd → 0
        assert near_abs(gaunt_real(1, 0, 1, 0, 1, 0), 0.0)

    def test_S10sq_S20(self):
        # ∫ S(1,0)^2 S(2,0) dΩ = 1 / sqrt(5 π)
        assert near(gaunt_real(1, 0, 1, 0, 2, 0), 1.0 / math.sqrt(5 * math.pi))

    def test_S1p1sq_S00(self):
        # ∫ S(1,+1)^2 S(0,0) dΩ = 1 / (2 sqrt π)
        assert near(gaunt_real(1, 1, 1, 1, 0, 0), 1.0 / (2 * math.sqrt(math.pi)))

    def test_S1m1sq_S00(self):
        # symmetry m → -m at n_- = 2
        assert near(gaunt_real(1, -1, 1, -1, 0, 0),
                    1.0 / (2 * math.sqrt(math.pi)))

    def test_cos_sin_vanishes(self):
        # ∫ cos(φ) sin(φ) dφ = 0  ⇒  cos·sin·1 vanishes
        assert near_abs(gaunt_real(1, 1, 1, -1, 0, 0), 0.0)

    def test_n_minus_odd(self):
        # one m_i < 0, two > 0: imaginary contribution, real Gaunt = 0
        assert near_abs(gaunt_real(1, -1, 1, 1, 1, 1), 0.0)

    def test_permutation_symmetry(self):
        a = gaunt_real(1, 0, 2, 1, 2, -1)
        b = gaunt_real(2, 1, 1, 0, 2, -1)
        c = gaunt_real(2, -1, 2, 1, 1, 0)
        assert near(a, b)
        assert near(a, c)

    def test_m0_matches_complex_gaunt(self):
        # When all m_i = 0, S_{l,0} = Y_l^0 so the real and complex Gaunts agree.
        for l1, l2, l3 in [(1, 1, 0), (2, 2, 0), (2, 2, 4),
                           (3, 3, 0), (3, 3, 4), (4, 4, 8)]:
            assert near(gaunt_real(l1, 0, l2, 0, l3, 0),
                        gaunt(l1, 0, l2, 0, l3, 0))

    def test_higher_l(self):
        # Spot-check at moderate l (5, 5, 6) using sympy convention.
        # Computed in sympy via the unitary transform: ~3.7e-2
        v = gaunt_real(5, 3, 5, -3, 6, 0)
        # Cross-check against permuted argument list (must match exactly):
        assert near(v, gaunt_real(5, -3, 5, 3, 6, 0))
        assert near(v, gaunt_real(6, 0, 5, 3, 5, -3))

    def test_precision_keyword(self):
        ref = 1.0 / math.sqrt(5 * math.pi)
        assert near(gaunt_real(1, 0, 1, 0, 2, 0,
                               precision='float'),  ref, rtol=1e-6)
        assert near(gaunt_real(1, 0, 1, 0, 2, 0,
                               precision='double'), ref)
        assert near(gaunt_real(1, 0, 1, 0, 2, 0,
                               precision='longdouble'), ref)


# ── real_ylm_in_complex_ylm ──────────────────────────────────────────────────

class TestRealYlmInComplexYlm:
    def test_l0_identity(self):
        C = real_ylm_in_complex_ylm(0)
        assert C == [[complex(1.0, 0.0)]]

    def test_l1_entries(self):
        # S_{1, 0}  = Y_1^0     -> C[1,1] = 1
        # S_{1,+1} = (1/sqrt2)( Y_1^{-1} - Y_1^{+1} )
        # S_{1,-1} = (i/sqrt2)( Y_1^{-1} + Y_1^{+1} )
        C = real_ylm_in_complex_ylm(1)
        s = 1.0 / math.sqrt(2.0)
        # Index is C[m_r + l][m_c + l].
        assert near_abs(C[1][1].real, 1.0)
        assert near_abs(C[1][1].imag, 0.0)
        assert near_abs(C[2][0].real, +s)         # (m_r=+1, m_c=-1)
        assert near_abs(C[2][2].real, -s)         # (m_r=+1, m_c=+1)
        assert near_abs(C[0][0].imag, +s)         # (m_r=-1, m_c=-1)
        assert near_abs(C[0][2].imag, +s)         # (m_r=-1, m_c=+1)

    def test_unitarity(self):
        # C C^H = I for l = 0..5.
        for l in range(6):
            n = 2 * l + 1
            C = real_ylm_in_complex_ylm(l)
            for i in range(n):
                for j in range(n):
                    acc = sum(C[i][k] * C[j][k].conjugate() for k in range(n))
                    expected = 1.0 if i == j else 0.0
                    assert near_abs(acc.real, expected, atol=1e-14)
                    assert near_abs(acc.imag, 0.0,      atol=1e-14)

    def test_precision_keyword(self):
        # All three precisions return the same Python-complex matrix.
        C_d = real_ylm_in_complex_ylm(1, precision='double')
        C_f = real_ylm_in_complex_ylm(1, precision='float')
        C_l = real_ylm_in_complex_ylm(1, precision='longdouble')
        for i in range(3):
            for j in range(3):
                assert near_abs(C_f[i][j].real, C_d[i][j].real, atol=5e-7)
                assert near_abs(C_f[i][j].imag, C_d[i][j].imag, atol=5e-7)
                assert near_abs(C_l[i][j].real, C_d[i][j].real, atol=1e-15)
                assert near_abs(C_l[i][j].imag, C_d[i][j].imag, atol=1e-15)

    def test_negative_l_raises(self):
        with pytest.raises(ValueError):
            real_ylm_in_complex_ylm(-1)
