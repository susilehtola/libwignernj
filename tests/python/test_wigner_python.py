# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# pytest suite for the Python wigner extension.
# Run from the build directory after: pip install -e . --no-build-isolation
# or after cmake --build build with PYTHONPATH set.

import math
import pytest

try:
    from wigner import (wigner3j, wigner6j, wigner9j,
                        clebsch_gordan, racah_w, gaunt)
except ImportError:
    pytest.skip("wigner extension not installed", allow_module_level=True)

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
