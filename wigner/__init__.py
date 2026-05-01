# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""
wigner — exact Wigner 3j/6j/9j symbols and related coefficients.

All functions accept angular momentum arguments as integers, floats
(half-integers), or fractions.Fraction objects.  An optional keyword
argument ``precision`` selects the output type:

  'float'      -- single precision (~7 digits)
  'double'     -- double precision (~15 digits, default)
  'longdouble' -- extended/quad precision where available

Returns 0.0 if selection rules are violated (not an error).

Example::

    import wigner
    wigner.wigner3j(1, 1, 0, 0, 0, 0)          # -1/sqrt(3) ≈ -0.5774
    wigner.wigner3j(0.5, 0.5, 1, 0.5, -0.5, 0) # 1/sqrt(6) ≈ 0.4082
"""
__version__ = "0.1.0"

from ._wigner import (
    wigner3j,
    wigner6j,
    wigner9j,
    clebsch_gordan,
    racah_w,
    gaunt,
)

__all__ = ["wigner3j", "wigner6j", "wigner9j",
           "clebsch_gordan", "racah_w", "gaunt"]
