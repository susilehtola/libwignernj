# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""
wignernj — exact Wigner 3j/6j/9j symbols and related coefficients via prime factorization.

All functions accept angular momentum arguments as integers, floats
(half-integers), or fractions.Fraction objects.  An optional keyword
argument ``precision`` selects the output type:

  'float'      -- single precision (~7 digits)
  'double'     -- double precision (~15 digits, default)
  'longdouble' -- extended/quad precision where available

Returns 0.0 if selection rules are violated (not an error).

Phase conventions
-----------------
The Wigner 3j, 6j, and 9j symbols, the Clebsch-Gordan coefficient,
the Racah W coefficient, and the Fano X-coefficient are pure SU(2)
algebraic objects -- their values
are fixed by the Racah/Wigner combinatorial formulas alone, with no
spherical-harmonic phase convention entering anywhere.  The
Clebsch-Gordan coefficient uses the Condon-Shortley sign convention of
Edmonds (1957) and Varshalovich et al. (1988).  The ``gaunt`` and
``gaunt_real`` routines assume the Condon-Shortley phase for Y_l^m;
``gaunt_real`` further fixes the real spherical harmonics by the
Wikipedia/Condon-Shortley construction.  See the header comment of the
underlying C library (``include/wignernj.h``) for the explicit formulas.

Example::

    import wignernj
    wignernj.wigner3j(1, 1, 0, 0, 0, 0)          # -1/sqrt(3) ≈ -0.5774
    wignernj.wigner3j(0.5, 0.5, 1, 0.5, -0.5, 0) # 1/sqrt(6) ≈ 0.4082
    wignernj.gaunt(2, 1, 2, -1, 2, 0)            # complex Y_l^m Gaunt
    wignernj.gaunt_real(2, 1, 2, -1, 0, 0)       # real-Y Gaunt

Functions exposed: ``wigner3j``, ``wigner6j``, ``wigner9j``,
``clebsch_gordan``, ``racah_w``, ``fano_x``, ``gaunt``, ``gaunt_real``.
"""
__version__ = "0.5.0"

from ._wignernj import (
    wigner3j,
    wigner6j,
    wigner9j,
    clebsch_gordan,
    racah_w,
    fano_x,
    gaunt,
    gaunt_real,
)

__all__ = ["wigner3j", "wigner6j", "wigner9j",
           "clebsch_gordan", "racah_w", "fano_x", "gaunt", "gaunt_real"]
