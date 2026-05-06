#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""libwignernj Python API demonstration.

Calls every public symbol family exposed by the ``wignernj`` package
once with small, textbook-scale arguments and prints the result
alongside an analytic reference value.  Exits 0 on success, 1 if any
computed value disagrees with its reference by more than 1e-14.

The Python wrapper accepts angular-momentum arguments directly as
integers, floats (half-integers), or ``fractions.Fraction`` objects;
no factor-of-two encoding is needed.  An optional ``precision=``
keyword selects ``'float'``, ``'double'`` (default), or ``'longdouble'``
output.

Run (against an installed ``wignernj`` package):

    pip install -e .   # from the source tree
    python examples/python/all_symbols.py
"""
import math
import sys

import wignernj


TOL = 1e-14
PI  = math.pi


def check(label, computed, expected, tol=TOL):
    print(f"  {label:<38s} = {computed:+.15f}   (expected {expected:+.15f})")
    if abs(computed - expected) > tol:
        sys.stderr.write(
            f"  FAIL: |{computed:g} - {expected:g}| = "
            f"{abs(computed - expected):g} exceeds tolerance {tol:g}\n"
        )
        return 1
    return 0


def main():
    failed = 0

    print("libwignernj Python API demonstration -- one call per symbol family")
    print("------------------------------------------------------------------")

    # 1. Wigner 3j symbol.   ( 1 1 0 )
    #                        ( 0 0 0 )  =  -1/sqrt(3)
    failed += check(
        "wigner3j(1,1,0; 0,0,0)",
        wignernj.wigner3j(1, 1, 0,  0, 0, 0),
        -1.0 / math.sqrt(3.0),
    )

    # 2. Wigner 6j symbol.   { 1 1 1 }
    #                        { 1 1 1 }  =  1/6
    failed += check(
        "wigner6j{1,1,1; 1,1,1}",
        wignernj.wigner6j(1, 1, 1,  1, 1, 1),
        1.0 / 6.0,
    )

    # 3. Wigner 9j symbol.   { 1 1 0 }
    #                        { 1 1 0 }  =  1/3
    #                        { 0 0 0 }
    failed += check(
        "wigner9j{1,1,0; 1,1,0; 0,0,0}",
        wignernj.wigner9j(1, 1, 0,  1, 1, 0,  0, 0, 0),
        1.0 / 3.0,
    )

    # 4. Clebsch-Gordan coefficient.   <1,0; 1,0 | 2,0> = sqrt(2/3)
    failed += check(
        "clebsch_gordan(1,0; 1,0 | 2,0)",
        wignernj.clebsch_gordan(1, 0,  1, 0,  2, 0),
        math.sqrt(2.0 / 3.0),
    )

    # 5. Racah W coefficient.   W(1,1,1,1; 1,1) = 1/6
    failed += check(
        "racah_w(1,1,1,1; 1,1)",
        wignernj.racah_w(1, 1, 1,  1, 1, 1),
        1.0 / 6.0,
    )

    # 6. Fano X coefficient.   X(1,1,1; 1,1,1; 1,1,2) = 1/2
    failed += check(
        "fano_x(1,1,1; 1,1,1; 1,1,2)",
        wignernj.fano_x(1, 1, 1,  1, 1, 1,  1, 1, 2),
        0.5,
    )

    # 7. Gaunt coefficient (complex Y_l^m).
    #    integral(Y_1^0 Y_1^0 Y_2^0) dOmega = 1/sqrt(5*pi)
    failed += check(
        "gaunt(1,0; 1,0; 2,0)",
        wignernj.gaunt(1, 0,  1, 0,  2, 0),
        1.0 / math.sqrt(5.0 * PI),
    )

    # 8. Gaunt coefficient (real spherical harmonics; m=0 -> same as complex).
    failed += check(
        "gaunt_real(1,0; 1,0; 2,0)",
        wignernj.gaunt_real(1, 0,  1, 0,  2, 0),
        1.0 / math.sqrt(5.0 * PI),
    )

    if not failed:
        print("\nAll symbols agree with their analytic references.")
    return failed


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
