#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# Generate C test files with sympy reference values for libwignernj.
# Run from the repository root:  python3 tests/gen_refs.py
#
# Requires: pip install sympy

import random
import os
import sys
import time

random.seed(42)

try:
    import sympy
    from sympy.physics.wigner import (
        wigner_3j, wigner_6j, wigner_9j, clebsch_gordan,
        gaunt as _sympy_gaunt, racah as _sympy_racah
    )
    from sympy import Rational, N
except ImportError:
    print("sympy is required: pip install sympy", file=sys.stderr)
    sys.exit(1)

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ── selection-rule helpers (all args in 2*j integer convention) ─────────────

def triangle_ok(ta, tb, tc):
    if (ta + tb + tc) & 1:
        return False
    if tc < abs(ta - tb):
        return False
    if tc > ta + tb:
        return False
    return True

def sel_3j(tj1, tj2, tj3, tm1, tm2, tm3):
    if tm1 + tm2 + tm3 != 0:
        return False
    if abs(tm1) > tj1 or abs(tm2) > tj2 or abs(tm3) > tj3:
        return False
    if (tj1 - tm1) & 1 or (tj2 - tm2) & 1 or (tj3 - tm3) & 1:
        return False
    return triangle_ok(tj1, tj2, tj3)

def sel_6j(tj1, tj2, tj3, tj4, tj5, tj6):
    return (triangle_ok(tj1, tj2, tj3) and
            triangle_ok(tj1, tj5, tj6) and
            triangle_ok(tj4, tj2, tj6) and
            triangle_ok(tj4, tj5, tj3))

def sel_9j(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33):
    return (triangle_ok(tj11, tj12, tj13) and
            triangle_ok(tj21, tj22, tj23) and
            triangle_ok(tj31, tj32, tj33) and
            triangle_ok(tj11, tj21, tj31) and
            triangle_ok(tj12, tj22, tj32) and
            triangle_ok(tj13, tj23, tj33))

def sel_cg(tj1, tm1, tj2, tm2, tJ, tM):
    if tm1 + tm2 != tM:
        return False
    if abs(tm1) > tj1 or abs(tm2) > tj2 or abs(tM) > tJ:
        return False
    if (tj1 - tm1) & 1 or (tj2 - tm2) & 1 or (tJ - tM) & 1:
        return False
    return triangle_ok(tj1, tj2, tJ)

def sel_racah(tj1, tj2, tJ, tj3, tj12, tj23):
    # W(j1,j2,J,j3;j12,j23) nonzero only if {j1,j2,j12;j3,J,j23} is valid
    return sel_6j(tj1, tj2, tj12, tj3, tJ, tj23)

def sel_gaunt(tl1, tm1, tl2, tm2, tl3, tm3):
    if tl1 & 1 or tl2 & 1 or tl3 & 1:   # l must be integer
        return False
    if tm1 + tm2 + tm3 != 0:
        return False
    if abs(tm1) > tl1 or abs(tm2) > tl2 or abs(tm3) > tl3:
        return False
    if (tl1 - tm1) & 1 or (tl2 - tm2) & 1 or (tl3 - tm3) & 1:
        return False
    if not triangle_ok(tl1, tl2, tl3):
        return False
    # parity: l1+l2+l3 must be even for 3j(l;0,0,0) to be non-zero
    if ((tl1 + tl2 + tl3) // 2) & 1:
        return False
    return True

def sel_gaunt_real(tl1, tm1, tl2, tm2, tl3, tm3):
    """Selection rules for the real-spherical-harmonic Gaunt coefficient.
    The l constraints are the same as for the complex Gaunt; m constraints
    are |m_i| <= l_i and integer m, but no fixed sum rule (the implicit
    rule comes from the absolute-value sign tuple algebra).  Both signs of
    m_i are allowed here (cosine vs sine real harmonics)."""
    if tl1 & 1 or tl2 & 1 or tl3 & 1:    # l must be integer
        return False
    if tm1 & 1 or tm2 & 1 or tm3 & 1:    # m must be integer
        return False
    if abs(tm1) > tl1 or abs(tm2) > tl2 or abs(tm3) > tl3:
        return False
    if not triangle_ok(tl1, tl2, tl3):
        return False
    # parity: l1+l2+l3 must be even for 3j(l;0,0,0) to be non-zero
    if ((tl1 + tl2 + tl3) // 2) & 1:
        return False
    return True

# ── sympy wrappers ──────────────────────────────────────────────────────────

def R(tx):
    return Rational(tx, 2)

def sym_3j(tj1, tj2, tj3, tm1, tm2, tm3):
    return wigner_3j(R(tj1), R(tj2), R(tj3), R(tm1), R(tm2), R(tm3))

def sym_6j(tj1, tj2, tj3, tj4, tj5, tj6):
    return wigner_6j(R(tj1), R(tj2), R(tj3), R(tj4), R(tj5), R(tj6))

def sym_9j(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33):
    return wigner_9j(R(tj11), R(tj12), R(tj13),
                     R(tj21), R(tj22), R(tj23),
                     R(tj31), R(tj32), R(tj33))

def sym_cg(tj1, tm1, tj2, tm2, tJ, tM):
    # sympy clebsch_gordan(j1, j2, j3, m1, m2, m3) = <j1 m1; j2 m2 | j3 m3>
    return clebsch_gordan(R(tj1), R(tj2), R(tJ), R(tm1), R(tm2), R(tM))

def sym_racah(tj1, tj2, tJ, tj3, tj12, tj23):
    # sympy racah(aa,bb,cc,dd,ee,ff) = W(aa,bb,cc,dd;ee,ff)
    return _sympy_racah(R(tj1), R(tj2), R(tJ), R(tj3), R(tj12), R(tj23))

def sym_gaunt(tl1, tm1, tl2, tm2, tl3, tm3):
    # sympy gaunt(l1,l2,l3,m1,m2,m3) with integer l values
    return _sympy_gaunt(tl1 // 2, tl2 // 2, tl3 // 2,
                        tm1 // 2, tm2 // 2, tm3 // 2)

def sym_gaunt_real(tl1, tm1, tl2, tm2, tl3, tm3):
    """Real-spherical-harmonic Gaunt coefficient via the unitary
    transform between real and complex harmonics:

        S_{l, 0}   = Y_l^0
        S_{l,m>0}  = (1/sqrt(2)) (Y_l^{-m} + (-1)^m Y_l^m)
        S_{l,m<0}  = (i/sqrt(2)) (Y_l^{ m} - (-1)^|m| Y_l^{-m})

    expanded into a sum (≤ 2 terms) of complex Gaunts.  Returns a
    sympy expression (which simplifies to a real number)."""
    from sympy import sqrt as ssqrt, I, S as _S
    ls = [tl1 // 2, tl2 // 2, tl3 // 2]
    ms = [tm1 // 2, tm2 // 2, tm3 // 2]
    a = [abs(m) for m in ms]
    sigma = [(1 if m > 0 else (-1 if m < 0 else 0)) for m in ms]

    def Tsym(sigma_i, s_i, a_i):
        if a_i == 0:
            return _S(1)
        sign_a = -_S(1) if (a_i & 1) else _S(1)
        if sigma_i > 0:
            return (sign_a if s_i > 0 else _S(1)) / ssqrt(2)
        else:
            return (-I * sign_a if s_i > 0 else I) / ssqrt(2)

    total = _S(0)
    for s1 in ([1, -1] if a[0] else [1]):
        for s2 in ([1, -1] if a[1] else [1]):
            for s3 in ([1, -1] if a[2] else [1]):
                if s1*a[0] + s2*a[1] + s3*a[2] != 0:
                    continue
                c = Tsym(sigma[0], s1, a[0]) * \
                    Tsym(sigma[1], s2, a[1]) * \
                    Tsym(sigma[2], s3, a[2])
                gc = _sympy_gaunt(ls[0], ls[1], ls[2],
                                   s1*a[0], s2*a[1], s3*a[2])
                total += c * gc
    # Result is mathematically real; take the real part to drop any
    # leftover imaginary symbolic noise.
    from sympy import re
    return re(total)

# ── random valid-tuple generators ──────────────────────────────────────────

def rand_3j_jvalues(tj_max):
    """Return (tj1,tj2,tj3) satisfying triangle + same parity."""
    tj1 = random.randrange(0, tj_max + 1)
    tj2 = random.randrange(0, tj_max + 1)
    lo = abs(tj1 - tj2); hi = min(tj1 + tj2, tj_max)
    choices = [t for t in range(lo, hi + 1) if (t % 2) == (lo % 2)]
    if not choices:
        return None
    tj3 = random.choice(choices)
    return tj1, tj2, tj3

def rand_3j_mvalues(tj1, tj2, tj3):
    """Return (tm1,tm2,tm3) satisfying 3j selection rules."""
    tm1 = random.choice(list(range(-tj1, tj1 + 1, 2)))
    tm2 = random.choice(list(range(-tj2, tj2 + 1, 2)))
    tm3 = -tm1 - tm2
    if not sel_3j(tj1, tj2, tj3, tm1, tm2, tm3):
        return None
    return tm1, tm2, tm3

def rand_6j_jvalues(tj_max):
    """Return (tj1,...,tj6) satisfying all four triangle conditions."""
    # Build valid tuple by sampling j-values that satisfy all 4 triangles
    tj1 = random.randrange(0, tj_max + 1)
    tj2 = random.randrange(0, tj_max + 1)
    # tj3: valid for (j1,j2,j3)
    lo = abs(tj1 - tj2); hi = min(tj1 + tj2, tj_max)
    c3 = [t for t in range(lo, hi + 1) if (t % 2) == (lo % 2)]
    if not c3: return None
    tj3 = random.choice(c3)
    # tj4, tj5: freely chosen
    tj4 = random.randrange(0, tj_max + 1)
    # tj5: valid for (j1,j5,j6) and (j4,j5,j3)
    tj5 = random.randrange(0, tj_max + 1)
    # tj6: must satisfy (j1,j5,j6), (j4,j2,j6), (j4,j5,j3) already fixed
    # Try to find valid tj6
    lo6a = abs(tj1 - tj5); hi6a = tj1 + tj5
    lo6b = abs(tj4 - tj2); hi6b = tj4 + tj2
    lo6 = max(lo6a, lo6b); hi6 = min(hi6a, hi6b, tj_max)
    # also parity constraints
    c6 = [t for t in range(lo6, hi6 + 1)
          if (t % 2) == (lo6a % 2) and (t % 2) == (lo6b % 2)
          and triangle_ok(tj4, tj5, tj3)]
    if not c6: return None
    # filter: (j4,j5,j3) must also hold (already checked via tj3 choice)
    if not triangle_ok(tj4, tj5, tj3): return None
    c6 = [t for t in c6 if triangle_ok(tj1, tj5, t) and triangle_ok(tj4, tj2, t)]
    if not c6: return None
    tj6 = random.choice(c6)
    if not sel_6j(tj1, tj2, tj3, tj4, tj5, tj6): return None
    return tj1, tj2, tj3, tj4, tj5, tj6

def rand_9j_jvalues(tj_max):
    """Return (tj11,...,tj33) satisfying all six triangle conditions."""
    for _ in range(50):  # limited attempts per call
        args = [random.randrange(0, tj_max + 1) for _ in range(9)]
        if sel_9j(*args):
            return tuple(args)
    return None

# ── case generators ─────────────────────────────────────────────────────────

def gen_3j(target_small, target_large, tj_small=10, tj_large=200):
    """Enumerate ALL non-zero 3j for j≤5 then sample larger j."""
    print("  3j: enumerating j ≤ 5 fully...", flush=True)
    small_all = []
    for tj1 in range(0, tj_small + 1):
        for tj2 in range(0, tj_small + 1):
            lo = abs(tj1 - tj2)
            for tj3 in range(lo, min(tj1 + tj2, tj_small) + 1):
                if not triangle_ok(tj1, tj2, tj3):
                    continue
                for tm1 in range(-tj1, tj1 + 1, 2):
                    for tm2 in range(-tj2, tj2 + 1, 2):
                        tm3 = -tm1 - tm2
                        if not sel_3j(tj1, tj2, tj3, tm1, tm2, tm3):
                            continue
                        v = sym_3j(tj1, tj2, tj3, tm1, tm2, tm3)
                        if v == 0:
                            continue
                        small_all.append(
                            (tj1, tj2, tj3, tm1, tm2, tm3, float(N(v, 34))))
    print(f"    found {len(small_all)} non-zero small-j cases", flush=True)
    random.shuffle(small_all)
    small = small_all[:target_small]

    print("    collecting large-j cases...", flush=True)
    seen = set((c[:6] for c in small_all))
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 50:
        attempts += 1
        jv = rand_3j_jvalues(tj_large)
        if jv is None:
            continue
        if max(jv) <= tj_small:
            continue
        mv = rand_3j_mvalues(*jv)
        if mv is None:
            continue
        args = jv + mv
        if args in seen:
            continue
        seen.add(args)
        v = sym_3j(*args)
        if v == 0:
            continue
        large.append(args + (float(N(v, 34)),))
    print(f"    → {len(small)} small + {len(large)} large = {len(small)+len(large)} total",
          flush=True)
    return small + large


def gen_6j(target_small, target_large, tj_small=10, tj_large=100):
    print("  6j: random sampling...", flush=True)
    seen = set()
    small = []
    attempts = 0
    while len(small) < target_small and attempts < target_small * 500:
        attempts += 1
        jv = rand_6j_jvalues(tj_small)
        if jv is None or jv in seen:
            continue
        seen.add(jv)
        v = sym_6j(*jv)
        if v == 0:
            continue
        small.append(jv + (float(N(v, 34)),))
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 500:
        attempts += 1
        jv = rand_6j_jvalues(tj_large)
        if jv is None or jv in seen or max(jv) <= tj_small:
            continue
        seen.add(jv)
        v = sym_6j(*jv)
        if v == 0:
            continue
        large.append(jv + (float(N(v, 34)),))
    total = len(small) + len(large)
    print(f"    → {len(small)} small + {len(large)} large = {total} total",
          flush=True)
    return small + large


def gen_9j(target_small, target_large, tj_small=10, tj_large=30):
    print("  9j: random sampling...", flush=True)
    seen = set()
    small = []
    attempts = 0
    while len(small) < target_small and attempts < target_small * 2000:
        attempts += 1
        jv = rand_9j_jvalues(tj_small)
        if jv is None or jv in seen:
            continue
        seen.add(jv)
        v = sym_9j(*jv)
        if v == 0:
            continue
        small.append(jv + (float(N(v, 34)),))
        if attempts % 5000 == 0:
            print(f"    ... {len(small)} so far, {attempts} attempts",
                  flush=True)
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 2000:
        attempts += 1
        jv = rand_9j_jvalues(tj_large)
        if jv is None or jv in seen or max(jv) <= tj_small:
            continue
        seen.add(jv)
        v = sym_9j(*jv)
        if v == 0:
            continue
        large.append(jv + (float(N(v, 34)),))
    total = len(small) + len(large)
    print(f"    → {len(small)} small + {len(large)} large = {total} total",
          flush=True)
    return small + large


def gen_cg(target_small, target_large, tj_small=10, tj_large=200):
    """Enumerate all non-zero CG for j≤5, then sample larger j."""
    print("  CG: enumerating j ≤ 5 fully...", flush=True)
    small_all = []
    for tj1 in range(0, tj_small + 1):
        for tj2 in range(0, tj_small + 1):
            lo = abs(tj1 - tj2)
            for tJ in range(lo, min(tj1 + tj2, tj_small) + 1):
                if not triangle_ok(tj1, tj2, tJ):
                    continue
                for tm1 in range(-tj1, tj1 + 1, 2):
                    for tm2 in range(-tj2, tj2 + 1, 2):
                        tM = tm1 + tm2
                        if not sel_cg(tj1, tm1, tj2, tm2, tJ, tM):
                            continue
                        v = sym_cg(tj1, tm1, tj2, tm2, tJ, tM)
                        if v == 0:
                            continue
                        small_all.append(
                            (tj1, tm1, tj2, tm2, tJ, tM, float(N(v, 34))))
    print(f"    found {len(small_all)} non-zero small-j cases", flush=True)
    random.shuffle(small_all)
    small = small_all[:target_small]

    print("    collecting large-j cases...", flush=True)
    seen = set(c[:6] for c in small_all)
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 50:
        attempts += 1
        jv = rand_3j_jvalues(tj_large)
        if jv is None:
            continue
        tj1, tj2, tJ = jv
        if max(jv) <= tj_small:
            continue
        tm1_choices = list(range(-tj1, tj1 + 1, 2))
        tm2_choices = list(range(-tj2, tj2 + 1, 2))
        if not tm1_choices or not tm2_choices:
            continue
        tm1 = random.choice(tm1_choices)
        tm2 = random.choice(tm2_choices)
        tM = tm1 + tm2
        args = (tj1, tm1, tj2, tm2, tJ, tM)
        if not sel_cg(*args) or args in seen:
            continue
        seen.add(args)
        v = sym_cg(*args)
        if v == 0:
            continue
        large.append(args + (float(N(v, 34)),))
    print(f"    → {len(small)} small + {len(large)} large = {len(small)+len(large)} total",
          flush=True)
    return small + large


def gen_racah(target_small, target_large, tj_small=10, tj_large=60):
    print("  Racah W: random sampling...", flush=True)
    seen = set()
    small = []
    attempts = 0
    while len(small) < target_small and attempts < target_small * 500:
        attempts += 1
        # W(j1,j2,J,j3;j12,j23): generate via 6j structure
        jv = rand_6j_jvalues(tj_small)
        if jv is None:
            continue
        tj1, tj2, tj12, tj3, tJ, tj23 = jv
        args = (tj1, tj2, tJ, tj3, tj12, tj23)
        if args in seen or not sel_racah(*args):
            continue
        seen.add(args)
        v = sym_racah(*args)
        if v == 0:
            continue
        small.append(args + (float(N(v, 34)),))
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 500:
        attempts += 1
        jv = rand_6j_jvalues(tj_large)
        if jv is None:
            continue
        tj1, tj2, tj12, tj3, tJ, tj23 = jv
        args = (tj1, tj2, tJ, tj3, tj12, tj23)
        if args in seen or max(args) <= tj_small or not sel_racah(*args):
            continue
        seen.add(args)
        v = sym_racah(*args)
        if v == 0:
            continue
        large.append(args + (float(N(v, 34)),))
    total = len(small) + len(large)
    print(f"    → {len(small)} small + {len(large)} large = {total} total",
          flush=True)
    return small + large


def gen_gaunt(target_small, target_large, tl_small=10, tl_large=40):
    """Enumerate non-zero Gaunt for l≤5, then sample larger l."""
    print("  Gaunt: enumerating l ≤ 5 fully...", flush=True)
    small_all = []
    for tl1 in range(0, tl_small + 1, 2):
        for tl2 in range(0, tl_small + 1, 2):
            lo = abs(tl1 - tl2)
            for tl3 in range(lo, min(tl1 + tl2, tl_small) + 1, 2):
                if not triangle_ok(tl1, tl2, tl3):
                    continue
                if ((tl1 + tl2 + tl3) // 2) & 1:   # l1+l2+l3 odd → zero
                    continue
                for tm1 in range(-tl1, tl1 + 1, 2):
                    for tm2 in range(-tl2, tl2 + 1, 2):
                        tm3 = -tm1 - tm2
                        if not sel_gaunt(tl1, tm1, tl2, tm2, tl3, tm3):
                            continue
                        v = sym_gaunt(tl1, tm1, tl2, tm2, tl3, tm3)
                        if v == 0:
                            continue
                        small_all.append(
                            (tl1, tm1, tl2, tm2, tl3, tm3, float(N(v, 34))))
    print(f"    found {len(small_all)} non-zero small-l cases", flush=True)
    random.shuffle(small_all)
    small = small_all[:target_small]

    print("    collecting large-l cases...", flush=True)
    seen = set(c[:6] for c in small_all)
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 200:
        attempts += 1
        tl1 = random.randrange(0, tl_large // 2 + 1) * 2
        tl2 = random.randrange(0, tl_large // 2 + 1) * 2
        lo = abs(tl1 - tl2); hi = min(tl1 + tl2, tl_large)
        # step 4 to keep l1+l2+l3 even: lo has same parity as tl1+tl2 (even)
        choices = [t for t in range(lo, hi + 1, 2)
                   if not (((tl1 + tl2 + t) // 2) & 1)]
        if not choices:
            continue
        tl3 = random.choice(choices)
        if max(tl1, tl2, tl3) <= tl_small:
            continue
        if not triangle_ok(tl1, tl2, tl3):
            continue
        tm1_ch = list(range(-tl1, tl1 + 1, 2))
        tm2_ch = list(range(-tl2, tl2 + 1, 2))
        if not tm1_ch or not tm2_ch:
            continue
        tm1 = random.choice(tm1_ch)
        tm2 = random.choice(tm2_ch)
        tm3 = -tm1 - tm2
        args = (tl1, tm1, tl2, tm2, tl3, tm3)
        if not sel_gaunt(*args) or args in seen:
            continue
        seen.add(args)
        v = sym_gaunt(*args)
        if v == 0:
            continue
        large.append(args + (float(N(v, 34)),))
    print(f"    → {len(small)} small + {len(large)} large = {len(small)+len(large)} total",
          flush=True)
    return small + large


def gen_gaunt_real(target_small, target_large, tl_small=10, tl_large=40):
    """Reference values for the real-spherical-harmonic Gaunt coefficient.

    For l ≤ 5 we enumerate every (l1,m1,l2,m2,l3,m3) that satisfies the
    selection rules and pick a random subset of non-zero values.  For
    l up to tl_large we sample randomly.  The reference value is computed
    from sympy via sym_gaunt_real (decomposition into complex Gaunts).
    """
    print("  Gaunt-real: enumerating l ≤ 5 fully...", flush=True)
    small_all = []
    for tl1 in range(0, tl_small + 1, 2):
        for tl2 in range(0, tl_small + 1, 2):
            lo = abs(tl1 - tl2)
            for tl3 in range(lo, min(tl1 + tl2, tl_small) + 1, 2):
                if not triangle_ok(tl1, tl2, tl3):
                    continue
                if ((tl1 + tl2 + tl3) // 2) & 1:    # l1+l2+l3 odd → zero
                    continue
                for tm1 in range(-tl1, tl1 + 1, 2):
                    for tm2 in range(-tl2, tl2 + 1, 2):
                        for tm3 in range(-tl3, tl3 + 1, 2):
                            if not sel_gaunt_real(tl1, tm1, tl2, tm2, tl3, tm3):
                                continue
                            v = sym_gaunt_real(tl1, tm1, tl2, tm2, tl3, tm3)
                            if v == 0:
                                continue
                            small_all.append(
                                (tl1, tm1, tl2, tm2, tl3, tm3, float(N(v, 34))))
    print(f"    found {len(small_all)} non-zero small-l cases", flush=True)
    random.shuffle(small_all)
    small = small_all[:target_small]

    print("    collecting large-l cases...", flush=True)
    seen = set(c[:6] for c in small_all)
    large = []
    attempts = 0
    while len(large) < target_large and attempts < target_large * 300:
        attempts += 1
        tl1 = random.randrange(0, tl_large // 2 + 1) * 2
        tl2 = random.randrange(0, tl_large // 2 + 1) * 2
        lo = abs(tl1 - tl2); hi = min(tl1 + tl2, tl_large)
        choices = [t for t in range(lo, hi + 1, 2)
                   if not (((tl1 + tl2 + t) // 2) & 1)]
        if not choices:
            continue
        tl3 = random.choice(choices)
        if max(tl1, tl2, tl3) <= tl_small:
            continue
        if not triangle_ok(tl1, tl2, tl3):
            continue
        tm1 = random.choice(list(range(-tl1, tl1 + 1, 2)))
        tm2 = random.choice(list(range(-tl2, tl2 + 1, 2)))
        tm3 = random.choice(list(range(-tl3, tl3 + 1, 2)))
        args = (tl1, tm1, tl2, tm2, tl3, tm3)
        if not sel_gaunt_real(*args) or args in seen:
            continue
        seen.add(args)
        v = sym_gaunt_real(*args)
        if v == 0:
            continue
        large.append(args + (float(N(v, 34)),))
    print(f"    → {len(small)} small + {len(large)} large = "
          f"{len(small)+len(large)} total", flush=True)
    return small + large


# ── C file emitters ─────────────────────────────────────────────────────────

SPDX = """\
/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * AUTO-GENERATED by tests/gen_refs.py -- DO NOT EDIT
 * Reference values computed with sympy {ver}.
 */
"""

def emit_3j(cases):
    path = os.path.join(OUTDIR, "test_3j.c")
    with open(path, "w") as f:
        f.write(SPDX.format(ver=sympy.__version__))
        f.write('#include "run_tests.h"\n')
        f.write('#include "../include/wignernj.h"\n\n')
        f.write("typedef struct { int tj1,tj2,tj3,tm1,tm2,tm3; double val; } w3j_t;\n")
        f.write("static const w3j_t g_3j[] = {\n")
        for c in cases:
            f.write("    {%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")
        f.write("int main(void)\n{\n    int i;\n")
        f.write("    for (i=0;i<(int)(sizeof(g_3j)/sizeof(g_3j[0]));i++) {\n")
        f.write("        const w3j_t *c=&g_3j[i];\n")
        f.write("        double got=wigner3j(c->tj1,c->tj2,c->tj3,c->tm1,c->tm2,c->tm3);\n")
        f.write("        TEST_NEAR(got,c->val,1e-12);\n    }\n")
        f.write("    SUMMARY();\n}\n")
    print(f"  Wrote {path} ({len(cases)} cases)")


def emit_6j(cases):
    path = os.path.join(OUTDIR, "test_6j.c")
    with open(path, "w") as f:
        f.write(SPDX.format(ver=sympy.__version__))
        f.write('#include "run_tests.h"\n')
        f.write('#include "../include/wignernj.h"\n\n')
        f.write("typedef struct { int tj1,tj2,tj3,tj4,tj5,tj6; double val; } w6j_t;\n")
        f.write("static const w6j_t g_6j[] = {\n")
        for c in cases:
            f.write("    {%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")
        f.write("int main(void)\n{\n    int i;\n")
        f.write("    for (i=0;i<(int)(sizeof(g_6j)/sizeof(g_6j[0]));i++) {\n")
        f.write("        const w6j_t *c=&g_6j[i];\n")
        f.write("        double got=wigner6j(c->tj1,c->tj2,c->tj3,c->tj4,c->tj5,c->tj6);\n")
        f.write("        TEST_NEAR(got,c->val,1e-12);\n    }\n")
        f.write("    SUMMARY();\n}\n")
    print(f"  Wrote {path} ({len(cases)} cases)")


def emit_9j(cases):
    path = os.path.join(OUTDIR, "test_9j.c")
    with open(path, "w") as f:
        f.write(SPDX.format(ver=sympy.__version__))
        f.write('#include "run_tests.h"\n')
        f.write('#include "../include/wignernj.h"\n\n')
        f.write("typedef struct {\n")
        f.write("    int tj11,tj12,tj13,tj21,tj22,tj23,tj31,tj32,tj33;\n")
        f.write("    double val;\n} w9j_t;\n")
        f.write("static const w9j_t g_9j[] = {\n")
        for c in cases:
            f.write("    {%d,%d,%d,%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")
        f.write("int main(void)\n{\n    int i;\n")
        f.write("    for (i=0;i<(int)(sizeof(g_9j)/sizeof(g_9j[0]));i++) {\n")
        f.write("        const w9j_t *c=&g_9j[i];\n")
        f.write("        double got=wigner9j(c->tj11,c->tj12,c->tj13,\n")
        f.write("                            c->tj21,c->tj22,c->tj23,\n")
        f.write("                            c->tj31,c->tj32,c->tj33);\n")
        f.write("        TEST_NEAR(got,c->val,1e-12);\n    }\n")
        f.write("    SUMMARY();\n}\n")
    print(f"  Wrote {path} ({len(cases)} cases)")


def emit_derived(cg_cases, racah_cases, gaunt_cases):
    path = os.path.join(OUTDIR, "test_derived.c")
    with open(path, "w") as f:
        f.write(SPDX.format(ver=sympy.__version__))
        f.write('#include "run_tests.h"\n')
        f.write('#include "../include/wignernj.h"\n\n')
        f.write("/* Clebsch-Gordan */\n")
        f.write("typedef struct{int tj1,tm1,tj2,tm2,tJ,tM;double val;}cg_t;\n")
        f.write("static const cg_t g_cg[]={\n")
        for c in cg_cases:
            f.write("    {%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")

        f.write("/* Racah W */\n")
        f.write("typedef struct{int tj1,tj2,tJ,tj3,tj12,tj23;double val;}racah_t;\n")
        f.write("static const racah_t g_racah[]={\n")
        for c in racah_cases:
            f.write("    {%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")

        f.write("/* Gaunt */\n")
        f.write("typedef struct{int tl1,tm1,tl2,tm2,tl3,tm3;double val;}gaunt_t;\n")
        f.write("static const gaunt_t g_gaunt[]={\n")
        for c in gaunt_cases:
            f.write("    {%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")

        f.write("int main(void)\n{\n    int i;\n\n")

        f.write("    /* Clebsch-Gordan */\n")
        f.write("    for(i=0;i<(int)(sizeof(g_cg)/sizeof(g_cg[0]));i++){\n")
        f.write("        const cg_t *c=&g_cg[i];\n")
        f.write("        double got=clebsch_gordan(c->tj1,c->tm1,c->tj2,c->tm2,c->tJ,c->tM);\n")
        f.write("        TEST_NEAR(got,c->val,1e-12);\n    }\n\n")

        f.write("    /* Racah W */\n")
        f.write("    for(i=0;i<(int)(sizeof(g_racah)/sizeof(g_racah[0]));i++){\n")
        f.write("        const racah_t *c=&g_racah[i];\n")
        f.write("        double got=racah_w(c->tj1,c->tj2,c->tJ,c->tj3,c->tj12,c->tj23);\n")
        f.write("        TEST_NEAR(got,c->val,1e-12);\n    }\n\n")

        f.write("    /* Gaunt */\n")
        f.write("    for(i=0;i<(int)(sizeof(g_gaunt)/sizeof(g_gaunt[0]));i++){\n")
        f.write("        const gaunt_t *c=&g_gaunt[i];\n")
        f.write("        double got=gaunt(c->tl1,c->tm1,c->tl2,c->tm2,c->tl3,c->tm3);\n")
        f.write("        TEST_NEAR(got,c->val,1e-12);\n    }\n\n")

        f.write("    SUMMARY();\n}\n")
    total = len(cg_cases) + len(racah_cases) + len(gaunt_cases)
    print(f"  Wrote {path} ({len(cg_cases)} CG + {len(racah_cases)} Racah + "
          f"{len(gaunt_cases)} Gaunt = {total})")


def emit_gaunt_real(cases):
    path = os.path.join(OUTDIR, "test_gaunt_real.c")
    with open(path, "w") as f:
        f.write(SPDX.format(ver=sympy.__version__))
        f.write('#include "run_tests.h"\n')
        f.write('#include "../include/wignernj.h"\n')
        f.write('#include <math.h>\n\n')

        f.write("/* Hand-derived sanity checks (independent of sympy). */\n"
                "#ifndef M_PI\n#define M_PI 3.14159265358979323846\n#endif\n\n")

        f.write("/* sympy-derived reference values for ∫ S_{l1,m1} S_{l2,m2} S_{l3,m3} dΩ */\n")
        f.write("typedef struct{int tl1,tm1,tl2,tm2,tl3,tm3;double val;}gaunt_real_t;\n")
        f.write("static const gaunt_real_t g_gaunt_real[]={\n")
        for c in cases:
            f.write("    {%d,%d,%d,%d,%d,%d, %s},\n" % c)
        f.write("};\n\n")

        f.write("int main(void)\n{\n    int i;\n    double sqrt_pi = sqrt(M_PI);\n\n")

        f.write("    /* (a) Hand-derived reference cases. */\n")
        f.write("    TEST_NEAR(gaunt_real(0,0,0,0,0,0), 1.0/(2*sqrt_pi), 1e-15);\n")
        f.write("    TEST_NEAR(gaunt_real(2,0,2,0,4,0), 1.0/sqrt(5*M_PI), 1e-14);\n")
        f.write("    TEST_NEAR(gaunt_real(2,2,2,2,0,0), 1.0/(2*sqrt_pi), 1e-14);\n")
        f.write("    TEST_NEAR(gaunt_real(2,-2,2,-2,0,0), 1.0/(2*sqrt_pi), 1e-14);\n")
        f.write("    TEST_NEAR(gaunt_real(4,4,4,4,0,0), 1.0/(2*sqrt_pi), 1e-14);\n")
        f.write("    /* parity / triangle violations and pairwise cancellations: */\n")
        f.write("    TEST_ABS(gaunt_real(2,0,2,0,2,0), 0.0, 1e-15);   /* l1+l2+l3 odd */\n")
        f.write("    TEST_ABS(gaunt_real(2,2,2,-2,0,0), 0.0, 1e-15);  /* cos*sin*1 vanishes */\n")
        f.write("    TEST_ABS(gaunt_real(2,2,2,2,0,0)-gaunt_real(2,2,2,2,0,0), 0.0, 1e-15);\n")
        f.write("    TEST_ABS(gaunt_real(2,-2,2,2,2,2), 0.0, 1e-15);  /* odd # of m_i<0 */\n\n")

        f.write("    /* (b) Permutation symmetry: G^R is symmetric in the 3 (l_i, m_i). */\n")
        f.write("    {\n")
        f.write("        double a = gaunt_real(2,0,4,2,4,-2);\n")
        f.write("        double b = gaunt_real(4,2,2,0,4,-2);\n")
        f.write("        double c = gaunt_real(4,-2,4,2,2,0);\n")
        f.write("        TEST_NEAR(a,b,1e-14);\n")
        f.write("        TEST_NEAR(a,c,1e-14);\n")
        f.write("    }\n\n")

        f.write("    /* (c) m=0: real Gaunt agrees with the complex Gaunt. */\n")
        f.write("    {\n")
        f.write("        const int cases[][3] = {\n")
        f.write("            {0,0,0},{2,2,0},{4,4,0},{2,2,4},{4,2,4},\n")
        f.write("            {6,4,4},{6,6,0},{6,6,4},{6,6,8},{8,8,8},\n")
        f.write("            {10,8,4},{12,12,0},{12,10,6},{20,16,4}\n")
        f.write("        };\n")
        f.write("        const int n = (int)(sizeof(cases)/sizeof(cases[0]));\n")
        f.write("        for (int j = 0; j < n; j++) {\n")
        f.write("            int l1 = cases[j][0], l2 = cases[j][1], l3 = cases[j][2];\n")
        f.write("            double r = gaunt_real(l1,0,l2,0,l3,0);\n")
        f.write("            double c = gaunt(l1,0,l2,0,l3,0);\n")
        f.write("            TEST_NEAR(r,c,1e-14);\n")
        f.write("        }\n")
        f.write("    }\n\n")

        f.write("    /* (d) sympy reference table */\n")
        f.write("    for (i = 0; i < (int)(sizeof(g_gaunt_real)/sizeof(g_gaunt_real[0])); i++) {\n")
        f.write("        const gaunt_real_t *c = &g_gaunt_real[i];\n")
        f.write("        double got = gaunt_real(c->tl1,c->tm1,c->tl2,c->tm2,c->tl3,c->tm3);\n")
        f.write("        TEST_NEAR(got, c->val, 1e-12);\n")
        f.write("    }\n\n")
        f.write("    SUMMARY();\n}\n")
    print(f"  Wrote {path} ({len(cases)} sympy cases + hand-derived/symmetry checks)")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()

    print("Generating 3j references...")
    c3j = gen_3j(target_small=800, target_large=200, tj_small=10, tj_large=200)
    emit_3j(c3j)

    print("Generating 6j references...")
    c6j = gen_6j(target_small=800, target_large=200, tj_small=10, tj_large=100)
    emit_6j(c6j)

    print("Generating 9j references...")
    c9j = gen_9j(target_small=700, target_large=300, tj_small=10, tj_large=30)
    emit_9j(c9j)

    print("Generating CG references...")
    ccg = gen_cg(target_small=350, target_large=150, tj_small=10, tj_large=200)

    print("Generating Racah W references...")
    crac = gen_racah(target_small=350, target_large=150, tj_small=10, tj_large=60)

    print("Generating Gaunt references...")
    cgnt = gen_gaunt(target_small=250, target_large=100, tl_small=10, tl_large=40)

    emit_derived(ccg, crac, cgnt)

    print("Generating real-Gaunt references...")
    cgr = gen_gaunt_real(target_small=200, target_large=60,
                          tl_small=10, tl_large=40)
    emit_gaunt_real(cgr)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s")
    print(f"  3j:{len(c3j)}  6j:{len(c6j)}  9j:{len(c9j)}")
    print(f"  CG:{len(ccg)}  Racah:{len(crac)}  Gaunt:{len(cgnt)}  GauntReal:{len(cgr)}")


if __name__ == "__main__":
    main()
