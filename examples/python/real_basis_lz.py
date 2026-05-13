#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
"""Build the orbital angular momentum operator l_z in the real
spherical harmonic basis from its diagonal complex-basis form, using
libwignernj's ``real_ylm_in_complex_ylm`` as the change of basis.  Python
parallel to ``examples/c/real_basis_lz.c``.

In the complex Y_l^m basis, l_z is diagonal with entry
``m * delta_{m m'}`` (in units of hbar).  libwignernj returns the
unitary C that encodes the basis-vector relation

    S_{l, m_r}  =  sum_{m_c}  C[m_r, m_c]  Y_l^{m_c}.

The corresponding operator-matrix similarity transform is then

    l_z_real  =  conj(C)  @  l_z_complex  @  transpose(C),

i.e. the textbook ``<S_{m_r}|O|S_{m_r'}>`` expanded out with the ket
using ``C[m_r', m_c']`` and the bra using ``conj(C[m_r, m_c])``.  The
result is purely imaginary and antisymmetric: l_z couples each
cosine-type real harmonic to its sin-type partner with strength
``+/- i`` times the diagonal m_c entry.

For l = 1 this gives the textbook

    l_z = (  0   0   i )
          (  0   0   0 )
          ( -i   0   0 )

in the basis ordering m_r = -1, 0, +1.

Run (against an installed ``wignernj`` package):

    python examples/python/real_basis_lz.py
"""
import sys

import wignernj


def matmul(A, B, n):
    return [[sum(A[r][k] * B[k][c] for k in range(n))
             for c in range(n)] for r in range(n)]


def main():
    l = 1
    n = 2 * l + 1

    # C: basis-vector matrix, S = C Y.  The Python binding returns
    # a list-of-lists indexed as C[m_r + l][m_c + l].
    C = wignernj.real_ylm_in_complex_ylm(l)

    # Cstar: elementwise conjugate.  CT: transpose (not conjugated).
    Cstar = [[C[r][c].conjugate() for c in range(n)] for r in range(n)]
    CT    = [[C[c][r]             for c in range(n)] for r in range(n)]

    # Diagonal complex-basis l_z, m index = m + l.
    Lzc = [[complex(0, 0)] * n for _ in range(n)]
    for m in range(-l, l + 1):
        Lzc[m + l][m + l] = complex(m, 0)

    # Lzr = Cstar @ Lzc @ CT.
    Lzr = matmul(Cstar, matmul(Lzc, CT, n), n)

    def show(name, M):
        print(f"{name} =")
        for row in M:
            print("  " + "  ".join(
                f"{(z.real if abs(z.real) > 1e-15 else 0.0):+.3f}"
                f"{(z.imag if abs(z.imag) > 1e-15 else 0.0):+.3f}i"
                for z in row))
        print()

    print(f"l = {l} real-basis representation of l_z")
    print("Basis ordering: m_r = -l, -l+1, ..., l-1, l")
    print("(real-Y convention matches libwignernj's gaunt_real)")
    print()
    show("C (basis-vector matrix, S = C Y)",       C)
    show("l_z_complex (diagonal, entry = m)",      Lzc)
    show("l_z_real = conj(C) @ l_z_complex @ C^T", Lzr)

    max_err = max(abs(Lzr[i][j] - Lzr[j][i].conjugate())
                  for i in range(n) for j in range(n))
    print(f"Hermiticity residual max |L - L^H| = {max_err:.3e}")
    return 0 if max_err < 1e-13 else 1


if __name__ == "__main__":
    sys.exit(main())
