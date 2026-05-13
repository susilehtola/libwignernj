// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// Build the orbital angular momentum operator l_z in the real
// spherical harmonic basis from its diagonal complex-basis form,
// using libwignernj's wignernj::real_ylm_in_complex_ylm as the change of
// basis.  C++ parallel to examples/c/real_basis_lz.c.
//
// In the complex Y_l^m basis, l_z is diagonal with entry m delta_{m m'}
// (in units of hbar).  libwignernj returns the unitary C that encodes
// the basis-vector relation
//
//     S_{l, m_r}  =  sum_{m_c}  C[m_r, m_c]  Y_l^{m_c}.
//
// The corresponding operator-matrix similarity transform is then
//
//     l_z_real  =  conj(C)  @  l_z_complex  @  transpose(C),
//
// which is the textbook  <S_{m_r}|O|S_{m_r'}>  expanded out with the
// ket using C[m_r', m_c'] and the bra using bar(C[m_r, m_c]).  The
// result is purely imaginary and antisymmetric: l_z couples each
// cosine-type real harmonic to its sin-type partner with strength
// +/- i times the diagonal m_c entry.
//
// For l = 1 this gives the textbook
//
//   l_z = (  0   0   i )
//         (  0   0   0 )
//         ( -i   0   0 )
//
// in the basis ordering m_r = -1, 0, +1.
//
// Build (out-of-tree, against an installed libwignernj):
//     c++ -std=c++11 -o real_basis_lz_cpp real_basis_lz.cpp -lwignernj -lm

#include "wignernj.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

using cdouble = std::complex<double>;

// Column-major (i, j) element of a dim x dim matrix stored as a flat
// std::vector<cdouble> with leading dimension dim.
static cdouble &at(std::vector<cdouble> &M, int dim, int i, int j)
{
    return M[static_cast<std::size_t>(j) * static_cast<std::size_t>(dim)
             + static_cast<std::size_t>(i)];
}
static cdouble at(const std::vector<cdouble> &M, int dim, int i, int j)
{
    return M[static_cast<std::size_t>(j) * static_cast<std::size_t>(dim)
             + static_cast<std::size_t>(i)];
}

// dst += A @ B for dim x dim column-major matrices.
static void cgemm_add(int dim,
                      const std::vector<cdouble> &A,
                      const std::vector<cdouble> &B,
                      std::vector<cdouble> &dst)
{
    for (int j = 0; j < dim; ++j) {
        for (int i = 0; i < dim; ++i) {
            cdouble acc(0.0, 0.0);
            for (int k = 0; k < dim; ++k)
                acc += at(A, dim, i, k) * at(B, dim, k, j);
            at(dst, dim, i, j) += acc;
        }
    }
}

static void print_matrix(const char *name, int dim,
                         const std::vector<cdouble> &M)
{
    std::printf("%s =\n", name);
    for (int i = 0; i < dim; ++i) {
        std::printf("  ");
        for (int j = 0; j < dim; ++j) {
            double re = at(M, dim, i, j).real();
            double im = at(M, dim, i, j).imag();
            if (std::fabs(re) < 1e-15) re = 0.0;
            if (std::fabs(im) < 1e-15) im = 0.0;
            std::printf("  %+.3f%+.3fi", re, im);
        }
        std::printf("\n");
    }
    std::printf("\n");
}

int main()
{
    const int l   = 1;
    const int dim = 2 * l + 1;
    const std::size_t n = static_cast<std::size_t>(dim) * dim;

    // C: basis-vector matrix, S = C Y.
    std::vector<cdouble> C = wignernj::real_ylm_in_complex_ylm<double>(l);
    // Cstar: elementwise conjugate.  CT: transpose (not conjugated).
    std::vector<cdouble> Cstar(n), CT(n), Lzc(n, cdouble(0.0, 0.0)),
                          tmp(n, cdouble(0.0, 0.0)),
                          Lzr(n, cdouble(0.0, 0.0));

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            at(Cstar, dim, i, j) = std::conj(at(C, dim, i, j));
            at(CT,    dim, i, j) = at(C, dim, j, i);
        }
    }

    for (int m = -l; m <= l; ++m)
        at(Lzc, dim, m + l, m + l) = cdouble(static_cast<double>(m), 0.0);

    cgemm_add(dim, Lzc,   CT,  tmp);
    cgemm_add(dim, Cstar, tmp, Lzr);

    std::printf("l = %d real-basis representation of l_z\n", l);
    std::printf("Basis ordering: m_r = -l, -l+1, ..., l-1, l\n");
    std::printf("(real-Y convention matches libwignernj's gaunt_real)\n\n");

    print_matrix("C (basis-vector matrix, S = C Y)",       dim, C);
    print_matrix("l_z_complex (diagonal, entry = m)",      dim, Lzc);
    print_matrix("l_z_real = conj(C) @ l_z_complex @ C^T", dim, Lzr);

    double max_err = 0.0;
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            cdouble d = at(Lzr, dim, i, j) - std::conj(at(Lzr, dim, j, i));
            max_err = std::max(max_err, std::abs(d));
        }
    }
    std::printf("Hermiticity residual max |L - L^H| = %.3e\n", max_err);
    return (max_err < 1e-13) ? 0 : 1;
}
