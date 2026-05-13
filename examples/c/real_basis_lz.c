/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Build the orbital angular momentum operator l_z in the real
 * spherical harmonic basis from its diagonal complex-basis form,
 * using libwignernj's wignernj_real_ylm_in_complex_ylm as the change of
 * basis.
 *
 * In the complex Y_l^m basis, l_z is diagonal with entry m delta_{m m'}
 * (in units of hbar).  libwignernj returns the unitary C that encodes
 * the basis-vector relation
 *
 *     S_{l, m_r}  =  sum_{m_c}  C[m_r, m_c]  Y_l^{m_c}.
 *
 * The corresponding operator-matrix similarity transform is then
 *
 *     l_z_real  =  conj(C)  @  l_z_complex  @  transpose(C),
 *
 * which is the textbook  <S_{m_r}|O|S_{m_r'}>  expanded out with the
 * ket using C[m_r', m_c'] and the bra using bar(C[m_r, m_c]).  The
 * result is purely imaginary and antisymmetric: l_z couples each
 * cosine-type real harmonic to its sin-type partner with strength
 * +/- i times the diagonal m_c entry.
 *
 * For l = 1 this gives the textbook
 *
 *   l_z = (  0   0   i )
 *         (  0   0   0 )
 *         ( -i   0   0 )
 *
 * in the basis ordering m_r = -1, 0, +1 (Wikipedia / Condon-Shortley
 * real-Y convention; sin partner first, then m=0, then cos partner).
 */
#include "wignernj.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Column-major (re, im)-interleaved buffer. */
#define RE(buf, dim, r, c) ((buf)[2u * ((size_t)(c) * (size_t)(dim) + (size_t)(r))    ])
#define IM(buf, dim, r, c) ((buf)[2u * ((size_t)(c) * (size_t)(dim) + (size_t)(r)) + 1u])

/* dst[i, j] += sum_k A[i, k] * B[k, j], complex.  All matrices column-major,
 * dimension dim x dim. */
static void cgemm_add(int dim, const double *A, const double *B, double *dst)
{
    for (int j = 0; j < dim; j++) {
        for (int i = 0; i < dim; i++) {
            double acc_re = 0.0, acc_im = 0.0;
            for (int k = 0; k < dim; k++) {
                double a_re = RE(A, dim, i, k);
                double a_im = IM(A, dim, i, k);
                double b_re = RE(B, dim, k, j);
                double b_im = IM(B, dim, k, j);
                acc_re += a_re * b_re - a_im * b_im;
                acc_im += a_re * b_im + a_im * b_re;
            }
            RE(dst, dim, i, j) += acc_re;
            IM(dst, dim, i, j) += acc_im;
        }
    }
}

static void print_matrix(const char *name, int dim, const double *M)
{
    printf("%s =\n", name);
    for (int i = 0; i < dim; i++) {
        printf("  ");
        for (int j = 0; j < dim; j++) {
            double re = RE(M, dim, i, j);
            double im = IM(M, dim, i, j);
            /* Clean -0.0 to +0.0 for legible output. */
            if (fabs(re) < 1e-15) re = 0.0;
            if (fabs(im) < 1e-15) im = 0.0;
            printf("  % .3f%+.3fi", re, im);
        }
        printf("\n");
    }
    printf("\n");
}

int main(void)
{
    /* Demo for l = 1.  Generalises to any l by changing this constant. */
    const int l = 1;
    const int dim = 2 * l + 1;
    const size_t nflat = (size_t)dim * (size_t)dim;
    const size_t nscalar = 2u * nflat;

    /* C: basis-vector matrix, S = C Y. */
    double *C = (double *)calloc(nscalar, sizeof(double));
    /* C_star: elementwise complex conjugate of C. */
    double *Cstar = (double *)calloc(nscalar, sizeof(double));
    /* C_T: transpose of C (NOT conjugated). */
    double *CT = (double *)calloc(nscalar, sizeof(double));
    /* l_z_complex: diagonal with entry m delta_{m m'}; m index = m + l. */
    double *Lzc = (double *)calloc(nscalar, sizeof(double));
    /* Workspace: l_z_complex @ transpose(C). */
    double *tmp = (double *)calloc(nscalar, sizeof(double));
    /* Output: l_z in real basis. */
    double *Lzr = (double *)calloc(nscalar, sizeof(double));
    if (!C || !Cstar || !CT || !Lzc || !tmp || !Lzr) {
        fprintf(stderr, "alloc failed\n");
        return 1;
    }

    /* Fill C.  The library entry point takes a typed `wignernj_cdouble_t *`;
     * we keep the local buffer as a flat (re, im)-interleaved `double *` and
     * cast at the call, well-defined because `wignernj_cdouble_t` is layout-
     * compatible with `double[2]` on every supported toolchain. */
    wignernj_real_ylm_in_complex_ylm(l, (wignernj_cdouble_t *)C);

    /* Build C_star (elementwise conjugate) and C_T (transpose). */
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            RE(Cstar, dim, i, j) =  RE(C, dim, i, j);
            IM(Cstar, dim, i, j) = -IM(C, dim, i, j);
            RE(CT,    dim, i, j) =  RE(C, dim, j, i);
            IM(CT,    dim, i, j) =  IM(C, dim, j, i);
        }
    }

    /* Build the diagonal complex-basis l_z. */
    for (int m = -l; m <= l; m++) {
        RE(Lzc, dim, m + l, m + l) = (double)m;
    }

    /* tmp = l_z_complex @ transpose(C) */
    cgemm_add(dim, Lzc, CT, tmp);
    /* Lzr = conj(C) @ tmp */
    cgemm_add(dim, Cstar, tmp, Lzr);

    printf("l = %d real-basis representation of l_z\n", l);
    printf("Basis ordering: m_r = -l, -l+1, ..., l-1, l\n");
    printf("(real-Y convention matches libwignernj's gaunt_real)\n\n");

    print_matrix("C (basis-vector matrix, S = C Y)",         dim, C);
    print_matrix("l_z_complex (diagonal, entry = m)",        dim, Lzc);
    print_matrix("l_z_real = conj(C) @ l_z_complex @ C^T",   dim, Lzr);

    /* Sanity: the result must be Hermitian (and in fact purely imaginary
     * and antisymmetric for orbital l_z). */
    double max_err = 0.0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            double d_re = RE(Lzr, dim, i, j) - RE(Lzr, dim, j, i);
            double d_im = IM(Lzr, dim, i, j) + IM(Lzr, dim, j, i);
            if (fabs(d_re) > max_err) max_err = fabs(d_re);
            if (fabs(d_im) > max_err) max_err = fabs(d_im);
        }
    }
    printf("Hermiticity residual max |L - L^H| = %.3e\n", max_err);

    free(C); free(Cstar); free(CT); free(Lzc); free(tmp); free(Lzr);
    return (max_err < 1e-13) ? 0 : 1;
}
