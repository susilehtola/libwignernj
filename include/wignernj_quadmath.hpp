// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// C++11 quadmath (__float128) extension for the libwignernj header-only
// wrapper.  Brings the C++ wrapper to parity with wignernj_quadmath.h:
// the same wignernj::symbol3j<T>() / wignernj::cg<T>() / ... template
// surface, now usable with T = __float128.
//
// Citation: if libwignernj contributes to published work, please cite
//   S. Lehtola, "libwignernj: a reusable C/C++/Fortran/Python library
//   for exact Wigner symbols and related coefficients", arXiv:2605.06634
//   (2026), doi:10.48550/arXiv.2605.06634.
//
// Usage:
//   #include "wignernj_quadmath.hpp"   // pulls in wignernj.hpp + wignernj_quadmath.h
//   __float128 v = wignernj::symbol3j<__float128>(2, 2, 0, 0, 0, 0);
//
// Requires libwignernj built with -DBUILD_QUADMATH=ON so that the
// wigner3j_q / wigner6j_q / ... C entry points this header forwards to
// are exported.  Link -lwignernj -lquadmath -lm.
//
// std::complex<__float128> is a libstdc++ extension and not portable
// across standard libraries (libc++ does not provide it); the
// real-to-complex Y_lm basis-transformation overload therefore takes
// wignernj_cfloat128_t * (a layout-compatible {__float128 _pair[2]}
// struct in C++ mode) rather than std::complex<__float128> *.

#ifndef WIGNERNJ_QUADMATH_HPP
#define WIGNERNJ_QUADMATH_HPP

#include "wignernj.hpp"
#include "wignernj_quadmath.h"

#include <stdexcept>
#include <vector>

namespace wignernj {

// ── integer-argument template specialisations at __float128 ───────────────

template<> inline __float128
symbol3j<__float128>(int a,int b,int c,int d,int e,int f)
{ return ::wigner3j_q(a,b,c,d,e,f); }

template<> inline __float128
symbol6j<__float128>(int a,int b,int c,int d,int e,int f)
{ return ::wigner6j_q(a,b,c,d,e,f); }

template<> inline __float128
symbol9j<__float128>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return ::wigner9j_q(a,b,c,d,e,f,g,h,i); }

template<> inline __float128
cg<__float128>(int a,int b,int c,int d,int e,int f)
{ return ::clebsch_gordan_q(a,b,c,d,e,f); }

template<> inline __float128
racahw<__float128>(int a,int b,int c,int d,int e,int f)
{ return ::racah_w_q(a,b,c,d,e,f); }

template<> inline __float128
fanox<__float128>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return ::fano_x_q(a,b,c,d,e,f,g,h,i); }

template<> inline __float128
gaunt<__float128>(int a,int b,int c,int d,int e,int f)
{ return ::gaunt_q(a,b,c,d,e,f); }

template<> inline __float128
gauntreal<__float128>(int a,int b,int c,int d,int e,int f)
{ return ::gaunt_real_q(a,b,c,d,e,f); }

// Once the integer-API specialisations above are visible, the
// real-valued (double-argument) convenience overloads in wignernj.hpp
// pick up __float128 via their `template<typename T = double>` form
// with no further work, e.g.
//     wignernj::symbol3j<__float128>(0.5, 0.5, 0.0, 0.5, -0.5, 0.0);

// ── real <-> complex Y_lm basis matrix at __float128 ──────────────────────
//
// Non-template overload of wignernj::real_ylm_in_complex_ylm taking the
// C-side wignernj_cfloat128_t * buffer (std::complex<__float128> is not
// portably available).  Selected by ordinary overload resolution; no
// explicit template argument needed at the call site:
//     wignernj_cfloat128_t C[(2*l+1)*(2*l+1)];
//     wignernj::real_ylm_in_complex_ylm(l, C);

inline void real_ylm_in_complex_ylm(int l, wignernj_cfloat128_t *C_out)
{ ::wignernj_real_ylm_in_complex_ylm_q(l, C_out); }

// Vector-returning convenience form, parallel to the float/double/
// long-double overloads in wignernj.hpp.  Named with the `_q` suffix
// because overloading purely on return type is not permitted in C++,
// and the float/double/long-double versions overload on the bare name
// via the std::complex<T> template parameter, which __float128 cannot
// share.
inline std::vector<wignernj_cfloat128_t> real_ylm_in_complex_ylm_q(int l)
{
    if (l < 0)
        throw std::invalid_argument(
            "wignernj::real_ylm_in_complex_ylm_q: l must be >= 0");
    const std::size_t dim = static_cast<std::size_t>(2 * l + 1);
    std::vector<wignernj_cfloat128_t> C(dim * dim);
    ::wignernj_real_ylm_in_complex_ylm_q(l, C.data());
    return C;
}

} // namespace wignernj

#endif /* WIGNERNJ_QUADMATH_HPP */
