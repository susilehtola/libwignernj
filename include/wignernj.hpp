// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// C++11 header-only wrapper around the libwignernj C API.
//
// "Header-only" means the C++ side has no separate .cpp -- every template and
// overload below is inline.  It does NOT mean stand-alone: each specialisation
// forwards to a C function (wigner3j, wigner6j_l, ...) defined in libwignernj,
// so you still need to link -lwignernj -lm (and -lmpfr if MPFR support was
// compiled in) for the resulting program to link.
//
// Usage:
//   #include "wignernj.hpp"
//   double v = wignernj::symbol3j<double>(2, 2, 0, 0, 0, 0);
//   double v = wignernj::symbol3j(1.0, 1.0, 0.0, 0.0, 0.0, 0.0); // real-valued
//   float  v = wignernj::symbol3j<float>(0.5, 0.5, 1.0, 0.5, -0.5, 0.0);
//
// Throws std::invalid_argument if a real-valued argument is not a half-integer.
// Returns 0 silently for symbols that vanish by selection rules.
//
// Phase conventions: identical to those of the underlying C library.  3j, 6j,
// 9j, Clebsch-Gordan, and Racah W are pure SU(2) algebraic objects and carry
// no spherical-harmonic phase convention; the Clebsch-Gordan coefficient uses
// the Condon-Shortley sign of Edmonds and Varshalovich.  The Gaunt and real-
// Gaunt routines assume the Condon-Shortley phase for Y_l^m, with the real
// spherical harmonics defined by the Wikipedia/Condon-Shortley construction.
// See include/wignernj.h for the explicit formulas.

#ifndef WIGNERNJ_HPP
#define WIGNERNJ_HPP

#include "wignernj.h"
#include <stdexcept>
#include <cmath>
#include <string>

namespace wignernj {

// ── integer-argument API (primary templates, explicit specialisations) ─────

template<typename T> T symbol3j(int,int,int,int,int,int);
template<typename T> T symbol6j(int,int,int,int,int,int);
template<typename T> T symbol9j(int,int,int,int,int,int,int,int,int);
template<typename T> T cg      (int,int,int,int,int,int);
template<typename T> T racahw  (int,int,int,int,int,int);
template<typename T> T fanox   (int,int,int,int,int,int,int,int,int);
template<typename T> T gaunt   (int,int,int,int,int,int);
template<typename T> T gauntreal(int,int,int,int,int,int);

template<> inline float
symbol3j<float>(int a,int b,int c,int d,int e,int f)
{ return wigner3j_f(a,b,c,d,e,f); }
template<> inline double
symbol3j<double>(int a,int b,int c,int d,int e,int f)
{ return wigner3j(a,b,c,d,e,f); }
template<> inline long double
symbol3j<long double>(int a,int b,int c,int d,int e,int f)
{ return wigner3j_l(a,b,c,d,e,f); }

template<> inline float
symbol6j<float>(int a,int b,int c,int d,int e,int f)
{ return wigner6j_f(a,b,c,d,e,f); }
template<> inline double
symbol6j<double>(int a,int b,int c,int d,int e,int f)
{ return wigner6j(a,b,c,d,e,f); }
template<> inline long double
symbol6j<long double>(int a,int b,int c,int d,int e,int f)
{ return wigner6j_l(a,b,c,d,e,f); }

template<> inline float
symbol9j<float>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return wigner9j_f(a,b,c,d,e,f,g,h,i); }
template<> inline double
symbol9j<double>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return wigner9j(a,b,c,d,e,f,g,h,i); }
template<> inline long double
symbol9j<long double>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return wigner9j_l(a,b,c,d,e,f,g,h,i); }

template<> inline float
cg<float>(int a,int b,int c,int d,int e,int f)
{ return clebsch_gordan_f(a,b,c,d,e,f); }
template<> inline double
cg<double>(int a,int b,int c,int d,int e,int f)
{ return clebsch_gordan(a,b,c,d,e,f); }
template<> inline long double
cg<long double>(int a,int b,int c,int d,int e,int f)
{ return clebsch_gordan_l(a,b,c,d,e,f); }

template<> inline float
racahw<float>(int a,int b,int c,int d,int e,int f)
{ return racah_w_f(a,b,c,d,e,f); }
template<> inline double
racahw<double>(int a,int b,int c,int d,int e,int f)
{ return racah_w(a,b,c,d,e,f); }
template<> inline long double
racahw<long double>(int a,int b,int c,int d,int e,int f)
{ return racah_w_l(a,b,c,d,e,f); }

template<> inline float
fanox<float>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return fano_x_f(a,b,c,d,e,f,g,h,i); }
template<> inline double
fanox<double>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return fano_x(a,b,c,d,e,f,g,h,i); }
template<> inline long double
fanox<long double>(int a,int b,int c,int d,int e,int f,int g,int h,int i)
{ return fano_x_l(a,b,c,d,e,f,g,h,i); }

template<> inline float
gaunt<float>(int a,int b,int c,int d,int e,int f)
{ return gaunt_f(a,b,c,d,e,f); }
template<> inline double
gaunt<double>(int a,int b,int c,int d,int e,int f)
{ return ::gaunt(a,b,c,d,e,f); }
template<> inline long double
gaunt<long double>(int a,int b,int c,int d,int e,int f)
{ return gaunt_l(a,b,c,d,e,f); }

template<> inline float
gauntreal<float>(int a,int b,int c,int d,int e,int f)
{ return gaunt_real_f(a,b,c,d,e,f); }
template<> inline double
gauntreal<double>(int a,int b,int c,int d,int e,int f)
{ return ::gaunt_real(a,b,c,d,e,f); }
template<> inline long double
gauntreal<long double>(int a,int b,int c,int d,int e,int f)
{ return gaunt_real_l(a,b,c,d,e,f); }

// ── convenience helpers ───────────────────────────────────────────────────

namespace detail {

inline int to_two_j(double v, const char *name)
{
    int tv = static_cast<int>(2.0 * v + (v >= 0 ? 0.5 : -0.5));
    if (std::abs(2.0 * v - tv) > 1e-9)
        throw std::invalid_argument(
            std::string("wigner: ") + name + " is not a half-integer");
    return tv;
}

} // namespace detail

// ── real-valued overloads (accept integer or half-integer doubles) ─────────

template<typename T = double>
inline T symbol3j(double j1, double j2, double j3,
                  double m1, double m2, double m3)
{
    return symbol3j<T>(detail::to_two_j(j1,"j1"), detail::to_two_j(j2,"j2"),
                       detail::to_two_j(j3,"j3"), detail::to_two_j(m1,"m1"),
                       detail::to_two_j(m2,"m2"), detail::to_two_j(m3,"m3"));
}

template<typename T = double>
inline T symbol6j(double j1, double j2, double j3,
                  double j4, double j5, double j6)
{
    return symbol6j<T>(detail::to_two_j(j1,"j1"), detail::to_two_j(j2,"j2"),
                       detail::to_two_j(j3,"j3"), detail::to_two_j(j4,"j4"),
                       detail::to_two_j(j5,"j5"), detail::to_two_j(j6,"j6"));
}

template<typename T = double>
inline T symbol9j(double j11, double j12, double j13,
                  double j21, double j22, double j23,
                  double j31, double j32, double j33)
{
    return symbol9j<T>(
        detail::to_two_j(j11,"j11"), detail::to_two_j(j12,"j12"),
        detail::to_two_j(j13,"j13"), detail::to_two_j(j21,"j21"),
        detail::to_two_j(j22,"j22"), detail::to_two_j(j23,"j23"),
        detail::to_two_j(j31,"j31"), detail::to_two_j(j32,"j32"),
        detail::to_two_j(j33,"j33"));
}

template<typename T = double>
inline T cg(double j1, double m1, double j2, double m2, double J, double M)
{
    return cg<T>(detail::to_two_j(j1,"j1"), detail::to_two_j(m1,"m1"),
                 detail::to_two_j(j2,"j2"), detail::to_two_j(m2,"m2"),
                 detail::to_two_j(J,"J"),   detail::to_two_j(M,"M"));
}

template<typename T = double>
inline T racahw(double j1, double j2, double J, double j3,
                double j12, double j23)
{
    return racahw<T>(detail::to_two_j(j1,"j1"),   detail::to_two_j(j2,"j2"),
                     detail::to_two_j(J,"J"),     detail::to_two_j(j3,"j3"),
                     detail::to_two_j(j12,"j12"), detail::to_two_j(j23,"j23"));
}

template<typename T = double>
inline T fanox(double j1, double j2, double j12,
               double j3, double j4, double j34,
               double j13, double j24, double J)
{
    return fanox<T>(detail::to_two_j(j1,"j1"),   detail::to_two_j(j2,"j2"),
                    detail::to_two_j(j12,"j12"), detail::to_two_j(j3,"j3"),
                    detail::to_two_j(j4,"j4"),   detail::to_two_j(j34,"j34"),
                    detail::to_two_j(j13,"j13"), detail::to_two_j(j24,"j24"),
                    detail::to_two_j(J,"J"));
}

template<typename T = double>
inline T gaunt(double l1, double m1, double l2, double m2,
               double l3, double m3)
{
    return gaunt<T>(detail::to_two_j(l1,"l1"), detail::to_two_j(m1,"m1"),
                    detail::to_two_j(l2,"l2"), detail::to_two_j(m2,"m2"),
                    detail::to_two_j(l3,"l3"), detail::to_two_j(m3,"m3"));
}

template<typename T = double>
inline T gauntreal(double l1, double m1, double l2, double m2,
                   double l3, double m3)
{
    return gauntreal<T>(detail::to_two_j(l1,"l1"), detail::to_two_j(m1,"m1"),
                        detail::to_two_j(l2,"l2"), detail::to_two_j(m2,"m2"),
                        detail::to_two_j(l3,"l3"), detail::to_two_j(m3,"m3"));
}

} // namespace wignernj

#endif /* WIGNERNJ_HPP */
