! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! Fortran 90 interface to libwignernj using iso_c_binding.
!
! The module provides:
!   1. Raw C-interop interfaces using 2*j integer arguments for every
!      coupling-coefficient routine (3j, 6j, 9j, Clebsch-Gordan,
!      Racah W, Fano X, Gaunt, real-Gaunt) -- same convention as the C
!      API.  The single exception is wignernj_real_ylm_in_complex_ylm,
!      which takes a plain integer orbital angular momentum l.
!   2. Real-valued convenience wrappers w3j, w6j, w9j, wcg, wracahw, wfanox,
!      wgaunt, wgaunt_real that accept double-precision real j/m arguments
!      and convert internally; plus wreal_ylm_in_complex_ylm(l, c_out)
!      filling the (2l+1) x (2l+1) real/complex Y_lm overlap matrix as a
!      complex(c_double_complex) array.
!
! Phase conventions: identical to those of the underlying C library.  The
! Wigner 3j/6j/9j symbols, Clebsch-Gordan coefficient, Racah W, and
! Fano X-coefficient are pure SU(2) algebraic objects and carry no
! spherical-harmonic phase convention;
! the Clebsch-Gordan sign uses the Condon-Shortley convention of Edmonds
! and Varshalovich.  The Gaunt, real-Gaunt and real-to-complex Y_lm
! basis-overlap routines assume the Condon-Shortley phase for Y_l^m,
! with the real spherical harmonics defined by the Wikipedia/
! Condon-Shortley construction.  See the header comment of
! include/wignernj.h for the explicit formulas.
!
! Compile with preprocessing enabled (gfortran uses .F90 extension automatically,
! or pass -cpp).  The preprocessor guard on c_long_double handles platforms
! where 'long double' has no Fortran equivalent.

module wignernj
  use iso_c_binding, only: c_int, c_float, c_double, c_float_complex, &
                            c_double_complex
#if defined(__GFORTRAN__) || defined(_CRAYFTN)
  use iso_c_binding, only: c_long_double, c_long_double_complex
#endif
#ifdef WIGNERNJ_HAVE_QUADMATH
  ! c_float128 is the gfortran/ifx extension that maps directly to C's
  ! __float128 and produces no -Wc-binding-type warning under bind(c);
  ! real128 from iso_fortran_env is the same physical kind on both
  ! compilers but is technically not a C-interoperable kind, which
  ! gfortran 16 now diagnoses on every bind(c) function returning it.
  use iso_c_binding, only: c_float128, c_float128_complex
#endif
  use iso_fortran_env, only: error_unit
  implicit none

  ! -------------------------------------------------------------------------
  ! Raw C interfaces: 2*j integer arguments
  ! -------------------------------------------------------------------------
  interface

    ! --- 3j ---
    function wigner3j_f(tj1,tj2,tj3,tm1,tm2,tm3) bind(c,name='wigner3j_f')
      import c_int, c_float
      integer(c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
      real(c_float)         :: wigner3j_f
    end function

    function wigner3j(tj1,tj2,tj3,tm1,tm2,tm3) bind(c,name='wigner3j')
      import c_int, c_double
      integer(c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
      real(c_double)        :: wigner3j
    end function

#if c_long_double > 0
    function wigner3j_l(tj1,tj2,tj3,tm1,tm2,tm3) bind(c,name='wigner3j_l')
      import c_int, c_long_double
      integer(c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
      real(c_long_double)   :: wigner3j_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function wigner3j_q(tj1,tj2,tj3,tm1,tm2,tm3) bind(c,name='wigner3j_q')
      import c_int, c_float128
      integer(c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
      real(c_float128)         :: wigner3j_q
    end function
#endif

    ! --- 6j ---
    function wigner6j_f(tj1,tj2,tj3,tj4,tj5,tj6) bind(c,name='wigner6j_f')
      import c_int, c_float
      integer(c_int), value :: tj1, tj2, tj3, tj4, tj5, tj6
      real(c_float)         :: wigner6j_f
    end function

    function wigner6j(tj1,tj2,tj3,tj4,tj5,tj6) bind(c,name='wigner6j')
      import c_int, c_double
      integer(c_int), value :: tj1, tj2, tj3, tj4, tj5, tj6
      real(c_double)        :: wigner6j
    end function

#if c_long_double > 0
    function wigner6j_l(tj1,tj2,tj3,tj4,tj5,tj6) bind(c,name='wigner6j_l')
      import c_int, c_long_double
      integer(c_int), value :: tj1, tj2, tj3, tj4, tj5, tj6
      real(c_long_double)   :: wigner6j_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function wigner6j_q(tj1,tj2,tj3,tj4,tj5,tj6) bind(c,name='wigner6j_q')
      import c_int, c_float128
      integer(c_int), value :: tj1, tj2, tj3, tj4, tj5, tj6
      real(c_float128)         :: wigner6j_q
    end function
#endif

    ! --- 9j ---
    function wigner9j_f(t11,t12,t13,t21,t22,t23,t31,t32,t33) &
        bind(c,name='wigner9j_f')
      import c_int, c_float
      integer(c_int), value :: t11,t12,t13,t21,t22,t23,t31,t32,t33
      real(c_float)         :: wigner9j_f
    end function

    function wigner9j(t11,t12,t13,t21,t22,t23,t31,t32,t33) &
        bind(c,name='wigner9j')
      import c_int, c_double
      integer(c_int), value :: t11,t12,t13,t21,t22,t23,t31,t32,t33
      real(c_double)        :: wigner9j
    end function

#if c_long_double > 0
    function wigner9j_l(t11,t12,t13,t21,t22,t23,t31,t32,t33) &
        bind(c,name='wigner9j_l')
      import c_int, c_long_double
      integer(c_int), value :: t11,t12,t13,t21,t22,t23,t31,t32,t33
      real(c_long_double)   :: wigner9j_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function wigner9j_q(t11,t12,t13,t21,t22,t23,t31,t32,t33) &
        bind(c,name='wigner9j_q')
      import c_int, c_float128
      integer(c_int), value :: t11,t12,t13,t21,t22,t23,t31,t32,t33
      real(c_float128)         :: wigner9j_q
    end function
#endif

    ! --- Clebsch-Gordan ---
    function clebsch_gordan_f(tj1,tm1,tj2,tm2,tJ,tM) &
        bind(c,name='clebsch_gordan_f')
      import c_int, c_float
      integer(c_int), value :: tj1, tm1, tj2, tm2, tJ, tM
      real(c_float)         :: clebsch_gordan_f
    end function

    function clebsch_gordan(tj1,tm1,tj2,tm2,tJ,tM) &
        bind(c,name='clebsch_gordan')
      import c_int, c_double
      integer(c_int), value :: tj1, tm1, tj2, tm2, tJ, tM
      real(c_double)        :: clebsch_gordan
    end function

#if c_long_double > 0
    function clebsch_gordan_l(tj1,tm1,tj2,tm2,tJ,tM) &
        bind(c,name='clebsch_gordan_l')
      import c_int, c_long_double
      integer(c_int), value :: tj1, tm1, tj2, tm2, tJ, tM
      real(c_long_double)   :: clebsch_gordan_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function clebsch_gordan_q(tj1,tm1,tj2,tm2,tJ,tM) &
        bind(c,name='clebsch_gordan_q')
      import c_int, c_float128
      integer(c_int), value :: tj1, tm1, tj2, tm2, tJ, tM
      real(c_float128)         :: clebsch_gordan_q
    end function
#endif

    ! --- Racah W ---
    function racah_w_f(tj1,tj2,tJ,tj3,tj12,tj23) bind(c,name='racah_w_f')
      import c_int, c_float
      integer(c_int), value :: tj1, tj2, tJ, tj3, tj12, tj23
      real(c_float)         :: racah_w_f
    end function

    function racah_w(tj1,tj2,tJ,tj3,tj12,tj23) bind(c,name='racah_w')
      import c_int, c_double
      integer(c_int), value :: tj1, tj2, tJ, tj3, tj12, tj23
      real(c_double)        :: racah_w
    end function

#if c_long_double > 0
    function racah_w_l(tj1,tj2,tJ,tj3,tj12,tj23) bind(c,name='racah_w_l')
      import c_int, c_long_double
      integer(c_int), value :: tj1, tj2, tJ, tj3, tj12, tj23
      real(c_long_double)   :: racah_w_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function racah_w_q(tj1,tj2,tJ,tj3,tj12,tj23) bind(c,name='racah_w_q')
      import c_int, c_float128
      integer(c_int), value :: tj1, tj2, tJ, tj3, tj12, tj23
      real(c_float128)         :: racah_w_q
    end function
#endif

    ! --- Fano X-coefficient ---
    function fano_x_f(tj1,tj2,tj12,tj3,tj4,tj34,tj13,tj24,tJ) &
        bind(c,name='fano_x_f')
      import c_int, c_float
      integer(c_int), value :: tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ
      real(c_float)         :: fano_x_f
    end function

    function fano_x(tj1,tj2,tj12,tj3,tj4,tj34,tj13,tj24,tJ) &
        bind(c,name='fano_x')
      import c_int, c_double
      integer(c_int), value :: tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ
      real(c_double)        :: fano_x
    end function

#if c_long_double > 0
    function fano_x_l(tj1,tj2,tj12,tj3,tj4,tj34,tj13,tj24,tJ) &
        bind(c,name='fano_x_l')
      import c_int, c_long_double
      integer(c_int), value :: tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ
      real(c_long_double)   :: fano_x_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function fano_x_q(tj1,tj2,tj12,tj3,tj4,tj34,tj13,tj24,tJ) &
        bind(c,name='fano_x_q')
      import c_int, c_float128
      integer(c_int), value :: tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ
      real(c_float128)         :: fano_x_q
    end function
#endif

    ! --- Gaunt ---
    function gaunt_f(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_f')
      import c_int, c_float
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_float)         :: gaunt_f
    end function

    function gaunt(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt')
      import c_int, c_double
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_double)        :: gaunt
    end function

#if c_long_double > 0
    function gaunt_l(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_l')
      import c_int, c_long_double
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_long_double)   :: gaunt_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function gaunt_q(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_q')
      import c_int, c_float128
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_float128)         :: gaunt_q
    end function
#endif

    ! --- Real-spherical-harmonic Gaunt ---
    function gaunt_real_f(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_real_f')
      import c_int, c_float
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_float)         :: gaunt_real_f
    end function

    function gaunt_real(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_real')
      import c_int, c_double
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_double)        :: gaunt_real
    end function

#if c_long_double > 0
    function gaunt_real_l(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_real_l')
      import c_int, c_long_double
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_long_double)   :: gaunt_real_l
    end function
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    function gaunt_real_q(tl1,tm1,tl2,tm2,tl3,tm3) bind(c,name='gaunt_real_q')
      import c_int, c_float128
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      real(c_float128)         :: gaunt_real_q
    end function
#endif

    ! --- real <-> complex spherical-harmonic basis transformation ---
    ! The C entry points fill a column-major (2l+1) x (2l+1) complex
    ! matrix into a flat buffer with leading dimension 2l+1.  Fortran
    ! column-major and the C-side column-major layout describe the
    ! same memory, so the natural Fortran array C(0:2l, 0:2l) (or
    ! equivalently C(2l+1, 2l+1)) maps directly: C(m_r+l+1, m_c+l+1)
    ! is the (m_r, m_c) entry.
    subroutine wignernj_real_ylm_in_complex_ylm_f(l, c_out) &
        bind(c,name='wignernj_real_ylm_in_complex_ylm_f')
      import c_int, c_float_complex
      integer(c_int), value :: l
      complex(c_float_complex), intent(out) :: c_out(*)
    end subroutine

    subroutine wignernj_real_ylm_in_complex_ylm(l, c_out) &
        bind(c,name='wignernj_real_ylm_in_complex_ylm')
      import c_int, c_double_complex
      integer(c_int), value :: l
      complex(c_double_complex), intent(out) :: c_out(*)
    end subroutine

#if c_long_double > 0
    subroutine wignernj_real_ylm_in_complex_ylm_l(l, c_out) &
        bind(c,name='wignernj_real_ylm_in_complex_ylm_l')
      import c_int, c_long_double_complex
      integer(c_int), value :: l
      complex(c_long_double_complex), intent(out) :: c_out(*)
    end subroutine
#endif

#ifdef WIGNERNJ_HAVE_QUADMATH
    subroutine wignernj_real_ylm_in_complex_ylm_q(l, c_out) &
        bind(c,name='wignernj_real_ylm_in_complex_ylm_q')
      import c_int, c_float128_complex
      integer(c_int), value :: l
      complex(c_float128_complex), intent(out) :: c_out(*)
    end subroutine
#endif

    ! --- per-thread cache control ---
    subroutine wignernj_warmup_to(n_max) bind(c,name='wignernj_warmup_to')
      import c_int
      integer(c_int), value :: n_max
    end subroutine

    subroutine wignernj_thread_cleanup() bind(c,name='wignernj_thread_cleanup')
    end subroutine

    function wignernj_max_factorial_arg() &
        bind(c,name='wignernj_max_factorial_arg')
      import c_int
      integer(c_int) :: wignernj_max_factorial_arg
    end function

    ! --- per-symbol max-factorial helpers ---
    function wigner3j_max_factorial(tj1,tj2,tj3,tm1,tm2,tm3) &
        bind(c,name='wigner3j_max_factorial')
      import c_int
      integer(c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
      integer(c_int) :: wigner3j_max_factorial
    end function

    function wigner6j_max_factorial(tj1,tj2,tj3,tj4,tj5,tj6) &
        bind(c,name='wigner6j_max_factorial')
      import c_int
      integer(c_int), value :: tj1, tj2, tj3, tj4, tj5, tj6
      integer(c_int) :: wigner6j_max_factorial
    end function

    function wigner9j_max_factorial(tj11,tj12,tj13,tj21,tj22,tj23, &
                                    tj31,tj32,tj33) &
        bind(c,name='wigner9j_max_factorial')
      import c_int
      integer(c_int), value :: tj11, tj12, tj13, tj21, tj22, tj23, &
                               tj31, tj32, tj33
      integer(c_int) :: wigner9j_max_factorial
    end function

    function clebsch_gordan_max_factorial(tj1,tm1,tj2,tm2,tJ,tM) &
        bind(c,name='clebsch_gordan_max_factorial')
      import c_int
      integer(c_int), value :: tj1, tm1, tj2, tm2, tJ, tM
      integer(c_int) :: clebsch_gordan_max_factorial
    end function

    function racah_w_max_factorial(tj1,tj2,tJ,tj3,tj12,tj23) &
        bind(c,name='racah_w_max_factorial')
      import c_int
      integer(c_int), value :: tj1, tj2, tJ, tj3, tj12, tj23
      integer(c_int) :: racah_w_max_factorial
    end function

    function fano_x_max_factorial(tj1,tj2,tj12,tj3,tj4,tj34, &
                                  tj13,tj24,tJ) &
        bind(c,name='fano_x_max_factorial')
      import c_int
      integer(c_int), value :: tj1, tj2, tj12, tj3, tj4, tj34, &
                               tj13, tj24, tJ
      integer(c_int) :: fano_x_max_factorial
    end function

    function gaunt_max_factorial(tl1,tm1,tl2,tm2,tl3,tm3) &
        bind(c,name='gaunt_max_factorial')
      import c_int
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      integer(c_int) :: gaunt_max_factorial
    end function

    function gaunt_real_max_factorial(tl1,tm1,tl2,tm2,tl3,tm3) &
        bind(c,name='gaunt_real_max_factorial')
      import c_int
      integer(c_int), value :: tl1, tm1, tl2, tm2, tl3, tm3
      integer(c_int) :: gaunt_real_max_factorial
    end function

  end interface

contains

  ! -------------------------------------------------------------------------
  ! Convenience wrappers: real-valued double-precision j/m arguments.
  ! Use nint(2*j) to convert; validate half-integer by checking round-trip.
  ! On error, write to stderr and return 0.
  ! -------------------------------------------------------------------------

  pure function is_half_int(v) result(ok)
    real(c_double), intent(in) :: v
    logical :: ok
    integer(c_int) :: tv
    tv = nint(2.0_c_double * v)
    ok = abs(2.0_c_double * v - tv) <= 1.0e-9_c_double
  end function is_half_int

  pure function to_tj(v) result(tv)
    real(c_double), intent(in) :: v
    integer(c_int) :: tv
    tv = nint(2.0_c_double * v)
  end function to_tj

  function w3j(j1,j2,j3,m1,m2,m3) result(val)
    real(c_double), intent(in) :: j1, j2, j3, m1, m2, m3
    real(c_double) :: val
    if (.not. (is_half_int(j1) .and. is_half_int(j2) .and. is_half_int(j3) .and. &
               is_half_int(m1) .and. is_half_int(m2) .and. is_half_int(m3))) then
      write(error_unit,'(A)') 'wigner(Fortran): w3j: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = wigner3j(to_tj(j1), to_tj(j2), to_tj(j3), to_tj(m1), to_tj(m2), to_tj(m3))
  end function w3j

  function w6j(j1,j2,j3,j4,j5,j6) result(val)
    real(c_double), intent(in) :: j1, j2, j3, j4, j5, j6
    real(c_double) :: val
    if (.not. (is_half_int(j1) .and. is_half_int(j2) .and. is_half_int(j3) .and. &
               is_half_int(j4) .and. is_half_int(j5) .and. is_half_int(j6))) then
      write(error_unit,'(A)') 'wigner(Fortran): w6j: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = wigner6j(to_tj(j1), to_tj(j2), to_tj(j3), to_tj(j4), to_tj(j5), to_tj(j6))
  end function w6j

  function w9j(j11,j12,j13,j21,j22,j23,j31,j32,j33) result(val)
    real(c_double), intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    real(c_double) :: val
    if (.not. (is_half_int(j11) .and. is_half_int(j12) .and. is_half_int(j13) .and. &
               is_half_int(j21) .and. is_half_int(j22) .and. is_half_int(j23) .and. &
               is_half_int(j31) .and. is_half_int(j32) .and. is_half_int(j33))) then
      write(error_unit,'(A)') 'wigner(Fortran): w9j: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = wigner9j(to_tj(j11),to_tj(j12),to_tj(j13), &
                   to_tj(j21),to_tj(j22),to_tj(j23), &
                   to_tj(j31),to_tj(j32),to_tj(j33))
  end function w9j

  function wcg(j1,m1,j2,m2,J,M) result(val)
    real(c_double), intent(in) :: j1, m1, j2, m2, J, M
    real(c_double) :: val
    if (.not. (is_half_int(j1) .and. is_half_int(m1) .and. is_half_int(j2) .and. &
               is_half_int(m2) .and. is_half_int(J)  .and. is_half_int(M))) then
      write(error_unit,'(A)') 'wigner(Fortran): wcg: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = clebsch_gordan(to_tj(j1), to_tj(m1), to_tj(j2), to_tj(m2), to_tj(J), to_tj(M))
  end function wcg

  function wracahw(j1,j2,J,j3,j12,j23) result(val)
    real(c_double), intent(in) :: j1, j2, J, j3, j12, j23
    real(c_double) :: val
    if (.not. (is_half_int(j1) .and. is_half_int(j2) .and. is_half_int(J)   .and. &
               is_half_int(j3) .and. is_half_int(j12) .and. is_half_int(j23))) then
      write(error_unit,'(A)') 'wigner(Fortran): wracahw: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = racah_w(to_tj(j1), to_tj(j2), to_tj(J), to_tj(j3), to_tj(j12), to_tj(j23))
  end function wracahw

  function wfanox(j1,j2,j12,j3,j4,j34,j13,j24,J) result(val)
    real(c_double), intent(in) :: j1, j2, j12, j3, j4, j34, j13, j24, J
    real(c_double) :: val
    if (.not. (is_half_int(j1) .and. is_half_int(j2) .and. is_half_int(j12) .and. &
               is_half_int(j3) .and. is_half_int(j4) .and. is_half_int(j34) .and. &
               is_half_int(j13) .and. is_half_int(j24) .and. is_half_int(J))) then
      write(error_unit,'(A)') 'wigner(Fortran): wfanox: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = fano_x(to_tj(j1),  to_tj(j2),  to_tj(j12), &
                 to_tj(j3),  to_tj(j4),  to_tj(j34), &
                 to_tj(j13), to_tj(j24), to_tj(J))
  end function wfanox

  function wgaunt(l1,m1,l2,m2,l3,m3) result(val)
    real(c_double), intent(in) :: l1, m1, l2, m2, l3, m3
    real(c_double) :: val
    if (.not. (is_half_int(l1) .and. is_half_int(m1) .and. is_half_int(l2) .and. &
               is_half_int(m2) .and. is_half_int(l3) .and. is_half_int(m3))) then
      write(error_unit,'(A)') 'wigner(Fortran): wgaunt: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = gaunt(to_tj(l1), to_tj(m1), to_tj(l2), to_tj(m2), to_tj(l3), to_tj(m3))
  end function wgaunt

  function wgaunt_real(l1,m1,l2,m2,l3,m3) result(val)
    real(c_double), intent(in) :: l1, m1, l2, m2, l3, m3
    real(c_double) :: val
    if (.not. (is_half_int(l1) .and. is_half_int(m1) .and. is_half_int(l2) .and. &
               is_half_int(m2) .and. is_half_int(l3) .and. is_half_int(m3))) then
      write(error_unit,'(A)') 'wigner(Fortran): wgaunt_real: argument is not a half-integer'
      val = 0.0_c_double; return
    end if
    val = gaunt_real(to_tj(l1), to_tj(m1), to_tj(l2), to_tj(m2), to_tj(l3), to_tj(m3))
  end function wgaunt_real

  ! -------------------------------------------------------------------------
  ! Real <-> complex spherical-harmonic basis transformation.  Fills a
  ! (2l+1) x (2l+1) column-major complex(c_double_complex) matrix C where
  !     S_{l,m_r} = sum_{m_c} C(m_r+l+1, m_c+l+1) Y_l^{m_c}
  ! using the same real-Y convention as wgaunt_real / gaunt_real.
  ! -------------------------------------------------------------------------
  subroutine wreal_ylm_in_complex_ylm(l, c_out)
    integer, intent(in) :: l
    complex(c_double_complex), intent(out) :: c_out(2*l+1, 2*l+1)
    call wignernj_real_ylm_in_complex_ylm(int(l, c_int), c_out)
  end subroutine wreal_ylm_in_complex_ylm

#ifdef WIGNERNJ_HAVE_QUADMATH
  ! -------------------------------------------------------------------------
  ! Convenience wrappers at quadruple precision.  Arguments are
  ! real(c_float128); validation reuses the double-precision is_half_int by
  ! comparing the round-trip through nint(2*j).  On error, write to stderr
  ! and return 0.
  ! -------------------------------------------------------------------------

  pure function is_half_int_q(v) result(ok)
    real(c_float128), intent(in) :: v
    logical :: ok
    integer(c_int) :: tv
    tv = nint(2.0_c_float128 * v)
    ok = abs(2.0_c_float128 * v - tv) <= 1.0e-30_c_float128
  end function is_half_int_q

  pure function to_tj_q(v) result(tv)
    real(c_float128), intent(in) :: v
    integer(c_int) :: tv
    tv = nint(2.0_c_float128 * v)
  end function to_tj_q

  function w3jq(j1,j2,j3,m1,m2,m3) result(val)
    real(c_float128), intent(in) :: j1, j2, j3, m1, m2, m3
    real(c_float128) :: val
    if (.not. (is_half_int_q(j1) .and. is_half_int_q(j2) .and. is_half_int_q(j3) .and. &
               is_half_int_q(m1) .and. is_half_int_q(m2) .and. is_half_int_q(m3))) then
      write(error_unit,'(A)') 'wigner(Fortran): w3jq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = wigner3j_q(to_tj_q(j1), to_tj_q(j2), to_tj_q(j3), &
                     to_tj_q(m1), to_tj_q(m2), to_tj_q(m3))
  end function w3jq

  function w6jq(j1,j2,j3,j4,j5,j6) result(val)
    real(c_float128), intent(in) :: j1, j2, j3, j4, j5, j6
    real(c_float128) :: val
    if (.not. (is_half_int_q(j1) .and. is_half_int_q(j2) .and. is_half_int_q(j3) .and. &
               is_half_int_q(j4) .and. is_half_int_q(j5) .and. is_half_int_q(j6))) then
      write(error_unit,'(A)') 'wigner(Fortran): w6jq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = wigner6j_q(to_tj_q(j1), to_tj_q(j2), to_tj_q(j3), &
                     to_tj_q(j4), to_tj_q(j5), to_tj_q(j6))
  end function w6jq

  function w9jq(j11,j12,j13,j21,j22,j23,j31,j32,j33) result(val)
    real(c_float128), intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    real(c_float128) :: val
    if (.not. (is_half_int_q(j11) .and. is_half_int_q(j12) .and. is_half_int_q(j13) .and. &
               is_half_int_q(j21) .and. is_half_int_q(j22) .and. is_half_int_q(j23) .and. &
               is_half_int_q(j31) .and. is_half_int_q(j32) .and. is_half_int_q(j33))) then
      write(error_unit,'(A)') 'wigner(Fortran): w9jq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = wigner9j_q(to_tj_q(j11),to_tj_q(j12),to_tj_q(j13), &
                     to_tj_q(j21),to_tj_q(j22),to_tj_q(j23), &
                     to_tj_q(j31),to_tj_q(j32),to_tj_q(j33))
  end function w9jq

  function wcgq(j1,m1,j2,m2,J,M) result(val)
    real(c_float128), intent(in) :: j1, m1, j2, m2, J, M
    real(c_float128) :: val
    if (.not. (is_half_int_q(j1) .and. is_half_int_q(m1) .and. is_half_int_q(j2) .and. &
               is_half_int_q(m2) .and. is_half_int_q(J)  .and. is_half_int_q(M))) then
      write(error_unit,'(A)') 'wigner(Fortran): wcgq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = clebsch_gordan_q(to_tj_q(j1), to_tj_q(m1), to_tj_q(j2), to_tj_q(m2), &
                           to_tj_q(J),  to_tj_q(M))
  end function wcgq

  function wracahwq(j1,j2,J,j3,j12,j23) result(val)
    real(c_float128), intent(in) :: j1, j2, J, j3, j12, j23
    real(c_float128) :: val
    if (.not. (is_half_int_q(j1) .and. is_half_int_q(j2) .and. is_half_int_q(J)   .and. &
               is_half_int_q(j3) .and. is_half_int_q(j12) .and. is_half_int_q(j23))) then
      write(error_unit,'(A)') 'wigner(Fortran): wracahwq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = racah_w_q(to_tj_q(j1), to_tj_q(j2), to_tj_q(J), &
                    to_tj_q(j3), to_tj_q(j12), to_tj_q(j23))
  end function wracahwq

  function wfanoxq(j1,j2,j12,j3,j4,j34,j13,j24,J) result(val)
    real(c_float128), intent(in) :: j1, j2, j12, j3, j4, j34, j13, j24, J
    real(c_float128) :: val
    if (.not. (is_half_int_q(j1)  .and. is_half_int_q(j2)  .and. is_half_int_q(j12) .and. &
               is_half_int_q(j3)  .and. is_half_int_q(j4)  .and. is_half_int_q(j34) .and. &
               is_half_int_q(j13) .and. is_half_int_q(j24) .and. is_half_int_q(J))) then
      write(error_unit,'(A)') 'wigner(Fortran): wfanoxq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = fano_x_q(to_tj_q(j1),  to_tj_q(j2),  to_tj_q(j12), &
                   to_tj_q(j3),  to_tj_q(j4),  to_tj_q(j34), &
                   to_tj_q(j13), to_tj_q(j24), to_tj_q(J))
  end function wfanoxq

  function wgauntq(l1,m1,l2,m2,l3,m3) result(val)
    real(c_float128), intent(in) :: l1, m1, l2, m2, l3, m3
    real(c_float128) :: val
    if (.not. (is_half_int_q(l1) .and. is_half_int_q(m1) .and. is_half_int_q(l2) .and. &
               is_half_int_q(m2) .and. is_half_int_q(l3) .and. is_half_int_q(m3))) then
      write(error_unit,'(A)') 'wigner(Fortran): wgauntq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = gaunt_q(to_tj_q(l1), to_tj_q(m1), to_tj_q(l2), to_tj_q(m2), to_tj_q(l3), to_tj_q(m3))
  end function wgauntq

  function wgaunt_realq(l1,m1,l2,m2,l3,m3) result(val)
    real(c_float128), intent(in) :: l1, m1, l2, m2, l3, m3
    real(c_float128) :: val
    if (.not. (is_half_int_q(l1) .and. is_half_int_q(m1) .and. is_half_int_q(l2) .and. &
               is_half_int_q(m2) .and. is_half_int_q(l3) .and. is_half_int_q(m3))) then
      write(error_unit,'(A)') 'wigner(Fortran): wgaunt_realq: argument is not a half-integer'
      val = 0.0_c_float128; return
    end if
    val = gaunt_real_q(to_tj_q(l1), to_tj_q(m1), to_tj_q(l2), to_tj_q(m2), &
                       to_tj_q(l3), to_tj_q(m3))
  end function wgaunt_realq

  subroutine wreal_ylm_in_complex_ylmq(l, c_out)
    integer, intent(in) :: l
    complex(c_float128_complex), intent(out) :: c_out(2*l+1, 2*l+1)
    call wignernj_real_ylm_in_complex_ylm_q(int(l, c_int), c_out)
  end subroutine wreal_ylm_in_complex_ylmq
#endif

end module wignernj
