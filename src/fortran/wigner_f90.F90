! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! Fortran 90 interface to libwignernj using iso_c_binding.
!
! The module provides:
!   1. Raw C-interop interfaces using 2*j integer arguments (same convention
!      as the C API).
!   2. Real-valued convenience wrappers w3j, w6j, w9j, wcg, wracahw, wgaunt
!      that accept double-precision real j/m arguments and convert internally.
!
! Compile with preprocessing enabled (gfortran uses .F90 extension automatically,
! or pass -cpp).  The preprocessor guard on c_long_double handles platforms
! where 'long double' has no Fortran equivalent.

module wigner
  use iso_c_binding, only: c_int, c_float, c_double
#if defined(__GFORTRAN__) || defined(_CRAYFTN)
  use iso_c_binding, only: c_long_double
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

end module wigner
