! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! libwignernj Fortran API demonstration.
!
! Calls every public symbol family exposed by the wignernj module once
! with small, textbook-scale arguments and prints the result alongside
! an analytic reference value.  The program exits with a non-zero
! status if any computed value disagrees with its reference by more
! than 1e-14, otherwise zero.
!
! The Fortran wrappers w3j, w6j, w9j, wcg, wracahw, wfanox, wgaunt,
! wgaunt_real accept j as a real(c_double) (integer or half-integer),
! mirroring Edmonds 1957 / Varshalovich et al. 1988 conventions.  The
! library doubles each argument internally before forwarding to the C
! API, which is why all arguments are passed as ordinary reals.
!
! Build (out-of-tree, against an installed libwignernj):
!     gfortran -o all_symbols_f all_symbols.f90 \
!         -I /path/to/installed/include \
!         -L /path/to/installed/lib -lwignernj_f03 -lwignernj
program all_symbols
  use, intrinsic :: iso_c_binding, only: c_double
  use, intrinsic :: iso_fortran_env, only: error_unit
  use wignernj
  implicit none

  real(c_double), parameter :: pi  = 3.14159265358979323846264338327950288_c_double
  real(c_double), parameter :: tol = 1.0e-14_c_double
  integer :: failed
  real(c_double) :: v, ref

  failed = 0

  print '(A)', 'libwignernj Fortran API demonstration -- one call per symbol family'
  print '(A)', '-------------------------------------------------------------------'

  ! 1. Wigner 3j symbol.   ( 1 1 0 )
  !                        ( 0 0 0 )  =  -1/sqrt(3)
  v   = w3j(1.0_c_double, 1.0_c_double, 0.0_c_double, &
            0.0_c_double, 0.0_c_double, 0.0_c_double)
  ref = -1.0_c_double / sqrt(3.0_c_double)
  call check('w3j(1,1,0; 0,0,0)', v, ref, tol, failed)

  ! 2. Wigner 6j symbol.   { 1 1 1 }
  !                        { 1 1 1 }  =  1/6
  v   = w6j(1.0_c_double, 1.0_c_double, 1.0_c_double, &
            1.0_c_double, 1.0_c_double, 1.0_c_double)
  ref = 1.0_c_double / 6.0_c_double
  call check('w6j{1,1,1; 1,1,1}', v, ref, tol, failed)

  ! 3. Wigner 9j symbol.   { 1 1 0 }
  !                        { 1 1 0 }  =  1/3
  !                        { 0 0 0 }
  v   = w9j(1.0_c_double, 1.0_c_double, 0.0_c_double, &
            1.0_c_double, 1.0_c_double, 0.0_c_double, &
            0.0_c_double, 0.0_c_double, 0.0_c_double)
  ref = 1.0_c_double / 3.0_c_double
  call check('w9j{1,1,0; 1,1,0; 0,0,0}', v, ref, tol, failed)

  ! 4. Clebsch-Gordan coefficient.   <1,0; 1,0 | 2,0> = sqrt(2/3)
  v   = wcg(1.0_c_double, 0.0_c_double, &
            1.0_c_double, 0.0_c_double, &
            2.0_c_double, 0.0_c_double)
  ref = sqrt(2.0_c_double / 3.0_c_double)
  call check('wcg(1,0; 1,0 | 2,0)', v, ref, tol, failed)

  ! 5. Racah W coefficient.   W(1,1,1,1; 1,1) = 1/6
  v   = wracahw(1.0_c_double, 1.0_c_double, 1.0_c_double, &
                1.0_c_double, 1.0_c_double, 1.0_c_double)
  ref = 1.0_c_double / 6.0_c_double
  call check('wracahw(1,1,1,1; 1,1)', v, ref, tol, failed)

  ! 6. Fano X coefficient.   X(1,1,1; 1,1,1; 1,1,2) = 1/2
  v   = wfanox(1.0_c_double, 1.0_c_double, 1.0_c_double, &
               1.0_c_double, 1.0_c_double, 1.0_c_double, &
               1.0_c_double, 1.0_c_double, 2.0_c_double)
  ref = 0.5_c_double
  call check('wfanox(1,1,1; 1,1,1; 1,1,2)', v, ref, tol, failed)

  ! 7. Gaunt coefficient over complex Y_l^m.
  !    integral(Y_1^0 Y_1^0 Y_2^0) dOmega = 1/sqrt(5*pi)
  v   = wgaunt(1.0_c_double, 0.0_c_double, &
               1.0_c_double, 0.0_c_double, &
               2.0_c_double, 0.0_c_double)
  ref = 1.0_c_double / sqrt(5.0_c_double * pi)
  call check('wgaunt(1,0; 1,0; 2,0)', v, ref, tol, failed)

  ! 8. Gaunt over real spherical harmonics (m=0 -> same as complex).
  v   = wgaunt_real(1.0_c_double, 0.0_c_double, &
                    1.0_c_double, 0.0_c_double, &
                    2.0_c_double, 0.0_c_double)
  ref = 1.0_c_double / sqrt(5.0_c_double * pi)
  call check('wgaunt_real(1,0; 1,0; 2,0)', v, ref, tol, failed)

  if (failed == 0) then
     print '(A)', ''
     print '(A)', 'All symbols agree with their analytic references.'
     stop 0
  else
     stop 1
  end if

contains

  subroutine check(label, computed, expected, tol, failed)
    character(*),  intent(in)    :: label
    real(c_double), intent(in)    :: computed, expected, tol
    integer,        intent(inout) :: failed
    print '("  ", A, T42, " = ", SP, F19.15, "   (expected ", SP, F19.15, ")")', &
        label, computed, expected
    if (abs(computed - expected) > tol) then
       write(error_unit, '("  FAIL: |", G0, " - ", G0, "| = ", G0, &
                          &" exceeds tolerance ", G0)') &
           computed, expected, abs(computed - expected), tol
       failed = failed + 1
    end if
  end subroutine check

end program all_symbols
