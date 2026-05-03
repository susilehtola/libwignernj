! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! Fortran test program for the wigner module's quadruple-precision
! convenience wrappers (real(real128)).  Built only when the C library
! is configured with -DBUILD_QUADMATH=ON.
!
! Compile with -cpp; .F90 extension already enables the preprocessor in
! gfortran by default.

program test_wigner_quadmath
#ifdef WIGNERNJ_HAVE_QUADMATH
  use wigner
  use iso_fortran_env, only: real128
  implicit none

  integer        :: npass, nfail, total
  real(real128)  :: val, exp_val
  ! 4 ulp at quad precision (FLT128_EPSILON ~= 1.93e-34)
  real(real128), parameter :: q_tol = 4.0_real128 * 1.92592994438723585305597794258493e-34_real128

  npass = 0
  nfail = 0

  ! 3j(1,1,0; 0,0,0) = -1/sqrt(3)
  val = w3jq(1.0_real128, 1.0_real128, 0.0_real128, &
             0.0_real128, 0.0_real128, 0.0_real128)
  exp_val = -1.0_real128 / sqrt(3.0_real128)
  call check_q(val, exp_val, q_tol, 'w3jq(1,1,0;0,0,0)', npass, nfail)

  ! 6j{1,1,1;1,1,1} = 1/6
  val = w6jq(1.0_real128, 1.0_real128, 1.0_real128, &
             1.0_real128, 1.0_real128, 1.0_real128)
  exp_val = 1.0_real128 / 6.0_real128
  call check_q(val, exp_val, q_tol, 'w6jq(1,1,1;1,1,1)', npass, nfail)

  ! Selection-rule zero
  val = w3jq(1.0_real128, 1.0_real128, 1.0_real128, &
             0.0_real128, 0.0_real128, 0.0_real128)
  call check_q_abs(val, 0.0_real128, q_tol, 'w3jq parity zero', npass, nfail)

  ! Large j: 3j(50 50 0; 0 0 0) = 1/sqrt(101) -- exercises multi-word bigint
  val = w3jq(50.0_real128, 50.0_real128, 0.0_real128, &
             0.0_real128, 0.0_real128, 0.0_real128)
  exp_val = 1.0_real128 / sqrt(101.0_real128)
  call check_q(val, exp_val, q_tol, 'w3jq large-j', npass, nfail)

  ! CG <1/2 1/2; 1/2 -1/2 | 1 0> = 1/sqrt(2)
  val = wcgq(0.5_real128,  0.5_real128, &
             0.5_real128, -0.5_real128, &
             1.0_real128,  0.0_real128)
  exp_val = 1.0_real128 / sqrt(2.0_real128)
  call check_q(val, exp_val, q_tol, 'wcgq(1/2,1/2;1/2,-1/2|1,0)', npass, nfail)

  ! Raw integer API: wigner3j_q(2,2,0, 0,0,0) = -1/sqrt(3)
  val = wigner3j_q(2, 2, 0, 0, 0, 0)
  exp_val = -1.0_real128 / sqrt(3.0_real128)
  call check_q(val, exp_val, q_tol, 'wigner3j_q raw', npass, nfail)

  ! Cross-precision: real128 must agree with double to double precision
  block
    real(real128) :: vq
    real(8)       :: vd
    vq = w9jq(0.5_real128, 0.5_real128, 1.0_real128, &
              0.5_real128, 0.5_real128, 0.0_real128, &
              1.0_real128, 1.0_real128, 1.0_real128)
    vd = w9j(0.5d0, 0.5d0, 1.0d0,  0.5d0, 0.5d0, 0.0d0,  1.0d0, 1.0d0, 1.0d0)
    call check_q(vq, real(vd, kind=real128), 1.0e-14_real128, &
                 'w9jq agrees with w9j (double)', npass, nfail)
  end block

  ! Racah W: W(j1,j2,J,j3;j12,j23) = (-1)^(j1+j2+J+j3) * 6j{j1,j2,j12;j3,J,j23}
  ! W(1,1,0,1;0,1) phase = (-1)^3 = -1
  block
    real(real128) :: vq, vw
    vq = wracahwq(1.0_real128, 1.0_real128, 0.0_real128, &
                  1.0_real128, 0.0_real128, 1.0_real128)
    vw = -w6jq(1.0_real128, 1.0_real128, 0.0_real128, &
               1.0_real128, 0.0_real128, 1.0_real128)
    call check_q(vq, vw, q_tol, 'wracahwq vs -w6jq', npass, nfail)
  end block

  ! Fano X = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * 9j.
  ! At all-equal-j=2 the norm is 5^2 = 25.
  block
    real(real128) :: vq, vw
    vq = wfanoxq(2.0_real128, 2.0_real128, 2.0_real128, &
                 2.0_real128, 2.0_real128, 2.0_real128, &
                 2.0_real128, 2.0_real128, 2.0_real128)
    vw = 25.0_real128 * w9jq(2.0_real128, 2.0_real128, 2.0_real128, &
                             2.0_real128, 2.0_real128, 2.0_real128, &
                             2.0_real128, 2.0_real128, 2.0_real128)
    call check_q(vq, vw, q_tol, 'wfanoxq = 25 * w9jq', npass, nfail)
  end block

  ! Gaunt: complex-spherical-harmonic Gaunt at l=0,m=0 → 1/(2*sqrt(pi))
  block
    real(real128) :: vq, expected
    real(real128), parameter :: piq = 3.14159265358979323846264338327950288_real128
    vq = wgauntq(0.0_real128, 0.0_real128, 0.0_real128, &
                 0.0_real128, 0.0_real128, 0.0_real128)
    expected = 0.5_real128 / sqrt(piq)
    call check_q(vq, expected, q_tol*8, 'wgauntq(0,0,0,0,0,0)', npass, nfail)
  end block

  ! Real-spherical-harmonic Gaunt: at all-m=0 it equals the complex Gaunt.
  block
    real(real128) :: vqr, vqc
    vqr = wgaunt_realq(1.0_real128, 0.0_real128, &
                       1.0_real128, 0.0_real128, &
                       2.0_real128, 0.0_real128)
    vqc = wgauntq(1.0_real128, 0.0_real128, &
                  1.0_real128, 0.0_real128, &
                  2.0_real128, 0.0_real128)
    call check_q(vqr, vqc, q_tol*8, 'wgaunt_realq m=0 = wgauntq', npass, nfail)
  end block

  total = npass + nfail
  print '(I0,"/",I0," passed")', npass, total
  if (nfail > 0) then
    error stop 1
  end if

contains

  subroutine check_q(actual, expected, tol, label, npass, nfail)
    real(real128),    intent(in)    :: actual, expected, tol
    character(len=*), intent(in)    :: label
    integer,          intent(inout) :: npass, nfail
    real(real128) :: scale, diff
    scale = max(abs(expected), 1.0_real128)
    diff  = abs(actual - expected)
    if (diff <= tol * scale) then
      npass = npass + 1
    else
      nfail = nfail + 1
      print '(A,A,A)', 'FAIL ', label, ':'
      print *,         '  got     =', actual
      print *,         '  expected=', expected
    end if
  end subroutine check_q

  subroutine check_q_abs(actual, expected, tol, label, npass, nfail)
    real(real128),    intent(in)    :: actual, expected, tol
    character(len=*), intent(in)    :: label
    integer,          intent(inout) :: npass, nfail
    if (abs(actual - expected) <= tol) then
      npass = npass + 1
    else
      nfail = nfail + 1
      print '(A,A,A)', 'FAIL ', label, ':'
      print *,         '  got     =', actual
      print *,         '  expected=', expected
    end if
  end subroutine check_q_abs

#else
  ! Quadmath wasn't enabled at build time.  This program is still
  ! compiled (so CMake can wire up the test target unconditionally) but
  ! becomes a no-op pass.
  implicit none
  print *, 'quadmath not enabled at build time -- nothing to test'
#endif
end program test_wigner_quadmath
