! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! libwignernj Fortran API quadmath (real(real128)) demonstration.
!
! Exercises the w3jq / w6jq / w9jq / wcgq / wracahwq / wfanoxq / wgauntq /
! wgaunt_realq / wreal_ylm_in_complex_ylmq convenience wrappers added by
! the wignernj module when libwignernj is built with -DBUILD_QUADMATH=ON.
!
! Built only when the C library is configured with -DBUILD_QUADMATH=ON
! (the WIGNERNJ_HAVE_QUADMATH macro is then defined for downstream
! preprocessing).  Compile with -cpp; .F90 extension enables the
! preprocessor in gfortran by default.

program quadmath_demo
#ifdef WIGNERNJ_HAVE_QUADMATH
  use wignernj
  use iso_c_binding,   only: c_float128_complex
  use iso_fortran_env, only: real128, error_unit
  implicit none

  integer :: failed
  real(real128) :: val, ref
  ! 4 ulp at quad precision (FLT128_EPSILON ~= 1.93e-34).
  real(real128), parameter :: tol = &
      4.0_real128 * 1.92592994438723585305597794258493e-34_real128
  complex(c_float128_complex) :: cmat(3, 3)

  failed = 0

  print '(A)', 'libwignernj Fortran quadmath API demonstration'
  print '(A)', '----------------------------------------------'

  ! 3j(1,1,0; 0,0,0) = -1/sqrt(3)
  val = w3jq(1.0_real128, 1.0_real128, 0.0_real128, &
             0.0_real128, 0.0_real128, 0.0_real128)
  ref = -1.0_real128 / sqrt(3.0_real128)
  call check('w3jq(1,1,0; 0,0,0)', val, ref)

  ! 6j{1,1,1; 1,1,1} = 1/6
  val = w6jq(1.0_real128, 1.0_real128, 1.0_real128, &
             1.0_real128, 1.0_real128, 1.0_real128)
  ref = 1.0_real128 / 6.0_real128
  call check('w6jq{1,1,1; 1,1,1}', val, ref)

  ! 9j{1/2,1/2,1; 1/2,1/2,0; 1,1,1} = sqrt(6)/18
  val = w9jq(0.5_real128, 0.5_real128, 1.0_real128, &
             0.5_real128, 0.5_real128, 0.0_real128, &
             1.0_real128, 1.0_real128, 1.0_real128)
  ref = sqrt(6.0_real128) / 18.0_real128
  call check('w9jq{1/2,1/2,1; 1/2,1/2,0; 1,1,1}', val, ref)

  ! CG(1/2 1/2; 1/2 -1/2 | 1 0) = 1/sqrt(2)
  val = wcgq(0.5_real128,  0.5_real128, &
             0.5_real128, -0.5_real128, &
             1.0_real128,  0.0_real128)
  ref = 1.0_real128 / sqrt(2.0_real128)
  call check('wcgq(1/2,1/2; 1/2,-1/2 | 1,0)', val, ref)

  ! Racah W(1,1,1,1; 1,1) = (-1)^4 * 6j{1,1,1; 1,1,1} = 1/6
  val = wracahwq(1.0_real128, 1.0_real128, 1.0_real128, &
                 1.0_real128, 1.0_real128, 1.0_real128)
  ref = 1.0_real128 / 6.0_real128
  call check('wracahwq(1,1,1,1; 1,1)', val, ref)

  ! Fano X at all-j=2 = 25 * 9j{2,2,2; 2,2,2; 2,2,2}
  val = wfanoxq(2.0_real128, 2.0_real128, 2.0_real128, &
                2.0_real128, 2.0_real128, 2.0_real128, &
                2.0_real128, 2.0_real128, 2.0_real128)
  ref = 25.0_real128 * w9jq(2.0_real128, 2.0_real128, 2.0_real128, &
                            2.0_real128, 2.0_real128, 2.0_real128, &
                            2.0_real128, 2.0_real128, 2.0_real128)
  call check('wfanoxq(2,2,2; 2,2,2; 2,2,2)', val, ref)

  ! Real Gaunt at all-m=0 equals complex Gaunt.
  val = wgaunt_realq(1.0_real128, 0.0_real128, &
                     1.0_real128, 0.0_real128, &
                     2.0_real128, 0.0_real128) &
      - wgauntq(1.0_real128, 0.0_real128, &
                1.0_real128, 0.0_real128, &
                2.0_real128, 0.0_real128)
  ref = 0.0_real128
  call check('wgaunt_realq - wgauntq at all-m=0', val, ref)

  ! Real <-> complex Y_lm at l=1: C(+1,-1) real part is 1/sqrt(2).
  call wreal_ylm_in_complex_ylmq(1, cmat)
  val = real(cmat(3, 1), kind=real128)
  ref = 1.0_real128 / sqrt(2.0_real128)
  call check('wreal_ylm_in_complex_ylmq[l=1] re(+1,-1)', val, ref)

  if (failed /= 0) then
    write(error_unit, '(I0, A)') failed, ' check(s) failed'
    stop 1
  end if
  print '(A)', ''
  print '(A)', 'All quadmath checks passed.'

contains

  subroutine check(label, v, e)
    character(len=*), intent(in) :: label
    real(real128),    intent(in) :: v, e
    real(real128) :: d, sc
    print '(A, A, ES45.34E3)',          '  ', label, v
    print '(A, T49, A, ES45.34E3, A)',  '  ', '(expected ', e, ')'
    d  = abs(v - e)
    sc = merge(abs(e), 1.0_real128, abs(e) > 1.0e-300_real128)
    if (d > tol * sc) then
      write(error_unit, '(A, A, A, ES12.3E3)') &
          '  FAIL: ', label, ' diff = ', d
      failed = failed + 1
    end if
  end subroutine check

#else
  implicit none
  print '(A)', 'libwignernj was built without -DBUILD_QUADMATH=ON; nothing to demo.'
#endif
end program quadmath_demo
