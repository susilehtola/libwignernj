! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! Fortran test program for the wigner module convenience wrappers.
program test_wigner_fortran
  use wignernj
  use iso_c_binding, only: c_double, c_int, c_double_complex
  implicit none

  integer  :: npass, nfail, total
  real(c_double) :: val, exp_val
  real(c_double) :: w, w6val, phase_val, g, w0, norm
  real(c_double), parameter :: pi = 3.14159265358979323846d0

  npass = 0
  nfail = 0

  ! ── 3j symbol ──────────────────────────────────────────────────────────────

  ! (1 1 0; 0 0 0) = -1/sqrt(3)
  val = w3j(1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0)
  exp_val = -1.0d0 / sqrt(3.0d0)
  call check_near(val, exp_val, 2.0d-14, 'w3j(1,1,0;0,0,0)', npass, nfail)

  ! (1/2 1/2 1; 1/2 -1/2 0) = 1/sqrt(6)
  val = w3j(0.5d0, 0.5d0, 1.0d0, 0.5d0, -0.5d0, 0.0d0)
  exp_val = 1.0d0 / sqrt(6.0d0)
  call check_near(val, exp_val, 2.0d-14, 'w3j(1/2,1/2,1;1/2,-1/2,0)', npass, nfail)

  ! Parity zero: (1 1 1; 0 0 0) = 0
  val = w3j(1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0)
  call check_abs(val, 0.0d0, 1.0d-15, 'w3j(1,1,1;0,0,0)=0', npass, nfail)

  ! (10 10 0; 0 0 0) = 1/sqrt(21)
  val = w3j(10.0d0, 10.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0)
  exp_val = 1.0d0 / sqrt(21.0d0)
  call check_near(val, exp_val, 2.0d-14, 'w3j(10,10,0;0,0,0)', npass, nfail)

  ! Raw integer API: wigner3j(2,2,0, 0,0,0) = -1/sqrt(3)
  val = wigner3j(2_c_int, 2_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int)
  call check_near(val, -1.0d0/sqrt(3.0d0), 2.0d-14, 'wigner3j raw', npass, nfail)

  ! ── 6j symbol ──────────────────────────────────────────────────────────────

  ! {1 1 1; 1 1 1} = 1/6
  val = w6j(1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0)
  call check_near(val, 1.0d0/6.0d0, 2.0d-14, 'w6j(1,1,1;1,1,1)', npass, nfail)

  ! {1/2 1/2 1; 1/2 1/2 0} = 1/2
  val = w6j(0.5d0, 0.5d0, 1.0d0, 0.5d0, 0.5d0, 0.0d0)
  call check_near(val, 0.5d0, 2.0d-14, 'w6j(1/2,1/2,1;1/2,1/2,0)', npass, nfail)

  ! Triangle violated: {0 0 1; 0 0 1} = 0
  val = w6j(0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0)
  call check_abs(val, 0.0d0, 1.0d-15, 'w6j triangle zero', npass, nfail)

  ! ── 9j symbol ──────────────────────────────────────────────────────────────

  ! {1/2 1/2 1; 1/2 1/2 0; 1 1 1} = sqrt(6)/18
  val = w9j(0.5d0, 0.5d0, 1.0d0,  0.5d0, 0.5d0, 0.0d0,  1.0d0, 1.0d0, 1.0d0)
  call check_near(val, sqrt(6.0d0)/18.0d0, 2.0d-14, 'w9j(sqrt6/18)', npass, nfail)

  ! Selection-rule zero
  val = w9j(0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0)
  call check_abs(val, 0.0d0, 1.0d-14, 'w9j triangle zero', npass, nfail)

  ! ── Clebsch-Gordan ─────────────────────────────────────────────────────────

  ! <1/2 1/2; 1/2 -1/2 | 1 0> = 1/sqrt(2)
  val = wcg(0.5d0, 0.5d0, 0.5d0, -0.5d0, 1.0d0, 0.0d0)
  call check_near(val, 1.0d0/sqrt(2.0d0), 2.0d-14, 'wcg(1/2,1/2;1/2,-1/2|1,0)', npass, nfail)

  ! <1/2 -1/2; 1/2 1/2 | 0 0> = -1/sqrt(2)
  val = wcg(0.5d0, -0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0)
  call check_near(val, -1.0d0/sqrt(2.0d0), 2.0d-14, 'wcg(1/2,-1/2;1/2,1/2|0,0)', npass, nfail)

  ! m conservation: m1+m2 != M → 0
  val = wcg(1.0d0, 1.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0)
  call check_abs(val, 0.0d0, 1.0d-15, 'wcg m-conservation zero', npass, nfail)

  ! ── Racah W ────────────────────────────────────────────────────────────────

  ! W(1,1,0,1;0,1) vs phase * 6j{1,1,0;1,0,1}
  ! phase = (-1)^((2+2+0+2)/2) = (-1)^3 = -1
  w     = wracahw(1.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0)
  w6val = w6j(1.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0)
  phase_val = -1.0d0
  call check_near(w, phase_val*w6val, 2.0d-14, 'wracahw vs w6j', npass, nfail)

  ! ── Gaunt ──────────────────────────────────────────────────────────────────

  ! G(1,0,1,0,0,0) = norm * 3j(2,2,0;0,0,0)^2
  g    = wgaunt(1.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0)
  w0   = w3j(1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0)
  norm = sqrt(3.0d0 * 3.0d0 * 1.0d0 / (4.0d0 * pi))
  call check_near(g, norm*w0*w0, 2.0d-14, 'wgaunt(1,0,1,0,0,0)', npass, nfail)

  ! m1+m2+m3 != 0 → 0
  val = wgaunt(1.0d0, 1.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0)
  call check_abs(val, 0.0d0, 1.0d-15, 'wgaunt m-sum zero', npass, nfail)

  ! ── Real-spherical-harmonic Gaunt ───────────────────────────────────────────

  ! G^R(0,0,0,0,0,0) = G(0,0,0,0,0,0) = 1/(2*sqrt(pi))
  val = wgaunt_real(0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0)
  exp_val = 0.5d0 / sqrt(pi)
  call check_near(val, exp_val, 2.0d-14, 'wgaunt_real(0,0,0,0,0,0)', npass, nfail)

  ! Real Gaunt agrees with complex Gaunt at all-m=0
  val = wgaunt_real(1.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0)
  exp_val = wgaunt(1.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0)
  call check_near(val, exp_val, 2.0d-14, 'wgaunt_real m=0 = wgaunt', npass, nfail)

  ! ── Fano X ─────────────────────────────────────────────────────────────────

  ! X(j1,j2,j12;j3,j4,j34;j13,j24,J)
  !   = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] * {9j}
  ! At all-equal-j=1, the four (2j12+1) factors are 3 each, so
  ! norm = sqrt(3^4) = 9; at all-equal-j=2 the norm is 5^2 = 25.
  val     = wfanox(1.0d0,1.0d0,1.0d0, 1.0d0,1.0d0,1.0d0, 1.0d0,1.0d0,1.0d0)
  exp_val = 9.0d0 * w9j(1.0d0,1.0d0,1.0d0, 1.0d0,1.0d0,1.0d0, 1.0d0,1.0d0,1.0d0)
  call check_near(val, exp_val, 2.0d-14, 'wfanox = 9 * 9j (all-j=1)', npass, nfail)

  val     = wfanox(2.0d0,2.0d0,2.0d0, 2.0d0,2.0d0,2.0d0, 2.0d0,2.0d0,2.0d0)
  exp_val = 25.0d0 * w9j(2.0d0,2.0d0,2.0d0, 2.0d0,2.0d0,2.0d0, 2.0d0,2.0d0,2.0d0)
  call check_near(val, exp_val, 2.0d-14, 'wfanox = 25 * 9j (all-j=2)', npass, nfail)

  ! Selection-rule zero (one triangle violated) → X = 0
  val = wfanox(1.0d0,1.0d0,1.0d0, 1.0d0,1.0d0,1.0d0, 1.0d0,1.0d0,3.0d0)
  call check_abs(val, 0.0d0, 1.0d-15, 'wfanox triangle zero', npass, nfail)

  ! ── Real <-> complex Y_lm basis-overlap matrix ─────────────────────────────
  block
    complex(c_double_complex) :: C(3, 3)
    real(c_double) :: s
    s = 1.0d0 / sqrt(2.0d0)
    call wreal_ylm_in_complex_ylm(1, C)
    ! Fortran column-major: C(row, col), row = m_r+l+1, col = m_c+l+1.
    ! (m_r=+1, m_c=-1) -> C(3, 1): real part = +1/sqrt(2)
    call check_near(real(C(3, 1)), s, 2.0d-14, &
                    'wreal_ylm_in_complex_ylm[+1,-1].re', npass, nfail)
    ! (m_r=-1, m_c=-1) -> C(1, 1): imag part = +1/sqrt(2)
    call check_near(aimag(C(1, 1)), s, 2.0d-14, &
                    'wreal_ylm_in_complex_ylm[-1,-1].im', npass, nfail)
    ! (m_r=0, m_c=0) -> C(2, 2): real part = 1
    call check_near(real(C(2, 2)), 1.0d0, 2.0d-14, &
                    'wreal_ylm_in_complex_ylm[ 0, 0].re', npass, nfail)
  end block

  ! ── Summary ────────────────────────────────────────────────────────────────

  total = npass + nfail
  if (nfail == 0) then
    write(*,'(A,I0,A,I0,A)') 'PASS: ', npass, '/', total, ' tests passed'
    stop 0
  else
    write(*,'(A,I0,A,I0,A)') 'FAIL: ', npass, '/', total, ' tests passed'
    stop 1
  end if

contains

  subroutine check_near(got, expected, tol, label, np, nf)
    real(c_double), intent(in)   :: got, expected, tol
    character(len=*), intent(in) :: label
    integer, intent(inout)       :: np, nf
    real(c_double) :: diff, ref
    diff = abs(got - expected)
    ref  = max(abs(expected), 1.0d-300)
    if (diff <= tol * ref) then
      np = np + 1
    else
      nf = nf + 1
      write(*,'(A,A,A,ES15.8,A,ES15.8,A,ES10.3)') &
        'FAIL: ', trim(label), '  got=', got, '  exp=', expected, '  diff=', diff
    end if
  end subroutine check_near

  subroutine check_abs(got, expected, tol, label, np, nf)
    real(c_double), intent(in)   :: got, expected, tol
    character(len=*), intent(in) :: label
    integer, intent(inout)       :: np, nf
    real(c_double) :: diff
    diff = abs(got - expected)
    if (diff <= tol) then
      np = np + 1
    else
      nf = nf + 1
      write(*,'(A,A,A,ES15.8,A,ES15.8,A,ES10.3)') &
        'FAIL: ', trim(label), '  got=', got, '  exp=', expected, '  diff=', diff
    end if
  end subroutine check_abs

end program test_wigner_fortran
