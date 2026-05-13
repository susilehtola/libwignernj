! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! Build the orbital angular momentum operator l_z in the real
! spherical harmonic basis from its diagonal complex-basis form,
! using libwignernj's wreal_ylm_in_complex_ylm as the change of basis.
! Fortran parallel to examples/c/real_basis_lz.c.
!
! In the complex Y_l^m basis, l_z is diagonal with entry m delta_{m m'}
! (in units of hbar).  libwignernj returns the unitary C that encodes
! the basis-vector relation
!
!     S_{l, m_r}  =  sum_{m_c}  C(m_r+l+1, m_c+l+1)  Y_l^{m_c}.
!
! The corresponding operator-matrix similarity transform is then
!
!     l_z_real  =  conjg(C)  @  l_z_complex  @  transpose(C),
!
! which expands the textbook  <S_{m_r}|O|S_{m_r'}>  with the ket
! using C(m_r'+l+1, m_c'+l+1) and the bra using conjg(C(m_r+l+1,
! m_c+l+1)).  The result is purely imaginary and antisymmetric: l_z
! couples each cosine-type real harmonic to its sin-type partner with
! strength +/- i times the diagonal m_c entry.
!
! For l = 1 this gives the textbook
!
!   l_z = (  0   0   i )
!         (  0   0   0 )
!         ( -i   0   0 )
!
! in the basis ordering m_r = -1, 0, +1.
!
! Build (out-of-tree, against an installed libwignernj):
!     gfortran -o real_basis_lz_f real_basis_lz.f90 \
!         -I /path/to/installed/include \
!         -L /path/to/installed/lib -lwignernj_f03 -lwignernj
program real_basis_lz
  use, intrinsic :: iso_c_binding, only: c_double, c_double_complex
  use wignernj
  implicit none

  integer, parameter :: l   = 1
  integer, parameter :: dim = 2 * l + 1
  complex(c_double_complex) :: C(dim, dim), Cstar(dim, dim), CT(dim, dim)
  complex(c_double_complex) :: Lzc(dim, dim), Lzr(dim, dim)
  integer :: i, j, m
  real(c_double) :: max_err

  ! Fill C: basis-vector matrix, S = C Y.
  call wreal_ylm_in_complex_ylm(l, C)

  ! Cstar: elementwise conjugate.  CT: transpose (not conjugated).
  Cstar = conjg(C)
  CT    = transpose(C)

  ! Diagonal complex-basis l_z, m index = m + l + 1.
  Lzc = (0.0_c_double, 0.0_c_double)
  do m = -l, l
     Lzc(m + l + 1, m + l + 1) = cmplx(real(m, c_double), 0.0_c_double, &
                                       c_double_complex)
  end do

  ! Lzr = Cstar @ Lzc @ CT, via the intrinsic matmul.
  Lzr = matmul(Cstar, matmul(Lzc, CT))

  print '(A,I0,A)', 'l = ', l, ' real-basis representation of l_z'
  print '(A)',      'Basis ordering: m_r = -l, -l+1, ..., l-1, l'
  print '(A)',      '(real-Y convention matches libwignernj''s gaunt_real)'
  print '(A)',      ''

  call print_matrix('C (basis-vector matrix, S = C Y)',       C)
  call print_matrix('l_z_complex (diagonal, entry = m)',     Lzc)
  call print_matrix('l_z_real = conj(C) @ l_z_complex @ C^T', Lzr)

  ! Hermiticity check.
  max_err = 0.0_c_double
  do j = 1, dim
     do i = 1, dim
        max_err = max(max_err, abs(Lzr(i, j) - conjg(Lzr(j, i))))
     end do
  end do
  print '(A, ES10.3)', 'Hermiticity residual max |L - L^H| = ', max_err

  if (max_err < 1.0e-13_c_double) then
     stop 0
  else
     stop 1
  end if

contains

  subroutine print_matrix(name, M)
    character(*),               intent(in) :: name
    complex(c_double_complex), intent(in) :: M(:, :)
    integer :: ii, jj, n
    real(c_double) :: re, im
    n = size(M, 1)
    print '(A, A)', name, ' ='
    do ii = 1, n
       write(*, '(A)', advance='no') '  '
       do jj = 1, n
          re = real(M(ii, jj))
          im = aimag(M(ii, jj))
          if (abs(re) < 1.0e-15_c_double) re = 0.0_c_double
          if (abs(im) < 1.0e-15_c_double) im = 0.0_c_double
          write(*, '("  ", SP, F6.3, SP, F6.3, "i")', advance='no') re, im
       end do
       print '(A)', ''
    end do
    print '(A)', ''
  end subroutine print_matrix

end program real_basis_lz
