! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Susi Lehtola
!
! Smoke test for the Fortran 90 interface of an installed libwignernj.
program test_fort
  use wignernj
  implicit none
  real(8) :: v, ref
  v   = w3j(1.0d0, 1.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0)
  ref = -1.0d0 / sqrt(3.0d0)
  if (abs(v - ref) > 1.0d-14) then
     write(0, '(A,ES23.15,A,ES23.15)') &
          "Fortran: w3j returned ", v, ", expected ", ref
     stop 1
  end if
  print '(A,F18.15)', "Fortran: w3j(1,1,0;0,0,0) = ", v
end program test_fort
