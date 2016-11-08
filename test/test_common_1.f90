! This file is part of RMPCDMD
! Copyright (c) 2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program test_common_1
  use common
  use tester
  implicit none

  type(tester_t) :: test

  double precision, dimension(3) :: x1, x2, x3

  call test% init()

  x1 = [1, 0, 0]
  x2 = [0, 1, 0]
  x3 = cross(x1, x2)

  call test% assert_close(x3, [0.d0, 0.d0, 1.d0])

  x1 = [1, 2, 3]
  x2 = [4, 5, 6]
  x3 = cross(x1, x2)

  call test% assert_close(x3, [-3.d0, 6.d0, -3.d0])

  call test% print()

end program test_common_1
