program test_common_0
  use common
  use tester
  implicit none

  type(tester_t) :: test

  type(profile_t) :: p

  double precision :: x0, x1
  integer :: n

  call test% init()

  x0 = 7.d0
  x1 = 10.d0
  n = 19

  call p% init(x0, x1, n)

  call test% assert_positive(p% dx)
  call test% assert_positive(p% count)

  call test% print()

end program test_common_0
