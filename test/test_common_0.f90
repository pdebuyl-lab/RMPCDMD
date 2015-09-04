program test_common_0
  use common
  use tester
  implicit none

  type(tester_t) :: test

  type(profile_t) :: p

  double precision :: x0, x1
  integer :: i
  integer :: n
  integer, parameter :: nloop = 100

  integer :: seed_size, clock
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  clock = 1002581
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call test% init()

  x0 = 7.d0
  x1 = 10.d0
  n = 19

  call p% init(x0, x1, n)

  call test% assert_positive(p% dx)
  call test% assert_positive(p% count)

  do i = 1, nloop
     call random_number(x0)
     x0 = 7.d0 + x0 * 3
     x1 = x0
     call p% bin(x0, x1)
  end do

  call test% assert_positive(p% count)
  call test% assert_equal(sum(p% count), nloop)
  call test% print()

end program test_common_0
