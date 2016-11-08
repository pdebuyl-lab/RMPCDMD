! This file is part of RMPCDMD
! Copyright (c) 2015-2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program test_common_0
  use common
  use tester
  implicit none

  type(tester_t) :: test

  type(profile_t) :: p

  double precision :: x0, x1
  double precision, allocatable :: x(:)
  integer :: i
  integer :: n
  integer, parameter :: nloop = 100

  integer :: seed_size, clock
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  clock = 71943989
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)


  x0 = 7.d0
  x1 = 10.d0
  n = 19

  call test% init(tolerance64=x1/sqrt(dble(nloop)))

  call p% init(x0, x1, n)

  allocate(x(size(p% data)))

  do i = 1, n
     x(i) = flin( (dble(i)-0.5d0)/dble(n) )
  end do

  call test% assert_positive(p% dx)
  call test% assert_positive(p% count)

  do i = 1, nloop
     call random_number(x0)
     x0 = flin(x0)
     call p% bin(x0, x0)
  end do

  call test% assert_positive(p% count)
  call test% assert_equal(sum(p% count), nloop)

  call p% norm()

  call test% assert_close(p% data, x)

  call p% reset()
  do i = 1, n
     x(i) = cos( flin((dble(i)-0.5d0)/dble(n)) )
  end do

  do i = 1, nloop
     call random_number(x0)
     x0 = flin(x0)
     x1 = cos(x0)
     call p% bin(x0, x1)
  end do

  call p% norm()

  write(*,*) p % data
  write(*,*) x

  call test% assert_close(p% data, x)

  call test% print()

contains

  double precision elemental function flin(x)
    double precision, intent(in) :: x

    flin = 7.d0 + x * 3.d0

  end function flin

end program test_common_0
