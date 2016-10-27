!> Test the elastic network force routine
!!
!! The geometry is a cubic box of side 5 where the atoms are placed following this diagram
!! <pre>
!!   3
!!   |
!! --1--2    4--
!! </pre>
!!
!! The link between '1' and '4' crosses the periodic box.

program test_elastic_network_0
  use particle_system
  use md
  use iso_c_binding
  use tester
  implicit none

  type(tester_t) :: test

  type(particle_system_t) :: p
  integer, parameter :: N = 4
  integer, parameter :: n_links = 3
  double precision, parameter :: k = 7
  integer :: links(2,n_links)
  double precision :: links_d(n_links)

  double precision :: L(3), e, expected_e

  call test%init()

  call p% init(N)
  L = 5.d0

  p%force = 0
  p%pos = 1
  p%pos(1,1) = 1

  expected_e = 0
  ! first link
  ! r0 = 1/2
  ! distance = 2
  p%pos(1,2) = 3
  links(:,1) = [1, 2]
  links_d(1) = 0.5d0
  expected_e = expected_e + 1.5d0**2
  ! second link
  ! r0 = 1
  ! distance = 1
  p%pos(2,3) = 2
  links(:,2) = [1, 3]
  links_d(2) = 1
  expected_e = expected_e + 0
  ! third link
  ! r0 = 1
  ! distance = 1.5
  p%pos(1,4) = 4.5d0
  links(:,3) = [1, 4]
  links_d(3) = 1
  expected_e = expected_e + 0.5d0**2

  expected_e = expected_e*k/2

  e = elastic_network(p, links, links_d, k, L)

  write(*,*) 'expected_e', expected_e
  write(*,*) 'e', e
  call test%assert_close(e, expected_e)

  call test%assert_close(p%force(:,1), [k, 0.d0, 0.d0])

  call test%assert_close(maxval(abs(sum(p%force, dim=2))), 0.d0)

  call test% print()

end program test_elastic_network_0
