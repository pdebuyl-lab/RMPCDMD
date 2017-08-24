! This file is part of RMPCDMD
! Copyright (c) 2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program test_md_0
  use rmpcdmd_module
  use threefry_module
  use iso_c_binding
  use tester
  implicit none

  integer, parameter :: N = 679

  type(tester_t) :: test

  type(threefry_rng_t) :: state(1)
  integer :: clock
  integer :: i, j
  double precision :: x(3,N), v(3,N)
  integer :: L(3)
  double precision :: x_old(3,N), v_old(3,N)
  logical :: will_collide(N)

  type(cell_system_t) :: cells
  type(particle_system_t) :: p

  call test% init()
  call system_clock(count=clock)
  call threefry_rng_init(state, int(clock, c_int64_t))

  L = [ 1, 20, 5 ]
  call cells%init(L, 1.d0)
  call p%init(N,1)

  cells%is_md = .false.

  cells%bc = [PERIODIC_BC, BOUNCE_BACK_BC, PERIODIC_BC]
  do j = 1, 10

     call p%random_placement(cells%edges, state=state(1))
     call p%sort(cells)

     do i = 1, N
        p%vel(1,:) = (-0.5d0+threefry_double(state(1)))*2*L(1)
        p%vel(2,:) = (-0.5d0+threefry_double(state(1)))*2*L(2)
        p%vel(3,:) = (-0.5d0+threefry_double(state(1)))*2*L(3)
     end do

     x_old = p%pos
     v_old = p%vel

     will_collide = ( (x_old(2,:)+v_old(2,:) < 0) .or. &
          (x_old(2,:)+v_old(2,:) > L(2)) )

     call cell_md_pos_ywall(cells, p, 1.d0, md_flag=.false.)

     x = p%pos
     v = p%vel

     do i = 1, N
        call test%assert_close(norm2(v(:,i)), norm2(v_old(:,i)))
        call test%assert_positive(x(2,i))
        call test%assert_positive(L(2)-x(2,i))
        if (will_collide(i)) call test%assert_close(v(:,i), -v_old(:,i))
     end do

  end do

  cells%bc = [PERIODIC_BC, SPECULAR_BC, PERIODIC_BC]
  do j = 1, 10

     call p%random_placement(cells%edges, state=state(1))
     call p%sort(cells)

     do i = 1, N
        p%vel(1,:) = (-0.5d0+threefry_double(state(1)))*2*L(1)
        p%vel(2,:) = (-0.5d0+threefry_double(state(1)))*2*L(2)
        p%vel(3,:) = (-0.5d0+threefry_double(state(1)))*2*L(3)
     end do

     x_old = p%pos
     v_old = p%vel

     will_collide = ( (x_old(2,:)+v_old(2,:) < 0) .or. &
          (x_old(2,:)+v_old(2,:) > L(2)) )

     call cell_md_pos_ywall(cells, p, 1.d0, md_flag=.false.)

     x = p%pos
     v = p%vel

     do i = 1, N
        call test%assert_close(norm2(v(:,i)), norm2(v_old(:,i)))
        call test%assert_positive(x(2,i))
        call test%assert_positive(L(2)-x(2,i))
        if (will_collide(i)) then
           call test%assert_close(v(1,i), v_old(1,i))
           call test%assert_close(v(2,i), -v_old(2,i))
           call test%assert_close(v(3,i), v_old(3,i))
        end if
     end do

  end do

  call test% print()

end program test_md_0
