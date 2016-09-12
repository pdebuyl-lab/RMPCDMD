program test_correlator_0
  use correlator
  use tester
  implicit none

  type(tester_t) :: test

  double precision, parameter :: dt = 0.2
  double precision, parameter :: gamma = atan(1.d0)*4
  double precision, parameter :: v0 = 0.03
  integer, parameter :: l = 8
  integer, parameter :: b = 4
  integer, parameter :: N=8*l**(b+1)
  type(correlator_t) :: acf1, acf2, acf3, acf4

  integer :: i, j, i_block
  double precision :: t

  double precision :: w
  double precision :: x(3)
  double precision :: y(4)
  double precision :: z(2, 5)

  double precision :: comp(l)

  call test%init(d_tol = 1d-8)

  call acf1%init(l, b)
  call acf2%init(l, b, dim=size(x))
  call acf3%init(l, b, n=size(y))
  call acf4%init(l, b, dim=size(z, dim=1), n=size(z, dim=2))

  do i = 0, N-1
     t = i*dt
     w = t*v0
     x = [-1, 0, 1]
     y = [-1, 0, 1, 2]
     do j = 1, size(z, dim=2)
        z(1, j) = cos(gamma*i)
        z(2, j) = 0
     end do
     call acf1%add(i, correlate_block_distsq, x=w)
     call acf2%add(i, correlate_block_dot, xvec=x)
     call acf3%add(i, correlate_block_dot, x_n=y)
     call acf4%add(i, correlate_block_dot, xvec_n=z)
  end do

  do i_block = 1, b
     comp = [ ((i*dt*v0*l**(i_block-1))**2, i=0, l-1) ]
     call test%assert_close(acf1%correlation(1,1,:,i_block) / acf1%count(i_block), comp)
  end do

  do i_block = 1, b
     comp = 1
     call test%assert_close(acf2%correlation(1,1,:,i_block) / acf2%count(i_block), comp)
     comp = 0
     call test%assert_close(acf2%correlation(2,1,:,i_block) / acf2%count(i_block), comp)
     comp = 1
     call test%assert_close(acf2%correlation(3,1,:,i_block) / acf2%count(i_block), comp)
  end do

  do i_block = 1, b
     comp = 1
     call test%assert_close(acf2%correlation(1,1,:,i_block) / acf2%count(i_block), comp)
     comp = 0
     call test%assert_close(acf2%correlation(2,1,:,i_block) / acf2%count(i_block), comp)
     comp = 1
     call test%assert_close(acf2%correlation(3,1,:,i_block) / acf2%count(i_block), comp)
  end do

  do i_block = 1, b
     comp = 1
     call test%assert_close(acf3%correlation(1,1,:,i_block) / acf3%count(i_block), comp)
     comp = 0
     call test%assert_close(acf3%correlation(1,2,:,i_block) / acf3%count(i_block), comp)
     comp = 1
     call test%assert_close(acf3%correlation(1,3,:,i_block) / acf3%count(i_block), comp)
     comp = 4
     call test%assert_close(acf3%correlation(1,4,:,i_block) / acf3%count(i_block), comp)
  end do

  do i_block = 1, b
     if (i_block==1) then
        comp = [ ((-1)**(modulo(i+2, 2)), i=0, l-1) ]
     else
        comp = 1
     end if
     call test%assert_close(acf4%correlation(1,1,:,i_block) / acf4%count(i_block), comp)
     write(15,*) acf4%correlation(1,1,:,i_block) / acf4%count(i_block)
  end do

  call test%print()

end program test_correlator_0
