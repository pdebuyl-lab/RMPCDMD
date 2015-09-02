module mpcd
  use particle_system
  use cell_system
  use mt19937ar_module
  implicit none

  private

  public :: compute_temperature, simple_mpcd_step

contains

  function rand_sphere(state) result(n)
    type(mt19937ar_t), intent(inout) :: state
    double precision :: n(3)

    logical :: s_lt_one
    double precision :: s, alpha

    s_lt_one = .false.
    do while (.not. s_lt_one)
       n(1) = genrand_real1(state)
       n(2) = genrand_real1(state)
       n(1:2) = 2*n(1:2) - 1
       s = n(1)**2 + n(2)**2
       if ( s<1.d0 ) s_lt_one = .true.
    end do
    alpha = 2.d0 * sqrt(1.d0 - s)
    n(1) = n(1)*alpha
    n(2) = n(2)*alpha
    n(3) = 1.d0 - 2.d0*s
  end function rand_sphere

  subroutine simple_mpcd_step(particles, cells, state, temperature)
    class(particle_system_t), intent(in) :: particles
    class(cell_system_t), intent(in) :: cells
    type(mt19937ar_t), intent(inout) :: state
    double precision, intent(in), optional :: temperature

    integer :: i, start, n
    integer :: cell_idx, count
    double precision :: local_v(3), omega(3,3), vec(3)
    logical :: thermostat

    thermostat = present(temperature)
    if (thermostat) error stop 'thermostatting not implemented'

    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)
       count = count + 1

       local_v = 0
       do i = start, start + n - 1
          local_v = local_v + particles% vel(:, i)
       end do
       local_v = local_v / n

       vec = rand_sphere(state)
       omega = &
            reshape( (/ &
            vec(1)**2, vec(1)*vec(2) + vec(3), vec(1)*vec(3) - vec(2) ,&
            vec(2)*vec(1) - vec(3) , vec(2)**2 , vec(2)*vec(3) + vec(1),&
            vec(3)*vec(1) + vec(2), vec(3)*vec(2) - vec(1), vec(3)**2 &
            /), (/3, 3/))

       do i = start, start + n - 1
          particles% vel(:, i) = local_v + matmul(omega, (particles% vel(:, i)-local_v))
       end do

    end do

  end subroutine simple_mpcd_step

  function compute_temperature(particles, cells) result(te)
    class(particle_system_t), intent(in) :: particles
    class(cell_system_t), intent(in) :: cells
    double precision :: te

    integer :: i, start, n
    integer :: cell_idx, count
    double precision :: local_v(3), local_kin

    cell_idx = 1
    count = 0
    te = 0
    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)
       count = count + 1

       local_v = 0
       local_kin = 0
       do i = start, start + n - 1
          local_v = local_v + particles% vel(:, i)
       end do
       local_v = local_v / n
       do i = start, start + n - 1
          local_kin = local_kin + sum((particles% vel(:, i)-local_v)**2)
       end do
       te = te + local_kin / dble(3*(n-1))

    end do

    te = te / count

  end function compute_temperature

end module mpcd
