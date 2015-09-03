module mpcd
  use particle_system
  use cell_system
  use mt19937ar_module
  implicit none

  private

  public :: compute_temperature, simple_mpcd_step
  public :: wall_mpcd_step
  public :: mpcd_stream

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

  !! Lamura, Gompper, Ihle and Kroll, EPL 56, 319-325 (2001)
  !! http://dx.doi.org/10.1209/epl/i2001-00522-9
  subroutine wall_mpcd_step(particles, cells, state, wall_temperature, wall_v, wall_n)
    use hilbert
    class(particle_system_t), intent(in) :: particles
    class(cell_system_t), intent(in) :: cells
    type(mt19937ar_t), intent(inout) :: state
    double precision, intent(in) :: wall_temperature(2)
    double precision, intent(in) :: wall_v(3,2)
    integer, intent(in) :: wall_n(2)

    integer :: i, start, n
    integer :: cell_idx, count
    double precision :: local_v(3), omega(3,3), vec(3)
    integer :: n_virtual
    integer :: cell(3)
    integer :: wall_idx
    double precision :: virtual_v(3)

    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)
       count = count + 1
       n_virtual = 0

       ! Find whether we are in a wall cell
       cell = compact_h_to_p(cell_idx, cells% M)
       if (cell(3) == 1) then
          wall_idx = 1
          vec = wall_v(:,wall_idx)
       else if (cell(3) == cells% L(3)) then
          wall_idx = 2
          vec = wall_v(:,wall_idx)
       else
          wall_idx = -1
          vec = local_v
       end if

       if ( (wall_idx > 0) .and. (n < wall_n(wall_idx)) ) then
          n_virtual = wall_n(wall_idx) - n
          call mt_normal_data(virtual_v, state)
          virtual_v = virtual_v * sqrt(n_virtual*wall_temperature(wall_idx))
       end if

       local_v = 0
       do i = start, start + n - 1
          local_v = local_v + particles% vel(:, i)
       end do
       if (n_virtual > 0) then
          local_v = (local_v + virtual_v) / wall_n(wall_idx)
       else
          local_v = local_v / n
       end if

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

  end subroutine wall_mpcd_step

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

  !! Advance mpcd particles
  !
  ! If the cell system has a wall in the z direction, a bounce-back collision is used.
  subroutine mpcd_stream(particles, cells, dt)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt

    integer :: i
    double precision :: pos_min(3), pos_max(3), delta

    pos_min(3) = 0
    pos_max(3) = cells% L(3) * cells% a

    do i = 1, particles% Nmax
       particles% pos(:,i) = particles% pos(:,i) + particles% vel(:,i)*dt
       particles% pos(1,i) = modulo( particles% pos(1,i) , cells% L(1)*cells% a )
       particles% pos(2,i) = modulo( particles% pos(2,i) , cells% L(2)*cells% a )
       if (cells% has_walls) then
          if (particles% pos(3,i) < pos_min(3)) then
             ! bounce position
             delta = pos_min(3) - particles% pos(3,i)
             particles% pos(3,i) = pos_min(3) + delta
             ! bounce velocity
             particles% vel(3,i) = -particles% vel(3,i)
          else if (particles% pos(3,i) > pos_max(3)) then
             ! bounce position
             delta = particles% pos(3,i) - pos_max(3)
             particles% pos(3,i) = pos_max(3) - delta
             ! bounce velocity
             particles% vel(3,i) = -particles% vel(3,i)
          end if
       else
          particles% pos(2,i) = modulo( particles% pos(2,i) , cells% L(2)*cells% a )
       end if
    end do

  end subroutine mpcd_stream

end module mpcd