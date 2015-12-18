module mpcd
  use common
  use particle_system
  use cell_system
  use mt19937ar_module
  implicit none

  private

  public :: compute_temperature, simple_mpcd_step
  public :: wall_mpcd_step
  public :: mpcd_stream_periodic, mpcd_stream_zwall
  public :: compute_rho, compute_vx

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
    integer :: cell_idx
    double precision :: local_v(3), omega(3,3), vec(3)
    logical :: thermostat

    thermostat = present(temperature)
    if (thermostat) error stop 'thermostatting not implemented'

    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)

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
  subroutine wall_mpcd_step(particles, cells, state, wall_temperature, wall_v, wall_n, bulk_temperature)
    use hilbert
    class(particle_system_t), intent(in) :: particles
    class(cell_system_t), intent(in) :: cells
    type(mt19937ar_t), intent(inout) :: state
    double precision, optional, intent(in) :: wall_temperature(2)
    double precision, optional, intent(in) :: wall_v(3,2)
    integer, optional, intent(in) :: wall_n(2)
    double precision, intent(in), optional :: bulk_temperature

    integer :: i, start, n
    integer :: cell_idx
    double precision :: local_v(3), omega(3,3), vec(3)
    integer :: n_virtual
    integer :: cell(3)
    integer :: wall_idx
    double precision :: virtual_v(3), t_factor
    logical :: all_present, all_absent
    logical :: bulk_thermostat

    all_present = present(wall_temperature) .and. present(wall_v) .and. present(wall_n)
    all_absent = .not. present(wall_temperature) .and. .not. present(wall_v) .and. .not. present(wall_n)
    if ( .not. (all_present .or. all_absent) ) &
         error stop 'wall parameters must be all present or all absent in wall_mpcd_step'

    if (present(bulk_temperature)) then
       bulk_thermostat = .true.
       t_factor = sqrt(bulk_temperature)
    else
       bulk_thermostat = .false.
    end if

    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)
       n_virtual = 0

       ! Find whether we are in a wall cell
       cell = compact_h_to_p(cell_idx - 1, cells% M) + 1
       if (all_present .and. (cell(3) == 1)) then
          wall_idx = 1
       else if (all_present .and. (cell(3) == cells% L(3))) then
          wall_idx = 2
       else
          wall_idx = -1
       end if

       if (wall_idx > 0) then
          if (n < wall_n(wall_idx)) then
             n_virtual = wall_n(wall_idx) - n
             call mt_normal_data(virtual_v, state)
             virtual_v = virtual_v * sqrt(n_virtual*wall_temperature(wall_idx))
          end if
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

       if (bulk_thermostat) then
          virtual_v = 0
          do i = start, start + n - 1
             call mt_normal_data(particles% vel(:, i), state)
             particles% vel(:, i) = particles% vel(:, i)*t_factor
             virtual_v = virtual_v + particles% vel(:, i)
          end do
          virtual_v = local_v - virtual_v / dble(n)
          do i = start, start + n - 1
             particles% vel(:, i) = particles% vel(:, i) + virtual_v
          end do
       else
          vec = rand_sphere(state)
          omega = &
               reshape( (/ &
               vec(1)**2, vec(1)*vec(2) + vec(3), vec(1)*vec(3) - vec(2) ,&
               vec(2)*vec(1) - vec(3) , vec(2)**2 , vec(2)*vec(3) + vec(1),&
               vec(3)*vec(1) + vec(2), vec(3)*vec(2) - vec(1), vec(3)**2 &
               /), (/3, 3/))
       end if
    end do

  end subroutine wall_mpcd_step

  function compute_temperature(particles, cells, tz) result(te)
    use hilbert, only : compact_h_to_p
    type(particle_system_t), intent(in) :: particles
    type(cell_system_t), intent(in) :: cells
    type(profile_t), intent(inout), optional :: tz

    double precision :: te

    integer :: i, start, n
    integer :: cell_idx, count
    integer :: cell(3)
    double precision :: local_v(3), local_kin
    logical :: do_tz

    if (present(tz)) then
       do_tz = .true.
       if (.not. (allocated(tz% data) .and. allocated(tz% count))) error stop 'profile_t: data not allocated'
    else
       do_tz = .false.
    end if

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
       local_kin = local_kin / dble(3*(n-1))
       te = te + local_kin

       if (do_tz) then
          cell = compact_h_to_p(cell_idx - 1, cells% M)
          call tz% bin(cells% origin(3) + (cell(3)+0.5d0)*cells% a, local_kin)
       end if

    end do

    te = te / count

  end function compute_temperature

  !! Compute density profile along z
  subroutine compute_rho(particles, rhoz)
    type(particle_system_t), intent(in) :: particles
    type(histogram_t), intent(inout) :: rhoz

    integer :: i

    if (.not. allocated(rhoz% data)) error stop 'histogram_t: data not allocated'

    do i = 1, particles% Nmax
       call rhoz% bin(particles% pos(3, i))
    end do

  end subroutine compute_rho
  
  subroutine compute_vx(particles, vx)
    type(particle_system_t), intent(in) :: particles
    type(profile_t), intent(inout) :: vx

    integer :: i

    if (.not. allocated(vx% data)) error stop 'histogram_t: data not allocated'

    do i = 1, particles% Nmax
       call vx% bin(particles% pos(3, i), particles% vel(1, i))
    end do

  end subroutine compute_vx

  !! Advance mpcd particles
  !!
  !! If the cell system has a wall in the z direction, a bounce-back collision is used.
  subroutine mpcd_stream_zwall(particles, cells, dt,g)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt

    integer :: i
    double precision :: pos_min(3), pos_max(3), delta
    double precision, dimension(3), intent(in):: g
    double precision, dimension(3) :: old_pos, old_vel
    double precision :: t_c, t_b, t_ab
    double precision :: time

    pos_min = 0
    pos_max = cells% edges

    do i = 1, particles% Nmax
       old_pos = particles% pos(:,i) 
       old_vel = particles% vel(:,i)
       particles% pos(:,i) = particles% pos(:,i) + particles% vel(:,i)*dt + g*dt**2/2
       particles% pos(1,i) = modulo( particles% pos(1,i) , cells% edges(1) )
       particles% pos(2,i) = modulo( particles% pos(2,i) , cells% edges(2) )
       particles% vel(:,i) = particles% vel(:,i) + g*dt
       if (cells% has_walls) then
          if (particles% pos(3,i) < pos_min(3)) then
             if (g(3)<0) then
                t_c = (-old_vel(3) - sqrt(old_vel(3)**2 - 2*g(3)*old_pos(3)))/g(3)
                t_b = -2*(old_pos(3) + g(3)*t_c)/g(3)
                t_ab = modulo(t_b,dt-t_c)
                ! bounce velocity
                particles% vel(3,i) = -(old_vel(3) + g(3)*t_c) + g(3)*t_ab
                particles% vel(2,i) = -particles% vel(2,i)
                particles% vel(1,i) = -particles% vel(1,i)
                ! bounce position
                particles% pos(3,i) = -(old_vel(3) + g(3)*t_c)*t_ab + g(3)*t_ab**2/2
                particles% pos(2,i) = old_pos(2) + old_vel(2)*t_c + particles% vel(2,i)*t_ab
                particles% pos(1,i) = old_pos(2) + old_vel(2)*t_c + particles% vel(1,i)*t_ab
             else if (g(3)>0) then
                t_c = (-old_vel(3) - sqrt(old_vel(3)**2 - 2*g(3)*old_pos(3)))/g(3)
                ! bounce velocity
                particles% vel(3,i) = -(old_vel(3) + g(3)*t_c) + g(3)*(dt - t_c)
                particles% vel(2,i) = -particles% vel(2,i)
                particles% vel(1,i) = -particles% vel(1,i)
                ! bounce position
                particles% pos(3,i) = -(old_vel(3) + g(3)*t_c)*(dt - t_c) + g(3)*(dt - t_c)**2/2
                particles% pos(2,i) = old_pos(2) + old_vel(2)*t_c + particles% vel(2,i)*(dt-t_c)
                particles% pos(1,i) = old_pos(2) + old_vel(2)*t_c + particles% vel(1,i)*(dt-t_c)
             else !no gravity in this direction
                ! bounce velocity
                particles% vel(:,i) = -particles% vel(:,i)
                ! bounce position
                t_c = abs(old_pos(3)/old_vel(3))
                particles% pos(:,i) = old_pos + old_vel*t_c + particles% vel(:,i)*(dt - t_c)
             end if
          else if (particles% pos(3,i) > pos_max(3)) then
             if (g(3)>0) then
                t_c = (-old_vel(3) + sqrt(old_vel(3)**2 - 2*g(3)*(old_pos(3)-pos_max(3))))/g(3)
                t_b = -2*(old_pos(3) + g(3)*t_c)/g(3)
                t_ab =  modulo(t_b,dt-t_c)
                ! bounce velocity
                particles% vel(3,i) = -(old_vel(3) + g(3)*t_c) + g(3)*(t_ab)
                particles% vel(2,i) = -particles% vel(2,i)
                particles% vel(1,i) = -particles% vel(1,i)
                ! bounce position
                particles% pos(3,i) = pos_max(3) -(old_vel(3) + g(3)*t_c)*(t_ab) + g(3)*(t_ab)**2/2
                particles% pos(2,i) = old_pos(2) + old_vel(2)*t_c + particles% vel(2,i)*t_ab
                particles% pos(1,i) = old_pos(1) + old_vel(1)*t_c + particles% vel(1,i)*t_ab
             else if (g(3)<0) then
                t_c = (-old_vel(3) + sqrt(old_vel(3)**2 - 2*g(3)*(old_pos(3)-pos_max(3))))/g(3)
                ! bounce velocity
                particles% vel(3,i) = -(old_vel(3) + g(3)*t_c) + g(3)*(dt - t_c)
                particles% vel(2,i) = -particles% vel(2,i)
                particles% vel(1,i) = -particles% vel(1,i)
                ! bounce position
                particles% pos(3,i) = pos_max(3) -(old_vel(3) + g(3)*t_c)*(dt - t_c) + g(3)*(dt - t_c)**2/2
                particles% pos(2,i) = old_pos(2) + old_vel(2)*t_c + particles% vel(2,i)*(dt - t_c)
                particles% pos(1,i) = old_pos(1) + old_vel(2)*t_c + particles% vel(1,i)*(dt - t_c)
             else ! no gravity in this direction
                ! bounce velocity
                particles% vel(:,i) = -particles% vel(:,i)
                ! particle position
                t_c = abs((pos_max(3) - old_pos(3))/old_vel(3)) 
                particles% pos(:,i) = old_pos + old_vel*t_c + particles% vel(:,i)*(dt - t_c)
             end if
          end if
       else
          particles% pos(3,i) = modulo( particles% pos(3,i) , cells% edges(3) )
       end if
    end do

  end subroutine mpcd_stream_zwall

  subroutine mpcd_stream_periodic(particles, cells, dt)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt

    integer :: i, jump(3)
    double precision :: edges(3)

    edges = cells% edges

    do i = 1, particles% Nmax
       particles% pos(:,i) = particles% pos(:,i) + particles% vel(:,i)*dt
       jump = floor(particles% pos(:,i) / edges)
       particles% image(:,i) = particles% image(:,i) + jump
       particles% pos(:,i) = particles% pos(:,i) - jump * edges
    end do

  end subroutine mpcd_stream_periodic

end module mpcd
