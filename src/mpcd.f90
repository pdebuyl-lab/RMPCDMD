module mpcd
  use common
  use particle_system
  use cell_system
  use threefry_module
  implicit none

  private

  public :: compute_temperature, simple_mpcd_step
  public :: wall_mpcd_step
  public :: mpcd_stream_periodic, mpcd_stream_zwall
  public :: mpcd_stream_xforce_yzwall
  public :: compute_rho, compute_vx

contains

  function rand_sphere(state) result(n)
    type(threefry_rng_t), intent(inout) :: state
    double precision :: n(3)

    logical :: s_lt_one
    double precision :: s, alpha

    s_lt_one = .false.
    do while (.not. s_lt_one)
       n(1) = threefry_double(state)
       n(2) = threefry_double(state)
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
    use omp_lib
    class(particle_system_t), intent(inout) :: particles
    class(cell_system_t), intent(in) :: cells
    type(threefry_rng_t), intent(inout) :: state(:)
    double precision, intent(in), optional :: temperature

    integer :: i, start, n
    integer :: cell_idx
    double precision :: local_v(3), omega(3,3), vec(3)
    logical :: thermostat
    integer :: thread_id

    thermostat = present(temperature)
    if (thermostat) error stop 'thermostatting not implemented'

    call particles%time_step%tic()
    !$omp parallel
    thread_id = omp_get_thread_num() + 1
    !$omp do private(start, n, local_v, i, vec, omega)
    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)

       local_v = 0
       do i = start, start + n - 1
          local_v = local_v + particles% vel(:, i)
       end do
       local_v = local_v / n

       vec = rand_sphere(state(thread_id))
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
    !$omp end do
    !$omp end parallel
    call particles%time_step%tac()

  end subroutine simple_mpcd_step

  !! Lamura, Gompper, Ihle and Kroll, EPL 56, 319-325 (2001)
  !! http://dx.doi.org/10.1209/epl/i2001-00522-9
  subroutine wall_mpcd_step(particles, cells, state, wall_temperature, wall_v, wall_n, thermostat, bulk_temperature)
    use hilbert
    use omp_lib
    class(particle_system_t), intent(inout) :: particles
    class(cell_system_t), intent(in) :: cells
    type(threefry_rng_t), intent(inout) :: state(:)
    double precision, optional, intent(in) :: wall_temperature(2)
    double precision, optional, intent(in) :: wall_v(3,2)
    integer, optional, intent(in) :: wall_n(2)
    logical, intent(in), optional :: thermostat
    double precision, intent(in), optional :: bulk_temperature

    integer :: i, start, n
    integer :: cell_idx
    double precision :: local_v(3), omega(3,3), vec(3)
    integer :: n_virtual
    integer :: cell(3)
    integer :: wall_idx
    double precision :: virtual_v(3), t_factor
    logical :: all_present, all_absent
    logical :: do_thermostat
    integer :: thread_id

    all_present = present(wall_temperature) .and. present(wall_v) .and. present(wall_n)
    all_absent = .not. present(wall_temperature) .and. .not. present(wall_v) .and. .not. present(wall_n)
    if ( .not. (all_present .or. all_absent) ) &
         error stop 'wall parameters must be all present or all absent in wall_mpcd_step'

    if (present(thermostat)) then
       do_thermostat = thermostat
       if (do_thermostat) then
          if (present(bulk_temperature)) then
             t_factor = sqrt(bulk_temperature)
          else
             error stop 'thermostat requested but no temperature given in wall_mpcd_step'
          end if
       end if
    else
       do_thermostat = .false.
    end if

    call particles%time_step%tic()
    !$omp parallel
    thread_id = omp_get_thread_num() + 1
    !$omp do private(start, n, n_virtual, virtual_v, cell, wall_idx, local_v, i, vec, omega)
    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)
       n_virtual = 0

       ! Find whether we are in a wall cell
       cell = compact_h_to_p(cell_idx - 1, cells% M) + 1
       if (all_present .and. (cell(3) == 1)) then
          wall_idx = 1
          do_thermostat = .false.
       else if (all_present .and. (cell(3) == cells% L(3))) then
          wall_idx = 2
          do_thermostat = .false.
       else
          wall_idx = -1
       end if

       if (wall_idx > 0) then
          if (n < wall_n(wall_idx)) then
             n_virtual = wall_n(wall_idx) - n
             virtual_v(1) = threefry_normal(state(thread_id))
             virtual_v(2) = threefry_normal(state(thread_id))
             virtual_v(3) = threefry_normal(state(thread_id))
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

       if (do_thermostat) then
          virtual_v = 0
          do i = start, start + n - 1
             vec(1) = threefry_normal(state(thread_id))
             vec(2) = threefry_normal(state(thread_id))
             vec(3) = threefry_normal(state(thread_id))
             particles% vel(:, i) = vec*t_factor
             virtual_v = virtual_v + particles% vel(:, i)
          end do
          virtual_v = local_v - virtual_v / dble(n)
          do i = start, start + n - 1
             particles% vel(:, i) = particles% vel(:, i) + virtual_v
          end do
       else
          vec = rand_sphere(state(thread_id))
          omega = &
               reshape( (/ &
               vec(1)**2, vec(1)*vec(2) + vec(3), vec(1)*vec(3) - vec(2) ,&
               vec(2)*vec(1) - vec(3) , vec(2)**2 , vec(2)*vec(3) + vec(1),&
               vec(3)*vec(1) + vec(2), vec(3)*vec(2) - vec(1), vec(3)**2 &
               /), (/3, 3/))
          do i = start, start + n - 1
             particles%vel(:,i) = local_v + matmul(omega, particles%vel(:,i)-local_v)
          end do
       end if
    end do
    !$omp end do
    !$omp end parallel
    call particles%time_step%tac()

  end subroutine wall_mpcd_step

  function compute_temperature(particles, cells, tz) result(te)
    use hilbert, only : compact_h_to_p
    type(particle_system_t), intent(inout) :: particles
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

    call particles%time_ct%tic()
    cell_idx = 1
    count = 0
    te = 0
    !$omp parallel do private(start, n, local_v, local_kin, i, cell) &
    !$omp& reduction(+:count) reduction(+:te)
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
    call particles%time_ct%tac()

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

    call particles%time_stream%tic()
    !$omp parallel do private(old_pos, old_vel, t_c, t_b, t_ab)
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
    call particles%time_stream%tac()

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

  !! Advance mpcd particles
  !!
  !! This routines allows a x-direction forcing and specular or bounce-back conditions in y
  !! and z
  subroutine mpcd_stream_xforce_yzwall(particles, cells, dt,g)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt
    double precision, intent(in):: g

    integer :: i, bc(3), im(3)
    double precision :: delta, L(3), gvec(3)
    double precision, dimension(3) :: old_pos, old_vel
    double precision, dimension(3) :: new_pos, new_vel
    double precision :: t_c, t_b, t_ab

    L = cells%edges
    bc = cells%bc
    gvec = 0
    gvec(1) = g

    call particles%time_stream%tic()
    !$omp parallel do private(old_pos, old_vel, new_pos, new_vel, im)
    do i = 1, particles% Nmax
       old_pos = particles% pos(:,i)
       old_vel = particles% vel(:,i)
       new_pos = old_pos + old_vel*dt + gvec*dt**2/2
       new_vel = old_vel + gvec*dt
       im = 0

       if ( (new_pos(2)<0) .or. (new_pos(2)>L(2)) .or. (new_pos(3)<0) .or. (new_pos(3)>L(3)) ) &
            call yzwall_collision(old_pos, old_vel, new_pos, new_vel, im, L, dt, bc, g)

       im = floor(new_pos/L)
       new_pos = new_pos - im*L
       particles%pos(:,i) = new_pos
       particles%vel(:,i) = new_vel
       particles%image(:,i) = particles%image(:,i) + im
    end do
    call particles%time_stream%tac()

  end subroutine mpcd_stream_xforce_yzwall

  !! Collide a particle in y and z with constant acceleration in x
  !!
  !! Periodic boundary conditions are applied in x
  subroutine yzwall_collision(x0, v0, x, v, im, L, t, bc, g)
    double precision, dimension(3), intent(inout) :: x0, v0, x, v
    integer, dimension(3), intent(inout) :: im
    double precision, intent(in) :: L(3), t
    integer, intent(in) :: bc(3)
    double precision, intent(in), optional :: g

    double precision, dimension(3) :: gvec
    double precision :: t_collision, t_remainder, tt
    integer :: i, coll_dim, actual_bc
    integer :: jump(3)

    gvec = 0
    if (present(g)) gvec(1) = g

    coll_dim = 0
    t_collision = huge(t_collision)
    do i = 2, 3
       if (v(i) > 0) then
          tt = (L(i)-x0(i))/v0(i)
          if ( (tt>0) .and. (tt<t_collision) ) then
             t_collision = tt
             coll_dim = i
          end if
       else
          tt = -x0(i)/v0(i)
          if ( (tt>0) .and. (tt<t_collision) ) then
             t_collision = tt
             coll_dim = i
          end if
       end if
    end do

    if (coll_dim==0) return

    t_remainder = t - t_collision

    x = x0 + v0*t_collision + gvec*t_collision**2 / 2
    v = v0 + gvec*t_collision

    if (bc(coll_dim) == BOUNCE_BACK_BC) then
       v = -v
    else if (bc(coll_dim) == SPECULAR_BC) then
       v(coll_dim) = -v(coll_dim)
    else if (bc(coll_dim) == PERIODIC_BC) then
       jump = floor(x/L)
       x(coll_dim) = x(coll_dim) - jump(coll_dim)*L(coll_dim)
    end if

    wall_loop: do while (.true.)
       coll_dim = change_23(coll_dim)

       if (v(coll_dim)>0) then
          t_collision = (L(coll_dim)-x(coll_dim))/v(coll_dim)
       else
          t_collision = -x(coll_dim)/v(coll_dim)
       end if
       if ((t_collision < 0) .or. (t_collision > t_remainder)) exit wall_loop
       x = x + v*t_collision + gvec*t_collision**2 / 2
       v = v + gvec*t_collision
       t_remainder = t_remainder - t_collision
       if (bc(coll_dim) == BOUNCE_BACK_BC) then
          v = -v
       else if (bc(coll_dim) == SPECULAR_BC) then
          v(coll_dim) = -v(coll_dim)
       else if (bc(coll_dim) == PERIODIC_BC) then
          jump = floor(x/L)
          x(coll_dim) = x(coll_dim) - jump(coll_dim)*L(coll_dim)
       end if
    end do wall_loop

    t_collision = t_remainder
    x = x + v*t_collision + gvec*t_collision**2 / 2
    v = v + gvec*t_collision

  end subroutine yzwall_collision

  pure function change_23(i) result(r)
    integer, intent(in) :: i
    integer :: r
    if (i==2) then
       r = 3
    else if (i==3) then
       r = 2
    end if
  end function change_23

end module mpcd
