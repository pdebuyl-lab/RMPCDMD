!> Routines to perform MPCD dynamics
!!
!! \manual{algorithms}.
!!
!!
!! MPCD collisions are implemented in simple_mpcd_step and wall_mpcd_step, that takes into
!! account a wall in the z-direction.
!!
!! Streaming is achieved by the mpcd_stream_periodic, mpcd_stream_zwall,
!! mpcd_stream_xforce_yzwall and mpcd_stream_nogravity_zwall routines. Only a single stream
!! routine should be used, depending on the simulation setup. A further call to md_vel is
!! needed in the presence of forces.

module mpcd
  use common
  use particle_system
  use cell_system
  use threefry_module
  implicit none

  private

  public :: compute_temperature, simple_mpcd_step
  public :: wall_mpcd_step
  public :: mpcd_at_step
  public :: mpcd_stream_periodic, mpcd_stream_zwall
  public :: mpcd_stream_xforce_yzwall
  public :: compute_rho, compute_vx
  public :: bulk_reaction
  public :: mpcd_stream_nogravity_zwall

contains

  !> Return random point on the surface of a sphere
  !!
  !! Ref. \cite marsaglia_random_sphere_1972
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

  !> Perform a collision.
  !!
  !! Use the rule defined in Ref. \cite malevanets_kapral_mpcd_1999 to collide the particles
  !! cell-wise. \manual{algorithms,mpcd}
  subroutine simple_mpcd_step(particles, cells, state, alpha)
    use omp_lib
    class(particle_system_t), intent(inout) :: particles
    class(cell_system_t), intent(in) :: cells
    type(threefry_rng_t), intent(inout) :: state(:)
    double precision, intent(in), optional :: alpha

    integer :: i, start, n
    integer :: cell_idx
    double precision :: local_v(3), omega(3,3), vec(3)
    double precision :: sin_alpha, cos_alpha
    integer :: thread_id

    if (present(alpha)) then
       sin_alpha = sin(alpha)
       cos_alpha = cos(alpha)
    else
       sin_alpha = 1
       cos_alpha = 0
    end if

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
            cos_alpha+vec(1)**2*(1-cos_alpha), vec(1)*vec(2)*(1-cos_alpha) - vec(3)*sin_alpha, &
            vec(1)*vec(3)*(1-cos_alpha) + vec(2)*sin_alpha ,&
            vec(2)*vec(1)*(1-cos_alpha) + vec(3)*sin_alpha , cos_alpha+vec(2)**2*(1-cos_alpha) ,&
            vec(2)*vec(3)*(1-cos_alpha) - vec(1)*sin_alpha,&
            vec(3)*vec(1)*(1-cos_alpha) - vec(2)*sin_alpha, vec(3)*vec(2)*(1-cos_alpha) + vec(1)*sin_alpha,&
            cos_alpha + vec(3)**2*(1-cos_alpha) &
            /), (/3, 3/))

       do i = start, start + n - 1
          particles% vel(:, i) = local_v + matmul(omega, (particles% vel(:, i)-local_v))
       end do

    end do
    !$omp end do
    !$omp end parallel
    call particles%time_step%tac()

  end subroutine simple_mpcd_step

  !> Perform a collision.
  !!
  !! MPCD with Anderson thermostat defined in Ref. \cite noguchi_epl_2007 to collide the
  !! particles cell-wise. \manual{algorithms,mpcd}
  subroutine mpcd_at_step(particles, cells, state, temperature)
    use omp_lib
    class(particle_system_t), intent(inout) :: particles
    class(cell_system_t), intent(in) :: cells
    type(threefry_rng_t), intent(inout) :: state(:)
    double precision, intent(in) :: temperature

    integer :: i, start, n
    integer :: cell_idx
    double precision :: local_v(3), vec(3), virtual_v(3)
    double precision :: t_factor
    integer :: thread_id

    t_factor = sqrt(temperature)

    call particles%time_step%tic()
    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(cell_idx, start, n, local_v, virtual_v, i, vec)
    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)

       local_v = 0
       do i = start, start + n - 1
          local_v = local_v + particles% vel(:, i)
       end do
       local_v = local_v / n

       virtual_v = 0
       do i = start, start + n - 1
          vec(1) = threefry_normal(state(thread_id))*t_factor
          vec(2) = threefry_normal(state(thread_id))*t_factor
          vec(3) = threefry_normal(state(thread_id))*t_factor
          particles% vel(:, i) = vec
          virtual_v = virtual_v + particles% vel(:, i)
       end do
       virtual_v = local_v - virtual_v / dble(n)
       do i = start, start + n - 1
          particles% vel(:, i) = particles% vel(:, i) + virtual_v
       end do

    end do
    !$omp end do
    !$omp end parallel
    call particles%time_step%tac()

  end subroutine mpcd_at_step

  !> Collisions in partially filled cells at the walls use the rule of
  !! Lamura, Gompper, Ihle and Kroll, EPL 56, 319-325 (2001)
  !! http://dx.doi.org/10.1209/epl/i2001-00522-9
  subroutine wall_mpcd_step(particles, cells, state, wall_temperature, wall_v, wall_n, &
       thermostat, bulk_temperature, alpha)
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
    double precision, intent(in), optional :: alpha

    integer :: i, start, n
    integer :: cell_idx
    double precision :: local_v(3), omega(3,3), vec(3)
    double precision :: sin_alpha, cos_alpha
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

    if (present(alpha)) then
       sin_alpha = sin(alpha)
       cos_alpha = cos(alpha)
    else
       sin_alpha = 1
       cos_alpha = 0
    end if

    call particles%time_step%tic()
    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(cell_idx, start, n, n_virtual, virtual_v, cell, wall_idx, local_v, i, vec, omega)
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

       if ( (do_thermostat) .and. (wall_idx<0) ) then
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
               cos_alpha+vec(1)**2*(1-cos_alpha), vec(1)*vec(2)*(1-cos_alpha) - vec(3)*sin_alpha, &
               vec(1)*vec(3)*(1-cos_alpha) + vec(2)*sin_alpha ,&
               vec(2)*vec(1)*(1-cos_alpha) + vec(3)*sin_alpha , cos_alpha+vec(2)**2*(1-cos_alpha) ,&
               vec(2)*vec(3)*(1-cos_alpha) - vec(1)*sin_alpha,&
               vec(3)*vec(1)*(1-cos_alpha) - vec(2)*sin_alpha, vec(3)*vec(2)*(1-cos_alpha) + vec(1)*sin_alpha,&
               cos_alpha + vec(3)**2*(1-cos_alpha) &
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

  !> Compute the temperature of a mpcd solvent
  !!
  !! \manual{algorithms,temperature-computation}
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

  !> Compute density profile along z
  subroutine compute_rho(particles, rhoz)
    type(particle_system_t), intent(in) :: particles
    type(histogram_t), intent(inout) :: rhoz

    integer :: i

    if (.not. allocated(rhoz% data)) error stop 'histogram_t: data not allocated'

    do i = 1, particles% Nmax
       call rhoz% bin(particles% pos(3, i))
    end do

  end subroutine compute_rho

  !> Compute x-velocity profile along z
  subroutine compute_vx(particles, vx)
    type(particle_system_t), intent(in) :: particles
    type(profile_t), intent(inout) :: vx

    integer :: i

    if (.not. allocated(vx% data)) error stop 'histogram_t: data not allocated'

    do i = 1, particles% Nmax
       call vx% bin(particles% pos(3, i), particles% vel(1, i))
    end do

  end subroutine compute_vx

  !> Advance mpcd particles
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

  !> Stream MPCD particles in a periodic system
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

  !> Stream MPCD particles with a force in the x-direction and specular or bounce-back
  !! conditions in y
  !!
  !! MPCD particles near the walls must not be in the interaction range of a colloid, this
  !! is not verified programmatically.
  subroutine mpcd_stream_xforce_yzwall(particles, cells, dt, g)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt
    double precision, intent(in):: g

    integer :: i, bc(3), im(3)
    double precision :: delta, L(3), gvec(3)
    double precision, dimension(3) :: old_pos, old_vel
    double precision, dimension(3) :: new_pos, new_vel
    double precision :: t_c, t_b, t_ab
    logical :: y_out, z_out

    L = cells%edges
    bc = cells%bc
    gvec = 0
    gvec(1) = g

    call particles%time_stream%tic()
    !$omp parallel do private(old_pos, old_vel, new_pos, new_vel, im, y_out, z_out)
    do i = 1, particles% Nmax
       old_pos = particles% pos(:,i)
       old_vel = particles% vel(:,i)
       new_pos = old_pos + old_vel*dt + particles%force(:,i)*dt**2/2
       new_vel = old_vel
       im = 0
       y_out = ((new_pos(2)<0) .or. (new_pos(2)>L(2)))
       z_out = ((new_pos(3)<0) .or. (new_pos(3)>L(3)))

       if ( &
            ((bc(2)==PERIODIC_BC) .and. z_out ) .or. &
            ((bc(2)/=PERIODIC_BC) .and. (y_out .or. z_out)) &
            ) then
          call yzwall_collision(old_pos, old_vel, new_pos, new_vel, im, L, dt, bc, g)
          particles%wall_flag(i) = 1
       end if

       particles%pos(:,i) = new_pos
       particles%vel(:,i) = new_vel
    end do
    call particles%time_stream%tac()

  end subroutine mpcd_stream_xforce_yzwall

  !> Collide a particle in y and z with constant acceleration in x
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
    im = 0

    coll_dim = 0
    t_collision = huge(t_collision)
    do i = 2, 3
       if (v0(i) > 0) then
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
       im(coll_dim) = im(coll_dim) + jump(coll_dim)
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
          im(coll_dim) = im(coll_dim) + jump(coll_dim)
       end if
    end do wall_loop

    t_collision = t_remainder
    x = x + v*t_collision + gvec*t_collision**2 / 2 + im*L
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

  subroutine bulk_reaction(p, c, from, to, rate, tau, state)
    use omp_lib
    type(particle_system_t), intent(inout) :: p
    type(cell_system_t), intent(inout) :: c
    integer, intent(in) :: from, to
    double precision, intent(in) :: rate, tau
    type(threefry_rng_t), intent(inout) :: state(:)

    integer :: cell_idx, start, n, i, pick, s
    double precision :: local_rate
    integer :: thread_id

    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(start, n, local_rate, i, pick, s)
    do cell_idx = 1, c%N
       if ( (c%cell_count(cell_idx) <= 1) .or. .not. c%is_reac(cell_idx) ) cycle

       start = c%cell_start(cell_idx)
       n = c%cell_count(cell_idx)

       local_rate = 0
       do i = start, start + n - 1
          s = p%species(i)
          if (s==from) then
             local_rate = local_rate + 1
             pick = i
          end if
       end do

       local_rate = local_rate*rate
       if (threefry_double(state(thread_id)) < (1 - exp(-local_rate*tau))) then
          p%species(pick) = to
       end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine bulk_reaction

  subroutine mpcd_stream_nogravity_zwall(particles, cells, dt)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt

    integer :: i
    double precision :: delta, L(3)
    double precision, dimension(3) :: old_pos, old_vel
    double precision, dimension(3) :: new_pos, new_vel
    double precision :: t_c, t_b, t_ab
    logical :: bounce

    L = cells%edges

    call particles%time_stream%tic()
    !$omp parallel do private(i, old_pos, old_vel, new_pos, new_vel, bounce, t_c)
    particles_loop: do i = 1, particles% Nmax
       old_pos = particles%pos(:,i)
       old_vel = particles%vel(:,i)
       new_pos = old_pos + old_vel*dt + particles%force(:,i)*dt**2 / 2
       bounce = .false.

       if (new_pos(3)<0) then
          t_c = -old_pos(3)/old_vel(3)
          bounce = .true.
       else if (new_pos(3)>L(3)) then
          t_c = (L(3)-old_pos(3))/old_vel(3)
          bounce = .true.
       else
          new_vel = old_vel
       end if

       if (bounce) then
          new_pos = new_pos - old_vel*2*(dt-t_c)
          new_vel = - old_vel
          particles%wall_flag(i) = 1
       end if

       particles%pos(:,i) = new_pos
       particles%vel(:,i) = new_vel
    end do particles_loop
    call particles%time_stream%tac()

  end subroutine mpcd_stream_nogravity_zwall

end module mpcd
