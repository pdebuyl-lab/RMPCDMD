!> Model a chemotactic experiment in a microfluidic channel
!!
!! In this simulation, an inlet (x=0) is fed with A and S fluid species in the lower and
!! upper halves in the y direction, respectively. A constant accerelation is applied in the
!! x direction and walls in the z direction confine the flow, leading to a Poiseuille
!! velocity profile.
!!
!! The colloid is a passive sphere, an active sphere or a dimer nanomotor.
!!
!! \param g                magnitude of acceleration
!! \param buffer_length    length of the inlet buffer
!! \param max_speed        maximum velocity of profile to initialize the velocities
!! \param prob             probability of reaction
!! \param alpha            angle of collision
!! \param store_rho_xy     store the xy density of solvent particles on a grid
!! \param dimer            simulate a dimer nanomotor (boolean, else it is a single sphere)
!! \param L                length of simulation box in the 3 dimensions
!! \param rho              fluid number density
!! \param T                Temperature. Used for setting initial velocities and for wall thermostatting.
!! \param d                length of rigid link
!! \param N_in_front       place N sphere in front (higher x), for the dimer nanomotor
!! \param tau              MPCD collision time
!! \param N_MD_steps       number MD steps occuring in tau
!! \param N_loop           number of MPCD timesteps
!! \param colloid_sampling interval (in MD steps) of sampling the colloid position and velocity
!! \param steps_fixed      number of steps during which the colloid is fixed (only when buffer_length>0)
!! \param equilibration_loops number of MPCD steps for equilibration (only when buffer_length=0)
!! \param sigma_C          radius of C sphere
!! \param sigma_N          radius of N sphere
!! \param track_y_shift    shift of the track in the y direction with respect to Ly/2
!! \param epsilon_C        interaction parameter of C sphere with both solvent species (2 elements)
!! \param epsilon_N        interaction parameter of N sphere with both solvent species (2 elements)

program chemotactic_cell
  use rmpcdmd_module
  use hdf5
  use h5md_module
  use threefry_module
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  type(threefry_rng_t), allocatable :: state(:)

  integer, parameter :: N_species = 3

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj
  type(lj_params_t) :: walls_colloid_lj

  type(profile_t) :: vx

  integer :: rho
  integer :: N
  integer :: error
  integer :: N_colloids
  integer, parameter :: n_bins_conc = 90
  double precision :: conc_z_cyl(n_bins_conc)

  double precision :: sigma_N, sigma_C, max_cut, alpha, sigma_sphere
  double precision :: shift
  double precision :: track_y_shift
  double precision :: sigma(3,2), sigma_cut(3,2), epsilon(3,2)
  double precision,allocatable :: mass(:)

  double precision :: v_com(3), wall_v(3,2), wall_t(2)
  double precision :: local_mass, total_mass

  double precision :: e1, e2, e_wall
  double precision :: tau, dt , T
  double precision :: d,prob
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: colloid_sampling
  integer :: n_extra_sorting
  double precision :: kin_e, temperature
  integer, dimension(N_species) :: n_solvent, catalytic_change, bulk_change
  type(h5md_element_t) :: n_solvent_el, catalytic_change_el, bulk_change_el
  type(h5md_element_t) :: omega_el

  double precision :: colloid_pos(3,2)
  double precision :: com_pos(3)
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group, params_group
  type(h5md_element_t) :: rho_xy_el
  type(thermo_t) :: thermo_data
  type(particle_system_io_t) :: dimer_io
  type(particle_system_io_t) :: solvent_io
  integer(HID_T) :: box_group
  type(h5md_element_t) :: elem_vx, elem_vx_count

  type(PTo) :: config
  integer :: i, L(3),  n_threads
  integer :: j, k, m
  integer :: i_release, i_block

  type(timer_t), target :: flag_timer, change_timer, buffer_timer, varia
  double precision :: total_time
  type(timer_list_t) :: timer_list
  integer(HID_T) :: timers_group

  integer, allocatable :: rho_xy(:,:,:)

  integer, parameter :: block_length = 8
  type(axial_correlator_t) :: axial_cf
  type(correlator_t) :: omega_acf
  integer(HID_T) :: correlator_group

  double precision :: unit_r(3), omega(3), rel_v(3), norm_xy

  double precision :: g(3) !gravity
  logical :: fixed, on_track, stopped, N_in_front, dimer, N_type
  logical :: store_rho_xy
  integer :: buffer_length
  logical :: sampling
  integer :: equilibration_loops
  double precision :: max_speed, z, Lz
  integer :: steps_fixed
  type(args_t) :: args

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  call flag_timer%init('flag')
  call change_timer%init('change')
  call buffer_timer%init('buffer')
  call varia%init('varia')

  call timer_list%init(13)
  call timer_list%append(flag_timer)
  call timer_list%append(change_timer)
  call timer_list%append(buffer_timer)
  call timer_list%append(varia)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call h5open_f(error)

  call hfile%create(args%output_file, 'RMPCDMD::chemotactic_cell', 'N/A', 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)

  g = 0
  g(1) = PTread_d(config, 'g', loc=params_group)
  buffer_length = PTread_i(config, 'buffer_length', loc=params_group)
  max_speed = PTread_d(config,'max_speed', loc=params_group)
  prob = PTread_d(config,'probability', loc=params_group)
  alpha = PTread_d(config,'alpha', loc=params_group)
  store_rho_xy = PTread_l(config, 'store_rho_xy', loc=params_group)
  dimer = PTread_l(config, 'dimer', loc=params_group)
  N_type = PTread_l(config, 'N_type', loc=params_group)
  if (dimer) then
     N_colloids = 2
  else
     N_colloids = 1
  end if
  L = PTread_ivec(config, 'L', 3, loc=params_group)
  L(1) = L(1) + buffer_length

  rho = PTread_i(config, 'rho', loc=params_group)
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T', loc=params_group)
  d = PTread_d(config, 'd', loc=params_group)
  N_in_front = PTread_l(config, 'N_in_front', loc=params_group)
  track_y_shift = PTread_d(config, 'track_y_shift', loc=params_group)

  wall_v = 0
  wall_t = [T, T]

  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  if (modulo(N_MD_steps, colloid_sampling) /= 0) then
     error stop 'colloid_sampling must divide N_MD with no remainder'
  end if
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)
  steps_fixed = PTread_i(config, 'steps_fixed', loc=params_group)
  equilibration_loops = PTread_i(config, 'equilibration_loops', loc=params_group)

  sigma_C = PTread_d(config, 'sigma_C', loc=params_group)
  sigma_N = PTread_d(config, 'sigma_N', loc=params_group)

  epsilon(:,1) = PTread_dvec(config, 'epsilon_C', N_species, loc=params_group)
  epsilon(:,2) = PTread_dvec(config, 'epsilon_N', N_species, loc=params_group)
  sigma(:,1) = sigma_C
  sigma(:,2) = sigma_N

  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)
  epsilon = 1.d0
  sigma(1,1) = 2*sigma_C
  sigma(1,2) = sigma_C + sigma_N
  sigma(2,1) = sigma_C + sigma_N
  sigma(2,2) = 2*sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)
  call colloid_lj% init(epsilon(1:2,:), sigma(1:2,:), sigma_cut(1:2,:))


  epsilon = 1.d0
  sigma(1,:) = [sigma_C, sigma_N]
  sigma_cut = sigma*3**(1.d0/6.d0)
  shift = max(sigma_C, sigma_N)*2**(1./6.) + 0.25
  call walls_colloid_lj% init(epsilon(1:1,:), sigma(1:1,:), sigma_cut(1:1,:), shift)

  if (N_type) then
     sigma_sphere = sigma_N
  else
     sigma_sphere = sigma_C
  end if

  allocate(mass(N_colloids))
  if (dimer) then
     mass(1) = rho * sigma_C**3 * 4 * 3.14159265/3
     mass(2) = rho * sigma_N**3 * 4 * 3.14159265/3
  else
     mass = rho * sigma_sphere**3 * 4 * 3.14159265/3
  end if

  call solvent% init(N,N_species, system_name='solvent')

  call colloids% init(N_colloids,2, mass, system_name='colloids') !there will be 2 species of colloids

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call PTkill(config)

  call axial_cf%init(block_length, N_loop*N_MD_steps/colloid_sampling, N_loop*N_MD_steps)
  call omega_acf%init(block_length, get_n_blocks(block_length, n_blocks_max=7, &
       n_samples=N_loop*N_MD_steps/colloid_sampling))

  if (dimer) then
     colloids% species(1) = 1
     colloids% species(2) = 2
  else
     if (N_type) then
        colloids% species = 2
     else
        colloids% species = 1
     end if
  end if
  colloids% vel = 0

  dimer_io%force_info%store = .false.
  dimer_io%id_info%store = .false.
  dimer_io%position_info%store = .true.
  dimer_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%position_info%step = colloid_sampling
  dimer_io%position_info%time = colloid_sampling*dt
  dimer_io%image_info%store = .true.
  dimer_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%image_info%step = colloid_sampling
  dimer_io%image_info%time = colloid_sampling*dt
  dimer_io%velocity_info%store = .true.
  dimer_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%velocity_info%step = colloid_sampling
  dimer_io%velocity_info%time = colloid_sampling*dt
  dimer_io%species_info%store = .true.
  dimer_io%species_info%mode = H5MD_FIXED
  call dimer_io%init(hfile, 'dimer', colloids)

  solvent_io%force_info%store = .false.
  solvent_io%id_info%store = .false.
  solvent_io%position_info%store = .true.
  solvent_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  solvent_io%position_info%step = N_loop*N_MD_steps
  solvent_io%position_info%time = N_loop*N_MD_steps*dt
  solvent_io%image_info%store = .true.
  solvent_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  solvent_io%image_info%step = N_loop*N_MD_steps
  solvent_io%image_info%time = N_loop*N_MD_steps*dt
  solvent_io%velocity_info%store = .true.
  solvent_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  solvent_io%velocity_info%step = N_loop*N_MD_steps
  solvent_io%velocity_info%time = N_loop*N_MD_steps*dt
  solvent_io%species_info%store = .true.
  solvent_io%species_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  solvent_io%species_info%step = N_loop*N_MD_steps
  solvent_io%species_info%time = N_loop*N_MD_steps*dt
  call solvent_io%init(hfile, 'solvent', solvent)

  solvent% force = 0
  solvent% species = 1
  call solvent_cells%init(L, 1.d0,has_walls = .true.)
  call vx% init(0.d0, solvent_cells% edges(3), L(3))

  if (store_rho_xy) allocate(rho_xy(N_species, L(2), L(1)))
  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  if (store_rho_xy)  call rho_xy_el%create_time(fields_group, 'rho_xy', rho_xy, ior(H5MD_LINEAR,H5MD_STORE_TIME), &
       step=N_MD_steps, time=N_MD_steps*dt, offset_by_one=.true.)
  call elem_vx% create_time(fields_group, 'vx', vx% data, ior(H5MD_TIME, H5MD_STORE_TIME), offset_by_one=.true.)
  call elem_vx_count% create_time(fields_group, 'vx_count', vx% count, ior(H5MD_TIME, H5MD_STORE_TIME), offset_by_one=.true.)

  call h5gclose_f(fields_group, error)

  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, ior(H5MD_LINEAR,H5MD_STORE_TIME), step=N_MD_steps, &
       time=N_MD_steps*dt, offset_by_one=.true.)
  call catalytic_change_el%create_time(hfile%observables, 'catalytic_change', &
       catalytic_change, ior(H5MD_LINEAR,H5MD_STORE_TIME), step=N_MD_steps, &
       time=N_MD_steps*dt, offset_by_one=.true.)
  call bulk_change_el%create_time(hfile%observables, 'bulk_change', &
       bulk_change, ior(H5MD_LINEAR,H5MD_STORE_TIME), step=N_MD_steps, &
       time=N_MD_steps*dt, offset_by_one=.true.)
  call omega_el%create_time(hfile%observables, 'omega', &
       omega(3), ior(H5MD_LINEAR,H5MD_STORE_TIME), step=colloid_sampling, &
       time=colloid_sampling*dt, offset_by_one=.true.)

  if (dimer) then
     colloids% pos(3,:) = solvent_cells% edges(3)/2.d0
     if (N_in_front) then
        colloids% pos(1,1) = sigma_C*2**(1.d0/6.d0) + 1
        colloids% pos(1,2) = colloids% pos(1,1) + d
        colloids% pos(2,:) = solvent_cells% edges(2)/2.d0 + track_y_shift
     else
        colloids% pos(1,2) = sigma_N*2**(1.d0/6.d0) + 1
        colloids% pos(1,1) = colloids% pos(1,2) + d
        colloids% pos(2,:) = solvent_cells% edges(2)/2.d0 + track_y_shift
     end if
  else
     colloids% pos(3,:) = solvent_cells% edges(3)/2.d0
     colloids% pos(2,:) = solvent_cells% edges(2)/2.d0 + track_y_shift
     colloids% pos(1,:) = sigma_sphere*2**(1.d0/6.d0) + 1
  end if

  call h5gcreate_f(dimer_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  Lz = solvent_cells%edges(3)
  do i=1, solvent% Nmax
     z = solvent%pos(3,i)
     solvent% vel(1,i) = threefry_normal(state(1))*sqrt(T) + &
          max_speed*z*(Lz-z)/(Lz/2)**2
     solvent% vel(2,i) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,i) = threefry_normal(state(1))*sqrt(T)
  end do

  do m = 1, solvent% Nmax
     if (solvent% pos(2,m) < (L(2)/2.d0)) then
        solvent% species(m) = 1
     else
        solvent% species(m) = 3
     end if
  end do

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, 10*int(300*max(sigma_C,sigma_N)**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

  solvent% force(2:3,:) = 0
  solvent% force(1,:) = g(1)
  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj, 3)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force
  catalytic_change = 0

  i = 0
  if (buffer_length>0) then
     fixed = .true.
     on_track = .true.
  else
     fixed = .false.
     on_track = .false.
  end if
  stopped = .false.

  solvent_cells%bc = [PERIODIC_BC, SPECULAR_BC, BOUNCE_BACK_BC]

  sampling = .false.
  i_release = 0
  i = 0
  i_block = 0
  write(*,*) 'Running for', N_loop, 'loops'
  !start RMPCDMD
  setup: do while (.not. stopped)
     if (modulo(i,20) == 0) write(*,'(i09)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        call mpcd_stream_xforce_yzwall(solvent, solvent_cells, dt, g(1))

        colloids% pos_rattle = colloids% pos

        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(k))
        end do

        if (dimer) then
           call rattle_dimer_pos(colloids, d, dt, solvent_cells% edges)
        end if

        if ((.not. on_track) .and. (buffer_length/=0))then
           do k=1, colloids% Nmax
              if (colloids% pos(1,k) > solvent_cells% edges(1)) then
                 stopped = .true.
              end if
           end do
        end if

        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1) .or. (so_max >= skin*0.9) ) then
           call varia%tic()
           call apply_pbc(colloids, solvent_cells% edges)
           call apply_pbc(solvent, solvent_cells% edges)
           call varia%tac()
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, max_cut + skin, solvent_cells)
           call varia%tic()
           solvent% pos_old = solvent% pos
           colloids% pos_old = colloids% pos
           call varia%tac()
           n_extra_sorting = n_extra_sorting + 1
        end if

        call buffer_particles(solvent,solvent_cells%edges)

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)

        !$omp parallel do
        do k = 1, solvent%Nmax
           solvent% force(1,k) = g(1)
           solvent% force(2:3,k) = 0
        end do
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
        if (.not. on_track) then
           e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj, 3)
        end if
        if (on_track) then
           colloids% force(2,:) = 0
           colloids% force(3,:) = 0
           if (fixed) then
              colloids% force(1,:) = 0
           end if
        end if

        call md_vel(solvent, dt)

        do k=1, colloids% Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
                dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(k))
        end do
        if (dimer) then
           call rattle_dimer_vel(colloids, d, dt, solvent_cells% edges)
        end if

        if (.not.fixed) then
           ! this should be solved for a single passive colloid
           if ((.not. N_type) .or. (dimer)) then
              if (.not. on_track) then
                 call flag_particles
                 call change_species
              end if
           end if
        end if

        com_pos = (colloids%pos(:,1)+colloids%image(:,1)*solvent_cells%edges + &
             colloids%pos(:,2)+colloids%image(:,2)*solvent_cells%edges)
        unit_r = rel_pos(colloids%pos(:,1), colloids%pos(:,2), solvent_cells%edges)
        norm_xy = norm2(unit_r(1:2))
        rel_v = colloids%vel(:,1)-colloids%vel(:,2)
        omega = cross(unit_r, rel_v) / norm_xy**2
        unit_r = unit_r / norm2(unit_r)

        if (sampling) then
           v_com = sum(colloids%vel, dim=2)/2
           call axial_cf%add_fast((i-i_release)*N_MD_steps+j-1, v_com, unit_r)
        end if

        if ((sampling) .and. (modulo(j, colloid_sampling)==0)) then
           ! correlators
           call axial_cf%add(i_block, com_pos, unit_r)
           call omega_acf%add(i_block, correlate_block_dot, x=omega(3))
           i_block = i_block + 1
           ! colloid trajectory
           call dimer_io%position%append(colloids%pos)
           call dimer_io%velocity%append(colloids%vel)
           call dimer_io%image%append(colloids%image)
           call omega_el%append(omega(3))
        end if

     end do md_loop

     call solvent_cells%random_shift(state(1))

     call apply_pbc(colloids, solvent_cells% edges)
     call apply_pbc(solvent, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)
     solvent% pos_old = solvent% pos
     colloids% pos_old = colloids% pos

     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[rho, rho], alpha=alpha)

     call compute_vx(solvent, vx)
     if ((sampling) .and. (modulo(i-i_release+1, 50) == 0)) then
        call vx% norm()
        call elem_vx% append(vx% data, i-i_release+1, (i-i_release+1)*tau)
        call elem_vx_count% append(vx% count, i-i_release+1, (i-i_release+1)*tau)
        call vx% reset()
     end if

     call varia%tic()
     if ((sampling) .and. (store_rho_xy)) then
        call compute_rho_xy
        call rho_xy_el%append(rho_xy)
     end if
     call varia%tac()

     temperature = compute_temperature(solvent, solvent_cells)
     total_mass = 0
     kin_e = sum(solvent% vel**2)/2
     v_com = sum(solvent% vel, dim=2)
     do k = 1, colloids% Nmax
        m = colloids%species(k)
        if (m==0) cycle
        local_mass = colloids%mass(m)
        kin_e = kin_e + local_mass*sum(colloids% vel(:,k)**2)/2
        v_com = v_com + local_mass*colloids%vel(:,k)
        total_mass = total_mass + local_mass
     end do
     v_com = v_com / (solvent%Nmax + total_mass)

     n_solvent = 0
     do k = 1, solvent%Nmax
        m = solvent%species(k)
        if (m <= 0) cycle
        n_solvent(m) = n_solvent(m) + 1
     end do

     if (sampling) then
        call thermo_data%append(hfile, temperature, e1+e2+e_wall, kin_e, e1+e2+e_wall+kin_e, v_com)
        call n_solvent_el%append(n_solvent)
        call catalytic_change_el%append(catalytic_change)
        call bulk_change_el%append(bulk_change)
     end if

     call h5fflush_f(hfile%id, H5F_SCOPE_GLOBAL_F, error)

     if (fixed) then
        if (i >= steps_fixed) then
           write(*,*) 'fixed', fixed
           fixed = .false.
        end if
     end if

     if (on_track) then
        if (dimer) then
           if ((colloids% pos(1,1) > (buffer_length+sigma_C)) &
                .and. (colloids% pos(1,2) > (buffer_length+sigma_N))) then
              on_track = .false.
           end if
        else
           if (colloids% pos(1,1) > (buffer_length+sigma_sphere)) then
              on_track = .false.
           end if
        end if
     end if

     i = i+1
     if ( &
          ((.not. sampling) .and. (buffer_length == 0) .and. (i >= equilibration_loops)) .or. &
          ((i_release==0) .and. (buffer_length > 0) .and. (.not. on_track)) ) then
        i_release = i
        sampling = .true.
        write(*,*) 'i_release =', i_release
     end if
     if (i-i_release > N_loop) exit setup
  end do setup

  call thermo_data%append(hfile, temperature, e1+e2+e_wall, kin_e, e1+e2+e_wall+kin_e, v_com, add=.false., force=.true.)

  write(*,*) 'n extra sorting', n_extra_sorting

  ! create a group for block correlators and write the data

  call h5gcreate_f(hfile%id, 'block_correlators', correlator_group, error)

  call axial_cf%write(correlator_group, colloid_sampling, colloid_sampling*dt, 1, dt)

  call write_correlator_block(correlator_group, 'planar_angular_velocity_autocorrelation', &
       omega_acf, colloid_sampling, colloid_sampling*dt)

  call h5gclose_f(correlator_group, error)

  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

  call h5gcreate_f(hfile%id, 'timers', timers_group, error)
  call timer_list%append(solvent%time_stream)
  call timer_list%append(solvent%time_md_vel)
  call timer_list%append(solvent%time_step)
  call timer_list%append(solvent%time_count)
  call timer_list%append(solvent%time_sort)
  call timer_list%append(solvent%time_ct)
  call timer_list%append(solvent%time_max_disp)
  call timer_list%append(neigh%time_update)
  call timer_list%append(neigh%time_force)

  call timer_list%write(timers_group, total_time)

  call h5md_write_dataset(timers_group, 'total', total_time)

  call h5gclose_f(timers_group, error)

  if (store_rho_xy) call rho_xy_el%close()
  call elem_vx% close()
  call elem_vx_count% close()
  call dimer_io%close()
  call hfile%close()
  call h5close_f(error)

contains

  subroutine flag_particles
    double precision :: dist_to_C_sq
    integer :: r, s
    double precision :: x(3)

    call flag_timer%tic()
    do s = 1,neigh% n(1)
       r = neigh%list(s,1)
       if (solvent% species(r) == 1) then
          x = rel_pos(colloids% pos(:,1),solvent% pos(:,r),solvent_cells% edges)
          dist_to_C_sq = dot_product(x, x)
          if (dist_to_C_sq < solvent_colloid_lj%cut_sq(1,1)) then
             if (threefry_double(state(1)) <= prob) then
                solvent% flag(r) = 1
             end if
          end if
       end if
    end do
    call flag_timer%tac()

  end subroutine flag_particles

  subroutine change_species
    double precision :: dist_to_C_sq
    double precision :: dist_to_N_sq
    integer :: m
    double precision :: x(3)

    call change_timer%tic()
    catalytic_change = 0
    !$omp parallel do private(x, dist_to_C_sq, dist_to_N_sq) reduction(+:catalytic_change)
    do m = 1, solvent% Nmax
       if (solvent% flag(m) == 1) then
          if (dimer) then
             x = rel_pos(colloids% pos(:,1), solvent% pos(:,m), solvent_cells% edges)
             dist_to_C_sq = dot_product(x, x)
             x = rel_pos(colloids% pos(:,2), solvent% pos(:,m), solvent_cells% edges)
             dist_to_N_sq = dot_product(x, x)
             if ( &
                (dist_to_C_sq > solvent_colloid_lj%cut_sq(1,1)) &
                .and. &
                (dist_to_N_sq > solvent_colloid_lj%cut_sq(1,2)) &
                ) &
                then
                solvent% species(m) = 2
                solvent% flag(m) = 0
                catalytic_change(1) = catalytic_change(1) - 1
                catalytic_change(2) = catalytic_change(2) + 1
             end if
          else
             x = rel_pos(colloids% pos(:,1), solvent% pos(:,m), solvent_cells% edges)
             dist_to_C_sq = dot_product(x, x)
             if (dist_to_C_sq > solvent_colloid_lj%cut_sq(1,1)) then
                solvent% species(m) = 2
                solvent% flag(m) = 0
                catalytic_change(1) = catalytic_change(1) - 1
                catalytic_change(2) = catalytic_change(2) + 1
             end if
          end if
       end if
    end do
    call change_timer%tac()

  end subroutine change_species

  subroutine concentration_field_cylindrical
    double precision :: dimer_orient(3),x(3),y(3),z(3)
    double precision :: solvent_pos(3,solvent% Nmax)
    double precision :: dz,r,theta,x_pos,y_pos,z_pos
    integer :: o
    integer :: check
    logical :: far_enough_from_wall
    double precision :: range_min1(3),range_min2(3),range_max1(3),range_max2(3)

    dz = 2.d0*d/n_bins_conc
    dimer_orient = colloids% pos(:,2) - colloids% pos(:,1)
    z = dimer_orient/sqrt(dot_product(dimer_orient,dimer_orient))

    x = (/0.d0, 1.d0, -dimer_orient(2)/dimer_orient(3)/)
    x = x/sqrt(dot_product(x,x))
    y = (/z(2)*x(3)-z(3)*x(2),z(3)*x(1)-z(1)*x(3),z(1)*x(2)-z(2)*x(1)/)
    conc_z_cyl = 0

    range_min1 = colloids%pos(:,1) - d/2.0*z - (/0.d0,0.d0,1.d0/)*2*max_cut
    range_min2 = colloids%pos(:,1) - d/2.0*z + (/0.d0,0.d0,1.d0/)*2*max_cut
    range_max1 = colloids%pos(:,1) + 3.d0*d/2.0*z - (/0.d0,0.d0,1.d0/)*2*max_cut
    range_max2 = colloids%pos(:,1) - 3.d0*d/2.0*z + (/0.d0,0.d0,1.d0/)*2*max_cut

    if ( (range_min1(3)<solvent_cells%edges(3)).and.(range_min1(3)>0).and. &
       (range_max1(3)<solvent_cells%edges(3)).and.(range_max1(3)>0).and. &
       (range_min2(3)<solvent_cells%edges(3)).and.(range_min2(3)>0).and. &
       (range_max2(3)<solvent_cells%edges(3)).and.(range_max2(3)>0) ) then
       far_enough_from_wall = .true.
    else
       far_enough_from_wall = .false.
    end if
    if (far_enough_from_wall) then
       do o = 1, solvent% Nmax
          solvent_pos(:,o) = solvent% pos(:,o) - colloids% pos(:,1)
          x_pos = dot_product(x,solvent_pos(:,o))
          y_pos = dot_product(y, solvent_pos(:,o))
          z_pos = dot_product(z, solvent_pos(:,o))
          solvent_pos(:,o) = (/x_pos,y_pos,z_pos/)
       end do
       do o = 1, solvent% Nmax
          r = sqrt(solvent_pos(1,o)**2 + solvent_pos(2,o)**2)
          theta = atan(solvent_pos(2,o)/solvent_pos(1,o))
          solvent_pos(1,o) = r
          solvent_pos(2,o) = theta
          if ((solvent_pos(1,o) < 2*max_cut).and.(solvent_pos(3,o)<1.5d0*d).and.(solvent_pos(3,o)>-0.5d0*d)) then
             if (solvent% species(o)==2) then
                check = floor((solvent_pos(3,o)+0.5d0*d)/dz)
                check = check+1
                conc_z_cyl(check) = conc_z_cyl(check) + 1
             end if
          end if
       end do
       colloid_pos(:,1) = 0
       colloid_pos(3,1) = colloids% pos(3,1)
       colloid_pos(:,2) = 0
       colloid_pos(3,2) = d + colloids% pos(3,1)
    else
       conc_z_cyl = 0
       colloid_pos = 0
    end if
  end subroutine concentration_field_cylindrical

  subroutine buffer_particles(particles,edges)
     type(particle_system_t), intent(inout) :: particles
     double precision, intent(in) :: edges(3)

     integer :: k, s

     call buffer_timer%tic()
     bulk_change = 0
     !$omp parallel do private(s) reduction(+:bulk_change)
     do k = 1, particles% Nmax
        s = particles% species(k)
        if (s <= 0) continue
        if ((particles% pos(1,k) > 0) .and. (particles% pos(1,k) < buffer_length)) then
           if (particles% pos(2,k) < edges(2)/2.d0) then
              bulk_change(s) = bulk_change(s) - 1
              particles% species(k) = 1
              s = 1
              bulk_change(s) = bulk_change(s) + 1
           else
              bulk_change(s) = bulk_change(s) - 1
              particles% species(k) = 3
              s = 3
              bulk_change(s) = bulk_change(s) + 1
           end if
        end if
     end do
     call buffer_timer%tac()

  end subroutine buffer_particles

  subroutine compute_rho_xy
    integer :: i, s, ix, iy

    rho_xy = 0
    do i = 1, solvent%Nmax
       s = solvent%species(i)
       if (s <= 0) continue
       ix = modulo(floor(solvent%pos(1,i)/solvent_cells%a), L(1)) + 1
       iy = modulo(floor(solvent%pos(2,i)/solvent_cells%a), L(2)) + 1
       rho_xy(s, iy, ix) = rho_xy(s, iy, ix) + 1
    end do

  end subroutine compute_rho_xy

end program chemotactic_cell
