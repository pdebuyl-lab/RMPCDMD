! This file is part of RMPCDMD
! Copyright (c) 2015-2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Simulate a single colloidal rigid-body particle
!!
!! This simulation models a chemically active colloid particle in either a periodic
!! simulation box or with confinement in the y direction.
!!
!! The coordinates of the colloid particle's beads must be provided in a H5MD file, as a
!! "fixed-in-time" dataset. The body of the particle can operate as a rigid-body (with
!! either RATTLE or quaternion dynamics) or as an elastic network.
!!
!! \param L                     length of simulation box in the 3 dimensions
!! \param rho                   fluid number density
!! \param T                     Temperature. Used for setting initial velocities and (if enabled) bulk thermostatting.
!! \param tau                   MPCD collision time
!! \param alpha                 angle of collision
!! \param probability           probability to change A to B upon collision
!! \param bulk_rate             rate of B->A reaction
!! \param N_MD                  number MD steps occuring in tau
!! \param N_loop                number of MPCD timesteps
!! \param colloid_sampling      interval (in MD steps) of sampling the colloid position and velocity
!! \param do_solvent_io         if true (T), a snapshot of the solvent in the final step is dump to the datafile
!! \param equilibration_loops   number of MPCD steps for equilibration
!! \param epsilon_C             interaction parameter of C sphere with both solvent species (2 elements)
!! \param epsilon_N             interaction parameter of N sphere with both solvent species (2 elements)
!! \param data_filename         filename for input Janus coordinates
!! \param data_group            particles group in the input file
!! \param epsilon_colloid       interaction parameter for colloid-colloid interactions
!! \param reaction_radius       radius for the reaction around the Janus particle
!! \param link_treshold         distance criterion for finding rigid-body links
!! \param do_read_links         read link information from a file
!! \param links_file            filename for the links data
!! \param do_rattle             perform RATTLE
!! \param do_lennard_jones      compute colloid-colloid Lennard-Jones forces
!! \param do_elastic            compute colloid-colloid elastic network forces
!! \param elastic_k             elastic constant for the elastic network
!! \param rattle_pos_tolerance  absolute tolerance for Rattle (position part)
!! \param rattle_vel_tolerance  absolute tolerance for Rattle (velocity part)
!! \param do_quaternion         perform quaternion velocity Verlet
!! \param quaternion_treshold   treshold for the iterative procedure for the quaternion integrator
!! \param sigma                 radius of the colloidal beads for colloid-solvent interactions
!! \param sigma_colloid         radius of the colloidal beads for colloid-colloid interactions
!! \param polar_r_max           maximal radius for the polar fields
!! \param do_ywall              use a confining potential in the y direction, 9-3 Lennard-Jones
!! \param wall_sigma            wall LJ sigma
!! \param wall_epsilon          wall LJ epsilon
!! \param wall_shift            wall shift
!! \param fluid_wall            boundary condition for the fluid

program single_body
  use rmpcdmd_module
  use hdf5
  use h5md_module
  use threefry_module
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj
  type(lj_params_t) :: walls_colloid_lj

  integer, parameter :: N_species = 2

  integer :: rho
  integer :: N
  integer :: error

  double precision :: max_cut
  double precision :: epsilon(2,2)
  double precision :: sigma, sigma_v(2,2), sigma_cut(2,2)
  double precision :: mass(2)
  double precision :: wall_sigma(3, N_species), wall_epsilon(3, N_species), wall_shift(3)

  double precision :: e1, e2, e3, e_wall
  double precision :: tau, dt , T
  double precision :: mpcd_alpha
  double precision :: prob, reaction_radius
  double precision :: bulk_rate
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting
  double precision :: loop_i_last_sort

  type(PTo) :: config

  integer :: i, L(3)
  integer :: j, k, cell_idx

  double precision :: total_time
  type(timer_t), target :: varia, main, time_flag, time_refuel, time_change
  type(timer_t), target :: so_timer
  type(timer_t), target :: bulk_reac_timer
  type(timer_t), target :: q_timer
  type(timer_list_t) :: timer_list
  integer(HID_T) :: timers_group

  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads
  type(h5md_file_t) :: hfile
  type(h5md_file_t) :: input_data_file
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group, params_group
  integer(HID_T) :: connectivity_group
  integer(HID_T) :: box_group
  integer(HID_T) :: tmp_id
  integer(HID_T) :: pos_dset
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  double precision :: v_com(3), x(3)
  double precision :: wall_v(3,2)
  double precision :: tmp_mass
  type(particle_system_io_t) :: janus_io
  type(particle_system_io_t) :: solvent_io

  integer, dimension(N_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el

  type(h5md_element_t) :: q_el, omega_body_el, janus_pos_el, janus_vel_el, u_el
  double precision, dimension(3) :: one_x, one_y, one_z

  double precision :: polar_r_max
  type(polar_fields_t) :: polar
  type(planar_fields_t) :: planar
  integer(HID_T) :: polar_id
  integer(HID_T) :: planar_id

  logical :: do_rattle, do_read_links, do_lennard_jones, do_elastic, do_quaternion, do_solvent_io
  logical :: do_ywall
  integer, allocatable :: links(:,:)
  double precision, allocatable :: links_d(:)
  double precision :: link_treshold
  integer :: i_link
  double precision :: dist, rattle_pos_tolerance, rattle_vel_tolerance
  double precision :: elastic_k
  double precision :: unit_r(3)
  double precision :: com_pos(3)
  double precision :: alpha, beta, z0
  type(rigid_body_t) :: rigid_janus
  double precision :: quaternion_treshold
  character(len=144) :: fluid_wall

  integer, parameter :: block_length = 8
  type(axial_correlator_t) :: axial_cf
  type(correlator_t) :: omega_cf
  integer(HID_T) :: correlator_group

  integer :: equilibration_loops
  integer :: colloid_sampling, coordinates_sampling
  logical :: sampling

  type(args_t) :: args
  character(len=144) :: data_filename
  character(len=144) :: links_file
  character(len=144) :: data_group
  logical :: attr_exists

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call main%init('main')
  call varia%init('varia')
  call so_timer%init('so_timer')
  call time_flag%init('flag')
  call time_refuel%init('refuel')
  call time_change%init('change')
  call bulk_reac_timer%init('bulk_reac')
  call q_timer%init('quaternion_vv')

  call timer_list%init(26)
  call timer_list%append(varia)
  call timer_list%append(so_timer)
  call timer_list%append(time_flag)
  call timer_list%append(time_refuel)
  call timer_list%append(time_change)
  call timer_list%append(bulk_reac_timer)
  call timer_list%append(q_timer)

  call h5open_f(error)
  call hfile%create(args%output_file, 'RMPCDMD::single_body', &
       RMPCDMD_REVISION, 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)


  prob = PTread_d(config,'probability', loc=params_group)
  bulk_rate = PTread_d(config,'bulk_rate', loc=params_group)

  L = PTread_ivec(config, 'L', 3, loc=params_group)
  rho = PTread_i(config, 'rho', loc=params_group)
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T', loc=params_group)
  
  tau = PTread_d(config, 'tau', loc=params_group)
  mpcd_alpha = PTread_d(config,'alpha', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  coordinates_sampling = PTread_i(config, 'coordinates_sampling', loc=params_group)
  do_solvent_io = PTread_l(config, 'do_solvent_io', loc=params_group)
  if (modulo(N_MD_steps, colloid_sampling) /= 0) then
     error stop 'colloid_sampling must divide N_MD with no remainder'
  end if
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)
  equilibration_loops = PTread_i(config, 'equilibration_loops', loc=params_group)

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  do_rattle = PTread_l(config, 'do_rattle', loc=params_group)
  do_lennard_jones = PTread_l(config, 'do_lennard_jones', loc=params_group)
  rattle_pos_tolerance = PTread_d(config, 'rattle_pos_tolerance', loc=params_group)
  rattle_vel_tolerance = PTread_d(config, 'rattle_vel_tolerance', loc=params_group)

  do_quaternion  = PTread_l(config, 'do_quaternion', loc=params_group)
  quaternion_treshold = PTread_d(config, 'quaternion_treshold', loc=params_group)

  do_elastic = PTread_l(config, 'do_elastic', loc=params_group)
  if (do_elastic) elastic_k = PTread_d(config, 'elastic_k', loc=params_group)

  sigma = PTread_d(config, 'sigma_colloid', loc=params_group)
  sigma_v = sigma
  sigma_cut = sigma_v*2**(1.d0/6.d0)
  mass = rho * sigma**3 * 4 * 3.14159265/3

  epsilon = PTread_d(config, 'epsilon_colloid', loc=params_group)

  call colloid_lj% init(epsilon, sigma_v, sigma_cut)

  ! no wall in x and z
  wall_sigma = -1
  wall_shift = 0

  do_ywall = PTread_l(config, 'do_ywall', loc=params_group)
  fluid_wall = 'PERIODIC'
  if (do_ywall) then
     wall_sigma(3,:) = PTread_d(config, 'wall_sigma', loc=params_group)
     wall_shift(3) = PTread_d(config, 'wall_shift', loc=params_group)
     wall_epsilon = PTread_d(config, 'wall_epsilon', loc=params_group)
     call walls_colloid_lj% init(wall_epsilon, &
          wall_sigma, 3.d0**(1.d0/6.d0)*wall_sigma, wall_shift)
     fluid_wall = PTread_s(config, 'fluid_wall', loc=params_group)
  end if

  ! solvent index first, colloid index second, in solvent_colloid_lj
  sigma = PTread_d(config, 'sigma', loc=params_group)
  sigma_v = sigma
  sigma_cut = sigma_v*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  epsilon(:,1) = PTread_dvec(config, 'epsilon_C', 2, loc=params_group)
  epsilon(:,2) = PTread_dvec(config, 'epsilon_N', 2, loc=params_group)

  call solvent_colloid_lj% init(epsilon, sigma_v, sigma_cut)

  call solvent% init(N, N_species, system_name='solvent') !there will be 2 species of solvent particles

  data_filename = PTread_s(config, 'data_filename', loc=params_group)
  data_group = PTread_s(config, 'data_group', loc=params_group)
  link_treshold = PTread_d(config, 'link_treshold', loc=params_group)
  do_read_links = PTread_l(config, 'do_read_links', loc=params_group)
  if (do_read_links) links_file = PTread_s(config, 'links_file', loc=params_group)

  polar_r_max = PTread_d(config, 'polar_r_max', loc=params_group)

  reaction_radius = PTread_d(config, 'reaction_radius', loc=params_group)

  call h5gclose_f(params_group, error)
  call PTkill(config)

  call solvent_cells%init(L, 1.d0, &
       has_walls= ( (trim(fluid_wall) == 'SPECULAR') &
       .or. (trim(fluid_wall) == 'BOUNCE_BACK')) )

  wall_v = 0
  solvent_cells%bc = PERIODIC_BC
  if (trim(fluid_wall) == 'PERIODIC') then
     solvent_cells%bc(3) = PERIODIC_BC
  else if (trim(fluid_wall) == 'SPECULAR') then
     solvent_cells%bc(3) = SPECULAR_BC
  else if (trim(fluid_wall) == 'BOUNCE_BACK') then
     solvent_cells%bc(3) = BOUNCE_BACK_BC
  else
     error stop 'unknown value for parameter fluid_wall'
  end if

  call axial_cf%init(block_length, N_loop, N_loop*N_MD_steps)

  call colloids%init_from_file(data_filename, data_group, H5MD_FIXED)

  ! center colloid in simulation box
  com_pos = sum(colloids%pos, dim=2) / size(colloids%pos, dim=2)
  do i = 1, size(colloids%pos, dim=2)
     colloids%pos(:,i) = colloids%pos(:,i) + (solvent_cells%edges/2 - com_pos)
  end do

  call input_data_file%open(data_filename, H5F_ACC_RDONLY_F)

  call h5aexists_by_name_f(input_data_file%particles, trim(data_group)//'/position', 'alpha', attr_exists, error)
  ! If the parameters alpha,beta,gamma are present, read them
  ! Only alpha is checked, if someone puts in alpha and not the rest, long hdf5 error log will show
  if (attr_exists) then
     call h5oopen_f(input_data_file%particles, trim(data_group)//'/position', pos_dset, error)
     call h5md_read_attribute(pos_dset, 'alpha', alpha)
     call h5md_read_attribute(pos_dset, 'beta', beta)
     call h5md_read_attribute(pos_dset, 'z0', z0)
     call h5oclose_f(pos_dset, error)
  end if

  call input_data_file%close()

  colloids%image = 0
  colloids%vel = 0
  colloids%force = 0
  colloids%mass = mass

  if (do_read_links) then
     call read_links(links_file, i_link, links, links_d)
  else
     i_link = 0
     do i = 1, colloids%Nmax
        do j = i+1, colloids%Nmax
           x = rel_pos(colloids%pos(:,i), colloids%pos(:,j), solvent_cells%edges)
           dist = norm2(x)
           if (dist < link_treshold) then
              ! count link
              i_link = i_link + 1
           end if
        end do
     end do
     allocate(links(2,i_link), links_d(i_link))

     i_link = 0
     do i = 1, colloids%Nmax
        do j = i+1, colloids%Nmax
           x = rel_pos(colloids%pos(:,i), colloids%pos(:,j), solvent_cells%edges)
           dist = norm2(x)
           if (dist < link_treshold) then
              ! add link
              i_link = i_link + 1
              links(1, i_link) = i
              links(2, i_link) = j
              links_d(i_link) = dist
           end if
        end do
     end do
  end if

  call rigid_janus%init(colloids, 1, colloids%Nmax, solvent_cells%edges)

  write(*,*) 'number of links:', i_link

  janus_io%force_info%store = .true.
  janus_io%force_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%force_info%step = N_MD_steps*coordinates_sampling
  janus_io%force_info%time = N_MD_steps*coordinates_sampling*dt
  janus_io%id_info%store = .false.
  janus_io%position_info%store = .true.
  janus_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%position_info%step = janus_io%force_info%step
  janus_io%position_info%time = janus_io%force_info%time
  janus_io%image_info%store = .true.
  janus_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%image_info%step = janus_io%force_info%step
  janus_io%image_info%time = janus_io%force_info%time
  janus_io%velocity_info%store = .true.
  janus_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%velocity_info%step = janus_io%force_info%step
  janus_io%velocity_info%time = janus_io%force_info%time
  janus_io%species_info%store = .true.
  janus_io%species_info%mode = H5MD_FIXED
  call janus_io%init(hfile, 'janus', colloids)

  if (do_solvent_io) then
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
  end if

  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, H5MD_LINEAR, step=N_MD_steps, &
       time=N_MD_steps*dt)

  if (do_quaternion) then
     call q_el%create_time(hfile%observables, 'q', &
          rigid_janus%q, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=colloid_sampling, time=colloid_sampling*dt)
     call omega_body_el%create_time(hfile%observables, 'omega_body', &
          rigid_janus%omega_body, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=colloid_sampling, time=colloid_sampling*dt)
     call omega_cf%init(block_length, get_n_blocks(block_length, 8, N_loop), dim=3)
  end if

  call u_el%create_time(hfile%observables, 'u', &
       unit_r, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=colloid_sampling, &
       time=colloid_sampling*dt)

  call janus_pos_el%create_time(hfile%observables, 'janus_pos', &
       rigid_janus%pos, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=colloid_sampling, &
       time=colloid_sampling*dt)

  call janus_vel_el%create_time(hfile%observables, 'janus_vel', &
       rigid_janus%vel, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=colloid_sampling, &
       time=colloid_sampling*dt)

  call h5gcreate_f(janus_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  if (do_solvent_io) then
     call h5gcreate_f(solvent_io%group, 'box', box_group, error)
     call h5md_write_attribute(box_group, 'dimension', 3)
     call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
     call h5gclose_f(box_group, error)
  end if

  call h5gcreate_f(hfile%id, 'fields', fields_group, error)

  call h5gcreate_f(hfile%id, 'connectivity', connectivity_group, error)
  call h5md_write_dataset(connectivity_group, 'janus_links', links-1)
  call h5dopen_f(connectivity_group, 'janus_links', tmp_id, error)
  call h5md_write_attribute(tmp_id, 'particles_group', 'janus')
  call h5dclose_f(tmp_id, error)
  call h5gclose_f(connectivity_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  call polar%init(N_species, 64, sigma, polar_r_max, 64)
  call timer_list%append(polar%time_polar_update)
  call planar%init(N_species, 64, -12.d0, 12.d0, 64, -12.d0, 12.d0, 1.d0)
  call timer_list%append(planar%timer)
  call neigh% init(colloids% Nmax, int(500*max(sigma,1.15d0)**3))

  skin = 1.5
  n_extra_sorting = 0
  loop_i_last_sort = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  if (do_lennard_jones) then
     e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  else
     e2 = 0
  end if
  if (do_elastic) then
     e3 = elastic_network(colloids, links, links_d, elastic_k, solvent_cells%edges)
  else
     e3 = 0
  end if
  if (do_ywall) then
     e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj)
  else
     e_wall = 0
  end if

  solvent% force_old = solvent% force
  colloids% force_old = colloids% force
  if (do_quaternion) call rigid_janus%compute_force_torque(colloids)

  write(*,*) 'Running for', equilibration_loops, '+', N_loop, 'loops'
  call main%tic()
  sampling = .false.
  do i = 0, N_loop+equilibration_loops
     if (i==equilibration_loops) sampling = .true.
     if (modulo(i,32) == 0) write(*,'(i08)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        if ((do_ywall) .and. (solvent_cells%bc(3)/=PERIODIC_BC)) then
           call cell_md_pos_zwall(solvent_cells, solvent, dt, md_flag=.true.)
        else
           call cell_md_pos(solvent_cells, solvent, dt, md_flag=.true.)
        end if

        if (do_rattle) then
           call varia%tic()
           ! Extra copy for rattle
           colloids% pos_rattle = colloids% pos
           do k=1, colloids% Nmax
              colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                   dt**2 * colloids% force(:,k) / (2 * colloids% mass(colloids%species(k)))
           end do
           call varia%tac()

           call rattle_body_pos(colloids, links, links_d, dt, solvent_cells% edges, rattle_pos_tolerance)

           com_pos = modulo(sum((colloids%pos(:,:)+ &
                colloids%image(:,:)*spread(solvent_cells%edges, dim=2, ncopies=colloids%Nmax) ) &
                , dim=2)/ colloids%Nmax, solvent_cells%edges)

        else if (do_quaternion) then
           call q_timer%tic()
           call rigid_janus%vv1(colloids, dt, quaternion_treshold, solvent_cells%edges)
           com_pos = modulo(rigid_janus%pos, solvent_cells%edges)
           call q_timer%tac()
        end if

        so_max = cell_maximum_displacement(solvent_cells, solvent, delta_t=dt*(N_MD_steps*i+j - loop_i_last_sort))
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1d0) .or. (so_max >= skin*0.9d0) ) then
           if ((do_ywall) .and. (solvent_cells%bc(3)/=PERIODIC_BC)) then
              call cell_md_pos_zwall(solvent_cells, solvent, (N_MD_steps*i+j - loop_i_last_sort)*dt, md_flag=.false.)
           else
              call cell_md_pos(solvent_cells, solvent, (N_MD_steps*i+j - loop_i_last_sort)*dt, md_flag=.false.)
           end if
           call cell_md_vel(solvent_cells, solvent, (N_MD_steps*i+j - loop_i_last_sort)*dt, md_flag=.false.)

           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call solvent% sort(solvent_cells)
           loop_i_last_sort = N_MD_steps*i + j
           call neigh% update_list(colloids, solvent, max_cut + skin, solvent_cells, solvent_colloid_lj)
           call so_timer%tic()
           !$omp parallel do
           do k = 1, solvent%Nmax
              solvent% pos_old(:,k) = solvent% pos(:,k)
           end do
           call compute_cell_wise_max_v
           call so_timer%tac()
           call varia%tic()
           colloids% pos_old = colloids% pos
           n_extra_sorting = n_extra_sorting + 1
           call varia%tac()
        end if

        call varia%tic()
        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)
        call varia%tac()

        call so_timer%tic()
        !$omp parallel do private(cell_idx, k)
        do cell_idx = 1, solvent_cells%N
           if (solvent_cells%is_md(cell_idx)) then
              do k = solvent_cells%cell_start(cell_idx), &
                   solvent_cells%cell_start(cell_idx) + solvent_cells%cell_count(cell_idx) - 1
                 solvent%force(:,k) = 0
              end do
           end if
        end do
        call so_timer%tac()
        colloids% force = 0

        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        if (do_lennard_jones) e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
        if (do_elastic) e3 = elastic_network(colloids, links, links_d, elastic_k, solvent_cells%edges)
        if (do_ywall) e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj)

        call cell_md_vel(solvent_cells, solvent, dt, md_flag=.true.)

        if (do_quaternion) then
           call q_timer%tic()
           call rigid_janus%compute_force_torque(colloids)
           call rigid_janus%vv2(colloids, dt)
           call q_timer%tac()
        else if (do_rattle) then
           call varia%tic()
           do k = 1, colloids%Nmax
              colloids% vel(:,k) = colloids% vel(:,k) + &
                   dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids%mass(colloids%species(k)))
           end do
           call varia%tac()
        end if

        if (do_rattle) then
           call rattle_body_vel(colloids, links, links_d, dt, solvent_cells% edges, rattle_vel_tolerance)
        end if

        if (sampling) then
           v_com = sum(colloids%vel, dim=2)/colloids%Nmax
           com_pos = sum((colloids%pos(:,:)+ &
                colloids%image(:,:)*spread(solvent_cells%edges, dim=2, ncopies=colloids%Nmax) ) &
                , dim=2)/ colloids%Nmax
           if (do_rattle .and. (.not. do_quaternion)) then
              rigid_janus%pos = com_pos
              rigid_janus%vel = v_com
           end if
           unit_r = get_unit_r()
           call axial_cf%add_fast((i-equilibration_loops)*N_MD_steps+j-1, v_com, unit_r)
        end if

        if ((sampling) .and. (modulo(j, colloid_sampling)==0)) then
           call u_el%append(unit_r)
           call janus_pos_el%append(rigid_janus%pos)
           call janus_vel_el%append(rigid_janus%vel)
           if (do_quaternion) then
              call q_el%append(rigid_janus%q)
              call omega_body_el%append(rigid_janus%omega_body)
           end if
        end if

        call time_flag%tic()
        call flag_particles
        call time_flag%tac()
        call time_change%tic()
        call change_species
        call time_change%tac()

     end do md_loop

     if (sampling) then
        temperature = compute_temperature(solvent, solvent_cells)
        kin_e = 0
        v_com = 0
        call so_timer%tic()
        !$omp parallel do private(k) reduction(+:kin_e) reduction(+:v_com)
        do k = 1, solvent%Nmax
           kin_e = kin_e + sum(solvent%vel(:,k)**2)
           v_com = v_com + solvent%vel(:,k)
        end do
        call so_timer%tac()
        call varia%tic()
        tmp_mass = solvent%Nmax
        do k = 1, colloids%Nmax
           v_com = v_com + colloids%mass(colloids%species(k)) * colloids%vel(:,k)
           tmp_mass = tmp_mass + colloids%mass(colloids%species(k))
           kin_e = kin_e + colloids%mass(colloids%species(k)) * sum(colloids%vel(:,k)**2)
        end do
        v_com = v_com / tmp_mass
        kin_e = kin_e / 2

        call thermo_data%append(hfile, temperature, e1+e2+e3+e_wall, kin_e, e1+e2+e3+e_wall+kin_e, v_com)
        call axial_cf%add(i-equilibration_loops, com_pos, unit_r)
        call varia%tac()
     end if

     call solvent_cells%random_shift(state(1))

     call cell_md_pos(solvent_cells, solvent, ((i+1)*N_MD_steps - loop_i_last_sort)*dt, md_flag=.false.)
     call cell_md_vel(solvent_cells, solvent, ((i+1)*N_MD_steps - loop_i_last_sort)*dt, md_flag=.false.)

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     loop_i_last_sort = N_MD_steps*(i+1)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)
     !$omp parallel do
     do k = 1, solvent%Nmax
        solvent% pos_old(:,k) = solvent% pos(:,k)
     end do
     colloids% pos_old = colloids% pos

     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=[T, T], wall_v=wall_v, wall_n=[rho, rho], &
          alpha=mpcd_alpha)
     call compute_cell_wise_max_v

     call bulk_reac_timer%tic()
     call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate, tau, state)
     call bulk_reac_timer%tac()

     if (sampling) then
        call so_timer%tic()
        n_solvent = 0
        !$omp parallel do private(k,j) reduction(+:n_solvent)
        do k = 1, solvent%Nmax
           j = solvent%species(k)
           if (j <= 0) cycle
           n_solvent(j) = n_solvent(j) + 1
        end do
        call so_timer%tac()
        call n_solvent_el%append(n_solvent)

        call varia%tic()
        v_com = sum(colloids%vel, dim=2)/colloids%Nmax
        call varia%tac()
        call polar%update(com_pos, v_com, unit_r/norm2(unit_r), solvent, solvent_cells)
        one_x = qrot(rigid_janus%q, [1.d0, 0.d0, 0.d0])
        one_y = qrot(rigid_janus%q, [0.d0, 1.d0, 0.d0])
        one_z = qrot(rigid_janus%q, [0.d0, 0.d0, 1.d0])
        call planar%update(com_pos, v_com, one_x, one_y, one_z, qrot(rigid_janus%q, rigid_janus%omega_body), &
             solvent, solvent_cells)

        if (do_quaternion) call omega_cf%add(i-equilibration_loops, correlate_block_dot, xvec=rigid_janus%omega_body)

        if ((i > equilibration_loops) .and. (modulo(i-equilibration_loops, coordinates_sampling)==0)) then
           call janus_io%position%append(colloids%pos)
           call janus_io%force%append(colloids%force)
           call janus_io%velocity%append(colloids%vel)
           call janus_io%image%append(colloids%image)
        end if
    end if

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  ! create a group for block correlators and write the data

  call h5gcreate_f(hfile%id, 'block_correlators', correlator_group, error)
  call axial_cf%write(correlator_group, N_MD_steps, N_MD_steps*dt, 1, dt)
  if (do_quaternion) &
       call write_correlator_block(correlator_group, 'omega_body_autocorrelation', omega_cf, N_MD_steps, N_MD_steps*dt)
  call h5gclose_f(correlator_group, error)

  ! write solvent coordinates for last step

  call thermo_data%append(hfile, temperature, e1+e2+e3+e_wall, kin_e, e1+e2+e3+e_wall+kin_e, v_com, add=.false., force=.true.)
  if (do_solvent_io) then
     call solvent_io%position%append(solvent%pos)
     call solvent_io%velocity%append(solvent%vel)
     call solvent_io%image%append(solvent%image)
     call solvent_io%species%append(solvent%species)
     call solvent_io%close()
  end if

  ! store polar fields

  where (polar%count>0)
     polar%v(1,:,:,:) = polar%v(1,:,:,:) / polar%count
     polar%v(2,:,:,:) = polar%v(2,:,:,:) / polar%count
  end where
  call h5md_write_dataset(fields_group, 'polar_concentration', dble(polar%count)/N_loop)
  call h5oopen_f(fields_group, 'polar_concentration', polar_id, error)
  call h5md_write_attribute(polar_id, 'r_min', polar%r_min)
  call h5md_write_attribute(polar_id, 'r_max', polar%r_max)
  call h5md_write_attribute(polar_id, 'dr', polar%dr)
  call h5md_write_attribute(polar_id, 'dtheta', polar%dtheta)
  call h5oclose_f(polar_id, error)

  call h5md_write_dataset(fields_group, 'polar_velocity', polar%v)
  call h5oopen_f(fields_group, 'polar_velocity', polar_id, error)
  call h5md_write_attribute(polar_id, 'r_min', polar%r_min)
  call h5md_write_attribute(polar_id, 'r_max', polar%r_max)
  call h5md_write_attribute(polar_id, 'dr', polar%dr)
  call h5md_write_attribute(polar_id, 'dtheta', polar%dtheta)
  call h5oclose_f(polar_id, error)

  ! store planar fields

  where (planar%count>0)
     planar%v(1,:,:,:) = planar%v(1,:,:,:) / planar%count
     planar%v(2,:,:,:) = planar%v(2,:,:,:) / planar%count
  end where
  call h5md_write_dataset(fields_group, 'planar_concentration', dble(planar%count)/N_loop)
  call h5oopen_f(fields_group, 'planar_concentration', planar_id, error)
  call h5md_write_attribute(planar_id, 'x_min', planar%x_min)
  call h5md_write_attribute(planar_id, 'dx', planar%dx)
  call h5md_write_attribute(planar_id, 'y_min', planar%y_min)
  call h5md_write_attribute(planar_id, 'dy', planar%dy)
  call h5md_write_attribute(planar_id, 'thickness', planar%thickness)
  call h5oclose_f(planar_id, error)

  call h5md_write_dataset(fields_group, 'planar_velocity', planar%v)
  call h5oopen_f(fields_group, 'planar_velocity', planar_id, error)
  call h5md_write_attribute(planar_id, 'x_min', planar%x_min)
  call h5md_write_attribute(planar_id, 'dx', planar%dx)
  call h5md_write_attribute(planar_id, 'y_min', planar%y_min)
  call h5md_write_attribute(planar_id, 'dy', planar%dy)
  call h5md_write_attribute(planar_id, 'thickness', planar%thickness)
  call h5oclose_f(planar_id, error)

  ! Store timing data
  call timer_list%append(solvent%time_stream)
  call timer_list%append(solvent%time_step)
  call timer_list%append(solvent%time_count)
  call timer_list%append(solvent%time_sort)
  call timer_list%append(solvent%time_ct)
  call timer_list%append(solvent%time_md_pos)
  call timer_list%append(solvent%time_md_vel)
  call timer_list%append(solvent%time_max_disp)
  call timer_list%append(solvent%time_apply_pbc)
  call timer_list%append(colloids%time_self_force)
  call timer_list%append(colloids%time_rattle_pos)
  call timer_list%append(colloids%time_rattle_vel)
  call timer_list%append(colloids%time_elastic)
  call timer_list%append(colloids%time_apply_pbc)
  call timer_list%append(neigh%time_update)
  call timer_list%append(neigh%time_force)
  call timer_list%append(colloids%time_max_disp)

  call h5gcreate_f(hfile%id, 'timers', timers_group, error)
  call timer_list%write(timers_group, total_time)
  call h5md_write_dataset(timers_group, 'total', total_time)
  call h5md_write_dataset(timers_group, main%name, main%total)
  call h5gclose_f(timers_group, error)

  call h5gclose_f(fields_group, error)
  call janus_io%close()
  call n_solvent_el%close()
  if (do_quaternion) then
     call q_el%close()
     call omega_body_el%close()
  end if
  call u_el%close()
  call janus_pos_el%close()
  call janus_vel_el%close()
  call hfile%close()
  call h5close_f(error)

contains

  subroutine flag_particles
    double precision :: dist_to_C_sq
    integer :: i, r, s
    double precision :: x(3)

    !$omp parallel
    !$omp do private(i, s, r, x, dist_to_C_sq)
    do i = 1, colloids%Nmax
       if (colloids%species(i)==1) then
          do s = 1,neigh% n(i)
             r = neigh%list(s,i)
             if (solvent% species(r) == 1) then
                x = rel_pos(colloids% pos(:,i),solvent% pos(:,r),solvent_cells% edges)
                dist_to_C_sq = dot_product(x, x)
                if (dist_to_C_sq < solvent_colloid_lj%cut_sq(1,1)) then
                   !$omp atomic
                   solvent%flags(r) = ior(solvent%flags(r), REAC_MASK)
                end if
             end if
          end do
       end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine flag_particles

  subroutine change_species
    integer :: m, thread_id, i
    double precision :: dist

    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(i, m, dist)
    do i = 1, solvent_cells%N
       change_loop: do m = solvent_cells%cell_start(i), solvent_cells%cell_start(i) + solvent_cells%cell_count(i) - 1
          if (solvent%species(m) /= 1) cycle change_loop
          if ((btest(solvent%flags(m), REAC_BIT)) .and. (.not. btest(solvent%flags(m), MD_BIT))) then
             dist = norm2(rel_pos(modulo(rigid_janus%pos, solvent_cells%edges), solvent%pos(:,m), solvent_cells%edges))
             if (dist > reaction_radius) then
                if (threefry_double(state(thread_id)) < prob) solvent%species(m) = 2
                solvent%flags(m) = ibclr(solvent%flags(m), REAC_BIT)
             end if
          end if
       end do change_loop
    end do
    !$omp end do
    !$omp end parallel

  end subroutine change_species
  
  subroutine refuel
    double precision :: dist_to_C_sq
    double precision :: dist_to_N_sq
    double precision :: far
    double precision :: x(3)
    integer :: n

    far = (L(1)*0.45)**2

    !$omp parallel do private(x, dist_to_C_sq, dist_to_N_sq)
    do n = 1,solvent% Nmax
       if (solvent% species(n) == 2) then
          x = rel_pos(colloids% pos(:,1), solvent% pos(:,n), solvent_cells% edges)
          dist_to_C_sq = dot_product(x, x)
          x= rel_pos(colloids% pos(:,2), solvent% pos(:,n), solvent_cells% edges)
          dist_to_N_sq = dot_product(x, x)
          if ((dist_to_C_sq > far) .and. (dist_to_N_sq > far)) then
             solvent% species(n) = 1
          end if
       end if
    end do
  end subroutine refuel

  subroutine read_links(filename, n_links, links, links_d)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_links
    integer, allocatable, intent(out) :: links(:,:)
    double precision, allocatable, intent(out) :: links_d(:)

    integer :: i, iostat
    integer :: i1, i2
    double precision :: l

    n_links = 0
    open(11, file=filename)
    count_loop: do while (.true.)
       read(11, *, iostat=iostat) i1, i2, l
       if (iostat.lt.0) exit count_loop
       n_links = n_links + 1
    end do count_loop
    close(11)

    allocate(links(2, n_links), links_d(n_links))
    open(11, file=filename)
    do i = 1, n_links
       read(11, *, iostat=iostat) i1, i2, l
       links(1, i) = i1+1
       links(2, i) = i2+1
       links_d(i) = l
    end do
    close(11)

  end subroutine read_links

  function get_unit_r() result(u_r)
    double precision :: u_r(3)
    double precision :: r1(3), r2(3), r3(3)
    double precision, parameter :: one_z(3) = [0.d0, 0.d0, 1.d0]

    if (do_quaternion) then
       u_r = qrot(rigid_janus%q, one_z)
    else if (attr_exists) then
       r1 = rel_pos(colloids%pos(:,1), com_pos, solvent_cells%edges)
       r2 = rel_pos(colloids%pos(:,2), com_pos, solvent_cells%edges)
       r3 = rel_pos(colloids%pos(:,3), com_pos, solvent_cells%edges)
       u_r = (r1 + alpha*r2 + beta*r3) / z0
    else
       u_r = one_z
    end if

  end function get_unit_r

  subroutine compute_cell_wise_max_v
    integer :: i, j
    double precision :: local_max_v

    local_max_v = 0
    !$omp parallel do private(i, j) reduction(max:local_max_v)
    do i = 1, solvent_cells%N
       local_max_v = 0
       do j = solvent_cells%cell_start(i), solvent_cells%cell_start(i) + solvent_cells%cell_count(i) - 1
          local_max_v = max(local_max_v, norm2(solvent%vel(:,j)))
       end do
    end do
    solvent_cells%max_v = local_max_v

  end subroutine compute_cell_wise_max_v

end program single_body
