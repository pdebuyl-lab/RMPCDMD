! This file is part of RMPCDMD
! Copyright (c) 2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Simulate a single Janus particle
!!
!! This simulations models a chemically active Janus particle in a periodic simulation box.
!!
!! The coordinates of the Janus particle's beads must be provided in a H5MD file, as a
!! "fixed-in-time" dataset. The body of the particle can operate as a rigid-body (RATTLE) or
!! an elastic network.
!!
!! \param L                     length of simulation box in the 3 dimensions
!! \param rho                   fluid number density
!! \param T                     Temperature. Used for setting initial velocities and (if enabled) bulk thermostatting.
!! \param tau                   MPCD collision time
!! \param probability           probability to change A to B upon collision
!! \param bulk_rate             rate of B->A reaction
!! \param N_MD                  number MD steps occuring in tau
!! \param N_loop                number of MPCD timesteps
!! \param colloid_sampling      interval (in MD steps) of sampling the colloid position and velocity
!! \param equilibration_loops   number of MPCD steps for equilibration
!! \param epsilon_C             interaction parameter of C sphere with both solvent species (2 elements)
!! \param epsilon_N             interaction parameter of N sphere with both solvent species (2 elements)
!! \param data_filename         filename for input Janus coordinates
!! \param epsilon_colloid       interaction parameter for colloid-colloid interactions
!! \param link_treshold         distance criterion for finding rigid-body links
!! \param do_read_links         read link information from a file
!! \param links_file            filename for the links data
!! \param do_rattle             perform RATTLE
!! \param do_lennard_jones      compute colloid-colloid Lennard-Jones forces
!! \param do_elastic            compute colloid-colloid elastic network forces
!! \param elastic_k             elastic constant for the elastic network
!! \param rattle_pos_tolerance  absolute tolerance for Rattle (position part)
!! \param rattle_vel_tolerance  absolute tolerance for Rattle (velocity part)
!! \param sigma                 radius of the colloidal beads for colloid-solvent interactions
!! \param sigma_colloid         radius of the colloidal beads for colloid-colloid interactions
!! \param polar_r_max           maximal radius for the polar fields

program single_janus_pbc
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
  type(particle_system_t) :: colloids_com
  type(neighbor_list_t) :: neigh
  type(neighbor_list_t) :: neigh_com
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer, parameter :: N_species = 2

  integer :: rho
  integer :: N
  integer :: error

  double precision :: max_cut
  double precision :: epsilon(2,2)
  double precision :: sigma, sigma_v(2,2), sigma_cut(2,2)
  double precision :: mass(2)

  double precision :: e1, e2, e3
  double precision :: tau, dt , T
  double precision :: prob
  double precision :: bulk_rate
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting

  type(PTo) :: config

  integer :: i, L(3)
  integer :: j, k

  double precision :: total_time
  type(timer_t), target :: varia, main, time_flag, time_refuel, time_change
  type(timer_t), target :: bulk_reac_timer
  type(timer_list_t) :: timer_list
  integer(HID_T) :: timers_group

  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group, params_group
  integer(HID_T) :: connectivity_group
  integer(HID_T) :: box_group
  integer(HID_T) :: tmp_id
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  double precision :: v_com(3), x(3)
  double precision :: tmp_mass
  type(particle_system_io_t) :: janus_io
  type(particle_system_io_t) :: solvent_io

  integer, dimension(N_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el

  type(histogram_t) :: radial_hist
  type(h5md_element_t) :: radial_hist_el
  double precision :: polar_r_max
  type(polar_fields_t) :: polar
  integer(HID_T) :: polar_id

  logical :: do_rattle, do_read_links, do_lennard_jones, do_elastic
  integer, allocatable :: links(:,:)
  double precision, allocatable :: links_d(:)
  double precision :: link_treshold
  integer :: i_link
  double precision :: dist, rattle_pos_tolerance, rattle_vel_tolerance
  double precision :: elastic_k
  double precision :: unit_r(3)
  double precision :: com_pos(3), head_pos(3)
  integer :: head

  integer, parameter :: block_length = 8
  type(axial_correlator_t) :: axial_cf
  integer(HID_T) :: correlator_group

  integer :: equilibration_loops
  integer :: colloid_sampling
  logical :: sampling

  type(args_t) :: args
  character(len=144) :: data_filename
  character(len=144) :: links_file

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call main%init('main')
  call varia%init('varia')
  call time_flag%init('flag')
  call time_refuel%init('refuel')
  call time_change%init('change')
  call bulk_reac_timer%init('bulk_reac')

  call timer_list%init(20)
  call timer_list%append(varia)
  call timer_list%append(time_flag)
  call timer_list%append(time_refuel)
  call timer_list%append(time_change)
  call timer_list%append(bulk_reac_timer)

  call h5open_f(error)
  call hfile%create(args%output_file, 'RMPCDMD::single_janus_pbc', &
       'N/A', 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)


  prob = PTread_d(config,'probability', loc=params_group)
  bulk_rate = PTread_d(config,'bulk_rate', loc=params_group)

  L = PTread_ivec(config, 'L', 3, loc=params_group)
  rho = PTread_i(config, 'rho', loc=params_group)
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T', loc=params_group)
  
  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
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

  do_elastic = PTread_l(config, 'do_elastic', loc=params_group)
  if (do_elastic) elastic_k = PTread_d(config, 'elastic_k', loc=params_group)

  sigma = PTread_d(config, 'sigma_colloid', loc=params_group)
  sigma_v = sigma
  sigma_cut = sigma_v*2**(1.d0/6.d0)
  mass = rho * sigma**3 * 4 * 3.14159265/3

  epsilon = PTread_d(config, 'epsilon_colloid', loc=params_group)

  call colloid_lj% init(epsilon, sigma_v, sigma_cut)

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
  link_treshold = PTread_d(config, 'link_treshold', loc=params_group)
  do_read_links = PTread_l(config, 'do_read_links', loc=params_group)
  if (do_read_links) links_file = PTread_s(config, 'links_file', loc=params_group)

  polar_r_max = PTread_d(config, 'polar_r_max', loc=params_group)

  call h5gclose_f(params_group, error)
  call PTkill(config)

  call axial_cf%init(block_length, N_loop, N_loop*N_MD_steps)

  call colloids%init_from_file(data_filename, 'janus', H5MD_FIXED)
  head = get_head_idx(data_filename)
  colloids%image = 0
  colloids%vel = 0
  colloids%force = 0
  colloids%mass = mass
  call colloids_com%init(1)

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

  write(*,*) 'number of links:', i_link

  janus_io%force_info%store = .false.
  janus_io%id_info%store = .false.
  janus_io%position_info%store = .true.
  janus_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%position_info%step = colloid_sampling
  janus_io%position_info%time = colloid_sampling*dt
  janus_io%image_info%store = .true.
  janus_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%image_info%step = colloid_sampling
  janus_io%image_info%time = colloid_sampling*dt
  janus_io%velocity_info%store = .true.
  janus_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%velocity_info%step = colloid_sampling
  janus_io%velocity_info%time = colloid_sampling*dt
  janus_io%species_info%store = .true.
  janus_io%species_info%mode = H5MD_FIXED
  call janus_io%init(hfile, 'janus', colloids)
  ! The index in the H5MD file is 0-based.
  call h5md_write_attribute(janus_io%position%id, 'head', head-1)

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

  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)

  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, H5MD_LINEAR, step=N_MD_steps, &
       time=N_MD_steps*dt)

  call h5gcreate_f(janus_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  call radial_hist%init(0.d0, 5*sigma, 100, 2)
  call radial_hist_el%create_time(fields_group, 'radial_histogram', radial_hist%data, &
       H5MD_LINEAR, step=N_MD_steps, time=N_MD_steps*dt)
  call h5md_write_attribute(radial_hist_el%id, 'xmin', radial_hist%xmin)
  call h5md_write_attribute(radial_hist_el%id, 'dx', radial_hist%dx)

  call h5gcreate_f(hfile%id, 'connectivity', connectivity_group, error)
  call h5md_write_dataset(connectivity_group, 'janus_links', links-1)
  call h5dopen_f(connectivity_group, 'janus_links', tmp_id, error)
  call h5md_write_attribute(tmp_id, 'particles_group', 'janus')
  call h5dclose_f(tmp_id, error)
  call h5gclose_f(connectivity_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  call polar%init(N_species, 64, sigma, polar_r_max, 64)
  call neigh% init(colloids% Nmax, int(500*max(sigma,1.d0)**3))
  call neigh_com%init(colloids_com%Nmax, int(rho*16*polar_r_max**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)
  call neigh_com%make_stencil(solvent_cells, polar_r_max)

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
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  write(*,*) 'Running for', equilibration_loops, '+', N_loop, 'loops'
  call main%tic()
  sampling = .false.
  do i = 0, N_loop+equilibration_loops
     if (i==equilibration_loops) sampling = .true.
     if (modulo(i,32) == 0) write(*,'(i05)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        call md_pos(solvent, dt)

        call varia%tic()
        ! Extra copy for rattle
        colloids% pos_rattle = colloids% pos
        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(colloids%species(k)))
        end do
        call varia%tac()

        if (do_rattle) call rattle_body_pos(colloids, links, links_d, dt, solvent_cells% edges, rattle_pos_tolerance)

        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin/2) .or. (so_max >= skin/2) ) then
           call varia%tic()
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call varia%tac()
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, max_cut + skin, solvent_cells, solvent_colloid_lj)
           call varia%tic()
           !$omp parallel do
           do k = 1, solvent%Nmax
              solvent% pos_old(:,k) = solvent% pos(:,k)
           end do
           colloids% pos_old = colloids% pos
           n_extra_sorting = n_extra_sorting + 1
           call varia%tac()
        end if

        call varia%tic()
        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)

        !$omp parallel do
        do k = 1, solvent%Nmax
           solvent% force(:,k) = 0
        end do
        colloids% force = 0
        call varia%tac()
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        if (do_lennard_jones) e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
        if (do_elastic) e3 = elastic_network(colloids, links, links_d, elastic_k, solvent_cells%edges)

        call md_vel(solvent, dt)

        call varia%tic()
        do k = 1, colloids%Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids%mass(colloids%species(k)))
        end do
        call varia%tac()

        if (do_rattle) call rattle_body_vel(colloids, links, links_d, dt, solvent_cells% edges, rattle_pos_tolerance)

        if (sampling) then
           v_com = sum(colloids%vel, dim=2)/colloids%Nmax
           com_pos = sum((colloids%pos(:,:)+ &
                colloids%image(:,:)*spread(solvent_cells%edges, dim=2, ncopies=colloids%Nmax) ) &
                , dim=2)/ colloids%Nmax
           head_pos = colloids%pos(:,head)+colloids%image(:,head)*solvent_cells%edges
           unit_r = rel_pos(head_pos, com_pos, solvent_cells%edges)
           unit_r = unit_r / norm2(unit_r)
           call axial_cf%add_fast((i-equilibration_loops)*N_MD_steps+j-1, v_com, unit_r)
        end if

        if ((sampling) .and. (modulo(j, colloid_sampling)==0)) then
           call janus_io%position%append(colloids%pos)
           call janus_io%velocity%append(colloids%vel)
           call janus_io%image%append(colloids%image)
        end if

        call time_flag%tic()
        call flag_particles
        call time_flag%tac()
        call time_change%tic()
        call change_species
        call time_change%tac()

     end do md_loop

     call varia%tic()

     if (sampling) then
        temperature = compute_temperature(solvent, solvent_cells)
        kin_e = sum(solvent% vel**2)
        v_com = sum(solvent%vel, dim=2)
        tmp_mass = solvent%Nmax
        do k = 1, colloids%Nmax
           v_com = v_com + colloids%mass(colloids%species(k)) * colloids%vel(:,k)
           tmp_mass = tmp_mass + colloids%mass(colloids%species(k))
           kin_e = kin_e + colloids%mass(colloids%species(k)) * sum(colloids%vel(:,k)**2)
        end do
        v_com = v_com / tmp_mass
        kin_e = kin_e / 2

        call thermo_data%append(hfile, temperature, e1+e2+e3, kin_e, e1+e2+e3+kin_e, v_com)
        call axial_cf%add(i-equilibration_loops, com_pos, unit_r)

     end if

     call solvent_cells%random_shift(state(1))

     call varia%tac()

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

     call simple_mpcd_step(solvent, solvent_cells, state)

     call bulk_reac_timer%tic()
     call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate, tau, state)
     call bulk_reac_timer%tac()

     if (sampling) then
        n_solvent = 0
        do k = 1, solvent%Nmax
           j = solvent%species(k)
           if (j <= 0) continue
           n_solvent(j) = n_solvent(j) + 1
        end do
        call n_solvent_el%append(n_solvent)

        radial_hist%data = 0
        call compute_radial_histogram(radial_hist, colloids%pos(:,1), &
             solvent_cells%edges, solvent)
        call radial_hist_el%append(radial_hist%data)

        call varia%tic()
        colloids_com%pos(:,1) = modulo(com_pos, solvent_cells%edges)
        call neigh_com%update_list(colloids_com, solvent, polar_r_max, solvent_cells)
        v_com = sum(colloids%vel, dim=2)/colloids%Nmax
        call polar%update(colloids_com%pos(:,1), v_com, unit_r, solvent, neigh_com, solvent_cells)
        call varia%tac()
     end if

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  ! create a group for block correlators and write the data

  call h5gcreate_f(hfile%id, 'block_correlators', correlator_group, error)
  call axial_cf%write(correlator_group, N_MD_steps, N_MD_steps*dt, 1, dt)
  call h5gclose_f(correlator_group, error)

  ! write solvent coordinates for last step

  call thermo_data%append(hfile, temperature, e1+e2+e3, kin_e, e1+e2+e3+kin_e, v_com, add=.false., force=.true.)
  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

  ! store polar fields

  polar%c = polar%c / N_loop
  where (polar%count>0)
     polar%v(1,:,:,:) = polar%v(1,:,:,:) / polar%count
     polar%v(2,:,:,:) = polar%v(2,:,:,:) / polar%count
  end where
  call h5md_write_dataset(fields_group, 'polar_concentration', polar%c)
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

  ! Store timing data
  call timer_list%append(solvent%time_stream)
  call timer_list%append(solvent%time_step)
  call timer_list%append(solvent%time_count)
  call timer_list%append(solvent%time_sort)
  call timer_list%append(solvent%time_ct)
  call timer_list%append(solvent%time_md_pos)
  call timer_list%append(solvent%time_md_vel)
  call timer_list%append(solvent%time_max_disp)
  call timer_list%append(colloids%time_self_force)
  call timer_list%append(colloids%time_rattle_pos)
  call timer_list%append(colloids%time_rattle_vel)
  call timer_list%append(colloids%time_elastic)
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
  call solvent_io%close()
  call radial_hist_el%close()
  call n_solvent_el%close()
  call hfile%close()
  call h5close_f(error)

contains

  subroutine flag_particles
    double precision :: dist_to_C_sq
    integer :: i, r, s, thread_id
    double precision :: x(3)

    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(i, s, r, x, dist_to_C_sq)
    do i = 1, colloids%Nmax
       if (colloids%species(i)==1) then
          do s = 1,neigh% n(i)
             r = neigh%list(s,i)
             if (solvent% species(r) == 1) then
                x = rel_pos(colloids% pos(:,i),solvent% pos(:,r),solvent_cells% edges)
                dist_to_C_sq = dot_product(x, x)
                if (dist_to_C_sq < solvent_colloid_lj%cut_sq(1,1)) then
                   if (threefry_double(state(thread_id)) <= prob) then
                      !$omp atomic write
                      solvent% flag(r) = 1
                   end if
                end if
             end if
          end do
       end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine flag_particles

  subroutine change_species
    double precision :: dist_to_colloid_sq
    integer :: m, m_colloid
    double precision :: x(3)
    logical :: do_change
    integer :: s_m, s_colloid

    !$omp parallel do private(m, s_m, x, s_colloid, dist_to_colloid_sq, m_colloid, do_change)
    do m = 1, solvent% Nmax
       s_m = solvent%species(m)
       if (solvent% flag(m) == 1) then
          do_change = .true.
          colloid_dist_loop: do m_colloid = 1, colloids%Nmax
             s_colloid = colloids%species(m_colloid)
             x = rel_pos(colloids% pos(:,m_colloid), solvent% pos(:,m), solvent_cells% edges)
             dist_to_colloid_sq = dot_product(x, x)
             if (dist_to_colloid_sq < solvent_colloid_lj%cut_sq(s_m,s_colloid)) then
                do_change = .false.
                exit colloid_dist_loop
             end if
          end do colloid_dist_loop
          if (do_change) then
             solvent% species(m) = 2
             solvent% flag(m) = 0
          end if
       end if
    end do

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

  !> Get index of the head bead in the Janus particle
  function get_head_idx(filename) result(head)
    character(len=*), intent(in) :: filename
    integer :: head

    integer(HID_T) :: file_id, obj_id, attr_id, space_id
    integer :: error
    integer(HSIZE_T) :: dims(1)
    logical :: is_simple

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    call h5oopen_f(file_id, 'particles/janus/position', obj_id, error)
    call h5aopen_f(obj_id, 'head', attr_id, error)
    call h5aget_space_f(attr_id, space_id, error)
    call h5sis_simple_f(space_id, is_simple, error)
    if (.not. is_simple) error stop 'head attribute not simple'
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, head, dims, error)
    call h5sclose_f(space_id, error)
    call h5aclose_f(attr_id, error)
    call h5oclose_f(obj_id, error)
    call h5fclose_f(file_id, error)

    ! The index in the H5MD file is 0-based.
    head = head + 1

  end function get_head_idx

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

end program single_janus_pbc
