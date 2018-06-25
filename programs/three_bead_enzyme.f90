! This file is part of RMPCDMD
! Copyright (c) 2015-2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Simulate a three bead enzyme model
!!
!! \param L                       length of simulation box in the 3 dimensions
!! \param rho                     fluid number density
!! \param T                       Temperature. Used for setting initial velocities and (if enabled) bulk thermostatting.
!! \param tau                     MPCD collision time
!! \param bulk_rmpcd              use bulkd rmpcd reaction for B->A instead of resetting
!! \param bulk_rate               rates for the A->B and B->A bulk reaction
!! \param N_MD                    number MD steps occuring in tau
!! \param N_loop                  number of MPCD timesteps
!! \param colloid_sampling        interval (in MD steps) of sampling the colloid position and velocity
!! \param equilibration_loops     number of MPCD steps for equilibration
!! \param sigma_E                 radius of enzymatic_site
!! \param sigma_N                 radius of N sphere
!! \param link_d                  length of rigid link
!! \param link_angles             angle of the model at rest
!! \param elastic_k               stiffness of the link
!! \param angle_k                 stiffness of the angular link
!! \param epsilon_C               interaction parameter of C sphere with both solvent species (2 elements)
!! \param epsilon_N               interaction parameter of N sphere with both solvent species (2 elements)
!! \param epsilon_colloid         interaction parameter among C spheres

program three_bead_enzyme
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

  integer, parameter :: N_species = 3
  integer, parameter :: N_colloids = 3
  integer, parameter :: N_species_colloids = 2

  integer :: rho
  integer :: N
  integer :: error

  double precision :: sigma_N, sigma_E, max_cut
  double precision :: epsilon(N_species,N_species_colloids)
  double precision :: sigma(N_species,N_species_colloids), sigma_cut(N_species,N_species_colloids)
  double precision :: mass(2)

  double precision :: elastic_k, link_d(2)
  double precision :: angle_k, link_angle, link_angles(2)

  double precision :: e1, e2, e3
  double precision :: tau, dt , T
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting

  type(PTo) :: config

  integer :: i, L(3)
  double precision :: current_time
  integer :: j, k
  type(timer_t), target :: varia, main, time_flag, time_change
  double precision :: total_time
  type(timer_list_t) :: timer_list
  integer(HID_T) :: timers_group

  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group, params_group
  integer(HID_T) :: box_group
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  double precision :: v_com(3)
  double precision :: com_pos(3)
  double precision :: unit_r(3)
  type(particle_system_io_t) :: dimer_io
  type(particle_system_io_t) :: solvent_io
  double precision :: bulk_rate(2)
  logical :: bulk_rmpcd

  logical :: enzyme_bound
  integer :: bound_molecule_id
  double precision :: enzyme_capture_radius
  double precision :: proba_p, proba_s
  double precision :: rate_release_p, rate_release_s, total_rate, next_reaction_time
  double precision :: substrate_fraction, product_relative_fraction

  integer, dimension(N_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el

  integer, parameter :: block_length = 8
  type(axial_correlator_t) :: axial_cf
  integer(HID_T) :: correlator_group

  integer :: equilibration_loops
  integer :: colloid_sampling
  logical :: sampling
  type(args_t) :: args

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call main%init('main')
  call varia%init('varia')
  call time_flag%init('flag')
  call time_change%init('change')

  call timer_list%init(15)
  call timer_list%append(varia)
  call timer_list%append(time_flag)
  call timer_list%append(time_change)

  call h5open_f(error)
  call hfile%create(args%output_file, 'RMPCDMD::three_bead_enzyme', &
       RMPCDMD_REVISION, 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)

  bulk_rmpcd = PTread_l(config, 'bulk_rmpcd', loc=params_group)
  bulk_rate = PTread_d(config, 'bulk_rate', loc=params_group)

  proba_s = PTread_d(config, 'proba_s', loc=params_group)
  proba_p = PTread_d(config, 'proba_p', loc=params_group)
  rate_release_s = PTread_d(config, 'rate_release_s', loc=params_group)
  rate_release_p = PTread_d(config, 'rate_release_p', loc=params_group)
  total_rate = rate_release_s + rate_release_p
  enzyme_capture_radius = PTread_d(config, 'enzyme_capture_radius', loc=params_group)
  substrate_fraction = PTread_d(config, 'substrate_fraction', loc=params_group)
  product_relative_fraction = PTread_d(config, 'product_relative_fraction', loc=params_group)

  L = PTread_ivec(config, 'L', 3, loc=params_group)
  rho = PTread_i(config, 'rho', loc=params_group)
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T', loc=params_group)
  link_d = PTread_dvec(config, 'link_d', 2, loc=params_group)
  link_angles = PTread_dvec(config, 'link_angles', 2, loc=params_group)
  link_angle = link_angles(1)
  elastic_k = PTread_d(config, 'elastic_k', loc=params_group)
  angle_k = PTread_d(config, 'angle_k', loc=params_group)

  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  if (modulo(N_MD_steps, colloid_sampling) /= 0) then
     error stop 'colloid_sampling must divide N_MD with no remainder'
  end if
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)
  equilibration_loops = PTread_i(config, 'equilibration_loops', loc=params_group)

  sigma_E = PTread_d(config, 'sigma_E', loc=params_group)
  sigma_N = PTread_d(config, 'sigma_N', loc=params_group)

  ! solvent index first, colloid index second, in solvent_colloid_lj
  epsilon(:,1) = PTread_dvec(config, 'epsilon_E', 3, loc=params_group)
  epsilon(:,2) = PTread_dvec(config, 'epsilon_N', 3, loc=params_group)

  sigma(:,1) = sigma_E
  sigma(:,2) = sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)

  epsilon(:,:) = PTread_d(config, 'epsilon_colloid', loc=params_group)

  sigma(1,1) = 2*sigma_E
  sigma(1,2) = sigma_E + sigma_N
  sigma(2,1) = sigma_E + sigma_N
  sigma(2,2) = 2*sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init(epsilon(1:N_species_colloids,1:N_species_colloids), &
       sigma(1:N_species_colloids,1:N_species_colloids), &
       sigma_cut(1:N_species_colloids,1:N_species_colloids))

  mass(1) = rho * sigma_E**3 * 4 * 3.14159265/3
  mass(2) = rho * sigma_N**3 * 4 * 3.14159265/3

  call solvent% init(N,2, system_name='solvent') !there will be 2 species of solvent particles

  call colloids% init(N_colloids, 2, mass, system_name='colloids') !there will be 2 species of colloids

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call h5gclose_f(params_group, error)
  call PTkill(config)

  call axial_cf%init(block_length, N_loop, N_loop*N_MD_steps)

  colloids% species(1) = 2
  colloids% species(2) = 1
  colloids% species(3) = 2

  colloids% vel = 0
  colloids% force = 0

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

  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
     if (threefry_double(state(1)) < substrate_fraction) then
        ! substrate
        if (threefry_double(state(1)) < product_relative_fraction) then
           solvent%species(k) = 2
        else
           solvent%species(k) = 1
        end if
     else
        ! neutral species
        solvent%species(k) = 3
     end if
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0

  call solvent_cells%init(L, 1.d0)

  do i = 1, N_colloids
     colloids% pos(:,i) = solvent_cells%edges / 2
  end do
  colloids%pos(1,2) = colloids%pos(1,2) + link_d(1)
  colloids%pos(1,3) = colloids%pos(1,3) + link_d(2)*cos(link_angle)
  colloids%pos(2,3) = colloids%pos(2,3) + link_d(2)*sin(link_angle)
  
  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=N_MD_steps, &
       time=N_MD_steps*dt)

  call h5gcreate_f(dimer_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, int(rho*30*max(sigma_E,sigma_N)**3))

  skin = 2
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  call h5gclose_f(fields_group, error)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  e3 = compute_bead_force()

  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  call h5fflush_f(hfile%id, H5F_SCOPE_GLOBAL_F, error)
  write(*,*) 'Running for', equilibration_loops, '+', N_loop, 'loops'
  write(*,*) 'mass', mass 
  call main%tic()
  sampling = .false.
  enzyme_bound = .false.
  next_reaction_time = huge(next_reaction_time)

  do i = 0, N_loop+equilibration_loops
     if (i==equilibration_loops) sampling = .true.
     if (modulo(i,20) == 0) write(*,'(i05)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        current_time = (i*N_MD_steps + (j-1))*dt
        call md_pos(solvent, dt)
        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(colloids%species(k)))
        end do
        
        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1d0) .or. (so_max >= skin*0.89d0) ) then
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
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
        e3 = compute_bead_force()

        call md_vel(solvent, dt)

        call varia%tic()
        do k=1, colloids% Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(colloids%species(k)))
        end do
        call varia%tac()

        if (sampling) then
           v_com = sum(colloids%vel, dim=2)/2
           unit_r = rel_pos(colloids%pos(:,1), colloids%pos(:,2), solvent_cells%edges)
           unit_r = unit_r / norm2(unit_r)
           call axial_cf%add_fast((i-equilibration_loops)*N_MD_steps+j-1, v_com, unit_r)
        end if

        if ((sampling) .and. (modulo(j, colloid_sampling)==0)) then
           call dimer_io%position%append(colloids%pos)
           call dimer_io%velocity%append(colloids%vel)
           call dimer_io%image%append(colloids%image)
        end if

     end do md_loop

     if (i>equilibration_loops) call select_substrate

     if (i*N_MD_steps*dt > next_reaction_time) then
        if (threefry_double(state(1)) < rate_release_s / total_rate) then
           ! release substrate
           call unbind_molecule(1)
        else
           ! release product
           call unbind_molecule(2)
        end if
        next_reaction_time = huge(next_reaction_time)
     end if

     call varia%tic()

     if (sampling) then
        temperature = compute_temperature(solvent, solvent_cells)
        kin_e = 0
        v_com = 0
        !$omp parallel do private(k) reduction(+:kin_e) reduction(+:v_com)
        do k = 1, solvent%Nmax
           if (solvent%species(k)>0) then
              kin_e = kin_e + sum(solvent%vel(:,k)**2) / 2
              v_com = v_com + solvent%vel(:,k)
           end if
        end do
        do k = 1, colloids%Nmax
           v_com = v_com + colloids%mass(colloids%species(k)) * colloids%vel(:,k)
           kin_e = kin_e + colloids%mass(colloids%species(k)) * sum(colloids%vel(:,k)**2) / 2
        end do

        call thermo_data%append(hfile, temperature, e1+e2+e3, kin_e, e1+e2+e3+kin_e, v_com)

        com_pos = colloids%pos(:,2)
        call axial_cf%add(i-equilibration_loops, com_pos, unit_r)

     end if

     call solvent_cells%random_shift(state(1))
     call varia%tac()

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)
     call varia%tic()
     !$omp parallel do
     do k = 1, solvent%Nmax
        solvent% pos_old(:,k) = solvent% pos(:,k)
     end do
     colloids% pos_old = colloids% pos
     call varia%tac()

     call simple_mpcd_step(solvent, solvent_cells, state)
     
     if (bulk_rmpcd) then
        call bulk_reaction(solvent, solvent_cells, 1, 2, bulk_rate(1), tau, state)
        call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate(2), tau, state)
     end if

     if (sampling) then
        n_solvent = 0
        do k = 1, solvent%Nmax
           j = solvent%species(k)
           if (j <= 0) cycle
           n_solvent(j) = n_solvent(j) + 1
        end do
        call n_solvent_el%append(n_solvent)
     end if

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  call thermo_data%append(hfile, temperature, e1+e2+e3, kin_e, e1+e2+e3+kin_e, &
       v_com, add=.false., force=.true.)

  ! create a group for block correlators and write the data

  call h5gcreate_f(hfile%id, 'block_correlators', correlator_group, error)
  call axial_cf%write(correlator_group, N_MD_steps, N_MD_steps*dt, 1, dt)
  call h5gclose_f(correlator_group, error)

  ! write solvent coordinates for last step

  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

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
  call timer_list%append(neigh%time_update)
  call timer_list%append(neigh%time_force)
  call timer_list%append(colloids%time_max_disp)

  call h5gcreate_f(hfile%id, 'timers', timers_group, error)
  call timer_list%write(timers_group, total_time)
  call h5md_write_dataset(timers_group, 'total', total_time)
  call h5md_write_dataset(timers_group, main%name, main%total)
  call h5gclose_f(timers_group, error)

  call dimer_io%close()
  call solvent_io%close()
  call n_solvent_el%close()
  call hfile%close()
  call h5close_f(error)

contains

  function compute_bead_force() result(en)
    double precision :: en

    integer :: i
    double precision :: f(3), r12(3), r
    double precision :: r32(3), d12, d32, costheta

    en = 0
    do i = 1, 2
       r12 = rel_pos(colloids%pos(:,i), colloids%pos(:,i+1), solvent_cells%edges)
       r = norm2(r12)
       en = en + elastic_k*(r-link_d(i))**2/2
       f = -elastic_k*(r-link_d(i))*r12/r
       colloids%force(:,i) = colloids%force(:,i) + f
       colloids%force(:,i+1) = colloids%force(:,i+1) - f
    end do

    r12 = rel_pos(colloids%pos(:,1), colloids%pos(:,2), solvent_cells%edges)
    d12 = norm2(r12)
    r32 = rel_pos(colloids%pos(:,3), colloids%pos(:,2), solvent_cells%edges)
    d32 = norm2(r32)
    costheta = dot_product(r12, r32)/(d12*d32)
    f = -fprime(costheta) * (r32/(d32*d12) - costheta * r12/d12**2)
    colloids%force(:,1) = colloids%force(:,1) + f
    colloids%force(:,2) = colloids%force(:,2) - f
    f = -fprime(costheta) * (r12/(d12*d32) - costheta * r32/d32**2)
    colloids%force(:,3) = colloids%force(:,3) + f
    colloids%force(:,2) = colloids%force(:,2) - f

    en = en + angle_k * (acos(costheta)-link_angle)**2 / 2

  end function compute_bead_force

  function compute_bead_energy() result(en)
    double precision :: en

    double precision :: r12(3), d12, r32(3), d32, costheta

    r12 = rel_pos(colloids%pos(:,1), colloids%pos(:,2), solvent_cells%edges)
    d12 = norm2(r12)
    r32 = rel_pos(colloids%pos(:,3), colloids%pos(:,2), solvent_cells%edges)
    d32 = norm2(r32)
    costheta = dot_product(r12, r32)/(d12*d32)
    en = angle_k * (acos(costheta)-link_angle)**2 / 2
  end function compute_bead_energy

  !> derivative of the angular harmonic arccos term
  function fprime(x) result(r)
    double precision, intent(in) :: x
    double precision :: r

    r = - angle_k * (acos(x) - link_angle) / sqrt(1-x**2)

  end function fprime


  !> Select substrate molecules for binding to enzyme
  !!
  !! All substrate molecules are tested for their distance to the enzyme active site and are
  !! allowed to react once per stay in the enzyme region.
  subroutine select_substrate

    integer :: m, thread_id, i, s_sp
    double precision :: dist, x_enzyme(3)
    double precision :: total_p, xi

    integer, parameter :: list_size = 256
    integer :: n_s, n_p
    integer :: idx
    integer :: list_s(list_size), list_p(list_size)
    logical :: reaction_happens

    ! Do nothing if enzyme is already bound
    if (enzyme_bound) return

    n_s = 0
    n_p = 0

    x_enzyme = modulo(colloids%pos(:,2), solvent_cells%edges)

    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(i, m, s_sp, dist)
    do i = 1, solvent_cells%N
       select_loop: do m = solvent_cells%cell_start(i), solvent_cells%cell_start(i) + solvent_cells%cell_count(i) - 1
          s_sp = solvent%species(m)
          if (s_sp == 3) cycle select_loop
          dist = norm2(rel_pos(x_enzyme, solvent%pos(:,m), solvent_cells%edges))
          if (.not. btest(solvent%flags(m), ENZYME_REGION_BIT)) then
             if ( (dist <= enzyme_capture_radius) .and. (dist > solvent_colloid_lj%cut(s_sp, colloids%species(2)) ) ) then
                ! select for current round
                solvent%flags(m) = ibset(solvent%flags(m), ENZYME_REGION_BIT)
                if (s_sp==1) then
                   !$omp atomic
                   n_s = n_s + 1
                   if (n_s > list_size) then
                      stop 'exceed size of list_s in select_substrate'
                   end if
                   list_s(n_s) = m
                else if (s_sp==2) then
                   !$omp atomic
                   n_p = n_p + 1
                   if (n_p > list_size) then
                      stop 'exceed size of list_p in select_substrate'
                   end if
                   list_p(n_p) = m
                end if
             end if
          else
             if (dist > enzyme_capture_radius) then
                solvent%flags(m) = ibclr(solvent%flags(m), ENZYME_REGION_BIT)
             end if
          end if
       end do select_loop
    end do
    !$omp end do
    !$omp end parallel

    total_p = n_s*proba_s + n_p*proba_p

    reaction_happens = (threefry_double(state(1)) < total_p)

    ! Pick either a substrate or a product if a reaction happens

    if (reaction_happens) then
       xi = threefry_double(state(1))*total_p

       if (xi < n_s*proba_s) then
          ! pick in substrates
          idx = floor(xi * n_s / total_p) + 1
          m = list_s(idx)
          ! use list_s(idx)
       else
          ! pick in products
          idx = floor(xi * n_p / total_p) + 1
          m = list_p(idx)
       end if
       call bind_molecule(m)
    end if

  end subroutine select_substrate


  !> Bind solvent molecule idx to enzyme
  subroutine bind_molecule(idx)
    integer, intent(in) :: idx
    double precision :: com_v(3), w_ab(3), excess_kinetic_energy, mu_ab
    double precision :: en1, en2
    integer :: p(3), cell_idx

    ! adjust velocity
    ! compute excess kinetic energy

    com_v = (colloids%vel(:,2)*colloids%mass(colloids%species(2)) + solvent%vel(:,idx)) &
         / (colloids%mass(colloids%species(2)) + 1)
    w_ab = colloids%vel(:,2) - solvent%vel(:,idx)

    mu_ab = 1/(1+1/colloids%mass(colloids%species(2)))

    excess_kinetic_energy = mu_ab * dot_product(w_ab, w_ab) / 2

    colloids%vel(:,2) = com_v

    ! todo: deposit extra energy

    p = floor( (solvent%pos(:, idx) / solvent_cells%a) - solvent_cells%origin )
    cell_idx = compact_p_to_h(p, solvent_cells%M) + 1

    ! transfer mass
    colloids%mass(1) = colloids%mass(1) + 1

    bound_molecule_id = solvent%id(idx)
    solvent%species(idx) = 0

    en1 = compute_bead_energy()

    link_angle = link_angles(2)

    en2 = compute_bead_energy()

    call add_energy_to_cell(cell_idx, excess_kinetic_energy + en1 - en2)

    e3 = e3 + en2 - en1

    enzyme_bound = .true.

    next_reaction_time = current_time - log(threefry_double(state(1)))/total_rate ! sample from Poisson process

  end subroutine bind_molecule

  !> Bind solvent molecule idx to enzyme
  subroutine unbind_molecule(to_species)
    integer, intent(in) :: to_species

    integer :: i, p(3), cell_idx
    double precision :: x_new(3), dist, en1, en2
    logical :: too_close

    ! transfer mass
    colloids%mass(1) = colloids%mass(1) - 1

    ! Place the molecule outside of the interaction range of all colloids
    too_close = .true.
    do while (too_close)
       x_new = colloids%pos(:,2) + rand_sphere(state(1))*solvent_colloid_lj%cut(to_species, colloids%species(2))*1.001d0
       too_close = .false.
       do i = 1, 3
          dist = norm2(rel_pos(colloids%pos(:,i), x_new, solvent_cells%edges))
          too_close = too_close .or. &
               (dist < solvent_colloid_lj%cut(to_species, colloids%species(i)))
       end do
    end do

    i = solvent%id_to_idx(bound_molecule_id)

    solvent%species(i) = to_species
    solvent%pos(:,i) = x_new

    ! Use c.o.m. velocity so that no kinetic energy exchange must take place
    solvent%vel(:,i) = colloids%vel(:,2)

    en1 = compute_bead_energy()

    link_angle = link_angles(1)

    en2 = compute_bead_energy()

    p = floor( (solvent%pos(:, i) / solvent_cells%a) - solvent_cells%origin )
    cell_idx = compact_p_to_h(p, solvent_cells%M) + 1
    call add_energy_to_cell(cell_idx, en1 - en2)

    e3 = e3 + en2 - en1

    enzyme_bound = .false.

  end subroutine unbind_molecule

  subroutine add_energy_to_cell(cell_idx, energy)
    integer, intent(in) :: cell_idx
    double precision, intent(in) :: energy

    integer :: i, start, n, n_effective
    double precision :: com_v(3), kin, factor

    start = solvent_cells% cell_start(cell_idx)
    n = solvent_cells% cell_count(cell_idx)

    n_effective = 0
    com_v = 0
    do i = start, start + n - 1
       if (solvent%species(i) > 0) then
          com_v = com_v + solvent% vel(:, i)
          n_effective = n_effective + 1
       end if
    end do

    write(31,*) 'time', current_time, 'n_effective', n_effective, 'energy', energy
    if (n_effective == 0) then
       stop 'n_effective == 0 in add_energy_to_cell'
    end if
    if (n_effective == 1) then
       stop 'n_effective == 1 in add_energy_to_cell'
    end if
    com_v = com_v / n_effective

    kin = 0
    do i = start, start + n - 1
       if (solvent%species(i) > 0) then
          kin = kin + sum((solvent%vel(:,i)-com_v)**2)
       end if
    end do
    kin = kin / 2

    write(31,*) 'com_v', com_v, 'kin', kin

    if ((kin + energy) <= 0) then
       stop 'negative energy in add_energy_to_cell'
    end if
    factor = sqrt((kin + energy)/kin)
    do i = start, start + n - 1
       solvent%vel(:,i) = com_v + factor*(solvent%vel(:,i) - com_v)
    end do

  end subroutine add_energy_to_cell


end program three_bead_enzyme
