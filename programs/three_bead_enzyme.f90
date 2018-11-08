! This file is part of RMPCDMD
! Copyright (c) 2018 Pierre de Buyl
! License: BSD 3-clause (see file LICENSE)

!> Simulate a three bead enzyme model
!!
!! \param L                          length of simulation box in the 3 dimensions
!! \param rho                        fluid number density
!! \param T                          Temperature. Used for setting initial velocities and (if enabled) bulk thermostatting.
!! \param tau                        MPCD collision time
!! \param N_enzymes                  Number of enzyme molecules
!! \param substrate_fraction         initial fraction of non-inert (substrate and product) fluid particles
!! \param product_relative_fraction  inital relative fraction of product among non-inert fluid particles
!! \param N_MD                       number MD steps occuring in tau
!! \param N_loop                     number of MPCD timesteps
!! \param colloid_sampling           interval (in MD steps) of sampling the colloid position and velocity
!! \param equilibration_loops        number of MPCD steps for equilibration
!! \param sigma_E                    radius of enzymatic_site
!! \param sigma_N                    radius of N sphere
!! \param link_d                     length of rigid link
!! \param link_angles                angle of the model at rest
!! \param elastic_k                  stiffness of the link
!! \param angle_k                    stiffness of the angular link
!! \param epsilon_E                  interaction parameter of E sphere with solvent species (3 elements)
!! \param epsilon_N                  interaction parameter of N sphere with solvent species (3 elements)
!! \param epsilon_colloid            interaction parameter among colloids
!! \param enzyme_capture_radius      capture and release radius for substrate/product
!! \param rate_release_s             rate of release of substrate from enzyme
!! \param rate_release_p             rate of release of product from enzyme
!! \param proba_s                    probability of capture of substrate
!! \param proba_p                    probability of capture of product
!! \param bulk_rmpcd                 use bulkd rmpcd reaction for B->A instead of resetting
!! \param bulk_rate                  rates for the A->B and B->A bulk reaction

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
  integer, parameter :: N_species_colloids = 2
  integer :: N_colloids, N_enzymes

  integer :: rho
  integer :: N
  integer :: error

  double precision :: sigma_N, sigma_E, max_cut
  double precision :: epsilon(N_species,N_species_colloids)
  double precision :: sigma(N_species,N_species_colloids), sigma_cut(N_species,N_species_colloids)
  double precision, allocatable :: mass(:)

  double precision :: elastic_k, link_d(2)
  double precision :: angle_k, link_angles(2)

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
  integer(HID_T) :: dummy_id
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  double precision :: v_com(3)
  double precision :: tmp_vec(3), normal_vec(3), tmp_q(4)
  type(particle_system_io_t) :: dimer_io
  double precision :: bulk_rate(2)
  logical :: bulk_rmpcd

  integer :: enzyme_i
  integer :: enz_1, enz_2, enz_3
  double precision :: dist
  double precision :: enzyme_capture_radius
  double precision :: proba_p, proba_s
  double precision :: rate_release_p, rate_release_s, total_rate
  double precision :: substrate_fraction, product_relative_fraction
  logical :: placement_too_close

  !! State variables for the enzyme kinetics
  integer, allocatable :: bound_molecule_id(:)
  double precision, allocatable :: next_reaction_time(:)
  logical, allocatable :: enzyme_bound(:)
  double precision, allocatable :: link_angle(:)

  integer, dimension(N_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el
  type(h5md_element_t) :: link_angle_el

  !! Histogram variable
  type(histogram_t), allocatable :: radial_hist(:)
  double precision, allocatable :: dummy_hist(:,:,:)

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

  call timer_list%init(17)
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
  elastic_k = PTread_d(config, 'elastic_k', loc=params_group)
  angle_k = PTread_d(config, 'angle_k', loc=params_group)

  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  if (colloid_sampling <= 0) then
     error stop 'colloid_sampling must be positive'
  end if
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)
  equilibration_loops = PTread_i(config, 'equilibration_loops', loc=params_group)

  N_enzymes = PTread_i(config, 'N_enzymes', loc=params_group)
  N_colloids = 3*N_enzymes

  allocate(link_angle(N_enzymes))
  allocate(next_reaction_time(N_enzymes))
  allocate(bound_molecule_id(N_enzymes))
  allocate(enzyme_bound(N_enzymes))

  link_angle = link_angles(1)

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
  sigma = sigma*1.1d0
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init(epsilon(1:N_species_colloids,1:N_species_colloids), &
       sigma(1:N_species_colloids,1:N_species_colloids), &
       sigma_cut(1:N_species_colloids,1:N_species_colloids))

  allocate(mass(3*N_enzymes))
  do i = 1, N_enzymes
     mass(3*(i-1)+1) = rho * sigma_N**3 * 4 * pi/3
     mass(3*(i-1)+2) = rho * sigma_E**3 * 4 * pi/3
     mass(3*(i-1)+3) = rho * sigma_N**3 * 4 * pi/3
  end do


  allocate(radial_hist(N_colloids))
  do i = 1, N_colloids
     call radial_hist(i)%init(min(sigma_E, sigma_N)/4, max(sigma_E, sigma_N)*3, 100, n_species=N_species)
  end do

  call solvent% init(N, N_species, system_name='solvent')

  call colloids% init(N_colloids, N_species_colloids, mass, system_name='colloids')
  deallocate(mass)

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call h5gclose_f(params_group, error)
  call PTkill(config)

  do enzyme_i = 1, N_enzymes
     colloids% species(3*(enzyme_i-1)+1) = 2
     colloids% species(3*(enzyme_i-1)+2) = 1
     colloids% species(3*(enzyme_i-1)+3) = 2
  end do

  colloids% vel = 0
  colloids% force = 0

  dimer_io%force_info%store = .false.
  dimer_io%id_info%store = .false.
  dimer_io%position_info%store = .true.
  dimer_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%position_info%step = colloid_sampling
  dimer_io%position_info%time = colloid_sampling*N_MD_steps*dt
  dimer_io%image_info%store = .true.
  dimer_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%image_info%step = dimer_io%position_info%step
  dimer_io%image_info%time = dimer_io%position_info%time
  dimer_io%velocity_info%store = .true.
  dimer_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%velocity_info%step = dimer_io%position_info%step
  dimer_io%velocity_info%time = dimer_io%position_info%time
  dimer_io%species_info%store = .true.
  dimer_io%species_info%mode = H5MD_FIXED
  call dimer_io%init(hfile, 'dimer', colloids)

  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
     if (dble(k)/dble(solvent%Nmax) <= substrate_fraction) then
        if (dble(k)/dble(solvent%Nmax)/substrate_fraction <= product_relative_fraction) then
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

  ! place enzymes without overlap
  do enzyme_i = 1, N_enzymes
     enz_1 = 3*(enzyme_i-1)+1
     enz_2 = enz_1+1
     enz_3 = enz_1+2

     placement_too_close = .true.
     do while (placement_too_close)

        colloids% pos(:,enz_2) = [ threefry_double(state(1))*solvent_cells%edges(1), &
             threefry_double(state(1))*solvent_cells%edges(2), &
             threefry_double(state(1))*solvent_cells%edges(3) ]

        ! Pick one vector to place bead 1 and another one normal to the first
        tmp_vec = rand_sphere(state(1))
        normal_vec = rand_sphere(state(1))
        normal_vec = normal_vec - tmp_vec * dot_product(normal_vec, tmp_vec)
        ! Place normal_vec at an angle link_angles(1) w.r.t. tmp_vec
        tmp_q = qnew(s=cos(link_angles(1)/2), v=sin(link_angles(1)/2)*normal_vec/norm2(normal_vec))
        normal_vec = qvector(qmul(tmp_q, qmul(qnew(v=tmp_vec), qconj(tmp_q))))

        colloids%pos(:,enz_1) = colloids%pos(:,enz_2) + tmp_vec*link_d(1)
        colloids%pos(:,enz_3) = colloids%pos(:,enz_2) + normal_vec*link_d(2)

        placement_too_close = .false.

        do j = 1, 3*(enzyme_i-1)
           tmp_vec = colloids%pos(:,j)
           do i = 1, 3
              dist = norm2(rel_pos(colloids%pos(:,(enzyme_i-1)*3+i), &
                   tmp_vec, solvent_cells%edges))
              if ( dist < &
                   colloid_lj%sigma(colloids%species((enzyme_i-1)*3+i), colloids%species(j)) &
                   ) then
                 placement_too_close = .true.
              end if
           end do
        end do
     end do
  end do

  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=N_MD_steps, &
       time=N_MD_steps*dt)

  call link_angle_el%create_time(hfile%observables, 'link_angle', &
       link_angle, ior(H5MD_LINEAR, H5MD_STORE_TIME), step=N_MD_steps, &
       time=N_MD_steps*dt)

  call h5gcreate_f(dimer_io%group, 'box', box_group, error)
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

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  e3 = compute_bead_force()

  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  call h5fflush_f(hfile%id, H5F_SCOPE_GLOBAL_F, error)
  write(*,*) 'Running for', equilibration_loops, '+', N_loop, 'loops'

  call main%tic()
  sampling = .false.
  enzyme_bound = .false.
  next_reaction_time = huge(next_reaction_time(1))

  do i = 0, N_loop+equilibration_loops
     if (i==equilibration_loops) sampling = .true.
     if (modulo(i,64) == 0) write(*,'(i05)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        current_time = (i*N_MD_steps + (j-1))*dt
        call md_pos(solvent, dt)
        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids%mass(k))
        end do
        
        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1d0) .or. (so_max >= skin*0.89d0) ) then
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
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
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids%mass(k))
        end do
        call varia%tac()

     end do md_loop

     call solvent_cells%random_shift(state(1))
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

     !! TODO: add thermostat and hydro option
     call simple_mpcd_step(solvent, solvent_cells, state)

     if (bulk_rmpcd) then
        call bulk_reaction(solvent, solvent_cells, 1, 2, bulk_rate(1), tau, state)
        call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate(2), tau, state)
     end if

     call neigh% make_stencil(solvent_cells, enzyme_capture_radius)
     call neigh% update_list(colloids, solvent, enzyme_capture_radius, solvent_cells)
     if (i>equilibration_loops) then
        call reset_enzyme_region_bit
        call select_substrate
     end if

     ! Check if unbinding should occur
     do enzyme_i = 1, N_enzymes
        if (i*N_MD_steps*dt > next_reaction_time(enzyme_i)) then
           if (threefry_double(state(1)) < rate_release_s / total_rate) then
              ! release substrate
              call unbind_molecule(enzyme_i, 1)
           else
              ! release product
              call unbind_molecule(enzyme_i, 2)
           end if
           next_reaction_time(enzyme_i) = huge(next_reaction_time(enzyme_i))
        end if
     end do

     ! Rebuild the neighbor list after possible unbindings
     call neigh% make_stencil(solvent_cells, max_cut+skin)
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
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

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
           v_com = v_com + colloids%mass(k) * colloids%vel(:,k)
           kin_e = kin_e + colloids%mass(k) * sum(colloids%vel(:,k)**2) / 2
        end do

        call thermo_data%append(hfile, temperature, e1+e2+e3, kin_e, e1+e2+e3+kin_e, v_com)

        ! update radial histogram
        do k = 1, colloids%Nmax
           call compute_radial_histogram(radial_hist(k), colloids%pos(:,k), solvent_cells%edges, solvent)
        end do
        n_solvent = 0
        !$omp parallel do private(k,j) reduction(+:n_solvent)
        do k = 1, solvent%Nmax
           j = solvent%species(k)
           if (j <= 0) cycle
           n_solvent(j) = n_solvent(j) + 1
        end do
        call n_solvent_el%append(n_solvent)
        call link_angle_el%append(link_angle)
        if (modulo(i, colloid_sampling)==0) then
           call dimer_io%position%append(colloids%pos)
           call dimer_io%velocity%append(colloids%vel)
           call dimer_io%image%append(colloids%image)
        end if

     end if

     call varia%tac()

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  call thermo_data%append(hfile, temperature, e1+e2+e3, kin_e, e1+e2+e3+kin_e, &
       v_com, add=.false., force=.true.)

  ! write solvent coordinates for last step

  call h5gcreate_f(hfile%id, 'fields', fields_group, error)

  ! copy radial histogram data
  allocate(dummy_hist(size(radial_hist(1)%data, 1), &
       size(radial_hist(1)%data, 2), size(radial_hist)))
  do i = 1, size(radial_hist)
     call correct_radial_histogram(radial_hist(i))
     dummy_hist(:,:,i) = radial_hist(i)%data / dble(N_loop)
  end do

  call dummy_element%create_fixed(fields_group, 'radial_histogram', dummy_hist)

  call h5oopen_f(fields_group, 'radial_histogram', dummy_id, error)
  call h5md_write_attribute(dummy_id, 'xmin', radial_hist(1)%xmin)
  call h5md_write_attribute(dummy_id, 'dx', radial_hist(1)%dx)
  call h5oclose_f(dummy_id, error)

  call h5gclose_f(fields_group, error)


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
  call timer_list%append(colloids%time_apply_pbc)
  call timer_list%append(neigh%time_update)
  call timer_list%append(neigh%time_force)
  call timer_list%append(colloids%time_max_disp)

  call h5gcreate_f(hfile%id, 'timers', timers_group, error)
  call timer_list%write(timers_group, total_time)
  call h5md_write_dataset(timers_group, 'total', total_time)
  call h5md_write_dataset(timers_group, main%name, main%total)
  call h5gclose_f(timers_group, error)

  call dimer_io%close()
  call n_solvent_el%close()
  call link_angle_el%close()
  call hfile%close()
  call h5close_f(error)

contains

  function compute_bead_force() result(en)
    double precision :: en

    integer :: i, i_enzyme
    double precision :: f(3), r12(3), r
    double precision :: r32(3), d12, d32, costheta

    en = 0

    do i_enzyme = 1, N_enzymes
       do i = 1, 2
          r12 = rel_pos(colloids%pos(:,(i_enzyme-1)*3+i), colloids%pos(:,(i_enzyme-1)*3+i+1), solvent_cells%edges)
          r = norm2(r12)
          en = en + elastic_k*(r-link_d(i))**2/2
          f = -elastic_k*(r-link_d(i))*r12/r
          colloids%force(:,(i_enzyme-1)*3+i) = colloids%force(:,(i_enzyme-1)*3+i) + f
          colloids%force(:,(i_enzyme-1)*3+i+1) = colloids%force(:,(i_enzyme-1)*3+i+1) - f
       end do

       r12 = rel_pos(colloids%pos(:,(i_enzyme-1)*3+1), colloids%pos(:,(i_enzyme-1)*3+2), solvent_cells%edges)
       d12 = norm2(r12)
       r32 = rel_pos(colloids%pos(:,(i_enzyme-1)*3+3), colloids%pos(:,(i_enzyme-1)*3+2), solvent_cells%edges)
       d32 = norm2(r32)
       costheta = dot_product(r12, r32)/(d12*d32)
       f = -fprime(costheta, link_angle(i_enzyme)) * (r32/(d32*d12) - costheta * r12/d12**2)
       colloids%force(:,(i_enzyme-1)*3+1) = colloids%force(:,(i_enzyme-1)*3+1) + f
       colloids%force(:,(i_enzyme-1)*3+2) = colloids%force(:,(i_enzyme-1)*3+2) - f
       f = -fprime(costheta, link_angle(i_enzyme)) * (r12/(d12*d32) - costheta * r32/d32**2)
       colloids%force(:,(i_enzyme-1)*3+3) = colloids%force(:,(i_enzyme-1)*3+3) + f
       colloids%force(:,(i_enzyme-1)*3+2) = colloids%force(:,(i_enzyme-1)*3+2) - f

       en = en + angle_k * (acos(costheta)-link_angle(i_enzyme))**2 / 2
    end do

  end function compute_bead_force

  function compute_bead_energy(enzyme_idx) result(en)
    integer, intent(in) :: enzyme_idx
    double precision :: en

    double precision :: r12(3), d12, r32(3), d32, costheta

    r12 = rel_pos(colloids%pos(:,(enzyme_idx-1)*3+1), colloids%pos(:,(enzyme_idx-1)*3+2), solvent_cells%edges)
    d12 = norm2(r12)
    r32 = rel_pos(colloids%pos(:,(enzyme_idx-1)*3+3), colloids%pos(:,(enzyme_idx-1)*3+2), solvent_cells%edges)
    d32 = norm2(r32)
    costheta = dot_product(r12, r32)/(d12*d32)
    en = angle_k * (acos(costheta)-link_angle(enzyme_idx))**2 / 2

  end function compute_bead_energy

  !> derivative of the angular harmonic arccos term
  function fprime(x, theta_0) result(r)
    double precision, intent(in) :: x
    double precision, intent(in) :: theta_0
    double precision :: r

    r = - angle_k * (acos(x) - theta_0) / sqrt(1-x**2)

  end function fprime


  !> Select substrate molecules for binding to enzyme
  !!
  !! All substrate molecules are tested for their distance to the enzyme active site and are
  !! allowed to react once per stay in the enzyme region.
  subroutine select_substrate

    integer :: i, m, s_sp
    double precision :: dist, x_enzyme(3)
    double precision :: total_p, xi

    integer, parameter :: list_size = 256
    integer :: n_s, n_p
    integer :: idx
    integer :: list_s(list_size), list_p(list_size)
    logical :: reaction_happens

    integer :: enzyme_i, enz_2


    select_substrate_loop: do enzyme_i = 1, N_enzymes

       if (enzyme_bound(enzyme_i)) cycle select_substrate_loop

       enz_2 = (enzyme_i-1)*3 + 2

       n_s = 0
       n_p = 0

       x_enzyme = modulo(colloids%pos(:,enz_2), solvent_cells%edges)

       select_loop: do i = 1, neigh%n(enz_2)
          m = neigh%list(i, enz_2)

          ! skip substrates undergoing MD
          if (btest(solvent%flags(m), MD_BIT)) cycle select_loop

          s_sp = solvent%species(m)

          ! skip neutral fluid particles
          if (s_sp == 3) cycle select_loop

          dist = norm2(rel_pos(x_enzyme, solvent%pos(:,m), solvent_cells%edges))

          if (.not. btest(solvent%flags(m), ENZYME_REGION_BIT)) then
             if ( (dist <= enzyme_capture_radius) .and. (dist > solvent_colloid_lj%cut(s_sp, colloids%species(enz_2)) ) ) then
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
          call bind_molecule(enzyme_i, m)
       end if

    end do select_substrate_loop

  end subroutine select_substrate


  !> Bind solvent molecule idx to enzyme
  subroutine bind_molecule(enzyme_idx, idx)
    integer, intent(in) :: enzyme_idx
    integer, intent(in) :: idx
    double precision :: com_v(3), w_ab(3), excess_kinetic_energy, mu_ab
    double precision :: en1, en2
    integer :: p(3), cell_idx
    integer :: enz_2

    ! adjust velocity
    ! compute excess kinetic energy

    enz_2 = 3*(enzyme_idx-1)+2

    com_v = (colloids%vel(:,enz_2)*colloids%mass(enz_2) + solvent%vel(:,idx)) &
         / (colloids%mass(enz_2) + 1)
    w_ab = colloids%vel(:,enz_2) - solvent%vel(:,idx)

    mu_ab = 1/(1+1/colloids%mass(enz_2))

    excess_kinetic_energy = mu_ab * dot_product(w_ab, w_ab) / 2

    colloids%vel(:,enz_2) = com_v

    p = solvent_cells%cartesian_indices(solvent%pos(:, idx))

    cell_idx = compact_p_to_h(p, solvent_cells%M) + 1

    ! transfer mass
    colloids%mass(enz_2) = colloids%mass(enz_2) + 1

    bound_molecule_id(enzyme_idx) = solvent%id(idx)
    solvent%species(idx) = 0

    ! compute conformational energy difference when changing the angle
    en1 = compute_bead_energy(enzyme_idx)
    link_angle(enzyme_idx) = link_angles(2)
    en2 = compute_bead_energy(enzyme_idx)

    call add_energy_to_cell(cell_idx, excess_kinetic_energy + en1 - en2)

    e3 = e3 + en2 - en1

    enzyme_bound(enzyme_idx) = .true.

    ! sample from Poisson process
    next_reaction_time(enzyme_idx) = current_time - log(threefry_double(state(1)))/total_rate

  end subroutine bind_molecule

  !> Bind solvent molecule idx to enzyme
  subroutine unbind_molecule(enzyme_idx, to_species)
    integer, intent(in) :: enzyme_idx
    integer, intent(in) :: to_species

    integer :: i, p(3), cell_idx
    double precision :: x_new(3), dist, en1, en2
    logical :: too_close
    integer :: enz_1, enz_2, enz_3

    enz_1 = 3*(enzyme_idx-1)+1
    enz_2 = enz_1+1
    enz_3 = enz_1+2

    ! transfer mass
    colloids%mass(enz_2) = colloids%mass(enz_2) - 1

    ! Place the molecule outside of the interaction range of all colloids
    too_close = .true.
    placement_loop: do while (too_close)
       x_new = colloids%pos(:,enz_2) + rand_sphere(state(1))*enzyme_capture_radius
       x_new = modulo(x_new, solvent_cells%edges)
       p = solvent_cells%cartesian_indices(x_new)
       if (solvent_cells%cell_count(compact_p_to_h(p, solvent_cells%M)+1) < 3) cycle placement_loop
       too_close = .false.
       do i = 1, 3*N_enzymes
          dist = norm2(rel_pos(colloids%pos(:,i), x_new, solvent_cells%edges))
          too_close = too_close .or. &
               (dist < solvent_colloid_lj%cut(to_species, colloids%species(i)))
       end do
    end do placement_loop

    ! Retrieve location in particle array
    i = solvent%id_to_idx(bound_molecule_id(enzyme_idx))

    solvent%species(i) = to_species
    solvent%pos(:,i) = x_new

    ! Use c.o.m. velocity so that no kinetic energy exchange must take place
    solvent%vel(:,i) = colloids%vel(:,enz_2)

    ! compute conformational energy difference when changing the angle
    en1 = compute_bead_energy(enzyme_idx)
    link_angle(enzyme_idx) = link_angles(1)
    en2 = compute_bead_energy(enzyme_idx)

    p = solvent_cells%cartesian_indices(solvent%pos(:, i))

    cell_idx = compact_p_to_h(p, solvent_cells%M) + 1
    call add_energy_to_cell(cell_idx, en1 - en2)

    e3 = e3 + en2 - en1

    enzyme_bound(enzyme_idx) = .false.

  end subroutine unbind_molecule

  subroutine add_energy_to_cell(cell_idx, energy)
    integer, intent(in) :: cell_idx
    double precision, intent(in) :: energy

    integer :: i, start, n, n_effective
    integer :: cell(3), actual_cell(3), actual_idx, cell_shift(3)
    integer :: counter
    double precision :: com_v(3), kin, factor
    double precision :: remaining_energy, xi(3), change

    integer, parameter :: max_counter = 20

    start = solvent_cells% cell_start(cell_idx)
    n = solvent_cells% cell_count(cell_idx)

    remaining_energy = energy
    cell = compact_h_to_p(cell_idx - 1, solvent_cells%M) + 1

    ! If energy is positive, go for one cell.
    ! todo: always go around several cells

    if (energy > 0) then
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

       factor = sqrt((kin + energy)/kin)
       do i = start, start + n - 1
          solvent%vel(:,i) = com_v + factor*(solvent%vel(:,i) - com_v)
       end do

    else

       do counter = 1, max_counter
          xi(1) = -1 + 3*threefry_double(state(1))
          xi(2) = -1 + 3*threefry_double(state(1))
          xi(3) = -1 + 3*threefry_double(state(1))

          cell_shift = floor(xi)
          actual_cell = modulo(cell + cell_shift , solvent_cells% L)
          actual_idx = compact_p_to_h(actual_cell-1, solvent_cells%M) + 1

          start = solvent_cells% cell_start(actual_idx)
          n = solvent_cells% cell_count(actual_idx)

          n_effective = 0
          com_v = 0
          do i = start, start + n - 1
             if (solvent%species(i) > 0) then
                com_v = com_v + solvent% vel(:, i)
                n_effective = n_effective + 1
             end if
          end do

          if (n_effective <= 1) cycle

          com_v = com_v / n_effective

          kin = 0
          do i = start, start + n - 1
             if (solvent%species(i) > 0) then
                kin = kin + sum((solvent%vel(:,i)-com_v)**2)
             end if
          end do
          kin = kin / 2

          change = max(-kin/2, remaining_energy)

          remaining_energy = remaining_energy - change

          factor = sqrt((kin + change)/kin)
          do i = start, start + n - 1
             solvent%vel(:,i) = com_v + factor*(solvent%vel(:,i) - com_v)
          end do

          if (abs(remaining_energy) < 1d-9) exit
       end do

       if (abs(remaining_energy) > 2d-1) stop 'error in add_energy_to_cell'

    end if

  end subroutine add_energy_to_cell

  subroutine reset_enzyme_region_bit

    integer :: cell_idx
    integer :: i, start, n

    !$omp parallel do private(cell_idx, i, start, n)
    do cell_idx = 1, solvent_cells%N

       start = solvent_cells%cell_start(cell_idx)
       n = solvent_cells%cell_count(cell_idx)

       if (.not. solvent_cells%is_md(cell_idx)) then
          do i = start, start + n - 1
             solvent%flags(i) = ibclr(solvent%flags(i), ENZYME_REGION_BIT)
          end do
       end if

    end do

  end subroutine reset_enzyme_region_bit

end program three_bead_enzyme
