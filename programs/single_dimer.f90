! This file is part of RMPCDMD
! Copyright (c) 2015-2018 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Simulate a single dimer colloid
!!
!! Consider a dimer in a simulation box, either periodic or with wall in the z direction.
!!
!! The dimer bond length is maintained either by a RATTLE constraint or by an harmonic bond.
!!
!! \param L           length of simulation box in the 3 dimensions
!! \param rho         fluid number density
!! \param T           Temperature. Used for setting initial velocities and (if enabled) bulk thermostatting.
!! \param tau         MPCD collision time
!! \param N_MD        number MD steps occuring in tau
!! \param N_loop      number of MPCD timesteps
!! \param colloid_sampling interval (in MD steps) of sampling the colloid position and velocity
!! \param equilibration_loops number of MPCD steps for equilibration
!! \param sigma       LJ sigma (2 elements)
!! \param d           length of rigid link
!! \param do_rattle   apply RATTLE constraint to the dimer bond
!! \param do_harmonic compute harmonic for on the dimer bond
!! \param harmonic_k  stiffness of the dimer bond
!! \param epsilon     interaction parameter of collodis with the solvent (2 elements)
!! \param do_zwall              use a confining potential in the z direction, 9-3 Lennard-Jones
!! \param wall_sigma            wall LJ sigma
!! \param wall_epsilon          wall LJ epsilon
!! \param wall_shift            wall shift
!! \param fluid_wall            boundary condition for the fluid
!! \param initial_condition     initial condition for the dimer. CENTER (fixed orientation dimer at center of box) or PLANAR_RANDOM (dimer at center, random orientation in xy plane)

program single_dimer
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

  integer, parameter :: N_species = 1
  integer, parameter :: N_colloid_species = 2

  integer :: rho
  integer :: N
  integer :: error

  double precision :: max_cut
  double precision :: epsilon(N_species, N_colloid_species)
  double precision :: sigma(N_species, N_colloid_species), sigma_cut(N_species, N_colloid_species)
  double precision :: mass(N_colloid_species)
  logical :: do_zwall
  logical :: do_rattle, do_harmonic
  double precision :: harmonic_k
  double precision :: wall_sigma(3, N_species), wall_epsilon(3, N_species), wall_shift(3)

  double precision :: e1, e2, e3, e_wall
  double precision :: tau, dt , T
  double precision :: d
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting
  integer :: loop_i_last_sort

  type(PTo) :: config

  integer :: i, L(3)
  integer :: j, k
  type(timer_t), target :: varia, main
  double precision :: total_time
  double precision :: theta
  type(timer_list_t) :: timer_list
  integer(HID_T) :: timers_group

  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: params_group
  integer(HID_T) :: box_group
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  double precision :: v_com(3)
  double precision :: com_pos(3)
  double precision :: unit_r(3)
  type(particle_system_io_t) :: dimer_io

  integer, parameter :: block_length = 8
  type(axial_correlator_t) :: axial_cf
  integer(HID_T) :: correlator_group

  integer :: equilibration_loops
  integer :: colloid_sampling
  logical :: sampling
  type(args_t) :: args

  character(len=144) :: initial_condition, fluid_wall

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call main%init('main')
  call varia%init('varia')

  call timer_list%init(13)
  call timer_list%append(varia)

  call h5open_f(error)
  call hfile%create(args%output_file, 'RMPCDMD::single_dimer', &
       RMPCDMD_REVISION, 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)

  L = PTread_ivec(config, 'L', 3, loc=params_group)
  rho = PTread_i(config, 'rho', loc=params_group)
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T', loc=params_group)
  d = PTread_d(config, 'd', loc=params_group)
  do_rattle = PTread_l(config, 'do_rattle', loc=params_group)
  do_harmonic = PTread_l(config, 'do_harmonic', loc=params_group)
  harmonic_k = PTread_d(config, 'harmonic_k', loc=params_group)
  
  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  if (modulo(N_MD_steps, colloid_sampling) /= 0) then
     error stop 'colloid_sampling must divide N_MD with no remainder'
  end if
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)
  equilibration_loops = PTread_i(config, 'equilibration_loops', loc=params_group)

  sigma(1,:) = PTread_dvec(config, 'sigma', N_colloid_species, loc=params_group)

  ! solvent index first, colloid index second, in solvent_colloid_lj
  epsilon(1,:) = PTread_dvec(config, 'epsilon', N_colloid_species, loc=params_group)

  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)

  call colloid_lj%init(&
       reshape([2*sigma(1,1), sigma(1,2)+sigma(1,2), sigma(1,2)+sigma(1,2), 2*sigma(1,2)], [2, 2]),&
       reshape([2*epsilon(1,1), epsilon(1,1)+epsilon(1,2), epsilon(1,1)+epsilon(1,2), 2*epsilon(1,2)], [2, 2]),&
       reshape([2*sigma_cut(1,1), sigma_cut(1,1)+sigma_cut(1,2), sigma(1,1)+sigma_cut(1,2), 2*sigma_cut(1,2)], [2, 2]))

  ! no wall in x and y
  wall_sigma = -1
  wall_shift = 0

  do_zwall = PTread_l(config, 'do_zwall', loc=params_group)
  fluid_wall = 'PERIODIC'
  if (do_zwall) then
     wall_sigma(3,:) = PTread_d(config, 'wall_sigma', loc=params_group)
     wall_shift(3) = PTread_d(config, 'wall_shift', loc=params_group)
     wall_epsilon = PTread_d(config, 'wall_epsilon', loc=params_group)
     call walls_colloid_lj% init(wall_epsilon, &
          wall_sigma, 3.d0**(1.d0/6.d0)*wall_sigma, wall_shift)
     fluid_wall = PTread_s(config, 'fluid_wall', loc=params_group)
  end if

  initial_condition = PTread_s(config, 'initial_condition', loc=params_group)

  mass(1) = rho * sigma(1,1)**3 * 4 * 3.14159265/3
  mass(2) = rho * sigma(1,2)**3 * 4 * 3.14159265/3

  call solvent% init(N,2, system_name='solvent') !there will be 2 species of solvent particles

  call colloids% init(2,2, mass, system_name='colloids') !there will be 2 species of colloids

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call h5gclose_f(params_group, error)
  call PTkill(config)

  call axial_cf%init(block_length, N_loop, N_loop*N_MD_steps)
  
  colloids% species(1) = 1
  colloids% species(2) = 2
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

  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)

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

  if (initial_condition == 'CENTER') then
     colloids% pos(:,1) = solvent_cells% edges/2.0
     colloids% pos(:,2) = solvent_cells% edges/2.0
     colloids% pos(1,2) = colloids% pos(1,2) + d
  else if (initial_condition == 'PLANAR_RANDOM') then
     theta = threefry_double(state(1))*2*pi
     colloids%pos(3,:) = solvent_cells%edges(3)/2
     colloids% pos(1,1) = solvent_cells%edges(1)/2.0 + d*cos(theta)/2
     colloids% pos(2,1) = solvent_cells%edges(2)/2.0 + d*sin(theta)/2
     colloids% pos(1,2) = solvent_cells%edges(1)/2.0 - d*cos(theta)/2
     colloids% pos(2,2) = solvent_cells%edges(2)/2.0 - d*sin(theta)/2
  else
     error stop 'Unknown initial condition IC'
  end if

  call h5gcreate_f(dimer_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, int(300*maxval(sigma)**3))

  skin = 1.5
  n_extra_sorting = 0
  loop_i_last_sort = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  e3 = compute_harmonic()
  if (do_zwall) then
     e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj)
  else
     e_wall = 0
  end if
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  call h5fflush_f(hfile%id, H5F_SCOPE_GLOBAL_F, error)
  write(*,*) 'Running for', equilibration_loops, '+', N_loop, 'loops'
  write(*,*) 'mass', mass 
  call main%tic()
  sampling = .false.
  do i = 0, N_loop+equilibration_loops
     if (i==equilibration_loops) sampling = .true.
     if (modulo(i,20) == 0) write(*,'(i05)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        if ((do_zwall) .and. (solvent_cells%bc(3)/=PERIODIC_BC)) then
           call cell_md_pos_zwall(solvent_cells, solvent, dt, md_flag=.true.)
        else
           call cell_md_pos(solvent_cells, solvent, dt, md_flag=.true.)
        end if

        call varia%tic()
        ! Extra copy for rattle
        if (do_rattle) &
             colloids% pos_rattle = colloids% pos
        if (do_harmonic) &
             e3 = compute_harmonic()
        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) &
                + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(k))
        end do
        if (do_rattle) &
             call rattle_dimer_pos(colloids, d, dt, solvent_cells% edges)

        call varia%tac()

        so_max = cell_maximum_displacement(solvent_cells, solvent, delta_t=dt*(N_MD_steps*i+j - loop_i_last_sort))
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1d0) .or. (so_max >= skin*0.89d0) ) then
           call varia%tic()

           if ((do_zwall) .and. (solvent_cells%bc(2)/=PERIODIC_BC)) then
              call cell_md_pos_zwall(solvent_cells, solvent, &
                   (N_MD_steps*i+j - loop_i_last_sort)*dt, md_flag=.false.)
           else
              call cell_md_pos(solvent_cells, solvent, &
                   (N_MD_steps*i+j - loop_i_last_sort)*dt, md_flag=.false.)
           end if
           call cell_md_vel(solvent_cells, solvent, &
                (N_MD_steps*i+j - loop_i_last_sort)*dt, md_flag=.false.)

           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call varia%tac()
           call solvent% sort(solvent_cells)
           loop_i_last_sort = N_MD_steps*i + j
           call neigh% update_list(colloids, solvent, max_cut + skin, solvent_cells, solvent_colloid_lj)
           call varia%tic()
           !$omp parallel do
           do k = 1, solvent%Nmax
              solvent% pos_old(:,k) = solvent% pos(:,k)
           end do
           call compute_cell_wise_max_v

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
        e3 = compute_harmonic()
        if (do_zwall) e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj)

        call cell_md_vel(solvent_cells, solvent, dt, md_flag=.true.)

        call varia%tic()
        do k=1, colloids% Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(k))
        end do

        if (do_rattle) &
             call rattle_dimer_vel(colloids, d, dt, solvent_cells% edges)
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

     call varia%tic()

     if (sampling) then
        temperature = compute_temperature(solvent, solvent_cells)
        kin_e = (colloids% mass(1)*sum(colloids% vel(:,1)**2) + &
             colloids% mass(2)*sum(colloids% vel(:,2)**2))/2 + &
             sum(solvent% vel**2)/2
        v_com = (sum(solvent% vel, dim=2) + &
             mass(1)*colloids%vel(:,1) + mass(2)*colloids%vel(:,2)) / &
             (solvent%Nmax + mass(1) + mass(2))

        call thermo_data%append(hfile, temperature, e1+e2+e3+e_wall, kin_e, e1+e2+e3+kin_e+e_wall, v_com)

        com_pos = ( colloids%pos(:,1)+colloids%image(:,1)*solvent_cells%edges + &
             colloids%pos(:,2)+colloids%image(:,2)*solvent_cells%edges)/2
        call axial_cf%add(i-equilibration_loops, com_pos, unit_r)

     end if

     call solvent_cells%random_shift(state(1))
     call varia%tac()

     call cell_md_pos(solvent_cells, solvent, ((i+1)*N_MD_steps - loop_i_last_sort)*dt, md_flag=.false.)
     call cell_md_vel(solvent_cells, solvent, ((i+1)*N_MD_steps - loop_i_last_sort)*dt, md_flag=.false.)

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     loop_i_last_sort = N_MD_steps*(i+1)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)
     call varia%tic()
     !$omp parallel do
     do k = 1, solvent%Nmax
        solvent% pos_old(:,k) = solvent% pos(:,k)
     end do
     colloids% pos_old = colloids% pos
     call varia%tac()

     call simple_mpcd_step(solvent, solvent_cells, state)
     call compute_cell_wise_max_v

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  ! create a group for block correlators and write the data

  call h5gcreate_f(hfile%id, 'block_correlators', correlator_group, error)
  call axial_cf%write(correlator_group, N_MD_steps, N_MD_steps*dt, 1, dt)
  call h5gclose_f(correlator_group, error)

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
  call hfile%close()
  call h5close_f(error)

contains

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

  function compute_harmonic() result(e)
    double precision :: e

    double precision :: f(3), x12(3), dist

    e = 0
    f = 0

    x12 = rel_pos(colloids%pos(:,1), colloids%pos(:,2), solvent_cells%edges)
    dist = norm2(x12)

    f = harmonic_k * (dist - d) * x12 / dist
    colloids%force(:,1) = colloids%force(:,1) - f
    colloids%force(:,2) = colloids%force(:,2) + f

    ! U = k/2 * (dist - d)**2
    ! f = k * (dist-d) * d dist / d x

  end function compute_harmonic

end program single_dimer
