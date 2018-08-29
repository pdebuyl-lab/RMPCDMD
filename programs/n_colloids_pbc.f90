! This file is part of RMPCDMD
! Copyright (c) 2015-2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Simulate an ensemble of spherical colloids
!!
!! The periodic simulation box is filled with a number of spherical colloids, that interact
!! with an attractive Lennard-Jones potential, and with solvent particles.
!! The temperature can be controlled with the MPCD Anderson thermostat.
!!
!! \param L                length of simulation box in the 3 dimensions
!! \param rho              fluid number density
!! \param T                Temperature. Used for setting initial velocities and for thermostatting.
!! \param T_final          Target temperature. Used for thermostatting with temperature program from T to T_final.
!! \param tau              MPCD collision time
!! \param do_thermostat    enable MPCD-AT thermostat
!! \param do_hydro         conserve cell-wise momentum (can be turned off only with thermostat enabled)
!! \param N_MD             number MD steps occuring in tau
!! \param colloid_sampling interval (in MD steps) of sampling the colloid position and velocity
!! \param N_loop           number of MPCD timesteps
!! \param N_colloids       number of colloids
!! \param epsilon          solvent-colloid epsilon
!! \param sigma            radius of the colloids
!! \param epsilon_colloids colloid-colloid epsilon

program n_colloids_pbc
  use rmpcdmd_module
  use hdf5
  use h5md_module
  use threefry_module
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  type(PTo) :: config

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer :: rho
  double precision :: sigma, sigma_cut, epsilon
  integer :: error

  double precision :: mass
  double precision :: so_max, co_max

  double precision :: e1, e2
  double precision :: temperature, local_mass, total_mass, kin_e, v_com(3)
  double precision :: tau, dt, T, T_init, T_final
  integer :: N_MD_steps, N_loop
  integer :: colloid_sampling
  integer :: n_extra_sorting
  integer :: n_threads
  logical :: do_hydro, do_thermostat

  type(threefry_rng_t), allocatable :: state(:)
  type(h5md_file_t) :: hfile
  integer(HID_T) :: params_group
  type(particle_system_io_t) :: colloids_io
  type(thermo_t) :: thermo_data
  integer(HID_T) :: box_group
  type(h5md_element_t) :: dummy_element

  double precision :: skin
  double precision :: rsq
  logical :: tooclose

  double precision :: total_time
  type(timer_list_t) :: timer_list
  integer(HID_T) :: timers_group

  integer :: i, L(3), N, N_colloids
  integer :: j, k
  type(args_t) :: args

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call timer_list%init(11)

  call h5open_f(error)
  call hfile%create(args%output_file, 'RMPCDMD::n_colloids_pbc', &
       RMPCDMD_REVISION, 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)

  L = PTread_ivec(config, 'L', 3, loc=params_group)
  rho = PTread_i(config, 'rho', loc=params_group)
  T_init = PTread_d(config, 'T', loc=params_group)
  T = T_init
  T_final = PTread_d(config, 'T_final', loc=params_group)
  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)

  N_colloids = PTread_i(config, 'N_colloids', loc=params_group)
  epsilon = PTread_d(config, 'epsilon', loc=params_group)
  sigma = PTread_d(config, 'sigma', loc=params_group)
  sigma_cut = sigma*2.d0**(1.d0/6.d0)

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  do_hydro = PTread_l(config, 'do_hydro', loc=params_group)
  do_thermostat = PTread_l(config, 'do_thermostat', loc=params_group)

  mass = rho * sigma**3 * 4 * 3.141/3
  call colloids%init(N_colloids, 1, [mass])
  colloids%species = 1
  colloids%vel = 0
  colloids%force = 0

  colloids_io%force_info%store = .false.
  colloids_io%id_info%store = .false.
  colloids_io%velocity_info%store = .true.
  colloids_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  colloids_io%velocity_info%step = N_MD_steps*colloid_sampling
  colloids_io%velocity_info%time = N_MD_steps*colloid_sampling*dt
  colloids_io%position_info%store = .true.
  colloids_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  colloids_io%position_info%step = N_MD_steps*colloid_sampling
  colloids_io%position_info%time = N_MD_steps*colloid_sampling*dt
  colloids_io%image_info%store = .true.
  colloids_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  colloids_io%image_info%step = N_MD_steps*colloid_sampling
  colloids_io%image_info%time = N_MD_steps*colloid_sampling*dt
  colloids_io%species_info%store = .true.
  colloids_io%species_info%mode = H5MD_FIXED
  call colloids_io%init(hfile, 'colloids', colloids)

  N = rho*L(1)*L(2)*L(3) - int(rho*4*3.142/3 * sigma**3*colloids%Nmax)

  write(*,*) N, 'solvent particles'

  if (N <= 0) error stop 'Not enough volume available for solvent'

  call solvent_colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  epsilon = PTread_d(config, 'epsilon_colloids', loc=params_group)
  sigma_cut = 3*sigma

  call h5gclose_f(params_group, error)
  call PTkill(config)

  call colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ 2*sigma ], [1, 1] ), reshape( [ 2*sigma_cut ], [1, 1] ) )

  call solvent% init(N)
  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)
  call h5gcreate_f(colloids_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  i = 1
  place_colloids: do while (i<=N_colloids)
     colloids%pos(1,i) = threefry_double(state(1))*L(1)
     colloids%pos(2,i) = threefry_double(state(1))*L(2)
     colloids%pos(3,i) = threefry_double(state(1))*L(3)

     tooclose = .false.
     check_distance: do j = 1, i-1
        rsq = sum(rel_pos(colloids%pos(:,i), colloids%pos(:,j), solvent_cells%edges)**2)
        if (rsq < colloid_lj%sigma(1,1)**2) then
           tooclose = .true.
           exit check_distance
        end if
     end do check_distance

     if (.not. tooclose) i=i+1
  end do place_colloids

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  sigma_cut = sigma*2.d0**(1.d0/6.d0) ! restore because colloid sigma_cut is larger
  skin = 0.8
  call neigh% init(colloids% Nmax, int(200*(sigma+skin)**2))
  call neigh% make_stencil(solvent_cells, sigma_cut+skin)
  call neigh% update_list(colloids, solvent, sigma_cut+skin, solvent_cells)

  n_extra_sorting = 0

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  solvent% pos_old = solvent% pos
  colloids% pos_old = colloids% pos
  write(*,*) 'Running for', N_loop, 'loops'
  do i = 1, N_loop
     if (modulo(i, 100)==0) write(*, '(i09)', advance='no') i
     md_loop: do j = 1, N_MD_steps
        call md_pos(solvent, dt)
        do k = 1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + dt**2 * colloids% force(:,k) / (2 * mass)
        end do
        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1) .or. (so_max >= skin*0.9) ) then
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, sigma_cut+skin, solvent_cells)
           solvent% pos_old = solvent% pos
           colloids% pos_old = colloids% pos
           n_extra_sorting = n_extra_sorting + 1
        end if

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)
        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)

        call md_vel(solvent, dt)
        colloids% vel = colloids% vel + dt * ( colloids% force + colloids% force_old ) / (2 * mass)

     end do md_loop

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, sigma_cut + skin, solvent_cells)
     solvent% pos_old = solvent% pos
     colloids% pos_old = colloids% pos

     T = T_init + (i-1)*(T_final-T_init)/(N_loop-1)
     call simple_mpcd_step(solvent, solvent_cells, state, &
          thermostat=do_thermostat, T=T, hydro=do_hydro)

     if (modulo(i,colloid_sampling)==0) then
        call colloids_io%velocity%append(colloids%vel)
        call colloids_io%position%append(colloids%pos)
        call colloids_io%image%append(colloids%image)
     end if

     temperature = compute_temperature(solvent, solvent_cells)
     total_mass = 0
     kin_e = sum(solvent% vel**2)/2
     v_com = sum(solvent% vel, dim=2)
     do k = 1, colloids% Nmax
        j = colloids%species(k)
        if (j==0) cycle
        local_mass = colloids%mass(j)
        kin_e = kin_e + local_mass*sum(colloids% vel(:,k)**2)/2
        v_com = v_com + local_mass*colloids%vel(:,k)
        total_mass = total_mass + local_mass
     end do
     v_com = v_com / (solvent%Nmax + total_mass)

     call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e, v_com)

  end do

  call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e, &
       v_com, add=.false., force=.true.)

  write(*,*) 'n extra sorting', n_extra_sorting

  call h5gcreate_f(hfile%id, 'timers', timers_group, error)
  call timer_list%append(solvent%time_stream)
  call timer_list%append(solvent%time_md_vel)
  call timer_list%append(solvent%time_step)
  call timer_list%append(solvent%time_count)
  call timer_list%append(solvent%time_sort)
  call timer_list%append(solvent%time_ct)
  call timer_list%append(solvent%time_max_disp)
  call timer_list%append(solvent%time_apply_pbc)
  call timer_list%append(neigh%time_update)
  call timer_list%append(neigh%time_force)
  call timer_list%append(colloids%time_self_force)

  call timer_list%write(timers_group, total_time)

  call h5md_write_dataset(timers_group, 'total', total_time)

  call h5gclose_f(timers_group, error)

  call colloids_io%close()
  call hfile%close()
  call h5close_f(error)

end program n_colloids_pbc
