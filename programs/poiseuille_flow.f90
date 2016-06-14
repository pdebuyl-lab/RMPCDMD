!> Simulate a forced flow between two plates
!!
!! Consider a pure fluid under a constant acceleration in the x-direction. Bounce-back
!! boundary conditions are used in the z-direction in addition to ghost cells for the
!! collisions near the walls.
!!
!! \param L           length of simulation box in the 3 dimensions
!! \param g           strength of the constant acceleration in x
!! \param rho         fluid number density
!! \param T           Temperature. Used for setting initial velocities, for wall
!!                    thermostatting and (if enabled) bulk thermostatting.
!! \param tau         MPCD collision time
!! \param alpha       MPCD collision angle
!! \param thermostat  whether to enable bulk thermostatting
!! \param N_therm     number of unsampled thermalization MPCD timesteps
!! \param N_loop      number of MPCD timesteps

program poiseuille_flow
  use common
  use cell_system
  use particle_system
  use hilbert
  use interaction
  use hdf5
  use h5md_module
  use particle_system_io
  use mpcd
  use md
  use threefry_module
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  type(threefry_rng_t), allocatable :: state(:)
  type(PTo) :: config

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent

  type(profile_t) :: tz
  type(histogram_t) :: rhoz
  type(profile_t) :: vx

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: elem_tz, elem_tz_count, elem_vx_count
  type(h5md_element_t) :: elem_rhoz
  type(h5md_element_t) :: elem_vx
  type(h5md_element_t) :: elem_T
  type(h5md_element_t) :: elem_v_com
  integer(HID_T) :: box_group, solvent_group
  integer(HID_T) :: fields_group
  type(particle_system_io_t) :: solvent_io

  integer :: i, L(3), error, N, n_threads
  integer :: rho
  integer :: N_loop, N_therm

  double precision :: v_com(3), wall_v(3,2), wall_t(2)
  double precision :: gravity_field(3)
  double precision :: T, set_temperature, tau
  double precision :: alpha
  logical :: thermostat
  type(args_t) :: args

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call h5open_f(error)

  L = PTread_ivec(config, 'L', 3)
  rho = PTread_i(config, 'rho')
  N = rho *L(1)*L(2)*L(3)

  set_temperature = PTread_d(config, 'T')
  thermostat = PTread_l(config, 'thermostat')
  
  tau = PTread_d(config, 'tau')
  alpha = PTread_d(config, 'alpha')
  N_therm = PTread_i(config, 'N_therm')
  N_loop = PTread_i(config, 'N_loop')
  gravity_field = 0
  gravity_field(1) = PTread_d(config, 'g')

  call solvent% init(N)

  do i=1, solvent% Nmax
     solvent% vel(1,i) = threefry_normal(state(1))
     solvent% vel(2,i) = threefry_normal(state(1))
     solvent% vel(3,i) = threefry_normal(state(1))
  end do
  solvent%vel = solvent%vel*sqrt(set_temperature)
  v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
  solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))

  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0)

  call solvent_cells%init(L, 1.d0, has_walls=.true.)
  solvent_cells% origin(3) = -0.5d0

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create(args%output_file, 'RMPCDMD:poiseuille_flow', &
       'N/A', 'Pierre de Buyl')

  call PTkill(config)

  call tz% init(0.d0, solvent_cells% edges(3), L(3))
  call rhoz% init(0.d0, solvent_cells% edges(3), L(3))
  call vx% init(0.d0, solvent_cells% edges(3), L(3))

  solvent_io%force_info%store = .false.
  solvent_io%species_info%store = .false.
  solvent_io%position_info%store = .true.
  solvent_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  solvent_io%position_info%step = 1
  solvent_io%position_info%time = tau
  solvent_io%image_info = solvent_io%position_info
  solvent_io%velocity_info = solvent_io%position_info
  solvent_io%id_info = solvent_io%position_info
  call solvent_io%init(datafile, 'solvent', solvent)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call h5md_write_attribute(box_group, 'boundary', ['periodic', 'periodic', 'periodic'])
  call elem% create_fixed(box_group, 'edges', L*1.d0)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(datafile%id, 'fields', fields_group, error)
  call elem_tz% create_time(fields_group, 'tz', tz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_tz_count% create_time(fields_group, 'tz_count', tz% count, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_vx% create_time(fields_group, 'vx', tz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_vx_count% create_time(fields_group, 'vx_count', vx% count, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_rhoz% create_time(fields_group, 'rhoz', rhoz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_T% create_time(datafile% observables, 'temperature', T, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_v_com% create_time(datafile% observables, 'center_of_mass_velocity', v_com, ior(H5MD_TIME, H5MD_STORE_TIME))
  call h5gclose_f(fields_group, error)

  call solvent% sort(solvent_cells)

  solvent%force(1,:) = gravity_field(1)
  solvent%force_old(1,:) = gravity_field(1)

  wall_v = 0
  wall_t = set_temperature
  solvent_cells%bc = [ PERIODIC_BC, PERIODIC_BC, BOUNCE_BACK_BC ]
  do i = 1, N_therm
     call mpcd_stream_xforce_yzwall(solvent, solvent_cells, tau, gravity_field(1))
     call md_vel(solvent, tau)
     call apply_pbc(solvent, solvent_cells%edges)
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1
     call solvent% sort(solvent_cells)
     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[rho, rho], thermostat=thermostat, &
          bulk_temperature=set_temperature, alpha=alpha)
  end do

  do i = 1, N_loop
     if (mod(i,100)==0) then
        write(*,*) i 
     end if
     call mpcd_stream_xforce_yzwall(solvent, solvent_cells, tau, gravity_field(1))
     call md_vel(solvent, tau)
     call apply_pbc(solvent, solvent_cells%edges)
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1
     call solvent% sort(solvent_cells)
     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[rho, rho], thermostat=thermostat, &
          bulk_temperature=set_temperature, alpha=alpha)
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)


     T = compute_temperature(solvent, solvent_cells, tz)
     call elem_v_com%append(v_com, i, i*tau)
     call elem_T% append(T, i, i*tau)

     call compute_rho(solvent, rhoz)
     call compute_vx(solvent, vx)

     if (modulo(i, 50) == 0) then
        call tz% norm()
        call elem_tz% append(tz% data, i, i*tau)
        call elem_tz_count% append(tz% count, i, i*tau)
        call tz% reset()
        call vx% norm()
        call elem_vx% append(vx% data, i, i*tau)
        call elem_vx_count% append(vx% count, i, i*tau)
        call vx% reset()
        rhoz% data = rhoz% data / (50.d0 * rhoz% dx)
        call elem_rhoz% append(rhoz% data, i, i*tau)
        rhoz% data = 0
     end if

     if (modulo(i,100)==0) then
        call solvent_io%position%append(solvent% pos)
        call solvent_io%image%append(solvent% image)
        call solvent_io%id%append(solvent% id)
        call solvent_io%velocity%append(solvent% vel)
        call h5fflush_f(datafile%id, H5F_SCOPE_GLOBAL_F, error)
     end if

  end do

  call elem_tz% close()
  call elem_tz_count% close()
  call elem_rhoz% close()
  call elem_vx% close()
  call elem_vx_count% close()
  call elem_T% close()
  call elem_v_com% close()

  call datafile% close()

  call h5close_f(error)

  write(*,'(a16,f8.3)') solvent%time_stream%name, solvent%time_stream%total
  write(*,'(a16,f8.3)') solvent%time_step%name, solvent%time_step%total
  write(*,'(a16,f8.3)') solvent%time_count%name, solvent%time_count%total
  write(*,'(a16,f8.3)') solvent%time_sort%name, solvent%time_sort%total
  write(*,'(a16,f8.3)') solvent%time_ct%name, solvent%time_ct%total
  write(*,'(a16,f8.3)') 'total                          ', &
       solvent%time_stream%total + solvent%time_step%total + solvent%time_count%total +&
       solvent%time_sort%total + solvent%time_ct%total

end program poiseuille_flow
