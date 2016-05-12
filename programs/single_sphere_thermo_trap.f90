program setup_sphere_thermo_trap
  use md
  use neighbor_list
  use common
  use cell_system
  use particle_system
  use particle_system_io
  use hilbert
  use interaction
  use hdf5
  use h5md_module
  use particle_system_io
  use mpcd
  use threefry_module
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  integer, parameter :: N_species = 1
  integer, parameter :: N_species_colloids = 1

  type(threefry_rng_t), allocatable :: state(:)

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj

  type(profile_t) :: tz
  type(histogram_t) :: rhoz
  type(profile_t) :: vx
  type(h5md_element_t) :: elem_tz, elem_tz_count
  type(h5md_element_t) :: elem_rhoz
  double precision, allocatable :: v_xz(:,:,:)
  integer, allocatable :: v_xz_count(:,:)
  type(h5md_element_t) :: v_xz_el

  integer :: rho
  integer :: N
  integer :: error

  double precision :: epsilon(N_species,N_species_colloids)
  double precision :: sigma(N_species,N_species_colloids), sigma_cut(N_species,N_species_colloids)
  double precision :: mass(N_species_colloids)
  double precision :: max_cut

  double precision :: v_com(3), wall_v(3,2), wall_t(2)

  double precision :: e1, e2
  double precision :: tau, dt , T, alpha
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: N_therm
  integer :: n_extra_sorting
  double precision :: kin_e, temperature

  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group
  type(h5md_element_t) :: vx_el
  type(thermo_t) :: thermo_data
  type(particle_system_io_t) :: sphere_io
  type(particle_system_io_t) :: solvent_io
  integer(HID_T) :: box_group

  type(PTo) :: config
  integer :: L(3), n_threads
  integer :: i, j

  type(timer_t) :: flag_timer, change_timer, varia
  integer(HID_T) :: timers_group

  !gravity
  double precision :: g
  ! trap parameters
  double precision :: k, trap_center(3)

  call PTparse(config,get_input_filename(),11)

  call flag_timer%init('flag')
  call change_timer%init('change')
  call varia%init('varia')

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, PTread_c_int64(config, 'seed'))

  call h5open_f(error)

  L = PTread_ivec(config, 'L', 3)
  if (modulo(L(2),2) /= 0) error stop 'non-even Ly is not supported'
  rho = PTread_i(config, 'rho')
  N = rho *L(1)*L(2)*L(3)
  tau =PTread_d(config, 'tau')
  alpha = PTread_d(config,'alpha')

  N_MD_steps = PTread_i(config, 'N_MD')
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop')
  N_therm = PTread_i(config, 'N_therm')

  wall_t = PTread_dvec(config, 'wall_T', 2)
  T = PTread_d(config, 'T')

  g = PTread_d(config, 'g')
  k = PTread_d(config, 'k')
  
  sigma = PTread_d(config, 'sigma')
  epsilon = PTread_d(config, 'epsilon')

  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)

  mass(1) = rho * sigma(1,1)**3 * 4 * 3.14159265/3

  call solvent% init(N,N_species)

  call colloids% init(1, N_species_colloids, mass)
  colloids% species(1) = 1
  colloids% vel = 0

  call hfile%create(PTread_s(config, 'h5md_file'), 'RMPCDMD::single_sphere_thermo_trap', &
       'N/A', 'Pierre de Buyl')
  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call PTkill(config)

  sphere_io%id_info%store = .false.
  sphere_io%force_info%store = .true.
  sphere_io%force_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  sphere_io%force_info%step = N_MD_steps
  sphere_io%force_info%time = N_MD_steps*dt
  sphere_io%position_info%store = .true.
  sphere_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  sphere_io%position_info%step = N_MD_steps
  sphere_io%position_info%time = N_MD_steps*dt
  sphere_io%image_info%store = .true.
  sphere_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  sphere_io%image_info%step = N_MD_steps
  sphere_io%image_info%time = N_MD_steps*dt
  sphere_io%velocity_info%store = .true.
  sphere_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  sphere_io%velocity_info%step = N_MD_steps
  sphere_io%velocity_info%time = N_MD_steps*dt
  sphere_io%species_info%store = .true.
  sphere_io%species_info%mode = H5MD_FIXED
  call sphere_io%init(hfile, 'sphere', colloids)

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
  call solvent_cells%init(L, 1.d0, has_walls=.true.)

  trap_center = solvent_cells%edges/2
  colloids%pos(:,1) = trap_center
  colloids%vel = 0
  colloids%force = 0

  call vx% init(0.d0, solvent_cells% edges(3), L(3))
  call tz% init(0.d0, solvent_cells% edges(3), L(3))
  call rhoz% init(0.d0, solvent_cells% edges(3), L(3))

  allocate(v_xz_count(L(3), L(1)))
  allocate(v_xz(2, L(3), L(1)))

  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  call vx_el%create_time(fields_group, 'vx', vx%data, ior(H5MD_LINEAR,H5MD_STORE_TIME), &
       step=N_MD_steps, time=N_MD_steps*dt)
  call elem_tz% create_time(fields_group, 'tz', tz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_tz_count% create_time(fields_group, 'tz_count', tz% count, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_rhoz% create_time(fields_group, 'rhoz', rhoz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call v_xz_el%create_time(fields_group, 'v_xz', v_xz, ior(H5MD_TIME, H5MD_STORE_TIME))
  call h5gclose_f(fields_group, error)

  call h5gcreate_f(sphere_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)
  solvent% pos_old = solvent% pos

  do i=1, solvent% Nmax
     solvent% vel(1,i) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,i) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,i) = threefry_normal(state(1))*sqrt(T)
  end do

  call solvent% sort(solvent_cells)

  skin = 1
  call neigh% init(colloids% Nmax, 25*floor((max_cut+skin)**2)*rho)

  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin+0.1)
  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

  solvent% force = 0
  colloids% force = 0
  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_harmonic_trap(colloids, k, trap_center)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  i = 0
  wall_v = 0

  write(*,*) 'Running for', N_loop, 'loops'
  !start RMPCDMD
  setup: do i = 1, N_loop
     if (modulo(i,50) == 0) write(*,'(i09)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        call mpcd_stream_nogravity_zwall(solvent, solvent_cells, dt)

        colloids% pos = colloids% pos + dt * colloids% vel + &
             dt**2 * colloids% force / (2*colloids% mass(1))

        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1) .or. (so_max >= skin*0.9) ) then
           call apply_pbc(colloids, solvent_cells% edges)
           call apply_pbc(solvent, solvent_cells% edges)
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)
           call varia%tic()
           solvent% pos_old = solvent% pos
           colloids% pos_old = colloids% pos
           call varia%tac()
           n_extra_sorting = n_extra_sorting + 1
        end if

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)

        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_harmonic_trap(colloids, k, trap_center)

        call md_vel(solvent, dt)
        colloids% vel = colloids% vel + &
             dt * ( colloids% force + colloids% force_old ) / (2 * colloids% mass(1))

     end do md_loop

     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1

     call apply_pbc(colloids, solvent_cells% edges)
     call apply_pbc(solvent, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)
     solvent%pos_old = solvent% pos
     colloids%pos_old = colloids% pos

     call rescale_at_walls
     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[rho, rho], alpha=alpha)

     if (i>N_therm) then
        temperature = compute_temperature(solvent, solvent_cells, tz)
        kin_e = colloids% mass(1)*sum(colloids% vel(:,1)**2)/2 + &
             sum(solvent% vel**2)/2
        v_com = (sum(solvent% vel, dim=2) + mass(1)*colloids%vel(:,1)) / (solvent%Nmax + mass(1))
        call thermo_data%append(hfile, temperature, e1, kin_e, e1+kin_e, v_com)
        call compute_rho(solvent, rhoz)

        call varia%tic()
        call compute_vx(solvent, vx)
        call vx% norm()
        call vx_el%append(vx%data)
        call vx% reset()
        call tz% norm()
        call elem_tz% append(tz% data, i, i*tau)
        call elem_tz_count% append(tz% count, i, i*tau)
        call tz% reset()
        rhoz% data = rhoz% data / rhoz% dx
        call elem_rhoz% append(rhoz% data, i, i*tau)
        rhoz% data = 0
        call compute_vxz
        call v_xz_el%append(v_xz, i, i*tau)
        call varia%tac()

        call sphere_io%position%append(colloids%pos)
        call sphere_io%velocity%append(colloids%vel)
        call sphere_io%force%append(colloids%force)
        call sphere_io%image%append(colloids%image)
     end if

  end do setup

  call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e, v_com, add=.false., force=.true.)

  write(*,*) ''
  write(*,*) 'n extra sorting', n_extra_sorting

  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

  call h5gcreate_f(hfile%id, 'timers', timers_group, error)
  call h5md_write_dataset(timers_group, solvent%time_stream%name, solvent%time_stream%total)
  call h5md_write_dataset(timers_group, solvent%time_md_vel%name, solvent%time_md_vel%total)
  call h5md_write_dataset(timers_group, solvent%time_step%name, solvent%time_step%total)
  call h5md_write_dataset(timers_group, solvent%time_count%name, solvent%time_count%total)
  call h5md_write_dataset(timers_group, solvent%time_sort%name, solvent%time_sort%total)
  call h5md_write_dataset(timers_group, solvent%time_ct%name, solvent%time_ct%total)
  call h5md_write_dataset(timers_group, solvent%time_max_disp%name, solvent%time_max_disp%total)
  call h5md_write_dataset(timers_group, solvent%time_apply_pbc%name, solvent%time_apply_pbc%total)
  call h5md_write_dataset(timers_group, neigh%time_update%name, neigh%time_update%total)
  call h5md_write_dataset(timers_group, varia%name, varia%total)
  call h5md_write_dataset(timers_group, neigh%time_force%name, neigh%time_force%total)

  call h5md_write_dataset(timers_group, 'total', solvent%time_stream%total + &
       solvent%time_step%total + solvent%time_count%total + solvent%time_sort%total + &
       solvent%time_ct%total + solvent%time_md_vel%total + solvent%time_max_disp%total + &
       flag_timer%total + change_timer%total + solvent%time_apply_pbc%total+ &
       neigh%time_update%total + varia%total + neigh%time_force%total)

  call h5gclose_f(timers_group, error)

  call elem_tz%close()
  call elem_tz_count%close()
  call elem_rhoz%close()
  call sphere_io%close()
  call hfile%close()
  call h5close_f(error)

contains

  function compute_force_harmonic_trap(p, k, x0) result(e)
    type(particle_system_t), intent(inout) :: p
    double precision, intent(in) :: k, x0(3)
    double precision :: e

    integer :: i

    e = 0
    do i = 1, p%Nmax
       p%force(:,i) = p%force(:,i) - k * (p%pos(:,i)-x0)
       e = e + k * sum( (p%pos(:,i)-x0)**2 ) / 2
    end do

  end function compute_force_harmonic_trap

  subroutine rescale_at_walls

    integer :: i
    integer :: cell_idx, n, start, wall_idx, cell(3)
    double precision :: local_v(3), local_k, local_T, factor

    do cell_idx = 1, solvent_cells% N
       if (solvent_cells% cell_count(cell_idx) <= 1) cycle

       start = solvent_cells% cell_start(cell_idx)
       n = solvent_cells% cell_count(cell_idx)

       ! Find whether we are in a wall cell
       cell = compact_h_to_p(cell_idx - 1, solvent_cells% M) + 1
       if (cell(3) == 1) then
          wall_idx = 1
       else if (cell(3) == solvent_cells% L(3)) then
          wall_idx = 2
          local_T = 1.1
       else
          wall_idx = -1
       end if
       if (wall_idx==-1) cycle
       local_T = wall_t(wall_idx)

       local_v = 0
       do i = start, start + n - 1
          local_v = local_v + solvent% vel(:, i)
       end do
       local_v = local_v / n

       local_k = 0
       do i = start, start + n - 1
          local_k = local_k + sum((solvent% vel(:, i)-local_v)**2)/2
       end do
       local_k = local_k/(3*(n-1))
       factor = sqrt(local_T/(2*local_k))

       do i = start, start + n - 1
          solvent% vel(:, i) = local_v + factor*(solvent% vel(:, i)-local_v)
       end do

    end do

  end subroutine rescale_at_walls

  subroutine compute_vxz
    integer :: i, j, s, ix, iy, iz

    v_xz_count = 0
    v_xz = 0
    do i = 1, solvent%Nmax
       s = solvent%species(i)
       if (s <= 0) cycle
       ix = modulo(floor(solvent%pos(1,i)/solvent_cells%a), L(1)) + 1
       iy = modulo(floor(solvent%pos(2,i)/solvent_cells%a), L(2)) + 1
       iz = modulo(floor(solvent%pos(3,i)/solvent_cells%a), L(3)) + 1
       if ( (iy == L(2)/2) .or. (iy == 1+L(2)/2) ) then
          v_xz_count(iz, ix) = v_xz_count(iz, ix) + 1
          v_xz(1, iz, ix) = v_xz(1, iz, ix) + solvent%vel(1, i)
          v_xz(2, iz, ix) = v_xz(2, iz, ix) + solvent%vel(3, i)
       end if
    end do

    do i = 1, L(1)
       do j = 1, L(3)
          if (v_xz_count(j, i) > 0) then
             v_xz(:, j, i) = v_xz(:, j, i) / v_xz_count(j, i)
          end if
       end do
    end do

  end subroutine compute_vxz

end program setup_sphere_thermo_trap
