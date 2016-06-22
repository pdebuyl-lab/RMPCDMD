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
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer, parameter :: N_species = 2
  integer, parameter :: n_colloids = 36

  integer :: rho
  integer :: N
  integer :: error

  double precision :: max_cut
  double precision :: epsilon(2,2)
  double precision :: sigma, sigma_v(2,2), sigma_cut(2,2)
  double precision :: mass(2)

  double precision :: e1, e2
  double precision :: tau, dt , T
  double precision :: prob
  double precision :: bulk_rate
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting

  type(PTo) :: config

  integer :: i, L(3)
  integer :: j, k
  type(timer_t) :: varia, main
  type(timer_t) :: time_flag, time_refuel, time_change
  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group
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

  integer, allocatable :: links(:,:)
  double precision, allocatable :: links_d(:)
  integer :: i_link
  double precision :: dist, rattle_pos_tolerance, rattle_vel_tolerance

  type(args_t) :: args

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

  call h5open_f(error)

  prob = PTread_d(config,'probability')
  bulk_rate = PTread_d(config,'bulk_rate')

  L = PTread_ivec(config, 'L', 3)
  rho = PTread_i(config, 'rho')
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T')
  
  tau = PTread_d(config, 'tau')
  N_MD_steps = PTread_i(config, 'N_MD')
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop')

  rattle_pos_tolerance = PTread_d(config, 'rattle_pos_tolerance')
  rattle_vel_tolerance = PTread_d(config, 'rattle_vel_tolerance')

  sigma = PTread_d(config, 'sigma_colloid')
  sigma_v = sigma
  sigma_cut = sigma_v*2**(1.d0/6.d0)
  mass = rho * sigma**3 * 4 * 3.14159265/3

  epsilon = PTread_d(config, 'epsilon_colloid')

  call colloid_lj% init(epsilon, sigma_v, sigma_cut)

  ! solvent index first, colloid index second, in solvent_colloid_lj
  sigma = PTread_d(config, 'sigma')
  sigma_v = sigma
  sigma_cut = sigma_v*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  epsilon(:,1) = PTread_dvec(config, 'epsilon_C', 2)
  epsilon(:,2) = PTread_dvec(config, 'epsilon_N', 2)

  call solvent_colloid_lj% init(epsilon, sigma_v, sigma_cut)

  call solvent% init(N,2) !there will be 2 species of solvent particles

  call colloids% init(n_colloids, 2, mass) !there will be 2 species of colloids

  call hfile%create(args%output_file, 'RMPCDMD::single_janus_pbc', &
       'N/A', 'Pierre de Buyl')
  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call PTkill(config)
  
  ! init Janus particle
  colloids%species = 1
  colloids%species(1:n_colloids/2) = 1
  colloids%species(n_colloids/2+1:n_colloids) = 2
  colloids%vel = 0
  colloids%force = 0
  colloids%pos = reshape([13.3745, 8.98251, 4.69069, 8.5896, 7.72357, 4.80796, 9.47026, &
8.68444, 3.05306, 11.5306, 8.71526, 3.21078, 12.0145, 9.40707, 6.13052, 8.17221, 9.88827, &
4.48735, 10.0301, 10.5029, 3.98571, 12.1461, 10.4592, 4.24321, 10.9462, 12.322, 4.66248, &
12.3154, 11.6388, 5.8488, 10.2087, 6.32609, 5.95645, 10.4149, 6.98389, 4.00391, 11.8483, &
7.4706, 5.50245, 10.3726, 8.77586, 4.90677, 10.4163, 8.16153, 6.93092, 13.7445, 10.3013, &
6.58686, 11.7959, 6.50942, 7.31858, 13.4875, 7.93538, 6.6553, 12.7693, 9.16432, 8.15863, &
8.87267, 11.5327, 5.35836, 8.89651, 11.1582, 7.45123, 10.5823, 10.5821, 5.87234, 8.52313, &
7.48241, 6.88098, 7.04389, 8.55093, 5.87505, 11.2238, 7.98061, 8.73575, 7.65272, 9.31661, &
8.08054, 9.40208, 10.2729, 9.34594, 10.6022, 9.87845, 7.7796, 7.34913, 10.6154, 6.3113, &
10.4922, 12.3975, 6.83583, 10.438, 11.8593, 8.7102, 12.2607, 11.1839, 7.82265, 9.29334, &
8.24435, 8.62586, 9.87488, 6.45058, 7.9423, 8.9374, 9.41897, 6.36459, 11.4998, 9.95905, &
9.53615], [3, n_colloids])

  i_link = 0
  do i = 1, colloids%Nmax
     do j = i+1, colloids%Nmax
        x = rel_pos(colloids%pos(:,i), colloids%pos(:,j), solvent_cells%edges)
        dist = norm2(x)
        if (dist < 2.8) then
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
        if (dist < 2.8) then
           ! add link
           i_link = i_link + 1
           links(1, i_link) = i
           links(2, i_link) = j
           links_d(i_link) = dist
        end if
     end do
  end do

  write(*,*) 'number of links:', i_link

  janus_io%force_info%store = .false.
  janus_io%id_info%store = .false.
  janus_io%position_info%store = .true.
  janus_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%position_info%step = N_MD_steps
  janus_io%position_info%time = N_MD_steps*dt
  janus_io%image_info%store = .true.
  janus_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%image_info%step = N_MD_steps
  janus_io%image_info%time = N_MD_steps*dt
  janus_io%velocity_info%store = .true.
  janus_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  janus_io%velocity_info%step = N_MD_steps
  janus_io%velocity_info%time = N_MD_steps*dt
  janus_io%species_info%store = .true.
  janus_io%species_info%mode = H5MD_FIXED
  call janus_io%init(hfile, 'janus', colloids)

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

  call random_number(solvent% vel(:, :))
  solvent% vel = (solvent% vel - 0.5d0)*sqrt(12*T)
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
  call h5gclose_f(fields_group, error)

  call h5gcreate_f(hfile%id, 'connectivity', connectivity_group, error)
  call h5md_write_dataset(connectivity_group, 'janus_links', links-1)
  call h5dopen_f(connectivity_group, 'janus_links', tmp_id, error)
  call h5md_write_attribute(tmp_id, 'particles_group', 'janus')
  call h5dclose_f(tmp_id, error)
  call h5gclose_f(connectivity_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, int(300*sigma**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

  ! insert polar histogram

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  write(*,*) 'Running for', N_loop, 'loops'
  call main%tic()
  do i = 1, N_loop
     if (modulo(i,5) == 0) write(*,'(i05)',advance='no') i
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

        call rattle_body_pos(colloids, links, links_d, dt, solvent_cells% edges, rattle_pos_tolerance)

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
           solvent% pos_old = solvent% pos
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

        call md_vel(solvent, dt)

        call varia%tic()
        do k = 1, colloids%Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids%mass(colloids%species(k)))
        end do
        call varia%tac()

        call rattle_body_vel(colloids, links, links_d, dt, solvent_cells% edges, rattle_pos_tolerance)

        call time_flag%tic()
        call flag_particles
        call time_flag%tac()
        call time_change%tic()
        call change_species
        call time_change%tac()

     end do md_loop

     call varia%tic()

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

     call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e, v_com)

     call janus_io%position%append(colloids%pos)
     call janus_io%velocity%append(colloids%vel)
     call janus_io%image%append(colloids%image)

     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call varia%tac()

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells, solvent_colloid_lj)

     call simple_mpcd_step(solvent, solvent_cells, state)

     call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate, tau, state)

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

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e, v_com, add=.false., force=.true.)
  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

  call janus_io%close()
  call solvent_io%close()
  call radial_hist_el%close()
  call n_solvent_el%close()
  call hfile%close()
  call h5close_f(error)

  write(*,'(a16,f8.3)') solvent%time_stream%name, solvent%time_stream%total
  write(*,'(a16,f8.3)') solvent%time_step%name, solvent%time_step%total
  write(*,'(a16,f8.3)') solvent%time_count%name, solvent%time_count%total
  write(*,'(a16,f8.3)') solvent%time_sort%name, solvent%time_sort%total
  write(*,'(a16,f8.3)') solvent%time_ct%name, solvent%time_ct%total
  write(*,'(a16,f8.3)') solvent%time_md_pos%name, solvent%time_md_pos%total
  write(*,'(a16,f8.3)') solvent%time_md_vel%name, solvent%time_md_vel%total
  write(*,'(a16,f8.3)') solvent%time_max_disp%name, solvent%time_max_disp%total
  write(*,'(a16,f8.3)') colloids%time_self_force%name, colloids%time_self_force%total
  write(*,'(a16,f8.3)') neigh%time_update%name, neigh%time_update%total
  write(*,'(a16,f8.3)') neigh%time_force%name, neigh%time_force%total
  write(*,'(a16,f8.3)') colloids%time_max_disp%name, colloids%time_max_disp%total
  write(*,'(a16,f8.3)') time_flag%name, time_flag%total
  write(*,'(a16,f8.3)') time_change%name, time_change%total
  write(*,'(a16,f8.3)') time_refuel%name, time_refuel%total
  write(*,'(a16,f8.3)') colloids%time_rattle_pos%name, colloids%time_rattle_pos%total
  write(*,'(a16,f8.3)') colloids%time_rattle_vel%name, colloids%time_rattle_vel%total

  write(*,'(a16,f8.3)') 'total', solvent%time_stream%total + solvent%time_step%total + &
       solvent%time_count%total + solvent%time_sort%total + solvent%time_ct%total + &
       solvent%time_md_pos%total + solvent%time_md_vel%total + &
       solvent%time_max_disp%total + solvent%time_self_force%total + &
       neigh%time_update%total + neigh%time_force%total + colloids%time_max_disp%total + &
       time_flag%total + time_change%total + time_refuel%total

  write(*,'(a16,f8.3)') varia%name, varia%total
  write(*,'(a16,f8.3)') main%name, main%total

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

    !$omp parallel do private(x, dist_to_colloid_sq, m_colloid, do_change)
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

end program single_janus_pbc
