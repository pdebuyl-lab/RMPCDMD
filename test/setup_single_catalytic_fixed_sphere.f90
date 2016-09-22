program setup_single_catalytic_fixed_sphere
  use common
  use cell_system
  use particle_system
  use particle_system_io
  use hilbert
  use neighbor_list
  use hdf5
  use h5md_module
  use interaction
  use threefry_module
  use mpcd
  use md
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj

  integer, parameter :: N_species = 2

  integer :: rho
  integer :: N
  integer :: error

  double precision :: epsilon, sigma, sigma_cut
  double precision :: mass

  double precision :: e1
  double precision :: tau, dt , T
  double precision :: prob
  double precision :: skin, so_max
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
  integer(HID_T) :: box_group
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  double precision :: v_com(3)
  type(particle_system_io_t) :: solvent_io

  integer, dimension(N_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el

  type(histogram_t) :: polar_hist
  type(h5md_element_t) :: polar_hist_el
  double precision, parameter :: pi = 4.d0*atan(1.d0)

  logical :: rmpcd_refuel
  double precision :: bulk_rate

  call PTparse(config,get_input_filename(),11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, PTread_c_int64(config, 'seed'))

  call main%init('main')
  call varia%init('varia')
  call time_flag%init('flag')
  call time_refuel%init('refuel')
  call time_change%init('change')

  call h5open_f(error)

  prob = PTread_d(config,'probability')
  rmpcd_refuel = PTread_l(config, 'rmpcd_refuel')
  bulk_rate = PTread_d(config, 'bulk_rate')

  L = PTread_ivec(config, 'L', 3)
  rho = PTread_i(config, 'rho')
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T')
  
  tau = PTread_d(config, 'tau')
  N_MD_steps = PTread_i(config, 'N_MD')
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop')

  epsilon = PTread_d(config, 'epsilon')
  sigma = PTread_d(config, 'sigma')
  sigma_cut = sigma*2**(1.d0/6.d0)
  mass = rho * sigma**3 * 4 * 3.14159265/3

  call solvent_colloid_lj% init(&
       reshape([epsilon, epsilon], [2, 1]),&
       reshape([sigma, sigma], [2, 1]),&
       reshape([sigma_cut, sigma_cut], [2, 1]))

  call solvent% init(N,2) !there will be 2 species of solvent particles

  call colloids% init(1,1, [mass]) !there will be 1 species of colloids

  call hfile%create(PTread_s(config, 'h5md_file'), &
       'RMPCDMD::setup_single_catalytic_fixed_sphere', 'N/A', 'Pierre de Buyl')
  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call PTkill(config)
  
  colloids% species(1) = 1
  colloids% vel = 0

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
  colloids% pos(:,1) = solvent_cells% edges/2.0

  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, H5MD_LINEAR, step=N_MD_steps, &
       time=N_MD_steps*dt)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  call polar_hist%init(0.d0, 3*sigma, 100, 2)
  call polar_hist_el%create_time(fields_group, 'polar_histogram', polar_hist%data, &
       H5MD_LINEAR, step=N_MD_steps, time=N_MD_steps*dt)
  call h5md_write_attribute(polar_hist_el%id, 'xmin', polar_hist%xmin)
  call h5md_write_attribute(polar_hist_el%id, 'dx', polar_hist%dx)
  call h5gclose_f(fields_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, int(300*sigma**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, sigma_cut+skin)

  call neigh% update_list(colloids, solvent, sigma_cut+skin, solvent_cells)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  write(*,*) 'Running for', N_loop, 'loops'
  call main%tic()
  do i = 1, N_loop
     if (modulo(i,5) == 0) write(*,'(i05)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        call md_pos(solvent, dt)

        so_max = solvent% maximum_displacement()

        if ( so_max >= skin ) then
           call varia%tic()
           call apply_pbc(solvent, solvent_cells% edges)
           call varia%tac()
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, sigma_cut + skin, solvent_cells)
           call varia%tic()
           solvent% pos_old = solvent% pos
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

        call varia%tac()
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)

        call md_vel(solvent, dt)

        call time_flag%tic()
        call flag_particles_nl
        call time_flag%tac()
        call time_change%tic()
        call change_species
        call time_change%tac()

     end do md_loop

     call varia%tic()

     temperature = compute_temperature(solvent, solvent_cells)
     kin_e = sum(solvent% vel**2)/2
     v_com = sum(solvent% vel, dim=2) / solvent%Nmax

     call thermo_data%append(hfile, temperature, e1, kin_e, e1+kin_e, v_com)

     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call varia%tac()

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, sigma_cut+skin, solvent_cells)

     call simple_mpcd_step(solvent, solvent_cells, state)

     ! either refuel or bulk reaction
     call time_refuel%tic()
     if (rmpcd_refuel) then
        call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate, tau, state)
     else
        call refuel
     end if
     call time_refuel%tac()

     n_solvent = 0
     do k = 1, solvent%Nmax
        j = solvent%species(k)
        if (j <= 0) continue
        n_solvent(j) = n_solvent(j) + 1
     end do
     call n_solvent_el%append(n_solvent)

     call update_polar_hist

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

  call solvent_io%close()
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
  write(*,'(a16,f8.3)') colloids%time_self_force%name, solvent%time_self_force%total
  write(*,'(a16,f8.3)') neigh%time_update%name, neigh%time_update%total
  write(*,'(a16,f8.3)') neigh%time_force%name, neigh%time_force%total
  write(*,'(a16,f8.3)') colloids%time_max_disp%name, colloids%time_max_disp%total
  write(*,'(a16,f8.3)') time_flag%name, time_flag%total
  write(*,'(a16,f8.3)') time_change%name, time_change%total
  write(*,'(a16,f8.3)') time_refuel%name, time_refuel%total

  write(*,'(a16,f8.3)') 'total', solvent%time_stream%total + solvent%time_step%total + &
       solvent%time_count%total + solvent%time_sort%total + solvent%time_ct%total + &
       solvent%time_md_pos%total + solvent%time_md_vel%total + &
       solvent%time_max_disp%total + solvent%time_self_force%total + &
       neigh%time_update%total + neigh%time_force%total + colloids%time_max_disp%total + &
       time_flag%total + time_change%total + time_refuel%total

  write(*,'(a16,f8.3)') varia%name, varia%total
  write(*,'(a16,f8.3)') main%name, main%total

contains

  subroutine flag_particles_nl
    double precision :: dist_to_C_sq
    integer :: r, s
    double precision :: x(3)

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

  end subroutine flag_particles_nl

  subroutine change_species
    double precision :: dist_to_sphere_sq
    integer :: m
    double precision :: x(3)

    !$omp parallel do private(x, dist_to_sphere_sq)
    do m = 1, solvent% Nmax
       if (solvent% flag(m) == 1) then
          x = rel_pos(colloids% pos(:,1), solvent% pos(:,m), solvent_cells% edges)
          dist_to_sphere_sq = dot_product(x, x)
          if (dist_to_sphere_sq > solvent_colloid_lj%cut_sq(1,1)) then
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

  subroutine update_polar_hist
    integer :: i, s
    double precision :: r, x(3), r_up

    polar_hist%data = 0
    do i = 1, solvent%Nmax
       s = solvent%species(i)
       if (s<=0) cycle
       x = rel_pos(colloids% pos(:,1), solvent% pos(:,i), solvent_cells% edges)
       r = sqrt(sum(x**2))
       call polar_hist%bin(r, s)
    end do

    do i = 1, polar_hist%n
       r = (i-1)*polar_hist%dx
       r_up = i*polar_hist%dx
       polar_hist%data(:,i) = polar_hist%data(:,i) / (4.d0/3.d0*pi*(r_up**3-r**3))
    end do

    call polar_hist_el%append(polar_hist%data)

  end subroutine update_polar_hist

end program setup_single_catalytic_fixed_sphere
