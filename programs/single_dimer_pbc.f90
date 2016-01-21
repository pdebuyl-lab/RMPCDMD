program setup_single_dimer
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
  type(lj_params_t) :: colloid_lj

  integer :: rho
  integer :: N
  integer :: error

  double precision :: sigma_N, sigma_C, max_cut
  double precision :: epsilon(2,2)
  double precision :: sigma(2,2), sigma_cut(2,2)
  double precision :: mass(2)

  double precision :: e1, e2
  double precision :: tau, dt , T
  double precision :: d,prob
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting

  double precision :: conc_z(400)
  double precision :: colloid_pos(3,2)

  type(PTo) :: config

  integer :: i, L(3)
  integer :: j, k
  type(timer_t) :: varia, main
  type(timer_t) :: time_flag, time_refuel, time_change
  type(threefry_rng_t), allocatable :: state(:)
  integer(c_int64_t) :: seed
  integer :: n_threads
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: box_group
  type(thermo_t) :: thermo_data
  double precision :: temperature, kin_e
  type(particle_system_io_t) :: dimer_io

  call PTparse(config,get_input_filename(),11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  seed = PTread_c_int64(config, 'seed')
  do i = 1, n_threads
     state(i)%counter%c0 = 0
     state(i)%counter%c1 = 0
     state(i)%key%c0 = 0
     state(i)%key%c1 = seed
  end do

  call main%init('main')
  call varia%init('varia')
  call time_flag%init('flag')
  call time_refuel%init('refuel')
  call time_change%init('change')

  call h5open_f(error)

  prob = PTread_d(config,'probability')

  L = PTread_ivec(config, 'L', 3)
  rho = PTread_i(config, 'rho')
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T')
  d = PTread_d(config, 'd')
  
  tau = PTread_d(config, 'tau')
  N_MD_steps = PTread_i(config, 'N_MD')
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop')

  sigma_C = PTread_d(config, 'sigma_C')
  sigma_N = PTread_d(config, 'sigma_N')
  
  epsilon(1,:) = PTread_dvec(config, 'epsilon_C', 2)
  epsilon(2,:) = PTread_dvec(config, 'epsilon_N', 2)

  sigma(1,:) = sigma_C
  sigma(2,:) = sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)

  epsilon(1,1) = PTread_d(config, 'epsilon_C_C')
  epsilon(1,2) = PTread_d(config, 'epsilon_N_C')
  epsilon(2,1) = PTread_d(config, 'epsilon_N_C')
  epsilon(2,2) = PTread_d(config, 'epsilon_N_N')

  sigma(1,1) = 2*sigma_C
  sigma(1,2) = sigma_C + sigma_N
  sigma(2,1) = sigma_C + sigma_N
  sigma(2,2) = 2*sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init(epsilon, sigma, sigma_cut)

  mass(1) = rho * sigma_C**3 * 4 * 3.14159265/3
  mass(2) = rho * sigma_N**3 * 4 * 3.14159265/3

  call solvent% init(N,2) !there will be 2 species of solvent particles

  call colloids% init(2,2, mass) !there will be 2 species of colloids

  call hfile%create(PTread_s(config, 'h5md_file'), 'RMPCDMD::single_dimer_pbc', &
       'N/A', 'Pierre de Buyl')
  call thermo_data%init(hfile, n_buffer=50, step=N_loop, time=N_MD_steps*dt)

  call PTkill(config)
  
  open(16,file ='dimerdata_Conc.txt')
  
  colloids% species(1) = 1
  colloids% species(2) = 2
  colloids% vel = 0

  dimer_io%force_info%store = .false.
  dimer_io%id_info%store = .false.
  dimer_io%position_info%store = .true.
  dimer_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%position_info%step = N_loop
  dimer_io%position_info%time = N_loop*dt
  dimer_io%image_info%store = .true.
  dimer_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%image_info%step = N_loop
  dimer_io%image_info%time = N_loop*dt
  dimer_io%velocity_info%store = .true.
  dimer_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%velocity_info%step = N_loop
  dimer_io%velocity_info%time = N_loop*dt
  dimer_io%species_info%store = .true.
  dimer_io%species_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  call dimer_io%init(hfile, 'dimer', colloids)

  call random_number(solvent% vel(:, :))
  solvent% vel = (solvent% vel - 0.5d0)*sqrt(12*T)
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)
  colloids% pos(:,1) = solvent_cells% edges/2.0
  colloids% pos(:,2) = solvent_cells% edges/2.0 
  colloids% pos(1,2) = colloids% pos(1,2) + d

  call h5gcreate_f(dimer_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, int(300*max(sigma_C,sigma_N)**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

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
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(k))
        end do

        call rattle_dimer_pos(colloids, d, dt, solvent_cells% edges)
        call varia%tac()

        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin/2) .or. (so_max >= skin/2) ) then
           call varia%tic()
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call varia%tac()
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, max_cut + skin, solvent_cells)
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

        call md_vel(solvent, solvent_cells% edges, dt)

        call varia%tic()
        do k=1, colloids% Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(k))
        end do

        call rattle_dimer_vel(colloids, d, dt, solvent_cells% edges)
        call varia%tac()

        call time_flag%tic()
        call flag_particles_nl
        call time_flag%tac()
        call time_change%tic()
        call change_species
        call time_change%tac()

     end do md_loop

     call varia%tic()

     temperature = compute_temperature(solvent, solvent_cells)
     kin_e = (colloids% mass(1)*sum(colloids% vel(:,1)**2) + &
          colloids% mass(2)*sum(colloids% vel(:,2)**2))/2 + &
          sum(solvent% vel**2)/2

     call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e)

     call dimer_io%position%append(colloids%pos)
     call dimer_io%velocity%append(colloids%vel)
     call dimer_io%image%append(colloids%image)

     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call varia%tac()

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

     call simple_mpcd_step(solvent, solvent_cells, state)
     
     call time_refuel%tic()
     call refuel
     call time_refuel%tac()
     
     if (modulo(i,100)==0) then
        call concentration_field
        write(16,*) conc_z, colloid_pos
     end if

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

  call dimer_io%close()
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

  subroutine flag_particles
  double precision :: dist_to_C_sq
  integer :: r
  double precision :: x(3)
  integer :: thread_id
  
  !$omp parallel
  thread_id = omp_get_thread_num() + 1
  !$omp do private(x, dist_to_C_sq)
  do  r = 1,solvent% Nmax
     if (solvent% species(r) == 1) then
       x = rel_pos(colloids% pos(:,1),solvent% pos(:,r),solvent_cells% edges) 
       dist_to_C_sq = dot_product(x, x)
       if (dist_to_C_sq < solvent_colloid_lj%cut_sq(1,1)) then
         if (threefry_double(state(thread_id)) <= prob) then
           solvent% flag(r) = 1 
         end if
       end if
    end if 
  end do
  !$omp end do
  !$omp end parallel
  
  end subroutine flag_particles
  
  
  subroutine change_species
    double precision :: dist_to_C_sq
    double precision :: dist_to_N_sq
    integer :: m
    double precision :: x(3)

    !$omp parallel do private(x, dist_to_C_sq, dist_to_N_sq)
    do m = 1, solvent% Nmax
       if (solvent% flag(m) == 1) then
          x = rel_pos(colloids% pos(:,1), solvent% pos(:,m), solvent_cells% edges)
          dist_to_C_sq = dot_product(x, x)
          x = rel_pos(colloids% pos(:,2), solvent% pos(:,m), solvent_cells% edges)
          dist_to_N_sq = dot_product(x, x)
          if ( &
               (dist_to_C_sq > solvent_colloid_lj%cut_sq(1, 1)) &
               .and. &
               (dist_to_N_sq > solvent_colloid_lj%cut_sq(2, 1)) &
               ) &
               then
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

  subroutine concentration_field
    double precision :: dimer_orient(3),x(3),y(3),z(3)
    double precision :: solvent_pos(3,solvent% Nmax)
    double precision :: dz,r,theta,x_pos,y_pos,z_pos
    integer :: o
    integer :: check

    dz = solvent_cells% edges(3)/400.d0
    dimer_orient = colloids% pos(:,2) - colloids% pos(:,1)
    z = dimer_orient/sqrt(dot_product(dimer_orient,dimer_orient))
    x = (/0.d0, 1.d0, -dimer_orient(2)/dimer_orient(3)/)
    x = x/sqrt(dot_product(x,x))
    y = (/z(2)*x(3)-z(3)*x(2),z(3)*x(1)-z(1)*x(3),z(1)*x(2)-z(2)*x(1)/)
    conc_z = 0

    !$omp parallel do private(x_pos, y_pos, z_pos)
    do o = 1, solvent% Nmax
       solvent_pos(:,o) = solvent% pos(:,o) - colloids% pos(:,1)
       x_pos = dot_product(x,solvent_pos(:,o))
       y_pos = dot_product(y, solvent_pos(:,o))
       z_pos = dot_product(z, solvent_pos(:,o))
       solvent_pos(:,o) = (/x_pos,y_pos,z_pos/)
    end do
    !$omp parallel do private(r, theta, check)
    do o = 1, solvent% Nmax
       r = sqrt(solvent_pos(1,o)**2 + solvent_pos(2,o)**2)
       theta = atan(solvent_pos(2,o)/solvent_pos(1,o))
       solvent_pos(1,o) = r
       solvent_pos(2,o) = theta
       solvent_pos(3,o) = solvent_pos(3,o)+colloids% pos(3,1)
       if (solvent% species(o)==2) then
          check = floor(solvent_pos(3,o)/dz)
          if ( (check>0) .and. (check<=400) ) then
             !$omp atomic
             conc_z(check) = conc_z(check) + 1
          end if
       end if 
    end do

    colloid_pos(:,1) = 0
    colloid_pos(3,1) = colloids% pos(3,1)
    colloid_pos(:,2) = 0
    colloid_pos(3,2) = d + colloids% pos(3,1)
  end subroutine concentration_field

end program setup_single_dimer
