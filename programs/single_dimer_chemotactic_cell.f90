program setup_single_dimer
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

  type(threefry_rng_t), allocatable :: state(:)
  
  integer, parameter :: N_species = 3

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj
  type(lj_params_t) :: walls_colloid_lj

  type(profile_t) :: vx

  integer :: rho
  integer :: N
  integer :: error

  double precision :: sigma_N, sigma_C, max_cut
  double precision :: epsilon(3,2), shift
  double precision :: sigma(3,2), sigma_cut(3,2)
  double precision :: mass(2)

  double precision :: v_com(3), wall_v(3,2), wall_t(2)

  double precision :: e1, e2, e_wall
  double precision :: tau, dt , T
  double precision :: d,prob
  double precision :: skin, co_max, so_max
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting
  double precision :: kin_e, temperature

  double precision :: conc_z(400)
  double precision :: colloid_pos(3,2)
  type(h5md_file_t) :: hfile
  type(h5md_element_t) :: dummy_element
  integer(HID_T) :: fields_group
  type(h5md_element_t) :: rho_xy_el
  type(thermo_t) :: thermo_data
  type(particle_system_io_t) :: dimer_io
  type(particle_system_io_t) :: solvent_io
  integer(HID_T) :: box_group

  type(PTo) :: config
  integer :: i, L(3),  n_threads
  integer :: j, k, m
  
  type(timer_t) :: flag_timer, change_timer, buffer_timer

  integer, allocatable :: rho_xy(:,:,:)

  double precision :: g(3) !gravity
  logical :: fixed, on_track, stopped, order
  integer :: bufferlength
  integer :: steps_fixed
  fixed = .true.
  on_track = .true.
  stopped = .false.

  call PTparse(config,get_input_filename(),11)

  call flag_timer%init('flag')
  call change_timer%init('change')
  call buffer_timer%init('buffer')

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, PTread_c_int64(config, 'seed'))

  call h5open_f(error)

  g = 0
  g(1) = PTread_d(config, 'g')
  bufferlength = PTread_i(config, 'buffer_length')
  prob = PTread_d(config,'probability')

  L = PTread_ivec(config, 'L', 3)
  L(1) = L(1) + bufferlength
  
  rho = PTread_i(config, 'rho')
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T')
  d = PTread_d(config, 'd')
  order = PTread_l(config, 'order')

  wall_v = 0
  wall_t = [T, T]
  
  tau =PTread_d(config, 'tau')
  N_MD_steps = PTread_i(config, 'N_MD')
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop')
  steps_fixed = PTread_i(config, 'steps_fixed')
  
  sigma_C = PTread_d(config, 'sigma_C')
  sigma_N = PTread_d(config, 'sigma_N')

  epsilon(:,1) = PTread_dvec(config, 'epsilon_C', N_species)
  epsilon(:,2) = PTread_dvec(config, 'epsilon_N', N_species)

  sigma(:,1) = sigma_C
  sigma(:,2) = sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)

  epsilon = 1.d0

  sigma(1,1) = 2*sigma_C
  sigma(1,2) = sigma_C + sigma_N
  sigma(2,1) = sigma_C + sigma_N
  sigma(2,2) = 2*sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init(epsilon(1:2,:), sigma(1:2,:), sigma_cut(1:2,:))

  epsilon = 1.d0
  sigma(1,:) = [sigma_C, sigma_N]
  sigma_cut = sigma*3**(1.d0/6.d0)
  shift = max(sigma_C, sigma_N)*2**(1./6.) + 0.25
  call walls_colloid_lj% init(epsilon(1:1,:), sigma(1:1,:), sigma_cut(1:1,:), shift)
  write(*,*) epsilon(1:2,:), sigma(1:2,:), sigma_cut(1:2,:), shift

  mass(1) = rho * sigma_C**3 * 4 * 3.14159265/3
  mass(2) = rho * sigma_N**3 * 4 * 3.14159265/3
  write(*,*) 'mass =', mass

  call solvent% init(N,N_species)

  call colloids% init(2,2, mass) !there will be 2 species of colloids

  call hfile%create(PTread_s(config, 'h5md_file'), 'RMPCDMD::single_dimer_chemotactic_cell', &
       'N/A', 'Pierre de Buyl')
  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call PTkill(config)

  dimer_io%force_info%store = .false.
  dimer_io%id_info%store = .false.
  dimer_io%position_info%store = .true.
  dimer_io%position_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%position_info%step = N_MD_steps
  dimer_io%position_info%time = N_MD_steps*dt
  dimer_io%image_info%store = .true.
  dimer_io%image_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%image_info%step = N_MD_steps
  dimer_io%image_info%time = N_MD_steps*dt
  dimer_io%velocity_info%store = .true.
  dimer_io%velocity_info%mode = ior(H5MD_LINEAR,H5MD_STORE_TIME)
  dimer_io%velocity_info%step = N_MD_steps
  dimer_io%velocity_info%time = N_MD_steps*dt
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

  open(17,file ='dimerdata_FullExp_1.txt')
  open(19,file ='dimerdata_vx_flow_wall.txt')  

  colloids% species(1) = 1
  colloids% species(2) = 2
  colloids% vel = 0

  do i=1, solvent% Nmax
     solvent% vel(1,i) = threefry_normal(state(1))
     solvent% vel(2,i) = threefry_normal(state(1))
     solvent% vel(3,i) = threefry_normal(state(1))
  end do
  solvent%vel = solvent%vel*sqrt(T)
  v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
  solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))

  solvent% force = 0

  do m = 1, solvent% Nmax
     if (solvent% pos(2,m) < (L(2)/2.d0)) then
        solvent% species(m) = 1
     else
        solvent% species(m) = 3
     end if
  end do

  call solvent_cells%init(L, 1.d0,has_walls = .true.)

  allocate(rho_xy(N_species, L(2), L(1)))
  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  call rho_xy_el%create_time(fields_group, 'rho_xy', rho_xy, ior(H5MD_LINEAR,H5MD_STORE_TIME), &
       step=N_loop*N_MD_steps, time=N_loop*N_MD_steps*dt)
  call h5gclose_f(fields_group, error)

  colloids% pos(3,:) = solvent_cells% edges(3)/2.d0
  if (order) then
     colloids% pos(1,1) = sigma_C*2**(1.d0/6.d0) + 1
     colloids% pos(1,2) = colloids% pos(1,1) + d
     colloids% pos(2,:) = solvent_cells% edges(2)/2.d0 + 1.5d0*maxval([sigma_C,sigma_N])
  else
     colloids% pos(1,2) = sigma_N*2**(1.d0/6.d0) + 1
     colloids% pos(1,1) = colloids% pos(1,2) + d
     colloids% pos(2,:) = solvent_cells% edges(2)/2.d0 + 1.5d0*maxval([sigma_C,sigma_N])
  end if
  call h5gcreate_f(dimer_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call h5gcreate_f(solvent_io%group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call dummy_element%create_fixed(box_group, 'edges', solvent_cells%edges)
  call h5gclose_f(box_group, error)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, 10*int(300*max(sigma_C,sigma_N)**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  call vx% init(0.d0, solvent_cells% edges(3), L(3))

  i = 0

  write(*,*) 'Running for', N_loop, 'loops'
  !start RMPCDMD
  setup: do i = 1, N_loop
     if (modulo(i,20) == 0) write(*,'(i09)',advance='no') i
     md_loop: do j = 1, N_MD_steps
        call mpcd_stream_zwall_light(solvent, solvent_cells, dt,g)

        colloids% pos_rattle = colloids% pos
        
        if (.not. fixed) then
           if (on_track) then
              !only update the flow direction
              do k=1, colloids% Nmax
                 colloids% pos(1,k) = colloids% pos(1,k) + dt * colloids% vel(1,k) + &
                      dt**2 * colloids% force(1,k) / (2 * colloids% mass(k))
              end do
           else
              do k=1, colloids% Nmax
                 colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                      dt**2 * colloids% force(:,k) / (2 * colloids% mass(k))
              end do
           end if 
           call rattle_dimer_pos(colloids, d, dt, solvent_cells% edges)
        end if  
   
        if (on_track) then
           do k=1, colloids% Nmax 
              if (colloids% pos(1,k) > bufferlength) then
                 on_track = .false.
                 write(*,*) 'on_track', on_track
              end if 
           end do
        end if
        
        if (.not. on_track) then
           do k=1, colloids% Nmax 
              if (colloids% pos(1,k) > solvent_cells% edges(1)) then
                 stopped = .true.
                 write(*,*) 'stopped', stopped
              end if 
           end do
        end if

        if (stopped) exit setup

        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin/2) .or. (so_max >= skin/2) ) then
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, max_cut + skin, solvent_cells)
           solvent% pos_old = solvent% pos
           colloids% pos_old = colloids% pos
           n_extra_sorting = n_extra_sorting + 1
        end if

        call buffer_timer%tic()
        call buffer_particles(solvent,solvent_cells% edges(:), bufferlength)
        call buffer_timer%tac()

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)

        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
        if (.not. on_track) then
           e_wall = lj93_zwall(colloids, solvent_cells% edges, walls_colloid_lj)
        end if 
        if (on_track) then
           colloids% force(2,:) = 0
           colloids% force(3,:) = 0
           if (fixed) then
              colloids% force(1,:) = 0
           end if 
        end if 

        call md_vel_flow_partial(solvent, dt, g)
        if (.not. fixed) then
           if (on_track) then
              !only update in the direction of the flow
              do k=1, colloids% Nmax
                 colloids% vel(1,k) = colloids% vel(1,k) + &
                   dt * ( colloids% force(1,k) + colloids% force_old(1,k) ) / (2 * colloids% mass(k))
              end do
           else
              do k=1, colloids% Nmax
                 colloids% vel(:,k) = colloids% vel(:,k) + &
                   dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(k))
              end do
           end if
           call rattle_dimer_vel(colloids, d, dt, solvent_cells% edges)
        end if 
        if (.not.fixed) then
           call flag_timer%tic()
           call flag_particles
           call flag_timer%tac()
           call change_timer%tic()
           call change_species
           call change_timer%tac()
        end if

     end do md_loop

     

     write(17,*) colloids% pos + colloids% image * spread(solvent_cells% edges, dim=2, ncopies=colloids% Nmax), &
                 colloids% vel, e1+e2+e_wall+(colloids% mass(1)*sum(colloids% vel(:,1)**2) &
                 +colloids% mass(2)*sum(colloids% vel(:,2)**2))/2 &
                 +sum(solvent% vel**2)/2
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1

     call compute_vx(solvent, vx)
     if (modulo(i, 50) == 0) then
        call vx% norm()
        write(19,*) vx% data
        flush(19)
        call vx% reset()
     end if

     call compute_rho_xy
     call rho_xy_el%append(rho_xy)

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[10, 10], bulk_temperature = T)
     
     temperature = compute_temperature(solvent, solvent_cells)
     kin_e = (colloids% mass(1)*sum(colloids% vel(:,1)**2) + &
          colloids% mass(2)*sum(colloids% vel(:,2)**2))/2 + &
          sum(solvent% vel**2)/2
     call thermo_data%append(hfile, temperature, e1+e2+e_wall, kin_e, e1+e2+e_wall+kin_e)

     call dimer_io%position%append(colloids%pos)
     call dimer_io%velocity%append(colloids%vel)
     call dimer_io%image%append(colloids%image)

     if (i >= steps_fixed) then
        fixed = .false.
     end if
     
  end do setup

  call thermo_data%append(hfile, temperature, e1+e2+e_wall, kin_e, e1+e2+e_wall+kin_e, add=.false., force=.true.)

  write(*,*) 'n extra sorting', n_extra_sorting

  call solvent_io%position%append(solvent%pos)
  call solvent_io%velocity%append(solvent%vel)
  call solvent_io%image%append(solvent%image)
  call solvent_io%species%append(solvent%species)

  call rho_xy_el%close()
  call dimer_io%close()
  call hfile%close()
  call h5close_f(error)

  write(*,'(a16,f15.3)') solvent%time_stream%name, solvent%time_stream%total
  write(*,'(a16,f15.3)') solvent%time_md_vel%name, solvent%time_md_vel%total
  write(*,'(a16,f15.3)') solvent%time_step%name, solvent%time_step%total
  write(*,'(a16,f15.3)') solvent%time_count%name, solvent%time_count%total
  write(*,'(a16,f15.3)') solvent%time_sort%name, solvent%time_sort%total
  write(*,'(a16,f15.3)') solvent%time_ct%name, solvent%time_ct%total
  write(*,'(a16,f15.3)') solvent%time_max_disp%name, solvent%time_max_disp%total
  write(*,'(a16,f15.3)') flag_timer%name, flag_timer%total
  write(*,'(a16,f15.3)') change_timer%name, buffer_timer%total
  write(*,'(a16,f15.3)') buffer_timer%name, buffer_timer%total
  write(*,'(a16,f15.3)') neigh%time_update%name, neigh%time_update%total
  write(*,'(a16,f15.3)') neigh%time_force%name, neigh%time_force%total

  write(*,'(a16,f15.3)') 'total                          ', &
       solvent%time_stream%total + solvent%time_step%total + solvent%time_count%total +&
       solvent%time_sort%total + solvent%time_ct%total + solvent%time_md_vel%total + solvent%time_max_disp%total + &
       flag_timer%total + change_timer%total + buffer_timer%total + neigh%time_update%total + neigh%time_force%total
  
contains

  subroutine flag_particles
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
               (dist_to_C_sq > solvent_colloid_lj%cut_sq(1,1)) &
               .and. &
               (dist_to_N_sq > solvent_colloid_lj%cut_sq(1,2)) &
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

    do o = 1, solvent% Nmax
       solvent_pos(:,o) = solvent% pos(:,o) - colloids% pos(:,1)
       x_pos = dot_product(x,solvent_pos(:,o))
       y_pos = dot_product(y, solvent_pos(:,o))
       z_pos = dot_product(z, solvent_pos(:,o))
       solvent_pos(:,o) = (/x_pos,y_pos,z_pos/)
    end do
    do o = 1, solvent% Nmax
       r = sqrt(solvent_pos(1,o)**2 + solvent_pos(2,o)**2)
       theta = atan(solvent_pos(2,o)/solvent_pos(1,o))
       solvent_pos(1,o) = r
       solvent_pos(2,o) = theta
       solvent_pos(3,o) = solvent_pos(3,o)+colloids% pos(3,1)
       if (solvent% species(o)==2) then
          check = floor(solvent_pos(3,o)/dz)
          conc_z(check) = conc_z(check) + 1
       end if 
    end do
    colloid_pos(:,1) = 0
    colloid_pos(3,1) = colloids% pos(3,1)
    colloid_pos(:,2) = 0
    colloid_pos(3,2) = d + colloids% pos(3,1)
  end subroutine concentration_field

  subroutine mpcd_stream_zwall_light(particles, cells, dt,g)
    type(particle_system_t), intent(inout) :: particles
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: dt
    double precision, dimension(3), intent(in):: g

    integer :: i
    double precision :: pos_min(3), pos_max(3)
    double precision, dimension(3) :: old_pos, old_vel
    double precision :: t_c

    pos_min = 0
    pos_max = cells% edges

    call solvent%time_stream%tic()
    do i = 1, particles% Nmax
       old_pos = particles% pos(:,i) 
       old_vel = particles% vel(:,i)
       particles% pos(:,i) = particles% pos(:,i) + dt * particles% vel(:,i) + dt**2 * (particles% force(:,i) + g)/ 2
       !particles% pos(2,i) = modulo( particles% pos(2,i) , cells% edges(2) )
       !particles% pos(1,i) = modulo( particles% pos(1,i) , cells% edges(1) )
       if (cells% has_walls) then
          if (particles% pos(3,i) < pos_min(3)) then
             t_c = abs(old_pos(3)/old_vel(3))
             particles% vel(:,i) = -(old_vel + g*t_c) + g*(dt-t_c)
             particles% pos(:,i) = old_pos + old_vel*t_c + g*t_c**2/2 - (old_vel + g*t_c)*(dt-t_c)+(dt-t_c)**2*g/2
             particles% wall_flag(i) = 1
          else if (particles% pos(3,i) > pos_max(3)) then
             t_c = abs((pos_max(3)-old_pos(3))/old_vel(3))
             particles% vel(:,i) = -(old_vel + g*t_c) + g*(dt-t_c)
             particles% pos(:,i) = old_pos + old_vel*t_c + g*t_c**2/2 - (old_vel + g*t_c)*(dt-t_c)+(dt-t_c)**2*g/2
             particles% wall_flag(i) = 1
          end if
       !else
          !particles% pos(3,i) = modulo( particles% pos(3,i) , cells% edges(3) )
       end if 
    end do
    call solvent%time_stream%tac()
  end subroutine mpcd_stream_zwall_light
  
  subroutine md_vel_flow_partial(particles, dt, ext_force)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt
    double precision, intent(in) :: ext_force(3)

    integer :: k

    call solvent%time_md_vel%tic()
    !$omp parallel do
    do k = 1, particles% Nmax
       if (particles% wall_flag(k) == 0) then
          particles% vel(:,k) = particles% vel(:,k) + &
               dt * ( particles% force(:,k) + particles% force_old(:,k) ) / 2 &
               + dt*ext_force
       else
         particles% wall_flag(k) = 0
       end if
    end do
    call solvent%time_md_vel%tac()

  end subroutine md_vel_flow_partial

  subroutine buffer_particles(particles,edges, bufferlength)
     type(particle_system_t), intent(inout) :: particles
     double precision, intent(in) :: edges(3)
     integer, intent(in) :: bufferlength
  
     integer :: k  

     do k = 1, particles% Nmax
        if (particles% pos(1,k) < bufferlength) then
           if (particles% pos(2,k) < edges(2)/2.d0) then
              particles% species(k) = 1
           else
              particles% species(k) = 3
           end if  
        end if 
     end do
  end subroutine buffer_particles

  subroutine compute_rho_xy
    integer :: i, s, ix, iy

    rho_xy = 0
    do i = 1, solvent%Nmax
       s = solvent%species(i)
       if (s <= 0) continue
       ix = modulo(floor(solvent%pos(1,i)/solvent_cells%a), L(1)) + 1
       iy = modulo(floor(solvent%pos(2,i)/solvent_cells%a), L(2)) + 1
       rho_xy(s, iy, ix) = rho_xy(s, iy, ix) + 1
    end do

  end subroutine compute_rho_xy

end program setup_single_dimer
