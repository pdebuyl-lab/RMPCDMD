!> Simulate a single dimer nanomotor
!!
!! Consider a dimer nanomotor in a periodic simulation cell filled with A particles. After a
!! collision with the catalytic sphere of the dimer, a A particle is converted to B.
!!
!! \param L           length of simulation box in the 3 dimensions
!! \param rho         fluid number density
!! \param T           Temperature. Used for setting initial velocities and (if enabled) bulk thermostatting.
!! \param tau         MPCD collision time
!! \param probability probability to change A to B upon collision
!! \param bulk_rmpcd  use bulkd rmpcd reaction for B->A instead of resetting
!! \param bulk_rate   rate of B->A reaction
!! \param N_MD        number MD steps occuring in tau
!! \param N_loop      number of MPCD timesteps
!! \param colloid_sampling interval (in MD steps) of sampling the colloid position and velocity
!! \param equilibration_loops number of MPCD steps for equilibration
!! \param sigma_C     radius of C sphere
!! \param sigma_N     radius of N sphere
!! \param d           length of rigid link
!! \param epsilon_C   interaction parameter of C sphere with both solvent species (2 elements)
!! \param epsilon_N   interaction parameter of N sphere with both solvent species (2 elements)
!! \param epsilon_C_C interaction parameter among C spheres
!! \param epsilon_N_C interaction parameter among N and C spheres
!! \param epsilon_N_N interaction parameter among N spheres

program single_dimer_pbc
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

  type(PTo) :: config

  integer :: i, L(3)
  integer :: j, k
  type(timer_t), target :: varia, main, time_flag, time_refuel, time_change
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
  double precision :: bulk_rate
  logical :: bulk_rmpcd

  integer, dimension(N_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el

  integer, parameter :: block_length = 8
  type(axial_correlator_t) :: axial_cf
  integer(HID_T) :: correlator_group

  integer :: equilibration_loops
  integer :: colloid_sampling
  logical :: sampling
  type(histogram_t) :: z_hist
  type(h5md_element_t) :: z_hist_el
  double precision :: cyl_shell_rmin, cyl_shell_rmax
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

  call timer_list%init(16)
  call timer_list%append(varia)
  call timer_list%append(time_flag)
  call timer_list%append(time_refuel)
  call timer_list%append(time_change)

  call h5open_f(error)
  call hfile%create(args%output_file, 'RMPCDMD::single_dimer_pbc', &
       RMPCDMD_REVISION, 'Pierre de Buyl')
  call h5gcreate_f(hfile%id, 'parameters', params_group, error)
  call hdf5_util_write_dataset(params_group, 'seed', args%seed)
  prob = PTread_d(config,'probability', loc=params_group)
  bulk_rmpcd = PTread_l(config, 'bulk_rmpcd', loc=params_group)
  bulk_rate = PTread_d(config, 'bulk_rate', loc=params_group)

  L = PTread_ivec(config, 'L', 3, loc=params_group)
  rho = PTread_i(config, 'rho', loc=params_group)
  N = rho *L(1)*L(2)*L(3)

  T = PTread_d(config, 'T', loc=params_group)
  d = PTread_d(config, 'd', loc=params_group)
  
  tau = PTread_d(config, 'tau', loc=params_group)
  N_MD_steps = PTread_i(config, 'N_MD', loc=params_group)
  colloid_sampling = PTread_i(config, 'colloid_sampling', loc=params_group)
  if (modulo(N_MD_steps, colloid_sampling) /= 0) then
     error stop 'colloid_sampling must divide N_MD with no remainder'
  end if
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop', loc=params_group)
  equilibration_loops = PTread_i(config, 'equilibration_loops', loc=params_group)

  sigma_C = PTread_d(config, 'sigma_C', loc=params_group)
  sigma_N = PTread_d(config, 'sigma_N', loc=params_group)

  ! solvent index first, colloid index second, in solvent_colloid_lj
  epsilon(:,1) = PTread_dvec(config, 'epsilon_C', 2, loc=params_group)
  epsilon(:,2) = PTread_dvec(config, 'epsilon_N', 2, loc=params_group)

  sigma(:,1) = sigma_C
  sigma(:,2) = sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)
  max_cut = maxval(sigma_cut)

  call solvent_colloid_lj% init(epsilon, sigma, sigma_cut)

  epsilon(1,1) = PTread_d(config, 'epsilon_C_C', loc=params_group)
  epsilon(1,2) = PTread_d(config, 'epsilon_N_C', loc=params_group)
  epsilon(2,1) = epsilon(1,2)
  epsilon(2,2) = PTread_d(config, 'epsilon_N_N', loc=params_group)

  sigma(1,1) = 2*sigma_C
  sigma(1,2) = sigma_C + sigma_N
  sigma(2,1) = sigma_C + sigma_N
  sigma(2,2) = 2*sigma_N
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init(epsilon, sigma, sigma_cut)

  mass(1) = rho * sigma_C**3 * 4 * 3.14159265/3
  mass(2) = rho * sigma_N**3 * 4 * 3.14159265/3

  call solvent% init(N,2, system_name='solvent') !there will be 2 species of solvent particles

  call colloids% init(2,2, mass, system_name='colloids') !there will be 2 species of colloids

  call thermo_data%init(hfile, n_buffer=50, step=N_MD_steps, time=N_MD_steps*dt)

  call h5gclose_f(params_group, error)
  call PTkill(config)

  call axial_cf%init(block_length, N_loop, N_loop*N_MD_steps)
  
  colloids% species(1) = 1
  colloids% species(2) = 2
  colloids% vel = 0

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
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)
  colloids% pos(:,1) = solvent_cells% edges/2.0
  colloids% pos(:,2) = solvent_cells% edges/2.0 
  colloids% pos(1,2) = colloids% pos(1,2) + d

  call n_solvent_el%create_time(hfile%observables, 'n_solvent', &
       n_solvent, H5MD_LINEAR, step=N_MD_steps, &
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

  call neigh% init(colloids% Nmax, int(300*max(sigma_C,sigma_N)**3))

  skin = 1.5
  n_extra_sorting = 0

  call neigh% make_stencil(solvent_cells, max_cut+skin)

  call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

  cyl_shell_rmin = max(sigma_N, sigma_C)
  cyl_shell_rmax = cyl_shell_rmin + 2
  call h5gcreate_f(hfile%id, 'fields', fields_group, error)
  call z_hist%init(-5.d0, 5.d0+d, 100)
  call z_hist_el%create_time(fields_group, 'cylindrical_shell_histogram', z_hist%data, &
       H5MD_LINEAR, step=N_MD_steps, time=N_MD_steps*dt)
  call h5md_write_attribute(z_hist_el%id, 'r_min', cyl_shell_rmin)
  call h5md_write_attribute(z_hist_el%id, 'r_max', cyl_shell_rmax)
  call h5md_write_attribute(z_hist_el%id, 'xmin', z_hist%xmin)
  call h5md_write_attribute(z_hist_el%id, 'dx', z_hist%dx)
  call h5gclose_f(fields_group, error)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
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

        call md_vel(solvent, dt)

        call varia%tic()
        do k=1, colloids% Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(k))
        end do

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

        call time_flag%tic()
        call flag_particles_nl
        call time_flag%tac()
        call time_change%tic()
        call change_species
        call time_change%tac()

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

        call thermo_data%append(hfile, temperature, e1+e2, kin_e, e1+e2+kin_e, v_com)

        com_pos = ( colloids%pos(:,1)+colloids%image(:,1)*solvent_cells%edges + &
             colloids%pos(:,2)+colloids%image(:,2)*solvent_cells%edges)/2
        call axial_cf%add(i-equilibration_loops, com_pos, unit_r)

     end if

     call solvent_cells%random_shift(state(1))
     call varia%tac()

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)
     call varia%tic()
     !$omp parallel do
     do k = 1, solvent%Nmax
        solvent% pos_old(:,k) = solvent% pos(:,k)
     end do
     colloids% pos_old = colloids% pos
     call varia%tac()

     call simple_mpcd_step(solvent, solvent_cells, state)
     
     call time_refuel%tic()
     if (bulk_rmpcd) then
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

     z_hist%data = 0
     call compute_cylindrical_shell_histogram(z_hist, colloids%pos(:,1), colloids%pos(:,2), &
          solvent_cells%edges, 2, cyl_shell_rmin, cyl_shell_rmax, solvent)
     call z_hist_el%append(z_hist%data)

  end do
  call main%tac()
  write(*,*) ''

  write(*,*) 'n extra sorting', n_extra_sorting

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
  call z_hist_el%close()
  call n_solvent_el%close()
  call hfile%close()
  call h5close_f(error)

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
             solvent% flag(r) = 1
          end if
       end if
    end do

  end subroutine flag_particles_nl
  
  subroutine change_species
    double precision :: dist_to_C_sq
    double precision :: dist_to_N_sq
    integer :: m
    double precision :: x(3)
    integer :: thread_id

    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num() + 1
    !$omp do private(x, dist_to_C_sq, dist_to_N_sq)
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
             if (threefry_double(state(thread_id)) <= prob) then
                solvent% species(m) = 2
             end if
             solvent% flag(m) = 0
          end if
       end if
    end do
    !$omp end do
    !$omp end parallel
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

end program single_dimer_pbc
