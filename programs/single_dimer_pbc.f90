program setup_single_dimer
  use common
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use hdf5
  use h5md_module
  use interaction
  use mt19937ar_module
  use mpcd
  use md
  use ParseText
  use iso_c_binding
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
  double precision :: kin_co

  type(mt19937ar_t), target :: mt
  type(PTo) :: config

  integer :: i, L(3), seed_size, clock
  integer :: j, k
  integer, allocatable :: seed(:)

  call PTparse(config,get_input_filename(),11)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call system_clock(count=clock)
  call init_genrand(mt, int(clock, c_long))

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
  write(*,*) 'mass =', mass

  call solvent% init(N,2) !there will be 2 species of solvent particles

  call colloids% init(2,2, mass) !there will be 2 species of colloids

  call PTkill(config)
  
  open(15,file ='dimerdata_chemKapTEST.txt')
  
  write(*, *) colloids% pos
  colloids% species(1) = 1
  colloids% species(2) = 2
  colloids% vel = 0
  
  call random_number(solvent% vel(:, :))
  solvent% vel = (solvent% vel - 0.5d0)*sqrt(12*T)
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)
  colloids% pos(:,1) = solvent_cells% edges/2.0
  colloids% pos(:,2) = solvent_cells% edges/2.0 
  colloids% pos(1,2) = colloids% pos(1,2) + d
  
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

  write(*,*) ''
  write(*,*) '    i           |    e co so     |   e co co     |   kin co      |   kin so      |   total       |   temp        |'
  write(*,*) ''

  kin_co = (mass(1)*sum(colloids% vel(:,1)**2)+mass(2)*sum(colloids% vel(:,2)**2))/2
  call thermo_write

  do i = 1, N_loop
     md: do j = 1, N_MD_steps
        call md_pos(solvent, dt)

        ! Extra copy for rattle
        colloids% pos_rattle = colloids% pos
        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(k))
        end do

        call rattle_dimer_pos(colloids, d, dt, solvent_cells% edges)

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

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)

        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)

        call md_vel(solvent, solvent_cells% edges, dt)

        do k=1, colloids% Nmax
           colloids% vel(:,k) = colloids% vel(:,k) + &
             dt * ( colloids% force(:,k) + colloids% force_old(:,k) ) / (2 * colloids% mass(k))
        end do

        call rattle_dimer_vel(colloids, d, dt, solvent_cells% edges)

        call flag_particles
        call change_species


     end do md


     write(15,*) colloids% pos + colloids% image * spread(solvent_cells% edges, dim=2, ncopies=colloids% Nmax), &
                 colloids% vel, e1+e2+(colloids% mass(1)*sum(colloids% vel(:,1)**2) &
                 +colloids% mass(2)*sum(colloids% vel(:,2)**2))/2 &
                 +sum(solvent% vel**2)/2
     
     solvent_cells% origin(1) = genrand_real1(mt) - 1
     solvent_cells% origin(2) = genrand_real1(mt) - 1
     solvent_cells% origin(3) = genrand_real1(mt) - 1

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, max_cut+skin, solvent_cells)

     call simple_mpcd_step(solvent, solvent_cells, mt)

     kin_co = (colloids% mass(1)*sum(colloids% vel(:,1)**2)+ colloids% mass(2)*sum(colloids% vel(:,2)**2))/2
     write(*,'(1i16,6f16.3,1e16.8)') i, e1, e2, &
          kin_co, sum(solvent% vel**2)/2, &
          e1+e2+kin_co+sum(solvent% vel**2)/2, &
          compute_temperature(solvent, solvent_cells), &
          sqrt(dot_product(colloids% pos(:,1) - colloids% pos(:,2),colloids% pos(:,1) - colloids% pos(:,2))) - d


  end do

  write(*,*) 'n extra sorting', n_extra_sorting
  
  call h5close_f(error)
  
contains

  subroutine thermo_write
    write(*,'(1i16,6f16.3,1e16.8)') i, e1, e2, &
         kin_co, sum(solvent% vel**2)/2, &
         e1+e2+kin_co+sum(solvent% vel**2)/2, &
         compute_temperature(solvent, solvent_cells), &
         sqrt(dot_product(colloids% pos(:,1) - colloids% pos(:,2),colloids% pos(:,1) - colloids% pos(:,2))) - d
  end subroutine thermo_write


  subroutine flag_particles
  double precision :: dist_to_C_sq
  real :: rndnumbers(solvent% Nmax)
  integer :: r
  double precision :: x(3)
  
  call random_number(rndnumbers)
  
  do  r = 1,solvent% Nmax
     if (solvent% species(r) == 1) then
       x = rel_pos(colloids% pos(:,1),solvent% pos(:,r),solvent_cells% edges) 
       dist_to_C_sq = dot_product(x, x)
       if (dist_to_C_sq < solvent_colloid_lj%cut_sq(1,1)) then
         if (rndnumbers(r) <= prob) then
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

    far = (L(1)/2)**2


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

end program setup_single_dimer
