program setup_simple_colloid
  use common
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use hdf5
  use h5md_module
  use interaction
  use threefry_module
  use mpcd
  use iso_c_binding
  use omp_lib
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer, parameter :: rho = 10
  integer :: N
  integer :: error

  double precision :: epsilon, sigma, sigma_cut
  double precision :: mass
  double precision :: so_max, co_max

  double precision :: e1, e2
  double precision :: tau, dt
  integer :: N_MD_steps

  integer :: i, L(3)
  integer :: jump(3)
  integer :: j, k
  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))

  do i = 1, n_threads
     state(i)%counter%c0 = 0
     state(i)%counter%c1 = 0
     state(i)%key%c0 = 0
     state(i)%key%c1 = 719287321987291_c_long
  end do

  call h5open_f(error)

  L = [20, 20, 20]
  N = rho *L(1)*L(2)*L(3)

  epsilon = 1
  sigma = 1
  sigma_cut = sigma*2**(1.d0/6.d0)

  call solvent_colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  epsilon = 1
  sigma = 1
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  mass = rho * sigma**3 * 4 * 3.141/3
  write(*,*) 'mass =', mass

  call solvent% init(N)

  call colloids% init(1)
  
  open(15,file ='data56.txt')
  
  write(*, *) colloids% pos
  colloids% species = 1
  colloids% vel = 0
  
  call random_number(solvent% vel(:, :))
  solvent% vel = (solvent% vel - 0.5d0)*sqrt(6.d0*2)
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)
  colloids% pos(:,1) = solvent_cells% edges/2.0
  
  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, 300)
  call neigh% make_stencil(solvent_cells, 1.5d0)

  call neigh% update_list(colloids, solvent, 1.5d0, solvent_cells)

  tau = 0.1d0
  N_MD_steps = 100
  dt = tau / N_MD_steps

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  do i = 1, 100
     md_equil: do j = 1, N_MD_steps
        solvent% pos_old = solvent% pos + dt * solvent% vel + dt**2 * solvent% force / 2
        !$omp parallel do private(k)
        do k = 1, solvent% Nmax
           solvent% pos(:,k) = modulo(solvent% pos_old(:,k), solvent_cells% edges)
        end do
        colloids% pos_old = colloids% pos + dt * colloids% vel + dt**2 * colloids% force / (2 * mass)
        do k = 1, colloids% Nmax
           jump = floor(colloids% pos_old(:,k) / solvent_cells% edges)
           colloids% image(:,k) = colloids% image(:,k) + jump
           colloids% pos(:,k) = colloids% pos_old(:,k) - jump*solvent_cells% edges
        end do

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)
        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)

        !$omp parallel do private(k)
        do k = 1, solvent% Nmax
           solvent% vel(:,k) = solvent% vel(:,k) + dt * ( solvent% force(:,k) + solvent% force_old(:,k) ) / 2
        end do
        colloids% vel = colloids% vel + dt * ( colloids% force + colloids% force_old ) / (2 * mass)

     end do md_equil

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, 1.5d0, solvent_cells)
     call simple_mpcd_step(solvent, solvent_cells, state)

  end do

  write(*,*) ''
  write(*,*) '    e co so     |   e co co     |   kin co      |   kin so      |   total       |   temp        |'
  write(*,*) ''

  do i = 1, 750
     so_max = 0
     co_max = 0
     md: do j = 1, N_MD_steps
        solvent% pos_old = solvent% pos + dt * solvent% vel + dt**2 * solvent% force / 2
        !$omp parallel do private(k)
        do k = 1, solvent% Nmax
           solvent% pos(:,k) = modulo(solvent% pos_old(:,k), solvent_cells% edges)
        end do
        colloids% pos_old = colloids% pos + dt * colloids% vel + dt**2 * colloids% force / (2 * mass)
        do k = 1, colloids% Nmax
           jump = floor(colloids% pos_old(:,k) / solvent_cells% edges)
           colloids% image(:,k) = colloids% image(:,k) + jump
           colloids% pos(:,k) = colloids% pos_old(:,k) - jump*solvent_cells% edges
        end do
        so_max = max(maxval(sqrt(sum((solvent% pos - solvent% pos_old)**2, dim=1))), so_max)
        co_max = max(maxval(sqrt(sum((colloids% pos - colloids% pos_old)**2, dim=1))), co_max)

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)
        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)

        !$omp parallel do private(k)
        do k = 1, solvent% Nmax
           solvent% vel(:,k) = solvent% vel(:,k) + dt * ( solvent% force(:,k) + solvent% force_old(:,k) ) / 2
        end do
        colloids% vel = colloids% vel + dt * ( colloids% force + colloids% force_old ) / (2 * mass)

     end do md

     write(15,*) colloids% pos + colloids% image * spread(solvent_cells% edges, dim=2, ncopies=colloids% Nmax)

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, 1.5d0, solvent_cells)

     call simple_mpcd_step(solvent, solvent_cells, state)

     write(*,'(6f16.3)') e1, e2, mass*sum(colloids% vel**2)/2, sum(solvent% vel**2)/2, &
         e1+e2+mass*sum(colloids% vel**2)/2+sum(solvent% vel**2)/2, &
         compute_temperature(solvent, solvent_cells)

  end do


  call h5close_f(error)

end program setup_simple_colloid
