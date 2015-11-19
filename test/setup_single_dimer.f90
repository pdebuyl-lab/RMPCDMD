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
  use iso_c_binding
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer, parameter :: rho = 9
  integer :: N
  integer :: error
  double precision, allocatable :: pos_rattle(:,:)

  double precision :: sigma, sigma_cut, epsilon1
  double precision, allocatable :: epsilon(:,:)
  double precision :: mass
  double precision :: so_max, co_max

  double precision :: e1, e2
  double precision :: tau, dt
  integer :: N_MD_steps

  type(mt19937ar_t), target :: mt

  integer :: i, L(3), seed_size, clock
  integer :: jump(3)
  integer :: j, k
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call system_clock(count=clock)
  call init_genrand(mt, int(clock, c_long))

  call h5open_f(error)

  L = [32, 32, 32]
  N = rho *L(1)*L(2)*L(3)
  
  allocate(epsilon(2,2))

  epsilon = reshape((/ 1.d0, 1.d0, 1.d0, 0.5d0 /),shape(epsilon))
  sigma = 2.d0
  sigma_cut = sigma*2**(1.d0/6.d0)

  call solvent_colloid_lj% init( epsilon, &
       reshape( [ sigma, sigma, sigma, sigma ], [2, 2] ), reshape( [ sigma_cut,sigma_cut,sigma_cut,sigma_cut ], [2, 2] ) )
  
  epsilon1 = 1.d0
  sigma = 2.d0
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init( reshape( [ epsilon1,epsilon1,epsilon1,epsilon1 ], [2, 2] ), &
       reshape( [ sigma,sigma,sigma,sigma ], [2, 2] ), reshape( [ sigma_cut,sigma_cut,sigma_cut,sigma_cut ], [2, 2] ) )

  mass = rho * sigma**3 * 4 * 3.141/3
  write(*,*) 'mass =', mass

  call solvent% init(N,2) !there will be 2 species of solvent particles

  call colloids% init(2,2) !there will be 2 species of colloids
  allocate(pos_rattle(3, colloids% Nmax))
  
  open(15,file ='dimerdata_chem01.txt')
  
  write(*, *) colloids% pos
  colloids% species(1) = 1
  colloids% species(2) = 2
  colloids% vel = 0
  
  call random_number(solvent% vel(:, :))
  solvent% vel = (solvent% vel - 0.5d0)*sqrt(6.d0*2)
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)
  colloids% pos(:,1) = solvent_cells% edges/2.0
  colloids% pos(:,2) = solvent_cells% edges/2.0 
  colloids% pos(1,2) = colloids% pos(1,2) + 2.0d0*sigma + 0.5d0
  
  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, int(300*sigma**3))
  call neigh% make_stencil(solvent_cells, 1.5d0*sigma)

  call neigh% update_list(colloids, solvent, 1.5d0*sigma, solvent_cells)

  tau = 0.1d0 !was 0.1d0 but we changed it to Kaprals' value 
  N_MD_steps = 100
  dt = tau / N_MD_steps

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  write(*,*) ''
  write(*,*) '    i   |    e co so     |   e co co     |   kin co      |   kin so      |   total       |   temp        |'
  write(*,*) ''

  do i = 1, 500
     so_max = 0
     co_max = 0
     md: do j = 1, N_MD_steps
        solvent% pos_old = solvent% pos + dt * solvent% vel + dt**2 * solvent% force / 2
        !$omp parallel do private(k)
        do k = 1, solvent% Nmax
           solvent% pos(:,k) = modulo(solvent% pos_old(:,k), solvent_cells% edges)
        end do
        
        pos_rattle = colloids% pos
        colloids% pos_old = colloids% pos + dt * colloids% vel + dt**2 * colloids% force / (2 * mass)
        
        call rattle_pos
        
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
        
        call rattle_vel
        
        call flag_particles
        
        call change_species

     end do md

     write(15,*) colloids% pos + colloids% image * spread(solvent_cells% edges, dim=2, ncopies=colloids% Nmax), &
                 colloids% vel

     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, 1.5d0*sigma, solvent_cells)

     call simple_mpcd_step(solvent, solvent_cells, mt)

     write(*,'(7f16.3)') real(i),e1, e2, mass*sum(colloids% vel**2)/2, sum(solvent% vel**2)/2, &
         e1+e2+mass*sum(colloids% vel**2)/2+sum(solvent% vel**2)/2, &
         compute_temperature(solvent, solvent_cells)

  end do
  
  call h5close_f(error)
  
contains
  
  subroutine rattle_pos
  double precision :: g !correction factor 
  double precision :: s(3) !direction vector 
  double precision :: r(3) !old direction vector
  double precision :: d
  
  d  = 2.d0*sigma + 0.5d0 !distance between the colloids in the dimer

  s = colloids% pos_old(:,1) - colloids% pos_old(:,2)
  r = pos_rattle(:,1) - pos_rattle(:,2)
  g = (dot_product(s,s) - d**2)/(2.d0*dt*(dot_product(s,r)*(2.d0/mass)))
  
  colloids% pos_old(:,1) = colloids% pos_old(:,1) - dt*g*r/mass
  colloids% pos_old(:,2) = colloids% pos_old(:,2) + dt*g*r/mass
  
  end subroutine rattle_pos


  subroutine rattle_vel
  double precision :: g !correction factor 
  double precision :: k !second correction factor
  double precision :: r_new(3) !new direction vector
  double precision :: s(3) !direction vector 
  double precision :: r(3) !old direction vector
  double precision :: d

  d = 2.d0*sigma + 0.5d0 !distance between the colloids in the dimer
  r_new = colloids% pos(:,1) - colloids% pos(:,2)
  s = colloids% pos_old(:,1) - colloids% pos_old(:,2)
  r = pos_rattle(:,1) - pos_rattle(:,2)
  g = (dot_product(s,s) - d**2)/(2.d0*dt*(dot_product(s,r)*(2.d0/mass)))
  
  colloids% vel(:,1) = colloids% vel(:,1) - dt*g*r/mass
  colloids% vel(:,2) = colloids% vel(:,2) + dt*g*r/mass
  
  k = dot_product(r_new,colloids% vel(:,1) - colloids% vel(:,2))/(d**2*(2.d0/mass))
  
  colloids% vel(:,1) = colloids% vel(:,1) - k*r_new/mass
  colloids% vel(:,2) = colloids% vel(:,2) + k*r_new/mass
  
  end subroutine rattle_vel
  
  
  subroutine flag_particles
  double precision :: dist_to_C_sq
  real :: prob 
  real :: rndnumbers(solvent% Nmax)
  integer :: r
  prob = 1.0
  
  call random_number(prob)
  
  do  r = 1,solvent% Nmax
     if (solvent% species(r) == 1) then
       dist_to_C_sq = dot_product(colloids% pos(:,1) - solvent% pos(:,r),colloids% pos(:,1) - solvent% pos(:,r))
       if (dist_to_C_sq < (sigma_cut)**2) then
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
  
  do m = 1, solvent% Nmax
     if (solvent% flag(m) == 1) then
       dist_to_C_sq = dot_product(colloids% pos(:,1) - solvent% pos(:,m),colloids% pos(:,1) - solvent% pos(:,m))
       dist_to_N_sq = dot_product(colloids% pos(:,2) - solvent% pos(:,m),colloids% pos(:,2) - solvent% pos(:,m))
       if ((dist_to_C_sq > (sigma_cut)**2) .and. (dist_to_N_sq > (sigma_cut)**2)) then
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
  integer :: n
  
  far = 25.d0
  
  do n = 1,solvent% Nmax
    if (solvent% species(i) == 2) then
      dist_to_C_sq = dot_product(colloids% pos(:,1) - solvent% pos(:,n),colloids% pos(:,1) - solvent% pos(:,n))
      dist_to_N_sq = dot_product(colloids% pos(:,2) - solvent% pos(:,n),colloids% pos(:,2) - solvent% pos(:,n))
      if ((dist_to_C_sq > far) .and. (dist_to_N_sq > far)) then
        solvent% species(n) = 1
      end if
    end if 
  end do
  
  end subroutine refuel

end program setup_single_dimer
