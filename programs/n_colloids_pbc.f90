program n_colloids_pbc
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
  use md
  use ParseText
  use iso_c_binding
  use omp_lib
  implicit none

  type(PTo) :: config

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer :: rho
  double precision :: sigma, sigma_cut, epsilon
  integer :: error

  double precision :: mass
  double precision :: so_max, co_max

  double precision :: e1, e2
  double precision :: tau, dt, T, T_init, T_final
  integer :: N_MD_steps, N_loop
  integer :: n_extra_sorting
  integer :: n_threads

  type(threefry_rng_t), allocatable :: state(:)

  double precision :: skin
  double precision :: rsq
  logical :: tooclose

  integer :: i, L(3), N, N_colloids
  integer :: j, k
  type(args_t) :: args

  args = get_input_args()
  call PTparse(config, args%input_file, 11)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, args%seed)

  call h5open_f(error)

  L = PTread_ivec(config, 'L', 3)
  rho = PTread_i(config, 'rho')
  T_init = PTread_d(config, 'T')
  T = T_init
  T_final = PTread_d(config, 'T_final')
  tau = PTread_d(config, 'tau')
  N_MD_steps = PTread_i(config, 'N_MD')
  dt = tau / N_MD_steps
  N_loop = PTread_i(config, 'N_loop')

  N_colloids = PTread_i(config, 'N_colloids')
  epsilon = PTread_d(config, 'epsilon')
  sigma = PTread_d(config, 'sigma')
  sigma_cut = sigma*2.d0**(1.d0/6.d0)

  mass = rho * sigma**3 * 4 * 3.141/3
  call colloids%init(N_colloids, 1, [mass])
  colloids%species = 1

  N = rho*L(1)*L(2)*L(3) - int(rho*4*3.142/3 * sigma**3*colloids%Nmax)
  write(*,*) N, 'colloids'

  if (N <= 0) error stop 'Not enough volume available for solvent'

  call solvent_colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  epsilon = PTread_d(config, 'epsilon_colloids')
  sigma_cut = 3*sigma

  call colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ 2*sigma ], [1, 1] ), reshape( [ 2*sigma_cut ], [1, 1] ) )

  call solvent% init(N)
  do k = 1, solvent%Nmax
     solvent% vel(1,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(2,k) = threefry_normal(state(1))*sqrt(T)
     solvent% vel(3,k) = threefry_normal(state(1))*sqrt(T)
  end do
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)

  i = 1
  place_colloids: do while (i<=N_colloids)
     colloids%pos(1,i) = threefry_double(state(1))*L(1)
     colloids%pos(2,i) = threefry_double(state(1))*L(2)
     colloids%pos(3,i) = threefry_double(state(1))*L(3)

     tooclose = .false.
     check_distance: do j = 1, i-1
        rsq = sum(rel_pos(colloids%pos(:,i), colloids%pos(:,j), solvent_cells%edges)**2)
        if (rsq < colloid_lj%sigma(1,1)**2) then
           tooclose = .true.
           exit check_distance
        end if
     end do check_distance

     if (.not. tooclose) i=i+1
  end do place_colloids

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent% sort(solvent_cells)

  sigma_cut = sigma*2.d0**(1.d0/6.d0) ! restore because colloid sigma_cut is larger
  skin = 0.8
  call neigh% init(colloids% Nmax, int(200*(sigma+skin)**2))
  call neigh% make_stencil(solvent_cells, sigma_cut+skin)
  call neigh% update_list(colloids, solvent, sigma_cut+skin, solvent_cells)

  n_extra_sorting = 0

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  write(*,*) ''
  write(*,*) '    e co so     |   e co co     |   kin co      |   kin so      |   total       |   temp        |'
  write(*,*) ''

  solvent% pos_old = solvent% pos
  colloids% pos_old = colloids% pos
  do i = 1, N_loop
     md_loop: do j = 1, N_MD_steps
        call md_pos(solvent, dt)
        do k = 1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + dt**2 * colloids% force(:,k) / (2 * mass)
        end do
        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.1) .or. (so_max >= skin*0.9) ) then
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, sigma_cut+skin, solvent_cells)
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

        call md_vel(solvent, dt)
        colloids% vel = colloids% vel + dt * ( colloids% force + colloids% force_old ) / (2 * mass)

     end do md_loop

     call apply_pbc(solvent, solvent_cells% edges)
     call apply_pbc(colloids, solvent_cells% edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, sigma_cut + skin, solvent_cells)
     solvent% pos_old = solvent% pos
     colloids% pos_old = colloids% pos

     T = T_init + (i-1)*(T_final-T_init)/(N_loop-1)
     call mpcd_at_step(solvent, solvent_cells, state, T)

     if (modulo(i,100)==0) &
     write(15,*) colloids% pos + colloids% image * spread(solvent_cells% edges, dim=2, ncopies=colloids% Nmax)
     if (modulo(i,10)==0) &
     write(*,'(6f16.3)') e1, e2, mass*sum(colloids% vel**2)/2, sum(solvent% vel**2)/2, &
         e1+e2+mass*sum(colloids% vel**2)/2+sum(solvent% vel**2)/2, &
         compute_temperature(solvent, solvent_cells)

  end do

  write(*,*) 'n extra sorting', n_extra_sorting

  call h5close_f(error)

end program n_colloids_pbc
