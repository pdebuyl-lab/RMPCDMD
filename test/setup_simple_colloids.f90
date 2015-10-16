program setup_simple_colloids
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

  integer, parameter :: N = 3000
  integer :: error

  double precision :: epsilon, sigma, sigma_cut

  double precision :: e1, e2, kin
  double precision :: tau, dt
  integer :: N_MD_steps

  type(mt19937ar_t), target :: mt

  integer :: i, L(3), seed_size, clock
  integer :: j
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

  L = [12, 12, 12]

  epsilon = 1
  sigma = 1
  sigma_cut = sigma*2**(1.d0/6.d0)

  call solvent_colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  epsilon = 1
  sigma = 2
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  call solvent% init(N)

  call colloids% init_from_file('input_data.h5', 'colloids', H5MD_LINEAR, 4)
  write(*, *) colloids% pos
  colloids% species = 1

  call random_number(solvent% vel(:, :))
  solvent% vel(:, :) = solvent% vel(:, :) - 0.5d0
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj)

  call solvent_cells%count_particles(solvent% pos)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, 300)
  call neigh% make_stencil(solvent_cells, 1.5d0)

  call neigh% update_list(colloids, solvent, 1.5d0, solvent_cells)

  tau = 0.1d0
  N_MD_steps = 100
  dt = tau / N_MD_steps

  do i = 1, 10
     md: do j = 1, N_MD_steps
        solvent% pos = solvent% pos + dt * solvent% vel + dt**2 * solvent% force / 2
        solvent% pos = modulo(solvent% pos, spread(solvent_cells% edges, 2, colloids% Nmax))
        colloids% pos = colloids% pos + dt * colloids% vel + dt**2 * colloids% force / 2
        colloids% pos = modulo(colloids% pos, spread(solvent_cells% edges, 2, colloids% Nmax))

        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        !e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)

        solvent% vel = solvent% vel + dt * ( solvent% force + solvent% force_old ) / 2
        colloids% vel = colloids% vel + dt * ( colloids% force + colloids% force_old ) / 2

        call neigh% update_list(colloids, solvent, 1.5d0, solvent_cells)

     end do md

     call solvent% sort(solvent_cells)
     call simple_mpcd_step(solvent, solvent_cells, mt)

     write(*,*) 'co ', colloids% pos(:,1), colloids% force(:,1)

     !write(*,*) 'so ', solvent% pos(:,1), solvent% force(:,1)

     !write(*,*) compute_temperature(solvent, solvent_cells)

  end do


  call h5close_f(error)

end program setup_simple_colloids
