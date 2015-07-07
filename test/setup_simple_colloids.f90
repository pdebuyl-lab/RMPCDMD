program setup_simple_colloids
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use hdf5
  use interaction
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

  integer :: i, L(3), seed_size, clock
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

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

  call colloids% init_from_file('input_data.h5', 'colloids')

  call random_number(solvent% vel(:, :))
  solvent% vel(:, :) = solvent% vel(:, :) - 0.5d0
  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0, colloids, solvent_colloid_lj)

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(solvent% pos)

  print *, sum(solvent_cells%cell_count)
  print *, solvent_cells%cell_count
  print *, solvent_cells%cell_start(1), solvent_cells%cell_start(solvent_cells%N)

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, 300)

  call neigh% update_list(colloids, solvent, 1.4d0, solvent_cells)

  open(12, file='neigh_list')
  do i = 1, colloids% Nmax
     write(12, *) neigh% n(i)
     write(12, *) neigh% list(:, i)
  end do

  open(12, file='sorted_pos')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  open(12, file='colloids_pos')
  do i=1, colloids% Nmax
     write(12, *) colloids% pos(:,i)
  end do
  close(12)

  do i = 1, 1000

     call solvent% random_placement(L*1.d0)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, 1.4d0, solvent_cells)

  end do

  open(12, file='sorted_pos_2')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  call h5close_f(error)

end program setup_simple_colloids
