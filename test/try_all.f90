program try_all
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use interaction
  use hdf5
  use h5md_module
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  integer(HID_T) :: colloids_group, box_group

  integer, parameter :: N = 1000
  integer, parameter :: N_colloids = 3

  double precision :: epsilon, sigma, sigma_cut

  integer :: i, L(3), seed_size, clock, error
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call h5open_f(error)

  L = [8, 5, 5]

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
  call colloids% init(N_colloids)
  colloids% species = 1

  call colloids% random_placement(L*1.d0)

  call random_number(solvent% vel(:, :))
  solvent% vel(:, :) = solvent% vel(:, :) - 0.5d0
  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0, colloids, solvent_colloid_lj)

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create('test_data.h5', 'RMPCDMD:try_all', '0.0 dev', 'Pierre de Buyl')
  call h5gcreate_f(datafile% particles, 'colloids', colloids_group, error)

  call h5gcreate_f(colloids_group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call h5md_write_attribute(box_group, 'boundary', ['periodic', 'periodic', 'periodic'])
  call elem% create_fixed(box_group, 'edges', L*1.d0)
  call h5gclose_f(box_group, error)

  call solvent% sort(solvent_cells)

  call neigh% init(N_colloids, 300)

  call neigh% update_list(colloids, solvent, 1.4d0, solvent_cells)

  open(12, file='neigh_list')
  do i = 1, N_colloids
     write(12, *) neigh% n(i)
     write(12, *) neigh% list(:, i)
  end do

  open(12, file='sorted_pos')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  open(12, file='colloids_pos')
  do i=1, N_colloids
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

  call h5gclose_f(colloids_group, error)

  call h5close_f(error)

end program try_all
