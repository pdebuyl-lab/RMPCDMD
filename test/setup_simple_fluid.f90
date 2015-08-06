program setup_fluid
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
  type(neighbor_list_t) :: neigh

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: e_solvent, e_solvent_v
  integer(HID_T) :: box_group, solvent_group

  integer, parameter :: N = 1000

  integer :: i, L(3), seed_size, clock, error
  integer :: j
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call h5open_f(error)

  L = [8, 5, 5]

  call solvent% init(N)

  call random_number(solvent% vel(:, :))
  solvent% vel(:, :) = solvent% vel(:, :) - 0.5d0
  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0)

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create('data_setup_simple_fluid.h5', 'RMPCDMD:setup_simple_fluid', '0.0 dev', 'Pierre de Buyl')

  call h5gcreate_f(datafile% particles, 'solvent', solvent_group, error)
  call h5gcreate_f(solvent_group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call h5md_write_attribute(box_group, 'boundary', ['periodic', 'periodic', 'periodic'])
  call elem% create_fixed(box_group, 'edges', L*1.d0)
  call h5gclose_f(box_group, error)

  call e_solvent% create_time(solvent_group, 'position', solvent% pos, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_v% create_time(solvent_group, 'velocity', solvent% vel, ior(H5MD_TIME, H5MD_STORE_TIME))

  call solvent% sort(solvent_cells)

  call e_solvent% append(solvent% pos, 0, 0.d0)
  call e_solvent_v% append(solvent% vel, 0, 0.d0)

  do i = 1, 50
     call random_number(solvent% vel)
     solvent% vel = solvent% vel - 0.5d0
     solvent% pos = solvent% pos + 0.1 * solvent% vel
     do j= 1, 3
        solvent% pos(:, j) = modulo( solvent% pos(:, j), L(j)*1.d0 )
     end do

     call solvent% sort(solvent_cells)
     call e_solvent% append(solvent% pos, i, i*1.d0)
     call e_solvent_v% append(solvent% vel, i, i*1.d0)
  end do

  call e_solvent% close()
  call e_solvent_v% close()

  call h5gclose_f(solvent_group, error)

  call datafile% close()

  call h5close_f(error)

end program setup_fluid
