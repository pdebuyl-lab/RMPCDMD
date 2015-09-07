program setup_fluid
  use common
  use cell_system
  use particle_system
  use hilbert
  use interaction
  use hdf5
  use h5md_module
  use mpcd
  use mt19937ar_module
  use iso_c_binding
  implicit none

  type(mt19937ar_t), target :: mt

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent

  type(profile_t) :: tz

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: e_solvent, e_solvent_v
  integer(HID_T) :: box_group, solvent_group

  integer, parameter :: N = 2000

  integer :: i, L(3), seed_size, clock, error
  integer, allocatable :: seed(:)

  double precision :: v_com(3), wall_v(3,2)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call system_clock(count=clock)
  call init_genrand(mt, int(clock, c_long))

  call h5open_f(error)

  L = [8, 5, 5]

  call solvent% init(N)

  call mt_normal_data(solvent% vel, mt)
  v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
  solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))

  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0)

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create('data_setup_simple_fluid.h5', 'RMPCDMD:setup_simple_fluid', '0.0 dev', 'Pierre de Buyl')

  call tz% init(0.d0, solvent_cells% edges(3), L(3))
  call h5gcreate_f(datafile% id, 'observables', datafile% observables, error)

  call h5gcreate_f(datafile% particles, 'solvent', solvent_group, error)
  call h5gcreate_f(solvent_group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call h5md_write_attribute(box_group, 'boundary', ['periodic', 'periodic', 'periodic'])
  call elem% create_fixed(box_group, 'edges', L*1.d0)
  call h5gclose_f(box_group, error)

  call e_solvent% create_time(solvent_group, 'position', solvent% pos, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_v% create_time(solvent_group, 'velocity', solvent% vel, ior(H5MD_TIME, H5MD_STORE_TIME))

  call solvent% sort(solvent_cells)

  wall_v = 0
  do i = 1, 200
     call simple_mpcd_step(solvent, solvent_cells, mt)
     call mpcd_stream(solvent, solvent_cells, 0.1d0)
     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)
  end do

  do i = 1, 200
     call simple_mpcd_step(solvent, solvent_cells, mt)
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)

     call mpcd_stream(solvent, solvent_cells, 0.1d0)

     call solvent% sort(solvent_cells)
     call e_solvent% append(solvent% pos, i, i*1.d0)
     call e_solvent_v% append(solvent% vel, i, i*1.d0)

     call solvent_cells%count_particles(solvent% pos)
     write(13,*) compute_temperature(solvent, solvent_cells, tz), sum(solvent% vel**2)/(3*solvent% Nmax), v_com

  end do

  call e_solvent% close()
  call e_solvent_v% close()

  call h5gclose_f(solvent_group, error)

  call tz% norm()
  call elem% create_fixed(datafile% observables, 'tz', tz% data)
  call elem% create_fixed(datafile% observables, 'tz_count', tz% count)

  call h5gclose_f(datafile% observables, error)

  call datafile% close()

  call h5close_f(error)

end program setup_fluid
