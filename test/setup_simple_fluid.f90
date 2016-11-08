! This file is part of RMPCDMD
! Copyright (c) 2015-2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program setup_fluid
  use common
  use cell_system
  use particle_system
  use hilbert
  use interaction
  use hdf5
  use h5md_module
  use mpcd
  use threefry_module
  use iso_c_binding
  use omp_lib
  implicit none

  type(threefry_rng_t), target, allocatable :: state(:)

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent

  type(profile_t) :: tz

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: e_solvent, e_solvent_v, e_solvent_spec, e_solvent_id
  type(h5md_element_t) :: e_solvent_image
  integer(HID_T) :: box_group, solvent_group

  integer, parameter :: N = 2000

  integer :: i, L(3), error, clock, n_threads

  double precision, parameter :: tau=0.1d0
  double precision :: v_com(3), wall_v(3,2)

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, 2985348049158689438_c_int64_t)

  call h5open_f(error)

  L = [8, 5, 5]

  call solvent% init(N)

  do i=1, solvent% Nmax
     solvent% vel(1,i) = threefry_normal(state(1))
     solvent% vel(2,i) = threefry_normal(state(1))
     solvent% vel(3,i) = threefry_normal(state(1))
  end do
  v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
  solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))

  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0, state=state(1))

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create('data_setup_simple_fluid.h5', 'RMPCDMD:setup_simple_fluid', '0.0 dev', 'Pierre de Buyl')

  call tz% init(0.d0, solvent_cells% edges(3), L(3))

  call h5gcreate_f(datafile% particles, 'solvent', solvent_group, error)
  call h5gcreate_f(solvent_group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call h5md_write_attribute(box_group, 'boundary', ['periodic', 'periodic', 'periodic'])
  call elem% create_fixed(box_group, 'edges', L*1.d0)
  call h5gclose_f(box_group, error)

  call e_solvent% create_time(solvent_group, 'position', solvent% pos, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_v% create_time(solvent_group, 'velocity', solvent% vel, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_spec% create_time(solvent_group, 'species', solvent% species, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_id% create_time(solvent_group, 'id', solvent% id, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_image% create_time(solvent_group, 'image', solvent% image, ior(H5MD_TIME, H5MD_STORE_TIME))

  call solvent% sort(solvent_cells)
  call solvent_cells%count_particles(solvent% pos)
  solvent_cells% origin(1) = threefry_double(state(1)) - 1
  solvent_cells% origin(2) = threefry_double(state(1)) - 1
  solvent_cells% origin(3) = threefry_double(state(1)) - 1

  wall_v = 0
  do i = 1, 200
     call simple_mpcd_step(solvent, solvent_cells, state)
     call mpcd_stream_periodic(solvent, solvent_cells, tau)
     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)
  end do
  do i = 1,solvent% Nmax
     if (solvent% pos(3,i) <= solvent_cells% edges(3)/2) then
        solvent% species(i) = 1
     else 
        solvent% species(i) = 2
     end if 
  end do

  solvent% image = 0

  do i = 1, 200
     call simple_mpcd_step(solvent, solvent_cells, state)
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)

     call mpcd_stream_periodic(solvent, solvent_cells, tau)

     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)
     call e_solvent% append(solvent% pos, i, i*tau)
     call e_solvent_image% append(solvent% image, i, i*tau)
     call e_solvent_v% append(solvent% vel, i, i*tau)
     call e_solvent_spec% append(solvent% species, i, i*tau)
     call e_solvent_id% append(solvent% id, i, i*tau)

     write(13,*) compute_temperature(solvent, solvent_cells, tz), sum(solvent% vel**2)/(3*solvent% Nmax), v_com

  end do

  clock = 0
  do i = 1 , solvent_cells% N
     if ( solvent_cells% cell_count(i) > 0 ) clock = clock + 1
  end do
  write(*,*) clock, 'filled cells'
  write(*,*) L(1)*L(2)*L(3), 'actual cells'

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
