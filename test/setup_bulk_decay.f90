! This file is part of RMPCDMD
! Copyright (c) 2015-2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program setup_bulk_decay
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

  integer :: rho, N, L(3)
  integer :: i, j, k
  integer :: error, clock, n_threads

  double precision, parameter :: tau=1.d0
  double precision :: bulk_rate
  double precision :: v_com(3)

  integer, parameter :: n_species = 2
  integer, dimension(n_species) :: n_solvent
  type(h5md_element_t) :: n_solvent_el

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, -4881001983349061041_c_int64_t)

  call h5open_f(error)

  L = 30
  rho = 9
  bulk_rate = 0.01
  N = rho*L(1)*L(2)*L(3)

  call solvent% init(N, n_species)

  do i=1, solvent% Nmax
     solvent% vel(1,i) = threefry_normal(state(1))
     solvent% vel(2,i) = threefry_normal(state(1))
     solvent% vel(3,i) = threefry_normal(state(1))
  end do
  v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
  solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))

  solvent% force = 0
  solvent% species = 2
  call solvent% random_placement(L*1.d0, state=state(1))

  call solvent_cells%init(L, 1.d0)
  solvent_cells%is_reac = .true.

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create('bulk_decay.h5', 'RMPCDMD:setup_bulk_decay', 'N/A', 'Pierre de Buyl')

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

  call n_solvent_el%create_time(datafile%observables, 'n_solvent', &
       n_solvent, H5MD_LINEAR, step=1, time=tau)

  call solvent% sort(solvent_cells)
  solvent_cells% origin(1) = threefry_double(state(1)) - 1
  solvent_cells% origin(2) = threefry_double(state(1)) - 1
  solvent_cells% origin(3) = threefry_double(state(1)) - 1

  do i = 1, 10
     call simple_mpcd_step(solvent, solvent_cells, state)
     call mpcd_stream_periodic(solvent, solvent_cells, tau)
     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call solvent% sort(solvent_cells)
  end do

  solvent% image = 0

  call update_n_solvent

  do i = 1, 200
     call simple_mpcd_step(solvent, solvent_cells, state)
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)

     call mpcd_stream_periodic(solvent, solvent_cells, tau)

     solvent_cells% origin(1) = threefry_double(state(1)) - 1
     solvent_cells% origin(2) = threefry_double(state(1)) - 1
     solvent_cells% origin(3) = threefry_double(state(1)) - 1
     call solvent% sort(solvent_cells)

     call bulk_reaction(solvent, solvent_cells, 2, 1, bulk_rate, tau, state)

     call update_n_solvent

     call e_solvent% append(solvent% pos, i, i*tau)
     call e_solvent_image% append(solvent% image, i, i*tau)
     call e_solvent_v% append(solvent% vel, i, i*tau)
     call e_solvent_spec% append(solvent% species, i, i*tau)
     call e_solvent_id% append(solvent% id, i, i*tau)

     write(*,*) i, compute_temperature(solvent, solvent_cells)
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
  call n_solvent_el%close()

  call h5gclose_f(solvent_group, error)

  call tz% norm()
  call elem% create_fixed(datafile% observables, 'tz', tz% data)
  call elem% create_fixed(datafile% observables, 'tz_count', tz% count)

  call h5gclose_f(datafile% observables, error)

  call datafile% close()

  call h5close_f(error)

contains

  subroutine update_n_solvent
    n_solvent = 0
    do k = 1, solvent%Nmax
       j = solvent%species(k)
       if (j <= 0) continue
       n_solvent(j) = n_solvent(j) + 1
    end do
    call n_solvent_el%append(n_solvent)
  end subroutine update_n_solvent

end program setup_bulk_decay
