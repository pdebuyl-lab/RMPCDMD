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

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent

  type(profile_t) :: tz
  type(histogram_t) :: rhoz

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: e_solvent, e_solvent_v
  type(h5md_element_t) :: elem_tz, elem_tz_count
  type(h5md_element_t) :: elem_rhoz
  type(h5md_element_t) :: elem_T
  integer(HID_T) :: box_group, solvent_group

  integer, parameter :: N = 10800

  integer :: i, L(3), seed_size, clock, error
  integer, allocatable :: seed(:)

  double precision :: v_com(3), wall_v(3,2), wall_t(2)
  double precision :: T

  double precision, parameter :: tau = 0.1d0
  type(threefry_rng_t), allocatable :: state(:)
  integer :: n_threads

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))

  do i = 1, n_threads
     state(i)%counter%c0 = 0
     state(i)%counter%c1 = 0
     state(i)%key%c0 = 0
     state(i)%key%c1 = 194387282_c_long
  end do

  call h5open_f(error)

  L = [6, 6, 30]

  call solvent% init(N)

  do i=1, solvent% Nmax
     solvent%vel(1,i) = threefry_normal(state(1))
     solvent%vel(2,i) = threefry_normal(state(1))
     solvent%vel(3,i) = threefry_normal(state(1))
  end do
  v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
  solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))

  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0)

  call solvent_cells%init(L, 1.d0, has_walls=.true.)
  solvent_cells% origin(3) = -0.5d0

  call solvent_cells%count_particles(solvent% pos)

  call datafile% create('data_setup_simple_fluid.h5', 'RMPCDMD:setup_simple_fluid', '0.0 dev', 'Pierre de Buyl')

  call tz% init(0.d0, solvent_cells% edges(3), L(3))
  call rhoz% init(0.d0, solvent_cells% edges(3), L(3))
  call h5gcreate_f(datafile% id, 'observables', datafile% observables, error)

  call h5gcreate_f(datafile% particles, 'solvent', solvent_group, error)
  call h5gcreate_f(solvent_group, 'box', box_group, error)
  call h5md_write_attribute(box_group, 'dimension', 3)
  call h5md_write_attribute(box_group, 'boundary', ['periodic', 'periodic', 'periodic'])
  call elem% create_fixed(box_group, 'edges', L*1.d0)
  call h5gclose_f(box_group, error)

  call e_solvent% create_time(solvent_group, 'position', solvent% pos, ior(H5MD_TIME, H5MD_STORE_TIME))
  call e_solvent_v% create_time(solvent_group, 'velocity', solvent% vel, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_tz% create_time(datafile% observables, 'tz', tz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_tz_count% create_time(datafile% observables, 'tz_count', tz% count, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_rhoz% create_time(datafile% observables, 'rhoz', rhoz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_T% create_time(datafile% observables, 'temperature', T, ior(H5MD_TIME, H5MD_STORE_TIME))

  call solvent% sort(solvent_cells)
  call solvent_cells%count_particles(solvent% pos)

  wall_v = 0
  wall_t = [0.9d0, 1.1d0]

  do i = 1, 1000
     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=[0.9d0, 1.1d0], wall_v=wall_v, wall_n=[10, 10])
     call mpcd_stream_zwall(solvent, solvent_cells, tau,[0d0,0d0,0d0])
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1
     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)
  end do

  do i = 1, 1000
     call wall_mpcd_step(solvent, solvent_cells, state, &
          wall_temperature=[0.9d0, 1.1d0], wall_v=wall_v, wall_n=[10, 10])
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)

     call mpcd_stream_zwall(solvent, solvent_cells, tau, [0d0,0d0,0d0])
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1

     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)

     T = compute_temperature(solvent, solvent_cells, tz)
     write(13,*) T, sum(solvent% vel**2)/(3*solvent% Nmax), v_com
     call elem_T% append(T, i, i*tau)

     call compute_rho(solvent, rhoz)

     if (modulo(i, 50) == 0) then
        call tz% norm()
        call elem_tz% append(tz% data, i, i*tau)
        call elem_tz_count% append(tz% count, i, i*tau)
        call tz% reset()
        rhoz% data = rhoz% data / (50.d0 * rhoz% dx)
        call elem_rhoz% append(rhoz% data, i, i*tau)
        rhoz% data = 0
     end if

  end do

  call e_solvent% append(solvent% pos, i, i*tau)
  call e_solvent_v% append(solvent% vel, i, i*tau)

  clock = 0
  do i = 1 , solvent_cells% N
     if ( solvent_cells% cell_count(i) > 0 ) clock = clock + 1
  end do
  write(*,*) clock, 'filled cells'
  write(*,*) L(1)*L(2)*(L(3)+1), 'actual cells'

  call e_solvent% close()
  call e_solvent_v% close()
  call elem_tz% close()
  call elem_rhoz% close()
  call elem_T% close()

  call h5gclose_f(solvent_group, error)

  call h5gclose_f(datafile% observables, error)

  call datafile% close()

  call h5close_f(error)

end program setup_fluid
