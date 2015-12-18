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
  type(histogram_t) :: rhoz
  type(profile_t) :: vx

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: e_solvent, e_solvent_v
  type(h5md_element_t) :: elem_tz, elem_tz_count, elem_vx_count
  type(h5md_element_t) :: elem_rhoz
  type(h5md_element_t) :: elem_vx
  type(h5md_element_t) :: elem_T
  integer(HID_T) :: box_group, solvent_group

  integer, parameter :: N = 10800

  integer :: i, L(3), seed_size, clock, error
  integer, allocatable :: seed(:)

  double precision :: v_com(3), wall_v(3,2), wall_t(2)
  double precision, parameter :: gravity_field(3) = [ 0.001d0, 0.d0, 0.d0 ]
  double precision :: T

  double precision, parameter :: set_temperature = 1
  double precision, parameter :: tau = 0.1d0

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call system_clock(count=clock)
  call init_genrand(mt, int(clock, c_long))

  call h5open_f(error)

  L = [6, 6, 30]

  call solvent% init(N)

  call mt_normal_data(solvent% vel, mt)
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
  call vx% init(0.d0, solvent_cells% edges(3), L(3))
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
  call elem_vx% create_time(datafile% observables, 'vx', tz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_vx_count% create_time(datafile% observables, 'vx_count', vx% count, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_rhoz% create_time(datafile% observables, 'rhoz', rhoz% data, ior(H5MD_TIME, H5MD_STORE_TIME))
  call elem_T% create_time(datafile% observables, 'temperature', T, ior(H5MD_TIME, H5MD_STORE_TIME))

  call solvent% sort(solvent_cells)
  call solvent_cells%count_particles(solvent% pos)

  wall_v = 0
  wall_t = [1.0d0, 1.0d0]
  do i = 1, 1000
     call wall_mpcd_step(solvent, solvent_cells, mt, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[10, 10], bulk_temperature=set_temperature)
     call mpcd_stream_zwall(solvent, solvent_cells, tau, gravity_field)
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1
     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)
  end do

  do i = 1, 1000
     call wall_mpcd_step(solvent, solvent_cells, mt, &
          wall_temperature=wall_t, wall_v=wall_v, wall_n=[10, 10], bulk_temperature=set_temperature)
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)

     call mpcd_stream_zwall(solvent, solvent_cells, tau, gravity_field)
     call random_number(solvent_cells% origin)
     solvent_cells% origin = solvent_cells% origin - 1

     call solvent% sort(solvent_cells)
     call solvent_cells%count_particles(solvent% pos)

     T = compute_temperature(solvent, solvent_cells, tz)
     write(13,*) T, sum(solvent% vel**2)/(3*solvent% Nmax), v_com
     call elem_T% append(T, i, i*tau)

     call compute_rho(solvent, rhoz)
     call compute_vx(solvent, vx)

     if (modulo(i, 50) == 0) then
        call tz% norm()
        call elem_tz% append(tz% data, i, i*tau)
        call elem_tz_count% append(tz% count, i, i*tau)
        call tz% reset()
        call vx% norm()
        call elem_vx% append(vx% data, i, i*tau)
        call elem_vx_count% append(vx% count, i, i*tau)
        call vx% reset()
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
  call elem_tz_count% close()
  call elem_rhoz% close()
  call elem_vx% close()
  call elem_vx_count% close()
  call elem_T% close()

  call h5gclose_f(solvent_group, error)

  call h5gclose_f(datafile% observables, error)

  call datafile% close()

  call h5close_f(error)

end program setup_fluid
