program setup_fluid
  use cell_system
  use particle_system
  use hilbert
  use interaction
  use hdf5
  use h5md_module
  use mpcd
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent

  type(h5md_file_t) :: datafile
  type(h5md_element_t) :: elem
  type(h5md_element_t) :: e_solvent, e_solvent_v
  integer(HID_T) :: box_group, solvent_group

  integer, parameter :: N = 2000

  integer :: i, L(3), seed_size, clock, error
  integer :: j
  integer, allocatable :: seed(:)

  double precision :: v_com(3)

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

  do i = 1, 20
     do j=1, solvent% Nmax
        solvent% vel(1, j) = localnormal()
        solvent% vel(2, j) = localnormal()
        solvent% vel(3, j) = localnormal()
     end do
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)
     solvent% vel = solvent% vel - spread(v_com, dim=2, ncopies=size(solvent% vel, dim=2))
     v_com = sum(solvent% vel, dim=2) / size(solvent% vel, dim=2)

     solvent% pos = solvent% pos + 0.1 * solvent% vel
     do j= 1, 3
        solvent% pos(:, j) = modulo( solvent% pos(:, j), L(j)*1.d0 )
     end do

     call solvent% sort(solvent_cells)
     call e_solvent% append(solvent% pos, i, i*1.d0)
     call e_solvent_v% append(solvent% vel, i, i*1.d0)

     call solvent_cells%count_particles(solvent% pos)
     write(*,*) compute_temperature(solvent, solvent_cells), sum(solvent% vel**2)/(3*solvent% Nmax)
  end do

  call e_solvent% close()
  call e_solvent_v% close()

  call h5gclose_f(solvent_group, error)

  call datafile% close()

  call h5close_f(error)

contains

  !> Normal random number
  !!
  !! Use the method from Marsaglia & Bray (1964) to generate numbers from a
  !! normal distribution.
  function localnormal() result(r)
    double precision :: r

    double precision :: u1, u2, radius
    double precision, save :: other
    logical :: found
    logical, save :: oneleft=.false.

    if ( oneleft ) then
       r = other
       oneleft = .false.
    else
       found = .false.
       do while (.not. found)
          call random_number(u1)
          u1 = 2*u1 - 1
          call random_number(u2)
          u2 = 2*u2 - 1
          radius = (u1**2+u2**2)
          if ( ( radius < 1 ) .and. (radius > 0) ) found = .true.
       end do
       other = u1 * sqrt( -2 * log(radius)/radius )
       r = u2 * sqrt( -2 * log(radius)/radius )
    end if
  end function localnormal

end program setup_fluid
