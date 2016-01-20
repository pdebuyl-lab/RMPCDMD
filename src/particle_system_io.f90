module particle_system_io
  use hdf5
  use h5md_module
  use particle_system
  implicit none

  private

  public :: particle_system_io_t
  public :: thermo_t

  type thermo_t
     integer :: n_buffer
     integer :: idx
     double precision, allocatable, dimension(:) :: temperature
     double precision, allocatable, dimension(:) :: potential_energy
     double precision, allocatable, dimension(:) :: kinetic_energy
     double precision, allocatable, dimension(:) :: internal_energy
   contains
     procedure :: init => thermo_init
     procedure :: append => thermo_append
  end type thermo_t

  type particle_system_io_t
     integer :: Nmax
     integer :: error
     integer(HID_T) :: group
     type(h5md_element_t) :: box
     type(h5md_element_t) :: position
     type(h5md_element_t) :: velocity
     type(h5md_element_t) :: force
     type(h5md_element_t) :: id
     type(h5md_element_t) :: species
   contains
     procedure :: init => ps_init
     procedure :: close => ps_close
  end type particle_system_io_t

contains

  subroutine ps_init(this, h5md_file, name, ps, position, position_step, position_time)
    class(particle_system_io_t), intent(out) :: this
    type(h5md_file_t), intent(inout) :: h5md_file
    character(len=*), intent(in) :: name
    type(particle_system_t), intent(in) :: ps
    integer, intent(in), optional :: position
    integer, intent(in), optional :: position_step
    double precision, intent(in), optional :: position_time

    this% Nmax = ps% Nmax

    call h5gcreate_f(h5md_file% particles, name, this% group, this% error)

    if (present(position)) then
       if (iand(position, H5MD_TIME) == H5MD_TIME) then
          call this% position% create_time(this% group, 'position', ps% pos, position)
       else if (iand(position, H5MD_LINEAR) == H5MD_LINEAR) then
          call this% position% create_time(this% group, 'position', ps% pos, position, position_step, position_time)
       else if (iand(position, H5MD_FIXED) == H5MD_FIXED) then
          call this% position% create_fixed(this% group, 'position', ps% pos)
       else
          stop 'unknown storage for position in particle_system_io init'
       end if
    end if
          

  end subroutine ps_init

  subroutine ps_close(this)
    class(particle_system_io_t), intent(inout) :: this

    call this% position% close()
    call h5gclose_f(this% group, this% error)

  end subroutine ps_close

  subroutine thermo_init(this, datafile, n_buffer, step, time)
    class(thermo_t), intent(out) :: this
    type(h5md_file_t), intent(inout) :: datafile
    integer, intent(in) :: n_buffer
    integer, intent(in) :: step
    double precision, intent(in) :: time

    type(h5md_element_t) :: e
    double precision :: dummy

    if (n_buffer <= 0) error stop 'n_buffer non-positive in thermo_init'

    this% n_buffer = n_buffer
    this% idx = 0

    allocate(this% temperature(n_buffer))
    allocate(this% potential_energy(n_buffer))
    allocate(this% kinetic_energy(n_buffer))
    allocate(this% internal_energy(n_buffer))

    call e%create_time(datafile%observables, 'temperature', dummy, H5MD_LINEAR, step, time)
    call e%close()
    call e%create_time(datafile%observables, 'potential_energy', dummy, H5MD_LINEAR, step, time)
    call e%close()
    call e%create_time(datafile%observables, 'kinetic_energy', dummy, H5MD_LINEAR, step, time)
    call e%close()
    call e%create_time(datafile%observables, 'internal_energy', dummy, H5MD_LINEAR, step, time)
    call e%close()

  end subroutine thermo_init

  subroutine thermo_append(this, datafile, temperature, potential_energy, kinetic_energy, internal_energy)
    class(thermo_t), intent(inout) :: this
    type(h5md_file_t), intent(inout) :: datafile
    double precision, intent(in) :: temperature, potential_energy, kinetic_energy, internal_energy

    integer :: i
    type(h5md_element_t) :: e

    i = this%idx + 1
    this%idx = i

    this%temperature(i) = temperature
    this%potential_energy(i) = potential_energy
    this%kinetic_energy(i) = kinetic_energy
    this%internal_energy(i) = internal_energy

    if (i == this%n_buffer) then
       call e%open_time(datafile%observables, 'temperature')
       call e%append_buffer(this%temperature)
       call e%close()
       call e%open_time(datafile%observables, 'potential_energy')
       call e%append_buffer(this%potential_energy)
       call e%close()
       call e%open_time(datafile%observables, 'kinetic_energy')
       call e%append_buffer(this%kinetic_energy)
       call e%close()
       call e%open_time(datafile%observables, 'internal_energy')
       call e%append_buffer(this%internal_energy)
       call e%close()
       this%idx = 0
    end if

  end subroutine thermo_append

end module particle_system_io