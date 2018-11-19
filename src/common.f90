! This file is part of RMPCDMD
! Copyright (c) 2015-2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Utility routines
!!
!! This module contains routines of general use for RMPCDMD: minimum distance convention
!! routine, histogramming, command-line options processing, timers.
!!
!! A parameter constant for pi is also defined.

module common
  use iso_c_binding
  implicit none
  private

  public :: rel_pos
  public :: profile_t
  public :: histogram_t
  public :: switch
  public :: get_input_filename
  public :: get_input_args
  public :: timer_t
  public :: timer_list_t
  public :: args_t
  public :: cross
  public :: pi
  public :: alist_t
  public :: enzyme_kinetics_t
  public :: numbered_string

  ! bitmask data
  public :: REAC_BIT
  public :: MD_BIT
  public :: WALL_BIT
  public :: PAST_MD_BIT
  public :: OUTBOUND_BIT
  public :: CATALYZED_BIT
  public :: ENZYME_REGION_BIT
  public :: REAC_MASK
  public :: MD_MASK
  public :: WALL_MASK
  public :: CATALYZED_MASK
  public :: ENZYME_REGION_MASK

  integer, parameter :: max_path_length = 255
  double precision, parameter :: pi = 4*atan(1.d0)

  integer, parameter :: REAC_BIT = 1
  integer, parameter :: MD_BIT = 2
  integer, parameter :: WALL_BIT = 3
  integer, parameter :: PAST_MD_BIT = 4
  integer, parameter :: OUTBOUND_BIT = 5
  integer, parameter :: CATALYZED_BIT = 6
  integer, parameter :: ENZYME_REGION_BIT = 7

  integer, parameter :: REAC_MASK = ibset(0, REAC_BIT)
  integer, parameter :: MD_MASK = ibset(0, MD_BIT)
  integer, parameter :: WALL_MASK = ibset(0, WALL_BIT)
  integer, parameter :: CATALYZED_MASK = ibset(0, CATALYZED_BIT)
  integer, parameter :: ENZYME_REGION_MASK = ibset(0, ENZYME_REGION_BIT)

  !> Container for a profile, e.g. v(x)
  !!
  !! The result is \f$v(x) = \frac{\sum_i v_i \delta(x_i - x)}{\sum_i \delta(x_i - x)}\f$
  type profile_t
     double precision, allocatable :: data(:)
     integer, allocatable :: count(:)
     double precision :: xmin
     double precision :: dx
     integer :: n
   contains
     generic, public :: init => profile_init
     procedure, private :: profile_init
     generic, public :: bin => profile_bin
     procedure, private :: profile_bin
     generic, public :: reset => profile_reset
     procedure, private :: profile_reset
     generic, public :: norm => profile_norm
     procedure, private :: profile_norm
  end type profile_t

  !> Container for a histogram, e.g. p(x)
  !!
  !! The result is \f$p(x) = \frac{\sum_i \delta(x_i - x)}{N}\f$
  type histogram_t
     double precision, allocatable :: data(:,:)
     double precision :: xmin
     double precision :: dx
     integer :: n
     integer :: n_species
   contains
     generic, public :: init => histogram_init
     procedure, private :: histogram_init
     generic, public :: bin => histogram_bin
     procedure, private :: histogram_bin
  end type histogram_t

  type timer_t
     double precision :: tic_time
     double precision :: total
     character(len=32) :: name
   contains
     generic, public :: init => timer_init
     procedure, private :: timer_init
     procedure :: tic
     procedure :: tac
  end type timer_t

  type timer_pointer_t
     type(timer_t), pointer :: p
  end type timer_pointer_t

  type timer_list_t
     type(timer_pointer_t), allocatable :: timers(:)
     integer :: current_idx
   contains
     generic, public :: init => timer_list_init
     procedure, private :: timer_list_init
     generic, public :: append => timer_list_append
     procedure, private :: timer_list_append
     generic, public :: write => timer_list_write
     procedure, private :: timer_list_write
  end type timer_list_t

  interface switch
     module procedure :: switch_d2
     module procedure :: switch_i1
     module procedure :: switch_i2
  end interface switch

  !> Container for the standard command-line arguments to RMPCDMD
  type args_t
     character(len=max_path_length) :: input_file
     character(len=max_path_length) :: output_file
     integer(c_int64_t) :: seed
  end type args_t

  !> Appendable lists of double precision data
  type alist_t
     double precision, allocatable :: data(:)
     integer :: current_idx
     integer :: block_size
   contains
     generic, public :: init => alist_init
     procedure, private :: alist_init
     generic, public :: append => alist_append
     procedure, private :: alist_append
  end type alist_t

  !> Container for the list of times for enzymatic kinetics
  type enzyme_kinetics_t
     type(alist_t) :: bind_substrate
     type(alist_t) :: release_substrate
     type(alist_t) :: bind_product
     type(alist_t) :: release_product
  end type enzyme_kinetics_t

contains

  !> Return x-y distance with minimum image convention
  pure function rel_pos(x, y, L) result(r)
    double precision, intent(in) :: x(3), y(3), L(3)

    double precision :: r(3)
    integer :: dim

    r = x - y

    do dim=1,3
       if ( r(dim) < -0.5d0*L(dim) ) then
          r(dim) = r(dim) + L(dim)
       else if ( r(dim) > 0.5d0*L(dim) ) then
          r(dim) = r(dim) - L(dim)
       end if
    end do

  end function rel_pos

  subroutine profile_init(this, xmin, xmax, n)
    class(profile_t), intent(out) :: this
    double precision, intent(in) :: xmin, xmax
    integer, intent(in) :: n

    this% n = n
    this% xmin = xmin
    this% dx = (xmax - xmin) / n
    if (this% dx <= 0) error stop 'negative step in profile_init'
    allocate(this% data(n))
    this% data = 0
    allocate(this% count(n))
    this% count = 0

  end subroutine profile_init

  subroutine profile_bin(this, x, value)
    class(profile_t), intent(inout) :: this
    double precision, intent(in) :: x, value

    integer :: idx

    idx = floor( (x - this% xmin) / this% dx ) + 1
    if ( ( idx > 0 ) .and. ( idx <= this% n ) ) then
       this% data(idx) = this% data(idx) + value
       this% count(idx) = this% count(idx) + 1
    end if

  end subroutine profile_bin

  subroutine profile_reset(this)
    class(profile_t), intent(inout) :: this

    this% data = 0
    this% count = 0

  end subroutine profile_reset

  subroutine profile_norm(this)
    class(profile_t), intent(inout) :: this

    where (this% count > 0)
       this% data = this% data / this% count
    end where

  end subroutine profile_norm

  subroutine histogram_init(this, xmin, xmax, n, n_species)
    class(histogram_t), intent(out) :: this
    double precision, intent(in) :: xmin, xmax
    integer, intent(in) :: n
    integer, optional, intent(in) :: n_species

    if (present(n_species)) then
       this%n_species = n_species
    else
       this%n_species = 1
    end if

    this% n = n
    this% xmin = xmin
    this% dx = (xmax - xmin) / n
    if (this% dx <= 0) error stop 'negative step in histogram_init'
    allocate(this% data(this%n_species, n))
    this% data = 0

  end subroutine histogram_init

  subroutine histogram_bin(this, x, s)
    class(histogram_t), intent(inout) :: this
    double precision, intent(in) :: x
    integer, optional, intent(in) :: s

    integer :: idx, s_var

    if (present(s)) then
       s_var = s
    else
       s_var = 1
    end if

    idx = floor( (x - this% xmin) / this% dx ) + 1
    if ( ( idx > 0 ) .and. ( idx <= this% n ) ) then
       this% data(s_var, idx) = this% data(s_var, idx) + 1
    end if

  end subroutine histogram_bin

  subroutine switch_d2(p1, p2)
    double precision, pointer, dimension(:,:), intent(inout) :: p1, p2
    double precision, pointer, dimension(:,:) :: p

    p => p1
    p1 => p2
    p2 => p

  end subroutine switch_d2

  subroutine switch_i1(p1, p2)
    integer, pointer, dimension(:), intent(inout) :: p1, p2
    integer, pointer, dimension(:) :: p

    p => p1
    p1 => p2
    p2 => p

  end subroutine switch_i1

  subroutine switch_i2(p1, p2)
    integer, pointer, dimension(:,:), intent(inout) :: p1, p2
    integer, pointer, dimension(:,:) :: p

    p => p1
    p1 => p2
    p2 => p

  end subroutine switch_i2

  character(len=max_path_length) function get_input_filename() result(r)

    if (command_argument_count() < 1) then
       error stop 'missing argument for parameter file'
    end if

    call get_command_argument(1, r)

  end function get_input_filename

  type(args_t) function get_input_args() result(args)

    integer :: iostat
    character(max_path_length) :: r

    if (command_argument_count() /= 3) then
       write(*,*) 'Welcome to RMPCMD http://lab.pdebuyl.be/rmpcdmd/'
       write(*,*) 'Usage:'
       write(*,*) '    rmpcdmd run program_name input output seed'
       write(*,*) &
            '    where input is the filename of the parameters file, output is the name'
       write(*,*) &
            '    of the H5MD output file and seed is a signed 64-bit integer'
       error stop
    end if

    call get_command_argument(1, args%input_file)
    call get_command_argument(2, args%output_file)
    call get_command_argument(3, r)
    read(r, *, iostat=iostat) args%seed
    if (iostat /= 0) then
       write(*,*) 'Error when reading the seed value (third command-line argument)'
       error stop
    end if

  end function get_input_args

  subroutine timer_init(this, name, system_name)
    class(timer_t), intent(out) :: this
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(in) :: system_name

    if (present(system_name)) then
       if (len(system_name) > 0) then
          this%name = system_name//' '//name
       else
          this%name = name
       end if
    else
       this%name = name
    end if
    this%total = 0

  end subroutine timer_init

  subroutine tic(this)
    use omp_lib
    class(timer_t), intent(inout) :: this
    this%tic_time = omp_get_wtime()
  end subroutine tic

  subroutine tac(this)
    use omp_lib
    class(timer_t), intent(inout) :: this
    this%total = this%total + omp_get_wtime() - this%tic_time
  end subroutine tac

  subroutine timer_list_init(this, n)
    class(timer_list_t), intent(out) :: this
    integer, intent(in) :: n

    allocate(this%timers(n))
    this%current_idx = 0

  end subroutine timer_list_init

  subroutine timer_list_append(this, timer_target)
    class(timer_list_t), intent(inout) :: this
    type(timer_t), target, intent(in) :: timer_target

    this%current_idx = this%current_idx + 1
    if (this%current_idx > size(this%timers)) error stop 'exceeded timer_list_t size'
    this%timers(this%current_idx)%p => timer_target

  end subroutine timer_list_append

  subroutine timer_list_write(this, group, total_out)
    use hdf5, only: HID_T
    use h5md_module, only: h5md_write_dataset
    implicit none
    class(timer_list_t), intent(inout) :: this
    integer(HID_T), intent(inout) :: group
    double precision, optional, intent(out) :: total_out

    integer :: i

    total_out = 0
    do i = 1, size(this%timers)
       if (associated(this%timers(i)%p)) then
          call h5md_write_dataset(group, this%timers(i)%p%name, this%timers(i)%p%total)
          total_out = total_out + this%timers(i)%p%total
       end if
    end do

  end subroutine timer_list_write

  pure function cross(x1, x2) result(r)
    double precision, intent(in) :: x1(3), x2(3)
    double precision :: r(3)

    r(1) = x1(2)*x2(3) - x1(3)*x2(2)
    r(2) = x1(3)*x2(1) - x1(1)*x2(3)
    r(3) = x1(1)*x2(2) - x1(2)*x2(1)

  end function cross

  subroutine alist_init(this, block_size)
    class(alist_t), intent(out) :: this
    integer, intent(in) :: block_size

    allocate(this%data(block_size))
    this%current_idx = 0
    this%block_size = block_size

  end subroutine alist_init

  subroutine alist_append(this, value)
    class(alist_t), intent(inout) :: this
    double precision, intent(in) :: value

    integer :: idx, len
    double precision, allocatable :: tmp_data(:)

    idx = this%current_idx
    len = size(this%data)

    if (idx == len) then
       call move_alloc(this%data, tmp_data)
       allocate(this%data(len+this%block_size))
       this%data(1:len) = tmp_data
       deallocate(tmp_data)
    end if

    idx = idx + 1
    this%data(idx) = value
    this%current_idx = idx

  end subroutine alist_append

  function numbered_string(base_string, index, length) result(s)
    character(len=*), intent(in) :: base_string
    integer, intent(in) :: index
    integer, intent(in) :: length
    character(len=:), allocatable :: s

    character(len=12) format_string

    allocate(character(len=len(trim(base_string))+length) :: s)

    s(1:len(trim(base_string))) = base_string

    write(format_string, '(a,i3.3,a,i3.3,a)') '(a,i', length, '.', length, ')'
    write(s, format_string) trim(base_string), index

  end function numbered_string

end module common
