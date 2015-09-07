module common
  implicit none
  private

  public :: rel_pos
  public :: profile_t
  public :: histogram_t

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

  type histogram_t
     double precision, allocatable :: data(:)
     double precision :: xmin
     double precision :: dx
     integer :: n
   contains
     generic, public :: init => histogram_init
     procedure, private :: histogram_init
     generic, public :: bin => histogram_bin
     procedure, private :: histogram_bin
  end type histogram_t

contains

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

  subroutine histogram_init(this, xmin, xmax, n)
    class(histogram_t), intent(out) :: this
    double precision, intent(in) :: xmin, xmax
    integer, intent(in) :: n

    this% n = n
    this% xmin = xmin
    this% dx = (xmax - xmin) / n
    if (this% dx <= 0) error stop 'negative step in histogram_init'
    allocate(this% data(n))
    this% data = 0

  end subroutine histogram_init

  subroutine histogram_bin(this, x)
    class(histogram_t), intent(inout) :: this
    double precision, intent(in) :: x

    integer :: idx

    idx = floor( (x - this% xmin) / this% dx ) + 1
    if ( ( idx > 0 ) .and. ( idx <= this% n ) ) then
       this% data(idx) = this% data(idx) + 1
    end if

  end subroutine histogram_bin

end module common
