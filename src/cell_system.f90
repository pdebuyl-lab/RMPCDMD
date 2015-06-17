module cell_system
  implicit none

  private

  public :: cell_system_t

  type cell_system_t
     integer :: L(3)
     integer :: N
     double precision :: a(3)
     double precision :: origin(3)
     integer, allocatable :: cell_count(:)
   contains
     procedure :: init
     procedure :: del
  end type cell_system_t

contains

  subroutine init(this, L, a)
    class(cell_system_t), intent(out) :: this
    integer, intent(in) :: L(3)
    double precision, intent(in) :: a

    this%L = L
    this%N = L(1)*L(2)*L(3)
    this%a = a
    this%origin = 0.d0
    allocate(this%cell_count(this%N))

  end subroutine init

  subroutine del(this)
    class(cell_system_t), intent(inout) :: this

    if (allocated(this%cell_count)) then
       deallocate(this%cell_count)
    end if

  end subroutine del

end module cell_system
