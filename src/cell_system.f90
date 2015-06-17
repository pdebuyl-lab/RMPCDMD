module cell_system
  use hilbert
  implicit none

  private

  public :: cell_system_t

  type cell_system_t
     integer :: L(3)
     integer :: N
     double precision :: a(3)
     double precision :: origin(3)
     integer, allocatable :: cell_count(:)
     type(hilbert_t) :: hilbert
   contains
     procedure :: init
     procedure :: del
     procedure :: count_particles
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

  subroutine count_particles(this, position)
    class(cell_system_t), intent(inout) :: this
    double precision, intent(in) :: position(:, :)

    integer :: i, idx, N_particles
    integer :: p(3)

    N_particles = size(position, 2)

    this%cell_count = 0

    do i=1, N_particles
       p = floor( (position(:, i) - this%origin ) / this%a )
       idx = this%hilbert%p_to_h(p)
       this%cell_count(idx) = this%cell_count(idx) + 1
    end do

  end subroutine count_particles

end module cell_system
