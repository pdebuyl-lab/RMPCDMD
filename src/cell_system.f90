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
     integer, allocatable :: cell_count_sum(:)
     integer :: M(3)
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

    integer :: i

    this%L = L
    this%a = a
    this%M = 1
    this%origin = 0.d0

    do i=1, 3
       do while ( 2**this%M(i) < L(i) )
          this%M(i) = this%M(i)+1
       end do
    end do
    this%N = 2**this%M(1)*2**this%M(2)*2**this%M(3)
    allocate(this%cell_count(this%N))
    allocate(this%cell_count_sum(this%N))

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
       idx = compact_p_to_h(p, this%M) + 1
       this%cell_count(idx) = this%cell_count(idx) + 1
    end do

    this%cell_count_sum = 0
    this%cell_count_sum(1) = this%cell_count(1)
    do i=2, this%N
       this%cell_count_sum(i) = this%cell_count_sum(i-1) + this%cell_count(i)
    end do

  end subroutine count_particles

end module cell_system
