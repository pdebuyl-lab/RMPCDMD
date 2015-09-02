module cell_system
  use hilbert
  implicit none

  private

  public :: cell_system_t

  type cell_system_t
     integer :: L(3)
     integer :: N
     double precision :: a
     double precision :: origin(3)
     integer, allocatable :: cell_count(:)
     integer, allocatable :: cell_start(:)
     integer :: M(3)
     logical :: has_walls
   contains
     procedure :: init
     procedure :: count_particles
     procedure :: sort_particles
  end type cell_system_t

contains

  subroutine init(this, L, a, has_walls)
    class(cell_system_t), intent(out) :: this
    integer, intent(in) :: L(3)
    double precision, intent(in) :: a
    logical, intent(in), optional :: has_walls

    integer :: i

    this%L = L
    if (present(has_walls)) then
       this% has_walls = has_walls
    else
       this% has_walls = .false.
    end if
    if (this% has_walls) then
       this%L(3) = L(3) + 1
    end if
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
    allocate(this%cell_start(this%N))

  end subroutine init

  subroutine count_particles(this, position)
    class(cell_system_t), intent(inout) :: this
    double precision, intent(in) :: position(:, :)

    integer :: i, idx, N_particles
    integer :: p(3)

    N_particles = size(position, 2)

    this%cell_count = 0

    !$omp parallel do private(p, idx)
    do i=1, N_particles
       p = floor( (position(:, i) - this%origin ) / this%a )
       idx = compact_p_to_h(p, this%M) + 1
       !$omp atomic
       this%cell_count(idx) = this%cell_count(idx) + 1
    end do

    this%cell_start(1) = 1
    do i=2, this%N
       this%cell_start(i) = this%cell_start(i-1) + this%cell_count(i-1)
    end do

  end subroutine count_particles

  subroutine sort_particles(this, position_old, position_new)
    class(cell_system_t), intent(inout) :: this
    double precision, intent(in) :: position_old(:, :)
    double precision, intent(out) :: position_new(:, :)

    integer :: i, idx, N, start, p(3)

    N = size(position_old, 2)

    !$omp parallel do private(p, idx, start)
    do i=1, N
       p = floor( (position_old(:, i) - this%origin ) / this%a )
       idx = compact_p_to_h(p, this%M) + 1
       !$omp atomic capture
       start = this%cell_start(idx)
       this%cell_start(idx) = this%cell_start(idx) + 1
       !$omp end atomic
       position_new(:, start) = position_old(:, i)
    end do

  end subroutine sort_particles

end module cell_system
