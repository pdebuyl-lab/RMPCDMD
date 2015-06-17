module hilbert
  implicit none
  private

  public :: hilbert_t
  
  integer, parameter :: D = 3

  type hilbert_t
     integer :: M(D)
   contains
     procedure :: p_to_h
  end type hilbert_t

contains

  subroutine init(this, M)
    class(hilbert_t), intent(out) :: this
    integer, intent(in) :: M(D)

    this%M = M

  end subroutine init

  function p_to_h(this, p)
    class(hilbert_t), intent(in) :: this
    integer, intent(in) :: p(D)
    integer :: p_to_h

    p_to_h = 1

  end function p_to_h

  function pure rotate_right(x, d)
    integer, intent(in) :: x
    integer, intent(in) :: d
    integer :: rotate_right

  end function pure

end module hilbert
