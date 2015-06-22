module common
  implicit none
  private

  public :: rel_pos

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

end module common
