module interaction
  implicit none

  type lj_params_t
     integer :: n1, n2
     double precision, allocatable :: epsilon(:,:)
     double precision, allocatable :: sigma(:,:)
     double precision, allocatable :: sigma_sq(:,:)
     double precision, allocatable :: cut(:,:)
     double precision, allocatable :: cut_sq(:,:)
   contains
     procedure :: init => lj_params_init
  end type lj_params_t

contains

  subroutine lj_params_init(this, epsilon, sigma, cut)
    class(lj_params_t), intent(out) :: this
    double precision, intent(in) :: epsilon(:,:)
    double precision, intent(in) :: sigma(:,:)
    double precision, intent(in) :: cut(:,:)

    integer :: n1, n2

    if ( size(epsilon, 1) /= size(sigma, 1) .or. &
         size(epsilon, 1) /= size(cut, 1) .or. &
         size(epsilon, 2) /= size(sigma, 2) .or. &
         size(epsilon, 2) /= size(cut, 2) ) then
       stop 'unequal sizes for arguments to lj_params_init'
    end if

    n1 = size(epsilon, 1)
    n2 = size(epsilon, 2)
    this% n1 = n1
    this% n2 = n2

    allocate(this% epsilon(n1, n2))
    allocate(this% sigma(n1, n2))
    allocate(this% sigma_sq(n1, n2))
    allocate(this% cut(n1, n2))
    allocate(this% cut_sq(n1, n2))

    this% epsilon = epsilon
    this% sigma = sigma
    this% sigma_sq = sigma**2
    this% cut = cut
    this% cut_sq = cut**2

  end subroutine lj_params_init


  pure function lj_force(d, r_sq, epsilon, sigma) result(f)
    double precision, intent(in) :: d(3)
    double precision, intent(in) :: r_sq, epsilon, sigma
    double precision :: f(3)

    double precision :: sig6_o_r6, force_over_r

    sig6_o_r6 = sigma**6/r_sq**3

    force_over_r = 24.d0*epsilon* sig6_o_r6/r_sq * (2.d0*sig6_o_r6 - 1.d0)

    f = force_over_r*d

  end function lj_force

end module interaction
