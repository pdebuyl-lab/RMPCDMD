module interaction
  use common
  use particle_system
  use neighbor_list
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
    this% cut = cut**2

  end subroutine lj_params_init

  subroutine compute_force(ps1, ps2, n_list, L, lj_params)
    type(particle_system_t), intent(inout) :: ps1
    type(particle_system_t), intent(inout) :: ps2
    type(neighbor_list_t), intent(in) :: n_list
    double precision, intent(in) :: L(3)
    type(lj_params_t), intent(in) :: lj_params

    integer :: i, j
    integer :: n_species1, n_species2
    double precision :: x(3), d(3), f(3), f1(3)
    integer :: s1, s2
    double precision :: r_sq

    n_species1 = ps1% n_species
    n_species2 = ps2% n_species

    do i = 1, ps1% Nmax
       if (ps1% species(i) <= 0) continue
       x = ps1% pos(:, i)
       s1 = ps1% species(i)
       f1 = 0
       do j = 1, n_list% n(i)
          s2 = ps2% species(j)
          if (s2 <= 0) continue
          d = rel_pos(x, ps2% pos(:, j), L)
          r_sq = sum(d**2)
          if ( r_sq <= lj_params% cut_sq(s1, s2) ) then
             f = lj_force(d, r_sq, lj_params% epsilon(s1, s2), lj_params% sigma(s1, s2))
             f1 = f1 + f
             ps2% force(:, j) = ps2% force(:, j) - f
          end if
       end do
       ps1% force(:, i) = ps1% force(:, i) + f1
    end do

  end subroutine compute_force

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
