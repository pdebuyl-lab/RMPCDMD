module interaction
  implicit none

  private

  public :: lj_params_t
  public :: lj_force
  public :: lj_energy
  public :: lj_force_9_6
  public :: lj_energy_9_6
  public :: lj_force_9_3
  public :: lj_energy_9_3

  type lj_params_t
     integer :: n1, n2
     double precision, allocatable :: epsilon(:,:)
     double precision, allocatable :: sigma(:,:)
     double precision, allocatable :: sigma_sq(:,:)
     double precision, allocatable :: cut(:,:)
     double precision, allocatable :: cut_sq(:,:)
     double precision, allocatable :: cut_energy(:,:)
     double precision, allocatable :: shift
   contains
     procedure :: init => lj_params_init
  end type lj_params_t

contains

  subroutine lj_params_init(this, epsilon, sigma, cut, shift)
    class(lj_params_t), intent(out) :: this
    double precision, intent(in) :: epsilon(:,:)
    double precision, intent(in) :: sigma(:,:)
    double precision, intent(in) :: cut(:,:)
    double precision, intent(in), optional :: shift

    integer :: n1, n2

    if ( size(epsilon, 1) /= size(sigma, 1) .or. &
         size(epsilon, 1) /= size(cut, 1) .or. &
         size(epsilon, 2) /= size(sigma, 2) .or. &
         size(epsilon, 2) /= size(cut, 2) ) then
       error stop 'unequal sizes for arguments to lj_params_init'
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
    allocate(this% cut_energy(n1, n2))
    allocate(this% shift)

    this% epsilon = epsilon
    this% sigma = sigma
    this% sigma_sq = sigma**2
    this% cut = cut
    this% cut_sq = cut**2
    if (present(shift)) then
       this% shift = shift
    else
       this% shift = 0.d0
    end if

    this% cut_energy = lj_energy(this%cut_sq, this%epsilon, this%sigma, this%sigma*0)

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

  !interaction of colloid with the wall (9-6 LJ from SklogWiki)
  pure function lj_force_9_6(d, r_sq, epsilon, sigma) result(f)
     double precision, intent(in) :: d(3)
     double precision, intent(in) :: r_sq, epsilon, sigma
     double precision :: f(3)

     double precision :: sig6_o_r6, sig3_o_r3, force_o_r
     sig6_o_r6 = sigma**6/r_sq**3
     sig3_o_r3 = sigma**3/r_sq**(3.d0/2.d0)

     force_o_r = 6.75d0*3.d0*epsilon*sig6_o_r6/r_sq*(3.d0*sig3_o_r3-2.d0)

     f = force_o_r*d

  end function lj_force_9_6

  pure function lj_force_9_3(z, z_sq, epsilon, sigma) result(e)
     double precision, intent(in) ::z, z_sq, epsilon, sigma
     double precision :: e

     double precision :: sig3_o_r3, force_o_z
     
     sig3_o_r3 = sigma**3/z_sq**(3.d0/2.d0)

     force_o_z = 9*sqrt(3.d0)/2*epsilon/z_sq* (3.d0*sig3_o_r3**3 - sig3_o_r3)
     e = force_o_z*z
  end function lj_force_9_3

  pure elemental function lj_energy(r_sq, epsilon, sigma, cut_energy) result(e)
    double precision, intent(in) :: r_sq, epsilon, sigma, cut_energy
    double precision :: e

    double precision :: sig6_o_r6

    sig6_o_r6 = sigma**6/r_sq**3

    e = 4.d0*epsilon * (sig6_o_r6**2 - sig6_o_r6) - cut_energy

  end function lj_energy

  pure function lj_energy_9_6(r_sq, epsilon, sigma) result(e)
    double precision, intent(in) :: r_sq, epsilon, sigma
    double precision :: e

    double precision :: sig6_o_r6, sig3_o_r3

    sig6_o_r6 = sigma**6/r_sq**3
    sig3_o_r3 = sigma**3/r_sq**(3.d0/2.d0)

    e = 6.75d0*epsilon * (sig6_o_r6*(sig3_o_r3 - 1.d0) + 4.d0/27.d0)

  end function lj_energy_9_6

  pure function lj_energy_9_3(z_sq, epsilon, sigma) result(e)
     double precision, intent(in) :: z_sq, epsilon, sigma
     double precision :: e

     double precision :: sig3_o_r3
     
     sig3_o_r3 = sigma**3/z_sq**(3.d0/2.d0)

     e = 3*sqrt(3.d0)/2*epsilon*(sig3_o_r3**3 - sig3_o_r3) + epsilon
  end function lj_energy_9_3

end module interaction
