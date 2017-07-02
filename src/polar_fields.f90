! This file is part of RMPCDMD
! Copyright (c) 2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Routines to compute polar concentration and velocity profiles
!!
!! Polar data is stored according to a given orientation axis. The coordinates are defined
!! as
!!
!! \f[
!!   \left\{\begin{array}{ll}
!!   x &= r \cos\varphi \sin\theta\cr
!!   y &= r \sin\varphi \sin\theta\cr
!!   z &= r \cos\theta\cr
!!   \end{array}\right.
!! \f]
!!
!! Only \f$r\f$ and \f$\theta\f$ are considered as the axial symmetry prevents the
!! definition a third degree of freedom from the orientation vector.
!!
!! Knowing the orientation unit vector \f$\hat u\f$ and the relative coordinates \f$\vec
!! x\f$ of a point, the polar coordinates are
!!
!! \f[
!!   \left\{\begin{array}{ll}
!!   r &= |x|\cr
!!   \theta &= \arccos\left( \frac{\hat u \cdot \vec x}{|x|} \right)
!!   \end{array}\right.
!! \f]

module polar_fields
  use common
  use particle_system
  use cell_system
  use neighbor_list
  implicit none

  private

  public :: polar_fields_t

  type polar_fields_t
     double precision, allocatable :: c(:,:,:)
     integer, allocatable :: count(:,:,:)
     double precision, allocatable :: v(:,:,:,:)
     double precision :: r_min
     double precision :: r_max
     double precision :: dr
     double precision :: dtheta
     type(timer_t), pointer :: time_polar_update
   contains
     procedure :: init
     procedure :: update
  end type polar_fields_t

contains

  !> Initialize a polar_fields_t object
  subroutine init(this, N_species, N_r, r_min, r_max, N_theta)
    class(polar_fields_t), intent(out) :: this
    integer, intent(in) :: N_species
    integer, intent(in) :: N_r
    double precision, intent(in) :: r_min
    double precision, intent(in) :: r_max
    integer, intent(in) :: N_theta

    allocate(this%c(N_species, N_theta, N_r))
    allocate(this%count(N_species, N_theta, N_r))
    allocate(this%v(2, N_species, N_theta, N_r))

    this%c = 0
    this%count = 0
    this%v = 0

    this%r_min = r_min
    this%r_max = r_max
    this%dr = (r_max - r_min) / N_r
    this%dtheta = pi/N_theta

    allocate(this%time_polar_update)
    call this%time_polar_update%init('polar_update')

  end subroutine init

  !> Update the fields
  subroutine update(this, x, v, unit_r, solvent, cells)
    class(polar_fields_t), intent(inout) :: this
    double precision, intent(in) :: x(3)
    double precision, intent(in) :: v(3)
    double precision, intent(in) :: unit_r(3)
    type(particle_system_t), intent(in) :: solvent
    type(cell_system_t), intent(in) :: cells

    double precision :: r_min, r_min_sq, r_max, r_max_sq, dr, dtheta
    double precision :: r_sq, v_r, v_th, r, theta
    double precision, dimension(3) :: d, one_r, out_of_plane, one_theta, L
    double precision, dimension(3) :: solvent_v
    integer :: idx, s2, i_r, i_th

    integer :: cell_idx, start, n

    r_min = this%r_min
    r_min_sq = r_min**2
    r_max = this%r_max
    r_max_sq = r_max**2
    dr = this%dr
    dtheta = this%dtheta
    L = cells%edges

    call this%time_polar_update%tic()
    !$omp parallel do &
    !$omp private(cell_idx, idx, start, n, s2, d, r_sq, r, i_r, theta, i_th, one_r, &
    !$omp out_of_plane, one_theta, solvent_v, v_th)
    do cell_idx = 1, cells%N

       ! TODO: check if minimum distance from x to cell corners is below r_max

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)

       do idx = start, start + n - 1
          s2 = solvent%species(idx)
          if (s2 <= 0) cycle
          d = rel_pos(solvent%pos(:,idx), x, L)
          r_sq = sum(d**2)
          if (( r_sq >= r_min_sq ) .and. ( r_sq < r_max_sq )) then
             ! compute r, i_r, theta, i_theta
             r = sqrt(r_sq)
             i_r = floor((r-r_min)/dr) + 1
             theta = acos( dot_product(unit_r,d)/r )
             i_th = floor(theta/dtheta) + 1
             this%c(s2, i_th, i_r) = this%c(s2, i_th, i_r) + 1
             one_r = d/r
             ! compute v_r, v_theta
             out_of_plane = cross(unit_r, one_r)
             one_theta = cross(out_of_plane, one_r)
             one_theta = one_theta/norm2(one_theta)
             solvent_v = solvent%vel(:,idx) - v
             v_r = dot_product(solvent_v, one_r)
             v_th = dot_product(solvent_v, one_theta)
             !$omp atomic
             this%v(1, s2, i_th, i_r) = this%v(1, s2, i_th, i_r) + v_r
             !$omp atomic
             this%v(2, s2, i_th, i_r) = this%v(1, s2, i_th, i_r) + v_th
             !$omp atomic
             this%count(s2, i_th, i_r) = this%count(s2, i_th, i_r) + 1
          end if
       end do
    end do
    call this%time_polar_update%tac()

  end subroutine update

end module polar_fields
