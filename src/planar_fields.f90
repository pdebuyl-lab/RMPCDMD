! This file is part of RMPCDMD
! Copyright (c) 2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Routines to compute planar concentration and velocity profiles

module planar_fields
  use common
  use particle_system
  use cell_system
  use neighbor_list
  implicit none

  private

  public :: planar_fields_t

  type planar_fields_t
     double precision, allocatable :: count(:,:,:)
     double precision, allocatable :: v(:,:,:,:)
     double precision :: x_min, dx
     double precision :: y_min, dy
     double precision :: thickness
     type(timer_t), pointer :: timer
   contains
     procedure :: init
     procedure :: update
  end type planar_fields_t

contains

  !> Initialize a planar_fields_t object
  subroutine init(this, N_species, N_x, x_min, x_max, N_y, y_min, y_max, thickness)
    class(planar_fields_t), intent(out) :: this
    integer, intent(in) :: N_species
    integer, intent(in) :: N_x
    double precision, intent(in) :: x_min
    double precision, intent(in) :: x_max
    integer, intent(in) :: N_y
    double precision, intent(in) :: y_min
    double precision, intent(in) :: y_max
    double precision, intent(in) :: thickness

    allocate(this%count(N_species, N_y, N_x))
    allocate(this%v(2, N_species, N_y, N_x))

    this%count = 0
    this%v = 0

    this%x_min = x_min
    this%dx = (x_max-x_min) / N_x
    this%y_min = y_min
    this%dy = (y_max-y_min) / N_y
    this%thickness = thickness

    allocate(this%timer)
    call this%timer%init('planar_update')

  end subroutine init

  !> Update the fields
  subroutine update(this, x, v, one_x, one_y, one_z, omega, solvent, cells)
    class(planar_fields_t), intent(inout) :: this
    double precision, intent(in) :: x(3)
    double precision, intent(in) :: v(3)
    double precision, intent(in) :: one_x(3)
    double precision, intent(in) :: one_y(3)
    double precision, intent(in) :: one_z(3)
    double precision, intent(in) :: omega(3)
    type(particle_system_t), intent(in) :: solvent
    type(cell_system_t), intent(in) :: cells

    double precision :: x_min, x_max, dx, y_min, y_max, dy, z_max
    double precision :: body_x, body_y, body_z, v_x, v_y
    double precision, dimension(3) :: d, L
    integer :: n_x, n_y, n_species
    double precision, dimension(3) :: solvent_v
    integer :: idx, s2, i_x, i_y

    integer, allocatable :: this_count(:,:,:)
    double precision, allocatable :: this_v(:,:,:,:)

    integer :: cell_idx, start, n

    L = cells%edges

    n_species = size(this%count, dim=1)
    n_y = size(this%count, dim=2)
    n_x = size(this%count, dim=3)

    x_min = this%x_min
    dx = this%dx
    x_max = x_min + n_x*dx
    y_min = this%y_min
    dy = this%dy
    y_max = y_min + n_y*dy
    z_max = this%thickness/2

    allocate(this_v(2, n_species, n_y, n_x))
    allocate(this_count(n_species, n_y, n_x))
    this_v = 0
    this_count = 0

    call this%timer%tic()
    !$omp parallel do &
    !$omp private(cell_idx, idx, start, n, s2, d, solvent_v, &
    !$omp body_x, body_y, body_z, v_x, v_y, i_x, i_y) &
    !$omp reduction(+:this_v) reduction(+:this_count)
    do cell_idx = 1, cells%N

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)

       do idx = start, start + n - 1
          s2 = solvent%species(idx)
          if (s2 <= 0) cycle
          d = rel_pos(solvent%pos(:,idx), x, L)
          body_z = dot_product(d, one_z)
          if (abs(body_z)>z_max) cycle
          body_x = dot_product(d, one_x)
          body_y = dot_product(d, one_y)

          if ( (body_x > x_min) .and. (body_x < x_max) .and. &
               (body_y > y_min) .and. (body_y < y_max) ) then

             i_x = floor((body_x-x_min)/dx)+1
             i_y = floor((body_y-y_min)/dy)+1
             solvent_v = solvent%vel(:,idx) - v
             solvent_v = solvent_v - cross(omega, d)
             v_x = dot_product(solvent_v, one_x)
             v_y = dot_product(solvent_v, one_y)
             this_count(s2, i_y, i_x) = this_count(s2, i_y, i_x) + 1
             this_v(:, s2, i_y, i_x) = this_v(:,s2, i_y, i_x) + [v_x, v_y]

          end if
       end do
    end do

    this%v = this%v + this_v
    this%count = this%count + this_count

    deallocate(this_v)
    deallocate(this_count)

    call this%timer%tac()

  end subroutine update

end module planar_fields
