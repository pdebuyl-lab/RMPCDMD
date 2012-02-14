
!> The polar_fields module defines a datatype for storing a polar concentration
!! field, a polar temperature field and a velocity field.

module polar_fields
  use MD, only: rel_pos
  use MPCD, only: L, PI
  use h5md
  implicit none

  !> A polar_fields_t contains the parameters for the distributions and the
  !! distributions themselves.
  type polar_fields_t
     !> Number of binnings.
     integer :: count
     !> Radial size of the bins.
     double precision :: dr
     !> Number of radial bins.
     integer :: Nr
     !> Angular size of the bins.
     double precision :: dtheta
     !> Number of angular bins.
     integer :: Nth
     !> Number of species.
     integer :: Ns
     !> Concentration field.
     double precision, allocatable :: c(:,:,:)
     !> Temperature field.
     double precision, allocatable :: t(:,:)
     !> Count for the temperature field.
     integer, allocatable :: t_count(:,:)
     !> Velocity field.
     double precision, allocatable :: v(:,:,:,:)
     !> Count for the velocity field.
     integer, allocatable :: v_count(:,:,:)
     !> h5md_t variable for the concentration.
     type(h5md_t) :: c_ID
     !> h5md_t variable for the velocity.
     type(h5md_t) :: v_ID
     !> h5md_t variable for the temperature.
     type(h5md_t) :: t_ID
  end type polar_fields_t

contains

  !> The init_polar_fields subroutine initializes a polar_fields_t variable.
  !! @param pft polar_fields_t variable
  !! @param Ns The number of species to take into account.
  !! @param Nr Number of points in the radial direction.
  !! @param dr Step in the radial direction.
  !! @param Nth Number of points in the angular direction.
  subroutine init_polar_fields(pft, Ns, Nr, dr, Nth, file_id, name)
    implicit none
    type(polar_fields_t), intent(out) :: pft
    integer, intent(in) :: Ns, Nr, Nth
    double precision, intent(in) :: dr
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name

    integer(HID_T) :: g_id
    
    pft % Nr = Nr
    pft % Nth = Nth
    pft % Ns = Ns

    pft % dr = dr
    pft % dtheta = 4.d0*atan(1.d0) / dble(Nth)

    ! allocate fields
    allocate(pft % c(Ns, Nth, Nr))
    allocate(pft % t(Nth, Nr))
    allocate(pft % t_count(Nth, Nr))
    allocate(pft % v(3, Ns, Nth, Nr))
    allocate(pft % v_count(Ns, Nth, Nr))
    ! reset fields and counters
    pft % c = 0.d0
    pft % v = 0.d0
    pft % v_count = 0
    pft % t = 0.d0
    pft % t_count = 0
    pft % count = 0

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)
    call h5gclose_f(g_id, h5_error)
    call h5md_create_obs(file_id, name//'/concentration', pft % c_ID, pft % c)
    call h5md_create_obs(file_id, name//'/velocity', pft % v_ID, pft % v, link_from=name//'/concentration')
    call h5md_create_obs(file_id, name//'/temperature', pft % t_ID, pft % t, link_from=name//'/concentration')
    
  end subroutine init_polar_fields

  !> Bins the different fields in a polar_fields_t variable.
  !! @param pft polar_fields_t variable
  !! @param x0 Center of the coordinates.
  !! @param dir Axis for the coordinates.
  subroutine update_polar_fields(pft, x0, dir)
    use MPCD, only: par_list, Vcom, so_r, so_v, so_species, so_sys, indices, oo_a, N_cells, cell_center
    implicit none
    type(polar_fields_t), intent(inout) :: pft
    double precision, intent(in) :: x0(3), dir(3)

    integer :: ci, cj, ck, mi, mj, mk, i, part, part_s, cc(3), t_count
    double precision :: norm, r, th, t
    integer :: cr_idx, cth_idx, extent
    integer :: r_idx, th_idx, sr_idx, sth_idx, sp_idx
    double precision, dimension(3) :: unit_r, x, sr, sth, sp, v, v_pol, x_cell
    logical :: do_t

    pft % count = pft % count + 1

    ! making the r-vector of norm 1
    norm = sqrt( sum( dir**2 ) )
    unit_r = dir / norm
    
    extent = ceiling( (pft % dr * pft % Nr) * oo_a  )

    call indices(x0, cc)
    do ck= cc(3) - extent, cc(3) + extent
       do cj = cc(2) - extent, cc(2) + extent
          do ci = cc(1) - extent, cc(1) + extent
             mi = modulo(ci-1,N_cells(1)) + 1 ; mj = modulo(cj-1,N_cells(2)) + 1
             mk = modulo(ck-1,N_cells(3)) + 1 ;
             if (par_list(0,mi,mj,mk)>0) then
                x_cell = cell_center(mi,mj,mk)
                call rel_pos(x_cell, x0, L, x)
                if ( sum(x**2) .le. ((pft % dr * pft % Nr)**2) ) then
                   do_t = .true.
                   t_count = 0
                   t = 0.d0
                end if

                v = Vcom(1:3,mi,mj,mk)/Vcom(4,mi,mj,mk)
                do i=1,par_list(0,mi,mj,mk)
                   part = par_list(i,mi,mj,mk)
                   part_s = so_species(part)
                   ! obtain x as relative to x0
                   call rel_pos(so_r(:,part) , x0 , L, x)
                   ! computation of polar position
                   r = sqrt(sum(x**2))
                   th = asin( sum(unit_r*x)/r )
                   r_idx = floor( r / pft % dr ) +1
                   th_idx = floor( (th + PI*0.5d0) / pft % dtheta ) + 1
                   if ( (r_idx .le. pft % Nr) .and. &
                        (th_idx .le. pft % Nth) ) then
                      pft % c(part_s, th_idx, r_idx) = &
                           pft % c(part_s, th_idx, r_idx) + 1
                      ! computation of base vector for the velocity
                      sr = x / r
                      sp = cross_product(x,unit_r)
                      norm = sqrt(sum(sp**2))
                      sp = sp / norm
                      sth = cross_product(sp,sr)
                      v_pol(1) = sum(sr*so_v(:,part))
                      v_pol(2) = sum(sth*so_v(:,part))
                      v_pol(3) = sum(sp*so_v(:,part))
                      pft % v(:,part_s,th_idx,r_idx) = &
                           pft % v(:,part_s,th_idx,r_idx) + v_pol
                      pft % v_count(part_s,th_idx,r_idx) = &
                           pft % v_count(part_s,th_idx,r_idx) + 1
                      if (do_t) then
                         t = t + so_sys % mass(part_s) * sum( (so_v(:,part)-v)**2 )
                      end if
                   end if
                end do
                if (par_list(0,mi,mj,mk) > 1) then
                   x_cell = cell_center(mi,mj,mk)
                   call rel_pos(x_cell, x0, L, x)
                   r = sqrt(sum(x**2))
                   th = asin( sum(unit_r*x)/r )
                   r_idx = floor( r / pft % dr ) +1
                   th_idx = floor( (th + PI*0.5d0) / pft % dtheta ) + 1
                   if ( (r_idx .le. pft % Nr) .and. &
                        (th_idx .le. pft % Nth) ) then
                      pft % t(th_idx, r_idx) = pft % t(th_idx, r_idx) + &
                           t / (3.d0*dble(par_list(0,mi,mj,mk)-1) )
                      pft % t_count(th_idx, r_idx) = pft % t_count(th_idx, r_idx) + 1
                   end if
                end if
             end if
          end do
       end do
    end do

  end subroutine update_polar_fields

  !> Writes the collected data to H5MD data groups, one group per field is
  !! used. The subroutine also resets the distribution to allow further binning.
  !! @param pft polar_fields_t variable
  !! @param step The integer time.
  !! @param time The time.
  subroutine write_polar_fields(pft, step, time)
    type(polar_fields_t), intent(inout) :: pft
    integer, intent(in) :: step
    double precision, intent(in) :: time

    integer :: r_idx, th_idx, s_idx
    double precision :: r, theta

    pft % c = pft % c / (dble(pft % count) * 2.d0 * PI)
    do r_idx = 1, pft % Nr
       do th_idx = 1, pft % Nth
          r = (dble(r_idx)-0.5d0) * pft % dr
          theta = -PI*0.5d0 + (dble(th_idx)-0.5d0) * pft % dtheta
          pft % c(:,th_idx,r_idx) = pft % c(:,th_idx,r_idx) / &
               (r * pft % dtheta * pft % dr * r * cos(theta))
       end do
    end do

    pft % v_count = max(1, pft % v_count)
    forall (r_idx=1:pft % Nr, th_idx=1:pft % Nth, s_idx=1:pft % Ns)
       pft % v(:,s_idx,th_idx,r_idx) = pft % v(:,s_idx,th_idx,r_idx) / dble(pft % v_count(s_idx,th_idx,r_idx))
    end forall
    pft % t_count = max(1, pft % t_count)
    forall (r_idx=1:pft % Nr, th_idx=1:pft % Nth)
       pft % t(th_idx,r_idx) = pft % t(th_idx,r_idx) / dble(pft % t_count(th_idx,r_idx))
    end forall

    ! write the fields
    call h5md_write_obs(pft % c_ID, pft % c, step, time)
    call h5md_write_obs(pft % v_ID, pft % v, step, time)
    call h5md_write_obs(pft % t_ID, pft % t, step, time)

    ! reset fields and counters
    pft % c = 0.d0
    pft % v = 0.d0
    pft % v_count = 0
    pft % t = 0.d0
    pft % t_count = 0
    pft % count = 0

  end subroutine write_polar_fields

  !> Computes the cross product of two vectors in 3D cartesian coordinates.
  !! @param a First input vector.
  !! @param b Second input vector.
  !! @return cross_product The cross product of a and b.
  function cross_product(a,b)
    implicit none
    double precision, intent(in), dimension(3) :: a,b
    double precision, dimension(3) :: cross_product

    cross_product(1) = a(2)*b(3)-a(3)*b(2)
    cross_product(2) = a(3)*b(1)-a(1)*b(3)
    cross_product(3) = a(1)*b(2)-a(2)*b(1)

  end function cross_product

end module polar_fields
  
