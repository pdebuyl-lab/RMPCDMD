
!> The module polar_dist defines a type for the storage of symmetric polar
!! distributions and the associated subroutines.

module polar_dist
  use MD, only: rel_pos
  use MPCD, only: L
  implicit none

  !> A polar_dist_t variable contains the parameters to compute a polar
  !! distribution function as well as the polar distribution function itself.
  type polar_dist_t
     !> Data for the polar distribution function.
     integer, allocatable :: g(:,:,:)
     !> Number of times that the binning has been performed.
     integer :: t_count
     !> Radial size of the bins.
     double precision :: dr
     !> Angular size of the bins.
     double precision :: dtheta
  end type polar_dist_t

contains

  !> Defines the parameters of a polar_dist_t variable and allocates the 
  !! distribution array.
  !! 
  !! @param gor The polar_dist_t variable.
  !! @param N_s The number of species to take into account.
  !! @param N_r The number of bins of the distribution.
  !! @param dr The radial bin size.
  !! @param N_theta The number of angular bins.
  subroutine init_polar(gor, N_s, N_r, dr, N_theta)
    implicit none
    type(polar_dist_t), intent(out) :: gor
    integer, intent(in) :: N_s,N_r
    double precision, intent(in) :: dr
    integer, intent(in) :: N_theta

    allocate(gor%g(N_s,N_theta,N_r))
    gor % g = 0 
    gor % t_count = 0
    gor % dr = dr
    gor % dtheta = 4.d0*atan(1.d0)/dble(N_theta)

  end subroutine init_polar

  !> Bins the particles from the array positions with indices given by list,
  !! taking x_0 as the center of the coordinates and dir_in as the axis.
  !! The counter t_count of the polar_dist_t variable is increased by 1.
  !! @param gor The polar_dist_t variable.
  !! @param x_0 The center of coordinates for the binning.
  !! @param dir_in The direction along which theta=0 for the binning.
  !! @param positions Array of positions.
  !! @param species Array of species.
  !! @param list The indices of particles in positions to consider.
  subroutine update_polar(gor, x_0, dir_in, positions, species, list)
    implicit none
    type(polar_dist_t), intent(inout) :: gor
    double precision, intent(in) :: x_0(3), dir_in(3)
    double precision, intent(in) :: positions(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: list(:)

    integer :: i, r_idx, th_idx, N, part
    double precision :: x(3), r, dir(3)
    
    N = size(list)

    dir = dir_in / sqrt( dir_in(1)**2 + dir_in(2)**2 + dir_in(3)**2 )
    
    do i=1, N
       part = list(i)
       call rel_pos(x_0, positions(:,part), L, x)
       r = sqrt( sum( x**2 ) )
       r_idx  = floor(r / gor % dr) + 1
       th_idx = floor(acos( sum(x*dir)/r )/gor % dtheta)+1
       if (r_idx .le. size(gor % g, dim=3) .and. th_idx .le. size(gor % g, dim=2)) &
            gor % g(species(part),th_idx, r_idx) = gor % g(species(part),th_idx, r_idx) + 1
    end do
  
    gor % t_count = gor % t_count + 1
  
  end subroutine update_polar

  !> Dumps the distribution in a HDF5 file, following the H5MD organization.
  !! @param gor The polar_dist_t variable.
  !! @param fileID The ID of the HDF5 file.
  !! @param group_name The name of the trajectory group.
  subroutine write_polar(gor, fileID, group_name)
    use HDF5
    implicit none
    type(polar_dist_t), intent(in) :: gor
    integer(HID_T), intent(inout) :: fileID
    character(len=*), intent(in), optional :: group_name
  
    double precision, allocatable :: g_real(:,:,:)
    integer :: i,j
    double precision :: r, theta
    integer(HID_T) :: d_id, s_id, a_id
    integer(HSIZE_T) :: dims(3)
    character(len=128) :: path
    integer :: h5_error
  
    allocate( g_real( size(gor % g, dim=1), size(gor % g, dim=2), size(gor % g, dim=3) ) )
    g_real = dble(gor % g) / ( dble(gor % t_count) * 2.d0 * 4.d0*atan(1.d0) )
    do i=1,size(g_real,dim=3)
       do j=1,size(g_real,dim=2)
          r = (dble(i)-0.5d0) * gor % dr
          theta = (dble(j)-0.5d0) * gor % dtheta
          g_real(:,j,i) = g_real(:,j,i) / ( r*gor % dtheta * gor % dr * r * sin(theta) )
       end do
    end do
    ! open dataset
    dims(1) = size(g_real, dim=1)
    dims(2) = size(g_real, dim=2)
    dims(3) = size(g_real, dim=3)
    call h5screate_simple_f(3, dims, s_id, h5_error)
    if (present(group_name)) then
       path = 'trajectory/'//group_name//'/polar_distribution'
    else
       path = 'trajectory/radial_distribution'
    end if
    call h5dcreate_f(fileID, path, H5T_NATIVE_DOUBLE, s_id, d_id, h5_error)
    call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, g_real, dims, h5_error)
    call h5sclose_f(s_id, h5_error)
    deallocate( g_real )
    ! write attributes: dr
    call h5screate_f(H5S_SCALAR_F, s_id, h5_error)
    call h5acreate_f(d_id, 'dr', H5T_NATIVE_DOUBLE, s_id, a_id, h5_error)
    call h5awrite_f(a_id, H5T_NATIVE_DOUBLE, gor % dr, dims, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(s_id, h5_error)

    call h5screate_f(H5S_SCALAR_F, s_id, h5_error)
    call h5acreate_f(d_id, 'dtheta', H5T_NATIVE_DOUBLE, s_id, a_id, h5_error)
    call h5awrite_f(a_id, H5T_NATIVE_DOUBLE, gor % dtheta, dims, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(s_id, h5_error)

    call h5dclose_f(d_id, h5_error)

  end subroutine write_polar

end module polar_dist
