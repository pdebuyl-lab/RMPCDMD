
!> The module radial_dist provides the rad_dist_t type that is used to define
!! the parameters and to store a radial distribution function.
!!
!! The following example illustrates the usage of the module.
!! @code
!!type(rad_dist_t) :: rd
!!call init_rad(rd, 100, 0.1d0)
!! !...
!!call update_rad(rd, x_0, pos_array, list_of_neighbour_indices)
!! !!...
!! @endcode
module radial_dist
  use MD, only: rel_pos
  use MPCD, only: L
  implicit none

  !> A rad_dist_t variable contains the parameters to compute a radial
  !! distribution function as well as the radial distribution function itself.
  type rad_dist_t
     !> Data for the radial distribution function.
     integer, allocatable :: g(:)
     !> Number of times that the binning has been performed.
     integer :: t_count
     !> Radial size of the bins.
     double precision :: dr
  end type rad_dist_t

contains
  !> Defines the parameters of a rad_dist_t variable and allocates the
  !! distribution array.
  !!
  !! @param gor The rad_dist_t variable.
  !! @param N_gor The number of bins of the distribution.
  !! @param dr The bin size.
  subroutine init_rad(gor, N_gor, dr)
    implicit none
    type(rad_dist_t), intent(out) :: gor
    integer, intent(in) :: N_gor
    double precision, intent(in) :: dr

    allocate(gor%g(N_gor))
    gor % g = 0 
    gor % t_count = 0
    gor % dr = dr

  end subroutine init_rad
  
  !> Bins the particles from the array positions with indices given by list,
  !! taking x_0 as the center of the coordinates.
  !! The counter t_count of the rad_dist_t variable is increased by 1.
  !! @param gor The rad_dist_t variable.
  !! @param x_0 The center of coordinates for the binning.
  !! @param positions Array of positions.
  !! @param list The indices of particles in positions to consider.
  subroutine update_rad(gor, x_0, positions, list)
    implicit none
    type(rad_dist_t), intent(inout) :: gor
    double precision, intent(in) :: x_0(3)
    double precision, intent(in) :: positions(:,:)
    integer, intent(in) :: list(:)

    integer :: i, idx, N
    double precision :: x(3)
    
    N = size(list)

    do i=1, N
       call rel_pos(x_0, positions(:,list(i)), L, x)
       idx = floor(sqrt( sum( x**2 ) ) / gor % dr) + 1
       if (idx .le. size(gor % g)) gor % g(idx) = gor % g(idx) + 1
    end do
  
    gor % t_count = gor % t_count + 1
  
  end subroutine update_rad
  
  !> Dumps the distribution in a HDF5 file, following the H5MD organization.
  !! @param gor The rad_dist_t variable.
  !! @param fileID The ID of the HDF5 file.
  !! @param group_name The name of the trajectory group.
  subroutine write_rad(gor, fileID, group_name)
    use HDF5
    implicit none
    type(rad_dist_t), intent(in) :: gor
    integer(HID_T), intent(inout) :: fileID
    character(len=*), intent(in), optional :: group_name
  
    double precision, allocatable :: g_real(:)
    integer :: i
    double precision :: r
    integer(HID_T) :: d_id, s_id, a_id
    integer(HSIZE_T) :: dims(1)
    character(len=128) :: path
    integer :: h5_error
  
    allocate( g_real( size(gor % g) ) )
    g_real = dble(gor % g) / ( dble(gor % t_count) * 4.d0/3.d0 * 4.d0*atan(1.d0) )
    do i=1,size(g_real)
       r = dble(i-1) * gor % dr
       g_real(i) = g_real(i) / ( (r+gor % dr)**3-r**3 )
    end do
    ! open dataset
    dims(1) = size(g_real)
    call h5screate_simple_f(1, dims, s_id, h5_error)
    if (present(group_name)) then
       path = 'trajectory/'//group_name//'/radial_distribution'
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

    call h5dclose_f(d_id, h5_error)

  end subroutine write_rad

end module radial_dist
