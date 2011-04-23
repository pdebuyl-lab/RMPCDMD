module radial_dist
  implicit none

  type rad_dist_t
     integer, allocatable :: g(:)
     integer :: t_count
     double precision :: dr
  end type rad_dist_t

contains

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
  
  subroutine update_rad(gor, x_0, positions, list)
    implicit none
    type(rad_dist_t), intent(inout) :: gor
    double precision, intent(in) :: x_0(3)
    double precision, intent(in) :: positions(:,:)
    integer, intent(in) :: list(:)

    integer :: i, idx, N
    double precision :: r
    
    N = size(list)

    do i=1, N
       idx = floor(sqrt( sum( (positions(:,list(i)) - x_0)**2 ) ) / gor % dr) + 1
       if (idx .le. size(gor % g)) gor % g(idx) = gor % g(idx) + 1
    end do
  
    gor % t_count = gor % t_count + 1
  
  end subroutine update_rad
  
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
