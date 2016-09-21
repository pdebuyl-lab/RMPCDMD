module correlator
  implicit none

  private

  public :: correlator_t
  public :: axial_correlator_t
  public :: correlate_block_distsq
  public :: correlate_block_dot
  public :: get_n_blocks
  public :: write_correlator_block

  type correlator_t
     integer :: block_length
     integer :: n_blocks
     integer :: dim
     integer :: n
     double precision, allocatable :: value(:,:,:,:)
     double precision, allocatable :: correlation(:,:,:,:)
     integer, allocatable :: count(:)
   contains
     procedure, public :: init
     procedure, public :: add
  end type correlator_t

  type axial_correlator_t
     type(correlator_t) :: msd
     type(correlator_t) :: oacf
     type(correlator_t) :: vacf
     type(correlator_t) :: parallel_vacf
     type(correlator_t) :: transverse_vacf
   contains
     procedure :: init => axial_init
     procedure :: add => axial_add
     procedure :: add_fast => axial_add_fast
     procedure :: write => axial_write
  end type axial_correlator_t

contains

  subroutine init(this, l, b, dim, n)
    class(correlator_t), intent(out) :: this
    integer, intent(in) :: l, b
    integer, intent(in), optional :: dim
    integer, intent(in), optional :: n

    if (present(dim)) then
       this%dim = dim
    else
       this%dim = 1
    end if

    if (present(n)) then
       this%n = n
    else
       this%n = 1
    end if

    this%block_length = l
    this%n_blocks = b

    allocate(this%value(this%dim, this%n, l, b))
    allocate(this%correlation(this%dim, this%n, l, b))
    allocate(this%count(b))

    this%value(:,:,:,:) = 0
    this%correlation(:,:,:,:) = 0
    this%count(:) = 0

  end subroutine init

  subroutine add(this, i, block_operation, x, x_n, xvec, xvec_n)
    class(correlator_t), intent(inout) :: this
    integer, intent(in) :: i
    double precision, intent(in), optional :: x
    double precision, intent(in), optional :: x_n(:)
    double precision, intent(in), optional :: xvec(:)
    double precision, intent(in), optional :: xvec_n(:,:)
    interface
       subroutine block_operation(s, c, idx, l, n, dim)
         double precision, intent(in) :: s(:,:,:)
         double precision, intent(inout) :: c(:,:,:)
         integer, intent(in) :: idx
         integer, intent(in) :: l, n, dim
       end subroutine block_operation
    end interface

    integer :: l, b
    integer :: i_block, normed_idx, x_count
    double precision, allocatable :: x_gen(:,:)

    x_count = 0
    allocate(x_gen(this%dim, this%n))
    if (present(x)) then
       x_count = x_count+1
       x_gen = x
    else if (present(x_n)) then
       x_count = x_count+1
       x_gen(1,:) = x_n
    else if (present(xvec)) then
       x_count = x_count+1
       x_gen(:,1) = xvec
    else if (present(xvec_n)) then
       x_count = x_count+1
       x_gen = xvec_n
    end if

    if (x_count/=1) stop 'multiple (or no) arguments to correlator_t add'

    l = this%block_length
    b = this%n_blocks

    block_loop: do i_block = 0, b-1
       if (modulo(i, l**i_block) /= 0) exit block_loop
       normed_idx = modulo(i / l**i_block, l)
       this%value(:, :, normed_idx+1, i_block+1) = x_gen
       if (i > l**(i_block+1)) then
          ! correlate block b, for point normed_idx
          call block_operation(this%value(:,:,:,i_block+1), &
               this%correlation(:,:,:,i_block+1), normed_idx, l, this%n, this%dim)
          this%count(i_block+1) = this%count(i_block+1) + 1
       end if
    end do block_loop

    deallocate(x_gen)

  end subroutine add

  subroutine correlate_block_dot(s, c, idx, l, n, dim)
    double precision, intent(in) :: s(:,:,:)
    double precision, intent(inout) :: c(:,:,:)
    integer, intent(in) :: idx, l, n, dim

    integer :: i, j, n_idx
    double precision :: current(dim,n)

    current = s(:,:,idx+1)

    do i = 0, l-1
       j = modulo(idx-i, l) + 1
       do n_idx = 1, n
          c(:, n_idx, j) = c(:, n_idx, j) + (current(:,n_idx)*s(:,n_idx,i+1))
       end do
    end do

  end subroutine correlate_block_dot

  subroutine correlate_block_distsq(s, c, idx, l, n, dim)
    double precision, intent(in) :: s(:,:,:)
    double precision, intent(inout) :: c(:,:,:)
    integer, intent(in) :: idx, l, n, dim

    integer :: i, j, n_idx
    double precision :: current(dim, n)

    current = s(:,:,idx+1)

    do i = 0, l-1
       j = modulo(idx-i, l) + 1
       do n_idx = 1, n
          c(:, n_idx, j) = c(:, n_idx, j) + (current(:,n_idx)-s(:,n_idx,i+1))**2
       end do
    end do

  end subroutine correlate_block_distsq

  function get_n_blocks(l, n_blocks_max, n_samples) result(n_blocks)
    integer, intent(in) :: l, n_blocks_max, n_samples
    integer :: n_blocks

    do n_blocks = 1, 7
       if (l**n_blocks >= n_samples/l) exit
    end do

  end function get_n_blocks

  subroutine axial_init(this, block_length, n_samples, n_samples_fast)
    class(axial_correlator_t), intent(inout) :: this
    integer, intent(in) :: block_length, n_samples, n_samples_fast

    integer :: n_blocks

    n_blocks = get_n_blocks(block_length, 7, n_samples)
    call this%msd%init(block_length, n_blocks, dim=3)
    call this%oacf%init(block_length, n_blocks, dim=3)

    n_blocks = get_n_blocks(block_length, 8, n_samples_fast)
    call this%vacf%init(block_length, n_blocks, dim=3)
    call this%parallel_vacf%init(block_length, n_blocks, dim=3)
    call this%transverse_vacf%init(block_length, n_blocks, dim=3)

  end subroutine axial_init

  subroutine axial_add(this, i, com_pos, unit_r)
    class(axial_correlator_t), intent(inout) :: this
    integer, intent(in) :: i
    double precision, intent(in) :: com_pos(3), unit_r(3)

    call this%msd%add(i, correlate_block_distsq, xvec=com_pos)
    call this%oacf%add(i, correlate_block_dot, xvec=unit_r)

  end subroutine axial_add

  subroutine axial_add_fast(this, i, v_com, unit_r)
    class(axial_correlator_t), intent(inout) :: this
    integer, intent(in) :: i
    double precision, intent(in) :: v_com(3), unit_r(3)

    double precision :: parallel_v(3), transverse_v(3)

    parallel_v = dot_product(v_com, unit_r) * unit_r
    transverse_v = v_com - parallel_v
    call this%vacf%add(i, correlate_block_dot, xvec=v_com)
    call this%parallel_vacf%add(i, correlate_block_dot, xvec=parallel_v)
    call this%transverse_vacf%add(i, correlate_block_dot, xvec=transverse_v)

  end subroutine axial_add_fast

  subroutine axial_write(this, correlator_group, sampling, dt, fast_sampling, fast_dt)
    use hdf5
    implicit none
    class(axial_correlator_t), intent(inout) :: this
    integer(HID_T) :: correlator_group
    integer, intent(in) :: sampling, fast_sampling
    double precision, intent(in) :: dt, fast_dt

    call write_correlator_block(correlator_group, 'mean_square_displacement', &
         this%msd, sampling, dt)
    call write_correlator_block(correlator_group, 'orientation_autocorrelation', &
         this%oacf, sampling, dt)

    call write_correlator_block(correlator_group, 'velocity_autocorrelation', &
         this%vacf, fast_sampling, fast_dt)
    call write_correlator_block(correlator_group, 'parallel_velocity_autocorrelation', &
         this%parallel_vacf, fast_sampling, fast_dt)
    call write_correlator_block(correlator_group, 'transverse_velocity_autocorrelation', &
         this%transverse_vacf, fast_sampling, fast_dt)

  end subroutine axial_write

  subroutine write_correlator_block(loc, name, c, step, time)
    use hdf5
    use h5md_module
    implicit none
    integer(HID_T), intent(inout) :: loc
    character(len=*), intent(in) :: name
    type(correlator_t), intent(in) :: c
    integer, intent(in) :: step
    double precision, intent(in) :: time

    integer(HID_T) :: group
    integer :: error

    call h5gcreate_f(loc, name, group, error)
    call h5md_write_dataset(group, 'value', c%correlation)
    call h5md_write_dataset(group, 'count', c%count)
    call h5md_write_dataset(group, 'step', step)
    call h5md_write_dataset(group, 'time', time)
    call h5gclose_f(group, error)

  end subroutine write_correlator_block

end module correlator
