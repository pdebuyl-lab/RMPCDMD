module correlator
  implicit none

  private

  public :: correlator_t
  public :: correlate_block_distsq
  public :: correlate_block_dot

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

end module correlator
