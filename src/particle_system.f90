module particle_system
  implicit none
  private

  public :: particle_system_t

  type particle_system_t
     integer :: Nmax
     double precision, pointer :: pos1(:,:)
     double precision, pointer :: pos2(:,:)
     double precision, pointer :: pos(:,:)
     double precision, pointer :: pos_old(:,:)
     double precision, pointer :: pos_pointer(:,:)
     double precision, pointer :: vel1(:,:)
     double precision, pointer :: vel2(:,:)
     double precision, pointer :: vel(:,:)
     double precision, pointer :: vel_old(:,:)
     double precision, pointer :: vel_pointer(:,:)
     double precision, pointer :: force1(:,:)
     double precision, pointer :: force2(:,:)
     double precision, pointer :: force(:,:)
     double precision, pointer :: force_old(:,:)
     double precision, pointer :: force_pointer(:,:)
     integer, pointer :: id1(:)
     integer, pointer :: id2(:)
     integer, pointer :: id(:)
     integer, pointer :: id_old(:)
     integer, pointer :: id_pointer(:)
     integer, pointer :: species1(:)
     integer, pointer :: species2(:)
     integer, pointer :: species(:)
     integer, pointer :: species_old(:)
     integer, pointer :: species_pointer(:)
   contains
     procedure :: init
     procedure :: random_placement
  end type particle_system_t

contains

  subroutine init(this, Nmax)
    class(particle_system_t), intent(out) :: this
    integer, intent(in) :: Nmax

    this% Nmax = Nmax
    allocate(this% pos1(3, Nmax))
    allocate(this% pos2(3, Nmax))
    this% pos => this% pos1
    this% pos_old => this% pos2

    allocate(this% vel1(3, Nmax))
    allocate(this% vel2(3, Nmax))
    this% vel => this% vel1
    this% vel_old => this% vel2

    allocate(this% force1(3, Nmax))
    allocate(this% force2(3, Nmax))
    this% force => this% force1
    this% force_old => this% force2

    allocate(this% id1(Nmax))
    allocate(this% id2(Nmax))
    this% id => this% id1
    this% id_old => this% id2

    allocate(this% species1(Nmax))
    allocate(this% species2(Nmax))
    this% species => this% species1
    this% species_old => this% species2

  end subroutine init

  subroutine random_placement(this, L, obstacles, radius)
    class(particle_system_t), intent(inout) :: this
    double precision, intent(in) :: L(3)
    double precision, intent(in), optional :: obstacles(:,:)
    double precision, intent(in), optional :: radius

    integer :: i, j, N_obstacles
    double precision :: x(3), rsq, radiussq
    logical :: tooclose

    if (present(obstacles) .or. present(radius)) then
       if ( .not. ( present(obstacles) .and. present(radius) ) ) then
          stop 'obstacles and radius must be present/absent together'
       end if
       N_obstacles = size(obstacles, 2)
    else
       N_obstacles = 0
    end if

    if (N_obstacles .gt. 0) then
       radiussq = radius**2
       do i = 1, this% Nmax
          tooclose = .true.
          do while ( tooclose )
             call random_number(x)
             x = x * L
             tooclose = .false.
             do j = 1, N_obstacles
                rsq = sum(rel_pos(x, obstacles(:, j), L)**2)
                if ( rsq < radiussq ) then
                   tooclose = .true.
                   exit
                end if
             end do
             this% pos(:, i) = x
          end do
       end do
    else
       call random_number(this% pos(:, :))
       do j=1, 3
          this% pos(j, :) = this% pos(j, :) * L(j)
       end do
    end if

  end subroutine random_placement

  function rel_pos(x, y, L) result(r)
    double precision, intent(in) :: x(3), y(3), L(3)

    double precision :: r(3)
    integer :: i, dim

    r = x - y

    do dim=1,3
       if ( r(dim) < -0.5d0*r(dim) ) then
          r(dim) = r(dim) + L(dim)
       else if ( r(dim) > 0.5d0*L(dim) ) then
          r(dim) = r(dim) - L(dim)
       end if
    end do

  end function rel_pos

end module particle_system
