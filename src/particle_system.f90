module particle_system
  use common
  implicit none
  private

  public :: particle_system_t

  type particle_system_t
     integer :: Nmax
     integer :: n_species
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
     procedure :: init_from_file
     procedure :: random_placement
  end type particle_system_t

contains

  subroutine init(this, Nmax, n_species)
    class(particle_system_t), intent(out) :: this
    integer, intent(in) :: Nmax
    integer, intent(in), optional :: n_species

    this% Nmax = Nmax
    if (present(n_species)) then
       this% n_species = n_species
    else
       this% n_species = 1
    end if

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

  subroutine init_from_file(this, filename, group_name)
    use h5md_module
    use hdf5
    implicit none
    class(particle_system_t), intent(out) :: this
    character(len=*), intent(in) :: filename, group_name

    type(h5md_file_t) :: f
    type(h5md_element_t) :: e
    logical :: valid, link_exists
    integer(HID_T) :: g, mem_space, file_space
    integer(HSIZE_T) :: dims(3), maxdims(3), start(3)
    integer :: rank

    call f% open(filename, H5F_ACC_RDONLY_F)

    call h5iis_valid_f(f% particles, valid, f% error)
    if (.not. valid) then
       stop 'particles group not found'
    end if

    call h5lexists_f(f% particles, group_name, link_exists, e% error)
    if (.not. link_exists) then
       write(*,*) 'particles group ', group_name, ' not found'
       stop
    end if
    call h5gopen_f(f% particles, group_name, g, f% error)

    call e% open_time(g, 'position')
    call this% init(e% Nmax)

    dims = [3, e% Nmax, 1]
    call h5screate_simple_f(3, dims, mem_space, e% error)

    call h5dget_space_f(e% v, file_space, e% error)
    call h5sget_simple_extent_ndims_f(file_space, rank, e% error)
    if (rank /= 3) then
       stop 'invalid dataset rank'
    end if
    call h5sget_simple_extent_dims_f(file_space, dims, maxdims, e% error)
    start = [0, 0, 0]
    dims(3) = 1
    call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, start, dims, e% error)
    call h5dread_f(e% v, H5T_NATIVE_DOUBLE, this% pos, dims, e% error, mem_space, file_space)
    call h5sclose_f(file_space, e% error)
    call h5sclose_f(mem_space, e% error)
    call h5dclose_f(e% v, e% error)

    call f% close()

  end subroutine init_from_file

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

end module particle_system
