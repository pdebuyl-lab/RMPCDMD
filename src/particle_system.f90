module particle_system
  use common
  use interaction
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
     double precision, pointer :: force3(:,:)
     double precision, pointer :: force4(:,:)
     double precision, pointer :: force(:,:)
     double precision, pointer :: force_store(:,:)
     double precision, pointer :: force_old(:,:)
     double precision, pointer :: force_old_store(:,:)
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
     integer, pointer :: image1(:,:)
     integer, pointer :: image2(:,:)
     integer, pointer :: image(:,:)
     integer, pointer :: image_old(:,:)
     integer, pointer :: image_pointer(:,:)
     integer, pointer :: flag1(:)
     integer, pointer :: flag2(:)
     integer, pointer :: flag(:)
     integer, pointer :: flag_old(:)
     integer, pointer :: flag_pointer(:)
   contains
     procedure :: init
     procedure :: init_from_file
     procedure :: random_placement
     procedure :: sort
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
    allocate(this% force3(3, Nmax))
    allocate(this% force4(3, Nmax))
    this% force => this% force1
    this% force_store => this% force2
    this% force_old => this% force3
    this% force_old_store => this% force4

    allocate(this% id1(Nmax))
    allocate(this% id2(Nmax))
    this% id => this% id1
    this% id_old => this% id2

    allocate(this% species1(Nmax))
    allocate(this% species2(Nmax))
    this% species => this% species1
    this% species_old => this% species2

    allocate(this% image1(3, Nmax))
    allocate(this% image2(3, Nmax))
    this% image => this% image1
    this% image_old => this% image2

    this% image = 0
    
    allocate(this% flag1(Nmax))
    allocate(this% flag2(Nmax))
    this% flag => this% flag1
    this% flag_old => this% flag2
    
    this% flag = 0

  end subroutine init

  subroutine init_from_file(this, filename, group_name, mode, idx)
    use h5md_module
    use hdf5
    implicit none
    class(particle_system_t), intent(out) :: this
    character(len=*), intent(in) :: filename, group_name
    integer, intent(in) :: mode
    integer, intent(in), optional :: idx

    type(h5md_file_t) :: f
    type(h5md_element_t) :: e
    integer(HID_T) :: g, mem_space, file_space
    integer(HSIZE_T) :: dims(3), maxdims(3), start(3)
    integer :: rank
    double precision, allocatable :: pos(:,:)

    call f% open(filename, H5F_ACC_RDONLY_F)

    call h5md_check_valid(f% particles, 'particles group not found')

    call h5md_check_exists(f% particles, group_name, 'group '//group_name//' not found')

    call h5gopen_f(f% particles, group_name, g, f% error)

    if (mode == H5MD_FIXED) then
       call e% read_fixed(g, 'position', pos)
       call this% init(size(pos, 2))
       this% pos = pos
       deallocate(pos)
    else if ( (mode == H5MD_TIME) .or. (mode == H5MD_LINEAR) ) then
       if (.not. present(idx) ) error stop 'missing idx in init_from_file'
       call e% open_time(g, 'position')
       call this% init(e% Nmax)
       dims = [3, e% Nmax, 1]
       call h5screate_simple_f(3, dims, mem_space, e% error)

       call h5dget_space_f(e% v, file_space, e% error)
       call h5sget_simple_extent_ndims_f(file_space, rank, e% error)
       if (rank /= 3) then
          error stop 'invalid dataset rank'
       end if
       call h5sget_simple_extent_dims_f(file_space, dims, maxdims, e% error)
       start = [0, 0, idx]
       dims(3) = 1
       call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, start, dims, e% error)
       call h5dread_f(e% v, H5T_NATIVE_DOUBLE, this% pos, dims, e% error, mem_space, file_space)
       call h5sclose_f(file_space, e% error)
       call h5sclose_f(mem_space, e% error)
       call h5dclose_f(e% v, e% error)
    end if

    call h5gclose_f(g, f% error)
    call f% close()

  end subroutine init_from_file

  subroutine random_placement(this, L, other, lj_params)
    class(particle_system_t), intent(inout) :: this
    double precision, intent(in) :: L(3)
    type(particle_system_t), intent(inout), optional :: other
    type(lj_params_t), intent(in), optional :: lj_params

    integer :: i, j, N_obstacles
    integer :: s1, s2
    double precision :: x(3), rsq
    logical :: tooclose

    if (present(other) .or. present(lj_params)) then
       if ( .not. ( present(other) .and. present(lj_params) ) ) then
          error stop 'other and lj_params must be present/absent together'
       end if
       N_obstacles = other% Nmax
    else
       N_obstacles = 0
    end if

    if (N_obstacles .gt. 0) then
       do i = 1, this% Nmax
          s1 = this% species(i)
          tooclose = .true.
          do while ( tooclose )
             call random_number(x)
             x = x * L
             tooclose = .false.
             do j = 1, N_obstacles
                s2 = other% species(j)
                rsq = sum(rel_pos(x, other% pos(:, j), L)**2)
                if ( rsq < lj_params% cut_sq(s1, s2) ) then
                   tooclose = .true.
                   exit
                end if
             end do
          end do
          this% pos(:, i) = x
       end do
    else
       call random_number(this% pos(:, :))
       do j=1, 3
          this% pos(j, :) = this% pos(j, :) * L(j)
       end do
    end if

  end subroutine random_placement

  subroutine sort(this, cells)
    use cell_system
    use hilbert
    implicit none
    class(particle_system_t), intent(inout) :: this
    type(cell_system_t), intent(inout) :: cells

    integer :: i, idx, start, p(3)
    integer :: L(3)
    logical :: nowalls

    L = cells% L
    nowalls = .not. cells% has_walls

    call cells%count_particles(this% pos)

    !$omp parallel do private(p, idx, start)
    do i=1, this% Nmax
       if (this% species(i) == 0) continue
       p = floor( (this% pos(:, i) / cells% a ) - cells% origin )
       if ( p(1) == L(1) ) p(1) = 0
       if ( p(2) == L(2) ) p(2) = 0
       if (nowalls) then
          if ( p(3) == L(3) ) p(3) = 0
       end if
       idx = compact_p_to_h(p, cells% M) + 1
       !$omp atomic capture
       start = cells% cell_start(idx)
       cells% cell_start(idx) = cells% cell_start(idx) + 1
       !$omp end atomic
       this% pos_old(:, start) = this% pos(:, i)
       this% vel_old(:, start) = this% vel(:, i)
       this% force_store(:, start) = this% force(:, i)
       this% force_old_store(:, start) = this% force_old(:, i)
       this% id_old(start) = this% id(i)
       this% species_old(start) = this% species(i)
    end do

    cells% cell_start = cells% cell_start - cells% cell_count

    call switch(this% pos, this% pos_old)
    call switch(this% vel, this% vel_old)
    call switch(this% force, this% force_store)
    call switch(this% force_old, this% force_old_store)
    call switch(this% id, this% id_old)
    call switch(this% species, this% species_old)

  end subroutine sort


end module particle_system
