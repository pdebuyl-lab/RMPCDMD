module particle_system

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

end module particle_system
