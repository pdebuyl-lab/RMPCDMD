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
   contains
     procedure :: init
     procedure :: del
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

  end subroutine init

  subroutine del(this)
    class(particle_system_t), intent(inout) :: this

    deallocate(this% pos1)
    deallocate(this% pos2)
    this% pos => null()
    this% pos_old => null()
    this% pos_pointer => null()

  end subroutine del

end module particle_system
