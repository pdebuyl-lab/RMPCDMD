module group
  
  type group_t
     integer :: GID
     integer :: N
     integer :: istart
     integer :: g_type
     integer :: species1, species2
     double precision :: dimer_length
     double precision :: elast_k, elast_rmax
     integer :: elast_nlink
     integer, allocatable :: elast_index(:,:)
     double precision, allocatable :: elast_r0(:)
  end type group_t

  integer, parameter :: ATOM_G=1
  integer, parameter :: DIMER_G=2
  integer, parameter :: ELAST_G=3

contains

  subroutine config_group(group_var, group_i, istart, CF)
    use ParseText
    implicit none
    type(group_t), intent(out) :: group_var
    integer, intent(in) :: group_i, istart
    type(PTo), intent(in) :: CF

    character(len=2) :: group_index
    character(len=16) :: g_string
    integer :: s(2)

    write(group_index,'(i2.02)') group_i

    g_string = PTread_s(CF,'group'//group_index)
    select case (g_string)
    case('atom')
       group_var % g_type = ATOM_G
       group_var % N = PTread_i(CF,'group'//group_index//'N')
       group_var % species1 = PTread_i(CF,'group'//group_index//'species')
    case('dimer')
       group_var % g_type = DIMER_G
       group_var % N = 2
       group_var % dimer_length = PTread_d(CF,'group'//group_index//'length')
       s = PTread_ivec(CF,'group'//group_index//'species',2)
       group_var % species1 = s(1)
       group_var % species2 = s(2)
    case('elast')
       group_var % g_type = ELAST_G
       group_var % N = PTread_i(CF,'group'//group_index//'N')
       group_var % elast_k = PTread_d(CF,'group'//group_index//'k')
       group_var % elast_rmax = PTread_d(CF,'group'//group_index//'rmax')
       group_var % species1 = PTread_i(CF,'group'//group_index//'species')
    case default
       write(*,*) 'unknown type for', g_string
       stop
    end select

    group_var % istart = istart

  end subroutine config_group


end module group
