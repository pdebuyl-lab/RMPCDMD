module group
  
  type group_t
     integer :: GID
     integer :: N
     integer :: istart
     integer :: g_type
     double precision :: dimer_length
     double precision :: elast_k
  end type group_t

  integer, parameter :: ATOM_G=1
  integer, parameter :: DIMER_G=2
  integer, parameter :: ELAST_G=3

contains

  subroutine config_group(group_var, group_i, CF)
    use ParseText
    implicit none
    type(group_t), intent(out) :: group_var
    integer, intent(in) :: group_i
    type(PTo), intent(in) :: CF

    character(len=2) :: group_index
    character(len=16) :: g_string

    write(group_index,'(i2.02)') group_i

    g_string = PTread_s(CF,'group'//group_index)
    select case (g_string)
    case('atom')
       group_var % g_type = ATOM_G
       group_var % N = PTread_i(CF,'group'//group_index//'N')
    case('dimer')
       group_var % g_type = DIMER_G
       group_var % N = 2
       group_var % dimer_length = PTread_d(CF,'group'//group_index//'length')
    case('elast')
       group_var % g_type = ELAST_G
       group_var % N = PTread_i(CF,'group'//group_index//'N')
       group_var % elast_k = PTread_d(CF,'group'//group_index//'k')
    case default
       write(*,*) 'unknown type for', g_string
       stop
    end select

  end subroutine config_group


end module group
