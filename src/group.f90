
!> This module allows the description of colloid groups of the type "ATOM", "DIMER" and "ELAST".

module group
  
  !> Contains information about an atomic group.
  type group_t
     !> Global group ID. This needs to be set to pass around groups between CPUs.
     integer :: GID
     !> Number of atoms in the group.
     integer :: N
     !> First index at which the group is found in the "at_*" arrays.
     integer :: istart
     !> Type of the group. Needs to be set to one of the hardcoded values of this module.
     integer :: g_type
     !> Number of subgroups
     integer :: N_sub
     !> subgroup array: holds the starting index in (1,:) and the ending index of a subgroup
     !! in (2,:). Dimensions (2,N_sub) where N_sub is the number of subgroups.
     integer, allocatable :: subgroup(:,:)
     !> Species holder 1. For value of first atom in a dimer and/or configuration purposes.
     integer :: species1
     !> Species holder 2. For value of second atom in a dimer and/or configuration purposes.
     integer :: species2
     !> In the case of a dimer, length of the dimer.
     double precision :: dimer_length
     !> In the case of an elastic network, stiffness of the network.
     double precision :: elast_k
  end type group_t

  !> g_type parameter for unbound LJ atoms.
  integer, parameter :: ATOM_G=1
  !> g_type parameter for a dimer.
  integer, parameter :: DIMER_G=2
  !> g_type parameter for an elastic network.
  integer, parameter :: ELAST_G=3

contains

  !> Sets a group variable from a configuration file.
  !>
  !> @param group_var Variable that is set from the file.
  !> @param group_i Index of the group, used to parse the file.
  !> @param istart Index of the first atom in the "at_*" arrays
  !> @param CF Configuration file.
  subroutine config_group(group_var, group_i, istart, CF)
    use ParseText
    implicit none
    type(group_t), intent(out) :: group_var
    integer, intent(in) :: group_i, istart
    type(PTo), intent(in) :: CF

    character(len=2) :: group_index
    character(len=16) :: g_string
    integer :: s(2), i0, i
    integer, allocatable :: subg_vec(:)

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
    case default
       write(*,*) 'unknown type for', g_string
       stop
    end select

    group_var % N_sub = PTread_i(CF,'group'//group_index//'Nsub')
    if (group_var % N_sub > 0) then
       allocate( subg_vec(group_var % N_sub) )
       allocate( group_var % subgroup(2,group_var % N_sub) )
       subg_vec = PTread_ivec(CF,'group'//group_index//'sub', size(subg_vec))
       i0 = 1
       do i=1,group_var % N_sub
          group_var % subgroup(1,i) = i0
          i0 = i0 + subg_vec(i)
          group_var % subgroup(2,i) = i0 - 1
       end do
       deallocate( subg_vec )
    end if
    
    group_var % istart = istart

  end subroutine config_group


end module group
