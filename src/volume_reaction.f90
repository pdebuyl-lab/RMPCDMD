
!> This module defines the volume chemical reaction for a MPCD system.
module volume_reaction
  implicit none

  !> Definition of a bulk chemical reaction for a MPCD system.
  type vol_reac_t
     !> The kind of chemical reaction.
     integer :: r_kind
     !> The rate of the reaction.
     double precision :: rate
     !> The number of reacting particles.
     integer :: N_reac
     !> The number of product particles.
     integer :: N_prod
     !> The reacting species. Dimensions: 1 (N_reac). Ordered
     !! according to the reaction scheme.
     integer, allocatable :: reac(:)
     !> Same as reac for product species.
     integer, allocatable :: prod(:)
     !> The reacting species given explicitly. Dimensions: 1 (so_sys % N_species).
     integer, allocatable :: species_reac(:)
  end type vol_reac_t

  !> r_kind parameter for unimolecular conversion reaction.
  integer, parameter :: AtoB_R=1

contains

  !> This routine parses the reaction configuration and sets up a 
  !! vol_reac_t variable.
  !! @param CF The configuration file.
  !! @param vol_reac_var The vol_reac_t variable.
  !! @param reac_i The index of the reaction.
  !! @param N_species The number of solvent species.
  subroutine add_vol_reac(CF, vol_reac_var, reac_i, N_species)
    use ParseText
    implicit none
    type(PTo), intent(in) :: CF
    type(vol_reac_t), intent(out) :: vol_reac_var
    integer, intent(in) :: reac_i, N_species

    character(len=10) :: reac_s
    character(len=16) :: r_kind_string

    write(reac_s, '(a8,i2.02)') 'vol_reac', reac_i

    r_kind_string = PTread_s(CF, reac_s//'kind')
    select case(r_kind_string)
    case('AtoB')
       vol_reac_var % r_kind = AtoB_R
       vol_reac_var % N_reac = 1
       vol_reac_var % N_prod = 1
       allocate( vol_reac_var % reac(vol_reac_var % N_reac) )
       vol_reac_var % reac = PTread_ivec(CF, reac_s//'reac', vol_reac_var % N_reac)
       allocate( vol_reac_var % prod(vol_reac_var % N_prod) )
       vol_reac_var % prod = PTread_ivec(CF, reac_s//'prod',vol_reac_var % N_prod)
       allocate( vol_reac_var % species_reac(N_species) )
       vol_reac_var % species_reac = 0
       vol_reac_var % species_reac( vol_reac_var % reac(1) ) = 1
    case default
       write(*,*) 'unknown reaction kind for ', reac_s
    end select
    vol_reac_var % rate = PTread_d(CF, reac_s//'rate')
    
  end subroutine add_vol_reac

end module volume_reaction
