
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

  !> This routine performs the A to B bulk chemical reaction.
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param vr_var The vol_reac_t variable.
  subroutine vol_A_to_B(ci,cj,ck,vr_var)
    use MPCD, only: par_list, so_species
    implicit none
    integer, intent(in) :: ci,cj,ck
    type(vol_reac_t), intent(in) :: vr_var

    integer :: i, part
    logical :: found

    found = .false.
    do i=1,par_list(0,ci,cj,ck)
       part = par_list(i,ci,cj,ck)
       if (so_species(part) .eq. vr_var % reac(1)) then
          so_species(part) = vr_var % prod(1)
          found = .true.
          exit
       end if
    end do
    
    if (not(found)) stop 'reactant not found in vol_A_to_B'

  end subroutine vol_A_to_B

  !> This function computes the combinatorial term h_\mu^\xi defined
  !! in Rohlf et al, Comp. Phys. Comm. 179, pp 132 (2008), Eq. (14).
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param vr_var The vol_reac_t variable.
  function combi(ci,cj,ck,vr_var)
    use MPCD, only: par_list, so_species, so_sys
    implicit none
    integer, intent(in) :: ci,cj,ck
    type(vol_reac_t), intent(in) :: vr_var
    integer :: combi

    integer :: i, part
    integer, allocatable :: local_n(:)
    allocate(local_n(so_sys % N_species))
    local_n = 0

    do i=1,par_list(0,ci,cj,ck)
       part = par_list(i,ci,cj,ck)
       local_n(so_species(part)) = local_n(so_species(part)) + 1
    end do
    combi = 1
    do i=1,so_sys % N_species
       if (local_n(i) < vr_var % species_reac(i)) then
          combi = 0
          deallocate(local_n)
          return
       end if
       if (vr_var % species_reac(i).eq.0) cycle
       combi = combi * partial_factorial( local_n(i) , vr_var % species_reac(i) )
    end do
    
    deallocate(local_n)

  contains
    
    !> Computes n! / (n-nu)!
    !! @param n
    !! @param nu
    function partial_factorial(n, nu)
      implicit none
      integer, intent(in) :: n, nu
      integer :: partial_factorial
      
      integer :: i

      partial_factorial = 1
      do i=n-nu, n
         partial_factorial = partial_factorial * i
      end do

    end function partial_factorial
    
  end function combi

  !> This routine performs all the steps that are needed to apply the
  !! RMPCD algorithm to a given MCPD cell.
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param n_reac The number of reactions to consider.
  !! @param reac_list An array of vol_reac_t variables.
  !! @param interval The time interval for the reaction.
  subroutine cell_reaction(ci,cj,ck,n_reac,reac_list,interval)
    use MPCD, only: par_list, so_species, so_sys, ran_state
    use mtprng, only: mtprng_rand_real1
    implicit none
    integer, intent(in) :: ci,cj,ck, n_reac
    type(vol_reac_t), intent(in) :: reac_list(:)
    double precision, intent(in) :: interval

    integer :: i, idx_reac
    double precision :: a0, p_something, x
    double precision, allocatable :: amu(:)
    allocate(amu(n_reac))

    a0 = 0.d0
    do i=1,n_reac
       amu(i) = combi(ci,cj,ck,reac_list(i))*reac_list(i)%rate
       a0 = a0 + amu(i)
    end do
    p_something = a0*interval

    if (mtprng_rand_real1(ran_state) >= p_something) then
       deallocate(amu)
       return
    end if
    amu = amu/a0
    do i=2,n_reac
       amu(i) = amu(i)+amu(i-1)
    end do
    x = mtprng_rand_real1(ran_state)
    do i=1,n_reac
       if (x<amu(i)) then
          idx_reac = i
          exit
       end if
    end do

    if (idx_reac .eq. AtoB_R) then
       call vol_A_to_B(ci,cj,ck,reac_list(idx_reac))
    else
       write(*,*) 'unknown reaction type in cell_reaction'
    end if

  deallocate(amu)

  end subroutine cell_reaction


end module volume_reaction
