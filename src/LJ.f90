
!> The module LJ serves to describe the Lennard-Jones interaction potential and
!! provides routines for the computation of the forces and potential energy, with
!! and without smoothing.
!!
!! The potential can be set with a cut-off of \f$ \sigma\ 2^{1/6} \f$ or \f$ 2.5 \sigma \f$.
!! The potential is always shifted by its value at the cut-off to avoid discontinuities.
!! The function \f$ \frac{x^4}{1+x^4} \f$ can be applied to smooth the potential.

module LJ
  
  !> A LJdata variable contains all the parameters for Lennard-Jones interaction
  !! between two ensembles of particles. Each ensemble may possess several 
  !! species. The interaction may be smoothed or not.
  type LJdata
     !> Number of particles in ensemble a.
     integer :: Na
     !> Number of particles in ensemble a.
     integer :: Nb
     !> Epsilon parameter, dimensions (Na,Nb).
     double precision, allocatable :: eps(:,:)
     !> Sigma parameter, dimensions (Na,Nb).
     double precision, allocatable :: sig(:,:)
     !> Cut-off radius, dimensions (Na,Nb).
     double precision, allocatable :: cut(:,:)
     !> Value of the potential at cut-off radius, dimensions (Na,Nb).
     double precision, allocatable :: V_c(:,:)
     !> Neighbouring radius (cut-off + skin), dimensions (Na,Nb).
     double precision, allocatable :: neigh(:,:)
     !> Switch for the smoothing of the potential, dimensions (Na,Nb).
     logical, allocatable :: smooth(:,:)
  end type LJdata

  !> The inter-atoms LJ parameters.
  type(LJdata) :: at_at
  !> The atom-solvent LJ parameters.
  type(LJdata) :: at_so

contains

  !> Configures LJ parameters from a configuration file.
  !>
  !> @param CF the configuration file.
  !> @param Nat the number of atom species.
  !> @param Nso The number of solvent species.
  subroutine config_LJdata(CF, Nat, Nso)
    use ParseText
    implicit none
    type(PTo), intent(in) :: CF
    integer, intent(in) :: Nat, Nso

    integer :: i, j
    character(len=2) LJindex
    double precision :: skin_factor
    integer, allocatable :: cut(:,:)

    skin_factor = PTread_d(CF, 'LJskinfactor')

    at_at%Na = Nat
    at_at%Nb = Nat

    allocate(at_at%eps(at_at%Na,at_at%Nb))
    allocate(at_at%sig(at_at%Na,at_at%Nb))
    allocate(at_at%cut(at_at%Na,at_at%Nb))
    allocate(at_at%neigh(at_at%Na,at_at%Nb))
    allocate(at_at%smooth(at_at%Na,at_at%Nb))
    allocate(at_at%V_c(at_at%Na,at_at%Nb))

    allocate(cut(at_at%Na,at_at%Nb))

    do i=1,at_at%Na
       write(LJindex,'(i02.2)') i
       at_at%eps(i,:) = PTread_dvec(CF,'at_ateps'//LJindex,at_at%Nb)
       at_at%sig(i,:) = PTread_dvec(CF,'at_atsig'//LJindex,at_at%Nb)
       at_at%smooth(i,:) = PTread_lvec(CF,'at_atsmooth'//LJindex,at_at%Nb)
       cut(i,:) = PTread_ivec(CF,'at_atcut'//LJindex,at_at%Nb)
    end do

    where (cut .eq.1) 
       at_at%cut = at_at%sig * 2.d0**(1.d0/6.d0)
       at_at%V_c = -at_at%eps
    end where
    forall (i=1:at_at%Na, j=1:at_at%Nb, cut(i,j).eq.2)
       at_at%cut(i,j) = at_at%sig(i,j) * 2.5d0
       at_at%V_c(i,j) = LJ_V(at_at%eps(i,j),at_at%sig(i,j), 0.d0,(at_at%sig(i,j) * 2.5d0)**2)
    end forall

    deallocate(cut)

    at_at%neigh = at_at%cut * skin_factor

    at_so%Na = Nat
    at_so%Nb = Nso

    allocate(at_so%eps(at_so%Na,at_so%Nb))
    allocate(at_so%sig(at_so%Na,at_so%Nb))
    allocate(at_so%cut(at_so%Na,at_so%Nb))
    allocate(at_so%neigh(at_so%Na,at_so%Nb))
    allocate(at_so%smooth(at_so%Na,at_so%Nb))
    allocate(at_so%V_c(at_so%Na,at_so%Nb))

    allocate(cut(at_so%Na,at_so%Nb))

    do i=1,at_so%Na
       write(LJindex,'(i02.2)') i
       at_so%eps(i,:) = PTread_dvec(CF,'at_soeps'//LJindex,at_so%Nb)
       at_so%sig(i,:) = PTread_dvec(CF,'at_sosig'//LJindex,at_so%Nb)
       at_so % smooth(i,:) = PTread_lvec(CF,'at_sosmooth'//LJindex,at_so%Nb)
       cut(i,:) = PTread_ivec(CF,'at_socut'//LJindex,at_so%Nb)
    end do

    where (cut .eq.1)
       at_so%cut = at_so%sig * 2.d0**(1.d0/6.d0)
       at_so%V_c = -at_so%eps
    end where
    forall (i=1:at_so%Na, j=1:at_so%Nb, cut(i,j).eq.2)
       at_so%cut(i,j) = at_so%sig(i,j) * 2.5d0
       at_so%V_c(i,j) = LJ_V(at_so%eps(i,j),at_so%sig(i,j), 0.d0,(at_so%sig(i,j) * 2.5d0)**2)
    end forall

    deallocate(cut)    

    at_so%neigh = at_so%cut * skin_factor
    
  end subroutine config_LJdata

  !> Returns the magnitude of the LJ force over r.
  !>
  !> @param eps The LJ epsilon parameter.
  !> @param sigma The LJ sigma parameter.
  !> @param rsq The squared distance between particles.
  !  ! The corresponding potential is
  ! \begin{equation}
  ! V_{LJ}(r) = 4 \epsilon \left( \frac{\sigma^{12}}{r^{12}} - \frac{\sigma^{6}}{r^{6}} + \frac{1}{4} \right)
  ! \end{equation}
  ! The magnitude of the force is
  ! \begin{equation}
  ! F = 4 \epsilon \left( 12 \frac{\sigma^{12}}{r^{13}} - 6 \frac{\sigma^{6}}{r^{7}}  \right)
  ! \end{equation}
  pure function LJ_force_or(eps, sigma, rsq)
    implicit none
    double precision :: LJ_force_or
    double precision, intent(in) :: eps, sigma, rsq

    double precision :: sig6_o_r6

    sig6_o_r6 = sigma**6/rsq**3

    LJ_force_or = 24.d0*eps* sig6_o_r6/rsq * (2.d0*sig6_o_r6 - 1.d0)

  end function LJ_force_or

  !> Returns the magnitude of the LJ force over r, with smoothing.
  !>
  !> @param eps The LJ epsilon parameter.
  !> @param sigma The LJ sigma parameter.
  !> @param V_c The shift at the cut-off.
  !> @param rsq The squared distance between particles.
  !> @param rcut The cut-off radius.
  !> @param h The smoothing parameter.
  ! returns the magnitude over r of the smoothed LJ force.
  ! The smoothed potential is V_LJ * x**4/(1.+x**4) where x = (r-rcut)/(sigma*h) 
  ! where V_LJ is shifted to correct the discontinuity at cut-off.
  ! See P. H. Colberg & F. Hofling, Comp. Phys. Comm. 182, pp 1120-1129 (2011)
  ! 
  pure function LJ_force_smooth_or(eps, sigma, V_c, rsq, rcut, h)
    implicit none
    double precision :: LJ_force_smooth_or
    double precision, intent(in) :: eps, sigma, V_c, rsq, rcut, h

    double precision :: sig6_o_r6, x, x4

    sig6_o_r6 = sigma**6/rsq**3

    x = (rcut - sqrt(rsq))/(sigma*h)
    x4 = x**4

    LJ_force_smooth_or = 24.d0*eps* sig6_o_r6/rsq * (2.d0*sig6_o_r6 - 1.d0) * x4/(1.d0+x4)

    LJ_force_smooth_or = LJ_force_smooth_or + &
         4.d0 * (4.d0 * eps * sig6_o_r6 * ( sig6_o_r6 - 1.d0 ) - V_c)/sqrt(rsq) * x**3 / ( (1.d0 + x4)**2 * sigma * h)

  end function LJ_force_smooth_or

  !> Returns the LJ potential of two particles.
  !>
  !> @param eps The LJ epsilon parameter.
  !> @param sigma The LJ sigma parameter.
  !> @param V_c The shift at the cut-off.
  !> @param rsq The squared distance between particles.
  pure function LJ_V(eps, sigma, V_c, rsq)
    implicit none
    double precision :: LJ_V
    double precision, intent(in) :: eps, sigma, V_c, rsq

    double precision :: sig6_o_r6

    sig6_o_r6 = sigma**6/rsq**3
    
    LJ_V = 4.d0 * eps * (sig6_o_r6 * ( sig6_o_r6 - 1.d0 )) - V_c

  end function LJ_V

  !> Returns the LJ potential of two particles, with smoothing.
  !>
  !> @param eps The LJ epsilon parameter.
  !> @param sigma The LJ sigma parameter.
  !> @param V_c The shift at the cut-off.
  !> @param rsq The squared distance between particles.
  !> @param rcut The cut-off radius.
  !> @param h The smoothing parameter.
  pure function LJ_V_smooth(eps, sigma, V_c, rsq, rcut, h)
    implicit none
    double precision :: LJ_V_smooth
    double precision, intent(in) :: eps, sigma, V_c, rsq, rcut, h

    double precision :: sig6_o_r6, x, x4

    x = (rcut - sqrt(rsq))/(sigma*h)
    x4 = x**4

    sig6_o_r6 = sigma**6/rsq**3

    LJ_V_smooth = (eps * 4.d0 * sig6_o_r6 * ( sig6_o_r6 - 1.d0 ) - V_c )* x4/(1.d0+x4)
    
  end function LJ_V_smooth

end module LJ
