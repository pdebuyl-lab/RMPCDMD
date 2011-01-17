module LJ
  
  type LJdata
     integer :: Na, Nb
     double precision, allocatable :: eps(:,:)
     double precision, allocatable :: sig(:,:)
     double precision, allocatable :: cut(:,:)
     double precision, allocatable :: neigh(:,:)
  end type LJdata

  type(LJdata) :: at_at, at_so

contains

  subroutine config_LJdata(CF, Nat, Nso)
    use ParseText
    implicit none
    type(PTo), intent(in) :: CF
    integer, intent(in) :: Nat, Nso

    integer :: i
    character(len=2) LJindex
    double precision :: skin_factor

    skin_factor = PTread_d(CF, 'LJskinfactor')

    at_at%Na = Nat
    at_at%Nb = Nat

    allocate(at_at%eps(at_at%Na,at_at%Nb))
    allocate(at_at%sig(at_at%Na,at_at%Nb))
    allocate(at_at%cut(at_at%Na,at_at%Nb))
    allocate(at_at%neigh(at_at%Na,at_at%Nb))

    do i=1,at_at%Na
       write(LJindex,'(i02.2)') i
       at_at%eps(i,:) = PTread_dvec(CF,'at_ateps'//LJindex,at_at%Nb)
       at_at%sig(i,:) = PTread_dvec(CF,'at_atsig'//LJindex,at_at%Nb)
    end do

    at_at%cut = at_at%sig * 2.d0**(1.d0/6.d0)
    at_at%neigh = at_at%cut * skin_factor

    at_so%Na = Nat
    at_so%Nb = Nso

    allocate(at_so%eps(at_so%Na,at_so%Nb))
    allocate(at_so%sig(at_so%Na,at_so%Nb))
    allocate(at_so%cut(at_so%Na,at_so%Nb))
    allocate(at_so%neigh(at_so%Na,at_so%Nb))

    do i=1,at_so%Na
       write(LJindex,'(i02.2)') i
       at_so%eps(i,:) = PTread_dvec(CF,'at_soeps'//LJindex,at_so%Nb)
       at_so%sig(i,:) = PTread_dvec(CF,'at_sosig'//LJindex,at_so%Nb)
    end do
    
    at_so%cut = at_so%sig * 2.d0**(1.d0/6.d0)
    at_so%neigh = at_so%cut * skin_factor
    
  end subroutine config_LJdata


end module LJ
