program test
  use sys
  use MPCD
  use ParseText
  implicit none
  
  type(PTo) :: CF

  integer :: i_time, i_in
  double precision :: MPCD_kin, MPCD_kin_0, MPCD_mom(3), MPCD_mom_0(3)

  call PTparse(CF,'sample_MPCD',9)

  call config_sys(so_sys,'so',CF)

  call config_MPCD(CF)

  call homogeneous_solvent
  
  call PTkill(CF)

  tau=1.d0

  call compute_en_mom(MPCD_kin_0,MPCD_mom_0)
!  write(*,'(4e25.18)') MPCD_kin_0,MPCD_mom_0

  open(11,file='kin_mom')

  do i_time = 1,100
     
     do i_in = 1,10
        call place_in_cells
        call compute_v_com
        call generate_omega
        call simple_MPCD_step
        call MPCD_stream
        call correct_so

     end do

     call compute_en_mom(MPCD_kin,MPCD_mom)
     write(11,'(4e28.18)') MPCD_kin-MPCD_kin_0,MPCD_mom-MPCD_mom_0

  end do

  close(11)

contains
  
  subroutine compute_en_mom(kin, mom)
    double precision, intent(out) :: kin, mom(3)
    integer :: i

    kin = 0.d0 ; mom = 0.d0
    do i=1,so_sys%N(0)
       kin = kin + sum(so_v(:,i)**2)*0.5d0*so_sys%mass(so_species(i))
       mom = mom + so_v(:,i)*so_sys%mass(so_species(i))
    end do
    
  end subroutine compute_en_mom

  subroutine correct_so
    integer :: i, dim

    do i=1,so_sys%N(0)
       do dim=1,3
          if (so_r(dim,i) < 0.d0) so_r(dim,i) = so_r(dim,i) + L(dim)
          if (so_r(dim,i) >= L(dim)) so_r(dim,i) = so_r(dim,i) - L(dim)
       end do
    end do
  end subroutine correct_so

end program test

