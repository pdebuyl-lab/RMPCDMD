program test
  use sys
  use MPCD
  use ParseText
  use MPCDMD
  use volume_reaction
  implicit none
  
  type(PTo) :: CF

  integer :: i_time, i_in, N_outer, N_loop
  double precision :: MPCD_kin, MPCD_kin_0, MPCD_mom(3), MPCD_mom_0(3)
  integer :: seed
  integer :: i

  call MPCDMD_info
  call mtprng_info(short=.true.)
  call PTinfo(short=.true.)

  call PTparse(CF,'sample_MPCD',9)

  seed = PTread_i(CF,'seed')
  if (seed < 0) then
     seed = nint(100*secnds(0.))
  end if
  call mtprng_init(seed, ran_state)

  call config_sys(so_sys,'so',CF)

  call config_MPCD(CF)

  N_reactions = PTread_i(CF,'N_reactions')
  if (N_reactions > 0) allocate(reaction_list(N_reactions))

  if (N_reactions > 0) then
     do i=1,N_reactions
        call add_vol_reac(CF, reaction_list(i), i, so_sys % N_species)
     end do
  end if

  call homogeneous_solvent(PTread_d(CF,'so_T'))

  tau=PTread_d(CF,'tau')
  MPCD_tau = tau
  do_shifting = PTread_l(CF, 'shifting')
  N_outer = PTread_i(CF,'N_outer')
  N_loop = PTread_i(CF,'N_loop')
  
  call PTkill(CF)


  call compute_en_mom(MPCD_kin_0,MPCD_mom_0)
  write(*,'(4e28.18)') MPCD_kin_0, MPCD_mom_0

  open(11,file='kin_mom')
  open(12,file='N_of_t')
  write(12,*) so_sys % N

  shift = 0.d0
  
  do i_time = 1,N_outer
     
     do i_in = 1,N_loop
        call place_in_cells
        call compute_v_com
        call generate_omega
        call chem_MPCD_step
        call MPCD_stream

        if (do_shifting) then
           shift(1) = (mtprng_rand_real1(ran_state)-0.5d0)*a
           shift(2) = (mtprng_rand_real1(ran_state)-0.5d0)*a
           shift(3) = (mtprng_rand_real1(ran_state)-0.5d0)*a
        end if

        call correct_so

     end do

     call compute_en_mom(MPCD_kin,MPCD_mom)
     write(11,'(4e28.18)') MPCD_kin-MPCD_kin_0,MPCD_mom-MPCD_mom_0
     write(12,*) so_sys % N

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
          if (so_r(dim,i) < shift(dim)) so_r(dim,i) = so_r(dim,i) + L(dim)
          if (so_r(dim,i) >= L(dim)+shift(dim)) so_r(dim,i) = so_r(dim,i) - L(dim)
       end do
    end do
  end subroutine correct_so


end program test

