program test
  use sys
  use group
  use LJ
  use MPCD
  use MD
  use ParseText
  use MPCDMD
  use h5md
  implicit none
  
  type(PTo) :: CF

  integer :: i_time, i_in, i, istart, reneigh
  integer :: N_MD_loop, N_loop, en_unit, at_x_unit, at_v_unit
  double precision :: max_d, realtime
  character(len=10) :: at_format
  integer :: collect_atom
  integer :: seed
  double precision :: at_sol_en, at_at_en, sol_kin, at_kin, energy
  double precision :: v_com(3), mass_temp

  integer(HID_T) :: file_ID
  type(h5md_t) :: posID
  type(h5md_t) :: enID, so_kinID, at_kinID, at_soID, at_atID, v_com_ID

  call MPCDMD_info
  call mtprng_info(short=.true.)
  call PTinfo(short=.true.)

  call PTparse(CF,'sample_MPCDMD',9)

  seed = PTread_i(CF,'seed')
  if (seed < 0) then
     seed = nint(100*secnds(0.))
  end if
  call mtprng_init(seed, ran_state)

  call config_sys(so_sys,'so',CF)
  call config_sys(at_sys,'at',CF)

  N_groups = PTread_i(CF, 'N_groups')

  if (N_groups <= 0) stop 'Ngroups is not a positive integer'
  allocate(group_list(N_groups))
  istart = 1
  do i=1,N_groups
     call config_group(group_list(i),i,istart,CF)
     istart = istart + group_list(i)%N
  end do


  if (at_sys%N_max<sum(group_list(:)%N)) stop 'at_sys%N_max < # atoms from group_list'

  call config_LJdata(CF, at_sys%N_species, so_sys%N_species)


  call config_MPCD(CF)

  call config_MD
  
  do i=1,N_groups
     if (group_list(i)%g_type == ATOM_G) then
        call config_atom_group(group_list(i))
     else if (group_list(i)%g_type == DIMER_G) then
        call config_dimer_group(group_list(i))
     else
        stop 'unknown group type'
     end if
  end do

  call init_atoms(CF)
  at_v = 0.d0

  write(*,*) so_sys%N_species
  write(*,*) so_sys%N_max
  write(*,*) so_sys%N

  write(*,*) at_sys%N_species
  write(*,*) at_sys%N_max
  write(*,*) at_sys%N

  write(*,*) at_at%eps
  write(*,*) at_at%sig

  write(*,*) at_so%eps
  write(*,*) at_so%sig

  write(*,*) so_species(1:10)
  write(*,*) at_species

  write(*,*) at_so%smooth
  write(*,*) at_at%smooth

  call fill_with_solvent(PTread_d(CF,'so_T'))
  call place_in_cells
  call make_neigh_list

  N_loop = PTread_i(CF, 'N_loop')
  N_MD_loop = PTread_i(CF, 'N_MD_loop')
  DT = PTread_d(CF, 'DT')
  h = PTread_d(CF, 'h')
  collect_atom = PTread_i(CF,'collect_atom')

  call PTkill(CF)

  
  at_f => at_f1
  at_f_old => at_f2
  so_f => so_f1
  so_f_old => so_f2

  call compute_f

  en_unit = 11
  open(en_unit,file='energy')
  at_x_unit=12
  open(at_x_unit,file='at_x')
  at_v_unit=13
  open(at_v_unit,file='at_v')
  write(at_format,'(a1,i02.2,a7)') '(', 3*at_sys%N(0),'e20.10)'
  write(*,*) at_format
  
  i_time = 0
  realtime = 0.d0
  call compute_tot_mom_energy(en_unit, at_sol_en, at_at_en, sol_kin, at_kin, energy)

  call begin_h5md

  call h5md_write_obs(at_soID, at_sol_en, i_time, realtime)
  call h5md_write_obs(at_atID, at_at_en, i_time, realtime)
  call h5md_write_obs(at_kinID, at_kin, i_time, realtime)
  call h5md_write_obs(so_kinID, sol_kin, i_time, realtime)
  call h5md_write_obs(enID, energy, i_time, realtime)
  call h5md_write_obs(v_com_ID, v_com, i_time, realtime)

  reneigh = 0
  max_d = min( minval( at_at%neigh - at_at%cut ) , minval( at_so%neigh - at_so%cut ) ) * 0.5d0
  write(*,*) 'max_d = ', max_d

  !at_v(:,1) = (/ 0.02d0, 0.01d0, 0.015d0 /)

  do i_time = 1,N_loop
     
     do i_in = 1,N_MD_loop
        call MD_step1

        if ( (maxval( sum( (so_r - so_r_neigh)**2 , dim=1 ) ) > max_d**2) .or. (maxval( sum( (at_r - at_r_neigh)**2 , dim=1 ) ) > max_d**2)) then
           reneigh = reneigh + 1
           call correct_so
           call place_in_cells
           call make_neigh_list
        end if

        call compute_f
        call MD_step2

        if (collect_atom > 0) then
           write(at_x_unit,at_format) ( at_r(:,i) , i=1,at_sys%N(0) )
           write(at_v_unit,at_format) ( at_v(:,i) , i=1,at_sys%N(0) )
        end if

        realtime=realtime+DT

     end do

     call correct_at
     
     if (collect_atom == 0) then
        write(at_x_unit,at_format) ( at_r(:,i) , i=1,at_sys%N(0) )
        write(at_v_unit,at_format) ( at_v(:,i) , i=1,at_sys%N(0) )
     end if

     call correct_so
     call place_in_cells
     call compute_v_com
     call generate_omega
     call simple_MPCD_step

     call compute_tot_mom_energy(en_unit, at_sol_en, at_at_en, sol_kin, at_kin, energy)

     v_com = 0.d0
     mass_temp = 0.d0
     do i=1,at_sys%N(0)
        v_com = v_com + at_sys%mass(at_species(i)) * at_v(:,i)
        mass_temp = mass_temp + at_sys%mass(at_species(i))
     end do
     v_com = v_com / mass_temp

     call h5md_write_obs(at_soID, at_sol_en, i_time, realtime)
     call h5md_write_obs(at_atID, at_at_en, i_time, realtime)
     call h5md_write_obs(at_kinID, at_kin, i_time, realtime)
     call h5md_write_obs(so_kinID, sol_kin, i_time, realtime)
     call h5md_write_obs(enID, energy, i_time, realtime)
     call h5md_write_obs(v_com_ID, v_com, i_time, realtime)
     call h5md_write_trajectory_data_d(posID, at_r, i_time, realtime)
  end do

  i_time = i_time-1

  call end_h5md

  write(*,*) reneigh, ' extra reneighbourings for ', N_loop*N_MD_loop, ' total steps'

contains
  
  subroutine correct_so
    integer :: i, dim

    do i=1,so_sys%N(0)
       do dim=1,3
          if (so_r(dim,i) < 0.d0) so_r(dim,i) = so_r(dim,i) + L(dim)
          if (so_r(dim,i) >= L(dim)) so_r(dim,i) = so_r(dim,i) - L(dim)
       end do
    end do
  end subroutine correct_so

  subroutine correct_at
    integer :: i, dim

    do i=1,at_sys%N(0)
       do dim=1,3
          if (at_r(dim,i) < 0.d0) at_r(dim,i) = at_r(dim,i) + L(dim)
          if (at_r(dim,i) >= L(dim)) at_r(dim,i) = at_r(dim,i) - L(dim)
       end do
    end do
  end subroutine correct_at

  subroutine begin_h5md
    call h5open_f(h5_error)
    call h5md_create_file(file_ID, 'data.h5', 'MPCDMD')

    call h5md_add_trajectory_data(file_ID, 'position', at_sys% N_max, 3, posID)
    call h5md_create_obs(file_ID, 'energy', enID, energy)
    call h5md_create_obs(file_ID, 'at_at_int', at_atID, at_at_en, link_from='energy')
    call h5md_create_obs(file_ID, 'at_so_int', at_soID, at_sol_en, link_from='energy')
    call h5md_create_obs(file_ID, 'so_kin', so_kinID, sol_kin, link_from='energy')
    call h5md_create_obs(file_ID, 'at_kin', at_kinID, at_kin, link_from='energy')
    call h5md_create_obs(file_ID, 'v_com', v_com_ID, v_com, link_from='energy')

  end subroutine begin_h5md

  subroutine end_h5md
    call h5md_close_ID(posID)
    call h5md_close_ID(enID)
    call h5md_close_ID(at_atID)
    call h5md_close_ID(at_soID)
    call h5md_close_ID(so_kinID)
    call h5md_close_ID(at_kinID)

    call h5fclose_f(file_ID, h5_error)
    
    call h5close_f(h5_error)

  end subroutine end_h5md

end program test

