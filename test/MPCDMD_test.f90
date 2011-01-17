program test
  use sys
  use group
  use LJ
  use MPCD
  use MD
  use ParseText
  implicit none
  
  type(PTo) :: CF

  integer :: i_time, i_in, i, reneigh
  integer :: N_MD_loop, N_loop, en_unit, at_unit
  double precision :: max_d

  call init_random_seed()

  call PTparse(CF,'sample_MPCDMD',9)

  call config_sys(so_sys,'so',CF)
  call config_sys(at_sys,'at',CF)

  call config_MPCD(CF)
  call config_MD

  call config_LJdata(CF, at_sys%N_species, so_sys%N_species)

  call init_atoms(CF)

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

  call fill_with_solvent
  call place_in_cells
  call make_neigh_list

  N_loop = PTread_i(CF, 'N_loop')
  N_MD_loop = PTread_i(CF, 'N_MD_loop')
  DT = PTread_d(CF, 'DT')

  call PTkill(CF)

  
  at_f => at_f1
  at_f_old => at_f2
  so_f => so_f1
  so_f_old => so_f2

  call compute_f

  en_unit = 11
  open(en_unit,file='energy')
  at_unit=12
  open(at_unit,file='at_x')

  call compute_tot_mom_energy(en_unit)

  reneigh = 0
  max_d = min( minval( at_at%neigh - at_at%cut ) , minval( at_so%neigh - at_so%cut ) ) * 0.5d0
  write(*,*) 'max_d = ', max_d

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
     end do
     
     write(at_unit,*) at_r(:,1)

     call correct_so
     call place_in_cells
     call compute_v_com
     call generate_omega
     call simple_MPCD_step

     call compute_tot_mom_energy(en_unit)
  end do

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

end program test

