
!> This module holds the variable to embed molecular dynamics particles (called "atoms" in this program) in a MPCD solvent.
!!
!! A number of global variables are defined in this module.
module MD
  use mtprng
  use group
  use sys
  use LJ
  use MPCD
  use reaction
  use HDF5
  implicit none

  !> Maximum number of solvent neighbours to each atom.
  integer, parameter :: max_neigh=32768 !8192
  !> Maximum number of atom neighbours to each solvent particle.
  integer, parameter :: max_so_neigh=16

  !> Information about the "atoms" system (N, ...).
  type(sys_t) :: at_sys

  !> Number of atomic groups in the system.
  integer :: N_groups
  !> List of atomic groups in the system.
  type(group_t), allocatable :: group_list(:)
  
  double precision, allocatable :: at_r(:,:), at_r_old(:,:), at_v(:,:)
  double precision, allocatable, target :: at_f1(:,:), at_f2(:,:)
  double precision, pointer :: at_f(:,:), at_f_old(:,:), at_f_temp(:,:)
  double precision, allocatable :: at_r_neigh(:,:)
  integer, allocatable :: at_jumps(:,:)
  integer, allocatable :: at_neigh_list(:,:)
  integer, allocatable :: at_species(:)

  !> A reverse lookup list: for each solvent particle, a list of neighbouring atoms.
  integer, allocatable :: so_neigh_list(:,:)

  !> The timestep for the molecular dynamics integrator.
  double precision :: DT

  integer, allocatable :: reac_table(:,:)
  integer, allocatable :: reac_product(:,:)
  double precision, allocatable :: reac_rates(:,:)
  double precision :: excess, max_d
  type(reac_t), allocatable :: at_so_reac(:,:)
  logical(kind=1), allocatable :: so_do_reac(:)
  
  double precision :: h


contains


  !> Allocates arrays for the atom system.
  !! Should be called after configuration of at_sys and so_sys.
  subroutine config_MD
    implicit none

    integer :: i,j
    
    allocate(at_r(3,at_sys%N_max))
    allocate(at_r_old(3,at_sys%N_max))
    allocate(at_v(3,at_sys%N_max))
    allocate(at_f1(3,at_sys%N_max))
    allocate(at_f2(3,at_sys%N_max))
    allocate(at_r_neigh(3,at_sys%N_max))
    allocate(at_jumps(3,at_sys%N_max))
    allocate(at_species(at_sys%N_max))

    allocate(at_neigh_list(0:max_neigh, at_sys%N_max))

    allocate(so_neigh_list(0:max_so_neigh, so_sys%N_max) )
    allocate(so_do_reac( so_sys % N_max ) )
    allocate(at_so_reac( at_sys % N_species , so_sys % N_species ) )

    j=1
    do i=1,at_sys%N_species
       at_species(j:j-1+at_sys%N(i)) = i
       j = j+at_sys%N(i)
    end do

  end subroutine config_MD

  !> Sets the species for an atom group.
  !! @param g_var A group_t variable.
  subroutine config_atom_group(g_var)
    implicit none
    type(group_t), intent(in) :: g_var
    
    integer :: i

    do i=g_var%istart, g_var%istart + g_var%N - 1
       at_species(i) = g_var % species1
    end do
    
  end subroutine config_atom_group

  !> Sets the species for a dimer group.
  !! @param g_var A group_t variable.
  subroutine config_dimer_group(g_var)
    implicit none
    type(group_t), intent(in) :: g_var

    at_species(g_var%istart) = g_var % species1
    at_species(g_var%istart+1) = g_var % species2
    
  end subroutine config_dimer_group

  !> Places all atoms at random in the simulation box while avoiding overlapping.
  subroutine init_atoms_random
    
    double precision :: x(3), dsqr
    integer :: i,j, iter
    logical :: too_close
    
    i=1
    iter=1
    do while (.true.)
       iter=iter+1
       x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ; x(3) = mtprng_rand_real1(ran_state) ;
       at_r(:,i) = x*L
       too_close = .false.
       do j=1,i-1
          call rel_pos(at_r(:,i),at_r(:,j), L, x)
          dsqr = sum( x**2 )
          if (dsqr .lt. at_at%cut(at_species(j),at_species(i))**2) then
             too_close = .true.
             exit
          end if
       end do
       if (.not.too_close) then
          i=i+1
       end if
       if (i>at_sys%N(0)) exit
       if (iter>100*at_sys%N(0)) then
          write(*,*) 'tried to place atoms more than', 100*at_sys%N(0), 'times'
          write(*,*) 'failed to place all atoms'
          stop
       end if
    end do
    at_v = 0.d0

  end subroutine init_atoms_random

  !> Fills a simulation box with the solvent at the given temperature, but with a flat velocity profile.
  !! @param temperature The temperature for the solvent.
  subroutine fill_with_solvent(temperature)
    double precision, intent(in) :: temperature
    integer :: i, j, iter
    double precision :: x(3), dsqr, t_factor
    logical :: too_close
    double precision :: tot_m, tot_v(3)

    t_factor = sqrt(3.d0*temperature)
    tot_m = 0.d0
    tot_v = 0.d0

    i=1
    iter=1
    do while (.true.)
       iter = iter+1
       x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ;  x(3) = mtprng_rand_real1(ran_state) ;
       so_r(:,i) = x*L
       too_close = .false.
       do j=1,at_sys%N(0)
          call rel_pos(so_r(:,i), at_r(:,j), L, x)
          dsqr = sum( x**2 )
          if (dsqr .lt. 2.d0**(1.d0/3.d0)*at_so%sig( at_species(j) , so_species(i) )**2 ) then
             too_close = .true.
             exit
          end if
       end do
       if (.not. too_close) then
          x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ;  x(3) = mtprng_rand_real1(ran_state) ;
          x = x-0.5d0
          so_v(:,i) = x*2.d0 * t_factor/sqrt(so_sys%mass(so_species(i)))
          tot_m = tot_m + so_sys%mass(so_species(i))
          tot_v = tot_v + so_sys%mass(so_species(i)) * so_v(:,i)
          i=i+1
       end if
       if (i>so_sys%N(0)) exit
       if (iter>100*so_sys%N(0)) then
          write(*,*) 'tried to place solvent particles more than', 100*so_sys%N(0), 'times'
          write(*,*) 'failed to place all solvent particle'
          stop
       end if
    end do
    tot_v = tot_v / tot_m
    do i=1,so_sys%N(0)
       so_v(:,i) = so_v(:,i) - tot_v
    end do
    
  end subroutine fill_with_solvent

  subroutine init_atoms(CF)
    use ParseText
    implicit none
    type(PTo), intent(in) :: CF

    character(len=24) :: at_init_mode

    at_init_mode = PTread_s(CF,'atinit')

    if (at_init_mode.eq.'center') then
       if (at_sys%N(0).eq.1) then
          at_r(:,1) = L/2.d0
          at_v(:,1) = 0.d0
       else
          write(*,*) 'too many atoms for center init mode, stopping'
          stop
       end if
    else if (at_init_mode.eq.'random') then
       call init_atoms_random
    else
       write(*,*) 'unknown atinit mode ', at_init_mode, ' stopping'
       stop
    end if

  end subroutine init_atoms

  subroutine make_neigh_list

    integer :: at_i, i, part
    integer :: ci, cj, ck, mi, mj, mk
    integer :: Si, Sj, Sk
    integer :: extent
    double precision :: dist_sqr, neigh_sqr, x(3)

    is_MD = .false.

    do i=1,so_sys%N(0)
       N_MD(i) = max_d/(sqrt(sum(so_v(:,i)**2))*DT)
    end do

    N_MD_max = minval(N_MD(1:so_sys%N(0)))
    
    at_neigh_list = 0
    so_neigh_list = 0

    do at_i=1,at_sys%N(0)
       
       extent = ceiling( maxval(at_so%neigh(at_species(at_i),:) ) * oo_a )

       Si = floor( (at_r(1,at_i)-shift(1)) * oo_a ) + 1
       Sj = floor( (at_r(2,at_i)-shift(2)) * oo_a ) + 1
       Sk = floor( (at_r(3,at_i)-shift(3)) * oo_a ) + 1
       do ck= Sk - extent, Sk + extent
          do cj = Sj - extent, Sj + extent
             do ci = Si - extent, Si + extent
                mi = modulo(ci-1,N_cells(1)) + 1 ; mj = modulo(cj-1,N_cells(2)) + 1 ; mk = modulo(ck-1,N_cells(3)) + 1 ; 
                do i=1,par_list(0,mi,mj,mk)
                   call rel_pos(so_r(:,par_list(i,mi,mj,mk)) , at_r(:,at_i) , L, x)
                   dist_sqr = sum( x**2 )
                   neigh_sqr = at_so%neigh( at_species(at_i), so_species(par_list(i,mi,mj,mk)))**2
                   if ( dist_sqr .le. neigh_sqr ) then
                      at_neigh_list(0,at_i) = at_neigh_list(0,at_i) + 1
                      if (at_neigh_list(0,at_i) > max_neigh) then
                         write(*,*) 'too many neighbours for atom',at_i
                         stop
                      end if
                      at_neigh_list(at_neigh_list(0,at_i),at_i) = par_list(i,mi,mj,mk)
                      part = par_list(i,mi,mj,mk)
                      so_neigh_list(0,part) = so_neigh_list(0,part) + 1
                      if ( so_neigh_list(0,part) > max_so_neigh ) stop 'too many neighbours for solvent'
                      so_neigh_list(so_neigh_list(0,part),part) = at_i
                      is_MD(par_list(i,mi,mj,mk)) = .true.
                   end if
                end do
             end do
          end do
       end do
    end do

    so_r_neigh = so_r
    at_r_neigh = at_r

  end subroutine make_neigh_list

  subroutine compute_f
    implicit none

    integer :: at_i, at_j, j, part, at_si, at_g, at_h, at_j_1
    double precision :: x(3), dist_sqr, LJcut_sqr, LJsig, f_var(3)
    double precision :: dist_min, at_dist_min
    
    so_f_temp => so_f
    so_f => so_f_old
    so_f_old => so_f_temp
    
    at_f_temp => at_f
    at_f => at_f_old
    at_f_old => at_f_temp

    so_f = 0.d0
    at_f = 0.d0

    dist_min = ( L(1) + L(2) + L(3) ) **2

    do at_i=1,at_sys%N(0)
       at_si = at_species(at_i)

       do j=1, at_neigh_list(0,at_i)
          part = at_neigh_list(j, at_i)
          call rel_pos(so_r(:,part), at_r(:,at_i), L, x)
          dist_sqr = sum( x**2 )
          if (dist_sqr < dist_min) dist_min=dist_sqr
          if ( dist_sqr .le. at_so%cut(at_si, so_species(part))**2 ) then
             if (at_so%smooth(at_si,so_species(part))) then
                f_var = LJ_force_smooth_or( &
                at_so%eps( at_si,so_species(part) ), at_so%sig(at_si,so_species(part)), &
                at_so%V_c( at_si,so_species(part) ), &
                dist_sqr, at_so%cut(at_si,so_species(part)), h ) * x
             else
                f_var = LJ_force_or(at_so%eps( at_si,so_species(part) ), at_so%sig(at_si,so_species(part)), dist_sqr) * x
             end if
             so_f(:,part) = so_f(:,part) + f_var
             at_f(:,at_i) = at_f(:,at_i) - f_var
          end if
       end do
    end do

    at_dist_min = ( L(1) + L(2) + L(3) ) **2

    do at_g = 1, N_groups
       
       if (group_list(at_g)%g_type .eq. ELAST_G) then
          call compute_f_elast(at_g)
       end if

       do at_h = at_g, N_groups
          
          if ( (at_g.eq.at_h) .and. (group_list(at_g) % g_type .ne. ATOM_G) ) cycle

          do at_i = group_list(at_g) % istart, group_list(at_g) % istart + group_list(at_g) % N - 1
             if (at_g .eq. at_h) then
                at_j_1 = at_i+1
             else
                at_j_1 = group_list(at_h) % istart
             end if
             do at_j = at_j_1, group_list(at_h) % istart + group_list(at_h) % N - 1

                call rel_pos(at_r(:,at_i), at_r(:,at_j), L, x)
                LJsig = at_at%sig( at_species(at_i), at_species(at_j) )
                LJcut_sqr = at_at%cut( at_species(at_i), at_species(at_j) )**2
                dist_sqr = sum( x**2 )
                if (dist_sqr .lt. at_dist_min) at_dist_min = dist_sqr
                if ( dist_sqr .le. LJcut_sqr ) then
                   if (at_at%smooth(at_species(at_i), at_species(at_j))) then
                      f_var = LJ_force_smooth_or( &
                           at_at%eps( at_species(at_i),at_species(at_j) ) , LJsig, &
                           at_at%V_c( at_species(at_i),at_species(at_j) ), dist_sqr, &
                           at_at%cut(at_species(at_i),at_species(at_j)), h) * x
                   else
                      f_var = LJ_force_or(at_at%eps( at_species(at_i),at_species(at_j) ) , LJsig, dist_sqr) * x
                   end if
                   at_f(:, at_i) = at_f(:,at_i) + f_var
                   at_f(:, at_j) = at_f(:,at_j) - f_var
                end if
             end do
          end do
       end do
    end do

  end subroutine compute_f

  subroutine MD_step1

    integer :: at_i, i

    do i=1,so_sys%N(0)
       if (is_MD(i)) so_r(:,i) = so_r(:,i) + so_v(:,i) * DT + so_f(:,i) * DT**2 * 0.5d0 * so_sys % oo_mass(so_species(i))
    end do
    do at_i=1,at_sys%N(0)
       at_r(:,at_i) = at_r(:,at_i) + at_v(:,at_i) * DT + at_f(:,at_i) * DT**2 * 0.5d0 * at_sys % oo_mass( at_species(at_i) )
    end do
    
  end subroutine MD_step1

  subroutine MD_step2

    integer :: at_i, i

    do i=1,so_sys%N(0)
       if (is_MD(i)) so_v(:,i) = so_v(:,i) + 0.5d0 * DT * (so_f(:,i) + so_f_old(:,i)) * so_sys % oo_mass( so_species(i) )
    end do

    do at_i=1, at_sys%N(0)
       at_v(:,at_i) = at_v(:,at_i) + 0.5d0 * DT * (at_f(:,at_i) + at_f_old(:,at_i) ) * at_sys % oo_mass( at_species(at_i) )
    end do

  end subroutine MD_step2

  subroutine compute_tot_mom_energy(file_unit, at_sol_en, at_at_en, sol_kin, at_kin, energy, total_v)
    integer, intent(in) :: file_unit
    double precision, intent(out) :: at_sol_en, at_at_en, sol_kin, at_kin, energy, total_v(3)
    double precision :: mom(3), at_mom(3), mass, at_mass

    integer :: at_i, at_j, j, part, at_si
    integer :: at_g, at_h, at_j_1
    double precision :: LJcut_sqr, LJsig, x(3), dist_sqr

    at_sol_en = 0.d0 ; at_at_en = 0.d0 ; sol_kin = 0.d0 ; at_kin = 0.d0
    mom = 0.d0 ; at_mom = 0.d0 ; mass = 0.d0 ; at_mass = 0.d0

    do at_i=1,at_sys%N(0)
       do j=1, at_neigh_list(0,at_i)
          part = at_neigh_list(j, at_i)
          LJsig = at_so%sig( at_species( at_i ) , so_species(part) )
          LJcut_sqr = at_so%cut( at_species(at_i), so_species(part) )**2
          call rel_pos(so_r(:,part), at_r(:,at_i), L, x)
          dist_sqr = sum( x**2 )
          if ( dist_sqr .le. LJcut_sqr ) then
             if (at_so%smooth(at_species(at_i),so_species(part))) then
                at_sol_en = at_sol_en + LJ_V_smooth( &
                     at_so%eps( at_species(at_i), so_species(part)) , &
                     LJsig, &
                     at_so%V_c( at_species(at_i), so_species(part)) , &
                     dist_sqr , at_so%cut( at_species(at_i), so_species(part) ), h)
             else
                at_sol_en = at_sol_en + LJ_V( at_so%eps( at_species(at_i), so_species(part)) , LJsig, at_so%V_c( at_species(at_i), so_species(part)), dist_sqr )
             end if
          end if
       end do
    end do

    do part=1,so_sys%N(0)
       sol_kin = sol_kin + 0.5d0 * so_sys%mass( so_species(part) ) * sum( so_v(:,part)**2 )
       mass = mass + so_sys%mass(so_species(part))
       mom = mom + so_sys%mass(so_species(part)) * so_v(:,part)
    end do

    do at_g = 1, N_groups

       if (group_list(at_g)%g_type .eq. ELAST_G) then
          at_at_en = at_at_en + compute_pot_elast(at_g)
       end if

       do at_h = at_g, N_groups

          if ( (at_g.eq.at_h) .and. (group_list(at_g) % g_type .ne. ATOM_G) ) cycle
    
          do at_i = group_list(at_g) % istart, group_list(at_g) % istart + group_list(at_g) % N - 1
             if (at_g .eq. at_h) then
                at_j_1 = at_i+1
             else
                at_j_1 = group_list(at_h) % istart
             end if
             do at_j = at_j_1, group_list(at_h) % istart + group_list(at_h) % N - 1
                call rel_pos( at_r(:,at_i), at_r(:,at_j), L, x)
                LJsig = at_at%sig( at_species(at_i),  at_species(at_j) )
                LJcut_sqr = at_at%cut( at_species(at_i),  at_species(at_j) )**2
                dist_sqr = sum( x**2 )
                if ( dist_sqr .le. LJcut_sqr ) then
                   if (at_at%smooth(at_species(at_i),at_species(at_j))) then
                      at_at_en = at_at_en + LJ_V_smooth( &
                           at_at%eps( at_species(at_i), at_species(at_j) ) , &
                           LJsig, at_at%V_c( at_species(at_i), at_species(at_j) ), &
                           dist_sqr , at_at%cut( at_species(at_i),  at_species(at_j) ) , h)
                   else
                      at_at_en = at_at_en + LJ_V(at_at%eps( at_species(at_i), at_species(at_j) ) , LJsig, at_at%V_c( at_species(at_i), at_species(at_j) ), dist_sqr )
                   end if
                end if
             end do
          end do
    
       end do
       
    end do

    do at_i = 1,at_sys%N(0)
       at_si = at_species(at_i)
       at_kin = at_kin + 0.5d0 * at_sys%mass( at_si ) * sum( at_v(:,at_i)**2 )
       at_mass = at_mass + at_sys%mass( at_si )
       at_mom = at_mom + at_sys%mass( at_si ) * at_v(:,at_i)
    end do

    excess = 0.d0
    do j=1, so_sys % N_species
       excess = excess + u_int(j)*so_sys % N(j)
    end do
    
    if (file_unit > 0) write(file_unit,'(7e30.20)') at_sol_en, at_at_en, sol_kin, at_kin, &
         excess, at_sol_en+at_at_en+sol_kin+at_kin, at_sol_en+at_at_en+sol_kin+at_kin+excess
    energy = at_sol_en+at_at_en+sol_kin+at_kin+excess
    total_v = (mom + at_mom) / (mass + at_mass)

  end subroutine compute_tot_mom_energy

  subroutine rel_pos(r1, r2, Lvar, rvar)
    double precision, intent(in) :: r1(3), r2(3), Lvar(3)
    double precision, intent(out) :: rvar(3)

    integer :: dim

    rvar = r1-r2
    do dim=1,3
       if ( rvar(dim) < -0.5d0*Lvar(dim) ) then
          rvar(dim) = rvar(dim) + Lvar(dim)
       else if ( rvar(dim) > 0.5d0*Lvar(dim) ) then
          rvar(dim) = rvar(dim) - Lvar(dim)
       end if
    end do

  end subroutine rel_pos

  subroutine config_reac_MD(CF)
    use sys
    use ParseText
    implicit none
    type(PTo), intent(in) :: CF
    character(len=14) :: temp_name
    integer :: i

    allocate(reac_table(at_sys%N_species,so_sys%N_species))
    allocate(reac_product(at_sys%N_species,so_sys%N_species))
    allocate(reac_rates(at_sys%N_species,so_sys%N_species))

    reac_table = 0
    reac_product = 0
    reac_rates = 0.d0

    do i=1,at_sys%N_species
       write(temp_name,'(a10,i02.2)') 'reac_table', i
       reac_table(i,:) = PTread_ivec(CF,trim(temp_name),so_sys%N_species)
       write(temp_name,'(a12,i02.2)') 'reac_product', i
       reac_product(i,:) = PTread_ivec(CF,trim(temp_name),so_sys%N_species)
       write(temp_name,'(a10,i02.2)') 'reac_rates', i
       reac_rates(i,:) = PTread_dvec(CF,trim(temp_name),so_sys%N_species)
    end do
    
  end subroutine config_reac_MD

  subroutine reac_MD_do(at_i, part, delta_U)
    implicit none
    integer, intent(in) :: at_i, part
    double precision, intent(out) :: delta_U

    double precision :: LJsig, dist_sqr, x(3)

    call rel_pos(at_r(:,at_i),so_r(:,part),L,x)
    dist_sqr = sum(x**2)

    if (dist_sqr .le. at_so%cut(at_species(at_i),so_species(part))**2) then

       delta_U = 0.d0
       
       LJsig = at_so%sig(at_species(at_i),so_species(part))
       if (dist_sqr .le. at_so%cut(at_species(at_i),so_species(part))**2) then
          delta_U = delta_U + 4.d0 * at_so%eps( at_species(at_i), so_species(part)) * &
               ( (LJsig**2/dist_sqr)**6 - (LJsig**2/dist_sqr)**3 + 0.25d0 )
       end if
       delta_U = delta_U + u_int(so_species(part))
       so_sys%N(so_species(part)) = so_sys%N(so_species(part)) - 1
       
       so_species(part) = reac_product(at_species(at_i),so_species(part))
       
       LJsig = at_so%sig(at_species(at_i),so_species(part))
       if (dist_sqr .le. at_so%cut(at_species(at_i),so_species(part))**2) then
          delta_U = delta_U - 4.d0 * at_so%eps( at_species(at_i), so_species(part)) * &
               ( (LJsig**2/dist_sqr)**6 - (LJsig**2/dist_sqr)**3 + 0.25d0 )
       end if
       delta_U = delta_U + u_int(so_species(part))
       so_sys%N(so_species(part)) = so_sys%N(so_species(part)) + 1
       
       excess = excess + delta_U

    end if

  end subroutine reac_MD_do

  subroutine config_elast_group(g_var)
    use ParseText
    type(group_t), intent(inout) :: g_var

    integer :: i, g_n
    
    i = g_var % istart
    g_n = g_var % N
    at_species(i:i+g_n/2-1) = g_var % species1
    at_species(i+g_n/2:i+g_n-1) = g_var % species2

  end subroutine config_elast_group

  subroutine config_elast_group2(g_var)
    use ParseText
    type(group_t), intent(inout) :: g_var

    integer :: i,j, N_link, i_link
    double precision, allocatable :: dist_table(:,:) 
    double precision :: x(3)
    

    allocate(dist_table(g_var%N,g_var%N))

    do i=1,g_var%N
       do j=i+1,g_var%N
          call rel_pos(at_r(:,g_var%istart+i-1),at_r(:,g_var%istart+j-1),L,x)
          dist_table(i,j) = sqrt(sum( x**2 ))
          dist_table(j,i) = dist_table(i,j)
       end do
    end do

    N_link = 0
    do i=1,g_var%N
       do j=i+1,g_var%N
          if (dist_table(j,i) .le. g_var%elast_rmax) then
             N_link = N_link + 1
          end if
       end do
    end do

    if (N_link .eq. 0) then
       stop 'elast group with 0 links, stopping.'
    end if

    g_var%nlink = N_link

    allocate(g_var%index(2,N_link))
    allocate(g_var%r0(N_link))
    i_link = 1
    do i=1,g_var%N
       do j=i+1,g_var%N
          if (dist_table(j,i) .le. g_var%elast_rmax) then
             g_var%index(1,i_link) = g_var%istart+i-1
             g_var%index(2,i_link) = g_var%istart+j-1
             g_var%r0(i_link) = dist_table(j,i)
             i_link = i_link + 1
          end if
       end do
    end do

    deallocate(dist_table)

  end subroutine config_elast_group2

  subroutine compute_f_elast(g_i)
    integer, intent(in) :: g_i

    integer :: i, part1, part2
    double precision :: r, f(3), x(3)

    do i=1, group_list(g_i)%nlink

       part1 = group_list(g_i) % index(1,i)
       part2 = group_list(g_i) % index(2,i)

       call rel_pos( at_r(:,part1) , at_r(:,part2) , L, x)
       
       r = sqrt( sum( x**2 ) )

       f = - group_list(g_i) % elast_k * ( r - group_list(g_i) % r0(i) ) * x / r

       at_f(:,part1) = at_f(:,part1) + f
       at_f(:,part2) = at_f(:,part2) - f

    end do

  end subroutine compute_f_elast

  function compute_pot_elast(g_i)
    double precision :: compute_pot_elast
    integer, intent(in) :: g_i

    integer :: i, part1, part2
    double precision :: r, x(3), u

    u = 0.d0

    do i=1, group_list(g_i)%nlink

       part1 = group_list(g_i) % index(1,i)
       part2 = group_list(g_i) % index(2,i)

       call rel_pos( at_r(:,part1) , at_r(:,part2) , L, x)
       
       r = sqrt( sum( x**2 ) )

       u = u + 0.5d0 * group_list(g_i) % elast_k * ( r - group_list(g_i) % r0(i) )**2

    end do

    compute_pot_elast = u

  end function compute_pot_elast

  !> Computes the center of mass position of a group or subgroup.
  !! @param g_var The group to consider
  !! @param sub_g Optionally, the index of a subgroup.
  !! @return com_r The coordinates.
  function com_r(g_var, sub_g)
    double precision :: com_r(3)
    type(group_t), intent(inout) :: g_var
    integer, intent(in), optional :: sub_g
    
    integer :: i, i1,iN
    double precision :: mass
    
    if (present(sub_g)) then
       i1 = g_var % subgroup(1, sub_g)
       iN = g_var % subgroup(2, sub_g)
    else
       i1 = g_var % istart
       iN = g_var % istart + g_var % N - 1
    end if
    
    com_r = 0.d0 ; mass = 0.d0
    do i = i1, iN
       com_r = com_r + (at_r(:,i) + at_jumps(:,i) * L ) * at_sys % mass(at_species(i))
       mass = mass + at_sys % mass(at_species(i))
    end do
    
    com_r = com_r / mass

  end function com_r

  !> Computes the center of mass velocity of a group or subgroup.
  !! @param g_var The group to consider
  !! @param sub_g Optionally, the index of a subgroup.
  !! @return com_r The velocity.
  function com_v(g_var, sub_g)
    double precision :: com_v(3)
    type(group_t), intent(inout) :: g_var
    integer, intent(in), optional :: sub_g
    
    integer :: i, i1,iN
    double precision :: mass
    
    if (present(sub_g)) then
       i1 = g_var % subgroup(1, sub_g)
       iN = g_var % subgroup(2, sub_g)
    else
       i1 = g_var % istart
       iN = g_var % istart + g_var % N - 1
    end if
    
    com_v = 0.d0 ; mass = 0.d0
    do i = i1, iN
       com_v = com_v + at_v(:,i) * at_sys % mass(at_species(i))
       mass = mass + at_sys % mass(at_species(i))
    end do
    
    com_v = com_v / mass

  end function com_v

  !> Computes the force applied on the center of mass of a group or subgroup.
  !! @param g_var The group to consider
  !! @param sub_g Optionally, the index of a subgroup.
  !! @return com_f The velocity.
  function com_f(g_var, sub_g)
    double precision :: com_f(3)
    type(group_t), intent(inout) :: g_var
    integer, intent(in), optional :: sub_g
    
    integer :: i, i1,iN
    
    if (present(sub_g)) then
       i1 = g_var % subgroup(1, sub_g)
       iN = g_var % subgroup(2, sub_g)
    else
       i1 = g_var % istart
       iN = g_var % istart + g_var % N - 1
    end if
    
    com_f = 0.d0
    do i = i1, iN
       com_f = com_f + at_f(:,i)
    end do
    
  end function com_f


  !> Applies the shake algorithm to the particles in the group g_var.
  !! The list of constraints should have been initialized.
  !! @param g_var A group_t variable.
  subroutine shake(g_var)
    implicit none
    type(group_t), intent(in) :: g_var

    integer :: i, j, k, iter, at_si, at_sj
    double precision :: eps=1d-12, max_err, lagrange, x_old(3), x(3), dist

    do iter=1,50000
       max_err = 0.d0
       do k=1,g_var % nlink
          i = g_var % index(1, k)
          at_si = at_species(i)
          j = g_var % index(2, k)
          at_sj = at_species(j)
          call rel_pos(at_r_old(:,i),at_r_old(:,j), L, x_old)
          call rel_pos(at_r(:,i),at_r(:,j), L, x)
          lagrange = &
               ( g_var % r0(k)**2 - sum( x**2 ) ) &
               / (at_sys % oo_mass(at_si) + at_sys % oo_mass(at_sj) ) &
               / 4.d0 &
               / sum( x*x_old )
          at_r(:,i) = at_r(:,i) + x_old * lagrange * (at_sys % oo_mass(at_si) + at_sys % oo_mass(at_sj))
          at_r(:,j) = at_r(:,j) - x_old * lagrange * (at_sys % oo_mass(at_si) + at_sys % oo_mass(at_sj))
          call rel_pos(at_r(:,i),at_r(:,j), L, x)
          dist = abs( sqrt(sum( x**2 )) - g_var % r0(k) )
          if ( dist > max_err ) max_err = dist
       end do
       if ( max_err < eps ) exit
    end do
    if ( max_err > eps ) then
       write(*,*) max_err, dist, iter
       stop 'shake fails'
    end if

  end subroutine shake

  !> Applies the rattle algorithm to the particles in the group g_var.
  !! The list of constraints should have been initialized.
  !! @param g_var A group_t variable.
  subroutine rattle(g_var)
    implicit none
    type(group_t), intent(in) :: g_var

    integer :: i, j, k, iter, at_si, at_sj
    double precision :: eps=1d-8, max_err, lagrange, x(3), dist, x_old(3)

    do iter=1,50000
       max_err = 0.d0
       do k=1,g_var % nlink
          i = g_var % index(1, k)
          at_si = at_species(i)
          j = g_var % index(2, k)
          at_sj = at_species(j)
          call rel_pos(at_r(:,i),at_r(:,j),L,x)
          call rel_pos(at_r_old(:,i), at_r_old(:,j), L, x_old)
          x = (x+x_old)*0.5d0
          lagrange = &
               sum( x * ( at_v(:,i) - at_v(:,j) )) / &
               ( g_var % r0(k)**2 * &
               (at_sys%oo_mass(at_si) + at_sys%oo_mass(at_sj)) )
          at_v(:,i) = at_v(:,i) - x * lagrange * at_sys % oo_mass(at_si)
          at_v(:,j) = at_v(:,j) + x * lagrange * at_sys % oo_mass(at_sj)
          dist = sqrt( sum( ( at_v(:,i)-at_v(:,j) ) * x ) )
          if ( dist > max_err ) max_err = dist
       end do
       if ( max_err < eps ) exit
    end do
    if ( max_err > eps ) then
       write(*,*) max_err, dist, iter
       stop 'rattle fails'
    end if

  end subroutine rattle

  subroutine reac_loop
    implicit none

    integer :: at_i, j, part, at_si, i_neigh, neigh_idx
    integer :: g_i
    integer :: p1si, p2si
    double precision :: x(3)
    double precision :: dist_min
    double precision :: delta_u
    logical :: too_many_atoms, thermal
    double precision :: dist_sqr, neigh_sqr
    integer :: total_N, ci, cj, ck, cc(3)

    total_N = so_sys % N(0)

    dist_min = sum(L)
    do g_i=1,N_groups
       do at_i=group_list(g_i)%istart, group_list(g_i)%istart + group_list(g_i)%N - 1
          at_si = at_species(at_i)

          do j=1, at_neigh_list(0,at_i)
             part = at_neigh_list(j, at_i)
             call rel_pos(so_r(:,part), at_r(:,at_i), L, x)
             dist_sqr = sum( x**2 )
             if (dist_sqr < dist_min) dist_min=dist_sqr
             if ( dist_sqr .le. at_so%cut(at_si, so_species(part))**2 ) then
                if (at_so_reac(at_si,so_species(part)) % on) then
                   if (.not. so_do_reac(part) ) then
                      !eval rate
                      if ( at_so_reac(at_si, so_species(part)) % rate * DT > mtprng_rand_real1(ran_state) ) then
                         if (at_so_reac(at_si, so_species(part)) % at_exit) then
                            so_do_reac(part) = .true.
                         else
                            stop 'immediate reaction not implemented yet'
                         end if
                      end if
                   end if
                end if ! (enable_reaction)
             else
                if (at_so_reac(at_si,so_species(part)) % on) then
                   if (so_do_reac(part)) then
                      !check for neighbours!
                      too_many_atoms = .false.
                      if (count_atom_neighbours(part,.true.)>0) too_many_atoms = .true.
                      if ( ( .not. at_so_reac(at_si, so_species(part)) % two_products ) ) then
                         if (.not. too_many_atoms) then
                            so_sys % N(so_species(part)) = so_sys % N(so_species(part)) - 1
                            delta_u = u_int(so_species(part))
                            thermal = at_so_reac(at_si, so_species(part)) % thermal
                            so_species(part) = at_so_reac(at_si, so_species(part)) % product1
                            so_sys % N(so_species(part)) = so_sys % N(so_species(part)) + 1
                            delta_u = delta_u - u_int(so_species(part))
                            if (thermal) then
                               call add_kin_kick_g_so(group_list(g_i), at_i, part, delta_u)
                            end if
                            so_do_reac(part) = .false.
                         end if
                      else
                         if (.not. too_many_atoms) then
                            if (so_sys % N(0) .ge. so_sys % N_max) then
                               write(*,*) 'particle number exceeded in reac_loop'
                               stop
                            end if
                            p1si = at_so_reac(at_si, so_species(part))%product1
                            p2si = at_so_reac(at_si, so_species(part))%product2
                            so_sys % N(so_species(part)) = so_sys % N(so_species(part)) - 1
                            so_sys % N(0) = so_sys % N(0) + 1
                            so_sys % N(p1si) = so_sys % N(p1si) + 1
                            so_sys % N(p2si) = so_sys % N(p2si) + 1
                            call reac_A_C_to_B1_B2_C(group_list(g_i), at_i, part, so_sys%N(0), p1si, p2si)
                            ! The solvent particle needs to be included in the neighbour
                            ! list
                         end if
                      end if
                   end if
                end if ! (enable_reaction)
             end if
          end do
       end do
    end do
    
    if (so_sys % N(0) > total_N) then
       ! add new particles to cell lists
       do part = total_N + 1, so_sys%N(0)
          do j=1,3
             if (so_r(j,part) < shift(j)) so_r(j,part) = so_r(j,part) + L(j)
             if (so_r(j,part) >= L(j)+shift(j)) so_r(j,part) = so_r(j,part) - L(j)
          end do
          call indices(so_r(:,part), cc)
          ci = cc(1) ; cj = cc(2) ; ck = cc(3)
          if ( ( maxval( cc - N_cells ) .gt. 0) .or. ( minval( cc ) .le. 0) ) then
             write(*,*) 'particle', part, 'out of bounds'
             stop
          end if
          par_list(0,ci,cj,ck) = par_list(0,ci,cj,ck) + 1
          if (par_list(0,ci,cj,ck) .ge. max_per_cell) then
             write(*,*) 'too many particles in cell', cc, 'particle', part
             stop
          end if
          par_list(par_list(0,ci,cj,ck), ci, cj, ck) = part
       end do
       ! end add new particles to cell lists

       ! add new particles to the neighbour list
       do g_i = 1, N_groups
          do at_i = group_list(g_i) % istart , group_list(g_i) % istart + group_list(g_i) % N - 1
             do part = total_N + 1, so_sys%N(0)
                call rel_pos(so_r(:,part), at_r(:,at_i), L, x)
                dist_sqr = sum(x**2)
                neigh_sqr = at_so%neigh(at_species(at_i),so_species(part))**2
                if (dist_sqr < neigh_sqr) then
                   at_neigh_list(0,at_i) = at_neigh_list(0,at_i) + 1
                   if (at_neigh_list(0,at_i) > max_neigh) then
                      write(*,*) 'too many neighbours for atom',at_i, 'in reac_loop'
                      stop
                   end if
                   at_neigh_list(at_neigh_list(0,at_i),at_i) = part
                   so_neigh_list(0,part) = so_neigh_list(0,part) + 1
                   if (so_neigh_list(0,part) > max_so_neigh) &
                        stop 'too many neighbours for solvent in reac_loop'
                   so_neigh_list(so_neigh_list(0,part),part) = at_i
                   is_MD(part) = .true.
                end if
             end do
          end do
       end do
       ! end add new particles to the neighbour list

    end if

  end subroutine reac_loop

  !> Add the kinetic energy kin_add to a solvent particle and a group of atoms or a single atom while
  !! preserving the center of mass velocity.
  !!
  !! @param g_var Group considered.
  !! @param at_i Index of the atom considered.
  !! @param so_i Index of the solvent particle considered.
  !! @param kin_add Amount of kinetic energy to add.
  subroutine add_kin_g_so(g_var, at_i, so_i, kin_add)
    type(group_t), intent(inout) :: g_var
    integer, intent(in) :: at_i, so_i
    double precision, intent(in) :: kin_add

    double precision :: v_com(3), vr_at(3), vr_so(3), old_g_v(3)
    double precision :: kin_0, kin_r, scale
    integer :: i

    if (g_var % g_type .eq. SHAKE_G) then
       v_com = ( g_var % v * g_var % mass + so_v(:,so_i) * so_sys % mass(so_species(so_i)) ) / &
            ( g_var % mass + so_sys % mass(so_species(so_i)) )
       vr_at = g_var % v - v_com
       old_g_v = g_var % v
       kin_0 = g_var % mass * sum(g_var % v**2)
       kin_r = g_var % mass * sum(vr_at**2) * 0.5d0
    else
       v_com = ( at_v(:,at_i) * at_sys % mass(at_species(at_i)) + so_v(:,so_i) * so_sys % mass(so_species(so_i)) ) / &
            ( at_sys % mass(at_species(at_i)) + so_sys % mass(so_species(so_i)) )
       vr_at = at_v(:,at_i) - v_com
       kin_0 = at_sys % mass(at_species(at_i)) * sum(at_v(:,at_i)**2)
       kin_r = at_sys % mass(at_species(at_i)) * sum(vr_at**2) * 0.5d0
    end if
    vr_so = so_v(:,so_i) - v_com
    kin_0 = 0.5d0 * (kin_0 +  so_sys % mass(so_species(so_i)) * sum(so_v(:,so_i)**2) )
    kin_r = kin_r + so_sys % mass(so_species(so_i)) * sum(vr_so**2) * 0.5d0

    scale = sqrt( 1.d0 + kin_add/kin_r )
    write(27,*) kin_add,kin_r,scale
    vr_at = scale*vr_at
    vr_so = scale*vr_so
    if (g_var % g_type .eq. SHAKE_G) then
       g_var % v = v_com + vr_at
       do i=g_var % istart,g_var % istart + g_var % N - 1
          at_v(:,i) = at_v(:,i) + g_var % v - old_g_v
       end do
    else
       at_v(:,at_i) = v_com + vr_at
    end if
    so_v(:,so_i) = v_com + vr_so

  end subroutine add_kin_g_so

  !> Add the kinetic energy kin_add to a solvent particle and a group of atoms or a single atom while
  !! preserving the center of mass velocity.
  !!
  !! @param g_var Group considered.
  !! @param at_i Index of the atom considered.
  !! @param so_i Index of the solvent particle considered.
  !! @param kin_add Amount of kinetic energy to add.
  subroutine add_kin_kick_g_so(g_var, at_i, so_i, kin_add)
    type(group_t), intent(inout) :: g_var
    integer, intent(in) :: at_i, so_i
    double precision, intent(in) :: kin_add

    double precision :: waj(3), tot_mass, reduced_mass, raj(3), old_g_v(3)
    double precision :: kin_0, kin_r, delta
    integer :: i

    if (g_var % g_type .eq. SHAKE_G) then
       waj = so_v(:,so_i) - g_var % v
       tot_mass = g_var % mass + so_sys% mass(so_species(so_i))
       reduced_mass = g_var % mass * so_sys% mass(so_species(so_i)) / tot_mass
       call rel_pos(so_r(:,so_i),g_var% r,L,raj)
    else
       waj = so_v(:,so_i) - at_v(:,at_i)
       tot_mass = at_sys% mass(at_species(at_i)) + so_sys% mass(so_species(so_i))
       reduced_mass = at_sys% mass(at_species(at_i)) * so_sys% mass(so_species(so_i)) / tot_mass
       call rel_pos(so_r(:,so_i),at_r(:,at_i),L,raj)
    end if
    raj = raj/sqrt(sum(raj**2))
    
    delta = - sum(waj*raj) + sqrt( (sum(waj*raj))**2+ 2.d0*kin_add/reduced_mass )
    write(27,'(6e15.7)') delta, reduced_mass, tot_mass, raj

    if (g_var % g_type .eq. SHAKE_G) then
       g_var % v = g_var% v - so_sys% mass(so_species(so_i))*delta*raj/tot_mass
       do i=g_var % istart,g_var % istart + g_var % N - 1
          at_v(:,i) = at_v(:,i) - so_sys% mass(so_species(so_i))*delta*raj/tot_mass
       end do
       so_v(:,so_i) = so_v(:,so_i) + g_var% mass*delta*raj/tot_mass
    else
       at_v(:,at_i) = at_v(:,at_i) - so_sys% mass(so_species(so_i))*delta*raj/tot_mass
       so_v(:,so_i) = so_v(:,so_i) + at_sys% mass(at_species(at_i))*delta*raj/tot_mass
    end if

  end subroutine add_kin_kick_g_so

  !> This function counts the neighbours of a given solvent particle.
  !!
  !! It takes the neighbour candidates from so_neigh_list and compares the distance with the
  !! cut-off of the appropriate LJ interaction parameter.
  !! @param so_i The solvent particle to consider.
  !! @param one_is_enough If set to .true., the count stop whenever the first neighbouring atom
  !! is found.
  function count_atom_neighbours(so_i, one_is_enough)
    implicit none
    integer :: count_atom_neighbours
    integer, intent(in) :: so_i
    logical, optional, intent(in) :: one_is_enough
    
    integer :: at_i, i
    double precision :: x(3), d_sqr
    logical :: flag

    if (present(one_is_enough)) then
       flag = one_is_enough
    else
       flag = .false.
    end if

    count_atom_neighbours = 0
    do i=1,so_neigh_list(0,so_i)
       at_i = so_neigh_list(i,so_i)
       call rel_pos(so_r(:,so_i), at_r(:,at_i), L, x)
       d_sqr = sum(x**2)
       if ( d_sqr <= 1.d0*at_so % cut(at_species(at_i), so_species(so_i))**2 ) then
          if (flag) then
             count_atom_neighbours = 1
             return
          else
             count_atom_neighbours = count_atom_neighbours + 1
          end if
       end if
    end do

  end function count_atom_neighbours

  !> Converts a A solvent particle to two (B1 and B2) solvent particles upon activation by
  !! a catalytic compound C.
  !!
  !! @param g_var The group to which C belongs.
  !! @param at_i The index of the catalytic atom.
  !! @param so_1 The first solvent particle (corresponds to A and B1).
  !! @param so_2 The second solvent particle (corresponds to B2).
  !! @param so_s1 The species of the first product solvent particle B1.
  !! @param so_s2 The species of the second product solvent particle B2.
  subroutine reac_A_C_to_B1_B2_C(g_var, at_i, so_1, so_2, so_s1, so_s2)
    type(group_t), intent(inout) :: g_var
    integer, intent(in) :: at_i, so_1, so_2, so_s1, so_s2

    double precision :: vel(3), rel_v, x(3)
    double precision :: A_mass, C_mass, B1_mass, B2_mass
    integer :: i

    A_mass = so_sys % mass(so_species(so_1))
    B1_mass = so_sys % mass(so_s1)
    B2_mass = so_sys % mass(so_s2)

    if (g_var % g_type .eq. SHAKE_G) then
       C_mass = g_var % mass
       vel = ( C_mass * g_var % v + A_mass * so_v(:,so_1)) / &
            (C_mass + A_mass)
       rel_v = sqrt( &
            A_mass*C_mass/(A_mass+C_mass)*(B1_mass+B2_mass)/(B1_mass*B2_mass) * &
            sum( (g_var % v - so_v(:,so_1))**2 ) )
    else
       C_mass = at_sys%mass(at_species(at_i))
       vel = (C_mass * at_v(:,at_i) + A_mass * so_v(:,so_1)) / &
            ( C_mass + A_mass )
       rel_v = sqrt( &
            A_mass*C_mass/(A_mass+C_mass)*(B1_mass+B2_mass)/(B1_mass*B2_mass) * &
            sum( (at_v(:,at_i) - so_v(:,so_1))**2 ) )
    end if

    if (g_var % g_type .eq. SHAKE_G) then
       do i=g_var % istart,g_var % istart + g_var % N - 1
          at_v(:,i) = at_v(:,i) - g_var % v + vel
       end do
       g_var % v = vel
    else
       at_v(:,at_i) = vel
    end if
    x = rand_sphere()
    so_v(:,so_1) = vel + rel_v * B2_mass/(B1_mass+B2_mass) * x
    so_v(:,so_2) = vel - rel_v * B1_mass/(B1_mass+B2_mass) * x
    so_r(:,so_2) = so_r(:,so_1)

    so_species(so_1) = so_s1
    so_species(so_2) = so_s2

  end subroutine reac_A_C_to_B1_B2_C

  !> Converts 2 B particles to a single A particle in the reaction \f$ 2B \to A \f$.
  !!
  !! @param so1 Indice de la première particule B.
  !! @param so2 Indice de la deuxième particule B.
  !! @param s_product Species to convert so1 and so2 into.
  !! @param deltak The excess kinetic energy from the recombination.
  subroutine reac_2B_to_A(so1,so2,s_product,deltak)
    implicit none
    integer, intent(in) :: so1, so2, s_product
    double precision, intent(out) :: deltak

    double precision :: wbb(3), mubb
    double precision :: m1, m2

    wbb = so_v(:,so1) - so_v(:,so2)
    m1 = so_sys%mass(so_species(so1))
    m2 = so_sys%mass(so_species(so2))
    mubb = m1*m2/(m1+m2)
    deltak = 0.5d0 * mubb * sum(wbb**2)

    so_v(:,so1) = ( m1*so_v(:,so1) + m2*so_v(:,so2) ) / (m1+m2)
    so_species(so1) = s_product
    so_r(:,so1) = (m1*so_r(:,so1)+m2*so_r(:,so2)) / (m1+m2)

  end subroutine reac_2B_to_A

  !> Removes a particle from the system.
  !!
  !! This subroutine removes particle del_i from the system. It replaces its data by the data
  !! of the last particle in the system, then decreases the number of particles by one unit.
  !! @param del_i The index of the particle to remove.
  !! @param cc The three indices of the cell to which the particle belongs.
  subroutine del_particle(del_i,cc)
    integer, intent(in) :: del_i, cc(3)

    integer :: i, j, nlast, at_i, at_neigh

    nlast = so_sys%N(0)

    so_r(:,del_i) = so_r(:,nlast)
    so_v(:,del_i) = so_v(:,nlast)
    so_f1(:,del_i) = so_f1(:,nlast)
    so_f2(:,del_i) = so_f2(:,nlast)
    so_r_neigh(:,del_i) = so_r_neigh(:,nlast)
    so_species(del_i) = so_species(nlast)
    is_MD(del_i) = is_MD(nlast)
    N_MD(del_i) = N_MD(nlast)

    so_neigh_list(:,del_i) = so_neigh_list(:,nlast)

    ! loop over atoms that are close to del_i
    do at_neigh = 1,so_neigh_list(0,del_i)
       at_i = so_neigh_list(at_neigh,del_i)
       ! loop over the neighbours of those atoms
       do i=1,at_neigh_list(0,at_i)
          j = at_neigh_list(i,at_i)
          ! if our particle to delete is found, remove it
          if (j.eq.del_i) then
             ! replace "j" by the last neighbour of the atom
             at_neigh_list(i,del_i) = at_neigh_list(at_neigh_list(0,del_i),del_i)
             at_neigh_list(0,del_i) = at_neigh_list(0,del_i) - 1
          end if
          exit
       end do
    end do

    ! remove del_i from the particle list, and replace nlast by del_i when found.
    ! cell indices should be provided.
    do i=1,par_list(0,cc(1),cc(2),cc(3))
       if (par_list(i,cc(1),cc(2),cc(3)).eq.del_i) then
          ! remove del_i from list and decrease count by 1
          par_list(i,cc(1),cc(2),cc(3)) = &
               par_list(par_list(0,cc(1),cc(2),cc(3)),cc(1),cc(2),cc(3))
          par_list(0,cc(1),cc(2),cc(3)) = par_list(0,cc(1),cc(2),cc(3)) - 1
          ! get out of loop
          exit
       end if
    end do

  end subroutine del_particle

end module MD
