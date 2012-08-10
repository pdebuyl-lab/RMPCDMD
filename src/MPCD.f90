
!> This module holds the variables for a MPCD simulation.
!!
!! Most variables are global. The MPCD particles variables are prefixed by
!! "so_" (for solvent) so as not to be mixed up with the MD module.

module MPCD
  use sys
  use mtprng
  use volume_reaction
  implicit none

  !> Value of Pi computed via N[Pi,35] in Mathematica
  double precision, parameter :: PI = 3.1415926535897932384626433832795029d0
  !> The size of the list of particles for each MPCD cell
  integer, parameter :: max_per_cell=128

  ! MPCD box variables
  !> The number of cells in each direction
  integer :: N_cells(3)
  !> The list of particles belonging to each cell
  integer, allocatable :: par_list(:,:,:,:)
  !> The list of cell-wise center of mass momenta.
  !! The leftmost index is
  !! - 1:3 : components of the mometum.
  !! - 4   : total mass in the cell
  double precision, allocatable :: Vcom(:,:,:,:)
  !> rotation matrix for each cell.
  double precision, allocatable :: omega(:,:,:,:,:)
  !> cell_active controls wether a cell reacts or not.
  logical, allocatable :: cell_active(:,:,:)
  !> cell_collide controls wether a cell collides fully or partially.
  logical, allocatable :: cell_collide(:,:,:)
  !> The length of a cell.
  double precision :: a
  !> The inverse length of a cell.
  double precision :: oo_a
  !> Size of the system in each direction. The box is [0:L(d)] in the direction d.
  double precision :: L(3)
  !> Inverse of L
  double precision :: oo_L(3)
  !> The MPCD streaming step.
  double precision :: tau
  !> The MPCD time step. May be different than tau if streaming is performed for partial steps.
  double precision :: MPCD_tau
  !> Switch to activate gridshifting.
  logical :: do_shifting
  !> Value of the shift. Is kept between -a/2 and a/2.
  double precision :: shift(3)

  ! MPCD particles variables
  !> Position of the MPCD solvent, dimension (3,N)
  double precision, allocatable :: so_r(:,:)
  !> Velocity of the MPCD solvent, dimension (3,N)
  double precision, allocatable :: so_v(:,:)
  !> Force for the MPCD solvent, copy 1. It is used in alternance with copy 2.
  double precision, allocatable, target :: so_f1(:,:)
  !> Force for the MPCD solvent, copy 2. It is used in alternance with copy 1.
  double precision, allocatable, target :: so_f2(:,:)
  !> Pointer to the force for the MPCD solvent. so_f points to the last computed force.
  double precision, pointer :: so_f(:,:)
  !> Pointer to the force for the MPCD solvent. so_f_old points to the before-to-last computed force.
  double precision, pointer :: so_f_old(:,:)
  !> Pointer to the force for the MPCD solvent. so_f_temp holds the pointer while switching so_f and so_f_old.
  double precision, pointer :: so_f_temp(:,:)
  !> List of solvent particles that are neighbour to colloid particles.
  double precision, allocatable :: so_r_neigh(:,:)
  !> The species of each solvent particle.
  integer, allocatable :: so_species(:)
  !> The per-species internal energy. Allows to consider exo- and endo-thermic reactions.
  double precision, allocatable :: u_int(:)
  !> Flag that indicates that a particle is found on this CPU.
  logical(kind=1), allocatable :: is_local(:)
  !> Flag that indicates that a particle exists: it has not been "destroyed" by a chemical reaction.
  logical(kind=1), allocatable :: exists(:)
  !> Flag that indicates that the force on a particle is non-zero and that MD stepping should be 
  !! used instead of simple streaming.
  logical(kind=1), allocatable :: is_MD(:)
  !> Number of MD steps that should be taken into account when performing a streaming.
  integer, allocatable :: N_MD(:)
  !> Stores the maximum number of MD steps that can be done before streaming the non-MD MPCD particles.
  integer :: N_MD_max


  ! Chemical reaction variables

  !> List of volume chemical reactions.
  type(vol_reac_t), allocatable :: reaction_list(:)
  !> Number of volume chemical reactions.
  integer :: N_reactions


  !> Information on the solvent system, based on sys_t from the sys group.
  type(sys_t) :: so_sys

  !> State of the random number generator.
  type(mtprng_state), save :: ran_state

contains
  
  !> Reads the cell size and allocates all MPCD solvent arrays.
  !! Also reads the number of particles for each species and the internal energy
  !! of each species.
  !! @param CF Configuration file.
  subroutine config_MPCD(CF)
    use ParseText
    implicit none
    type(PTo), intent(inout) :: CF

    character(len=10) :: temp_name

    integer :: i,j

    N_cells = PTread_ivec(CF,'N_cells',3)
    a = PTread_d(CF,'cell_unit')
    oo_a = 1.d0/a
    L = N_cells * a
    oo_L = 1.d0 / L

    allocate(par_list(0:max_per_cell,N_cells(1),N_cells(2),N_cells(3)))
    allocate(Vcom(4,N_cells(1),N_cells(2),N_cells(3)))
    allocate(omega(3,3,N_cells(1),N_cells(2),N_cells(3)))
    allocate(cell_active(N_cells(1),N_cells(2),N_cells(3)))
    cell_active = .true.
    allocate(cell_collide(N_cells(1),N_cells(2),N_cells(3)))
    cell_collide = .true.
    
    allocate(so_r(3,so_sys%N_max))
    allocate(so_v(3,so_sys%N_max))
    allocate(so_f1(3,so_sys%N_max))
    allocate(so_f2(3,so_sys%N_max))
    allocate(so_r_neigh(3,so_sys%N_max))
    allocate(so_species(so_sys%N_max))
    allocate(is_MD(so_sys%N_max))
    allocate(N_MD(so_sys%N_max))

    allocate(u_int(so_sys%N_species))

    j=1
    do i=1,so_sys%N_species
       so_species(j:j-1+so_sys%N(i)) = i
       j = j+so_sys%N(i)
    end do

    do i=1,so_sys%N_species
       write(temp_name,'(a8,i02.2)') 'so_u_int',i
       u_int(i) = PTread_d(CF,trim(temp_name))
    end do

  end subroutine config_MPCD

  !> Initializes a solvent with random positions in the simulation box, regardless
  !! of any obstacle.
  !! Velocities taken from a flat velocity profile with requested temperature.
  !! @param temperature The temperature given to the solvent.
  subroutine homogeneous_solvent(temperature)
    double precision, intent(in) :: temperature
    integer :: i, Nloop
    double precision :: x(3), t_factor, tot_m

    t_factor = sqrt(3.d0*temperature)

    do i=1,so_sys%N(0)
       x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ; x(3) = mtprng_rand_real1(ran_state) ; 
       so_r(:,i) = x*L
       x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ; x(3) = mtprng_rand_real1(ran_state) ; 
       x = x-0.5d0
       so_v(:,i) = x*2.d0 * t_factor/sqrt(so_sys%mass(so_species(i)))
    end do

    x = 0.d0
    tot_m = 0.d0
    do i=1,so_sys%N(0)
       x = x + so_sys%mass(so_species(i)) * so_v(:,i)
       tot_m = tot_m + so_sys%mass(so_species(i))
    end do
    x = x/tot_m

    do i=1,so_sys%N(0)
       so_v(:,i) = so_v(:,i) - x
    end do

    Nloop = 1
    do i=1,so_sys%N_species
       so_species(Nloop:Nloop-1+so_sys%N(i)) = i
       Nloop = Nloop + so_sys%N(i)
    end do
   

  end subroutine homogeneous_solvent

  !> Computes the cell indices for a given position.
  !!
  !! The indices are 1-based. The routine returns a 3 elements vector of indices.
  !!
  !! @param x Position for which one wants to locate a cell.
  !! @param x_index The 3 indices for the cell.
  subroutine indices(x,x_index)
    implicit none
    double precision, intent(in) :: x(3)
    integer, intent(out) :: x_index(3)

    x_index = floor( (x - shift)* oo_a ) + 1

  end subroutine indices

  !> Returns the geometric center of a cell.
  !!
  !! @param ci x-index
  !! @param cj y-index
  !! @param ck z-index
  !! @return cell_center The center of the cell.
  function cell_center(ci,cj,ck)
    implicit none
    integer, intent(in) :: ci,cj,ck
    double precision :: cell_center(3)

    cell_center = shift + (dble( (/ ci, cj, ck /) ) -0.5d0) * a

  end function cell_center

  !> Builds "par_list" by filling for each cell the list of solvent particles
  !! that are explicitly found in that cell.
  subroutine place_in_cells
    
    integer :: i, ci, cj, ck
    integer :: cc(3)

    par_list = 0

    do i=1,so_sys%N(0)
       call indices(so_r(:,i), cc)
       ci = cc(1) ; cj = cc(2) ; ck = cc(3)
       
       if ( ( maxval( (cc-1)/N_cells ) .ge. 1) .or. ( minval( (cc-1)/N_cells ) .lt. 0) ) then
          write(*,*) 'particle', i, 'out of bounds'
       end if

       par_list(0,ci,cj,ck) = par_list(0,ci,cj,ck) + 1
       if (par_list(0,ci,cj,ck) .ge. max_per_cell) then
          write(*,*) 'too many particles in cell', cc, 'particle', i
          stop
       end if
       par_list(par_list(0,ci,cj,ck), ci, cj, ck) = i

    end do

  end subroutine place_in_cells
  
  !> Computes the center of mass velocity of each cell.
  !!
  !! The first three elements of Vcom(:,ci,cj,ck) represent the momenta while
  !! the fourth is the mass of that cell.
  subroutine compute_v_com

    integer :: i,j, ci, cj, ck
    double precision :: vv(4)

    do ck=1,N_cells(3)
       do cj=1,N_cells(2)
          do ci=1,N_cells(1)
             vv = 0.d0
             do i=1,par_list(0,ci,cj,ck)
                j = par_list(i,ci,cj,ck)
                vv(1:3) = vv(1:3) + so_v(:,j)*so_sys%mass(so_species(j))
                vv(4) = vv(4) + so_sys%mass(so_species(j))
             end do
             Vcom(:,ci,cj,ck) = vv
          end do
       end do
    end do
  end subroutine compute_v_com

  !> Generates, for each cell, a rotation matrix whose angle is
  !! Pi/2 and whose orientation vector is chosen at random.
  subroutine generate_omega
    
    integer :: ci,cj,ck
    double precision :: n(3)

    do ck=1,N_cells(3)
       do cj=1,N_cells(2)
          do ci=1,N_cells(1)
             n = rand_sphere()
             omega(:,:,ci,cj,ck) = &
                  reshape( (/ &
                  n(1)**2, n(1)*n(2) + n(3), n(1)*n(3) - n(2) ,&
                  n(2)*n(1) - n(3) , n(2)**2 , n(2)*n(3) + n(1),&
                  n(3)*n(1) + n(2), n(3)*n(2) - n(1), n(3)**2 &
                  /), (/3, 3/))
          end do
       end do
    end do
  end subroutine generate_omega

  !> Collides the particles whithin each cell, according to their precomputed
  !! rotation matrix.
  subroutine simple_MPCD_step

    integer :: i, j
    integer :: ci,cj,ck
    double precision :: om(3,3),vv(4)

    do ck=1,N_cells(3)
       do cj=1,N_cells(2)
          do ci=1,N_cells(1)
             vv = Vcom(:,ci,cj,ck)
             if (vv(4)>0.d0) vv(1:3) = vv(1:3)/vv(4)
             om = omega(:,:,ci,cj,ck)
             do i=1,par_list(0,ci,cj,ck)
                j = par_list(i,ci,cj,ck)
                so_v(:,j) = so_v(:,j) - vv(1:3)
                so_v(:,j) = matmul(om,so_v(:,j))
                so_v(:,j) = so_v(:,j) + vv(1:3)
             end do
          end do
       end do
    end do

  end subroutine simple_MPCD_step

  !> Collides the particles whithin each cell, according to their precomputed
  !! rotation matrix. Applies chemical reactions if appropriate.
  subroutine chem_MPCD_step

    integer :: i, j
    integer :: ci,cj,ck
    double precision :: om(3,3),vv(4)

    do ck=1,N_cells(3)
       do cj=1,N_cells(2)
          do ci=1,N_cells(1)
             vv = Vcom(:,ci,cj,ck)
             if (vv(4)>0.d0) vv(1:3) = vv(1:3)/vv(4)
             om = omega(:,:,ci,cj,ck)
             do i=1,par_list(0,ci,cj,ck)
                j = par_list(i,ci,cj,ck)
                so_v(:,j) = so_v(:,j) - vv(1:3)
                so_v(:,j) = matmul(om,so_v(:,j))
                so_v(:,j) = so_v(:,j) + vv(1:3)
             end do
             if (N_reactions > 0) then
                if (cell_active(ci,cj,ck)) &
                     call cell_reaction(ci,cj,ck,N_reactions,reaction_list,MPCD_tau)
             end if
          end do
       end do
    end do

  end subroutine chem_MPCD_step

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

  !> Advances all particles according to their velocity, for 
  !! a time of tau.
  subroutine MPCD_stream
    integer :: i

    do i=1,so_sys%N(0)
       if (.not. is_MD(i)) so_r(:,i) = so_r(:,i) + so_v(:,i) * tau
    end do
    
  end subroutine MPCD_stream

  !> Returns a list of indices from solvent particles in cells that are close enough to
  !! the position x0.
  !! @param x0 the position around which we look for solvent particles.
  !! @param radius the search radius.
  !! @param list a list of indices, that is allocated in this routine.
  subroutine list_idx_from_x0(x0, radius, list)
    implicit none
    double precision, intent(in) :: x0(3)
    double precision, intent(in) :: radius
    integer, intent(inout), allocatable :: list(:)

    integer :: cc(3)
    integer :: il, jl, kl, iu, ju, ku, mi, mj, mk
    integer :: i, j, k
    integer :: in_cell, part, idx

    integer :: list_len

    call indices(x0,cc)

    il = cc(1) - ceiling( radius/a ) - 1
    jl = cc(2) - ceiling( radius/a ) - 1
    kl = cc(3) - ceiling( radius/a ) - 1

    iu = cc(1) + ceiling( radius/a ) - 1
    ju = cc(2) + ceiling( radius/a ) - 1
    ku = cc(3) + ceiling( radius/a ) - 1
    
    list_len = 0    
    do k = kl, ku
       do j = jl, ju
          do i=il,iu
             mi = modulo(i-1,N_cells(1)) + 1 ; mj = modulo(j-1,N_cells(2)) + 1 ; mk = modulo(k-1,N_cells(3)) + 1 ; 
             list_len = list_len + par_list(0, mi,mj,mk)
          end do
       end do
    end do

    allocate(list(list_len))

    idx=1
    do k = kl, ku
       do j = jl, ju
          do i=il,iu
             mi = modulo(i-1,N_cells(1)) + 1 ; mj = modulo(j-1,N_cells(2)) + 1 ; mk = modulo(k-1,N_cells(3)) + 1 ; 
             do in_cell=1,par_list(0, mi,mj,mk)
                part = par_list(in_cell, mi,mj,mk)
                list(idx) = part
                idx = idx + 1
             end do
          end do
       end do
    end do

  end subroutine list_idx_from_x0

  !> Returns a unit vector whose direction is randomly chosen on a sphere.
  !!
  !! @todo add refs for algo.
  !!
  !! @return rand_sphere A 3-d random vector.
  function rand_sphere() result(n)
    double precision :: n(3)

    logical :: s_lt_one
    double precision :: s, alpha

    s_lt_one = .false.
    do while (.not. s_lt_one)
       n(1) = 2.d0*mtprng_rand_real1(ran_state)-1.d0
       n(2) = 2.d0*mtprng_rand_real1(ran_state)-1.d0
       s = n(1)**2 + n(2)**2
       if ( s<1.d0 ) s_lt_one = .true.
    end do
    alpha = 2.d0 * sqrt(1.d0 - s)
    n(1) = n(1)*alpha
    n(2) = n(2)*alpha
    n(3) = 1.d0 - 2.d0*s
  end function rand_sphere

  !> Adds a specified amount of kinetic energy in a cell.
  !!
  !! @param cc The cell indices.
  !! @param kin_add The amount of energy to add.
  subroutine add_energy(cc, kin_add)
    implicit none
    integer, intent(in) :: cc(3)
    double precision, intent(in) :: kin_add

    double precision :: gamma, kin, v(3)
    integer :: i,j

    v = Vcom(1:3,cc(1),cc(2),cc(3)) / Vcom(4,cc(1),cc(2),cc(3))

    kin = 0.d0
    do i=1,par_list(0,cc(1),cc(2),cc(3))
       j = par_list(i,cc(1),cc(2),cc(3))
       kin = kin + so_sys%mass(so_species(j))*sum( (so_v(:,j) - v)**2 )
    end do

    gamma = sqrt( 1.d0 + ( 2.d0*kin_add/kin ) )

    do i=1,par_list(0,cc(1),cc(2),cc(3))
       j = par_list(i,cc(1),cc(2),cc(3))
       so_v(:,j) = gamma*(so_v(:,j) - v) + v
    end do

  end subroutine add_energy

  !> This routine performs the A to B bulk chemical reaction.
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param vr_var The vol_reac_t variable.
  subroutine vol_A_to_B(ci,cj,ck,vr_var)
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
          so_sys % N(vr_var % reac(1)) = so_sys % N(vr_var % reac(1)) - 1
          so_sys % N(vr_var % prod(1)) = so_sys % N(vr_var % prod(1)) + 1
          found = .true.
          exit
       end if
    end do

    if (.not. found) stop 'reactant not found in vol_A_to_B'

  end subroutine vol_A_to_B

  !> This routine performs the A to B bulk endothermic chemical reaction.
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param vr_var The vol_reac_t variable.
  subroutine vol_A_to_B_endo(ci,cj,ck,vr_var)
    implicit none
    integer, intent(in) :: ci,cj,ck
    type(vol_reac_t), intent(in) :: vr_var

    integer :: i, part, i_found
    logical :: found
    double precision :: deltau, kin, v(3), gamma

    v = Vcom(1:3,ci,cj,ck) / Vcom(4,ci,cj,ck)
    deltau = u_int( vr_var % prod(1) ) - u_int( vr_var % reac(1) )
    found = .false.
    kin = 0.d0
    do i=1,par_list(0,ci,cj,ck)
       part = par_list(i,ci,cj,ck)
       if (so_species(part) .eq. vr_var % reac(1) .and. &
            .not. found) then
          i_found = part
          found = .true.
       end if
       kin = kin + 0.5d0 * so_sys%mass(so_species(part)) * &
            sum( (so_v(:,part) - v)**2 )
    end do

    if (kin .gt. deltau .and. found) then
       so_species(i_found) = vr_var % prod(1)
       so_sys % N(vr_var % reac(1)) = so_sys % N(vr_var % reac(1)) - 1
       so_sys % N(vr_var % prod(1)) = so_sys % N(vr_var % prod(1)) + 1
       gamma = sqrt( 1.d0 - deltau/kin )
       do i=1,par_list(0,ci,cj,ck)
          part = par_list(i,ci,cj,ck)
          so_v(:,part) = v + gamma*(so_v(:,part)-v)
       end do
    end if

    if (.not. found) stop 'reactant not found in vol_A_to_B'

  end subroutine vol_A_to_B_endo

  !> This routine performs the A to B bulk chemical reaction.
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param vr_var The vol_reac_t variable.
  subroutine vol_A_B_to_C(ci,cj,ck,vr_var)
    implicit none
    integer, intent(in) :: ci,cj,ck
    type(vol_reac_t), intent(in) :: vr_var

    integer :: i, part, cc(3)
    logical :: found_A, found_B
    integer :: i_A, i_B
    double precision :: wbb(3), mubb, deltak
    double precision :: m1, m2

    cc(1) = ci ; cc(2) = cj ; cc(3) = ck

    found_A = .false.
    found_B = .false.
    do i=1,par_list(0,ci,cj,ck)
       part = par_list(i,ci,cj,ck)
       if (so_species(part) .eq. vr_var % reac(1) .and. .not. found_A) then
          i_A = part
          found_A = .true.
       else if (so_species(part) .eq. vr_var % reac(2) .and. .not. found_B) then
          i_B = part
          found_B = .true.
       end if
       if (found_A .and. found_B) exit
    end do

    if (i_A.eq.i_B) stop 'i_A and i_B equal in vol_A_B_to_C'

    wbb = so_v(:,i_A) - so_v(:,i_B)
    m1 = so_sys%mass(so_species(i_A))
    m2 = so_sys%mass(so_species(i_B))
    mubb = m1*m2/(m1+m2)
    deltak = 0.5d0 * mubb * sum(wbb**2)

    so_v(:,i_A) = ( m1*so_v(:,i_A) + m2*so_v(:,i_B) ) / (m1+m2)
    so_species(i_A) = vr_var % prod(1)

    call del_particle_noMD(i_B, cc)

    so_sys % N(vr_var % reac(1)) = so_sys % N(vr_var % reac(1)) - 1
    so_sys % N(vr_var % reac(1)) = so_sys % N(vr_var % reac(2)) - 1
    so_sys % N(vr_var % prod(1)) = so_sys % N(vr_var % prod(1)) + 1
    so_sys % N(0) = so_sys % N(0) - 1

    call add_energy(cc,deltak)

    if (.not. found_A .or. .not. found_B) stop 'reactant not found in vol_A_B_to_C'

  end subroutine vol_A_B_to_C

  !> This function computes the combinatorial term h_\mu^\xi defined
  !! in Rohlf et al, Comp. Phys. Comm. 179, pp 132 (2008), Eq. (14).
  !! @param ci Index i of the cell in which the particle belongs.
  !! @param cj Index j of the cell in which the particle belongs.
  !! @param ck Index k of the cell in which the particle belongs.
  !! @param vr_var The vol_reac_t variable.
  function combi(ci,cj,ck,vr_var)
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
      do i=n-nu+1, n
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
    implicit none
    integer, intent(in) :: ci,cj,ck, n_reac
    type(vol_reac_t), intent(in) :: reac_list(:)
    double precision, intent(in) :: interval

    integer :: i, idx_reac, idx_kind
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

    idx_kind = reac_list(idx_reac) % r_kind

    if (idx_kind .eq. AtoB_R) then
       call vol_A_to_B(ci,cj,ck,reac_list(idx_reac))
    else if (idx_kind .eq. AtoB_endothermic_R) then
       call vol_A_to_B_endo(ci,cj,ck,reac_list(idx_reac))
    else if (idx_kind .eq. A_B_to_C_R) then
       if (par_list(0,ci,cj,ck) .gt. 2) call vol_A_B_to_C(ci,cj,ck,reac_list(idx_reac))
    else
       write(*,*) 'unknown reaction type in cell_reaction'
    end if

  deallocate(amu)

  end subroutine cell_reaction

  !> Removes a particle from the system.
  !!
  !! This subroutine removes particle del_i from the system. It replaces its data by the data
  !! of the last particle in the system, then decreases the number of particles by one unit.
  !! @param del_i The index of the particle to remove.
  !! @param cc The three indices of the cell to which the particle belongs.
  subroutine del_particle_noMD(del_i,cc)
    integer, intent(in) :: del_i, cc(3)

    integer :: i, j, nlast, at_i, at_neigh
    integer :: lastcc(3)

    nlast = so_sys%N(0)

    so_r(:,del_i) = so_r(:,nlast)
    so_v(:,del_i) = so_v(:,nlast)
    so_f1(:,del_i) = so_f1(:,nlast)
    so_f2(:,del_i) = so_f2(:,nlast)
    so_r_neigh(:,del_i) = so_r_neigh(:,nlast)
    so_species(del_i) = so_species(nlast)
    is_MD(del_i) = is_MD(nlast)
    N_MD(del_i) = N_MD(nlast)

    ! remove del_i from the particle list
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

    ! seek nlast in par_list and replace by del_i
    call indices(so_r(:,nlast),lastcc)
    do i=1,par_list(0,lastcc(1),lastcc(2),lastcc(3))
       if (par_list(i,lastcc(1),lastcc(2),lastcc(3)).eq.nlast) then
          ! remove nlast from list and decrease count by 1
          par_list(i,lastcc(1),lastcc(2),lastcc(3)) = del_i
          ! get out of loop
          exit
       end if
    end do

  end subroutine del_particle_noMD

end module MPCD
