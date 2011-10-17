
!> This module holds the variables for a MPCD simulation.
!!
!! Most variables are global. The MPCD particles variables are prefixed by
!! "so_" (for solvent) so as not to be mixed up with the MD module.

module MPCD
  use sys
  use mtprng
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
    double precision :: x(3), t_factor

    t_factor = sqrt(3.d0*temperature)

    do i=1,so_sys%N(0)
       x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ; x(3) = mtprng_rand_real1(ran_state) ; 
       so_r(:,i) = x*L
       x(1) = mtprng_rand_real1(ran_state) ; x(2) = mtprng_rand_real1(ran_state) ; x(3) = mtprng_rand_real1(ran_state) ; 
       x = x-0.5d0
       so_v(:,i) = x*2.d0 * t_factor/sqrt(so_sys%mass(so_species(i)))
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
       n(1) = mtprng_rand_real1(ran_state)
       n(2) = mtprng_rand_real1(ran_state)
       s = n(1)**2 + n(2)**2
       if ( s<1.d0 ) s_lt_one = .true.
    end do
    alpha = 2.d0 * sqrt(1.d0 - s)
    n(1) = n(1)*alpha
    n(2) = n(2)*alpha
    n(3) = 1.d0 - 2.d0*s
  end function rand_sphere

end module MPCD
