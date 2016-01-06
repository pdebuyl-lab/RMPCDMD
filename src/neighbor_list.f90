module neighbor_list
  use common
  use hilbert
  use cell_system
  use particle_system
  use interaction
  implicit none
  private

  public :: neighbor_list_t, compute_force
  public :: compute_force_n2

  type neighbor_list_t
     integer :: Npoints
     integer :: Nmax
     integer, allocatable :: list(:, :)
     integer, allocatable :: n(:)
     integer, allocatable :: stencil(:,:)
     type(timer_t) :: time_update, time_force
   contains
     procedure :: init
     procedure :: update_list
     procedure :: make_stencil
  end type neighbor_list_t

contains

  subroutine init(this, Npoints, Nmax)
    class(neighbor_list_t), intent(out) :: this
    integer, intent(in) :: Npoints
    integer, intent(in) :: Nmax

    this% Nmax = Nmax
    this% Npoints = Npoints
    allocate(this% list(Nmax, Npoints))
    allocate(this% n(Npoints))
    this% list = 0
    this% n = 0

    call this%time_update%init('neigh list update')
    call this%time_force%init('neigh list force')

  end subroutine init

  pure function next_cell(p) result(cell)
    integer, intent(in) :: p(3)
    integer :: cell(3)

    integer :: i

    cell = p
    cell(1) = cell(1) + 1

    do i = 1, 2
       if (cell(i) .eq. 3) then
          cell(i) = -2
          cell(i+1) = cell(i+1) + 1
       end if
    end do

  end function next_cell

  !! Update the neighbor list
  !!
  !! Loop over the stencil for neighboring cells.
  subroutine update_list(this, system1, system2, radius, cells)
    class(neighbor_list_t), intent(inout) :: this
    type(particle_system_t), intent(in) :: system1
    type(particle_system_t), intent(in) :: system2
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: radius

    integer :: cell(3), M(3), actual_cell(3)
    integer :: neigh_idx, i, cell_i, cell_n, cell_start, list_idx
    integer :: j, stencil_size
    double precision :: x(3), y(3), L(3)
    double precision :: rsq, radiussq

    L = cells% L * cells% a
    radiussq = radius**2
    M = cells% M
    stencil_size = size(this% stencil, dim=2)

    call this%time_update%tic()
    !$omp parallel do private(x, cell, list_idx, j, actual_cell, neigh_idx, cell_n, cell_start, cell_i, y, rsq)
    do i = 1, system1% Nmax
       x = system1% pos(:, i)
       cell = floor( (x - cells% origin) / cells% a ) + 1
       list_idx = 0

       stencil: do j = 1, stencil_size
          ! 0-based index for hilbert routine
          actual_cell = modulo(cell + this% stencil(:,j) - 1 , cells% L)
          neigh_idx = compact_p_to_h( actual_cell, M ) + 1
          cell_n = cells% cell_count(neigh_idx)
          cell_start = cells% cell_start(neigh_idx)

          do cell_i = cell_start, cell_start + cell_n - 1
             y = system2% pos(:, cell_i)
             rsq = sum(rel_pos(x, y, L)**2)

             if (rsq .lt. radiussq) then
                list_idx = list_idx + 1
                if (list_idx .gt. this% Nmax) then
                   error stop 'maximum of neighbor list reached'
                end if
                this% list(list_idx, i) = cell_i
             end if
          end do
       end do stencil

       this% n(i) = list_idx

    end do
    call this%time_update%tac()

  end subroutine update_list

  function compute_force(ps1, ps2, n_list, L, lj_params) result(e)
    type(particle_system_t), intent(inout) :: ps1
    type(particle_system_t), intent(inout) :: ps2
    type(neighbor_list_t), intent(inout) :: n_list
    double precision, intent(in) :: L(3)
    type(lj_params_t), intent(in) :: lj_params
    double precision :: e

    integer :: i, j
    integer :: idx
    integer :: n_species1, n_species2
    double precision :: x(3), d(3), f(3), f1(3)
    integer :: s1, s2
    double precision :: r_sq

    n_species1 = ps1% n_species
    n_species2 = ps2% n_species

    e = 0

    call n_list%time_force%tic()
    !$omp parallel do private(x, s1, f1, j, idx, s2, d, r_sq, f) &
    !$omp& reduction(+:e)
    do i = 1, ps1% Nmax
       if (ps1% species(i) <= 0) continue
       x = ps1% pos(:, i)
       s1 = ps1% species(i)
       f1 = 0
       do j = 1, n_list% n(i)
          idx = n_list% list(j, i)
          s2 = ps2% species(idx)
          if (s2 <= 0) continue
          d = rel_pos(x, ps2% pos(:, idx), L)
          r_sq = sum(d**2)
          if ( r_sq <= lj_params% cut_sq(s1, s2) ) then
             f = lj_force(d, r_sq, lj_params% epsilon(s1, s2), lj_params% sigma(s1, s2))
             e = e + lj_energy(r_sq, lj_params% epsilon(s1, s2), lj_params% sigma(s1, s2))
             f1 = f1 + f
             ps2% force(:, idx) = ps2% force(:, idx) - f
          end if
       end do
       ps1% force(:, i) = ps1% force(:, i) + f1
    end do
    call n_list%time_force%tac()

  end function compute_force

  function compute_force_n2(ps, L, lj_params) result(e)
    type(particle_system_t), intent(inout) :: ps
    double precision, intent(in) :: L(3)
    type(lj_params_t), intent(in) :: lj_params
    double precision :: e

    integer :: i, j
    integer :: n_species
    double precision :: x(3), d(3), f(3)
    integer :: s1, s2
    double precision :: r_sq

    n_species = ps% n_species

    e = 0

    call ps%time_self_force%tic()
    do i = 1, ps% Nmax
       s1 = ps% species(i)
       if (s1 <= 0) continue
       x = ps% pos(:, i)
       do j = i + 1, ps% Nmax
          s2 = ps% species(j)
          if (s2 <= 0) continue
          d = rel_pos(x, ps% pos(:, j), L)
          r_sq = sum(d**2)
          if ( r_sq <= lj_params% cut_sq(s1, s2) ) then
             f = lj_force(d, r_sq, lj_params% epsilon(s1, s2), lj_params% sigma(s1, s2))
             ps% force(:, i) = ps% force(:, i) + f
             ps% force(:, j) = ps% force(:, j) - f
             e = e + lj_energy(r_sq, lj_params% epsilon(s1, s2), lj_params% sigma(s1, s2))
          end if
       end do
    end do
    call ps%time_self_force%tac()

  end function compute_force_n2

  !! Prepare a stencil of neighboring cells
  !!
  !! Precomputing the stencil of cells for the neighbor list is described in P. J. in 't Veld,
  !! S. J. Plimpton and G. S. Grest. Comp. Phys. Commun. 179, 320-329 (2008) but for regular
  !! Molecular Dynamics. The present version is suited for MPCD solvent.
  subroutine make_stencil(this, cells, cut)
    class(neighbor_list_t), intent(inout) :: this
    type(cell_system_t), intent(in) :: cells
    double precision :: cut

    integer :: max_i, count
    integer :: i,j,k
    double precision :: a

    a = cells% a

    max_i = floor(cut / a) + 1

    count = 0
    do i = -max_i, max_i
       do j = -max_i, max_i
          do k = -max_i, max_i
             if (closest( [i,j,k], [0,0,0] ) < cut) count = count + 1
          end do
       end do
    end do

    if (allocated(this% stencil)) deallocate(this% stencil)
    allocate(this% stencil(3, count))

    count = 0
    do i = -max_i, max_i
       do j = -max_i, max_i
          do k = -max_i, max_i
             if (closest( [i,j,k], [0,0,0] ) < cut) then
                count = count + 1
                this% stencil(:,count) = [i,j,k]
             end if
          end do
       end do
    end do

  end subroutine make_stencil

  pure function closest(x1, x2) result(c)
    integer, dimension(3), intent(in) :: x1, x2
    double precision :: c

    integer :: i, dist_sq

    dist_sq = 0
    do i = 1, 3
       dist_sq = dist_sq + max( 0, abs(x1(i) - x2(i)) - 1)**2
    end do
    c = sqrt(dble(dist_sq))

  end function closest

end module neighbor_list
