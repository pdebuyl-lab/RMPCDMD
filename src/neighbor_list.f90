module neighbor_list
  use common
  use hilbert
  use cell_system
  use particle_system
  use interaction
  implicit none
  private

  public :: neighbor_list_t, compute_force

  type neighbor_list_t
     integer :: Npoints
     integer :: Nmax
     integer, allocatable :: list(:, :)
     integer, allocatable :: n(:)
   contains
     procedure :: init
     procedure :: update_list
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

  subroutine update_list(this, system1, system2, radius, cells)
    class(neighbor_list_t), intent(inout) :: this
    type(particle_system_t), intent(in) :: system1
    type(particle_system_t), intent(in) :: system2
    type(cell_system_t), intent(in) :: cells
    double precision, intent(in) :: radius
    
    integer :: cell(3), neigh_cell(3), M(3), actual_cell(3)
    integer :: neigh_idx, i, cell_i, cell_n, cell_start, list_idx
    double precision :: x(3), y(3), L(3)
    double precision :: rsq, radiussq

    L = cells% L * cells% a
    radiussq = radius**2
    M = cells% M

    do i = 1, system1% Nmax
       x = system1% pos(:, i)
       cell = floor( (x - cells% origin) / cells% a )
       neigh_cell = [ -2, -2, -2]
       neigh_idx = compact_p_to_h( cell + neigh_cell, M ) + 1
       cell_i = 1
       cell_n = cells% cell_count(neigh_idx)
       cell_start = cells% cell_start(neigh_idx)
       list_idx = 0

       if ( cell_n .eq. 0 ) continue

       do while ( .true. )

          y = system2% pos(:, cell_start+cell_i)
          rsq = sum(rel_pos(x, y, L)**2)
          if (rsq .lt. radiussq) then
             list_idx = list_idx + 1
             if (list_idx .gt. this% Nmax) then
                error stop 'maximum of neighbor list reached'
             end if
             this% list(list_idx, i) = cell_start+cell_i
          end if

          cell_i = cell_i + 1

          if ( cell_i .gt. cell_n ) then
             neigh_cell = next_cell(neigh_cell)
             actual_cell = modulo(cell + neigh_cell, cells% L)
             neigh_idx = compact_p_to_h( actual_cell, M ) + 1
             cell_i = 1
             cell_n = cells% cell_count(neigh_idx)
             cell_start = cells% cell_start(neigh_idx)
          end if

          if ( neigh_cell(3) .eq. 3) exit

       end do
       this% n(i) = list_idx

    end do

  end subroutine update_list

  subroutine compute_force(ps1, ps2, n_list, L, lj_params)
    type(particle_system_t), intent(inout) :: ps1
    type(particle_system_t), intent(inout) :: ps2
    type(neighbor_list_t), intent(in) :: n_list
    double precision, intent(in) :: L(3)
    type(lj_params_t), intent(in) :: lj_params

    integer :: i, j
    integer :: n_species1, n_species2
    double precision :: x(3), d(3), f(3), f1(3)
    integer :: s1, s2
    double precision :: r_sq

    n_species1 = ps1% n_species
    n_species2 = ps2% n_species

    do i = 1, ps1% Nmax
       if (ps1% species(i) <= 0) continue
       x = ps1% pos(:, i)
       s1 = ps1% species(i)
       f1 = 0
       do j = 1, n_list% n(i)
          s2 = ps2% species(j)
          if (s2 <= 0) continue
          d = rel_pos(x, ps2% pos(:, j), L)
          r_sq = sum(d**2)
          if ( r_sq <= lj_params% cut_sq(s1, s2) ) then
             f = lj_force(d, r_sq, lj_params% epsilon(s1, s2), lj_params% sigma(s1, s2))
             f1 = f1 + f
             ps2% force(:, j) = ps2% force(:, j) - f
          end if
       end do
       ps1% force(:, i) = ps1% force(:, i) + f1
    end do

  end subroutine compute_force

end module neighbor_list
