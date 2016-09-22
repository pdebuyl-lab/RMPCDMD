program test_neighbor_list
  use common
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use threefry_module
  use iso_c_binding
  use tester
  implicit none

  type(tester_t) :: test

  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(cell_system_t) :: solvent_cells
  type(neighbor_list_t) :: neigh

  integer, parameter :: rho = 20
  integer :: N_solvent, N_colloids

  integer :: L(3)
  integer :: i, j, n_neigh
  integer :: neighbor_list_size

  double precision :: xij(3), dist_sq, cutoff

  type(threefry_rng_t) :: state(1)
  integer :: seed_size, clock
  integer, allocatable :: seed(:)

  call test% init()

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call threefry_rng_init(state, int(clock, c_int64_t))

  L = [20, 10, 35]
  N_solvent = L(1)*L(2)*L(3)*rho
  N_colloids = 20

  call solvent% init(N_solvent)
  call colloids% init(N_colloids)

  call solvent_cells%init(L, 1.d0)

  solvent% species = 1
  colloids% species = 1
  call solvent% random_placement(L*1.d0, state=state(1))
  call colloids% random_placement(L*1.d0, state=state(1))

  call solvent% sort(solvent_cells)

  cutoff = 3.6
  neighbor_list_size = rho * 6 * int(cutoff**3)

  call neigh% init(colloids% Nmax, neighbor_list_size)

  call neigh% make_stencil(solvent_cells, cutoff)

  call neigh% update_list(colloids, solvent, cutoff*0.8, solvent_cells)

  do i = 1, colloids% Nmax
     n_neigh = neigh% n(i)
     do j = 1, solvent% Nmax
        xij = rel_pos(colloids% pos(:, i), solvent% pos(:, j), solvent_cells% edges)
        dist_sq = sum(xij**2)
        if (dist_sq < cutoff**2) then
           call test% assert_equal( is_in_list(j, neigh% list(1:n_neigh,i)), .true. )
        end if
     end do
  end do

  call test% print()

contains

  function is_in_list(idx, list) result(l)
    integer, intent(in) :: idx
    integer, intent(in) :: list(:)

    logical :: l

    l = ( minval( abs(list - idx) ) == 0 )

  end function is_in_list

end program test_neighbor_list
