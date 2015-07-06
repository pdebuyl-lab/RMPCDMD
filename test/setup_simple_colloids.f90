program setup_simple_colloids
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use hdf5
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh

  integer, parameter :: N = 3000
  integer :: error

  integer :: i, L(3), seed_size, clock
  integer, allocatable :: seed(:)

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  call h5open_f(error)

  L = [12, 12, 12]

  call solvent% init(N)

  call colloids% init_from_file('input_data.h5', 'colloids')

  call random_number(solvent% vel(:, :))
  solvent% vel(:, :) = solvent% vel(:, :) - 0.5d0
  solvent% force = 0
  solvent% species = 1
  call solvent% random_placement(L*1.d0, colloids% pos, 1.d0)

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(solvent% pos)

  print *, sum(solvent_cells%cell_count)
  print *, solvent_cells%cell_count
  print *, solvent_cells%cell_start(1), solvent_cells%cell_start(solvent_cells%N)

  call sort

  call neigh% init(colloids% Nmax, 300)

  call neigh% update_list(colloids, solvent, 1.4d0, solvent_cells)

  open(12, file='neigh_list')
  do i = 1, colloids% Nmax
     write(12, *) neigh% n(i)
     write(12, *) neigh% list(:, i)
  end do

  open(12, file='sorted_pos')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  open(12, file='colloids_pos')
  do i=1, colloids% Nmax
     write(12, *) colloids% pos(:,i)
  end do
  close(12)

  do i = 1, 1000

     call solvent% random_placement(L*1.d0)
     call sort
     call neigh% update_list(colloids, solvent, 1.4d0, solvent_cells)

  end do

  open(12, file='sorted_pos_2')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  call h5close_f(error)

contains

  subroutine sort
    integer :: i, idx, start, p(3)

    call solvent_cells%count_particles(solvent% pos)

    do i=1, N
       p = floor( (solvent% pos(:, i) - solvent_cells% origin ) / solvent_cells% a )
       idx = compact_p_to_h(p, solvent_cells% M) + 1
       start = solvent_cells% cell_start(idx)
       solvent% pos_old(:, start) = solvent% pos(:, i)
       solvent% vel_old(:, start) = solvent% vel(:, i)
       solvent% force_old(:, start) = solvent% force(:, i)
       solvent% id_old(start) = solvent% id(i)
       solvent% species_old(start) = solvent% species(i)
       solvent_cells% cell_start(idx) = start + 1
    end do

    solvent% pos_pointer => solvent% pos
    solvent% pos => solvent% pos_old
    solvent% pos_old => solvent% pos_pointer

    solvent% vel_pointer => solvent% vel
    solvent% vel => solvent% vel_old
    solvent% vel_old => solvent% vel_pointer

    solvent% force_pointer => solvent% force
    solvent% force => solvent% force_old
    solvent% force_old => solvent% force_pointer

    solvent% id_pointer => solvent% id
    solvent% id => solvent% id_old
    solvent% id_old => solvent% id_pointer

    solvent% species_pointer => solvent% species
    solvent% species => solvent% species_old
    solvent% species_old => solvent% species_pointer

  end subroutine sort

end program setup_simple_colloids
