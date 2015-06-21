program try_all
  use cell_system
  use particle_system
  use hilbert
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids

  integer, parameter :: N = 1000
  integer, parameter :: N_colloids = 5

  integer :: i, idx, j, L(3), p(3)

  L = [8, 3, 4]

  call solvent% init(N)
  call colloids% init(N_colloids)

  call colloids% random_placement(L*1.d0)

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

  open(12, file='sorted_pos')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  open(12, file='colloids_pos')
  do i=1, N_colloids
     write(12, *) colloids% pos(:,i)
  end do
  close(12)

  call solvent% random_placement(L*1.d0)

  call sort

  open(12, file='sorted_pos_2')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

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

end program try_all
