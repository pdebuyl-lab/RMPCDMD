program try_all
  use cell_system
  use particle_system
  use hilbert
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent

  integer, parameter :: N = 1000

  integer :: i, idx, j, L(3), p(3)

  L = [8, 3, 4]

  call solvent% init(N)

  call random_number(solvent% vel(:, i))
  solvent% force = 0
  solvent% species = 1
  do i = 1, N
     call random_number(solvent% pos(:, i))
     do j=1, 3
        solvent% pos(j, i) = solvent% pos(j, i) * L(j)
     end do
  end do

  call solvent_cells%init(L, 1.d0)

  print *, solvent_cells%M

  call solvent_cells%count_particles(solvent% pos)

  print *, sum(solvent_cells%cell_count)
  print *, solvent_cells%cell_count
  print *, solvent_cells%cell_start(1), solvent_cells%cell_start(solvent_cells%N)

  call solvent_cells%sort_particles(solvent% pos, solvent% pos_old)
  solvent% pos_pointer => solvent% pos
  solvent% pos => solvent% pos_old
  solvent% pos_old => solvent% pos_pointer

  open(12, file='sorted_pos')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  call random_number(solvent% pos(:, :))
  do j = 1, 3
     solvent% pos(j, :) = solvent% pos(j, :) * L(j)
  end do

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
