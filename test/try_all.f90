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
  call solvent_cells%count_particles(solvent% pos)
  call solvent_cells%sort_particles(solvent% pos, solvent% pos_old)
  solvent% pos_pointer => solvent% pos
  solvent% pos => solvent% pos_old
  solvent% pos_old => solvent% pos_pointer

  open(12, file='sorted_pos_2')
  do i=1, N
     write(12, *) solvent% pos(:,i)
  end do
  close(12)

  call solvent_cells%del()
  call solvent% del()

end program try_all
