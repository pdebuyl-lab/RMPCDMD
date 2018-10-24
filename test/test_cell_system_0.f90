! This file is part of RMPCDMD
! Copyright (c) 2015-2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program try_all
  use cell_system
  use hilbert
  use tester
  implicit none

  type(cell_system_t) :: solvent_cells
  type(tester_t) :: test

  integer, parameter :: N = 1000000
  double precision, target, allocatable :: pos1(:, :), pos2(:, :)
  double precision, pointer :: pos(:,:), pos_old(:,:), pos_pointer(:,:)

  integer :: i, L(3), seed_size, clock
  integer, allocatable :: seed(:)
  integer :: h, c1(3), c2(3)

  allocate(pos1(3, N))
  allocate(pos2(3, N))

  call test% init()

  call random_seed(size = seed_size)
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37 * [ (i - 1, i = 1, seed_size) ]
  call random_seed(put = seed)
  deallocate(seed)

  pos => pos1
  pos_old => pos2

  L = [120, 30, 30]

  call random_number(pos)
  do i = 1, 3
     pos(i, :) = pos(i, :)*L(i)
  end do

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%count_particles(pos)

  call solvent_cells%sort_particles(pos, pos_old)

  pos_pointer => pos
  pos => pos_old
  pos_old => pos_pointer

  call solvent_cells%count_particles(pos)

  call test% assert_equal(solvent_cells%cell_start(1), 1)
  call test% assert_equal(solvent_cells%cell_start(solvent_cells%N), N+1)
  call test% assert_equal(sum(solvent_cells%cell_count), N)

  h = 1
  do i = 1, N
     if (i .ge. solvent_cells% cell_start(h) + solvent_cells% cell_count(h)) h = h + 1
     do while (solvent_cells% cell_count(h) .eq. 0)
        h = h+1
     end do
     c1 = solvent_cells%cartesian_indices(pos(:,i))
     c2 = compact_h_to_p(h-1, solvent_cells% M)
     call test% assert_equal(c1, c2)
  end do

  call test% print()

  deallocate(pos1)
  deallocate(pos2)

end program try_all
