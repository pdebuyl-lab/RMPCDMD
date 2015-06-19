program try_all
  use cell_system
  use hilbert
  implicit none

  type(cell_system_t) :: solvent_cells

  integer :: i, idx, j, L(3), p(3)
  double precision :: r(3, 1000), r_new(3, 1000)

  L = [8, 3, 5]

  do i = 1, size(r, 2)
     call random_number(r(:, i))
     do j=1, 3
        r(j, i) = r(j, i) * L(j)
     end do
  end do

  call solvent_cells%init(L, 1.d0)

  print *, solvent_cells%M

  call solvent_cells%count_particles(r)

  print *, sum(solvent_cells%cell_count)
  print *, solvent_cells%cell_count
  print *, solvent_cells%cell_start(1), solvent_cells%cell_start(solvent_cells%N)

  call solvent_cells%sort_particles(r, r_new)

  do i=1, 1000
       p = floor( (r_new(:, i) - solvent_cells%origin ) / solvent_cells%a )
       idx = compact_p_to_h(p, solvent_cells%M) + 1
       write(*, *) idx
  end do

  call solvent_cells%del()

end program try_all
