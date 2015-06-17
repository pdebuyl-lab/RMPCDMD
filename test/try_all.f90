program try_all
  use cell_system
  implicit none

  type(cell_system_t) :: solvent_cells

  integer :: i
  double precision :: r(3, 10)

  do i = 1, size(r, 2)
     r(:, i) = i * [0.01d0, 0.02d0, 0.05d0]
  end do

  call solvent_cells%init([8, 4, 4], 1.d0)

  call solvent_cells%count_particles(r)

  print *, solvent_cells%cell_count(1)

  call solvent_cells%del()

end program try_all
