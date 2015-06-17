program try_all
  use cell_system
  implicit none

  type(cell_system_t) :: solvent_cells

  integer :: L(3)

  call solvent_cells%init(L, 1.d0)

  call solvent_cells%del()

end program try_all
