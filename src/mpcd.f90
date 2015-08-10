module mpcd
  use particle_system
  use cell_system
  implicit none

  private

  public :: compute_temperature

contains

  function compute_temperature(particles, cells) result(te)
    class(particle_system_t), intent(in) :: particles
    class(cell_system_t), intent(in) :: cells
    double precision :: te

    integer :: i, start, n
    integer :: cell_idx, count
    double precision :: local_v(3), local_kin

    cell_idx = 1
    count = 0
    te = 0
    do cell_idx = 1, cells% N
       if (cells% cell_count(cell_idx) <= 1) cycle

       start = cells% cell_start(cell_idx)
       n = cells% cell_count(cell_idx)
       count = count + 1

       local_v = 0
       local_kin = 0
       do i = start, start + n - 1
          local_v = local_v + particles% vel(:, i)
       end do
       local_v = local_v / n
       do i = start, start + n - 1
          local_kin = local_kin + sum((particles% vel(:, i)-local_v)**2)
       end do
       te = te + local_kin / dble(3*(n-1))

    end do

    te = te / count

  end function compute_temperature

end module mpcd
