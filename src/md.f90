module md
  use particle_system
  implicit none

  private

  public :: md_pos
  public :: md_vel

contains

  subroutine md_pos(particles, edges, dt)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: edges(3)
    double precision, intent(in) :: dt

    integer :: k, jump(3)
    double precision :: tmp_x(3)

    !$omp parallel do private(tmp_x, jump)
    do k = 1, particles% Nmax
       tmp_x = particles% pos(:,k) + dt * particles% vel(:,k) + dt**2 * particles% force(:,k) / 2

       jump = floor(particles% pos(:,k) / edges)
       particles% image(:,k) = particles% image(:,k) + jump

       particles% pos(:,k) = tmp_x - jump * edges
    end do

  end subroutine md_pos

  subroutine md_vel(particles, edges, dt)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: edges(3)
    double precision, intent(in) :: dt

    integer :: k, jump(3)
    double precision :: tmp_x(3)

    !$omp parallel do
    do k = 1, particles% Nmax
       particles% vel(:,k) = particles% vel(:,k) + dt * ( particles% force(:,k) + particles% force_old(:,k) ) / 2
    end do

  end subroutine md_vel

end module md
