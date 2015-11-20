module md
  use particle_system
  implicit none

  private

  public :: md_pos
  public :: md_vel
  public :: rattle_dimer_pos
  public :: rattle_dimer_vel

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

    integer :: k

    !$omp parallel do
    do k = 1, particles% Nmax
       particles% vel(:,k) = particles% vel(:,k) + &
            dt * ( particles% force(:,k) + particles% force_old(:,k) ) / 2
    end do

  end subroutine md_vel

  subroutine rattle_dimer_pos(p, d, dt)
    type(particle_system_t), intent(inout) :: p
    double precision, intent(in) :: d
    double precision, intent(in) :: dt

    double precision :: g
    double precision :: s(3) ! direction vector
    double precision :: r(3) ! old direction vector
    double precision :: mass1, mass2, inv_mass

    r = p% pos_rattle(:,1) - p% pos_rattle(:,2)
    s = p% pos(:,1) - p% pos(:,2)
    mass1 = p%mass(p%species(1))
    mass2 = p%mass(p%species(2))
    inv_mass = 1/mass1 + 1/mass2

    g = (dot_product(s, s) - d**2) / (2 * inv_mass * dot_product(r, s))

    p% pos_old(:,1) = p% pos_old(:,1) - g*r/mass1
    p% pos_old(:,2) = p% pos_old(:,2) + g*r/mass2

  end subroutine rattle_dimer_pos

  subroutine rattle_dimer_vel(p, d, dt)
    type(particle_system_t), intent(inout) :: p
    double precision, intent(in) :: d
    double precision, intent(in) :: dt

    double precision :: k !second correction factor
    double precision :: s(3) !direction vector
    double precision :: mass1, mass2, inv_mass

    mass1 = p%mass(p%species(1))
    mass2 = p%mass(p%species(2))
    inv_mass = 1/mass1 + 1/mass2

    s = p% pos(:,1) - p% pos(:,2)

    k = dot_product(p%vel(:,1)-p%vel(:,2), s) / (d**2*inv_mass)

    p% vel(:,1) = p% vel(:,1) - k*s/mass1
    p% vel(:,2) = p% vel(:,2) + k*s/mass2

  end subroutine rattle_dimer_vel

end module md
