! This file is part of RMPCDMD
! Copyright (c) 2015-2017 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Routines to perform Molecular Dynamics integration
!!
!! This module contains routines to perform velocity Verlet integration, Rattle for a dimer
!! and for other rigid bodies, and quaternion-based rotational velocity Verlet.
!!
!! Routines to compute colloid-wall (with Lennard-Jones 9-3) and elastic network
!! interactions.

module md
  use particle_system
  use cell_system
  use interaction
  use common
  use quaternion
  implicit none

  private

  public :: md_pos
  public :: apply_pbc
  public :: md_vel
  public :: rattle_dimer_pos
  public :: rattle_dimer_vel
  public :: rattle_body_pos, rattle_body_vel
  public :: lj93_zwall
  public :: elastic_network
  public :: rigid_body_t
  public :: cell_md_pos, cell_md_vel
  public :: cell_md_pos_ywall
  public :: cell_md_pos_zwall

  !> Container for rigid body colloid
  !!
  !! Stores the center-of-mass coordinates, the orientation, etc.
  type rigid_body_t
     integer :: i_start
     integer :: i_stop
     double precision :: mass
     double precision :: pos(3)
     double precision, allocatable :: pos_body(:,:)
     double precision :: vel(3)
     double precision :: force(3)
     double precision :: q(4)
     double precision :: L(3)
     double precision :: torque(3)
     double precision :: I_body(3)
     double precision :: omega_body(3)
   contains
     procedure :: init => rigid_body_init
     procedure :: compute_force_torque => rigid_body_compute_force_torque
     procedure :: vv1 => rigid_body_vv1
     procedure :: vv2 => rigid_body_vv2
  end type rigid_body_t

contains

  subroutine md_pos(particles, dt)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt

    double precision :: dt_sq
    integer :: k

    dt_sq = dt**2/2

    call particles%time_md_pos%tic()
    !$omp parallel do
    do k = 1, particles% Nmax
       particles% pos(:,k) = particles% pos(:,k) + dt * particles% vel(:,k) + dt_sq * particles% force(:,k)
    end do
    call particles%time_md_pos%tac()

  end subroutine md_pos

  subroutine cell_md_pos(cells, particles, dt, md_flag)
    type(cell_system_t), intent(inout) :: cells
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt
    logical, intent(in) :: md_flag

    double precision :: dt_sq
    integer :: i, j

    dt_sq = dt**2/2

    call particles%time_md_pos%tic()
    !$omp parallel do private(i, j)
    do i = 1, cells%N
       if (cells%is_md(i) .eqv. md_flag) then
          do j = cells%cell_start(i), cells%cell_start(i) + cells%cell_count(i) - 1
             particles% pos(:,j) = particles% pos(:,j) + dt * particles% vel(:,j) + dt_sq * particles% force(:,j)
          end do
       end if
    end do

    call particles%time_md_pos%tac()

  end subroutine cell_md_pos

  subroutine cell_md_pos_ywall(cells, particles, dt, md_flag)
    type(cell_system_t), intent(inout) :: cells
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt
    logical, intent(in) :: md_flag

    double precision :: dt_sq
    double precision :: Ly, x(3), v(3), overtime
    integer :: i, j, bc
    logical :: collide

    dt_sq = dt**2/2
    Ly = cells%edges(2)
    bc = cells%bc(2)

    call particles%time_md_pos%tic()
    !$omp parallel do private(i, j, x, v, collide, overtime)
    do i = 1, cells%N
       if (cells%is_md(i) .eqv. md_flag) then
          do j = cells%cell_start(i), cells%cell_start(i) + cells%cell_count(i) - 1
             v = particles%vel(:,j)
             x = particles%pos(:,j) + dt * v + dt_sq * particles%force(:,j)
             collide = .false.
             if (x(2)<0.d0) then
                collide = .true.
                overtime = x(2)/v(2)
             else if (x(2)>Ly) then
                collide = .true.
                overtime = (x(2)-Ly)/v(2)
             end if
             if (collide) then
                if (bc==BOUNCE_BACK_BC) then
                   v = -v
                   x = x+2*overtime*v
                else if (bc==SPECULAR_BC) then
                   v(2) = -v(2)
                   x(2) = x(2)+2*overtime*v(2)
                end if
                particles%vel(:,j) = v
             end if
             particles% pos(:,j) = x
          end do
       end if
    end do

    call particles%time_md_pos%tac()

  end subroutine cell_md_pos_ywall

  subroutine cell_md_pos_zwall(cells, particles, dt, md_flag)
    type(cell_system_t), intent(inout) :: cells
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt
    logical, intent(in) :: md_flag

    double precision :: dt_sq
    double precision :: Lz, x(3), v(3), overtime
    integer :: i, j, bc
    logical :: collide

    dt_sq = dt**2/2
    Lz = cells%edges(3)
    bc = cells%bc(3)

    call particles%time_md_pos%tic()
    !$omp parallel do private(i, j, x, v, collide, overtime)
    do i = 1, cells%N
       if (cells%is_md(i) .eqv. md_flag) then
          do j = cells%cell_start(i), cells%cell_start(i) + cells%cell_count(i) - 1
             v = particles%vel(:,j)
             x = particles%pos(:,j) + dt * v + dt_sq * particles%force(:,j)
             collide = .false.
             if (x(3)<0.d0) then
                collide = .true.
                overtime = x(3)/v(3)
             else if (x(3)>Lz) then
                collide = .true.
                overtime = (x(3)-Lz)/v(3)
             end if
             if (collide) then
                if (bc==BOUNCE_BACK_BC) then
                   v = -v
                   x = x+2*overtime*v
                else if (bc==SPECULAR_BC) then
                   v(3) = -v(3)
                   x(3) = x(3)+2*overtime*v(3)
                end if
                particles%vel(:,j) = v
             end if
             particles% pos(:,j) = x
          end do
       end if
    end do

    call particles%time_md_pos%tac()

  end subroutine cell_md_pos_zwall

  subroutine apply_pbc(particles, edges)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: edges(3)

    integer :: k, jump(3)

    call particles%time_apply_pbc%tic()
    !$omp parallel do private(jump)
    do k = 1, particles% Nmax
       jump = floor(particles% pos(:,k) / edges)
       particles% image(:,k) = particles% image(:,k) + jump
       particles% pos(:,k) = particles% pos(:,k) - jump*edges
    end do
    call particles%time_apply_pbc%tac()

  end subroutine apply_pbc

  subroutine md_vel(particles, dt)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt

    integer :: k

    call particles%time_md_vel%tic()
    !$omp parallel do
    do k = 1, particles% Nmax
       if (.not. btest(particles%flags(k), WALL_BIT)) then
          particles% vel(:,k) = particles% vel(:,k) + &
               dt * ( particles% force(:,k) + particles% force_old(:,k) ) / 2
       else
          particles%flags(k) = ibclr(particles%flags(k), WALL_BIT)
       end if
    end do
    call particles%time_md_vel%tac()

  end subroutine md_vel

  subroutine cell_md_vel(cells, particles, dt, md_flag)
    type(cell_system_t), intent(in) :: cells
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: dt
    logical, intent(in) :: md_flag

    integer :: i, j

    call particles%time_md_vel%tic()
    !$omp parallel do private(i, j)
    do i = 1, cells%N
       if (cells%is_md(i) .eqv. md_flag) then
          do j = cells%cell_start(i), cells%cell_start(i) + cells%cell_count(i) - 1
             if (.not. btest(particles%flags(j), WALL_BIT)) then
                particles% vel(:,j) = particles% vel(:,j) + &
                     dt * ( particles% force(:,j) + particles% force_old(:,j) ) / 2
             else
                particles%flags(j) = ibclr(particles%flags(j), WALL_BIT)
             end if
          end do
       end if
    end do
    call particles%time_md_vel%tac()

  end subroutine cell_md_vel

  subroutine rattle_dimer_pos(p, d, dt,edges)
    type(particle_system_t), intent(inout) :: p
    double precision, intent(in) :: d
    double precision, intent(in) :: dt
    double precision, intent(in) :: edges(3)

    double precision :: g
    double precision :: s(3) ! direction vector
    double precision :: r(3) ! old direction vector
    double precision :: rsq, ssq, rs
    double precision :: mass1, mass2, inv_mass

    r = rel_pos(p% pos_rattle(:,1),p% pos_rattle(:,2), edges)
    s = rel_pos(p% pos(:,1),p% pos(:,2), edges)
    mass1 = p%mass(p%species(1))
    mass2 = p%mass(p%species(2))
    inv_mass = 1/mass1 + 1/mass2

    rsq = dot_product(r,r)
    ssq = dot_product(s,s)
    rs = dot_product(r,s)

    g = rs - sqrt(rs**2 - rsq*(ssq-d**2))
    g = g / (dt * inv_mass * rsq)

    p% pos(:,1) = p% pos(:,1) - g*dt*r/mass1
    p% pos(:,2) = p% pos(:,2) + g*dt*r/mass2

    p% vel(:,1) = p% vel(:,1) - g*r/mass1
    p% vel(:,2) = p% vel(:,2) + g*r/mass2

  end subroutine rattle_dimer_pos

  subroutine rattle_dimer_vel(p, d, dt,edges)
    type(particle_system_t), intent(inout) :: p
    double precision, intent(in) :: d
    double precision, intent(in) :: dt
    double precision, intent(in) :: edges(3)

    double precision :: k !second correction factor
    double precision :: s(3) !direction vector
    double precision :: mass1, mass2, inv_mass

    mass1 = p%mass(p%species(1))
    mass2 = p%mass(p%species(2))
    inv_mass = 1/mass1 + 1/mass2

    s = rel_pos(p% pos(:,1), p% pos(:,2), edges)

    k = dot_product(p%vel(:,1)-p%vel(:,2), s) / (d**2*inv_mass)

    p% vel(:,1) = p% vel(:,1) - k*s/mass1
    p% vel(:,2) = p% vel(:,2) + k*s/mass2

  end subroutine rattle_dimer_vel

  subroutine rattle_body_pos(p, links, distances, dt, edges, precision)
    type(particle_system_t), intent(inout) :: p
    integer, intent(in) :: links(:,:)
    double precision, intent(in) :: distances(:)
    double precision, intent(in) :: dt, edges(3), precision

    double precision :: g, d, error
    double precision :: s(3) ! direction vector
    double precision :: r(3) ! old direction vector
    double precision :: rsq, ssq, rs
    double precision :: i_mass1, i_mass2, inv_mass
    double precision, allocatable :: i_mass(:)
    integer, allocatable :: im(:,:), im_r(:,:)

    integer :: rattle_i, rattle_max, i_link, n_link
    integer :: i1, i2

    n_link = size(links, dim=2)
    rattle_max = 1000

    allocate(i_mass(size(p%mass)))
    i_mass = 1/p%mass

    allocate(im(3,p%Nmax))
    allocate(im_r(3,p%Nmax))

    s = p%pos(:,1)
    do i1 = 2, p%Nmax
       im(:,i1) = floor(((p%pos(:,i1)-s)/edges) + 0.5d0)
       p%pos(:,i1) = p%pos(:,i1) - im(:,i1)*edges
    end do

    s = p%pos_rattle(:,1)
    do i1 = 2, p%Nmax
       im_r(:,i1) = floor(((p%pos_rattle(:,i1)-s)/edges) + 0.5d0)
       p%pos_rattle(:,i1) = p%pos_rattle(:,i1) - im_r(:,i1)*edges
    end do

    call p%time_rattle_pos%tic()
    rattle_loop: do rattle_i = 1, rattle_max
       error = 0
       do i_link = 1, n_link
          i1 = links(1,i_link)
          i2 = links(2,i_link)
          d = distances(i_link)

          r = p%pos_rattle(:,i1)-p%pos_rattle(:,i2)
          s = p%pos(:,i1)-p%pos(:,i2)
          i_mass1 = i_mass(p%species(i1))
          i_mass2 = i_mass(p%species(i2))
          inv_mass = i_mass1 + i_mass2

          rsq = dot_product(r,r)
          ssq = dot_product(s,s)
          rs = dot_product(r,s)

          g = rs - sqrt(rs**2 - rsq*(ssq-d**2))
          g = g / (dt * inv_mass * rsq)

          p% pos(:,i1) = p% pos(:,i1) - g*dt*r*i_mass1
          p% pos(:,i2) = p% pos(:,i2) + g*dt*r*i_mass2

          p% vel(:,i1) = p% vel(:,i1) - g*r*i_mass1
          p% vel(:,i2) = p% vel(:,i2) + g*r*i_mass2
       end do

       do i_link = 1, n_link
          i1 = links(1,i_link)
          i2 = links(2,i_link)
          d = distances(i_link)
          r = p%pos(:,i1)-p%pos(:,i2)
          g = abs(sqrt(dot_product(r,r)) - d)
          if (g > error) error = g
       end do
       if (error < precision) exit rattle_loop

    end do rattle_loop
    call p%time_rattle_pos%tac()
    deallocate(i_mass)

    do i1 = 2, p%Nmax
       p%pos(:,i1) = p%pos(:,i1) + im(:,i1)*edges
    end do
    deallocate(im, im_r)

    if (rattle_i==rattle_max) write(*,*) 'rattle_max reached in rattle_body_pos'

  end subroutine rattle_body_pos

  subroutine rattle_body_vel(p, links, distances, dt, edges, precision)
    type(particle_system_t), intent(inout) :: p
    integer, intent(in) :: links(:,:)
    double precision, intent(in) :: distances(:)
    double precision, intent(in) :: dt, edges(3), precision

    double precision :: g, d, error
    double precision :: s(3), k
    double precision :: i_mass1, i_mass2, inv_mass
    double precision, allocatable :: i_mass(:)
    integer, allocatable :: im(:,:)

    integer :: rattle_i, rattle_max, i_link, n_link
    integer :: i1, i2

    n_link = size(links, dim=2)
    rattle_max = 1000
    allocate(i_mass(size(p%mass)))
    i_mass = 1/p%mass

    allocate(im(3,p%Nmax))

    s = p%pos(:,1)
    do i1 = 2, p%Nmax
       im(:,i1) = floor(((p%pos(:,i1)-s)/edges) + 0.5d0)
       p%pos(:,i1) = p%pos(:,i1) - im(:,i1)*edges
    end do

    call p%time_rattle_vel%tic()
    rattle_loop: do rattle_i = 1, rattle_max
       error = 0
       do i_link = 1, n_link
          i1 = links(1,i_link)
          i2 = links(2,i_link)
          d = distances(i_link)
          i_mass1 = i_mass(p%species(i1))
          i_mass2 = i_mass(p%species(i2))
          inv_mass = i_mass1 + i_mass2

          s = p%pos(:,i1)-p%pos(:,i2)

          k = dot_product(p%vel(:,i1)-p%vel(:,i2), s) / (d**2*inv_mass)

          p% vel(:,i1) = p% vel(:,i1) - k*s*i_mass1
          p% vel(:,i2) = p% vel(:,i2) + k*s*i_mass2
       end do

       do i_link = 1, n_link
          i1 = links(1,i_link)
          i2 = links(2,i_link)
          d = distances(i_link)
          i_mass1 = i_mass(p%species(i1))
          i_mass2 = i_mass(p%species(i2))
          inv_mass = i_mass1 + i_mass2
          s = p%pos(:,i1)-p% pos(:,i2)
          k = abs(dot_product(p%vel(:,i1)-p%vel(:,i2), s) / (d**2*inv_mass))
          if (k>error) error = k
       end do
       if (error < precision) exit rattle_loop
    end do rattle_loop
    call p%time_rattle_vel%tac()

    do i1 = 2, p%Nmax
       p%pos(:,i1) = p%pos(:,i1) + im(:,i1)*edges
    end do

    deallocate(i_mass, im)

    if (rattle_i==rattle_max) write(*,*) 'rattle_max reached in rattle_body_vel'

  end subroutine rattle_body_vel

  function lj93_zwall(particles, edges, lj_params) result(e)
    type(particle_system_t), intent(inout) :: particles
    double precision, intent(in) :: edges(3)
    type(lj_params_t), intent(in) :: lj_params
    double precision :: e

    integer :: i, s, dim
    double precision :: z, z_sq, f, dir
    
    e = 0
    do i = 1, particles%Nmax
       s = particles%species(i)
       if (s<=0) cycle
       do dim = 1, 3
          if (lj_params%sigma(dim, s)<0) cycle
          z = particles%pos(dim,i)
          if (z > edges(dim)/2) then
             z = edges(dim) - z - lj_params%shift(dim)
             dir = -1
          else
             dir = 1
             z = z - lj_params%shift(dim)
          end if
          z_sq = z**2
          if (z_sq <= lj_params% cut_sq(dim,s)) then
             f = lj_force_9_3(z, z_sq, lj_params%epsilon(dim,s), lj_params%sigma(dim,s))
             particles%force(dim,i) = particles%force(dim,i) + dir*f
             e = e + lj_energy_9_3(z_sq, lj_params%epsilon(dim,s), lj_params%sigma(dim,s))
          end if
       end do
    end do

  end function lj93_zwall

  function elastic_network(p, links, distances, k, edges) result(e)
    type(particle_system_t), intent(inout) :: p
    integer, intent(in) :: links(:,:)
    double precision, intent(in) :: distances(:)
    double precision, intent(in) :: k
    double precision, intent(in) :: edges(3)
    double precision :: e

    integer :: i, i1, i2, n_links
    double precision :: x(3), f(3), r, r0

    n_links = size(distances)

    call p%time_elastic%tic()
    e = 0
    do i = 1, n_links
       i1 = links(1,i)
       i2 = links(2,i)
       r0 = distances(i)
       x = rel_pos(p%pos(:,i1), p%pos(:,i2), edges)
       r = norm2(x)
       e = e + k*(r-r0)**2/2
       f = - k * (r - r0) * x / r
       p%force(:,i1) = p%force(:,i1) + f
       p%force(:,i2) = p%force(:,i2) - f
    end do
    call p%time_elastic%tac()

  end function elastic_network

  !! Compute body frame positions and inertia tensor
  subroutine rigid_body_init(this, ps, i_start, i_stop, edges)
    class(rigid_body_t), intent(out) :: this
    type(particle_system_t), intent(inout) :: ps
    integer, intent(in) :: i_start
    integer, intent(in) :: i_stop
    double precision, intent(in) :: edges(3)

    double precision :: mass
    double precision :: pos(3)
    integer :: i

    this%i_start = i_start
    this%i_stop = i_stop

    this%mass = 0
    this%pos = 0
    this%vel = 0
    do i = i_start, i_stop
       mass = ps%mass(ps%species(i))
       this%mass = this%mass + mass
       this%pos = this%pos + mass*(ps%pos(:,i)+ps%image(:,i)*edges)
       this%vel = this%vel + mass*ps%vel(:,i)
    end do
    this%pos = this%pos/this%mass
    this%vel = this%vel/this%mass

    allocate(this%pos_body(3,i_stop-i_start+1))

    this%I_body = 0
    do i = i_start, i_stop
       mass = ps%mass(ps%species(i))
       pos = ps%pos(:,i)+ps%image(:,i)*edges - this%pos
       this%pos_body(:,i) = pos
       this%I_body = this%I_body + mass*(sum(pos**2)-pos**2)
    end do

    this%force = 0
    this%torque = 0
    this%q = qnew(s=1.d0)
    this%L = 0

  end subroutine rigid_body_init

  !! Compute torque and force on the rigid body in the lab frame
  subroutine rigid_body_compute_force_torque(this, ps)
    class(rigid_body_t), intent(inout) :: this
    type(particle_system_t), intent(in) :: ps

    double precision :: pos(3), f(3)
    integer :: i

    this%force = 0
    this%torque = 0
    do i = this%i_start, this%i_stop
       f = ps%force(:,i)
       this%force = this%force + f
       pos = qrot(this%q, this%pos_body(:,i))
       this%torque = this%torque + cross(pos, f)
    end do

  end subroutine rigid_body_compute_force_torque

  !> Perform first velocity-Verlet step for quaternion dynamics
  !!
  !! The algorithm is described in \cite rozmanov_rotational_2010
  !! \manual{algorithms/quaternions}
  subroutine rigid_body_vv1(this, ps, dt, treshold, edges)
    class(rigid_body_t), intent(inout) :: this
    type(particle_system_t), intent(inout) :: ps
    double precision, intent(in) :: dt
    double precision, intent(in) :: treshold
    double precision, intent(in) :: edges(3)

    double precision :: L(3), L_body(3), L_body_dot(3), omega_body(3), q(4), q_old(4), q_dot(4), torque_body(3)
    double precision :: I_inv(3)
    double precision :: old_pos(3), delta_pos(3)

    integer :: i, j

    I_inv = 1/this%I_body

    ! timestep 0
    L = this%L
    q = this%q

    L_body = qrot(qconj(q), L)
    torque_body = qrot(qconj(q), this%torque)
    omega_body = L_body*I_inv

    L_body_dot = torque_body - cross(omega_body, L_body)

    ! timestep 1/2
    L_body = L_body + dt*L_body_dot/2
    omega_body = L_body*I_inv
    q_dot = qmul(q, qnew(v=omega_body))/2
    q = q + q_dot*dt/2
    q = qnormalize(q)
    L =  this%L + this%torque*dt/2

    ! Iterative procedure

    q_old = q
    iterative_q: do while (.true.)
       L_body = qrot(qconj(q), L)
       omega_body = I_inv*L_body
       q_dot = qmul(q, qnew(v=omega_body))/2
       q = this%q + q_dot*dt/2
       q = qnormalize(q)
       if (qnorm(q-q_old)<treshold) exit iterative_q
       q_old = q
    end do iterative_q

    this%q = qnormalize(this%q + dt*q_dot)

    this%pos = this%pos + this%vel*dt + this%force*dt**2/(2*this%mass)

    do i = this%i_start, this%i_stop
       old_pos = ps%pos(:,i)
       ps%pos(:,i) = this%pos + qrot(this%q, this%pos_body(:,i))
       delta_pos = ps%pos(:,i) - old_pos + edges/2
       ps%pos(:,i) = ps%pos(:,i) - floor(delta_pos / edges)*edges
    end do

    this%L = this%L + this%torque*dt/2
    this%vel = this%vel + this%force*dt/(2*this%mass)

  end subroutine rigid_body_vv1

  !> Perform second velocity-Verlet step for quaternion dynamics
  !!
  !! The algorithm is described in \cite rozmanov_rotational_2010
  !! \manual{algorithms/quaternions}
  subroutine rigid_body_vv2(this, ps, dt)
    class(rigid_body_t), intent(inout) :: this
    type(particle_system_t), intent(inout) :: ps
    double precision, intent(in) :: dt

    double precision :: L_body(3), omega(3)
    double precision :: pos(3)
    integer :: i

    this%L = this%L + this%torque*dt/2
    L_body = qrot(qconj(this%q), this%L)
    this%omega_body = L_body/this%I_body
    omega = qrot(this%q, this%omega_body)

    this%vel = this%vel + this%force*dt/(2*this%mass)

    do i = this%i_start, this%i_stop
       pos = qrot(this%q, this%pos_body(:,i))
       ps%vel(:,i) = this%vel + cross(omega, pos)
    end do

  end subroutine rigid_body_vv2

end module md
