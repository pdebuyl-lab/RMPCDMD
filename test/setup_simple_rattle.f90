! This file is part of RMPCDMD
! Copyright (c) 2015-2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

program setup_simple_colloids
  use common
  use cell_system
  use particle_system
  use hilbert
  use neighbor_list
  use hdf5
  use h5md_module
  use interaction
  use threefry_module
  use mpcd
  use md
  use iso_c_binding
  use omp_lib
  implicit none

  type(cell_system_t) :: solvent_cells
  type(particle_system_t) :: solvent
  type(particle_system_t) :: colloids
  type(neighbor_list_t) :: neigh
  type(lj_params_t) :: solvent_colloid_lj
  type(lj_params_t) :: colloid_lj

  integer, parameter :: rho = 10
  integer, parameter :: L(3) = [12, 12, 12]
  integer, parameter :: N = rho*L(1)*L(2)*L(3)

  integer :: error

  integer, parameter :: n_link=6
  integer :: links(2,n_link)
  double precision:: links_dist(n_link)

  double precision :: epsilon, sigma, sigma_cut
  double precision :: so_max, co_max

  double precision :: e1, e2
  double precision :: tau, dt
  integer :: N_MD_steps
  integer :: n_extra_sorting
  integer :: n_threads

  type(threefry_rng_t), allocatable :: state(:)

  double precision :: tmp_x(3)
  double precision :: skin

  integer :: i
  integer :: j, k

  n_threads = omp_get_max_threads()
  allocate(state(n_threads))
  call threefry_rng_init(state, -4466704919147519266_c_int64_t)

  call h5open_f(error)

  epsilon = 1
  sigma = 1.
  sigma_cut = sigma*2**(1.d0/6.d0)

  call solvent_colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ sigma ], [1, 1] ), reshape( [ sigma_cut ], [1, 1] ) )

  epsilon = 1
  sigma = 1
  sigma_cut = sigma*2**(1.d0/6.d0)

  call colloid_lj% init( reshape( [ epsilon ], [1, 1] ), &
       reshape( [ 2*sigma ], [1, 1] ), reshape( [ 2*sigma_cut ], [1, 1] ) )

  call solvent% init(N)

  call colloids%init(4)
  colloids% species = 1
  colloids% vel = 0

  colloids%mass = rho*sigma**3*4*3.14159/3

  colloids%pos(:,1) = [1.d0, 0.d0, -1/sqrt(2.d0)]
  colloids%pos(:,2) = [1.d0, 0.d0, 1/sqrt(2.d0)]
  colloids%pos(:,3) = [0.d0, 1.d0, 1/sqrt(2.d0)]
  colloids%pos(:,4) = [0.d0, -1.d0, 1/sqrt(2.d0)]
  colloids%pos = colloids%pos + 3

  links(:,1) = [1,2]
  links(:,2) = [1,3]
  links(:,3) = [1,4]
  links(:,4) = [2,3]
  links(:,5) = [2,4]
  links(:,6) = [3,4]
  do i = 1, n_link
     links_dist(i) = sqrt(sum( &
     rel_pos(colloids%pos(:,links(1,i)), colloids%pos(:,links(2,i)), solvent_cells%edges)**2 &
     ))
  end do

  call random_number(solvent% vel(:, :))
  solvent% vel = (solvent% vel - 0.5d0)*sqrt(6.d0*2)
  solvent% vel = solvent% vel - spread(sum(solvent% vel, dim=2)/solvent% Nmax, 2, solvent% Nmax)
  solvent% force = 0
  solvent% species = 1

  call solvent_cells%init(L, 1.d0)

  call solvent% random_placement(solvent_cells% edges, colloids, solvent_colloid_lj, state(1))

  call solvent% sort(solvent_cells)

  call neigh% init(colloids% Nmax, 1000)

  skin = 1.
  call neigh% make_stencil(solvent_cells, sigma_cut+skin)

  tau = 1.
  N_MD_steps = 200
  dt = tau / N_MD_steps

  n_extra_sorting = 0

  call neigh% update_list(colloids, solvent, sigma_cut + skin, solvent_cells)

  e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
  e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)
  solvent% force_old = solvent% force
  colloids% force_old = colloids% force

  write(*,*) ''
  write(*,*) '    e co so     |   e co co     |   kin co      |   kin so      |   total       |   temp        |'
  write(*,*) ''

  solvent% pos_old = solvent% pos
  colloids% pos_old = colloids% pos
  do i = 1, 100
     so_max = 0
     co_max = 0
     md_loop: do j = 1, N_MD_steps
        !$omp parallel do private(k, tmp_x)
        do k = 1, solvent% Nmax
           solvent% pos(:,k) = solvent% pos(:,k) + dt * solvent% vel(:,k) + dt**2 * solvent% force(:,k) / 2
        end do

        ! Extra copy for rattle
        colloids% pos_rattle = colloids% pos
        do k=1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + &
                dt**2 * colloids% force(:,k) / (2 * colloids% mass(colloids%species(k)))
        end do
        call rattle_body_pos(colloids, links, links_dist, dt, solvent_cells% edges, 1d-12)


        do k = 1, colloids% Nmax
           colloids% pos(:,k) = colloids% pos(:,k) + dt * colloids% vel(:,k) + dt**2 * colloids% force(:,k) / (2 * colloids%mass(1))
        end do
        so_max = solvent% maximum_displacement()
        co_max = colloids% maximum_displacement()

        if ( (co_max >= skin*0.2) .or. (so_max >= skin*0.8) ) then
           call apply_pbc(solvent, solvent_cells% edges)
           call apply_pbc(colloids, solvent_cells% edges)
           call solvent% sort(solvent_cells)
           call neigh% update_list(colloids, solvent, sigma_cut + skin, solvent_cells)
           solvent% pos_old = solvent% pos
           colloids% pos_old = colloids% pos
           n_extra_sorting = n_extra_sorting + 1
        end if

        call switch(solvent% force, solvent% force_old)
        call switch(colloids% force, colloids% force_old)
        solvent% force = 0
        colloids% force = 0
        e1 = compute_force(colloids, solvent, neigh, solvent_cells% edges, solvent_colloid_lj)
        e2 = compute_force_n2(colloids, solvent_cells% edges, colloid_lj)

        !$omp parallel do private(k)
        do k = 1, solvent% Nmax
           solvent% vel(:,k) = solvent% vel(:,k) + dt * ( solvent% force(:,k) + solvent% force_old(:,k) ) / 2
        end do
        colloids% vel = colloids% vel + dt * ( colloids% force + colloids% force_old ) / (2 * colloids%mass(1))
        call rattle_body_vel(colloids, links, links_dist, dt, solvent_cells% edges, 1d-9)
        

        write(15,*) colloids% pos + colloids% image * spread(solvent_cells% edges, dim=2, ncopies=colloids% Nmax)

     end do md_loop

     call apply_pbc(solvent, solvent_cells%edges)
     call apply_pbc(colloids, solvent_cells%edges)
     call solvent% sort(solvent_cells)
     call neigh% update_list(colloids, solvent, sigma_cut + skin, solvent_cells)
     solvent% pos_old = solvent% pos
     colloids% pos_old = colloids% pos

     call simple_mpcd_step(solvent, solvent_cells, state)

     write(*,'(6f16.3)') e1, e2, colloids%mass(1)*sum(colloids%vel**2)/2, sum(solvent% vel**2)/2, &
         e1+e2+colloids%mass(1)*sum(colloids%vel**2)/2+sum(solvent% vel**2)/2, &
         compute_temperature(solvent, solvent_cells)

  end do

  write(*,*) 'n extra sorting', n_extra_sorting

  call h5close_f(error)

end program setup_simple_colloids
