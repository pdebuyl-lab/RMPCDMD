.. _tutorial:

Tutorial
========

This tutorial introduces the reader to particle-based simulations of nanomotors. The main
simulation method consists of a coupled scheme of Multiparticle Collision Dynamics (MPCD),
for the fluid, and Molecular Dynamics (MD), for the motor. Chemical activity is introduced
in the fluid via a surface-induced catalytic effect and bulk kinetics.

To start working with RMPCMD, visit the :ref:`install` section first.

Introduction
------------

Nano- to micro-meter scale devices that propel themselves in solution
are built and studied since about a decade. They represent a promise of
future applications at scales where usual control strategies reach their
limits and, ideally, autonomous action replaces manual control.

It is possible to compute, via phoretic theory, the stationary regime of
operation of nanomotors in simple geometries. Still, testing geometrical
effects or including fluctuating behaviour is best done using numerical
simulations. A successful modeling strategy was started by Rückner and
Kapral in 2007 :cite:`ruckner_kapral_prl_2007`.
It builds on a particle-based fluid, explicitly
the Multiparticle Collision Dynamics (MPCD) algorithm. The flexibility
of particle-based simulations allowed for numerous extensions of their
work to Janus particles, polymer nanomotors, various chemical kinetics,
thermally active motors, among others.

The principle of the hybrid scheme is very close to full Molecular
Dynamics (MD), with the major difference that solvent-solvent
interactions are not explicitly computed and are replaced by cell-wise
collisions at fixed time intervals. This saves computational time and
renders otherwise untractable problems feasible.

In Ref.:cite:`ruckner_kapral_prl_2007`, the authors introduce a computational model for a
dimer nanomotor
that is convenient thanks to its simple geometry. There are two spheres
making up the dimer, linked by a rigid bond, one of which being
chemically active and the other not. Solvent particles in contact with
the chemically active sphere are converted from product to fuel. The
active sphere thus acts as a sink for reagent particles (the "fuel") and
a source for product particles.

Preliminary remarks
-------------------

This tutorial does not intend to cover all *possible* manners to conduct
nanomotor simulations. Rather, it aims at presenting one strategy for
modeling chemically powered nanomotors, that is the combination of a
chemically active MPCD fluid coupled to possibly catalytic colloidal
beads.

While this tutorial relies on the RMPCDMD software, it is worth mentioning that Peter
Colberg has also made his software for dimer nanomotor simulations available openly
:cite:`colberg_nanodimer_web`. ``nano-dimer`` is based on OpenCL and can benefit from GPU
acceleration.

Literature
----------

The MPCD algorithm was introduced in :cite:`malevanets_kapral_mpcd_1999` and
:cite:`malevanets_kapral_mpcd_2000`. General reviews on the MPCD simulation method are
available in the literature.

- Raymond Kapral, *Multiparticle collision dynamics: simulation of complex
  systems on mesoscales*, Adv. Chem. Phys. **140**, 89 (2008).
  :cite:`kapral_adv_chem_phys_2008`

- G. Gompper, T. Ihle, D. M. Kroll and R. G. Winkler, *Multi-Particle Collision Dynamics: A
  Particle-Based Mesoscale Simulation Approach to the Hydrodynamics of Complex Fluids*,
  Adv. Polymer Sci. **221**, 1 (2008).
  :cite:`gompper_et_al_adv_polym_sci_2008`


An overview of chemically powered synthetic nanomotors has been published by Kapral
:cite:`kapral_perspective_jcp_2013`.

There is no literature on the practical conduct of nanomotor simulations, however.

The MPCD fluid and Molecular Dynamics
-------------------------------------

A MPCD fluid consists of point particles with a mass (set to unity here
for convenience), a position :math:`x` and a velocity :math:`v`. The
particles evolve in two step: (i) indepedent streaming of the particles
for a duration :math:`\tau` and (ii) cell-wise collision of the
particles’ velocities.

For particle :math:`i` this results in the following equations:

.. math::
   :label: stream

   x_i' = x_i + v_i \tau

and

.. math::
   :label: collide

   v_i' = v_\xi + \omega_\xi ( v_i - v_\xi )

where the prime denotes the quantities after the corresponding step, :math:`\xi` is a cell,
:math:`\omega_\xi` is a rotation operator and :math:`v_\xi` is the center-of-mass velocity
in the cell. The cell consists in a regular lattice of cubic cells in space. Equations
:eq:`stream` and :eq:`collide` conserve mass, energy and linear momentum.

The viscosity for a MPCD fluid can be computed from its microscopic
properties:

.. math::

     \eta = \frac{k_BT\tau\rho}{2m} \left( \frac{ 5\gamma - (\gamma-1+e^{-\gamma})(2-\cos\alpha-\cos 2\alpha) }{(\gamma - 1 + e^{-\gamma})(2-\cos\alpha-\cos 2\alpha)} \right) + \frac{m}{18 a \tau} (\gamma -1 + e^{-\gamma})(1-\cos\alpha)


One can embed a body in a MPCD fluid by using a explicit potential
energy. Then, the streaming step is replaced by the velocity-Verlet
integration scheme. Collision involve only fluid particles and not the
colloid.

The dimer nanomotor
-------------------

Physical setup
^^^^^^^^^^^^^^

In this section, we review the propulsion of the dimer nanomotor
presented by Rückner and Kapral. The geometry of the motor and the
chemical kinetics are presented in the figure below.

The solvent consists of particles of types A and B, initially all
particles are set to A (the fuel). Fuel particles that enter the
interaction range of the catalytic sphere are flagged for reaction but
the actual change of A to B only occurs when the solvent particle is
outside of any interaction range. Else, the change would generate a
discontinuous jump the in the potential energy and disrupt the
trajectory. This chemical activity generates an excess of product
particles "B" around the catalytic sphere and a gradient of solvent
concentration is established.

.. figure:: simple_dimer.png

   Geometry and chemistry for the dimer nanomotor. The graph sketched below represents the
   local excess of “B” particles that is asymmetric for the “N” sphere. Many more “A” and
   “B” particles not shown.

In this type of simulation, the total energy is conserved but the system is maintained in
nonequilibrium by *refueling*, that is by changing B particles to species A when they are
far enough from the colloid.

The solvent and colloids interact via a purely repulsive Lennard-Jones potential of the form

.. math:: V(r) = 4 \epsilon \left( \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} - 1 \right)

where :math:`\epsilon` and :math:`\sigma` can be different depending on the combination of
solvent and colloid species.

Simulation setup
^^^^^^^^^^^^^^^^

Within RMPCDMD, the simulation program for the dimer is called ``single_dimer_pbc``. This
program requires a configuration file that contains the physical parameters, an example of
which is given in the listing below.

::

    # physical parameters
    T = .16666666
    L = 32 32 32
    rho = 9
    tau = 1.0
    probability = 1.0

    # simulation parameters
    N_MD = 200
    N_loop = 50

    # interaction parameters
    sigma_N = 4.0
    sigma_C = 2.0

    d = 6.8
    epsilon_N = 1.0 0.1
    epsilon_C = 1.0 1.0

    epsilon_N_N = 1.0
    epsilon_N_C = 1.0
    epsilon_C_C = 1.0

The configuration allows one to set the size of both spheres in the dimer as well as the
interaction parameters. The setting ``epsilon_N`` contain the prefactor to the Lennard-Jones
potential for the "N" sphere and all solvent species on a single line. In this example, all
the interaction parameters are set to 1 except for the interaction between the "N" sphere
and "B" solvent particles, as was done in :cite:`ruckner_kapral_prl_2007`.

Running the simulations
^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Make sure that you have built the code properly (see :ref:`install`) and that the
          command-line tool ``rmpcdmd`` is available at your command-line prompt. You will
          also need a working scientific Python environment (see :ref:`install_python`).

An example simulation setup is provided in the directory ``experiments`` of RMPCDMD. There,
the sub-directory ``01-single-dimer`` contains a parameter file.

Review the parameters in the file ``dimer.parameters`` then execute the
code

.. code:: bash

    make dimer.h5

The actual commands that are executed will be shown in the terminal.

Analyzing the data
^^^^^^^^^^^^^^^^^^

The output of the simulation is stored in the file ``dimer.h5``, that follows the H5MD
convention for storing molecular data :cite:`h5md_cpc_2014`. H5MD files are regular HDF5
files and can be inspected using the programs distributed by the HDF Group. Issue the
following command and observe the output:

.. code:: bash

    h5ls dimer.h5

HDF5 files have an internal directory-like structure. In ``dimer.h5``
you should find

::

    fields                   Group
    h5md                     Group
    observables              Group
    parameters               Group
    particles                Group
    timers                   Group

The elements are called "groups" in HDF5 terminology. Here, there is data about the
particles (positions, velocities, etc), observables (e.g.  temperature) and fields (here,
the histogram of "B" particles). The ``h5md`` group contains metadata (simulation creator,
H5MD version, etc.), the ``timers`` group contains timing data that is collected during the
simulation and ``parameters`` contains all the parameters with which the simulation was run.

The command

.. code:: bash

    h5ls -r dimer.h5

will visit all groups recursively. The output is then rather large. Let
us focus first on the velocity of the dimer, it is located at
``/particles/dimer/velocity``, where it is stored in ``value`` and the
time step information of the dataset is stored in ``step`` and ``time``.
In the present case, the velocity is sampled at regular time interval of
100 timesteps or equivalently 1 in units of :math:`\tau`.

All the data analysis in this tutorial is done using the Python language
and a set of libraries: NumPy for storing and computing with array data,
h5py for reading HDF5 files, matplotlib for plotting and SciPy for some
numerical routines. For installation, see appendix [install-py]. Some
generic programs are provided with as an introduction to reading the
files, such as ``h5md_plot.py``. Its usage is

.. code:: bash

    python h5md_plot.py dimer.h5 --obs temperature

(the ``obs`` option is preceded by two dashes) to display the
temperature in the course of time. This program can also display the
trajectory of the dimer

.. code:: bash

    python h5md_plot.py dimer.h5 --traj dimer/position

TODO MSD

Nanodimer in a flow
-------------------

Let’s consider a nanodimer moves in a square channel, where periodic
boundary in the :math:`x` directions and real walls in the :math:`y` and
:math:`z` directions are used. The dimer motor interacts with the walls
through long-ranged soft potentials, which restrict the dimer motion to
occur largely along the :math:`x`-direction. The flow is generated in
the :math:`-x`-direction, which is in the opposite direction to the
motor moving direction, by imposing a constant external force with
strength :math:`g` on each solvent molecule. Since the motor is moving
against the flow, as expected, if :math:`g` increases the motor speed
:math:`V_z` decreases and starts to move backward when :math:`g` is
larger than the critical value :math:`g_c`.

An example simulation setup is provided in the directory
``02-chemotactic-cell`` in ``experiments``. The parameters is listed in
listing [dimer-in-a-flow]. To run the simulation, use
``make simulation``, and check the propulsion speed :math:`V_z` with

.. code:: bash

    python plot_velocity.py chemotactic_cell.h5 --directed

Try different values of :math:`g` to see how :math:`V_z` changes with
flow strength :math:`g`.

::

    h5md_file = chemotactic_cell.h5

    # simulation parameters
    N_MD = 50
    N_loop = 1000
    probability = 1

    # number of initialisation steps
    steps_fixed = 100

    # cell parameters
    g = 0.001
    buffer_length = 10
    randomisation_length = 5
    L = 50 50 15
    seed = 4519199302125082433

    max_speed =  0.090

    # fluid parameters
    rho = 10
    T = 1
    tau = 0.5

    # dimer parameters
    sigma_C = 2
    sigma_N = 2
    d = 4.5

    epsilon_C = 1 1 1
    epsilon_N = 1 0.5 1

    # order (T: CN ; F: NC)
    order = F

    # local concentration
    number_of_angles = 6

The Janus nanomotor
-------------------

Physical setup
^^^^^^^^^^^^^^

A Janus motor is a single sphere with an active hemisphere on one side
and an inactive part on the other side (see Fig. [fig:JP]). The
propulsion velocity along its axis :math:`\hat{z}` for
diffusion-controlled reaction is known to be

.. math::

   V_z = c_1 \frac{k_B T}{\eta} \frac{\rho}{3 R} \Lambda,
   \label{eq:Vz}

where :math:`k_B T` is the thermal energy of the system with temperature
:math:`T`, :math:`\eta` is the solvent viscosity, and :math:`\rho` is
the solvent density. The Janus motor has radius :math:`R`, and the
effects due to interactions with the fuel and waste molecules are taken
into account in the factor :math:`\Lambda`. The coefficient :math:`c_1`
depends on the steady-state concentration of product :math:`B`
particles, which is affected by the way of refueling.

Bulk reaction
^^^^^^^^^^^^^

To keep motor active, one needs to maintain the system in a
nonequilibrium state by removing product molecules from and adding fuel
molecules into the system. While in experiments this is achieved by
adding fuel molecules at distant boundaries, in cells waste molecules
may be converted back to fuel molecules through chemical reactions
carried out by proteins or enzymes (:math:`E`). Here we aim to model the
later. Let :math:`n_A` and :math:`n_B` be the concentration of :math:`A`
and :math:`B` molecules, respectively, and :math:`n_E` be the
concentration of the proteins that carry out the irreversible reaction
:math:`B + E \to A + E` with reaction rate :math:`k`. The rate equation
of :math:`A` molecules is

.. math::

   \frac{d n_A}{dt} = k n_E n_B = k_2 n_B.
   \label{eq:rate_eq_A}

The enzyme concentration :math:`n_E` is a constant since enzymes only
facilitate the reaction, therefore one can rewrite
Eq. ([eq:rate\ :sub:`e`\ q\ :sub:`A`]) as
:math:`B \stackrel{k_2}\rightarrow A` with an effective reaction rate
:math:`k_2 = k n_E`.

In reactive multiparticle collision dynamics (RMPCD), reactive and
non-reactive collisions occurs at discrete time interval :math:`\tau`.
In each collision step, the reaction
:math:`B \stackrel{k_2}\rightarrow A` is carried out locally within each
collision cell. Specifically, in cell :math:`\xi` a :math:`B` molecule
is randomly picked from the :math:`n_B^{\xi}` product molecules in the
cell, and is converted to :math:`A` particle with probability
:math:`p = 1-e^{-k_2 n_B^{\xi} \tau}`. The code for the bulk reaction is
shown in listing [bulk\ :sub:`r`\ eaction], which can be found in the
subroutine ``bulk_reaction`` in ``src/mpcd.f90``.

::

    do cell_idx = 1, c%N
       if ( (c%cell_count(cell_idx) <= 1) .or. .not. c%is_reac(cell_idx) ) cycle

       start = c%cell_start(cell_idx)
       n = c%cell_count(cell_idx)

       local_rate = 0
       do i = start, start + n - 1
          s = p%species(i)
          if (s==from) then
             local_rate = local_rate + 1
             pick = i
          end if
       end do

       local_rate = local_rate*rate
       if (threefry_double(state(thread_id)) < (1 - exp(-local_rate*tau))) then
          p%species(pick) = to
       end if
    end do

To run the simulation to test bulk reaction, use

.. code:: bash

    ./setup_bulk_decay

in the directory ``/build``, and an exponential fit to the data can be
done with

.. code:: bash

    python plot_species_evolution.py bulk_decay.h5 --tau 1.0 --species 1 --rate 0.01

Composite Janus motor
---------------------

The Janus motor can propel itself powered by chemical reactions on the
active hemisphere surface. However, in the presence of thermal noises
the Janus motor changes its moving direction by rotational Brownian
motion. It is not possible to simulate Janus particle as a single sphere
interacting with the surrounding solvent molecules only through central
potentials, :math:`V(r)`. It is because the collisions described by
:math:`V(r)` only exchange momentum in the radial direction giving rise
to a body force, but no momentum exchange in the tangential direction so
that the Janus particle can not rotate. In 2013, Pierre and Kapral
ntroduced a composite model for Janus motor, see Fig. [fig:JP](b) . The
active (blue, :math:`C`) and inactive (red, :math:`N`) parts are
composed of spheres linked by rigid bonds. These spheres have the same
radius :math:`1`, and interact with the surrounding solvent particles
through :math:`V_{\alpha C}` and :math:`V_{\alpha N}`.

An example simulation setup is provided in the directory
``03-single-janus`` in ``experiments``. The parameters is listed in
listing [janus-parameters]

::

    # physical parameters
    T = .333333333
    L = 32 32 32
    rho = 9
    tau = 1.0
    probability = 1

    # simulation parameters
    N_MD = 50
    N_loop = 50
    seed = -9223372036854775808
    h5md_file = janus.h5

    # interaction parameters
    sigma_colloid = 1
    epsilon_colloid = 1

    sigma = 3
    epsilon_N = 1.0 0.5
    epsilon_C = 1.0 0.5

    epsilon_N_N = 1.0
    epsilon_N_C = 1.0
    epsilon_C_C = 1.0
    bulk_rate = 0.001

To run the simulation, use ``make simulation``, and check the proplusion
speed :math:`V_z` with

.. code:: bash

    python plot_velocity.py janus.h5 --directed

Controls of motor speed
^^^^^^^^^^^^^^^^^^^^^^^

In Eq. ([eq:Vz]), one can see the propulsion speed is determined by
factors, such as system temperature, fluid properties (viscosity and
density). The effects from concentration gradient of product particle is
given in the coefficient :math:`c_1` which can be altered by the bulk
reaction rate :math:`k_2`. While the factors above affect propulsion
speed, the moving direction is only determined by the factor
:math:`\Lambda` that accounts for the effect from the interactions with
the solvent species. In this section, we will try to explore the effects
from bulk reaction rate :math:`k_2` and interaction with the solvent
:math:`\Lambda`.

Example :math:`1`, forward moving Janus motor.

::

    epsilon_N = 1.0 0.5
    epsilon_C = 1.0 0.5
    bulk_rate = 0.001

Example :math:`2`, backward moving Janus motor.

::

    epsilon_N = 0.5 1.0
    epsilon_C = 0.5 1.0
    bulk_rate = 0.001

Example :math:`3`, Bulk reaction rate.

::

    epsilon_N = 1.0 0.5
    epsilon_C = 1.0 0.5
    bulk_rate = 0.0001
