Algorithms
==========

.. toctree::
   :maxdepth: 1

   algorithms/quaternions

MPCD
----

The MPCD algorithm introduced in :cite:`malevanets_kapral_mpcd_1999` is implemented in
:doxytag:`simple_mpcd_step` for periodic boundary conditions.

The presence of walls is taken into account in :doxytag:`wall_mpcd_step` following
:cite:`lamura_mpcd_epl_2001`: ghost particles are placed in the boundary cells, during the
collision, to reach the average density of the solvent. Optionally, a bulk Anderson
thermostat can also be applied in this routine.

Molecular Dynamics
------------------

Molecular Dynamics (MD) is implemented using the Velocity verlet integration scheme
:cite:`malevanets_kapral_mpcd_2000`, where solvent particles and colloids evolve at the MD
timestep :math:`dt`.

The implementation in RMPCDMD is found in :doxytag:`md_pos` and :doxytag:`md_vel`. It is
important to know that particles close to walls have their velocities updated in the
*stream* routines and are skipped in :doxytag:`md_vel`. This is only correct if they are
outside of the interaction range of all colloids, which is the case in all simulations here.

Boundary conditions
-------------------

Bounce-back boundary conditions are used for the walls, in addition to the modified
collision rule. The are presented in :cite:`allahyarov_gompper_mpcd_flows_2002` or
:cite:`whitmer_luitjen_2010` and implemented in :doxytag:`mpcd_stream_xforce_yzwall`. There,
the bounce-back collision is computer for the parabolic trajectories relevant for forced
flows.

Convenience routines
--------------------

Temperature computation
^^^^^^^^^^^^^^^^^^^^^^^

A general routine to compute the temperature :doxytag:`compute_temperature` computes
cell-wise the kinetic energy relative to the cell's center-of-mass.

 .. math::
    T = \frac{1}{3N_c} \sum_\xi \frac{1}{N_\xi-1} \sum_i m_i \left( v_i - v_\xi \right)^2

where :math:`N_c` is the number of cells, the variable :math:`\xi` represents a cell,
:math:`N_\xi` the number of particles in a cell, :math:`v_\xi` the center-of-mass velocity
of the cell, and :math:`m_i` and :math:`v_i` the mass and velocity of particle.
