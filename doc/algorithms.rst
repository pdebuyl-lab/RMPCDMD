Algorithms
==========

MPCD
----

The MPCD algorithim introduced in :cite:`malevanets_kapral_mpcd_1999` is implemented in
:doxytag:`simple_mpcd_step`.

MD
--

Velocity verlet, see :cite:`malevanets_kapral_mpcd_2000`

See :doxytag:`md_pos`

Convenience routines
--------------------

Temperature computation
^^^^^^^^^^^^^^^^^^^^^^^

A general routine to compute the temperature `compute_temperature
<../../html/namespacempcd.html#a0dfdbae6c690cf5d011e7fd0ef78969e>`_ computes cell-wise the
kinetic energy relative to the cell's center-of-mass.

 .. math::
    T = \frac{1}{3N_c} \sum_\xi \frac{1}{N_\xi-1} \sum_i m_i \left( v_i - v_\xi \right)^2

where :math:`N_c` is the number of cells, the variable :math:`\xi` represents a cell,
:math:`N_\xi` the number of particles in a cell, :math:`v_\xi` the center-of-mass velocity
of the cell, and :math:`m_i` and :math:`v_i` the mass and velocity of particle.
