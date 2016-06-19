.. RMPCDMD documentation master file, created by
   sphinx-quickstart on Wed May  4 12:33:14 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. only:: html

    RMPCDMD
    =======

.. only:: latex

    About RMPCDMD
    =============

:Authors: Pierre de Buyl, Peter Colberg, Laurens Deprez
:License: BSD

RMPCDMD is a software for the simulation of colloids via Molecular Dynamics,
embedded in a MPCD fluid.

Ready-to-execute simulation programs are provided for the dimer nanomotor in Periodic
Boundary Conditions (PBC), the forced Poiseuille flow or for N colloids in PBC. These
programs only require the setting of parameters in the ad-hoc text file for execution.

The implementation of the algorithms is organized in *modules* and is available for users to
prepare their own simulation setup, then requiring programming in Fortran.

See :ref:`install` for obtaining and building RMPCDMD. The contact person for RMPCDMD is
`Pierre de Buyl <http://pdebuyl.be>`_.


Features:

  - MPCD collision rule for the solvent
  - Chemical activity (either catalytic at a colloid or in the bulk)
  - Rattle constrained dynamics for rigid bodies
  - Walls (specular, bounce-back, virtual particles)
  - OpenMP multithreaded operation

Development of the code takes place on `GitHub <https://github.com/>`_. The code repository
is `pdebuyl-lab/RMPCDMD <https://github.com/pdebuyl-lab/RMPCDMD>`_

The source code features inline comments, published with Doxygen: `api <api/index.html>`_.

