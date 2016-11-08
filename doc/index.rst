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

:Authors: Pierre de Buyl, Peter Colberg, Laurens Deprez, Mu-Jie Huang
:License: BSD 3-clause
:Website: http://lab.pdebuyl.be/rmpcdmd/
:Version: |rmpcdmd_version|



RMPCDMD is a software for the simulation of colloids via Molecular Dynamics, embedded in
a MPCD fluid.
Ready-to-execute simulation programs are provided for the dimer nanomotor in Periodic
Boundary Conditions (PBC), the forced Poiseuille flow or for N colloids in PBC. These
programs only require the setting of parameters in the ad-hoc text file for execution.

Highlights:

- The simulation of dimer nanomotors, reproducing the pioneering work of Rückner and
  Kapral :cite:`ruckner_kapral_prl_2007` is well tested.
- This code is a research code, so other features are probably *under development*. This
  should not prevent you from using it!
- We have a :ref:`tutorial` on nanomotor simulations.

Features:

- MPCD collision rule for the solvent
- Chemical activity (either catalytic at a colloid or in the bulk)
- Rattle constrained dynamics for rigid bodies
- Walls (specular, bounce-back, virtual particles)
- Hilbert curve based spatial sorting of solvent particles
- `H5MD <http://nongnu.org/h5md>`_ trajectory file output :cite:`h5md_cpc_2014`
- Fortran 2008 codebase using modules and *no global variable*
- OpenMP multithreaded operation


Development and contact information:

- Development of the code takes place on `GitHub <https://github.com/pdebuyl-lab/RMPCDMD>`_.
- The contact for RMPCDMD is the main author, `Pierre de Buyl <http://pdebuyl.be>`_.
- Bug reports are welcome either by email or via `GitHub issues
  <https://github.com/pdebuyl-lab/RMPCDMD/issues>`_

The source code features inline comments, published with Doxygen: `api <api/index.html>`_.

The use of the Hilbert curve sorting and of the Threefry Random Number Generator
(:cite:`random123`) is inspired by Peter Colberg's code `nano-dimer`
:cite:`colberg_nanodimer_web`.

.. only:: html

    Download the documentation as a `pdf file <_static/RMPCDMD.pdf>`_ (does not have api
    links).
