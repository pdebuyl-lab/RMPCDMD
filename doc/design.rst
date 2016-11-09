Design
======

General principles
------------------

General purpose routines and computational algorithms are implemented as Fortran modules in
`src/`. Parts of the code is tested with programs in `test/`.

Actual simulations are found in `programs/`, where every program corresponds to a given
setup. The programs read parameters from a text file and store trajectories in `H5MD
<http://nongnu.org/h5md/>`_ files.

The aim of the code organization is that one can write a program that call trusted routines
from the modules, giving a high-level view on the simulation. There is no global variable in
the code so that a routine may only act on variables that are passed to it.

Storage of cell and particle data
---------------------------------

The storage of particle and cell data is defined in the module `particle_system` and
`cell_system`, respectively.

The derived type `particle_system_t` defines arrays for the position, velocity, etc of every
particle.

The derived type `cell_system_t` defines a cell list with a one-dimensional index that is
mapped to a three-dimensional compact Hilbert index
:cite:`hamilton_compact_hilbert_tr`. Particles are sorted according to this cell system, as
in the nano-dimer software :cite:`colberg_nanodimer_web`.

Operations on particles can proceed in two different ways:

  1. Particle-wise: The loop goes over every particle, following the particle index of the
     particle data in memory.
  2. Cell-wise: The loop runs over successive cells, where an inner loop runs over all
     particles. Thanks to the sorting of particles, the indices of the particles in a given
     cell are contiguous.

Programming language
--------------------

RMPCDMD is written in Fortran 2008 and uses derived types objects to store simulation data.

RMPCDMD is built using `CMake <https://cmake.org/>`_, which facilitates the inclusion of
external libraries (`fortran_h5md <https://github.com/pdebuyl/fortran_h5md>`_,
`random_module <https://github.com/pdebuyl/random_module>`_, `ParseText
<https://github.com/pdebuyl/ParseText>`_ and `fortran_tester
<https://github.com/pdebuyl/fortran_tester>`_). The `HDF5 <https://www.hdfgroup.org/HDF5/>`_
library, with the Fortran 2003 interface, is also necessary to build RMPCDMD.

Random numbers
--------------

A C library implements the Threefry Random Number Generator (RNG) :cite:`random123`. This
RNG can be used in separate independent threads, provided that the thread's numerical ID is
part of the *key* of the RNG. The other part of the key, that is made of two 64-bit
integers, is the seed that is provided by the user.

The RNG is implemented in `random_module <https://github.com/pdebuyl/random_module>`_, that
is fetched as a git submodule in RMPCDMD. A Fortran wrapper, using the ``bind(c)``
attribute, manages the Fortran/C interface.
