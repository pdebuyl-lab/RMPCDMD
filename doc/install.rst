.. _install:

Install
=======

Requirements
------------

To compile RMPCDMD, the following software is required

  - A Fortran 2003 compiler (e.g. `gfortran <https://gcc.gnu.org/wiki/GFortran>`_ ≥ 4.7 with
    support for `OpenMP <https://gcc.gnu.org/wiki/openmp>`_ ≥ 3.1)
  - A Fortran enabled `HDF5 <https://www.hdfgroup.org/HDF5/>`_ installation
  - `CMake <http://cmake.org/>`_
  - `GNU Make <https://www.gnu.org/software/make/>`_
  - `git <http://git-scm.com/>`_

Other utilities, mostly Fortran code, will be downloaded automatically as part of the
installation process.

RMPCDMD has been developed on Linux, using both the GFortran and Intel Fortran
compilers. Compilation under OS X has been successfully achieved by users of RMPCMD. The
Windows platform is currently not tested.

Installation on Linux
---------------------

On a Debian distribution (or derivative, e.g. Ubuntu), as root::

    apt-get install gfortran libhdf5 cmake git

will install the dependencies.

At the command line, execute the following::

    git clone https://github.com/pdebuyl-lab/RMPCDMD
    cd RMPCDMD
    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake ..
    make VERBOSE=1

This will first fetch RMPCDMD, then the related modules `ParseText`, `fortran_tester`,
`fortran_h5md` and `random_module`. The compilation is prepared by the ``cmake`` command and
executed by the ``make`` command.

If your HDF5 installation is properly setup, CMake should find it automatically. If this is
not the case, you may pass the option ``-DHDF5_ROOT=/path/to/hdf5`` to cmake.

Installation on OS X
--------------------

TODO

Installation on Windows
-----------------------

Installation on Windows is not tested at this time. The major difficulty is likely related
to a functional installation of HDF5 for Fortran under Windows, that is only provided for
the Intel compiler. Having no Windows computer with the Intel Fortran compiler at our
disposal, user feedback is welcome.
