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

If you encounter a configuration issue, it is often useful to remove all CMake data by
executing::

    rm -r CMake*

in the build directory (*not in the source directory*).

Installation on Linux
---------------------

On a Debian distribution (or derivative, e.g. Ubuntu), as root::

    apt-get install gfortran libhdf5 cmake git

will install the dependencies. (See `install_hdf5`_ for Ubuntu 14.04 or problematic HDF5
installs).

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
not the case, you may define the environment variable ``HDF5_ROOT=/path/to/hdf5`` before
invoking cmake (see the OS X installation notes for an example).

Installation on OS X
--------------------

The first step is to install a development environment, and thus the `XCode
<https://developer.apple.com/xcode/>`_ software from Apple (we have tested RMPCDMD with
XCode 7.3.1 on OS X El Capitan).

We used `MacPorts <https://www.macports.org/>`_ for gcc, git, cmake and make. Install
MacPorts (a GUI installer is available at https://www.macports.org/install.php) and then
the dependencies::

    sudo port install gcc-5 cmake git

The rest is similar to Linux except that the name of the compiler has to be specified
manually at the command-line to prevent the automatic selection of Apple's provided
compiler. Also, HDF5 is built locally::

    git clone https://github.com/pdebuyl-lab/RMPCDMD
    cd RMPCDMD
    git submodule init
    git submodule update
    CC=gcc-mp-5 FC=gfortran-mp-5 ./scripts/download_and_build_hdf5.sh
    mkdir build
    cd build
    HDF5_ROOT=../_hdf5-1.8.17 cmake .. -DCMAKE_C_COMPILER=gcc-mp-5 -DCMAKE_Fortran_COMPILER=gfortran-mp-5
    make

The compiler names given here may vary depending on your setup.

Installation on Windows
-----------------------

Installation on Windows is not tested at this time. The major difficulty is likely related
to a functional installation of HDF5 for Fortran under Windows, that is only provided for
the Intel compiler. Having no Windows computer with the Intel Fortran compiler at our
disposal, user feedback is welcome.

.. _install_hdf5:

HDF5
----

On OS X or with older Linux systems (such as Ubuntu 14.04 still in wide usage), it is
necessary to build HDF5 to enable the Fortran 2003 interface.

A script is provided to download and build HDF5 that has been tested on OS X and on
Ubuntu. It must be run from the main RMPCDMD directory::

    ./scripts/download_and_build_hdf5.sh

with compiler definitions on OS X (see above). This script must be run only once and the
environment variable ``HDF5_ROOT`` must be set when invoking cmake (also see the OS X
installation notes).

The script downloads HDF5 1.8.17 and installs it under ``_hdf5-1.8.17``. You can remove the
directory ``hdf5-1.8.17`` (no leading underscore) after the execution of the script.
