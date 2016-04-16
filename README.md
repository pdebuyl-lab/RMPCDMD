RMPCDMD: Reactive MultiParticle Collision Dynamics - Molecular Dynamics
=======================================================================

**Author:** [Pierre de Buyl](http://pdebuyl.be/)

RMPCDMD is a collection of Fortran modules and programs for the
mesoscopic modeling of chemically active fluids with embedded colloids.

This software represents the development version of the author to
perform his research on nanomotor modeling.

Contributions from Peter Colberg, either in the code or via
discussions, are warmly acknowledged.

## Status

A former version of this code, still available in the branches `trs`
and `trs_two_prod`, was used to obtain the results presented in P. de
Buyl and R. Kapral [Nanoscale 5, 1337-1344
(2013)](http://dx.doi.org/10.1039/C2NR33711H) and P. de Buyl,
A. S. Mikhailov and R. Kapral [EPL 103, 60009
(2013)](http://dx.doi.org/10.1209/0295-5075/103/60009).

The current version is under total refactoring to remove the use of
global variables, enable testing and plan for parallel accelerations.

## Compile the code

RMPCDMD has the following requirements:

- A Fortran 2003 compiler (e.g. [gfortran](https://gcc.gnu.org/wiki/GFortran))
- A Fortran enabled [HDF5](https://www.hdfgroup.org/HDF5/) installation
- [CMake](http://cmake.org/)
- [GNU Make](https://www.gnu.org/software/make/)
- [git](http://git-scm.com/)

At the command line, execute the following

    git clone https://github.com/pdebuyl-lab/RMPCDMD
    cd RMPCDMD
    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake ..
    make VERBOSE=1

## Run the code

The most convenient manner to execute a simulation is to visit an "experiment"
directory. From the root of the RMPCDMD software directory

    cd experiments/01-single-dimer
    make simulation

## License

BSD 3-clause, see [LICENSE](LICENSE).

## Contributors

Peter Colberg: general programming improvements, OpenMP, debugging  
Laurens Deprez: single colloid setup, gravity field and corresponding bounce-back, shake/rattle for dimers
