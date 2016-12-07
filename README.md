RMPCDMD: Reactive MultiParticle Collision Dynamics - Molecular Dynamics
=======================================================================

**Homepage:** [RMPCDMD](http://lab.pdebuyl.be/rmpcdmd/) (includes documentation)  
**License:** BSD 3-clause, see [LICENSE](LICENSE).

[![Build Status](https://travis-ci.org/pdebuyl-lab/RMPCDMD.svg?branch=master)](https://travis-ci.org/pdebuyl-lab/RMPCDMD)

RMPCDMD is a collection of Fortran modules and programs for the
mesoscopic modeling of chemically active fluids with embedded colloids.

## Status

RMPCDMD is ready to use for chemically powered dimer nanomotor simulations.
A former version of this code, still available in the branches `trs`
and `trs_two_prod`, was used to obtain the results presented in P. de
Buyl and R. Kapral [Nanoscale 5, 1337-1344
(2013)](http://dx.doi.org/10.1039/C2NR33711H) and P. de Buyl,
A. S. Mikhailov and R. Kapral [EPL 103, 60009
(2013)](http://dx.doi.org/10.1209/0295-5075/103/60009).

The current version has been totally refactored to remove the use of global variables,
enable testing and enable OpenMP multithreaded operation.

## Citation

RMPCDMD is presented in the article *RMPCDMD: Simulations of colloids with coarse-grained
hydrodynamics, chemical reactions and external fields*
[[arXiv:1608.04904](https://arxiv.org/abs/1608.04904)].

Please cite this paper if you use RMPCDMD in a research work. A bibtex entry is provided in
the [CITATION](CITATION) file.

## Requirements

RMPCDMD has the following requirements:

- A Fortran 2008 compiler (e.g. [gfortran](https://gcc.gnu.org/wiki/GFortran) ≥ 4.7 with support for [OpenMP](https://gcc.gnu.org/wiki/openmp) ≥ 3.1)
- A Fortran enabled [HDF5](https://www.hdfgroup.org/HDF5/) installation
- [CMake](http://cmake.org/)
- [GNU Make](https://www.gnu.org/software/make/)
- [git](http://git-scm.com/)

See the documentation for
[installation instructions](http://lab.pdebuyl.be/rmpcdmd/install.html) for Linux and OS X.

## Contact

- The contact for RMPCDMD is the main author, [Pierre de Buyl](http://pdebuyl.be/).
- Bug reports are welcome either by email or via
  [GitHub issues](https://github.com/pdebuyl-lab/RMPCDMD/issues).

## Contributors

Peter Colberg: general programming improvements, OpenMP, debugging  
Laurens Deprez: single colloid setup, gravity field and corresponding bounce-back, shake/rattle for dimers  
Mu-Jie Huang: parts of the tutorial
