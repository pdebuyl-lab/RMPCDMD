## RMPCDMD Experiments

This directory contains example experiments that can be run with RMPCDMD. They
consist in part of the reproduction of published results and in part of new
developments.

The software must be compiled before running any experiment. See the README.md
file at the root of the RMPCDMD software distribution for more information.

Every experiment directory contains a parameter file and a Makefile to run
preset simulations available in RMPCDMD and example analysis code. For
instance, to run the first experiment use the following commands:

    cd 01-single-dimer
    make simulation

This results in an output file `dimer.h5`. The velocity of the dimer can be
displayed with the command

    python plot_velocity.py dimer.h5

