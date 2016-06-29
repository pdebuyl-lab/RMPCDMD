.. _run:

Run the code
============

The execution of RMPCDMD is controlled via a single command-line program, ``rmpcdmd``.
RMPCDMD must be built before using this tool, see the :ref:`install` documentation.

``rmpcdmd run``
---------------

Usage::

    rmpcdmd run program input output seed

Arguments

*  ``program`` one of the simulation program coming with RMPCMD (i.e. single_dimer_pbc)
* ``input`` text file with simulation parameters
* ``output`` filename for the output data
* ``seed`` a signed 64-bit integer value. ``seed`` can be set to the value ``auto``, a seed
  will the be generated from the ``/dev/urandom`` device.

When run with no argument, ``rmpcdmd run`` will list the parameters and the possible values
for ``program``.

``rmpcdmd plot``
----------------

Usage::

    rmpcdmd plot [-h] datafile [--obs OBSERVABLE] [--traj GROUP/TRAJECTORY]

Arguments

* ``-h`` display the the full command-line syntax and exit
*  ``datafile`` a datafile produced by one of the simulation programs

One of ``[--obs OBSERVABLE]`` or ``[--traj GROUP/TRAJECTORY]`` can be given to display an
observable or a trajectory from the file.

``rmpcdmd seeder``
------------------

Usage::

    rmpcdmd seeder

Returns a signed 64-bit integer seed.

``rmpcdmd timers``
------------------

Usage::

    rmpcdmd timers [-h] datafile [--plot]

Arguments

* ``-h`` display the the full command-line syntax and exit
* ``datafile`` a datafile produced by one of the simulation programs
* ``--plot`` plots the timers data as a bargraph instead of printing to the terminal.

Prints (or plot in the ``--plot`` option is given) the value of the timers in the simulation
file ``datafile``.

``experiments/`` directory
--------------------------

The execution of some RMPCDMD simulations is illustrated in the directory ``experiments/``,
using makefiles for simplicity. An example simulation session is given below

.. code-block:: console

    user@pc$~$ cd /tmp/RMPCDMD/
    user@pc$/tmp/RMPCDMD$ cd experiments/01-single-dimer/
    user@pc$/tmp/RMPCDMD/experiments/01-single-dimer$ ls
    dimer.parameters  Makefile  plot_histogram.py  plot_velocity.py
    ruckner-kapral.parameters
    user@pc$/tmp/RMPCDMD/experiments/01-single-dimer$ make simulation
    /tmp/RMPCDMD/experiments/01-single-dimer/../../build/rmpcdmd run single_dimer_pbc
    dimer.parameters dimer.h5 auto
    RMPCDMD running single_dimer_pbc
    OMP_NUM_THREADS not set
    Start time -- Thu Jun 16 13:40:08 CEST 2016
    single_dimer_pbc dimer.parameters dimer.h5 3589052620060159831

     Running for         100 loops
     mass   1130.9733867645264        1130.9733867645264     
	5   10   15   20   25   30   35   40   45   50   55   60   65   70   75   80   85
       90   95  100 
     n extra sorting         747

    real    3m25.006s
    user    10m13.496s
    sys     0m0.988s
    End time -- Thu Jun 16 13:43:33 CEST 2016
    205s elapsed
    user@pc$/tmp/RMPCDMD/experiments/01-single-dimer$ 
