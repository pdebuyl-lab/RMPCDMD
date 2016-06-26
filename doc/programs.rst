Programs
========

Configuration files
-------------------

The parameters listed for each program are read from text configuration files. The syntax is given by example for the arguments that follow.

- **T** Temperature
- **L** Box size
- **enable_thermostat** Activate the thermostat

.. code::

    T = 1.2
    L = 32 32 32
    enable_thermostat = F

Boolean values are input as T (True) or F (False). For vectors (such as the box size **L**),
several components are listed, separated by a space. For more examples, the subdirectories
in ``experiments/`` contain configuration files for some simulations.

.. doxyheader:: ../programs/single_dimer_pbc.f90

.. doxyheader:: ../programs/poiseuille_flow.f90

.. doxyheader:: ../programs/chemotactic_cell.f90

