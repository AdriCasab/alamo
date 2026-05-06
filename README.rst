MMW Thermal Spallation In ALAMO
===============================

This fork extends `ALAMO <https://github.com/solidsgroup/alamo>`_ with a new
research integrator for millimetre-wave (MMW) drilling and thermal spallation
of rock.  The model couples enthalpy heat conduction, a Gaussian/Beer-Lambert
MMW source, prescribed-temperature flame-spallation boundary conditions,
mineral microstructure, thermoelastic stress, damage, grain-boundary cohesive
diagnostics, and first-pass spall detachment/removal.

The main executable is built from:

- ``src/mmwspalling.cc`` - dedicated executable entry point.
- ``src/Integrator/MMWSpalling.H`` - top-level multiphysics integrator.

The generated binary is:

.. code-block:: bash

   bin/mmwspalling-3d-g++


What This Project Adds
----------------------

Thermal model
~~~~~~~~~~~~~

- Conserved enthalpy state ``H`` with solid/liquid/vapour phase fractions.
- Explicit and implicit thermal paths for material-property validation.
- Temperature-dependent material tables with reusable ``H <-> T`` inversion.
- Gaussian/Beer-Lambert MMW beam source with time-varying power schedules.
- Radiation and convection losses on surface cells.
- Hu-style prescribed-temperature circular surface patch for flame-spallation
  validation.

Relevant files:

- ``src/Integrator/MMWSpalling.H``
- ``src/Numeric/MMWBeam.H``
- ``src/Numeric/Material/Material.H``
- ``src/Numeric/Material/Constant.H``
- ``src/Numeric/Material/Zhang.H``
- ``src/Numeric/Material/Table.H``


Microstructure and AMR
~~~~~~~~~~~~~~~~~~~~~~

- Voronoi and expression-based mineral microstructures.
- Per-phase material fields for ``kappa``, ``rho``, ``Cp``, ``beta``, ``E``,
  ``mu``, and damage parameters.
- Explicit grain topology:
  ``grain_id``, ``is_grain_boundary``, ``is_phase_boundary``, and legacy
  ``is_gb`` compatibility output.
- Heterogeneous conductivity and damage-modified conductivity.
- AMR repair for discrete phase/grain fields and derived material fields.
- Quartz alpha-beta transformation strain on the microstructure mechanics path.

Primary implementation:

- ``src/Integrator/MMWSpalling.H``


Mechanics, Damage, and Spallation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Heterogeneous thermoelastic mechanics coupled to the temperature field.
- Top-surface prescribed-temperature eigenstrain coupling for Hu validation.
- Drucker-Prager failure criterion and continuous irreversible damage
  evolution.
- Hu breakage indicator ``f_b = sigma_v / sigma_s`` from thermoelastic stress.
- Bilinear grain-boundary cohesive-zone law with irreversible history and
  diagnostics.
- Level-set-style ``phi`` field for surface advancement and first-pass spall
  detachment/removal.

Relevant files:

- ``src/Integrator/MMWSpalling.H``
- ``src/Numeric/DruckerPrager.H``
- ``src/Numeric/CohesiveZone.H``


Validation and Regression Tests
-------------------------------

The MMWSpalling tests live in ``tests/MMWSpalling/``.  They are staged so each
piece of the model is verified before the fully coupled thermal-spallation
pipeline is used.

Implemented test directories include:

- ``skeleton`` - executable/integrator smoke test.
- ``stefan`` - enthalpy and phase-fraction verification.
- ``beam`` - MMW source energy deposition.
- ``equilibrium`` - radiation/convection equilibrium.
- ``zhang_oglesby`` - MMW granite-heating thermal validation.
- ``voronoi`` - mineral microstructure generation.
- ``heterogeneous_kappa`` - heterogeneous conductivity.
- ``heterogeneous_enthalpy`` - per-phase enthalpy/table consistency.
- ``grain_topology`` - grain vs phase boundary flags.
- ``hu_conduction`` - Hu prescribed-temperature conduction validation.
- ``thermal_stress`` - thermoelastic stress regression.
- ``hu_thermoelastic`` - Hu Granite 2 vs Sandstone 2 stress validation.
- ``hu_thermoelastic_amr`` - AMR version of Hu thermoelastic validation.
- ``amr_microstructure_regrid`` - AMR repair for microstructure fields.
- ``dp_yield`` - Drucker-Prager damage evolution.
- ``hu_breakage_index`` - Hu breakage indicator validation.
- ``alpha_beta_transition`` - quartz alpha-beta transformation strain.
- ``gb_cohesive`` - bilinear cohesive-zone law and irreversibility.
- ``spall_event`` - deterministic spall-detachment and surface-advance event.

Most validation tests write an ``output/comparison.png`` figure comparing the
simulation output with the analytical or experimental target.


Build
-----

From the repository root:

.. code-block:: bash

   EIGEN=$PWD/ext \
     CPLUS_INCLUDE_PATH=/opt/homebrew/include \
     LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
     make -j8

If starting from a clean upstream checkout, run ``./configure`` first using the
normal ALAMO workflow for your platform.


Run a Simulation
----------------

Use 4 MPI ranks by default for these test and validation cases on the current
development machine:

.. code-block:: bash

   mpirun --oversubscribe --bind-to none -np 4 \
     bin/mmwspalling-3d-g++ tests/MMWSpalling/<case>/input

Examples:

.. code-block:: bash

   mpirun --oversubscribe --bind-to none -np 4 \
     bin/mmwspalling-3d-g++ tests/MMWSpalling/hu_breakage_index/input_granite2_37

   mpirun --oversubscribe --bind-to none -np 4 \
     bin/mmwspalling-3d-g++ tests/MMWSpalling/spall_event/input


Run a Regression Test
---------------------

The Python postprocessors use ``yt`` and ``numpy``.  On the current development
machine, use:

.. code-block:: bash

   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/<case>/test

Examples:

.. code-block:: bash

   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/hu_breakage_index/test
   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/gb_cohesive/test
   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/spall_event/test


Project Context
---------------

This repository remains an ALAMO-based research code.  The upstream ALAMO
documentation is still the reference for the framework, AMReX setup, mechanics
operators, and general build system:

`ALAMO documentation <https://solidsgroup.github.io/alamo/docs/>`_

For this fork's staged implementation notes and task history, see:

- ``ROADMAP.md`` - short status map.
- ``ACTIVE_STEP.md`` - current coding packet.
- ``ARCHIVE_DONE.md`` - completed-step history.
- ``in-main-tex-you-will-quizzical-treasure.md`` - long historical plan.
- ``main.tex`` - physics proposal and equations.


Citation
--------

If you use the underlying ALAMO framework, cite the ALAMO paper:

.. code-block:: bibtex

   @article{runnels2025alamo,
     title={The Alamo multiphysics solver for phase field simulations with strong-form mechanics and block structured adaptive mesh refinement},
     author={Runnels, Brandon and Agrawal, Vinamra and Meier, Maycon},
     journal={Journal of Open Source Software},
     year={2025}
   }
