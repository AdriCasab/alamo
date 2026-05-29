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
- Convective flame-jet Robin surface BC ``q = h_fl (T_flame - T_surface)`` for
  Kant-style spallation validation.

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
- LEFM stress-intensity spall criterion ``Sp = K_I / K_Ic(T)`` using the Tada
  edge-crack weight function, a temperature-dependent fracture-toughness table,
  and a Weibull per-cell flaw distribution.
- Confining-pressure-dependent onset via the ``confining.p`` parser key.
- Lateral / zlo roller boundary conditions for 1-D-confined static
  thermoelastic validation.
- Level-set-style ``phi`` field for surface advancement and first-pass spall
  detachment/removal.
- Vaporisation removal and unified rate-of-penetration (RoP) regime selection
  (spall vs vapour).

Relevant files:

- ``src/Integrator/MMWSpalling.H``
- ``src/Numeric/DruckerPrager.H``
- ``src/Numeric/CohesiveZone.H``
- ``src/Numeric/SpCriterion.H``
- ``src/BC/Operator/Elastic/ZloRoller321.H``


Validation and Regression Tests
-------------------------------

The MMWSpalling tests live in ``tests/MMWSpalling/`` and are organised into
three tiers by purpose.  See ``tests/MMWSpalling/README.md`` for the full
per-test index, including what each test checks and the pass/fail status of
every validation.

.. code-block:: text

   tests/MMWSpalling/
     unit/          fast synthetic / analytic checks of one component each
     validation/    quantitative comparison against published data/theory
       zhang/         MMW source + thermal model         (Zhang & Oglesby)
       hu/            prescribed-T thermal spallation     (Hu, damage_law branch)
       kant/          LEFM Sp onset + confining pressure  (Kant 2017, sp_weibull)
       rossi/         spall depth distribution            (Rossi 2018, sp_weibull)
     extensions/    model-extension demos - NOT validations (no reference data)

- **unit/** - 19 component tests (e.g. ``stefan``, ``beam``,
  ``convective_patch``, ``dp_yield``, ``gb_cohesive``, ``spall_event``,
  ``sp_weibull_unit``).  These isolate one numerical/physical building block
  each and run in seconds; they are the first line of defence.
- **validation/** - literature comparisons.  Currently passing: the Hu chain
  (``hu_conduction``, ``hu_thermoelastic``, ``hu_breakage_index``,
  ``hu_spall_onset``, and the canonical ``hu_end_to_end``) and the Kant 2017
  LEFM cases (``sp_kant_onset``, ``sp_kant_pressure_sweep``).  The Rossi case
  ``sp_rossi_damage_profile`` is deferred (its mechanics gates pass, but the
  Rossi depth band is unmet at the working mesh resolution), and
  ``zhang_oglesby`` is not in the routine rotation.
- **extensions/** - ``hu_spall_onset_voronoi``, ``hu_spall_onset_mmwbeam``, and
  ``hu_spall_onset_mmwbeam_lefm`` are sanity/diagnostic runs of model
  configurations that have no published reference; a pass means "ran and
  evolved sensibly", not scientific validation.

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
     bin/mmwspalling-3d-g++ tests/MMWSpalling/<tier>/<case>/input

Examples (note the ``unit/``, ``validation/<source>/``, ``extensions/`` tiers):

.. code-block:: bash

   mpirun --oversubscribe --bind-to none -np 4 \
     bin/mmwspalling-3d-g++ tests/MMWSpalling/validation/hu/hu_breakage_index/input_granite2_37

   mpirun --oversubscribe --bind-to none -np 4 \
     bin/mmwspalling-3d-g++ tests/MMWSpalling/unit/spall_event/input


Run a Regression Test
---------------------

The Python postprocessors use ``yt`` and ``numpy``.  On the current development
machine, use:

.. code-block:: bash

   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/<tier>/<case>/test

Examples:

.. code-block:: bash

   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/validation/hu/hu_breakage_index/test
   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/gb_cohesive/test
   /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/spall_event/test

The full per-test index, with what each test checks and the pass/fail status of
every validation, is in ``tests/MMWSpalling/README.md``.


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
