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
three tiers by purpose.  See ``tests/MMWSpalling/README.md`` for run
instructions; the full per-test list follows below.

.. code-block:: text

   tests/MMWSpalling/
     unit/          fast synthetic / analytic checks of one component each
     validation/    quantitative comparison against published data/theory
       zhang/         MMW source + thermal model         (Zhang & Oglesby)
       hu/            prescribed-T thermal spallation     (Hu, damage_law branch)
       kant/          LEFM Sp onset + confining pressure  (Kant 2017, sp_weibull)
       rossi/         spall depth distribution            (Rossi 2018, sp_weibull)
     extensions/    model-extension demos - NOT validations (no reference data)

Most tests write an ``output/comparison.png`` figure comparing the simulation
against the analytical or experimental target.

Unit tests
~~~~~~~~~~~

Fast, synthetic or analytic; each isolates one building block and runs in
seconds (``tests/MMWSpalling/unit/``):

- ``stefan`` - enthalpy formulation + phase change vs the analytic Stefan
  moving-front solution.
- ``beam`` - Gaussian / Beer-Lambert MMW source: domain-integrated absorbed
  energy and in-plane beam centring.
- ``equilibrium`` - radiation + convection surface losses vs the analytic
  steady-state energy balance.
- ``convective_patch`` - convective flame-jet Robin BC vs the Carslaw & Jaeger
  semi-infinite-solid analytic solution.
- ``voronoi`` - Voronoi mineral-microstructure generation (phase fractions and
  topology).
- ``grain_topology`` - correctness of the ``grain_id``, ``is_grain_boundary``,
  and ``is_phase_boundary`` flags.
- ``heterogeneous_kappa`` - per-phase and damage-modified conductivity vs the
  series-resistance temperature profile.
- ``heterogeneous_enthalpy`` - per-cell ``H = rho*Cp*(T - T_ref)`` consistency
  on the microstructure path.
- ``thermal_stress`` - heterogeneous thermoelastic stress regression.
- ``dp_yield`` - Drucker-Prager yield, irreversible damage evolution, and
  stiffness/conductivity degradation.
- ``gb_cohesive`` - bilinear grain-boundary cohesive-zone law: envelope,
  irreversibility, and fracture energy.
- ``alpha_beta_transition`` - quartz alpha-beta transformation eigenstrain.
- ``amr_microstructure_regrid`` - AMR repair of discrete microstructure and
  derived fields across a regrid (raw per-level check).
- ``spall_event`` - deterministic spall detachment and surface advance: ``phi``
  shift, detached-cell reset, surface-mask migration, closed-form ``h_spall``.
- ``regime_low_high_power`` - spall-vs-vaporisation regime selection and the
  unified rate-of-penetration (RoP).
- ``sp_lefm`` - LEFM Sp machinery on three prescribed-stress inputs: the K_I
  Tada-weight integrator, the K_Ic(T) Nasseri table closed form, depth-resolved
  per-cell K_I, the Weibull flaw distribution (including the ``weibull.V0``
  mesh-objectivity size effect), per-cell Sp variation, determinism, sign
  convention, and AMR regrid-repair.
- ``sp_v_n_regime`` - per-cell normal velocity, the ``regime_field`` diagnostic,
  and the local-Q (Gaussian) RoP reconstruction.

Validation tests
~~~~~~~~~~~~~~~~~

Quantitative comparisons against published experiments/theory
(``tests/MMWSpalling/validation/``).  Each group names the paper it reproduces.

**Zhang et al. (2023) / Oglesby et al. (2014) - MMW granite-heating thermal
model.**  Reproduces the Oglesby flat-granite millimetre-wave heating
experiment as modelled by Zhang et al.: a time-varying MMW source on a granite
surface, comparing the surface temperature beneath the beam centre (Zhang
Fig. 7, measured with a 137 GHz radiometer and emissivity-corrected).
Exercises the temperature-dependent ``rho``/``kappa``/``Cp`` material model,
the beam source, surface losses, and the H<->T inversion.

- ``zhang_oglesby`` - **DEFERRED**: long-running and not in the routine
  rotation; run it to obtain a current verdict.

**Hu et al. (2019) - lowest required surface temperature (LRST) for thermal
spallation.**  Reproduces Hu et al., "Lowest Required Surface Temperature for
Thermal Spallation in Granite and Sandstone Specimens: Experiments and
Simulations" (*Rock Mechanics and Rock Engineering*, 2019), which measured the
surface temperature and onset time at which a heated rock surface first spalls,
for Granite 2 and Sandstone 2.  This chain builds the prescribed-temperature
thermal-spallation pipeline on the ``damage_law`` branch:

- ``hu_conduction`` - prescribed-temperature surface-patch conduction
  (Granite 2 / Sandstone 2).  **PASS**.
- ``hu_thermoelastic`` - Granite 2 vs Sandstone 2 thermoelastic surface stress,
  single level.  **PASS** (documented bottom-clamp / L3 ordering caveats).
- ``hu_thermoelastic_amr`` - AMR version of the above.  **PASS** (same caveats).
- ``hu_breakage_index`` - Hu breakage indicator ``f_b = sigma_v / sigma_s``.
  **PASS** (surface targets; deep L6 relaxed under the bottom clamp).
- ``hu_spall_onset`` - LRST/onset via the DP-damage threshold proxy
  (homogeneous, removal off).  **PASS** (v1; sandstone damage-depth
  warning-only).
- ``hu_end_to_end`` - the canonical, fully coupled Hu run (conduction ->
  thermoelastic -> DP damage -> spall removal).  **PASS** (granite onset
  40.0 s, sandstone 90.5 s; ordering correct).

**Kant (2017) - flame-jet spallation onset, Central Aare granite.**  Reproduces
Kant's Central Aare granite flame-spallation experiments, in which a convective
flame jet heats a laterally-confined granite surface to the spalling onset
temperature.  ALAMO reproduces Kant's Eq. 15 closed-form onset temperature (the
LEFM ``K_I = K_Ic`` criterion) and its confining-pressure dependence, on the
``sp_weibull`` branch:

- ``sp_kant_onset`` - onset at confining pressure ``p = 0``; FEM vs Kant Eq. 15
  plus Weibull patchiness.  **PASS** (6/6 targets; verification onset 461.3
  degrees C, inside Kant's 390-560 degrees C band).
- ``sp_kant_pressure_sweep`` - confining-pressure dependence
  ``p in {0, 27, 48} MPa``.  **PASS** (monotone decreasing; within 8% of
  Eq. 15 at each pressure).

**Rossi et al. (2018) - spall crack-depth distribution.**  Reproduces Rossi's
post-cooling measurement of the spatial distribution of distinct thermal-spall
cracks versus depth (a peak at 100-200 um, falling to baseline by ~520 um), on
the ``sp_weibull`` branch:

- ``sp_rossi_damage_profile`` - **FAIL / DEFERRED**: the Sp depth-scan
  mechanics gates (P1-P4) PASS, but the predicted crack-depth peak lands at
  ~705 um rather than Rossi's 100-200 um.  This is a mesh-resolution limit, not
  a code bug (see ``tests/MMWSpalling/validation/rossi/rossi-validation-diagnostic-design.md`` section 11).  The
  test is kept as an ongoing P1-P4 mechanics regression; the depth-band gate
  (R1) is reported but not enforced.

Extensions
~~~~~~~~~~

Model-extension demos with **no published reference**, so they are sanity /
diagnostic runs, not validations - a pass means "ran and evolved sensibly"
(``tests/MMWSpalling/extensions/``):

- ``hu_spall_onset_voronoi`` - Voronoi-microstructure granite spall-onset and
  drilling on the DP ``damage_law`` branch (heterogeneous counterpart to the
  homogeneous ``hu_spall_onset`` validation).
- ``hu_spall_onset_mmwbeam`` - Voronoi granite drilled by a 50 kW MMW beam
  (``damage_law``); prints a drilling-progress diagnostic.
- ``hu_spall_onset_mmwbeam_lefm`` - ``sp_weibull`` LEFM spall under an MMW beam
  (relaxes the prescribed-T gate so sigma_xx is read from ``stress_mf``);
  binding code-correctness gates only, drilling metrics informational.


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
- ``docs/project/ARCHIVE_DONE.md`` - completed-step history.
- ``docs/project/in-main-tex-you-will-quizzical-treasure.md`` - long historical plan.
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
