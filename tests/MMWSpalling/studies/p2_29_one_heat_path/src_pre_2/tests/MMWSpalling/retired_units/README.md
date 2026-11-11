# Retired unit tests

Tests of deleted features, kept as records (not run, not in any gate). One line each.

- `constraint_ref_lock` -- P2-16: locked the constraint depth d (spallation.sp.constraint_depth), deleted in P2-16
- `face_cref` -- P2-16: the face list's constraint reference (P2-14a), deleted with d in P2-16
- `sp_melt_skeleton` -- P2-26: A/B test of spall.melt_aware's skeleton f(Lambda_L) (P1), deleted in P2-26 (census C3, user Q6); last on-disk run sweep_run.log (2026-09-17) PASS, the census (2026-11-05) recorded rc 1 -- not re-run
- `dp_yield` -- P2-27 (census C1a): the Drucker-Prager yield / damage growth (Step 8); the damage lineage is deleted
  (decision 1 of 2026-11-05: LEFM is the model). Its keys (`damage.*`) now abort; not re-run
- `gb_cohesive` -- P2-27 (census C1a): the grain-boundary cohesive zone (Step 9), deleted with the lineage; `cohesive.*` aborts
- `spall_event` -- P2-27 (census C1a): the damage-law removal (D >= threshold, principal-stress ratio, h = C_h sqrt(alpha t),
  two cells per event); it was in the tree gate (8 paths, `hash.sh`) until P2-27 stage 2, its last hashes are in
  `witness/p2_27a.md5`
- `regime_low_high_power` -- P2-27 (census C1a): the damage-law spall vs bulk-averaged vapour regimes (prescribed ratio /
  prescribed RoP), both retired with the lineage
The decks above write to their old `tests/MMWSpalling/unit/<name>/output` paths: they are records, not runnable as moved.
- `heterogeneous_kappa` -- P2-27 stage 3 (census C1b): its heterogeneity WAS the damage IC (`damage.ic.*`, `kappa_eff =
  kappa exp(-ln2 D)`), retired with the lineage; D and kappa_eff are no longer registered

P2-28 (census H): THE COLUMN REMOVAL PATH is deleted (one removal path, the face list). Each test below needed it or one of its keys (which now abort); not re-run. Successor in brackets.
- `robin_pinned` -- P2-28: the pinned Robin closure (robin_form = pinned) and its pin rule [face_heat, face_list]
- `robin_jet` -- P2-28: the jet map on the column path's Robin forms [face_heat (v)/(vi), the face gate]
- `jet_enthalpy` -- P2-28: the jet enthalpy march's column flux (per-bin recursion vs the strip closed form) [face_heat (v)/(vi) + the face gate; the per-bin recursion is a named coverage gap]
- `robin_feet` -- P2-28: the burner on its feet on the column path [face_feet]
- `jet_d2c` -- P2-28: the D2c steady closures on the column path (feed cap, stall, foot_rule modes) [face_feet; the feet descent's feed cap / stall / foot_rule modes are a named coverage gap]
- `pad_clear_ledger` -- P2-28: the column path's pad-clearance ledger [face_feet (e)]
- `pin_clock` -- P2-28: the pinned closure's idle clock [none (robin_form = pinned is retired)]
- `removal_increment` -- P2-28: the column removal's whole-cell increment and its prescribed-T leg [face_removal]
- `scalar_flaw` -- P2-28: the deterministic flaw on the column criterion [face_removal, face_buckle]
- `beam_void_closure` -- P2-28: the beam path through the column path's voided columns [beam, beam_tau, face_beam; a beam path across an interior void is a named coverage gap]
- `seed_void` -- P2-28: seeded interior voids for the column scans [face_list (CAV)]
- `surface_T_recon` -- P2-28: the column reconstruction (spall.surface_T = cell_average | diagnostic | reconstructed) [face_removal, face_heat (the face list's R1)]
- `sp_v_n_regime` -- P2-28: the column path's spall / vapour regime field [face_vapor]
- `sp_lefm` -- P2-28: site A's LEFM checks (prescribed stress and T hooks) [validation/kant/sp_kant_onset/test_face, sp_kant_pressure_sweep/test_face]
- `robin_face_ledger` -- P2-28: robin_face's (b) legs, the ledger on the column path (`input_ledger`) [face_feet (e), face_vapor, face_heat (v)]; robin_face's (a) form legs stay
- `gradient_blind_column` -- P2-28: the uniform slab on the column path (removal_energy = cell) [unit/gradient_blind, now the face-list deck]
