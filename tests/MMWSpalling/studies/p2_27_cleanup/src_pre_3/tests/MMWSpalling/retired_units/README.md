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
