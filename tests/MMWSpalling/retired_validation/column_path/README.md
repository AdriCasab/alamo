# The column path's validation decks (retired in P2-28)

Site A -- the column path's no-removal onset criterion (`Sp_field`) -- is deleted in P2-28 (census H, decision 8 of 2026-11-05). These decks and tests ran it; an onset deck would now RUN silently as a thermal-only deck (its criterion field is gone, standing rule 1), so they are moved here as records, not run. Their last runs are P2-28 stage 1's (`studies/p2_28_column_path/K1_kant_siteA.log`, `K1_pressure_siteA.log`), on the stage-1 binary.

- `sp_kant_onset/` -- the five site-A onset decks (`input_verification`, `_hfl`, `input_model`, `input_weibull`, `_cluster`) and `test` (targets 1-6 and the 25-case sweep) [validation/kant/sp_kant_onset/test_face: targets 3, 5 on the face list's evaluate-only mode]
- `sp_kant_pressure_sweep/` -- the site-A base deck and `test` (R1-R4) [validation/kant/sp_kant_pressure_sweep/test_face]
- `sp_kant_onset_meiercfg/` -- the onset at Meier's configuration on the column path [validation/kant/grimsel_calibration (c6)]
- `p0_steady_column/` -- the steady-front reconstruction check on the column path [unit/face_removal --floor]
- `rossi/` -- the Rossi damage-profile deck (decision 3 of 2026-11-05: retire)
The decks write to their old `tests/MMWSpalling/validation/...` paths: records, not runnable as moved.
