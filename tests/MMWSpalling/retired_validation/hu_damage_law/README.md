# Hu validation on the damage law -- retired at P2-27

The user, 2026-11-09 (census Q5): Hu's spallation onset is re-expressed on the LEFM rule
(`validation/hu/hu_spall_onset/input_granite2_lefm` + `test_lefm`: the Grimsel criterion under Hu's granite-2 ramp, a
record), and everything that ran Hu on the Drucker-Prager damage law retires with the lineage (P2-27 stage 2, census
C1a; decision 1 of 2026-11-05: LEFM is the model). Records, not run; their decks name keys that now abort and write to
their old paths.

- `hu_spall_onset/` -- the damage-law onset decks `input_granite2`, `input_sandstone2`, their `test`, the `diag/` runs,
  the AMR / bus-error analyses, `plot_mesh_fields`, the last `output/` and `sweep_run.log`
- `hu_end_to_end/` -- the damage-law end-to-end removal decks (granite 2, sandstone 2), "real removal, never run"
- `hu_spall_onset_mmwbeam/` -- extension: the damage-law Hu onset under an MMW beam (the LEFM sibling
  `extensions/hu_spall_onset_mmwbeam_lefm` stays)
- `hu_spall_onset_voronoi/` -- extension: the damage-law Hu onset on a Voronoi microstructure
