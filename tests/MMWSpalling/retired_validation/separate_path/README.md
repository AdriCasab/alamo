# The separate heat kernels' validation decks (retired in P2-29)

P2-29 (census D) deleted `AdvanceConstAlpha`, `AdvanceMaterial` and `AdvanceMaterialImplicit`: one heat path,
`AdvanceMicrostructure`. Each test below kept its name and took a production-kernel deck (one constant phase, or
`microstructure.phase0.law = zhang`); the OLD decks, their logs and their last outputs are the record here. A deck here
would now ABORT (`material.*` / `phase.*` / `implicit.*` retired by presence, `microstructure.nphases` required).

- `zhang_oglesby/` -- the global Zhang material + the implicit solver (2 s) [validation/zhang/zhang_oglesby/input:
  `phase0.law = zhang` on a QUARTER block, explicit 0.05 s; P2-29 a3/a4: the production kernel = the old explicit
  kernel = the full block to 2.3e-13 K]
- `hu_conduction/`, `hu_thermoelastic/` -- the material kernel with the HALF-conductance prescribed term and the
  homogeneous mechanics; `output/` is the 09-25 record (P2-29 m4/m5's baseline) [the same tests on the phase decks]
- `hu_breakage_index/` -- the AMR decks (max_level 1, no output since 09-17) [the same test on single-level 40^3 phase
  decks: AMR + mechanics fails on the production path, the user 2026-11-10]
- `hu_thermoelastic_amr/` -- RETIRED whole (decks + test; the user, 2026-11-10): AMR is parked (SR 12);
  `hu_thermoelastic` covers the thermoelastic physics. Re-instate with the AMR-repair packet.
