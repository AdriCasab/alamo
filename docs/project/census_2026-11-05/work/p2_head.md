# Census part 2: MMWSpalling.H, function by function (2026-11-05)

Read-only map of `src/Integrator/MMWSpalling.H` as it stands on 2026-11-05: **9,375 lines,
6,292 code / 2,674 comment / 409 blank**. Re-anchored by function name against
`docs/project/CODE_CENSUS.md` (2026-10-08, 11,033 lines then). Nothing under `src/` or
`tests/` was edited. Every line number below was opened.

**How the numbers were made** (all in `work/`, re-runnable from the repo root):

| script | output | what it does |
|---|---|---|
| `split.py` | `split_raw.tsv` | brace-depth splitter (strings/comments stripped); every block opened at class-body depth, signature start, end, code/comment/blank, leading comment block. Boundaries spot-checked by reading (all 106 depth-2 blocks). |
| `callers.py` | `callers.tsv` | who calls each function, across `MMWSpalling.H`, `Removal.H`, `FaceRemoval.H` |
| `aborts.py` | `aborts.tsv` | all 185 `Util::Abort` sites with their guarding condition |
| `ranges.py` | `lineage_by_line.tsv` | assigns every one of the 9,375 lines a lineage and a CLOSED flag (function base + hand-checked sub-ranges); prints the budget of report 7 |
| `table.py`, `members.py` | `table.md`, `members.md` | the two tables below |

**What "production" means here.** The production configuration is now a file:
`tests/MMWSpalling/validation/shared/model_form_flame.inputs` + `rock_grimsel.inputs`
(both Meier and the Kant calibration `FILE =` them). It sets `surface_patch.face_list = 1`,
`robin_form = surface`, `spall.surface_T = reconstructed`, `removal_energy = flake`,
`stress_source = local_thermoelastic`, `stress_law = grimsel`, `kic.law = grimsel_table`,
`phase0.law = grimsel`, `weibull.enabled = 1`, `energy_ledger.enabled = 1`,
`spall.buckle_span = 0.012`, `surface.follow_mask = 1`, one phase in expression mode,
`melt_aware = 0`, `surface_normal = 0`, vapor/bit/damage off, and `max_level = 0`. The
face-list Parse block (`MMWSpalling.H:967-1112`) refuses most of the column path's options, so
"reached on production" below means "reached with that file".

**Lineage codes.** PROD = PROD-SHARED. The other nine are as in the brief. CLOSED is a separate
flag: a default-off key whose study closed per CODE_CENSUS C.7.3 (`surface_normal`,
`melt_aware`, `flake_length`/`flake_fraction`, `flake_coherence_length`, `pad_shield_open`,
`jet_above_h`/`jet_T_above`, `jet_T_ent_mode = rec_fixed`, `robin_form = face` as a control).
Two classification calls worth stating:
- `UpdateSpAfterMechanics` is COLUMN (the Sp_field criterion site) and only its
  `stress_mf` branch is MECH (`:6061-6062, :6125-6126, :6270-6276`).
- `face_heat = path` (`FaceListPath` and its helpers) is FLAME. It is a flame heat-map
  feature, and the P2-11b review recorded that the user OK'd deleting the path code.

---

## Table 1. Functions in file order

`lines` is the signature-to-closing-brace range. `+N` is the leading comment block attached
above it, which is counted in the lineage column. Lineage shows the main lineage, then the
other lineages it contains, by line count. Indented rows are lambdas or sub-blocks of the row
above.

