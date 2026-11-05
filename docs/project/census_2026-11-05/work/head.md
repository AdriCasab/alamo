# Census 2026-11-05, part 1: every input key the MMWSpalling integrator parses

Read-only census of the source as of the 2026-11-03 tree (`src/Integrator/MMWSpalling.H` 9,375
lines, `MMWSpalling/Removal.H` 2,580, `MMWSpalling/FaceRemoval.H` 937; mtimes Nov 3 08:29). No file
under `src/` or `tests/` was edited. Supersedes `docs/project/CODE_CENSUS.md` §C.7 (2026-10-08),
whose line numbers drifted and whose key set predates P2-10, P2-14a..P2-17.

Line references: `M` = `src/Integrator/MMWSpalling.H`, `R` = `src/Integrator/MMWSpalling/Removal.H`,
`FR` = `src/Integrator/MMWSpalling/FaceRemoval.H`; a bare `:NNNN` means `M:NNNN`.

## How this was produced (reproducible)

All under `docs/project/census_2026-11-05/work/`:

| script | does |
|---|---|
| `funcmap.py` → `funcmap.tsv` | maps every line of the three headers to its enclosing member function (heuristic; innermost wins) |
| `keys.py` | the key table, written by hand from every `pp.*` call: full dotted key, parse line(s), default, member(s) set, parse-branch guard |
| `readers.py` → `readers.json` | for each member, every read outside the nine parse functions (`Parse`, `ParseDamageSettings`, `ParseCohesiveSettings`, `ParseSpallSettings`, `ParseVaporSettings`, `ParseBitSettings`, `SurfacePatch::Parse`, `CheckPadRing`, `ParseSpallationModel`); comments and string literals stripped first |
| `setters.py` → `setters.json`, `files_scanned.txt` | who sets each key: (a) the two production legs with `FILE =` includes expanded and argv applied in AMReX order (last definition wins); (b) 84 files in `tests/MMWSpalling/unit/*/` (`input*`, `test`, `*.py`; `key=value` argv strings and `"key": value` dicts); (c) 38 `input*` files under the nine validation directories; (d) the 11 tree-gate tests' inputs and scripts |
| `annot.py` | the hand-made part: class, reader-path guard, production status, evidence |
| `render.py` → `tables.md` | the tables below |

The mechanical reader scan was then checked by hand against the code for every key; false
positives removed (e.g. `FrCell::mu`, `ElasticConsts::E`, the `ThermalFreeStrain(T_ref)` parameter,
`ScanTAtZv`'s local `alpha`) and indirect readers added (`beam.power_at`, `SurfacePatch::T_at`,
`ElasticConstants()`, the Grimsel law methods, the derived Drucker-Prager members).

### Scope of keys

- 197 keys parsed directly in `MMWSpalling.H` (top-level, `surface_patch.*`, `microstructure.phase<i>.*`),
  12 `beam.*` keys (`Numeric/MMWBeam.H:35-88`, `queryclass` at :138), 11
  `microstructure.phase<i>.grimsel.*` keys (`Numeric/Material/Grimsel.H:93-117`, `queryclass` at :372),
  22 retired keys (`contains()` only). One row per key; `phase<i>` rows stand for every phase index.
- Grouped rows (not counted as keys): `material.zhang.*` (20), `material.constant.*` (10),
  `damage.ic.{constant,expression}.*` (`select<>` sub-scopes).
- Out of scope, listed in the appendix: the inherited `hc.*` (HeatConduction) and `el.*` (Mechanics)
  base-class keys, plus `amr.*`, `geometry.*`, `timestep`, `plot_file`, `stop_time`, `alamo.program`.

### Two facts about ALAMO parsing that every row depends on

1. **Top-level keys are not strict.** `MMWSpalling::Parse` receives the root `IO::ParmParse`
   (`src/mmwspalling.cc`: `select_only<MMWSpalling>`), so an unread `spall.*`, `weibull.*`,
   `spallation.*`, `microstructure.*`, `damage.*`, `cohesive.*`, `vapor.*`, `bit.*`, `losses.*`,
   `phase.*`, `conductivity.*`, `surface.*`, `material.*`, `alpha` key is **silently ignored**
   (standing rule 1; comment `:89-96`).
2. **`queryclass` scopes are strict.** `IO::ParmParse::queryclass` (`src/IO/ParmParse.H:959-985`)
   aborts on any key under its prefix that the class's `Parse` did not read. That covers
   `surface_patch.*` (:190), `beam.*` (:138), `microstructure.phase<i>.grimsel.*` (:372), `hc.*`,
   `el.*`, and the `select<>` sub-scopes `material.<type>.*`, `damage.ic.<type>.*`. Consequence:
   deleting a parse inside a strict scope turns every deck that still sets the key into an abort
   (decks must be edited in the same step); deleting a top-level parse makes such decks silently
   ignore the key.

### The production path these rows are judged against

Two legs (brief): `studies/p2_16_span_rule/legs/HOS_R/argv.sh` (Höser lab deck
`validation/hoeser/input_lab` + 41 `key=value` argv tokens) and `legs/CAL/argv.sh` (Kant calibration deck
`validation/kant/grimsel_calibration/input`), both including `validation/shared/rock_grimsel.inputs`
and `validation/shared/model_form_flame.inputs`. Effective configuration: microstructure path, one
phase (`mode = expression`, `law = grimsel`), `el.type = disable`, beam power 0,
`surface_patch.face_list = 1`, `robin_form = surface`, `nozzle_descent = prescribed`,
`jet_closure = enthalpy` (HOS_R) / `none` with `T_flame_expr` (CAL), `spallation.model = sp_weibull`,
`stress_source = local_thermoelastic`, `stress_law = grimsel`, `kic.law = grimsel_table`,
`surface_T = reconstructed`, `removal_energy = flake`, Weibull per-face flaws,
`spall.buckle_span` 0.012 (HOS_R) / 0.0 (CAL: rock file's 0.012 overridden later in the Kant deck),
`amr.max_level = 0`. Functions that execute on it each step: `BuildSurfaceTColumns`,
`UpdatePropertyFields`, `UpdateSpAfterMechanics`, `UpdateSpClusters`, `UpdateRemovalAfterCohesive`
**only up to its face-list dispatch at `R:543-546`** then `RecedeAlongLines` (FR), and
`AdvanceMicrostructure` (`BuildFaceList`, `BuildFlameColumns`, `FaceListHeat`→`FaceH`/`FaceHr`,
`FrAssignTs`, `FaceListTg`, `JetEnthalpyMarch` [HOS_R], `FrFaceFlux`, `EnergyLedger`,
`PatchDiagnostics`, `FaceListDiagnostics`), `InvertHtoT`. **Everything in `Removal.H` after line 546
is the column path and never runs on production.** `UpdateModel` never runs (`el.type = disable`,
`src/Integrator/Base/Mechanics.H:167`).

## Classification rule

Each key gets exactly one class. Precedence: RETIRED-ABORT > DEAD > CLOSED > feature class.

| class | meaning used here |
|---|---|
| PROD | defines the production model: set on a production leg (explicitly or via the shared files) **and** read on that leg's path, or REQUIRED by the parser there. A PROD key whose value is inert on production says so in column (a). |
| FLAME-CAP | flame capability off the production legs but needed to run some flame (feet, pinned form, constant-coefficient flame, D2c free surface, X1 path rule, jet bin width) |
| BEAM | beam, vapour, melt, latent-heat and drawworks-bit keys (column path today; the face list refuses vapour/bit/melt) |
| ZHANG | global-material and const-α thermal paths: `material.type`, `solver`, `implicit.*`, `phase.*`, Hu's prescribed-T `T_expr` |
| MICRO | Voronoi / multi-grain / grain-boundary machinery |
| MECH | the mechanics stress source and criterion alternatives that are not on production: eigenstrain `alpha`, α–β quartz strain, the linear K_Ic law, prescribed σ/T for criterion site A, the cluster-detach rule |
| DAMAGE | Drucker-Prager damage law, cohesive zone, damage-law spall sampling (`spall.sample`, `prescribed_ratio`, `C_h`/`alpha`/`t_spall`, `xi_s`, `damage_threshold`), damage-modified conductivity |
| CLOSED | a study switch the record closed; evidence named in the row |
| DIAG | CSV / thermo / debug switches |
| DEAD | parsed, no reader (or reader unreachable) |
| RETIRED-ABORT | in P2-10's or P2-16's retired-key list (`:97-132`) |

## Counts

| class | keys | grouped rows |
|---|---|---|
| PROD | 80 | |
| FLAME-CAP | 17 | |
| BEAM | 32 | |
| ZHANG | 12 | +2 (`material.zhang.*`, `material.constant.*`) |
| MICRO | 5 | |
| MECH | 14 | |
| DAMAGE | 35 | +1 (`damage.ic.*` sub-keys) |
| CLOSED | 12 | |
| DIAG | 7 | |
| DEAD | 6 | |
| RETIRED-ABORT | 22 | |
| **total** | **242** | 3 |

Keys set by **no** production leg, unit test, validation deck or gate input (13):
`phase.T_v`, `phase.lh_v`, `implicit.tol_abs`, `microstructure.phase<i>.alpha_beta.linear_strain`,
`microstructure.phase<i>.grimsel.T_ref`, `surface.refine_depth`, `surface.refine_height`,
`cohesive.G_c`, `cohesive.f_t`, `spall.h_spall_n_segments`, `vapor.prescribed_rop.enabled`,
`vapor.prescribed_rop.expr`, `weibull.V_ref`.

## Column guide

- **Default**: the `query_default` fallback, or REQUIRED (with its condition), or ABORT.
- **Readers**: reads of the member(s) outside the parse functions, `n: Function×count (file:line,…)`.
  Parse-time cross-key checks are not counted here; they are in report (iii).
- **Guard**: the parse branch the key sits in (left of ‖), then the condition under which the
  readers execute (right of ‖).
- **(a) Production**: effective value on each leg, then whether the production path reads it and
  whether the value matters there.
- **(b) Kept unit tests**: number of `unit/` directories whose `input*`, `test` or `*.py` sets the key,
  then `dir(values)`. Values inside test scripts include deliberate abort legs (e.g. `foot_npads=0`,
  `foot_pad_quantile=1.5`, the `face_heat` R_* legs, `lateral_flux_lock`'s `size_volume=flux`) and
  f-string placeholders (`{...}`).
- **(c) Validation decks**: files that set the key **directly** (a deck that only includes the shared
  rock/model-form files is not counted), grouped by directory, `dir×files(values)`.
- **(d) Tree gate**: the 11 tests run by `studies/tree_sweep_0922/witness/post_all.sh` + `wit.sh`
  (`unit/{jet_d2c, robin_face, robin_pinned, robin_jet, jet_enthalpy, robin_feet, beam_void_closure,
  scalar_flaw, spall_event}`, `validation/meier/sp_meier_pilot/input_2d_dev`, `studies/s1_surface_resolution
  --smoke`), listing every gate input/script that sets a **non-default** value (for REQUIRED keys,
  every setter). `jet_d2c` runs `robin_feet/input_feet`, `jet_enthalpy/input_{strip,quarter}`,
  `robin_jet/input_jet` with argv; those inputs are listed under their own directories.

## Tables
