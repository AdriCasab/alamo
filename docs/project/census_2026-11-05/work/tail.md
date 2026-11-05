
## Report (i): keys parsed more than once

A `contains()` test followed by a query is listed when the two sit at different sites; a plain
`contains()`-then-`query` pair on adjacent lines is the standard "given?" idiom and is only counted.

| key | sites | what to know |
|---|---|---|
| `hc.heat.alpha` | `src/Integrator/HeatConduction.H:79` (via `queryclass("hc")` at :136) **and** `:153` | the only key parsed by two classes. Both REQUIRED. The second copy feeds `therm_alpha`, read at `AdvanceConstAlpha:1476` and `R:297` (the h_spall positivity abort that runs on every spall step, production included) |
| `material.type` | `:179` `contains()`, then `select_default` re-queries it (`src/IO/ParmParse.H:1073`) | without it the whole `material.*` block is skipped silently (report ii) |
| `damage.ic.type` | `:495` `contains()`, then `select` `query_required`s it (`ParmParse.H:1045`) | same pattern |
| `microstructure.nphases` | `:292` `contains()` (path selector), `:332` `query_required` | its presence gates ~45 keys (report ii) |
| `damage.phi_deg`, `damage.cohesion` | `:5021`/`:5022` REQUIRED under `damage.enabled`; `:5026`/`:5027` defaulted otherwise | mutually exclusive branches, different defaults (REQUIRED vs 0.0 / 1.0) |
| `surface_patch.T_expr` | `:8027` REQUIRED under `mode = prescribed_T`; `:8034` default `""` otherwise | under `convective_flame` it is consumed and ignored |
| `surface_patch.h_expr`, `T_flame_expr` | `:8049`/`:8050` `contains()`; `:8263`/`:8269` query | `has_hx`/`has_Tx` drive five cross-checks before the query (:8062-8068, :8243, :8253-8261) |
| `surface_patch.h_conv`, `T_flame` | `:8259` `contains()` (abort with the expressions); `:8342`/`:8343` query | |
| `surface_patch.nozzle_z0`, `nozzle_feed` | `:8276` `contains()` (abort under feet); `:8318`/`:8319` query | |
| `surface_patch.foot_body_clearance` | `:8248` `contains()` (abort without feet); `:8308` query | |
| `surface_patch.jet_De_ref` | `:8139` query (momentum); `:8144` `contains()` abort otherwise; `:8229` loop | |
| `surface_patch.jet_fs_aspect`, `jet_fs_exponent` | `:8156`/`:8157` query; `:8163` `contains()` abort; `:8229` loop | |
| `surface_patch.{jet_decay_diameter, jet_De_ref, jet_free_surface, jet_fs_aspect, jet_fs_exponent, jet_bin_update, jet_negative_flux, jet_above_h, jet_T_above}` | own query inside the enthalpy branch, plus `:8226-8232` `contains()` loop in the `else` branch | the else-loop is the strict-scope abort, restated with a clearer message |
| `contains()`+query "given?" pairs | `grain_expr` :321/:323; `cohesive.delta_c` :5289/:5291; `cohesive.delta_max` :5290/:5292; `spall.alpha` :5687/:5690; `spall.t_spall` :5692/:5695; `vapor.A_beam` :5820/:5823; `vapor.T0` :5828/:5831; `surface_patch.face_heat` :7912/:7918; `barrier_expr` :7924/:7930; `pad_shield_open` :8004/:8005; `jet_core_length` :8077/:8078; `jet_T_rec_fixed` :8110/:8113; `jet_T_above` :8187/:8188 | each `contains()` is there because `query_default` ADDS the key (comment :8076) |
| retired keys (22) | `contains()` only, `:119`, `:128` | never parsed |

Member-name collisions (not double parses, but they make greps lie): `losses.h_conv` →
`MMWSpalling::h_conv` and `surface_patch.h_conv` → `SurfacePatch::h_conv`; `T_ref` is the member
of `phase.T_ref` (class), `microstructure.phase<i>.T_ref` (`MineralPhase`), `...grimsel.T_ref`
(`Material` base) and the parameter of `ThermalFreeStrain`; three "alpha" keys (`alpha`,
`hc.heat.alpha`, `spall.alpha`) set three different members.

## Report (ii): keys parsed inside a keyed branch (silently inert on decks that do not take it)

Only **top-level** keys can be silently inert. Inside a strict `queryclass` scope an unread key
aborts, so the `surface_patch.*` branches are safe **except** for the consumed-but-ignored cases in
the last block.

**A. Behind `contains("microstructure.nphases")` (`:292-624`)** — a deck without `nphases` silently
ignores: `microstructure.mode` :299, `number_of_grains` :306, `seed` :307, `phase_expr` :312,
`grain_expr` :321-323, every `microstructure.phase<i>.*` key :338-411 (25 keys, including the
otherwise-strict `grimsel.*` sub-scope :372, which is never entered), `conductivity.kappa_damage_alpha`
:489, `conductivity.h_gb0` :491, `damage.ic.*` :495-498. Nested:
- `mode = voronoi` (:304): `phase_expr`, `grain_expr` inert. `mode = expression`: `number_of_grains`, `seed` inert.
- phase index: `microstructure.phase<k>.*` for k ≥ `nphases` inert (loop bound :335).
- `damage.ic.type` selects one IC sub-scope; the other IC type's keys are inert.

**B. Behind `contains("material.type")` (`:179`)** — every `material.zhang.*` / `material.constant.*`
key without `material.type`, and the non-selected type's keys with it.

**C. `ParseDamageSettings`** — `damage.prescribed_stress.{xx,yy,zz,xy,yz,zx}_expr` behind
`damage.prescribed_stress.enabled` (:5038). (`phi_deg`, `cohesion`, `psi0`..`stiffness_floor` are read
on both branches.)

**D. `ParseCohesiveSettings`** — every `cohesive.*` key behind `cohesive.enabled` (:5278);
`cohesive.G_c`, `cohesive.f_t` additionally behind "a delta is missing" (:5293): **a deck giving both
deltas silently ignores G_c and f_t**; `prescribed_delta.expr` behind its enable (:5316).

**E. `ParseSpallSettings` after `if (!value.spall_enabled) return;` (`:5681`)** —
`spall.damage_threshold`, `xi_s`, `C_h`, `alpha`, `t_spall`, `prescribed_ratio.enabled`, `.expr`,
`sample`, `shallow_depth`, **`h_spall_max`**, `h_spall_n_segments`, `flake_coherence_length`,
`seed_void`. ⚠ **A latent trap with a reader outside the branch:** `spall.surface_T`, `flake_c_f` and
`flake_h_max` are deliberately parsed BEFORE the return (onset legs run the reconstruction with
removal off), but the reconstruction `BuildSurfaceTColumns:2979,3019` also reads `spall_h_spall_max`,
which is parsed only AFTER it. An onset deck with `spall.enabled = 0`, `surface_T = reconstructed` and
its own `spall.h_spall_max` would reconstruct with the member default 0.010 m. No current deck does
this (every reconstructed deck under `unit/` and `validation/` sets `spall.enabled = 1`). The same
ordering breaks the parse check `:5610-5612` (`flake_h_max < h_spall_max`), which always compares
against 0.010 (see Surprises).

**F. `ParseVaporSettings` after `:5818`** — `vapor.A_beam`, `vapor.T0`, `vapor.prescribed_rop.*`.

**G. `ParseBitSettings` after `:5859`** — every `bit.*` key.

**H. `ParseSpallationModel` after `if (model != SpWeibull) return;` (`:9164`)** — a `damage_law` deck
silently ignores `spallation.sp.a0`, `n_segments`, `stress_relief`, `kic.*`, `prescribed_stress.*`,
`prescribed_T_surface.*`, **and every `weibull.*` key including `weibull.m`** (moved outside the
`weibull.enabled` guard by P1g-3, :9277-9288, but still inside this one). Nested:
- `kic.law` (:9221): `kic.scale` inert under `linear`; `kic.K_Ic0`, `kic.T_ref`, `kic.T_melt` inert under
  `nasseri_table` / `grimsel_table` — **a production deck that sets `kic.K_Ic0` is silently ignored**.
- `weibull.enabled` (:9292): `a0_gb`, `a0_ig`, `V0`, `size_volume`, `V_ref`, `seed`, `detach_mode`, `A_crit`.
- `detach_mode = connected_cluster` (:9361): `A_crit`.
- `prescribed_stress.enabled` (:9245) / `prescribed_T_surface.enabled` (:9257): their expressions.
- Deliberately NOT behind it: `confining.p` :9151 and `spall.buckle_span` :9160 (P2-15 comment :9155-9158).

**I. Strict scopes, consumed but ignored (silent despite strictness)** —
`surface_patch.T_expr` under `convective_flame` (:8034, read and discarded); `surface_patch.type`
under `convective_flame` (:7940, checked only under `prescribed_T`, :7952); `surface_patch.face`
(:8014, only `zhi` legal, never read); `surface_patch.T_ext` (:8041, never read);
`surface_patch.pinned_idle_cycles` under `robin_form = surface` (:7993, logged by the Message at
:770-774 only when ≠ 2.0); `microstructure.phase<i>.grimsel.T_ref` (Grimsel.H:115, overwritten at
:7060). And `surface_patch.enabled = 0` turns every other `surface_patch.*` key into an abort
(`hoeser/input_lab` alone therefore does not run: it sets `enabled = 0` while including
`model_form_flame.inputs`, which sets `robin_form` and `face_list`; HOS_R's argv sets `enabled = 1`).

## Report (iii): the three interlock sites (CODE_CENSUS B14), what each refuses today

B14's runtime half has shrunk: the P1g-0b assertions and the P1f parse-time `SurfaceCellFlux`
self-test are gone with P2-10 (no match for `P1g-0b` / `self-test` / `sentinel` assertions in
`src/Integrator/MMWSpalling*`).

### Site 1 — `SurfacePatch::Parse` (`:7897-8356`) and `CheckPadRing` (`:7887-7895`, called on both exits :8338, :8355)

Value checks aside, it refuses:
- `face_heat` without `face_list = 1` (:7914-7916); `barrier_expr` without `face_heat = path` (:7926-7928)
- `mode = prescribed_T` with `type` ≠ `prescribed_temperature` (:7952-7955)
- `robin_form = face` without `convective_flame` (:7969-7971); `pinned`/`surface` without `convective_flame` (:7976-7981)
- `jet_closure = enthalpy` without `h_expr` (:8063-8065) or with `T_flame_expr` (:8066-8068)
- `jet_T_ent_mode = exhaust|rec_fixed` without `jet_entrained_mass = 1` (:8104-8106); `rec_fixed` without `jet_T_rec_fixed` (:8110-8112) or with it below `jet_T_ent` (:8114-8116)
- `jet_decay_diameter = momentum` without `jet_stagnation = decay` (:8133-8135) or without an explicit `jet_core_length` (:8136-8138); `jet_De_ref` without momentum (:8144-8146)
- `jet_fs_aspect`/`jet_fs_exponent` without `jet_free_surface = 1` (:8163-8165)
- `jet_above_h = 1` without `jet_T_above` (:8189-8192); `jet_T_above` without `jet_above_h` (:8193-8195); `jet_T_above ≤ jet_T_ent` (:8196-8199)
- any of 9 enthalpy-only jet keys under `jet_closure = none` (:8226-8232)
- `nozzle_descent = feet` without `h_expr` (:8242-8244); `jet_free_surface = 1` without feet (:8245-8247); `foot_body_clearance` without feet (:8248-8250)
- `h_expr`/`T_flame_expr` without `convective_flame` (:8253-8255), half of the pair (:8256-8258), or with `h_conv`/`T_flame` (:8259-8261)
- `nozzle_z0`/`nozzle_feed` under feet (:8276-8278); `foot_body_clearance = 1` with `foot_rule` ≠ `pads` (:8312-8314)
- constant flame (`convective_flame`, no expressions) without `h_conv > 0` / `T_flame > 0` (:8344-8351)
- `pad_shield_open ≥ 0` without a pad ring `foot_r_outer > foot_r_inner` (`CheckPadRing`; so the shield needs a feet run)

### Site 2 — `MMWSpalling::Parse` cross-key block (`:626-1124`) and the sub-parsers' cross-key aborts

- microstructure required by `damage.enabled` :626-627, `cohesive.enabled` :629-630, `spall.enabled` :632-633, `vapor.enabled` :646-647, `sp_weibull` :649-653; `kic.law = grimsel_table` needs a phase :640-645
- `surface.follow_mask = 1` needs removal :665-667 and refuses `prescribed_T` :668-671; follow_mask/ledger need microstructure :673-676; ledger single-level :677-679
- `amr.max_level > 0` with removal needs `surface.refine_depth > 0` :690-693 and `amr.nsubsteps = 1` :694-699
- `robin_form = face`: microstructure only :701-705, monotone emissivity :710-716
- `robin_form = pinned`: microstructure :722-724, follow_mask :725-727, spall + sp_weibull :728-731, monotone emissivity :732-736
- `robin_form = surface`: microstructure :748-750, follow_mask :751-753, spall + sp_weibull :754-757, monotone emissivity :758-761, `surface_T = reconstructed` :762-767; `pinned_idle_cycles` ≠ 2.0 only **messaged** :770-774
- `removal_events_csv = 1` needs spall + sp_weibull :780-786
- `surface_T`/`removal_energy` non-default need sp_weibull :790-799; `removal_energy = flake` needs spall :801-803
- `h_expr` needs `follow_mask` :839-843
- `foot_body_clearance = 1` needs removal :894-897
- `spall.buckle_span > 0` needs `face_list = 1` :949-953 and `surface_T = reconstructed` :954-956
- `face_heat = path` refuses `pad_shield_open ≥ 0`, `jet_free_surface`, feet, `jet_above_h`, `robin_form` ≠ surface, `jet_closure` ≠ enthalpy :981-992
- `face_list = 1` requires `convective_flame` + `h_expr` :994-997, `robin_form = surface` :998-1000, microstructure :1001-1004, `follow_mask` :1005-1008; and refuses (`fr_abort` :1027-1061): vapour, bit, `foot_body_clearance` (regime 4, owner "P2-11"), `melt_aware`, `surface_normal = 1`, `stress_source ≠ local_thermoelastic`, `amr.max_level > 0`, a removal model other than sp_weibull spallation, `surface_T ≠ reconstructed`, `removal_energy ≠ flake`, `flake_length = flux` / `flake_fraction = weibull` / `size_volume = flux`, `flake_coherence_length > 0`, `connected_cluster` with `A_crit > 0`
- default `vapor.A_beam` non-positive :1120-1122
- in the sub-parsers: `damage.enabled` needs active mechanics unless prescribed stress `:5066-5070`; `flake_h_max < h_spall_max` `:5610-5619` (⚠ compares against the default, see Surprises); `removal_energy = flake` needs `surface_T = reconstructed` `:5651-5658`; reconstruction/flake need `local_thermoelastic` `:5670-5678`; `bit.enabled` needs spall or vapour `:5861-5863`; sp_weibull refuses `cohesive.enabled` `:9169-9171` and `damage.enabled` `:9172-9175`
- in `Initialize`: `material` + `microstructure` together `:1135-1137`

### Site 3 — runtime aborts that check configuration

- `BuildFaceList:3705-3706`: `face_list = 1` needs cubic cells
- `UpdateSpAfterMechanics:6023-6025`: `local_thermoelastic` needs phase 0 with E > 0 and mu > 0 (this is what keeps `phase0.E`/`mu` load-bearing under the law)
- `R:299-304`: h_spall needs `spall.alpha`/`hc.heat.alpha` > 0 and `t_spall` > 0 — on every spall step, face list included, though h_spall is unused there
- `R:688-690`: sp_weibull spall needs phase 0 (column path); `R:733-741`: sp_weibull spall without a T-derived stress needs an active mechanics solve (column path)
- surface-owning level only: `FeetDescent:2304-2305`, `JetEnthalpyMarch:2472-2473`, `R:439-440` (foot clearance)
- geometry/value at run time: `FeetDescent:2333`, `:2369` (no live column under the feet); `BuildFlameColumns:2143` (nozzle collision), `:2160`, `:2198` (h/T finite); `FaceHr:3683`; `JetEnthalpyMarch:2753` (free-surface column in the axis bin); `R:2365`, `:2375` (removal ran out of rock / voided > 1 cell, column path)

## Surprises (for the deletion packet)

1. **Six DEAD keys, four of them on the production or validation decks.** `microstructure.phase<i>.alpha_attenuation`
   (REQUIRED, set 2000 on both legs, zero readers), `microstructure.phase<i>.name` (set on both legs),
   `surface_patch.face` (set on both legs; only `zhi` is legal), `surface_patch.T_ext` (27 Hu decks + 4 units),
   `surface_patch.type` (17 Hu decks), `microstructure.phase<i>.grimsel.T_ref` (set by nothing; its only consumer
   receives an overwritten copy, :7060). The three `surface_patch.*` ones are in a strict scope: removing the parse
   without removing the deck lines turns those decks into aborts.
2. **`spall.flake_h_max`'s parse check reads `spall.h_spall_max` before it is parsed** (`:5610` vs `:5752`, and
   `:5752` is behind the `spall.enabled` return), so it always compares with the member default 0.010 m. Production
   (0.0015 < 0.015) is unaffected; a deck with `h_spall_max = 0.008`, `flake_h_max = 0.009` would pass, and one with
   `h_spall_max = 0.020`, `flake_h_max = 0.012` would abort wrongly. Same parse-order trap as the memory note.
3. **`hc.heat.alpha` is not inert on production** (the rock/Kant decks call it inert): `R:297-301` aborts every spall
   step if it is ≤ 0 and `spall.alpha` is unset, before the face-list dispatch. The damage-law h_spall it guards is
   computed and discarded every production step (`spall.C_h`, `spall.alpha`, `spall.t_spall` also read there).
4. **`spall.removal_energy` has no runtime reader on the face list**: all four readers are column path (`R:985-2382`);
   on production it is only the face-list interlock `:1050-1051`. Likewise `spall.sample`, `spall.flake_coherence_length`,
   `spall.melt_aware`, `spall.surface_normal` (column half) never run on production although the model-form file sets them.
5. **`phase0.E`, `mu`, `beta`, `kappa`, `rho`, `Cp` are REQUIRED but value-inert under the Grimsel law**; only
   `E > 0 && mu > 0` matters (the abort at :6023). Two Poisson ratios survive in the input language (`grimsel.nu`
   used; `phase0.E/(2 mu) − 1` used only without the law, `R:43`).
6. **The tree gate pins closed keys**: `spall.surface_normal = 1` and `spall.flake_coherence_length = 0.004` are set by
   `unit/beam_void_closure/input` and `meier/input_2d_dev`; `jet_T_rec_fixed` (closed D2e diagnostic) by `unit/jet_d2c`;
   the pinned form (`pinned_idle_cycles`) by four gate tests; `spall.melt_aware = 1` by five gate inputs. Deleting any of
   them moves the 2639-file tree and needs a re-pin.
7. **13 keys are set by nothing at all** (list under Counts), including `spall.h_spall_n_segments`, which is read on
   every production step at its default (0 = defer to `spallation.sp.n_segments`), and `surface.refine_depth`/`height`,
   the whole AMR surface band (only `studies/amr_band/run.py` ever set them; the face list refuses AMR).
8. **`surface_patch.side_power_interval`** (a CSV writer, set on HOS_R) is parsed inside the `jet_closure = enthalpy`
   branch (:8206), so a `jet_closure = none` deck such as CAL cannot set it without a strict-scope abort.
9. **`weibull.m` and `weibull.a0_gb`**: `m` was moved out of the `weibull.enabled` guard but is still behind
   `sp_weibull` (:9164); `a0_gb` is REQUIRED (:9294) on production yet inert on the grain-free rock.
10. **The retired-key lists are only load-bearing for `spall.*`**; the 16 `surface_patch.*` retired keys would abort as
    "specified but not used" anyway (strict scope). `unit/retired_keys` and `unit/lateral_flux_lock` assert the
    lists' source text, so editing them needs those tests edited too.

## Appendix: inherited base-class keys (not censused row by row)

- `hc.*` — `Integrator::HeatConduction::Parse` (`src/Integrator/HeatConduction.H:79-98`): `heat.alpha` (REQUIRED;
  re-read at :153), `heat.refinement_threshold`, `ic.*` (`select_default<IC::Constant, IC::Sphere, IC::Expression>`),
  `bc.temp.*` (`select_default<BC::Constant, BC::Expression>`), `method`. Strict scope. Production sets `heat.alpha`,
  `heat.refinement_threshold = 1.0`, `ic.type = constant`, `ic.constant.value`, `bc.temp.constant.type.*` and `.val.*`.
- `el.*` — `Integrator::Mechanics::Parse` (`src/Integrator/Mechanics.H:124-157`) and `Base::Mechanics::Parse`
  (`src/Integrator/Base/Mechanics.H:56-59`, returns early on `type = disable`): `type`, `nmodels`, `model<n>.*`,
  `eta_ref_threshold`, `ref_threshold`, `model_neumann_boundary`, `ic.*`, `eta.reset_on_regrid`, `psi.ic.*`.
  Production: `el.type = disable`, `nmodels = 1`, `model1.{E, nu, F0}`.
- Integrator/AMR keys (`amr.*`, `geometry.*`, `timestep`, `stop_time`, `plot_file`, `alamo.program`) belong to
  `Integrator::Integrator` and are out of scope.
