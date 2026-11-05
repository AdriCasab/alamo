
---

## Report 1. The four conduction kernels (B6)

`Advance` dispatches at `:1457-1460`: if `microstructure_enabled`, `AdvanceMicrostructure`;
else if `material` and `solver = implicit`, `AdvanceMaterialImplicit`; else if `material`,
`AdvanceMaterial`; else `AdvanceConstAlpha`. A run never mixes kernels. `Initialize:1136`
aborts on `material` together with `microstructure`.

| | `AdvanceConstAlpha` 1472-1511 | `AdvanceMaterial` 1513-1633 | `AdvanceMaterialImplicit` 1642-1771 | `AdvanceMicrostructure` 1788-2085 |
|---|---|---|---|---|
| runs when | no `material`, no `microstructure` | `material.type`, `solver = explicit` | `material.type`, `solver = implicit` | `microstructure.nphases` (**production**) |
| conductivity | constant `hc.heat.alpha` x `Numeric::Laplacian` (`:1495`) | `mat->kappa(T)` per neighbour, harmonic face (`:1582-1600`) | `mat->kappa(T_old)` harmonic B-coefficients (`:1715-1727`) | per-cell `kappa_mf` (Grimsel-refreshed by `UpdatePropertyFields`) x `exp(-a_D D)`, harmonic, GB contact resistance in series (`k_face` `:1930-1966`) |
| enthalpy units / inversion | H in kelvin (`H = T - T_ref`, `ConvertTtoH:6680-6690`). Adds `Q` [W/m^3] and `q_loss/dz` to `a_th lap T` [K/s] with no 1/(rho c_p) (`:1508`), so it is consistent only if the input takes rho c_p = 1 | global `material_table` | BE on T, then `H = material_table.H_of_T(T_new)` (`:1762-1770`) | per-phase `phase_tables` (Grimsel `Table::Build`, `:7051-7066`) |
| voids | none | void-aware (`removed` ghosts, `:1557-1593`) | ignores voids | void-aware |
| surface patch | none | prescribed_T or constant h on the **domain** top `k = bigEnd` (`:1609-1617`) | none (Neumann hard-coded on every face, `:1735-1741`) | `SurfaceCellFlux` at the mask column top (follow_mask) in four forms, plus face-list `FrFaceFlux` (`:2003-2033`) |
| losses | inline on `surface_mf` (`:1504-1507`) | inline on `surface_mf` (`:1624-1628`) | inline, explicit (`:1700-1704`) | via `SurfaceCellFlux` (cell form `:4035-4044` for out-of-patch tops) |
| beam path | `BeamSource(..., DX[2])`, geometric path z0 - z | same, geometric: **solid cells below a cavity are attenuated through the void** (the E1 leak survives on this path) | same, geometric | `SolidFacePath` solid count (`:2048-2061`) |
| AMR / level | any | any | asserts `finest_level == 0` (`:1645`) | owner-gated surface terms (`:1825`) |

**Shared copy-paste.** The beam-local block appears in all four kernels: `:1476-1483`,
`:1518-1525`, `:1650-1657` and `:1842-1849`. The loss-local block appears four times as well.
The inline loss expression appears three times: `:1506`, `:1627` and `:1702`. The fourth
kernel uses `SurfaceCellFlux` instead, and `SurfaceCellFlux` repeats the expression in its own
three forms (`:3979, :4027, :4041`). `kc = kappa exp(-a_D D)` appears 8 times: `:1923`,
`:1943`, `:2500`, `:3152`, `:3757`, `:4149`, `:5125` and `FaceRemoval.H:454`. The
`a_D == 0` exp-skip exists only at `:1923/:1943/:3757`. The march (`:2500`), PatchDiagnostics
(`:3152`) and the ledger (`:4149`) call `exp` unconditionally.

**Production** runs `AdvanceMicrostructure` and nothing else. `solver = implicit` combined with
a microstructure is **silently ignored**: the dispatch checks `microstructure_enabled` first,
and nothing aborts on it (see report 4). `AdvanceConstAlpha`, `AdvanceMaterial` and
`AdvanceMaterialImplicit` are reached only without a microstructure block, which means the
Step 2-4 regressions (const-alpha) and the Zhang validation (`material.type`, explicit or
implicit). Hu runs on the microstructure kernel through its `prescribed_T` branch (`:1996`).

**After a ZHANG merge onto the phase-law path** (Zhang as `microstructure.phase0.law = zhang`,
a single-phase expression microstructure, its `Material` going through the same polymorphic
`Table::Build` the Grimsel law uses at `:7051-7066`, and `UpdatePropertyFields` generalised
from "grimsel" to "any T-dependent law"), these become redundant:
- `AdvanceMaterial` (121 lines).
- `AdvanceConstAlpha` (40). The Step 2-4 regressions move to a constant one-phase
  microstructure. Note that their H is in kelvin today, so the references re-pin.
- `BuildEnthalpyTable` and `SetPhaseFractionsConstAlpha`.
- The global-material and closed-form branches of `ConvertTtoH` (`:6676-6690`) and
  `InvertHtoT` (`:6724-6757`).
- The `material` + `microstructure` abort (`:1136`).
- Members `material`, `material_table`, `SolverKind`, `implicit_*`, and `T_ref/T_m/T_v/lh_*`
  (the const-alpha `phase.*` keys, `:155-161`).
- Parse `:174-181` and `:252-267`.
- `UpdateModel`'s homogeneous `alpha[]/eta` branch (`:1361-1398`, MECH, reachable only
  without a microstructure).

`AdvanceMaterialImplicit` has no microstructure equivalent: it is single-level, Neumann-only and
void-blind. Either it is ported (an implicit microstructure kernel) or the Zhang implicit leg is
re-run explicit. That is a decision, not a mechanical merge.

---

## Report 2. The surface-temperature machinery (B3)

Every place a column-top or face surface temperature is formed, and whether production
(face list) reaches it:

| site | lines | what | production (face list) |
|---|---|---|---|
| `RawSurfaceX` | 9000-9010 | two-cell x = ln((T_top-T0)/(T_1-T0)) | **yes**: called by `FrAssignTs` (`FaceRemoval.H:235`) and `RecedeAlongLines` (`FaceRemoval.H:501`) per surface cell along its line |
| `ReconFromX` | 9020-9109 | T_surf, ell_eff, h_f (c_f, `h_spall_max` cap, `flake_h_max` bound), T_bar | **yes**: same two callers (`FaceRemoval.H:240, :504`) |
| `SurfaceGasTempCap` | 9117-9124 | ceiling max(T_flame, jet_T_nozzle) | **yes** (`FaceRemoval.H:229, :469`) |
| `FrAssignTs` (FaceRemoval.H) | 162-253 | per listed cell T_s (carry from removal, or fresh recon, or fallback flag) | **yes**, called from `BuildFlameColumns:2232` |
| `FrFaceFlux` (FaceRemoval.H) | 138-156 | every listed face's flux at its cell's T_s; fallback = face Newton on e* | **yes**: kernel `:2030`, march `:2624/:2666`, ledger `:4187`, `FaceListDiagnostics:3919` |
| `SurfaceCellFlux` surface form | 3974-3986 | flux + losses at a supplied T_surf | **yes**, through `FrFaceFlux` (`T_surf_in = Ts`) |
| `SurfaceCellFlux` face form (Newton) | 3987-4022 | infer T_f on a linear half-cell | **yes**, only as `FrFaceFlux`'s fallback (`FaceRemoval.H:148`) and `RecedeAlongLines:523` |
| `SurfaceCellFlux` pinned branch | 4023-4033 | T_s = T_pin | no (face list requires `robin_form = surface`, `:999`) |
| `SurfaceCellFlux` cell form | 4035-4044 | flux on the cell average + losses on `surf_c` | **yes, for out-of-patch column tops only** (in-patch tops are all listed and their result is zeroed at `:2020`) |
| `FaceH` / `FaceHr` / `FaceListHeat` | 3663-3691, 3848-3901 | per-face h (the heat rule; no temperature) | **yes** |
| `BuildSurfaceTColumns` | 2866-3033 | per-column recon into `surf_T_col` etc. via `GatherColumnSurface` + `ColumnSolidRunBelowTop` | **runs** (`:1434`, keyed on `surface_T`), but under the face list its output only feeds diagnostics: the Sp_field site, the ledger's `surf_q_mean` (`:4285-4292`), and `SurfaceRobinT`, whose consumers discard it |
| `SharedSurfaceRecon` | 3046-3062 | rebuild a SurfaceRecon from the column arrays | Sp_field site (`:6159`) and the column removal path (`Removal.H:1316`, after the dispatch). Diagnostic only on production |
| `SurfaceRobinT` | 3069-3074 | column T_surf for `robin_form = surface` | called by kernel `:2004-2009` (result zeroed for listed cells, `:2020`), march `q_col :2607-2610` (discarded when `fl_on_m`, `:2617-2629`), ledger `:4167-4174` (zeroed `:4179`), `PatchDiagnostics:3163-3166` (live, but a column-form diagnostic) |
| `BuildPinEffective` / `PinRule` | 3076-3093, 3098-3115 | the pinned T_pin | no (refused by the face list) |
| ell_flux q_net (flake_length = flux) | 2974-3014 | flux-based heated depth | no (`fr_abort :1052`); CLOSED |
| Sp_field `T_recon_at` | 6181-6187 | lerp of T_surf to the top-cell average over dz/2 | runs (plotted only) |
| `Removal.H` column scan `ScanTAtZv` / `T_at_didx` | Removal.H 68-103, ~1390-1425 | the scan profile | no (after the `face_list` return at `Removal.H:543-547`) |

**Column-top sweeps per production step** (B3 said "about eight"):
- `GatherColumnSurface` three times: `BuildSurfaceTColumns:2883`, `UpdateSpAfterMechanics:6073`
  and `UpdateSpClusters:6378`. All three feed diagnostics only.
- `ColumnSolidRunBelowTop` once (`:2884`).
- `Removal.H`'s inline `k_top_global` once (`Removal.H:512-537`, before the dispatch). Production
  needs it.
- `ColumnTopSolid` in the kernel once (`:1828-1830`). Production needs it.

That is **5 top sweeps + 1 run sweep, each with an MPI reduce**. Two are needed (before and
after removal). The three `GatherColumnSurface` calls and the run sweep exist only for the
COLUMN diagnostics.

---

## Report 3. The criterion sites (B4)

| site | where | reached on production? |
|---|---|---|
| A: Sp_field | `UpdateSpAfterMechanics` `:5992-6343`. Builds `sigma_at_depth` (`:6219-6276`: prescribed / law / scalar / `stress_mf`), `K_I` at `:6296`, `K_Ic` at the tip `T_recon_at(a*icth)` `:6303-6316`, writes `Sp_field` `:6319`. Called at `Advance:1446`. | **Runs every step** (guard: `sp_weibull && owner`). On the face list its output is **plotted and hashed only**: `Sp_field`, `T_surf_field`, `ell_eff_field` and `recon_*`. `FaceRemoval.H` never reads `Sp_field_mf`. The only other readers are `UpdateSpClusters` (`:6386`, diagnostic) and the column path's detach gate (`Removal.H:1100-1104`). It does **not** apply the P2-15 buckling condition or `melt_aware`. Without removal (spall.enabled = 0) it is the only criterion, but the face list requires `spall.enabled`, and the Kant calibration now runs the face list (`kant/grimsel_calibration/input` `FILE =`s `model_form_flame.inputs`). |
| B: column depth scan | `Removal.H:1225-1570` (`UpdateRemovalAfterCohesive`) through the shared pieces `ScanTAtZv` (`Removal.H:68`), `ScanKIcPa` (108), `ScanSigmaT` (132) and `DepthScanH` (158) | **No**. `Removal.H:543-547` returns into `RecedeAlongLines` when `face_list = 1` |
| C: line rule | `FaceRemoval.H:331-935` (`RecedeAlongLines`). Criterion at `FaceRemoval.H:560-591` (same `ScanSigmaT`/`DepthScanH`, with `f_skel = 1.0` and `sn_cth = sn_icth = 1.0` passed as literals). Buckling withdrawal at `:597-630` (`BuckleRatio :315`) | **Yes, this is the production criterion** |

Sites B and C share one copy of the arithmetic, the P2-9 pieces at the top of `Removal.H`. Site
A remains a second hand-written copy: `sigma_at_depth` `:6219-6276` repeats `ScanSigmaT`'s two
branches verbatim, and its `K_Ic` switch `:6309-6316` repeats `ScanKIcPa`. The "two-firing-site
invariant" comments at `:6052-6056` and `:6200-6213` now describe A against B, and B is not the
production path. A port or deletion of site A should either call `ScanSigmaT`/`ScanKIcPa` or go
entirely.

---

## Report 4. Feature interlocks (B14): every abort that refuses a key combination

Value-range checks (finite, > 0, enum spelling) are left out. These are the 185 `Util::Abort`s
in `work/aborts.tsv` filtered to combinations.

**(a) `MMWSpalling::Parse` (`:88-1125`)**

| line | refused combination |
|---|---|
| 120 | any of 21 P2-10 retired keys present (`:98-118`: lateral_*, per_cell_removal, side_face_*, step_diag_radius, jet_wall_*, jet_body*) |
| 130 | `spallation.sp.constraint_depth` present (P2-16) |
| 627 / 630 / 633 / 647 / 651 | damage / cohesive / spall / vapor / sp_weibull without `microstructure.nphases` |
| 641 | `kic.law = grimsel_table` with no phase |
| 666 | `surface.follow_mask` without spall or vapor |
| 670 | `follow_mask` with `surface_patch.mode = prescribed_T` |
| 675 | `follow_mask` or ledger without microstructure |
| 678 | ledger with `amr.max_level > 0` |
| 692 / 697 | AMR + removal without `surface.refine_depth > 0` / with `nsubsteps != 1` |
| 704 / 713 | `robin_form = face` without microstructure / with `epsilon_high < epsilon` |
| 723 / 726 / 730 / 733 | `robin_form = pinned` without microstructure / follow_mask / (spall && sp_weibull) / with `eps_hi < eps` |
| 749 / 752 / 756 / 759 / 763 | `robin_form = surface`: the same four, plus without `spall.surface_T = reconstructed` |
| 784 | `removal_events_csv` without (spall && sp_weibull) |
| 792 | `surface_T` != cell_average, or `removal_energy = flake`, with `spallation.model != sp_weibull` |
| 802 | `removal_energy = flake` without `spall.enabled` |
| 842 | `h_expr` without `follow_mask` |
| 896 | `foot_body_clearance` without removal |
| 950 / 955 | `spall.buckle_span > 0` without `face_list` / without `surface_T = reconstructed` |
| 986-991 | `face_heat = path` with `pad_shield_open >= 0`, `jet_free_surface`, `nozzle_descent = feet`, `jet_above_h`, `robin_form != surface`, or `jet_closure != enthalpy` |
| 995 / 999 / 1002 / 1006 | `face_list` without (convective_flame && h_expr) / robin_form surface / microstructure / follow_mask |
| 1033-1061 (`fr_abort`, defined 1028) | `face_list` with vapor, bit, feet + foot_body_clearance, `melt_aware`, `surface_normal`, `stress_source != local_thermoelastic`, `max_level > 0`, a model other than sp_weibull spall, `surface_T != reconstructed`, `removal_energy != flake`, `flake_length = flux` / `flake_fraction = weibull` / `size_volume = flux`, `flake_coherence_length > 0`, or connected_cluster with `A_crit > 0` |

**(b) `SurfacePatch::Parse` / `CheckPadRing` (`:7887-8356`)**

| line | refused combination |
|---|---|
| 7890 (`CheckPadRing`, called on both exits `:8337`, `:8355`) | `pad_shield_open >= 0` without `foot_r_outer > foot_r_inner` |
| 7915 / 7927 | `face_heat` without `face_list` / `barrier_expr` without `face_heat = path` |
| 7954 | `mode = prescribed_T` with legacy `type != prescribed_temperature` |
| 7970 / 7979 | `robin_form = face` / `pinned` / `surface` without `mode = convective_flame` |
| 8064 / 8067 | `jet_closure = enthalpy` without `h_expr` / with `T_flame_expr` |
| 8105 / 8111 | `jet_T_ent_mode` exhaust or rec_fixed without `jet_entrained_mass` / rec_fixed without `jet_T_rec_fixed` |
| 8134 / 8137 / 8145 | `jet_decay_diameter = momentum` without `jet_stagnation = decay` / without an explicit `jet_core_length`; `jet_De_ref` without momentum |
| 8164 | `jet_fs_aspect`/`exponent` without `jet_free_surface` |
| 8190 / 8194 / 8198 | `jet_above_h` without `jet_T_above` / the converse / `jet_T_above <= jet_T_ent` |
| 8230 | any of 9 D2c/D2v jet keys without `jet_closure = enthalpy` |
| 8243 / 8246 / 8249 | feet without h_expr / `jet_free_surface` without feet / `foot_body_clearance` without feet |
| 8254 / 8257 / 8260 | `h_expr` without convective_flame / without (T_flame_expr or enthalpy) / with `h_conv` or `T_flame` |
| 8277 / 8313 | feet with `nozzle_z0`/`nozzle_feed` / `foot_body_clearance` with `foot_rule != pads` |

**(c) Elsewhere: the `Parse*Settings` helpers and run time**

| line | refused combination |
|---|---|
| 1136 (`Initialize`) | `material` together with `microstructure` (run time, not Parse) |
| 2305 / 2473 (FeetDescent, JetEnthalpyMarch) | feet / jet enthalpy on a non-owner level (run time) |
| 3706 (`BuildFaceList`) | face_list with non-cubic cells (run time, every step) |
| 5069 (`ParseDamageSettings`) | `damage.enabled` with `el.type = disable` and no prescribed stress |
| 5612 (`ParseSpallSettings`) | `flake_h_max > 0 && flake_h_max >= h_spall_max`. **Parse-order defect**: `spall.h_spall_max` is not read until `:5752`, after the early return at `:5681`, so this compares against the member default 0.010 (`:8570`), not the deck's value. A deck with `h_spall_max = 0.015, flake_h_max = 0.012` aborts wrongly, and `h_spall_max = 0.001, flake_h_max = 0.005` passes. Production (0.0015 vs 0.015) passes either way. |
| 5652 / 5673 | `removal_energy = flake` without `surface_T = reconstructed` / (`surface_T` recon or flake) without `local_thermoelastic` |
| 5862 (`ParseBitSettings`) | `bit.enabled` without spall or vapor |
| 9170 / 9173 / 9333 (`ParseSpallationModel`) | sp_weibull with cohesive / with damage; `size_volume = flux` without `V_ref` |
| 6024 / 6031 (`UpdateSpAfterMechanics`, run time) | `local_thermoelastic` with phase 0 lacking E, mu, or physical nu |
| `Removal.H:735` (run time) | sp_weibull spall, neither prescribed_T nor local_thermoelastic, with mechanics disabled |

**Combinations nothing refuses** (silently inert or silently ignored on production):
- `solver = implicit` with a microstructure is ignored (dispatch `:1457`).
- `spallation.sp.prescribed_stress.*` and `prescribed_T_surface.*` only change `Sp_field`
  (`:6113-6126`); the face list never reads them.
- `spall.prescribed_ratio`, `spall.sample != top_cell`, `spall.h_col_events_csv` and
  `spall.clusters_csv` act on the column path only. The deck comment
  `model_form_flame.inputs:6` says so for `sample`.
- `pinned_idle_cycles` under `surface` is only a `Util::Message` (`:770-774`).
- A beam (`beam.P0 > 0`) with the face list is not refused, but the beam deposits through
  column tops only (report 5).

---

## Report 5. BeamSource and the beam energy balance

- `BeamSource` is at `:4908-4947`: its 20-line header sits in the gap at `:4888-4907` and its
  `AMREX_GPU_HOST_DEVICE` at `:4907`. It has a 10-line lead comment and a 30-line body
  (`:4918-4947`). `BeamIncidentIntensity` is `:4949-4980`, `SolidFacePath` `:4315-4331`,
  `ColumnSolidAboveCount` `:4613-4665`, and `BeamEnergyBalance` `:4716-4887` (32-line header
  + 140 body).
- The beam is parsed by `Numeric::MMWBeam` (`MMWBeam.H:72-77`: `beam.alpha_expr` as a parser
  of T only, default "10"). The power schedule is `MMWBeam::power_at` (`MMWBeam.H:93-100`).
- The beam's lineage in this file: the beam counters (Parse `:140-149`, members `:7583-7598`),
  the kernel block (`:1833-1841`, `:2035-2061`) and the vapour parse.

**Is alpha(T) energy-conserving? No, only for constant alpha.** `BeamSource` uses
`a = alpha_f(T_cell)` (`:4932`) both for the cell's own absorption and for the attenuation over
the whole material path above it, `path = SolidFacePath(n_above)` (`:4933-4943`). The analytic
per-cell integral `(e^{-a s_top} - e^{-a s_bot})/dz` telescopes to `1 - e^{-a s_N}` down a
column only if every cell has the same `a`. With alpha(T) each cell attenuates with its own
`a_k` instead of the cumulative optical depth `sum_j a_j dz`, so the column can over- or
under-deposit. The code says so itself (`BeamEnergyBalance` header `:4723-4727`, "needs the
cumulative optical depth ... required before real alpha(T) MMW work"). The closure check that
would catch it, `P_dep <= P_inc (1 + 1e-9)`, is applied to the **uniform** profile only
(`:4859-4861`). The Gaussian profile carries a further O((dz/z_R)^2) error, because `Pi` uses
the geometric depth z0 - z (`:4966-4972`).
- On the ZHANG kernels the path is geometric (`dz_abs = -1`), so void cavities attenuate the
  beam (report 1).
- On the face list the beam enters only through column tops (`ktc`, `sab`). Listed side, up
  and down faces receive nothing, and `FrFaceFlux` has no beam term. A beam port must decide
  this.

**Per-phase `alpha_attenuation`.** It is parsed with `query_required` at `:347-348` into
`MineralPhase::alpha_atten` (`:8375`), so every microstructure deck must carry it.
`grep -rn alpha_atten src/` finds only those two lines, so it is **never read**. The production
rock file says the same (`rock_grimsel.inputs`: "beam only; never read (MODEL_CARD §12 item
11)"). The beam's attenuation comes from `beam.alpha_expr` alone, with no phase dependence.

---

## Report 6. Everything DEAD (CODE_CENSUS C.7.4 updated)

**Computed, then discarded.**

On the production step:
- `AdvanceMicrostructure`: `SurfaceCellFlux` is evaluated for every in-patch top cell
  (`:2003-2009`), then its `q_in/q_out` are zeroed for every listed cell (`:2020`). Every
  in-patch top cell is listed, because `BuildFaceList` uses the same radius test.
- `JetEnthalpyMarch`, under the face list:
  - the `jTc`/`jkc` gather with its two `ReduceRealSum` (`:2485-2504`) feeds only the
    `SurfaceCellFlux` + `PinRule` call in `q_col` (`:2603-2611`);
  - that call's result is then replaced by the listed faces' sum (`:2617-2629`).
- `EnergyLedger`: `SurfaceCellFlux` per listed cell is zeroed at `:4179`.
- `RefreshKappaEff` called by `FaceRemoval.H:929` (and `Removal.H:2576`) writes `kappa_eff_mf`.
  The kernel then overwrites every cell of it in the same step at `:2072`. Only the
  init/regrid calls matter.
- The Sp_field site, `UpdateSpClusters`, `BuildSurfaceTColumns` and three `GatherColumnSurface`
  sweeps run every step, and their only consumers are plotted fields and thermo columns
  (report 3). They are diagnostic, not dead, but nothing in the production trajectory reads
  them.
- `PatchDiagnostics` (`:3122-3196`) runs under the face list (`:2083-2084`, `sp_on &&
  use_expr`) and reports `patch_P_robin` from the column-form `SurfaceCellFlux`. That is a
  flux the face-list kernel does not apply.

Locals:
- `pc` (`:1926`, `(void)` `:1928`), `p_n` (`:1945-1947`), `q_rob`/`q_los` (`:2063`) in the
  kernel.
- The `phi` patch bound and unread at `:1916-1917` (kernel), `:4779-4780`
  (`BeamEnergyBalance`) and `:5920` (`GatherColumnSurface`).
- `face_top` at `:3913`/`:3920` (`FaceListDiagnostics`).
- `z_lo` (`:6855-6856`, `InitPhi`), `DX` (`:5175`, `UpdateDamageAfterMechanics`).

Parameters:
- `dt` of `BuildSurfaceTColumns` (`(void) dt` `:2950`) and `z_n` of `FaceListPath`
  (`(void) z_n` `:3518`).
- `dt` of `UpdateSpAfterMechanics` / `UpdateSpClusters` / `UpdateCohesiveAfterMechanics`.
- The redundant `(void)a_time; (void)a_step;` in `TimeStepComplete` (`:1192-1193`), which
  uses both above.

Members written or parsed and never read:
- `SurfacePatch::T_ext` (`:7822`, parsed `:8041`).
- `fl_wtop` (`:3290`, written `:3808`, `:3816`).
- `sp_surface_T_mode` (`:8965`, written `:5555-5567`; only the two bools are read).
- `MineralPhase::alpha_atten` (report 5).
- The legacy `surface_patch.type` key (`:7938-7940`), checked only for equality under
  prescribed_T (`:7953-7955`).

Unreachable or inert:
- Invariant (b) of `BeamEnergyBalance` exists only as comment (`:4740-4742`, `:4822-4828`).
- The k_eff bulk diagnostic counts phases rank-locally with no MPI reduction
  (`:7368-7378`, "fine for Step 6 (no MPI tests)"). It is inert for one phase but
  decomposition-dependent for several. It is a defect, not dead code.

**Stranded or stale text.**
- `:622` "P1b: the acceptance test for the support rule is that" (unfinished; still there).
- `:1401-1404`: the Advance header ("no material -> AdvanceConstAlpha") sits above
  `SurfaceRecon` and omits the microstructure branch.
- `:2846-2853`: the C1 pin-rule header sits above `BuildSurfaceTColumns`; it belongs to
  `BuildPinEffective` at `:3076`.
- `:3035-3039`: `SurfaceRobinT`'s header sits above `SharedSurfaceRecon`.
- `:3064-3068`: the orphan "P1 -- THE REFERENCE TEMPERATURE" fragment (still there).
- `:3198-3233`: `SurfaceCellFlux`'s 36-line header sits **694 lines** above the function
  (`:3928`); `:3234` is the `RiserWeight` tombstone.
- `:3939-3955`: the "P1f ... the lateral criterion needs the face temperature" rationale. The
  lateral criterion was deleted in P2-10; `T_face_out` is now read by `FrFaceFlux`.
- `:4479-4484`: `ColumnTopSolid`'s header sits above `UpdatePropertyFields`, 105 lines before
  the body (`:4589`). The census said 500 lines; the distance shrank but the comment is still
  stranded.
- `:5588-5598`: the "P1 -- THE GRADIENT-AWARE DRIVING STRESS, `absolute` (default)" block
  describes a key with no parse (deleted at both sites, `:6209`).
- `:5899-5905`: `GatherColumnSurface`'s header says "phi >= 0". The rule has been
  removed-only since D2q-1 (`:5928`).
- `:5942-5949`: `UpdateSpAfterMechanics`'s header sits above `ColumnSurfaceCos`.
- `:6574-6578`: `LawSigma` says the Sp_field site sits in "`amrex::ParallelFor` under
  `[=, this]`". It is `amrex::LoopOnCpu` with `[&]` (`:6136`), and the comment points at the
  deleted lateral criterion.
- `:8360-8361`: the duplicated "Step 5: Voronoi mineral microstructure" line.
- `:8764-8774` and `:8796`: P1/P1b support-rule text with no members (`spall.support_rule`
  is parsed nowhere).
- Stale line references:
  - `:637-638` (":206", ":293")
  - `:787`
  - `:2888` (`kc_cell`)
  - `:2967` (`Removal.H:1503-1517`)
  - `:2991` (":7920-7925")
  - `:4500` (":1436")
  - `:5478-5481` (`Removal.H:691`, ":4668")
  - `:5669`
  - `:6584` (`Removal.H:1833`)
  - `:7215` (":5380")
  - `:7703` (`Removal.H:589-594`)
  - `:8932` (`Removal.H:846-858`)

**C.7.4, item by item**

| C.7.4 item | now |
|---|---|
| kernel `pc`, `p_n`, `q_rob`, `q_los` discarded | **still** (`:1926-1928`, `:1945-1947`, `:2063`) |
| `BuildSurfaceTColumns` `dt` unused | **still** (`:2950`) |
| Removal.H PASS 3 temperature array | Removal.H; not re-checked here (part 3) |
| `SurfacePatch::Parse` legacy `type` | **still** (`:7938-7955`) |
| `BeamEnergyBalance` invariant (b) comment only | **still** (`:4740-4742`, `:4822-4828`) |
| Removal.H P1b support-rule comment blocks | Removal.H; grep still finds P1b text at `Removal.H:640, 1448, 1611, 2058, 2213, 2321, 2490` |
| `SurfaceRobinT` orphaned comment | **still** (`:3064-3068`) |
| `ColumnTopSolid` comment 500 lines above | **still stranded**, now 105 lines above (`:4479-4484` vs `:4589`) |
| `SurfaceCellFlux` `AMREX_GPU_HOST_DEVICE` attached to nothing | **gone** (`SurfaceCellFlux :3928` has no macro); its header comment is now stranded instead |
| Parse "P1b: the acceptance test ..." | **still** (`:622`) |
| `BuildFlameColumns` pad-ring block pasted twice | **still** (`:2166-2184` and `:2205-2222`) plus a third code copy in `FaceHr :3686-3689` |
| `UpdateSpAfterMechanics` "READ 1 of 2" pasted twice | **gone** (once, `:6278`) |

---

## Report 7. Line budget

From `work/ranges.py`. Every line is counted once. CLOSED lines are counted inside their
lineage and also flagged.

| lineage | lines | code | comment | blank | of which CLOSED |
|---|---|---|---|---|---|
| PROD-SHARED | 2,969 | 1,827 | 1,003 | 139 | 0 |
| FLAME | 1,686 | 1,264 | 396 | 26 | 170 |
| BEAM | 526 | 307 | 200 | 19 | 45 |
| ZHANG | 414 | 313 | 72 | 29 | 0 |
| MICRO | 543 | 349 | 164 | 30 | 0 |
| MECH | 354 | 299 | 35 | 20 | 0 |
| DAMAGE | 640 | 523 | 63 | 54 | 0 |
| COLUMN | 1,634 | 992 | 573 | 69 | 223 |
| DIAG | 564 | 418 | 126 | 20 | 0 |
| DEAD (stranded comment) | 45 | 0 | 42 | 3 | 0 |
| **total** | **9,375** | 6,292 | 2,674 | 409 | 438 |

**Removal scenario** (DAMAGE + COLUMN-only + CLOSED switches + ZHANG + stranded text):

| removed | lines |
|---|---|
| DAMAGE | 640 |
| COLUMN (incl. its 223 CLOSED) | 1,634 |
| ZHANG | 414 |
| CLOSED outside COLUMN: FLAME 170 (pad shield, above-plane pair, rec_fixed, robin_form = face aborts) + BEAM 45 (`melt_aware`) | 215 |
| DEAD stranded text | 45 |
| **total removed** | **2,948** |
| **remaining** | **6,427** (6,472 if the stranded text is kept) |

Adjustments to that 6,427:
- The ZHANG merge also makes `UpdateModel`'s homogeneous branch redundant (`:1361-1398`, about
  39 MECH lines) and adds roughly 40-80 lines: a generic T-dependent phase law in
  `UpdatePropertyFields` and the `law = zhang` parse. The net is about **6,430 ± 20** after
  the merge.
- Further FLAME cuts that are column-only on production and are not in the brief's list:
  - the pinned closure (refused by the face list): about 135 lines (`:717-740`,
    `:828-838`, `:1891-1897`, `BuildPinEffective`, `PinRule`, `:4023-4033`,
    `:7614-7632`, `:7976-7998`);
  - `face_heat = path`, which the user OK'd deleting: about 385 lines (`FaceListPath`,
    `FrPathOffsets`, `FaceCentre`, `:3303-3309`, `:3855-3863`, `:975-992`, `:1096-1111`,
    `:7810-7815`, `:7911-7934`, the march's path branch `:2630-2672`);
  - `PatchDiagnostics` (80).
  - Together these bring the file to about **5,830**.
- Deleting COLUMN also deletes plotted fields: `Sp_field`, `Sp_cluster_id`, `T_surf_field`,
  `ell_eff_field`, `h_spall_field`, `crack_event_count` and the `recon_*` thermo columns. That
  re-pins the gate header. It also requires editing the production callers
  `FaceRemoval.H:928` (`UpdateSurfaceMaskFromPhi`) and `Advance:1434/:1446/:1450`.
- MECH (354) stays by the brief's rules but is unreachable on production (`el.type = disable`).

---

## Report 8. The two closed-study branches at the criterion/removal sites a beam port will edit

The face list refuses both (`fr_abort`, `:1037` and `:1038`). The production criterion
therefore never takes either branch, but both still thread through the shared criterion
pieces the face list calls, as literal 1.0 arguments. Every site follows.

**`spall.melt_aware` (f_skel / MeltSkeleton)**
- `MMWSpalling.H:5502-5526`: parse of `spall.melt_aware` and `melt.lambda_crit`, with the
  range abort at `:5525`.
- `MMWSpalling.H:8599-8605`: members `melt_aware` and `melt_lambda_crit`, and their comment.
- `MMWSpalling.H:8606-8618`: `MeltSkeleton` and its lead comment.
- `MMWSpalling.H:1037`: the face-list refusal.
- `Removal.H:128-146`: `ScanSigmaT`'s **`f_skel` parameter**, which multiplies both branches
  (`:138`, `:143`). The face list calls it with the literal `1.0` at `FaceRemoval.H:577`. A
  port deletes the parameter and that argument; multiplying by exactly 1.0 is an identity, so
  the deletion is bit-safe.
- `Removal.H:1017-1057`: the `LL_col_prof` gather and its `ReduceRealSum`, and `LL_col_p`.
- `Removal.H:1366-1386`: `LLcol_i` and the `LL_at_didx` lambda.
- `Removal.H:1464-1479`: `f_skel = MeltSkeleton(LL_at_didx(...))`.
- `Removal.H:1501-1502`: `f_skel` in the `stress_mf` branch.
- `Removal.H:1527`: `f_skel` passed to `ScanSigmaT`.
- Site A (`UpdateSpAfterMechanics`) **never applied** the skeleton factor, so the two sites
  disagree under the key.
- Inputs: 9 set it per C.7.3, and `model_form_flame.inputs` sets `spall.melt_aware = 0`. By the
  P2-10 pattern the key moves to the retired-key list at `:98-118`.

**`spall.surface_normal` (the criterion half; the recession half is already superseded)**
- `MMWSpalling.H:5448-5501`: parse plus the `Util::Message` at `:5493-5500`.
- `MMWSpalling.H:8593-8598`: member `spall_surface_normal`.
- `MMWSpalling.H:8619-8622`: `SP_SN_COS_MIN`.
- `MMWSpalling.H:8551-8557`: `sn_recession_superseded`.
- `MMWSpalling.H:938-940`: its thermo registration.
- `MMWSpalling.H:1038`: the face-list refusal.
- `MMWSpalling.H:5942-5990`: `ColumnSurfaceCos`, including the stranded Sp header.
- Site A:
  - `:6081-6090` (`sn_cos_col`, `sn_cos_p`)
  - `:6189-6194` (`icth`)
  - `:6220` (`dv = depth * icth`)
  - `:6300-6304` (`T_kic` at `a_loc * icth`)
- `Removal.H:548-558`: `sn_cos_col`/`sn_cos_p`, after the face-list return.
- `Removal.H:158-216`: `DepthScanH`'s **`sn_cth`/`sn_icth` parameters**, used at `:174-175`,
  `:179`, `:201-204`, `:209-211` and `:215`. The face list passes `1.0, 1.0` at
  `FaceRemoval.H:591`. A port deletes both parameters and the two literals; this is bit-safe.
- `Removal.H:1224-1231`: per-column `sn_cth`/`sn_icth`.
- `Removal.H:1460-1463`: `z_v = z_abs * sn_icth`.
- `Removal.H:1544`: the `DepthScanH` call.
- `Removal.H:1558-1569`: `h_vert = h_col * sn_icth`, the superseded recession.
- `Removal.H:1807-1810`: `sn_active` switches off the coherence cap.
- `Removal.H:2000-2009`: the comment.
- `Removal.H:2403-2424`: `sn_recession_superseded`.
- Inputs: 4 set it per C.7.3, two of them inside the 2,639-file sweep (`unit/beam_void_closure`
  and the Meier `input_2d_dev` leg). Removing it re-pins those.

---

## Appendix. What from CODE_CENSUS C.2 no longer exists

| C.2 entry | status now |
|---|---|
| C.2.4 `UpdateColumnPhase`, `AssertStepSplit`, `WallStepDiagnostic` | **deleted** (0 grep hits) |
| C.2.5 the P2-1 wall law: `WJ_*`, `WallMu`, `WallZFill`, `WallJetBranch`, `WallJetH`, `SideFaceH`, `ColumnWallCorner`, `PatchAzimuthFraction`, `UpdateWallLawState` | **deleted** (0 hits; `side_face_h_law` is a retired-key abort, `:107`) |
| C.2.5 the D2z helpers `ColPhi`, `PairDelta`, `StepCellWeight`, `StaircaseFactor`, `AreaOfDelta`, `StepConductance` | **deleted** (0 hits) |
| C.2.3 `RiserWeight` | **deleted** (tombstone comment `:3234`) |
| C.2.9 `PerCellRemoval` | **deleted** (`spall.per_cell_removal` retired, `:102`) |
| C.2.1 `ParseSpallSettings` P1f parse-time `SurfaceCellFlux` self-test (B14) | **deleted** |
| C.2.2 kernel `col_in_patch` side-face block; C.2.4 march `SideHit` gather; C.2.10 the ledger's side-face/area-law/D2x/D2z/P1e statistics | **deleted**; the face list (`FrFaceFlux`) replaces all four side-face copies (B2 is closed) |
| C.2.1 `SurfacePatch` side-face block, `side_face_closure` (`SideTForm`) | **deleted** (tombstone `:7807-7808`) |
| P2-2e constraint depth (`spallation.sp.constraint_depth`) | **deleted** in P2-16 (retired-key abort `:130`) |
| carried state `side_T_eff`, `lat_accum_col`, `col_t_last`/`col_dt`/`col_nvoid`, `wall_corner`, `side_cols_above`, `wall_step_*` | **gone** (0 hits). `T0_col_prof` survives only in `Removal.H` |
| C.2.3 claim "SurfaceCellFlux: 20 positional params, 8 call sites" | now 18 positional + 4 defaulted (`:3928-3937`), **7** call sites (kernel, march, PatchDiagnostics, ledger, `FrFaceFlux` x2, `RecedeAlongLines`) |
| sizes | `Parse` 1,100 → 1,039; `AdvanceMicrostructure` 583 → 298; `JetEnthalpyMarch` 569 → 376; `EnergyLedger` 699 → 250; `ParseSpallSettings` 508 → 412; `BuildFlameColumns` 170 → 137 (body); file 11,033 → 9,375 |
| new since C.2 (none of these are in CODE_CENSUS) | face list: `FaceEntry`, `BuildFaceList`, `FaceListHeat`, `FaceListTg`, `FaceListDiagnostics`, `FaceH`/`FaceHr`, `FaceWeight`, `FrGauss`/`FrSmoothAt`/`FrSmoothedGradient`, `FrInitFraction`, `FillFractionGhosts`; path: `FrPathOffsets`, `FaceListPath`, `FaceCentre`; P2-15 buckle members + Parse block `:943-965`; retired-key guards `:90-130` |
