# Census 2026-11-05, part 3: the two removal paths

Read-only census. Every line number was taken from the working tree on 2026-11-05:
`src/Integrator/MMWSpalling/FaceRemoval.H` (937 lines), `src/Integrator/MMWSpalling/Removal.H`
(2,580 lines, 1,008 of them comment-only), and `src/Integrator/MMWSpalling.H` (9,375 lines).
In this document `F:` means FaceRemoval.H, `R:` means Removal.H and `M:` means MMWSpalling.H.
The helper scripts are in `work/`:
- `gate_keys.sh` reads the removal keys of the tree-gate inputs.
- `test_paths.sh` tags each test by the column-path features it uses.
- `fn_end.py` finds function line ranges.

**The plan this document serves:**
1. Port the beam onto the face list.
2. Rebuild melt and vapour on the enthalpy path.
3. Port Meier's feet.
4. Delete `Removal.H`.

**How the step reaches removal.** `Advance` (M:1424–1470) runs, in this order:
1. `BuildSurfaceTColumns` (M:1434).
2. Mechanics, then `UpdateSpAfterMechanics` (M:1445), the Sp_field site, criterion site A.
3. `UpdateSpClusters` (M:1449).
4. `UpdateRemovalAfterCohesive` (M:1454), which is the whole of `Removal.H`.
5. The conduction advance.

Under `surface_patch.face_list = 1` the function runs R:270–537 and then hands over at R:543–547:

```
if (surface_patch.face_list) { RecedeAlongLines(lev, time, dt, k_top_global); return; }
```

So the face list still depends on the prologue of the column function (Part B, stages B0–B6).

---

## A. The face list (FaceRemoval.H)

### A.1 Rule and stage map

| # | Stage | F: lines | What it does |
|---|---|---|---|
| — | header R0–R3 | 9–36 | the rule statement |
| R0 | types and state | 38–93 | `FrCell` 39–47, `FrLine` 49–58, `FrCarry` 59, `FrClass`/`FrNdotE` 61–71. Witness members 74–85, `FrKey` 87–93 |
| R0 | `FrCells()` | 95–131 | one surface cell per run of `fl_faces` sharing (col, k). e* is the exposed face with the largest n·e; mu = n·e*. A cell with n·e* ≤ 1e-6 is a misaligned **sheet** (124); a sheet has mu = 1 |
| F4 | `FrFaceFlux` | 133–156 | the **heat-BC** consumer, not removal. Every listed face is heated at its cell's T_s. A fallback cell (`tsfb`) gets T_s from the face Newton on `fl_faces[qline]`. Callers: kernel M:2030, march M:2624/2666, `FaceListDiagnostics` M:3919, ledger M:4187 |
| R1 (BC side) | `FrAssignTs` | 158–253 | runs at list rebuild, from `BuildFlameColumns` M:2232. A carried cell keeps the T_s handed over by removal (189–197). A fresh cell gathers T(c), T(c−e*) and the solidity of c−e* (203–227), then calls `RawSurfaceX` → `ReconFromX` (234–249). x = 0 gives Ts = 0 and `tsfb = 1` |
| R3 | `RecedeAlongLine` | 255–302 | consumes h/dx of solid fraction from the face inward (279–283). The heat goes through `RemapProfile` over W = min(Dmax, run−1) (286–288), then the index shift by n_v (294–300). E_left is telescoped (289–293) |
| B | `BuckleRatio` | 304–326 | bk = C_L·x/L_span, with C_L = π/√(6(1+ν)ε_th). ε_th comes from `Grimsel::strain_of` under the law, else β(T_s − T0_s). It returns +inf for a non-finite x or for ε_th ≤ 0 |
| — | `RecedeAlongLines` | 328–935 | the dispatch target, and the whole removal step under the key |
| | per-step resets | 339–357 | fr_* and bk_* witnesses; `FrInitFraction` (M:3585–3597); fr_score and fr_mu set to 0 |
| R0 | lines | 367–398 | `hcap` = min(flake_h_max, h_spall_max) (368–369). Dmax = floor(h_spall_max/(mu·dx)) (377–378). Lg = Dmax + 2 + ceil(hcap/(mu·dx)) (379–380). A sheet gets one line per face with mu = 1 (388–395) |
| — | **the one gather** | 400–465 | NF = 5 columns per line cell, plus Gc (2 per line); see A.2 |
| R1 | surface T on the line | 491–506 | `run` = solid run (496–498). If run ≥ 2: `x = RawSurfaceX(Tc, T_next, T0r)` (501), and if x > 0, `ReconFromX(..., h_spall_max/mu, flake_h_max/mu)` (502–506) |
| R1 | **the fallback** | 507–550 | x = 0 or run < 2. Face Newton on e* with `ef.h`, `ef.Tg` (518–525) via `SurfaceCellFlux(..., &Tf)`. Then q_net = q_robin − q_loss (529). If q_net > 0, T_s > T0 and k_c > 0: ell_q = k_c(T_s−T0)/q_net and h_f = c_f·ell_q, capped by /mu (531–543). Otherwise **noskin**: h_f = 0 (544–549) |
| R2 | the scan | 551–600 | nprof = min(Dmax, run−1) (552). Lambdas: `T_at_d` (553–557), `T0_at_d` (558–564), `T_at_zv` = `ScanTAtZv` anchored on T_s (566–569), `sigma_at_depth` → `ScanSigmaT(scan_g, T_z, T0_z, 1.0, conf, E, β, ν)` (572–578), `K_I` (579), `kic_pa_at_z`/`kic_pa_at_cell` = `ScanKIcPa` (581–586). `DepthScanH(nprof, 0, -(Lg+1), dx, 1.0, 1.0, a_f, true, ...)` (590–593). fire = h_col > 0 (595). **through** = no fire, score ≥ 1, no failure, the walk reached nprof and the run ended inside the window (598–600) |
| B | decision | 601–635 | runs only if `bk_on` (`spall.buckle_span > 0`, 488) and the line fired. x = mu·L_av on a through line (610–615), else mu·c_f·ell_eff (616–622). `BuckleRatio(scan_g, ν, β, T_s, T0_at_d(0), x, span)` (623). If !(bk ≤ 1) the line is withdrawn: fire = 0, through = 0 (625–633) |
| B | counters | 638–648 | bk_n_eval, bk_n_block, bk_n_inf, bk_n_block_{floor,wall,roof}, bk_max_pass. They are added once, not rank-summed. Registered at M:947–965, only under the key |
| — | winner per cell | 650–658 | a sheet keeps its highest-score line |
| R3 | recession | 660–762 | skip on !fire, or on fire && noskin && !through (678–688). f, Ho and L_av from the gather (691–700). h = L_av if through, else `rec.h_f` clipped to L_av (701–711). ell = `rec.ell_eff` or ∞ (712–713). `RecedeAlongLine` (714–717). Per-cell deltas into `upd` keyed by `FrKey` (721–742). **phi follow** on +z lines only (743–750). T_s is carried to the new face cell (752–761) |
| — | carry for non-fired cells | 763–775 | T_s is carried; `fr_id_remap` = \|ΣdH + E_left\|/E_left (774–775) |
| — | overlap and ledger | 777–790 | `fr_overlap`, `E_drop` (a cell voided under another line's books), then `ledger_E_flake += (E_left_tot + E_drop)·Vcell` (790) |
| R3 | **one write** | 792–861 | the +z column phi shift (814–821). A void sets f = 0, removed = 1, regime = regime_field = 1, RoP, spall_event, spall_thickness, RoP_spall, D = 0, with H frozen (825–837). A partial cell gets f = fn and H += dH (838–842). Witness fields fr_score and fr_mu (845–855) |
| — | CSV writer | 863–914 | `<plot_file>_removal_events.csv`. It keeps the column path's leading header and appends the line columns (873–882), plus `,bk` under the key (881, 909). Fixed placeholders: regime is always 1, Sp_top = 0, d_ell_rel/delta/q_ratio = 0, and ell_raw repeats ell_eff (894–903). `k_top` is read from the column path's `k_top_global[c.col]` (895) |
| — | tail | 918–934 | `InvertHtoT`, `RealFillBoundary`, `RefillTempGhosts` (from R:58–62), `UpdateSurfaceMaskFromPhi`, `RefreshKappaEff`, `fr_t_removal` |

**Thermo counters.** These are registered only under the key, at M:1062–1090:

| counter group | members |
|---|---|
| per class (floor / wall / roof) | `fr_n_*`, `fr_v_*`, `fr_fb_*`, `fr_smax_*`, `fr_noskin_*` |
| scalars | `fr_clip`, `fr_cap`, `fr_sheet`, `fr_multi`, `fr_overlap`, `fr_id_remap`, `fr_ncell`, `fr_through`, `fr_E_drop`, `fr_fresh`, `fr_misaligned` |
| timing | `fr_t_removal`, `fr_t_step` |

Its plotfile fields are `fr_f` (3 ghosts), `fr_score` and `fr_mu`, registered at M:1093–1095.

### A.2 What each stage reads

The gather is F:400–465: **G[o+0] = T** (`temp_mf`), **G[o+1] = T0** (`T_init_mf`, NaN→0 if absent, 461),
**G[o+2] = H**, **G[o+3] = solid** (`removed_mf` ≤ 0.5), **G[o+4] = f** (`fr_f_mf`). Per line:
**Gc[2l] = a_f** (`flaw_a_mf` under `weibull.enabled`, else `sp_a0`) and **Gc[2l+1] = k_c** (`kappa_mf`·e^{−α_D·D}).

| stage | uses | not used |
|---|---|---|
| R1 recon (501–506) | T[0], T[1], T0r = `phases[0].T_ref`, gas cap | H, f |
| R1 fallback (507–550) | T[0], k_c, the face's `h` and `Tg` from the list (`ef`) | — |
| R2 scan (551–600) | T[0..nprof], T0[0..nprof], a_f, `RockLaw()`, `ElasticConstants()` | **H, f** (the scan sees every solid cell as a whole cell) |
| B (601–635) | f (through: L_av = Σf·dx), ell_eff, T_s, T0[0] | — |
| R3 (660–762) | f[0..run), H[0..run), ell_eff, h | T |
| FrAssignTs (203–227) | `temp_old_mf` at c and c−e*, `removed_mf` at c−e* | — |

**Keys that gate it.**
- `surface_patch.face_list` (dispatch R:543). It also requires `face_list`'s parse conditions (M:994–1050): `convective_flame` with `h_expr`, `robin_form = surface`, the microstructure path, `surface.follow_mask = 1`, `sp_weibull` with `spall.enabled`, `local_thermoelastic`, `surface_T = reconstructed`, `removal_energy = flake`, and `amr.max_level = 0`.
- It **aborts** on each of these: `vapor.enabled`, `bit.enabled`, `foot_body_clearance`, `melt_aware`, `surface_normal`, `flake_length = flux`, `flake_fraction = weibull`, `size_volume = flux`, `flake_coherence_length > 0`, and the `connected_cluster` labeller (M:1026–1060).
- Inside: `spall.h_spall_max`, `spall.flake_h_max`, `spall.flake_c_f`, `spall.buckle_span` (B), `spallation.sp.n_segments`, `weibull.enabled` (flaw), `spall.removal_events_csv`, `kic.*`, `spall.stress_law`, `confining.p`, `sp_relief`.

**What differs from the column path in what fires.** The face list fires on `h_col > 0` alone (F:595). The column path also requires the detach gate (`Sp_cluster_id > 0` or `Sp_field ≥ 1`, R:1162–1165, 1558). So under the key, criterion site A (`UpdateSpAfterMechanics`) and `UpdateSpClusters` are **pure diagnostics**.

### A.3 The flame-specific assumptions, stated exactly

1. **Where the flame's h and T_gas enter.**
   - (a) `FrFaceFlux` F:148/153: `el.h`, `e.h` and the caller's `Tg` set the heat BC on every listed face.
   - (b) The removal fallback F:518–525: `SurfaceCellFlux(true, true, 1.0, ef.h, ef.Tg, k_c, dx, Tc, …, &Tf)`. The T_s it returns is the root of G(T_f − T_c) = h(T_g − T_f) − losses, with G = 2k/dx (M:3986–4010). Then q_net = h(T_g − T_s) − losses (F:529), and the skin length ell_q = k_c(T_s − T0)/q_net (F:533). **The flake thickness is therefore h-and-T_gas-defined on every fallback line.**
   - (c) `SurfaceGasTempCap()` (M:9117–9124 = max(`T_flame`, `jet_T_nozzle`)) caps T_s in `ReconFromX` (F:229/240, F:469/504). It is 0, meaning no cap, when `surface_patch.enabled = 0`.
   - (d) `e.h` comes from `FaceListHeat`/`FaceH` (`h_expr` at the face centre, M:3848–3901). `e.Tg` comes from `FaceListTg` (`flame_Tg[col]` or the march, M:2243–2276).
   - (e) The list exists only inside `surface_patch.radius` (M:3734–3736). It is built only when `surface_patch.enabled && !prescribed_T` (M:1883–1885), and `FrAssignTs` is called only from `BuildFlameColumns` (M:2228–2233).
2. **Where a profile that decays inward is assumed.**
   - `RawSurfaceX` (M:9000–9010) returns x > 0 only if T_c > T_next > T0. Its callers are F:234–235 and F:501.
   - `ReconFromX` (M:9020–9109) extrapolates an exponential: T_s = T0 + (T_c − T0)·x/(1−e^{−x}), ell = dx/x, h_f = c_f·ell (M:9036–9053).
   - `RemapProfile`'s sub-cell blend b = (1−e^{−r/ell})/(1−e^{−dx/ell}) (R:240–247) is fed ell = `rec.ell_eff` (F:712–713). On a fallback line that is ell_q, a flux length rather than a measured decay.
   - `ScanTAtZv` lerps from T_s at z = 0 to the first cell centre (R:75–86).
   - `DepthScanH` walks from the surface and **stops at the first sample with K_I < K_Ic** (R:171–196). Its 8-step bisection assumes pass → fail with depth (R:197–225). A subsurface peak is never reached.
   - `BuckleRatio`'s strain is evaluated at the **surface** temperature T_s (F:623).
   - The through rule (F:598–600) assumes the walk only ends at the run's end when the whole run passed.
3. **When `RawSurfaceX` returns 0 ("flat").**
   - In `FrAssignTs`, `e.Ts = 0` and `tsfb = 1` (F:236–248). The BC then solves the face Newton per face at the flame's h and T_g (F:144–151).
   - In removal, the fallback replaces `ReconFromX` entirely (F:507–550). Its `rec.flat = true`, its ell_eff is ell_q or ∞, and ReconFromX's flat h_f = dz is **not** read (F:514).
   - If q_net ≤ 0 or T_s ≤ T0, the line is **noskin**: h_f = 0, and it recedes only if `through` (F:544–549, 678–683).
   - ⚠ Under a beam-only deck h = 0, so q_net = −losses < 0 on every flat line, and every flat line is noskin. Volumetric heating makes T_next > T_c the normal state (MODEL_CARD §8), so **on the beam today the face list would recede only through `through`.**
4. **mu and where /mu is applied.**
   - **mu** = n·e*. Here n is the unit outward normal of the smoothed solid fraction (`FaceEntry.nx/ny/nz`, set in `BuildFaceList` from `FrSmoothedGradient`/`FaceWeight`, M:3752–3755), and e* is the exposed face that maximises it (F:113–125). mu = 1 on sheets and misaligned cells (F:124, 394).
   - **/mu is applied to the two NORMAL lengths only:**
     - `h_spall_max/mu` and `flake_h_max/mu` in R1 (F:241, 505) and in the fallback (F:535–536);
     - Dmax and Lg (F:377–380).
   - **x·mu** converts the along-line thickness back to normal thickness in B (F:614, 620). The CSV prints h_n = mu·h_line (F:905).
   - mu is **not** applied to dx in the reconstruction: the two-cell ratio along the line already returns ell/mu (F:19–24).
   - That cancellation holds only if T depends on the normal distance alone, i.e. a surface flux. Beer–Lambert deposits along −z only, from the top face of each column's first solid cell (M:2033–2061). So a wall line (dir 0–3) sees a radial, Gaussian-shaped profile, not exp(−normal/ell).

### A.4 Lines a beam port must change or bypass

| # | Site | Lines | Why |
|---|---|---|---|
| 1 | parse gate | M:994–1000 | requires `convective_flame` with `h_expr` and `robin_form = surface` |
| 2 | parse aborts | M:1033 (`vapor.enabled`), M:1037 (`melt_aware`, which is deleted rather than lifted, see below) | vapour is part of the MMW step |
| 3 | list build gate | M:1883–1886 | `fl_on` needs `surface_patch.enabled` and `!prescribed_T`. `BuildFlameColumns` is the only caller of `FaceListHeat`/`FrAssignTs`/`FaceListTg` (M:2228–2235) |
| 4 | patch radius | M:3734–3736 | faces outside `surface_patch.radius` are not listed. A beam deck has no patch |
| 5 | BC face temperature | F:138–156, F:162–253 | T_s comes from the flame's h/T_g or from `ReconFromX`. Under a beam the losses alone set the face Newton |
| 6 | R1 reconstruction | F:229, 234–242, 469, 500–506; M:9000–9124 | the gas-temperature cap, the decaying-inward test, the exponential T_s and h_f |
| 7 | R1 fallback | F:507–550 | `ef.h`, `ef.Tg`, ell_q = k(T_s−T0)/q_net, noskin. This is the flame-defined skin length that the beam must replace with a volumetric one, e.g. built from 1/α and the in-depth profile |
| 8 | R2 scan anchor and walk | F:551–600; R:64–103, R:159–227 | T_s anchors the profile. The walk stops at the first failure. The search depth is h_spall_max/(mu·dx), so a deep subsurface peak needs Dmax ≥ the peak depth |
| 9 | B thickness and strain | F:607–635 | x = c_f·ell_eff; ε_th at T_s |
| 10 | R3 heat remap | F:701–717; R:236–260 | h = `rec.h_f` and the exponential sub-cell blend with ell_eff |
| 11 | window sizing | F:368–369, 377–380 | Dmax and Lg from h_spall_max and flake_h_max |
| 12 | the beam's sight of a partial cell | M:1924, M:2048–2061; `SolidFacePath` M:4328–4331 | the kernel and `BeamSource` test `removed` only. A cell with f < 1 (`fr_f_mf`) absorbs as a full cell over path (n_above + ½)·dz. The kernel does not read `fr_f` at all |

**Two closed-study branches that sit at those sites, to delete in place.** Neither appears in FaceRemoval.H (`grep` finds no `melt_aware`, `MeltSkeleton`, `f_skel`, `surface_normal` or `sn_c*` there). Both live in the shared helpers it calls.

- **melt-aware.**
  - **In the helper:** `ScanSigmaT`'s `f_skel` parameter (R:131, 133, 138, 143); F:577 passes 1.0.
  - **Column path:**

    | piece | lines |
    |---|---|
    | the Λ_L gather | R:1017–1057 |
    | `LL_at_didx` | R:1367–1385 |
    | the skeleton factor | R:1472–1479, 1501 |

  - **In M:**

    | piece | lines |
    |---|---|
    | members and `MeltSkeleton` | M:8599–8618 |
    | parse | M:5510–5526 |
    | face-list abort | M:1037 |

- **surface_normal.**
  - **In the helper:** `DepthScanH`'s `sn_cth`/`sn_icth` (R:160–161, 174, 179, 202–204, 211, 215); F:591 passes 1.0, 1.0.
  - **Column path:**

    | piece | lines |
    |---|---|
    | the slope cosines | R:549–558, 1223–1231 |
    | uses in the stress lambda | R:1456–1463 |
    | uses in the firing record | R:1560–1569 |
    | coherence-cap bypass | R:1797–1810 |
    | comments | R:2000–2009 |
    | `sn_recession_superseded` | R:2403–2424 |

  - **In M:**

    | piece | lines |
    |---|---|
    | registration | M:938–940 |
    | face-list abort | M:1038 |
    | parse and message | M:5470–5499 |
    | `ColumnSurfaceCos` | M:5958–5990 |
    | site A | M:6085–6090, 6189–6194 |
    | members | M:8552–8557, 8593–8598, 8619–8622 |

  - ⚠ Two **tree-gate** decks set `surface_normal = 1`, and both are beam decks: `sp_meier_pilot/input_2d_dev` (dev2d) and `unit/beam_void_closure` (`gate_keys.sh`). M:5486–5490 explains that the key is a message rather than an abort for that reason. Deleting the key means re-keying both decks.

---

## B. The column path (Removal.H)

CODE_CENSUS C.3 numbered stages against a 3,654-line snapshot. Since then:
- P2-10 deleted the lateral criterion (old stage 5), the lateral thermal shift (old stage 17) and its diagnostics (old stage 18).
- P2-16 deleted the constraint depth.

No `lat_`, `constraint` or `d_ref` remains (`grep`). The stages are re-anchored below.

**Disposition key:**
- **FL** means the face list already has it; the equivalent is named.
- **PORT** means the beam, vapour, melt or feet logic the face list lacks.
- **DEL** means column-only machinery, a dead remnant, or a column-only diagnostic.
- **MOVE** means the face list depends on it.

### B.1 Shared helpers (lines 1–260). These must MOVE, not be deleted

| Helper | R: lines | Callers in R: | Callers in F: | Other callers (M:) |
|---|---|---|---|---|
| `ElasticConsts` / `ElasticConstants()` | 23–47 | 692, 750 | 471 | 6022, 6059 (site A) |
| `RefillTempGhosts` | 49–62 | 2571 | 926 | — |
| `ScanTAtZv` | 64–103 | 1397 | 567 | — |
| `ScanKIcPa` | 105–127 | 1411, 1423 | 582, 585 | — |
| `ScanSigmaT` | 129–148 | 1527 | 577 | — (drop `f_skel` on the move) |
| `DepthScanH` | 150–227 | 1543 | 591 | — (drop `sn_cth`/`sn_icth` on the move) |
| `RemapProfile` | 229–260 | 2108 | 288 | — |

These already live outside Removal.H and stay where they are:
- `BuckleRatio` (F:315–326);
- `LawSigma`/`RockLaw` (M:6591–6616);
- `RawSurfaceX`/`ReconFromX`/`SurfaceGasTempCap` (M:9000–9124);
- `SurfaceCellFlux` (M:3928–4045);
- `FrInitFraction` (M:3585–3597);
- `InvertHtoT` (M:6698–6760);
- `UpdateSurfaceMaskFromPhi` (M:6793–6806).

### B.2 `UpdateRemovalAfterCohesive` (R:262–2578), stage by stage

Stages B0–B7 run **before** the dispatch, so they run under the face list too.

| Stage | R: lines | Old C.3 | What it does | Keys | Disposition |
|---|---|---|---|---|---|
| B0 guards | 262–287 | 0 | `removal_enabled`, the microstructure assert, `SurfaceOwner`, domain sizes | — | **MOVE** to the dispatch site (the face list needs the same guards) |
| B1 damage-law h_spall | 289–321 | 0 | h_spall = C_h√(α t). It **runs and can abort under the face list but is unused there**. Also `stress_mf` FillBoundary | `spall.C_h/alpha/t_spall`, `spallation.model` | **DEL** (damage-law removal; see B.4) |
| B2 vapour calibration | 323–343 | 0 | `A_beam`, `T0_v`, `P_now`, 1−R, ω0/x0/y0, the A_beam abort | `vapor.*`, `beam.*` | **PORT** (vapour; see C) |
| B3 per-step field reset | 345–368 | 0 | sets to 0: `regime`, `RoP`, `spall_event/thickness/RoP_spall`, `vapor_*`, `regime_field` | `spall.enabled`, `vapor.enabled` | **MOVE**. The face list writes these fields on voided cells (F:830–835) and relies on this reset |
| B4 locals | 369–385 | 0 | `sp_weibull_spall`, `ptabs`, `col_index` | — | DEL |
| B5 **feet mechanical clearance** (D2c) | 387–506 | 1 | an annulus column proud of its pad: cells above `foot_clear_k` go void with regime 4 and D = 0, H frozen; the column phi shifts; `foot_mech_E/vol/cols` booked | `surface_patch.feet` + `foot_body_clearance` (431–432) | **PORT**. Re-express on the list: each cell above the pad index gets f = 0, removed = 1, regime 4, D = 0, `foot_mech_E += H·V`, without touching H. ⚠ As written it would not zero `fr_f_mf`, which is why the face list aborts the key (M:1035–1036) even though this block sits *before* the dispatch |
| B6 PASS 1 `k_top_global` | 508–537 | 2 | the column-top reduce (B3 copy) | — | **MOVE** (F:895 reads it for the CSV `k_top`) or replace with `ColumnTopSolid` (M:4589–4611) |
| B7 **dispatch** | 539–547 | — | → `RecedeAlongLines` | `surface_patch.face_list` | MOVE (it becomes the only call, at M:1454) |
| B8 Step 20 cosines | 549–558 | 2 | `ColumnSurfaceCos` | `spall.surface_normal` | **DEL** |
| B9 PASS 1.5 D_max_col | 560–626 | 3 | per-column max D over `spall.sample` | damage law, `spall.sample` | **DEL** |
| B10 candidate buffers | 628–667 | 4 | spall/vapour candidate vectors. 639–658 is a P1b support-rule comment **with no code** | — | DEL (vapour buffers: PORT as per-cell state) |
| B11 Step 16d captures | 669–706 | 4 | phase-0 E, β, ν, T_ref, κ; `surface_patch.T_f` | `prescribed_T` | DEL (FL has `ElasticConstants`, F:471) |
| B12 Step 19 selector | 708–741 | 6 | T-face vs `stress_mf`; aborts | `stress_source`, `el.type` | DEL (FL requires `local_thermoelastic`) |
| B13 confinement offset | 743–750 | 6 | `ElasticConstants().conf` | `confining.p` | FL (F:472 `conf_off`) |
| B14 gates and buffers | 752–799 | 6 | `HColEvent`; `pin_on`, `rev_on`, `inc_on`, `inc_active`; `col_rho/Cp`; `rev_*` | `robin_form = pinned`, `removal_events_csv` | DEL |
| B15 recon gate and `rc_*` | 800–848 | 6 | `recon_on`; `rc_Tsurf/ell/ellq/h_f/Tbar/flag/qrat` | `spall.surface_T` | DEL (FL: `LineOut.rec`, F:482) |
| B16 PASS 1b T gather | 850–902 | 7 | `solid_run_col`, `Dmax_col`, `T_col_prof` | `local_thermoelastic` / `prescribed_T` | FL (F:400–465, column 0) |
| B17 T0 gather | 904–947 | 7 | `T0_col_prof` | — | FL (F:445–446) |
| B18 H and phi gather | 949–1015 | 7 | `H_col_prof`, `phi_top_col` (`remap_on`) | `surface_T` + `removal_energy = flake` | FL (F:447) |
| B19 Λ_L gather | 1017–1057 | 7 | melt profile | `spall.melt_aware` | **DEL** (closed study) |
| B20 PASS 2 per-column criterion (site B) | 1059–1733 | 8 | owner-rank loop over column tops; the parts follow | — | — |
| ↳ owner reads | 1065–1145 | | `rev_*`, `col_rho/Cp`, xc/yc/zc | | DEL |
| ↳ detach gate | 1162–1165 | | cluster id > 0, or `Sp_field ≥ 1` | `weibull.enabled` | **DEL** (FL fires on h_col alone) |
| ↳ flux size factor | 1197–1204 | | `WeibullSizeFactor(ell_q³)` | `weibull.size_volume = flux` | DEL (FL aborts it) |
| ↳ scan setup | 1211–1298 | | `max_steps` (`solid_run`), `T_at_didx`, `T0_at_didx` | | FL (F:551–564) |
| ↳ recon read | 1300–1365 | | `SharedSurfaceRecon`, `rc_*`, `q_fit/q_real` (flame h/T_g, 1348–1360) | | FL (F:500–550) |
| ↳ T_at_zv, K_Ic | 1387–1424 | | `ScanTAtZv`, `ScanKIcPa`; the prescribed-T face (1399–1400) | | FL (F:566–586) |
| ↳ `k_i_at_zcrack` | 1426–1534 | | melt skeleton (1472–1479); **`stress_mf` branch** (1480–1503); `ScanSigmaT` (1527) | | FL for the T branch (F:571–580). DEL the `stress_mf` branch and the melt branch |
| ↳ `DepthScanH` + `crack_event_count` | 1536–1556 | | edge-triggered `cec` | | FL (F:590–593). DEL `cec` |
| ↳ fire record | 1558–1608 | | h_vert = h_col/cos; `h_spall_field` running max; `HColEvent` | | DEL |
| ↳ P1b comment | 1610–1628 | | no code | | DEL |
| ↳ damage-law branch | 1632–1661 | | D ≥ D_th, principal-stress ratio ≥ ξ_s, h = h_spall | `spallation.model = damage_law`, `spall.prescribed_ratio` | **DEL** (decision: B.4) |
| ↳ **vapour branch** | 1664–1731 | | H_vap, L_tot, three rate forms, fire if H_top ≥ H_vap | `vapor.enabled`, `vapor.prescribed_rop` | **PORT** (redesign; see C) |
| B21 reductions | 1734–1751 | 9 | the candidate and `rev_*` sums | — | DEL |
| B22 recon reductions | 1753–1786 | 9 | `rc_*`; `rc_dell` from `rc_ell_prev` | — | DEL |
| B23 coherence cap | 1788–1842 | 10 | recede at most `coh_cells` below the shallowest neighbour | `spall.flake_coherence_length` | **DEL** (FL aborts it; set to 4 mm in dev2d/beam_void_closure) |
| B24 `h_col_events` CSV | 1844–1892 | 11 | rank-serialised first-firing log | `spall.h_col_events_csv` | DEL |
| B25 bit | 1894–1931 | 12 | regime 3 target per column | `bit.enabled` | DEL, or PORT if the bit is kept (not in the plan; FL aborts it at M:1034) |
| B26 winner loop | 1933–2044 | 13 | regime 1 / 2 / 3 priority; whole-cell h = dz; `flake` h = `rc_h_f` | `vapor`, `bit`, `removal_energy` | FL for regime 1 (F:650–658, 701–711). **PORT** the spall-versus-vapour priority |
| B27 remap build | 2046–2117 | 14 | `RemapProfile` on regime 1; `nvoid_pred`; `E_flake_step`. `E_sup_local` (2068) is **dead**: no stage assigns regime 5 | `remap_on` | FL (F:266–302, 714–717) |
| B28 pin record | 2119–2138 | 14 | `pin_T/t/rho/Cp` on regime 1 | `robin_form = pinned` | DEL (FL requires `robin_form = surface`) |
| B29 PASS 3 apply | 2140–2297 | 15 | column phi shift, phi-crossing void flip, diagnostics, H frozen (2242–2257), remapped H with E_bnd (2268–2295) | — | FL (F:792–861). **PORT** the `vapor_event/thickness/RoP_vap` writes (2167–2173, 2235–2240) with vapour |
| B30 invariants | 2299–2426 | 16 | voids = predicted, with an out-of-rock abort; `sn_recession_superseded` | `inc_active` | DEL (FL: `fr_id_remap`; the void count is exact by construction in `RecedeAlongLine`) |
| B31 ledger sink | 2428–2436 | 14/15 | `ledger_E_flake += E_flake + E_bnd + E_sup` | `remap_on` | FL (F:790) |
| B32 `removal_events` CSV | 2438–2547 | 19 | the column header, recon columns, `ell_q`, regime-4 feet rows | `spall.removal_events_csv` | FL (F:863–914). **PORT** the regime-4 row (2501–2507) with the feet and the regime-2 rows with vapour |
| B33 AMR coverage | 2549–2550 | 20 | `CheckSurfaceCovered` (M:4442–4477) | `amr.max_level > 0` | PORT when AMR comes to the list (FL aborts AMR, M:1042). Note that the dispatch **skips** it today |
| B34 tail | 2552–2577 | 20 | `InvertHtoT`, ghosts, mask, κ_eff | — | FL (F:918–934) |

### B.3 What the face list lacks, by plan item

| Plan item | Column-path logic | What it must become on lines |
|---|---|---|
| beam | none in Removal.H; the beam is the kernel's `BeamSource` | A.4 rows 1–12. Vapour's rate (R:1700–1710) is the only removal-side beam read, and it should not be ported as is (C.3) |
| vapour | B2, B3 (355–360), B20 vapour branch, the B26 priority, the B29 diagnostics | the face cell's regime-2 test on H (not on `H_top` of a column), a recession charged against H (C.3), and per-cell regime 2 in `regime`/`regime_field`/`vapor_*` |
| melt | B19, the skeleton in B20, `MeltSkeleton` (closed study) | nothing to port. **No melt removal exists** on either path (C.2). Building one is new work on the enthalpy path |
| Meier's feet | B5 (387–506), the regime-4 CSV row (2501–2507), `foot_mech_E` (M:892–903, ledger M:4310) | B.2 B5, plus the zeroing of `fr_f_mf` and the face-list carry (`fr_carry`) for the new face cell |
| bit (optional) | B25, regime 3 in B26 | not in the plan |

### B.4 Decision the deletion forces

The damage-law removal (B1, B9 and the damage branch at 1632–1661) is the Hu path:
- `unit/spall_event` (`el.type = static`, no `spallation.model` key);
- `validation/hu/hu_end_to_end`;
- `unit/regime_low_high_power` (damage law + vapour).

The face list aborts any model other than `sp_weibull` (M:1045–1046). Deleting Removal.H deletes the damage-law removal and these tests, unless someone first decides to port it.

---

## C. Vapour and melt

### C.1 Where the vapour front lives

| Piece | Location |
|---|---|
| keys | `vapor.enabled` M:5816–5817; `vapor.A_beam` M:5820–5826 (default 0.5πω0², M:1112–1123); `vapor.T0` M:5828–5833 (default `T_amb`, M:1124–1125); `vapor.prescribed_rop` M:5837–5845 |
| phase constants | `phaseN.L_m` M:390, `T_vap_lo/hi` M:391–392, `L_v` M:393. Passed into the per-phase H(T) tables at M:7039–7042 and 7063–7066 |
| documentation | M:8641–8662: RoP_vap = (1−R)P/(ρ L_tot A_beam), L_tot = c_p(T_vap_lo − T0) + L_m + L_v |
| fields | `regime`, `RoP` M:534–535; `vapor_event`, `vapor_thickness`, `RoP_vap` M:543–548; `regime_field` (per-cell regime, sp_weibull with spall or vapour) M:551–560 |
| calibration, reset | R:323–343, R:355–360, R:364–368 |
| front update | R:1664–1731. H_vap = `H_of_T(T_vap_lo)` (1667); L_tot (1675–1676). Three forms: prescribed (1678–1684); sp_weibull local Gaussian (1700–1710): Q = 2P/(πω0²)·e^{−2r²/ω0²} at the waist, rop = (1−R)Q/(ρ L_tot); damage-law bulk (1721–1722): (1−R)P/(ρ L_tot A_beam). Fires if `H_top ≥ H_vap && rop > 0` (1725–1730); h_v = rop·dt |
| regime choice | R:1941–1965: sp_weibull gives spall priority (1953–1954); damage law takes the larger h (1958–1964); regime 2 recorded (2041) |
| apply | R:2201–2240: phi −= h_v, so every cell whose phi crosses zero is voided as regime 2, with `vapor_*` written. **H is not touched** (2242–2257), and there is no remap for regime 2 (2268) |
| face list | none. It aborts `vapor.enabled` (M:1033) |

### C.2 Where melt lives

- **Latent heat is already on the enthalpy path.** `L_m` and `L_v` are plateaus of the per-phase H(T) table, and `InvertHtoT` (M:6698–6760) returns T and Λ_S/Λ_L/Λ_V (fields M:277–279). Every path's H update absorbs them with no removal involved.
- **The melt-aware skeleton** f(Λ_L) = clamp(1 − Λ_L/λ_crit) (M:8599–8618, `melt.lambda_crit` M:5521–5526). It multiplies the scan's driving stress (R:1472–1479, 1501, via `ScanSigmaT`'s `f_skel`). It only **stalls** spallation under a molten cap (M:5510–5518). It adds or removes no energy. Its test is `unit/sp_melt_skeleton`.
- **No melt removal or evacuation exists anywhere.** Neither path has a melt regime.

### C.3 Energy ledger terms and who books into them

With `energy_ledger.enabled` the following are registered at M:806–812:

| term | definition |
|---|---|
| `ledger_E_H` | Σ H·V over **all** cells, void cells included (their H is frozen) |
| `ledger_E_src` | += dt·(P_beam + P_robin − P_loss) (M:4295) |
| `ledger_err` | the residual of the two rows above |
| `ledger_P_beam` | `beam_P_dep` from `BeamEnergyBalance`, the same `BeamSource` calls (M:4227) |
| `ledger_P_robin`, `ledger_P_loss` | the Robin and loss powers |
| `ledger_surface_missing_cols` | see M:4218–4223 |

Under `spall.surface_T` (M:818–825):

| term | definition |
|---|---|
| `ledger_E_flake` | the sink |
| `ledger_E_H_solid` | solid cells only |
| `ledger_err_solid` | (ΔE_solid − (E_src − E_flake − foot_mech_E))/E_src (M:4310) |

Under feet clearance (M:892–903): `foot_mech_E` (plus `foot_mech_vol`, `foot_mech_cols`).

| Regime | Path | Books into | Lines |
|---|---|---|---|
| 1 spall, `remap_on` | column | `ledger_E_flake` (E_flake_step + E_bnd) | R:2108–2109, 2292, 2433–2435 |
| 1 spall, no remap (whole cell) | column | **nothing** (cell leaves the solid set unbooked) | R:2431 gate |
| 1 spall | face list | `ledger_E_flake` (E_left + E_drop) | F:777–790 |
| 2 vapour | column | **nothing** | B29 voids without a sink; R:2431–2435 counts regime 1 and 5 only |
| 3 bit | column | **nothing** | — |
| 4 feet clearance | column | `foot_mech_E` | R:486, 495–496 |
| 5 support (P1b) | column | `E_sup_local`, which is **dead**: regime 5 is never assigned | R:2068, 2217 |

So no ledger term exists for vaporisation or melt. `ledger_err` cannot see any removal, because void H stays in E_H. `ledger_err_solid` would drift by the H of every vapour-voided cell.

### C.4 Can the vapour front spend beam power that the enthalpy field also received? (MODEL_CARD §8, "unverified")

**Yes.** The evidence:

1. **The beam goes into H everywhere.** The kernel does `H = H_old + dt·(div_kgrad + Q − q_loss/dz)` with `Q = BeamSource(P_now, …)` on every non-void cell (M:2048–2066). `BeamSource` returns (1−R)·Pi·(e^{−a s_top} − e^{−a s_bot})/dz (M:4918–4948). So the absorbed (1−R)P that falls inside the domain is added to H. `ledger_P_beam` books it as a source (M:4227, 4295).
2. **The vapour front turns the same (1−R)P into a recession rate**: rop_v = (1−R)·Q_local/(ρ L_tot) (R:1706–1709), or (1−R)P/(ρ L_tot A_beam) (R:1722). It charges ρ·L_tot = ρ[c_p(T_vap_lo − T0) + L_m + L_v] per removed volume.
3. **Nothing subtracts that charge from H.**
   - PASS 3 voids the cells and freezes their H (R:2242–2257).
   - The remap and the flake sink apply to regime 1 only (R:2075, 2268, 2431–2435).
   - The removed rock's enthalpy stays in `ledger_E_H` (M:4134), and its sensible and latent heat up to T_vap_lo was itself delivered by the beam through Q.
   - So the beam pays twice: once into H, and once as ρ·L_tot·h_v of recession.
4. **The cells removed need not hold the charged enthalpy.** The gate tests only the column-top cell's H against H_vap (R:1725). h_v = rop·dt can exceed dz, and the phi shift then voids deeper cells whatever their H (R:2201–2212).
5. **There is no compensating term.** A grep for `ledger_` finds no vapour sink (C.3).

**Secondary defects in the same block:**
- `Q_local` uses the Gaussian waist formula even under `beam.profile = uniform` (R:1703–1708, against `BeamIncidentIntensity`'s uniform branch at M:4964–4970).
- It ignores Gaussian spreading at the actual surface depth.
- It ignores radiation losses at the vapour temperature (MODEL_CARD §8).
- The tests that exercise it are `unit/sp_v_n_regime` (sp_weibull + vapour + beam) and `unit/regime_low_high_power` (damage law + vapour + beam). Both run only on the column path.

---

## D. Column-top state in MMWSpalling.H

| State | M: lines | Producers | Consumers | Used by FaceRemoval.H? |
|---|---|---|---|---|
| `BuildSurfaceTColumns` | 2866–3033 (called M:1434) | gathers `Ttop`, `T1`, `Ktop` (locals, 2885–2922) | Removal.H 1200, 1316, 2077, 2525; site A M:6159, 6290; `SurfaceRobinT` (below); ledger `surf_q_mean` M:4288–4291 | **No.** F re-reads T(c), T(c−e*) itself (F:203–227, 444) |
| `surf_T/ell/ellq/ell_raw/hf/tbar/flag_col`, `surf_x_hold` | decl. 8834–8852; filled 2924–3032 | `BuildSurfaceTColumns` | as above; `surf_x_hold` is the cross-step flat-hold (2961) | No |
| `sp_len_fb_h0/q0` | 8794–8795; reg. 227–229 | 2872–2873, 2987, 3013 | thermo only | No (FL aborts `flake_length = flux`) |
| `SharedSurfaceRecon` | 3046–3062 | — | Removal.H 1316; site A M:6159 | No |
| `SurfaceRobinT` | 3069–3074 | `surf_T_col` | kernel M:2012, march M:2610, pin counters M:3159, ledger M:4174 | No. On every listed cell the kernel (M:2015–2020), ledger (M:4180–4182) and march (M:2618–2626) override it with `FrFaceFlux` |
| `rc_*` buffers | locals R:812–840; `rc_ell_prev` member M:8858 | Removal.H 1322–1364, 1756–1786 | Removal.H CSV 2515–2533 | No (the F CSV writes zeros or ell_eff, F:901–903) |
| `ColumnTopSolid` | 4589–4611 | removed_mf reduce | `AdvanceMicrostructure` M:1828–1830 → `BuildFaceList` (M:1885, 3709), `FeetDescent`, `CheckSurfaceCovered` M:4449, Removal.H 444 | Indirectly: `BuildFaceList` uses k_top only to class a +z face TOP or UP (M:3765–3766) |
| `k_top_global` | local R:512–537 | PASS 1 | passed into `RecedeAlongLines` | **Yes, for the CSV `k_top` column only** (F:895) |
| `FrCell` list | F:95–131 from `fl_faces` | `BuildFaceList` M:3694–3844 builds from **`removed_mf`** (open faces, M:3721–3760) and `fr_f_mf` (normal) | — | It is the mask, confirmed. No column-top dependence |
| `ColumnSolidRunBelowTop` | 4677–4714 | — | Removal.H 868; `BuildSurfaceTColumns` 2883 | No (F computes `run`, F:496–498) |
| `GatherColumnSurface` | 5907–5940 | a duplicate of `ColumnTopSolid` | `BuildSurfaceTColumns` 2882; site A 6074; clusters 6378 | No |
| `ColumnSurfaceCos` | 5958–5990 | — | Removal.H 555; site A 6087 | No (`surface_normal`, DEL) |
| `CheckSurfaceCovered` | 4442–4477 | — | Removal.H 2550 only | No (and the dispatch skips it) |
| **Sp_field site** (`UpdateSpAfterMechanics`) | 5992–6343; `Sp_field_mf` reset/write 6111–6117 | site A | Removal.H detach gate 1163–1164, `rev_Sp` 1126, `HColEvent` 1605; `UpdateSpClusters` 6386–6393 | **No** (F CSV `Sp_top = 0`, F:896; F fires on h_col alone, F:595) |
| `UpdateSpClusters` / `Sp_cluster_id_mf` | 6353–6557; reg. 576–577 | — | Removal.H detach gate 1100, 1163 only | No |
| `h_spall_field_mf`, `crack_event_count_mf` | reg. 583–591 | written only by Removal.H 1094–1095, 1551–1554, 1580–1586 | plotfiles | No |
| `h_col_events` CSV | Removal.H 1844–1892; flags M:5536, 8904, 8913 | — | — | No |
| pin state `pin_T/t/rho/Cp/T_eff` | 7621; `BuildPinEffective` 3076–3115 | Removal.H 2121–2138 | the kernel and march under `robin_form = pinned` | No (FL requires `surface`) |
| `foot_clear_k` | 7681 | `FeetDescent` M:2411–2417 | Removal.H 442–505 only | No |
| `phi_mf` | — | Removal.H B5/B29; F:814–821 (+z lines only) | `UpdateSurfaceMaskFromPhi` M:6793–6806, AMR phi rebuild M:4341–4401, `InitPhi` M:6849 | F writes it for the column shift only. Under the face list's required `follow_mask` the decisions read the mask |

**Net result.** Of the column-top state, the face list uses exactly three things: `k_top_global` (for one CSV column), `ColumnTopSolid` (through `BuildFaceList`'s TOP/UP class), and the B3 field reset. With Removal.H gone and the face list mandatory, `BuildSurfaceTColumns` and `SharedSurfaceRecon` feed only site A's diagnostics and `surf_q_mean`. So do `SurfaceRobinT` (overridden on every listed cell), the Sp_field/cluster pair, `h_spall_field`, `crack_event_count`, `rc_ell_prev`, the pin state and `ColumnSurfaceCos`.

---

## What deleting Removal.H requires first

**Helpers to move** (out of Removal.H, into FaceRemoval.H or a shared fragment):
- [ ] `ElasticConsts`/`ElasticConstants` (R:23–47). Site A also calls it (M:6022, 6059).
- [ ] `RefillTempGhosts` (R:49–62).
- [ ] `ScanTAtZv` (R:64–103) and `ScanKIcPa` (R:105–127).
- [ ] `ScanSigmaT` (R:129–148). Drop `f_skel` on the move.
- [ ] `DepthScanH` (R:150–227). Drop `sn_cth`/`sn_icth` on the move.
- [ ] `RemapProfile` (R:229–260).

**Prologue pieces the face list silently depends on:**
- [ ] Guards (R:272–276).
- [ ] The per-step reset of `regime`, `RoP`, `spall_*`, `vapor_*` and `regime_field` (R:345–368).
- [ ] `k_top_global` for the CSV (R:508–537), or `ColumnTopSolid`.
- [ ] Repoint `M:1454` to call `RecedeAlongLines` directly.

**State to re-home or delete with it:**
- [ ] `foot_clear_k` consumer (feet port).
- [ ] Pin state (B28 and `BuildPinEffective`) when `robin_form = pinned` is retired.
- [ ] `BuildSurfaceTColumns`, `SharedSurfaceRecon` and `SurfaceRobinT`, if no surviving deck uses the column BC. Repoint site A to the list's T_s, or retire site A.
- [ ] The `Sp_cluster_id`, `h_spall_field`, `crack_event_count` and `h_col_events` producers or fields.
- [ ] `ColumnSurfaceCos`, `ColumnSolidRunBelowTop`, `GatherColumnSurface`, `CheckSurfaceCovered`, `rc_ell_prev`.
- [ ] The dead `E_sup_local` (R:2068, 2217).

**Ports to complete:**
- [ ] **Beam onto the list:** every row of A.4 (parse gate, list build without a patch, the T_s/ell/h_f/remap volumetric meaning, the scan for a subsurface peak, the beam's sight of `fr_f`).
- [ ] **Vapour:** rebuilt on H, not on the column-top beam formula. It must not double-spend the beam (C.4), and it needs its own ledger sink.
- [ ] **Melt:** new work. Nothing to port; the skeleton is a closed study, so delete it.
- [ ] **Feet clearance:** regime 4 on the list (B.2 B5), including zeroing `fr_f`, a `foot_mech_E` sink and the CSV row. Lift the abort at M:1035–1036.
- [ ] **Decide:** damage-law removal (Hu: `spall_event`, `hu_end_to_end`, `regime_low_high_power`), bit (regime 3), and AMR coverage (`CheckSurfaceCovered`).

**Gate inputs that exercise only the column path** (`wit.sh` + `post_all.sh`; keys read by `work/gate_keys.sh`):

**No input in the 2639 tree gate sets `surface_patch.face_list = 1`.** The face list is covered only by the argv-keyed units `unit/face_list`, `face_removal`, `face_buckle` and `face_heat`.

| Gate leg | Column-only features it uses |
|---|---|
| dev2d `validation/meier/sp_meier_pilot/input_2d_dev` | beam 11.4 kW, `surface_normal = 1`, `flake_coherence_length = 4 mm`, no patch |
| `unit/beam_void_closure/input` | beam, `surface_normal = 1`, coherence 4 mm, no patch (exits 1 by design) |
| `unit/scalar_flaw/input` | beam, `melt_aware = 1`, no patch |
| `unit/robin_pinned`, `unit/robin_jet`, `unit/jet_enthalpy/input_quarter`, `unit/robin_feet` (both inputs) | `robin_form = pinned`, `melt_aware = 1` |
| `unit/jet_enthalpy/input_strip` | `robin_form = cell`, `melt_aware = 1` |
| `unit/robin_face/input_ledger` | `robin_form = face`. `input_form` (cell form) sets no removal key |
| `unit/spall_event` (both inputs) | damage-law removal, `el.type = static` |
| `studies/s1_surface_resolution` (smoke) | `melt_aware = 1`, no patch |
| `unit/jet_d2c` (run by `post_all.sh`) | legs CLD/CLD3 set `foot_body_clearance = 1` (its test, lines 26 and 136) |

Also on the "removal" row of WHICH_TESTS but not in `wit.sh`:
- `unit/removal_increment`: pinned, plus `input_presc` with `prescribed_T`.

**Column-feature unit tests outside the gate** (`work/test_paths.sh`, a grep-level tag):

| feature | tests |
|---|---|
| vapour | `sp_v_n_regime`, `regime_low_high_power` |
| melt skeleton | `sp_melt_skeleton` |
| feet | `pad_clear_ledger` |
| `prescribed_T` + `stress_mf` | `validation/rossi/sp_rossi_damage_profile` |
| damage law | `validation/hu/hu_end_to_end` |
| coherence | `surface_T_recon`, `gradient_blind`, `seed_void`, `validation/p0_steady_column`, `kant/sp_kant_onset_meiercfg` |

**Source-lock tests that open Removal.H by path, which break on deletion:**
- `unit/lateral_flux_lock/test` (:46, :547, :565, :640);
- `unit/buckle_lock/test` (:31–32, :49–50, :124–125, clause B2 "Removal.H never reads sp_buckle_span");
- `unit/retired_keys/test` (:41).

`face_removal`, `face_heat`, `spall_event`, `removal_increment`, `regime_low_high_power`, `pad_clear_ledger` and `sp_melt_skeleton` also mention `Removal.H` (some may be comments only).

**Order this implies.**
1. Move the helpers and the prologue (zero-behaviour, provable on the gate).
2. Re-key the gate's column-only decks onto the list, or accept their retirement, before the 2639 identity can be restated on the face list.
3. Port.
4. Delete.
