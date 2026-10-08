# ROADMAP.md

This is the short project map. Claude should read this at startup, then read
`ACTIVE_STEP.md` for the detailed coding packet. The full historical plan remains
in `in-main-tex-you-will-quizzical-treasure.md`.

## ⚠ Cap and debt

`CLAUDE.md` budgets this file by **section, not total** (rule set 2026-10-07, replacing the
300 → 400 total cap): the **generated** `## Current Status` list (112 lines at 97 archived
steps) is excluded; **hand-written content ≤ 350 lines** (287 today — check:
`awk '/^## Current Status/{g=1} /^Active step:/{g=0} !g{n++} END{print n}' ROADMAP.md`);
`Active step:` + `Next milestones` ≤ 100 lines (75 today); standing rules 10–15 bullets; notes
≤ 40 × 3 lines; total alarm **600**. ⚠ Release rule exercised 2026-10-07 — 11 notes bullets
plus `## Context Files` (which duplicated `CLAUDE.md`'s own list) went to `ARCHIVE_DONE.md`
**verbatim** first, heading *"ROADMAP notes released 2026-10-07"*, and the index and status
list were regenerated in the same pass.

## Standing rules

⚠⚠ **Programme invariants, NEVER archived. Hard cap: 15 bullets.** These are not step
caveats — they are traps that will fire again on steps not yet named, so the release rule
governing `## Known Stale/Important Notes` does not apply to them. ⚠ If this block needs a
16th bullet, one of these has stopped being an invariant; retire it explicitly.

1. ⚠⚠ **A FLAW SIZE AND THE DRIVE IT WAS FITTED AGAINST ARE ONE CALIBRATION, AND SEPARATING THEM
   IS SILENT.** Three instances: the 80 µm refit alone multiplies the column rate by **2.46×**;
   `spall.support_rule`'s 80 µm fired **240 K early** against full stress (retired, P2-2a); and the
   Kant harness carried a **×10** refit for a whole campaign, reading **520.68 K** once its partner
   was deleted. ⚠ `differential` is **gone** — the rule is about the *pairing*, not that form. Any
   step changing the drive must **refit the flaw against Kant in the same step**, and ⚠⚠ **only a
   run can retire a key: ALAMO does not abort on an unused top-level key.**
2. ⚠⚠ **Never score hole diameter at 100 mm depth** — that is the pause underream (Meier
   p. 223). The calibrated section is **150–450 mm**.
3. ⚠⚠ **Quarter domain.** Event-log volumes are one quadrant; **×4** for any Meier volume.
4. ⚠⚠ **The event log's `k_top` is the top BEFORE that removal.** The current top is
   `k_top − n_voided`. Reading it as the post-removal state makes recently-fired columns
   read one cell too tall, and has already made a published figure wrong.
5. ⚠⚠ **THE LZ TRAP: imported analyzers hard-code `LZ = 0.4` at module level** and will
   silently mis-reduce a run on a different domain height. Check it in any reduction that
   imports an older analyzer.
6. ⚠⚠ **`surface_patch.nozzle_collision_radius` defaults to `radius` = 0.2 m under
   `prescribed`** — effectively no guard. Set it explicitly on any prescribed-feed leg.
7. ⚠⚠ **READING A LOGGED OUTPUT: FOUR TRAPS, and the planner walked into (i) and (ii) this week.**
   (i) **NEVER SUM ANY per-step `thermo` column for a total.** `amr.thermo.plot_int` is **not 1** —
   it is **6** on the Meier legs — so a sum under-counts by that factor. Integrate from the **event
   log**, or read a **state** variable (a cumulative column). Recover the factor from the log when
   in doubt: `sum(lat_bank) x 6.03 = 39,795 firings`. This put a 15.2 % leak at 92 % and produced an
   inferred second bug that did not exist (10-06 §11).
   (ii) **CHECK WHETHER THE QUANTITY IS DIMENSIONLESS.** `ledger_err_solid` is divided by
   `ledger_E_src` (`MMWSpalling.H:4587`), so 9e-3 is **not** 9 mJ — it was **7.6 kJ, 0.89 % of the
   budget**. And a *constant relative* error over a *linearly growing* source **is** a steady leak,
   so "it saturates" inverts (10-06 §10).
   (iii) **The event log stamps a step's START and `thermo` its END**, so matching the two by equal
   timestamp pairs unrelated steps. Aggregates are unaffected; per-step labelling needs a one-step
   offset.
   (iv) **`thermo.dat` is written with `%g` — six significant digits.** Compare the **printed**
   value, not a recomputed one, or a pass turns red on formatting alone.
8. ⚠⚠ **A NUMBER IS BACKED ONLY IF BOTH THE FILE IT WAS READ FROM AND THE CODE THAT REDUCED
   IT ARE IN THE REPOSITORY.** A scratchpad log is working material; a **cited** artefact is
   evidence. P1g-0b cited three files that were never committed — seven citations resolving
   to nothing — and that is how a hard-coded `T_gas` = 1200 K against a logged 1604–1634 K
   survived into two interpretations and reversed a conclusion.
   ⚠⚠ **A LEG'S `argv` IS PART OF ITS PROVENANCE, NOT ONLY ITS BINARY HASH.** `PROVENANCE.txt`
   carries hashes and compilation stamps; `metadata` is a *resolved* dump including every
   default, so rebuilding a command line from it is **reconstruction, not replay**. P1g-2c's
   `(c1b)` was unscoreable for exactly this, and Leg B was unattributable in P1g-2b. ⇒ **record
   the argv beside the hash, or the comparison is not available later.**
9. ⚠ **RE-VALIDATING AN ANCHOR: archive first, then re-run — never `--check`.** `--check`
   re-reduces output written by the **old** binary, so it cannot validate a source change;
   re-run the cheapest leg and compare **bytes**. And archive the output directory first
   whenever a harness `rmtree`s its own legs and their hashes live in an earlier manifest —
   the toolchain is **not byte-reproducible**, so that provenance is otherwise
   unrecoverable, and never `rm` what a manifest cites.
10. ⚠⚠ **A Meier number is quotable only from a leg that passed its own gates and whose
    results file says so.** Nothing else is quotable, however often it is repeated. And a
    bore is reported **per azimuth sector**, never as a single diameter — the maximum
    radius over all azimuths picks the grid diagonal, which produced one retracted claim.
    *(Which legs currently qualify is a step caveat and lives in the notes, not here.)*
11. ⚠ **Every long leg under `caffeinate -dis`.** Runs are not checkpointed; host sleep
    loses them. Harnesses assert it themselves.
12. **[standing INSTRUCTION, not a trap — kept here for visibility]** Parked at the
    user's instruction: AMR, the third mesh, the 1 mm run. Remove when they un-park it.
13. ⚠⚠ **THE ROCK IS GRAIN-FREE IN EVERY INPUT.** `is_grain_boundary` reads 0 on all cells
    (1,080,000 checked), so `weibull.a0_gb` is **never drawn** and the max flaw anywhere is
    **8.99 µm**. The 8 µm calibration was made on this same rock, so the threshold is
    self-consistent — nothing is firing 170 K too hot. ⚠ **No packet may claim
    grain-boundary flaws lead the firing**, and enabling Voronoi multiplies the drive by
    **√5 = 2.24×** and reopens `a0` together with the drive fraction (rule 1).
14. ⚠⚠ **ONE CRITERION, ONE KERNEL — AND A COMMENT IS NOT A BAR.** Three times now a decision
    site has disagreed with the flux site it is supposed to mirror, each documented in a nearby
    comment and each shipped anyway: **Pa vs MPa** (`K_Ic` unconverted), **the Newton length**
    (`d_n` vs `DX[0]`), **the coefficient cap** (unconditional vs conditional). ⚠ When a quantity
    has several consumers, the bar is a **source assertion that they resolve it identically**,
    with every by-configuration equality **named and its condition stated**.
    ⚠⚠ **AND CHECKING THE SITE YOU WERE POINTED AT IS NOT CHECKING THE CLAIM.** 2026-10-07: the
    planner wrote "`weibull.m` is inert, and `MMWSpalling.H:8379` proves it". The implementer
    checked that line, confirmed the volume factor is exactly 1 — and inherited the conclusion.
    But `m` appears **twice** on that line, and the second use (`(-log u)^(1/m)`, the draw's own
    shape) does **not** cancel, so `m` had always set the flaw scatter. The pointer was right and
    the claim was false. ⇒ **resolve the claim against the code, not the citation.**
15. ⚠⚠ **A QUANTITY BOUNDED BY THE MECHANISM'S OWN ACTION CANNOT BE ITS DISCRIMINATOR.**
    `lat_Tface_max` is pinned at the firing threshold on any leg that fires, because removal
    consumes every face above it — so a window drawn around it passed by **1.8 K** on a
    window that had itself moved **1,253 K**. Before writing a bar, ask whether the
    mechanism's action bounds the quantity.

⚠ **Seeded 2026-10-02 during the mechanical document pass**, from invariants already in
daily use plus two (5 and 6) that lost prominence in the `Current Status` collapse. The
full triage of the 339-bullet notes section is the **two-reader pass** and is NOT done.

## Current Status

Completed and passing — **one bullet per archive record**, generated so two different
steps can never be merged under a shared label.
⚠ Regenerate with `python3 docs/project/tools/roadmap_status.py` after any archive
  append, and run `archive_index.py` first so the line numbers agree.

- Step 1 - Skeleton Integrator And Build Registration — `ARCHIVE_DONE.md` L147
- Step 2 - Enthalpy Heat Formulation — `ARCHIVE_DONE.md` L159
- Step 3 - MMW Beam Source — `ARCHIVE_DONE.md` L173
- Step 4 - Surface Losses — `ARCHIVE_DONE.md` L184
- Step 4b - Zhang/Oglesby Granite Heating Validation — `ARCHIVE_DONE.md` L195
- Step 5 - Voronoi Microstructure With Mineral Phases — `ARCHIVE_DONE.md` L220
- Step 6 - Heterogeneous And Damage-Modified Conductivity — `ARCHIVE_DONE.md` L265
- Step 6b - H/T Consistency For Heterogeneous Microstructure — `ARCHIVE_DONE.md` L309
- Step 6c - Grain Topology Correction — `ARCHIVE_DONE.md` L349
- Step 6d - Hu Prescribed-Surface-Temperature Conduction — `ARCHIVE_DONE.md` L392
- Step 7 - Heterogeneous Thermoelastic Mechanics — `ARCHIVE_DONE.md` L448
- Step 7b - Hu Thermoelastic Stress Validation, Single-Level — `ARCHIVE_DONE.md` L506
- Step 7c - Early AMR Correctness For Microstructure Regrid — `ARCHIVE_DONE.md` L581
- Step 7d - Hu Thermoelastic Stress Validation With AMR — `ARCHIVE_DONE.md` L630
- Step 8 - Drucker-Prager Failure Criterion And Continuous Damage Evolution — `ARCHIVE_DONE.md` L699
- Step 8b - Hu Breakage-Probability Indicator Without Damage Evolution — `ARCHIVE_DONE.md` L770
- Step 9 - Grain-Boundary Cohesive-Zone Unit Test — `ARCHIVE_DONE.md` L844
- Step 10 - Surface Advancement And Spall Detachment — `ARCHIVE_DONE.md` L898
- Step 11 - Vaporisation Removal And Unified RoP — `ARCHIVE_DONE.md` L963
- Step 11b - Hu LRST And Onset-Spallation With DP Damage — `ARCHIVE_DONE.md` L1031
- Step 12 - Spallation-Model Toggle + Sp Onset Criterion — `ARCHIVE_DONE.md` L1106
- Step 14 - Hu Flame-Jet End-to-End Validation (damage_law) — `ARCHIVE_DONE.md` L1287
- Step 14b - Damage-Law Spall Sampling Knob — `ARCHIVE_DONE.md` L1368
- Step 15a - Convective Robin Surface-Patch BC — `ARCHIVE_DONE.md` L1442
- Step 15b - Sp Mechanics: Weibull Flaw Distribution + Sign Fix + Depth-Resolved Sampling — `ARCHIVE_DONE.md` L1531
- Step 15b - Kant Onset Validation (sp_kant_onset) — `ARCHIVE_DONE.md` L1646
- Step 15c - Surface labeller, detach_mode A/B, AMR regrid-repair, sensitivity sweep (v2 dual-scoring) — `ARCHIVE_DONE.md` L1735
- Step 16 - Sp-vs-depth h_spall mechanics (Rossi validation deferred to 16b) — `ARCHIVE_DONE.md` L1859
- Step 16b - Rossi crack-event count + sub-cell bisection (depth bands not met at working mesh) — `ARCHIVE_DONE.md` L1961
- Step 16-series Rossi chain — closed, R1 deferred to Step 18+ scope — `ARCHIVE_DONE.md` L2068
- Step 16d - T-derived σ_xx with prescribed-T face BC (option (a)) — `ARCHIVE_DONE.md` L2108
- Step 17 - per-cell v_n with simultaneous spall + vap (mechanics half) — `ARCHIVE_DONE.md` L2218
- Step 16c - h_col CSV log + linear cell-to-cell σ_xx (R1 FAILs on (L) artefact) — `ARCHIVE_DONE.md` L2341
- Step 18 - Kant 2017 confining-pressure sweep validation (revised-approach Stage 5) — `ARCHIVE_DONE.md` L2496
- Step 19 - sp_weibull spall under MMW beam (relax Step 16d surface_patch gate) — `ARCHIVE_DONE.md` L2601
- Refactor R1 — extract UpdateRemovalAfterCohesive into a partial header (completed) — `ARCHIVE_DONE.md` L2669
- P1 — Melt-aware spallation (completed, review-accepted) — `ARCHIVE_DONE.md` L2722
- P4 — Cuttings/debris shielding (WITHDRAWN before implementation, 2026-09-15) — `ARCHIVE_DONE.md` L2789
- E1 — Beam/void gap energy leak fix + deposition closure counter (completed, review-accepted) — `ARCHIVE_DONE.md` L2812
- S1 — 1-D surface-resolution study (fixed-flux and Robin sources) (completed, review-accepted) — `ARCHIVE_DONE.md` L3243
- A1 — Removal fixes (scalar flaw, removal-event log) + S1b firing temperature under production Weibull (completed, review-accepted) — `ARCHIVE_DONE.md` L3810
- C1 — Pinned-surface Robin closure (`robin_form = pinned`) (completed, review-accepted) — `ARCHIVE_DONE.md` L4373
- D2a — Impinging-jet face source: spatial Robin h(r, s), T_gas(s) with a prescribed nozzle path (completed, review-accepted) — `ARCHIVE_DONE.md` L4962
- D2a2 — Energy-conserving wall-jet gas temperature T_gas(r, s); no jet flux above the nozzle plane (completed, review-accepted) — `ARCHIVE_DONE.md` L5601
- D2b — Burner-on-feet descent and the scored Meier ROP / hole-diameter test (completed, review-accepted) — `ARCHIVE_DONE.md` L6341
- D2c — Steady-state and no-funnel closures + pre-registration (code only; completed, review-accepted) — `ARCHIVE_DONE.md` L7286
- D2d — Meier campaign: scored run, far-field control, T_nozzle sweep, funnel isolation, full-domain mesh check (completed, review-accepted) — `ARCHIVE_DONE.md` L8251
- D2e — Speed-up acceptance and the mesh-seed diagnosis (completed; NOT /verify-ed) — `ARCHIVE_DONE.md` L8925
- D2f — Contact support rule as Meier ran it: q = 1.0 mesh gate and quantile sensitivity (completed, review-accepted) — `ARCHIVE_DONE.md` L9566
- D2e — retroactive review (2026-09-19; attaches to the D2e entry above) — `ARCHIVE_DONE.md` L10062
- D2g-1 — Long mesh check of the D2f support rule (q = 0.9) to the 0.40 m bottom stop (completed, review-accepted; GATE FAILED) — `ARCHIVE_DONE.md` L10167
- D2h — Key-only diagnostic pairs on the D2f support rule (2026-09-20; reviewer diagnostic, no packet) — `ARCHIVE_DONE.md` L10480
- D2i — Support-rule mesh gate with the idle clock off (COMPLETED 2026-09-21, review-accepted; GATE FAILED) — `ARCHIVE_DONE.md` L10517
- D2j-0 — The hot-disk edge test with today's model (COMPLETED 2026-09-21; its verdict was WITHDRAWN the same evening by D2j-0b) — `ARCHIVE_DONE.md` L10597
- D2j-0b — The hot disk under a descending exclusion plane (reviewer diagnostic, 2026-09-21 evening; no packet) — `ARCHIVE_DONE.md` L10649
- D2k — Side-face flux on exposed vertical faces (COMPLETED 2026-09-22, review-accepted; arbiter PASS, Meier gate FAIL and inconclusive) — `ARCHIVE_DONE.md` L10717
- D2l — Gate 0 scan: why the model drills 2.7× too fast (COMPLETED 2026-09-22; window EMPTY ⇒ pre-registered STOP, no code written) — `ARCHIVE_DONE.md` L10810
- D2m — Is the cut diameter set by the feet or by the flame? (COMPLETED 2026-09-22; MIXED verdict, reported as mixed) — `ARCHIVE_DONE.md` L10885
- D2o-0 — Prescribed floor/wall heating: does a cylindrical hole appear? (COMPLETED 2026-09-22; review-accepted; NEGATIVE, and the field had a defect found afterwards) — `ARCHIVE_DONE.md` L10963
- D2o-0b — Prescribed floor/wall heating with the field bounded in radius (COMPLETED 2026-09-22; REJECTED by review as incomplete; Q1 FAILS; **closes the prescribed-field line**) — `ARCHIVE_DONE.md` L11180
- D2p-0 — Lateral recession at exposed step faces (CLOSED 2026-09-22 at DESIGN TIME, zero run cost; premise refuted, pre-registered fallback triggered) — `ARCHIVE_DONE.md` L11338
- D2p-1 — Corner regularisation (CLOSED 2026-09-22 at DESIGN TIME, never hashed, never launched; refuted by arithmetic) — `ARCHIVE_DONE.md` L11449
- D2q-1 — Make `removed_mf` authoritative (per-cell void, stage 1: pure refactor) (COMPLETED 2026-09-22, review-accepted; GATE PASSED 2639/2639) — `ARCHIVE_DONE.md` L11560
- D2q-2b (checkpoint stage) — Per-cell removal: the mechanism, built and gated (2026-09-22; arbiter deliberately not launched) — `ARCHIVE_DONE.md` L11658
- D2q-2c — The arbiter: does the corner drain spall the wall? (2026-09-22/23, review-accepted) — `ARCHIVE_DONE.md` L11745
- D2r-0 (pre-flight stage) — closed at item 4, PLANNER ERROR, no code written (2026-09-23) — `ARCHIVE_DONE.md` L11890
- D2r-0b — SUPERSEDED BEFORE LAUNCH by a code fact found in review (2026-09-23) — `ARCHIVE_DONE.md` L12012
- D2r-0c — SUPERSEDED BEFORE LAUNCH, inert by construction (2026-09-23) — `ARCHIVE_DONE.md` L12096
- D2r-1 — Heating the SIDE faces of above-plane columns (COMPLETED 2026-09-23, review-accepted; (a) PASSES) — `ARCHIVE_DONE.md` L12169
- D2t — Retire the mask edge with a smooth radial taper, and score SHAPE (COMPLETED 2026-09-23; (a) PASSES) — `ARCHIVE_DONE.md` L12352
- D2u — Meier shape packet: WITHDRAWN BEFORE LAUNCH 2026-09-23 (fatal design flaw + a prediction refuted by data on disk) — `ARCHIVE_DONE.md` L12454
- D2v — the body-in-the-hole h map (COMPLETED 2026-09-24; (a) PASSES, but the ATTRIBUTION IS CONFOUNDED — planner error) — `ARCHIVE_DONE.md` L12551
- D2w — the 2x2 attribution: the MAP owns the funnel, the SIDE PATH owns the rate (COMPLETED 2026-09-24, review-accepted) — `ARCHIVE_DONE.md` L12635
- D2x — split `f`: the derived area is measured, and it does NOT reach Meier's hole (COMPLETED 2026-09-24) — `ARCHIVE_DONE.md` L12771
- D2y — conservative surface area: PARKED BEFORE LAUNCH 2026-09-24 (never ran, no code) — `ARCHIVE_DONE.md` L12860
- D2z — the step-height closure: PLANNED AND PARKED 2026-09-25, pending D2z-0 (never ran, no code) — `ARCHIVE_DONE.md` L12924
- D2z-0 — retire the half-cell removal increment (COMPLETED 2026-09-25; P1/P2 refuted by the METRIC, P3/P4 held) — `ARCHIVE_DONE.md` L12993
- D2z-1 — whole-cell removal unconditional, the key deleted (COMPLETED 2026-09-27, review-accepted) — `ARCHIVE_DONE.md` L13054
- P4 REFUTED, in the informative direction — `ARCHIVE_DONE.md` L13220
- P3a REFUTED — the planner's crater-wall assumption — `ARCHIVE_DONE.md` L13420
- P3 REFUTED — and the mechanism is the finding: THE SIDE TERM DEFENDS ITSELF — `ARCHIVE_DONE.md` L13585
- P7 CONFIRMED — `PD` jammed, and it was pre-registered — `ARCHIVE_DONE.md` L13618
- P0b — ONE surface temperature for the BC, the criterion and the removal (SOURCE, revision 7) — review-accepted 2026-09-30 — `ARCHIVE_DONE.md` L14012
- P1 — the gradient-aware firing criterion (SOURCE, revision 5) — review-accepted 2026-10-01 — `ARCHIVE_DONE.md` L14142
- P1b — unsupported hot rock comes off (SOURCE, revision 2) — implemented 2026-10-01, **(c) FAILED on a grid artefact** — `ARCHIVE_DONE.md` L14286
- P1c — a grid-isotropic support rule (revision 2) — implemented 2026-10-01, **FAILED both primary gates; the cause was found and P1d is WITHDRAWN** — `ARCHIVE_DONE.md` L14401
- P1e — the staircase side-face correction — implemented 2026-10-01, **correction RIGHT, stall CURED, both shape gates FAILED; the attribution then CONFIRMED by field diagnostic** — `ARCHIVE_DONE.md` L14489
- P1f — lateral flake drive + lateral removal (revision 11) — implemented 2026-10-02, **ACCEPTED, every gate green; the mode is CORRECT and nearly INCONSEQUENTIAL in scope** — `ARCHIVE_DONE.md` L15875
- P1g-0 — the corner-face temperature fix and the azimuthal witness — implemented 2026-10-02, **ACCEPTED; mechanism CONFIRMED, deep deficit only PARTLY closed, and a SECOND mismatch found** — `ARCHIVE_DONE.md` L16108
- P1g-0b — one criterion, one kernel: the coefficient cap and the four-site argument lock — implemented 2026-10-02, **ACCEPTED (revision 2); the third instance of one bug family closed, and (d1) REVERSED** — `ARCHIVE_DONE.md` L16407
- P1g-1 — the lateral flake sink (sub-cell lateral recession) — implemented 2026-10-03, **ACCEPTED AS BOOKKEEPING; PHYSICALLY INERT — the recession RATE is identical to whole-cell** — `ARCHIVE_DONE.md` L17919
- P1g-2 — the lateral flake sink, THERMAL half — implemented 2026-10-03, **ACCEPTED; the thermal half IS real and attributable — it delays the ONSET, and the long-run recession converges back toward whole-cell** — `ARCHIVE_DONE.md` L18190
- P1g-3 — the flake thickness from the rock — implemented 2026-10-05, **REJECTED ON PROVENANCE; the physics findings are reviewer-verified and load-bearing. Closure carried as ord 0 of P2-0** — `ARCHIVE_DONE.md` L18448
- P2-0 — the wall heat-transfer coefficient from the annulus flow — ord 0/1 run 2026-10-06, **SUPERSEDED: the premise was a duct law applied where the flow is a jet. ord 0 and ord 1 are KEPT; the source work was never built** — `ARCHIVE_DONE.md` L18509
- P1g-2b — the lateral sub-cell mode at Meier scale — run 2026-10-06, **ACCEPTED with `(b3)` FAILED on the thermal leg; three measurements nobody had made, and the widening attributed to the right cause at last** — `ARCHIVE_DONE.md` L18590
- P1g-2c — the thermal shift's energy leak, CLOSED; P3 refuted (review-accepted 10-06) — `ARCHIVE_DONE.md` L18657
- P2-1 — the wall heat-transfer coefficient, DERIVED; the typed constant was too LOW, and the wall self-limits (review-accepted 10-07) — `ARCHIVE_DONE.md` L18802
- P2-2a — the stress-free reference, and the retirement of `differential` and the support rule (2026-10-08, rev 7, ACCEPTED — all 11 criteria PASS) — `ARCHIVE_DONE.md` L18954
- D2z-2 — witness repair, derived conductance, Meier (COMPLETED 2026-09-28) — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)
- D2aa — heat delivery under the burner body (COMPLETED 2026-09-28) — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)
- D2ab — the face assignment and the gap-derived coefficient (COMPLETED — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)
- D2ac — Meier at Meier's kinematics: THE FLOOR-FLUX EXCESS (COMPLETED — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)
- D2ad — calibrate the floor flux on Fig. 8.8's standoff: THE WINDOW IS EMPTY — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)
- P0 — a real surface temperature, and a mesh-independent anchor (COMPLETED — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)
- D2s — WITHDRAWN BEFORE LAUNCH (2026-09-23); premise and arithmetic both — `ARCHIVE_DONE.md` L14575 (no record of its own; in the verbatim Current Status)

Active step:
- **P2-2b rev 2 — the measured Grimsel properties and the relief factor. READY.** The *rock* half of
  the 2.27x, on the stress expression P2-2a proved. rev 1 was **BLOCKED at `(a1)`, correctly**; ord 0
  and the rest of ord 1 PASSED and are carried. (`ACTIVE_STEP.md`;
  **`Claude_markdowns/2026-10-08.md`**.)
  ⚠⚠ **`(a1)` WAS A PLANNER DEFECT TWICE OVER.** It enumerated hard-coded `rho*c_p` by grepping the
  **product**, found five, and barred the row at *"five and no sixth"* — but two sites read `rho` and
  `Cp` as **separate scalars** and cannot match a product pattern, so the bar could only pass by
  repeating the planner's own search error. ⇒ **the whitelist is abolished at the user's decision:
  no site in the rock path reads `rho` or `Cp` as a per-phase constant, zero exceptions, so there is
  no list to keep correct.** Each site resolves either by reading the per-cell field or by building
  from the **law**; the style is to read the two **separately** at the point of use and never store a
  pre-multiplied product (`rev_rhoCp`/`pin_rhoCp` retire). ⚠ No site is in the hot loop, so the speed
  argument for keeping any fixed is worth nothing.
  ⚠⚠ **THE PROPERTY CHANGE REVERSES SIGN AT THE COLD END** — measured: `rho*c_p` is **0.9250** of
  today's at 293 K, crosses 1.0 at **325.7 K**, reaches **1.4899** at 1200 K. "~1.3x" is never a
  scaling, it is the integral over the span the removed rock spans ⇒ `(d3)` reports the factors apart.
  ⚠⚠ **rev 1's P1 is WITHDRAWN**: 0.930 -> ~1.25 needs x1.344 but its own factors give **1.330-1.451**
  — unreachable from itself. Only the enthalpy factor is measured (**1.298** at the firing
  temperature); the mean-T factor is now **measured, not predicted**, and no endpoint is asserted.
  ⚠ **The latent-heat hazard is real in code and INERT on both decks** (`vapor.enabled = 0`, no
  phase-change scalar set) ⇒ `(b3)` cannot see it; new `(b6)` switches melting **on** as its only
  witness. What `(b3)` is really about is `(a2)`'s **-204.70 K** divergence between the two enthalpy
  curves at equal enthalpy.
  ⚠⚠ **`lambda(600 degC)/lambda(20 degC) = 0.501`** ⇒ the flux skin thins **1.4 -> 0.70 mm**, **0.35 of
  a 2 mm cell**, below the surface reconstruction's design range. **Every 2 mm rate here is
  mesh-limited and provisional**; pre-registered that the 2 mm vs 1 mm gap **GROWS**. ⚠ 1 mm is
  parked (rule 12) — un-parking is the user's call and is not assumed.
  ⚠⚠ **FIT TARGET SETTLED (user, 10-08): option A, `R` fitted so Kant fires at 850 K.** rev 2
  contradicted itself — it fitted to Kant's **measured** 823-873 K while also requiring Meier's `T_s`
  *unchanged*. Today both decks fire at **~757 K**, so the fit **raises `T_s` by ~93 K on both**, and
  that rise **is** the mechanism: hotter rock leaves, hotter rock cost more, the price rises to a
  predicted **~1.37 GJ/m3** (x1.47 on 0.930, top of the independent 1.33-1.54 band) against Meier's
  **1.40**. ⚠ 850 K is a deliberate **under-fit** — 823 (plane) and 873 (surface) are two different
  measurements, not an error bar, and the criterion fires on a **surface** temperature — so `R` comes
  out slightly too large and a result that closes here closes harder at 873. ⚠⚠ The **~2 % residual
  is suspiciously good and is flagged**: it stands only with `R` inside `(c2)`'s band and `(c7)`'s two
  expansion fits agreeing. ⚠ `K_Ic(T)` and `R` are **exactly degenerate** on Kant's single firing
  temperature ⇒ `(c2)`'s band needs `kic.*` **pinned in the bar**.

Next milestones after the active step:
⚠⚠ **FIVE ABSORPTION PATHS STAND MEASURED** — firing rate, energy budget, wall surface temperature,
wall coefficient, standoff (the last only under `nozzle_descent = prescribed`, a modelling choice).
A 15.3 % leak was worth **0.16 mm** of bore; a **x1.430** rise **+9.2 %** of power. ⚠ **An open
research thread, deliberately NOT a standing rule** (§9.5). ⚠⚠ **The wall is closed as the width
lever, five ways**; P2-0's premise is **inverted** — typed 169.3 too low near the floor (`ARCHIVE_DONE.md` *"P1g-2b"*/*"P2-0"*/*"P2-1"*).

- **P2-2b — properties and the relief factor.** The active step; see above.
- **P2-2d — THE FLUX EXPONENT.** The Weibull size term out of `InitializeFlaws` to `ell_q^3` behind
  `size_volume`, swept over `m` = 20/10/5/4, against Kant's measured **0.229-0.306**. ⚠⚠ Genuinely
  open: either the floor tolerates a 1.5-1.8x onset spread or **no variate choice helps** — the
  scatter follows the **exponent** demanded, not the variate (`m' = 2m`, §9.3). The frozen-floor
  cliff is **measured**, not predicted; the 48 % datum measured a bug P1e repaired.
  ⚠⚠ `P2-2c` is **RETIRED** (distribute strength) and the name is **not reused**.
- **P2-2e — THE COMPRESSIVE FAILURE CRITERION.** Buckling and peeling of a compressed layer
  (Rauenzahn; Evans & Hutchinson), which the relief factor stands in for. ⚠ Its case must be argued
  from P2-2b's data — if `R` lands near `UCS/sigma` with a grain-scale flaw, that is the evidence.
  ⚠ `differential` is **gone**, so this builds from the literature. A temperature-dependent `R`
  belongs here (microcracking grows above 300 degC).
- **P2-3 — the floor: conductivity's consequences and total delivered power.** ⚠⚠ **The prescribed
  feed has masked about half the delivery error for two weeks**: D2ac measured **1.98-2.23x** the
  required floor flux at Meier's ~100 mm standoff against only **1.37-1.47x** at our realised
  179-224 mm; P2-2's `(d5)` turns that into a number. ⚠ Measured, **not** what the archive says: the
  floor is **flat 0.48-0.60 MW/m2** to 35 mm, **0.78-0.88** at 38-46, **2.40** at 48; recession
  **2.4-2.7 m/h** to 40 then **3.5** at 40-50 — **not a cone**. The rim's high flux-per-footprint is
  an **area** effect (~33 mm exposed wall per rim column under `area_law = wetted`), not a
  coefficient one. ⚠ D2ad binds: halving the level moved the axis **-5 %**, the bore edge **-89 %**,
  three legs aborted — **a level cut is not a translation**; and its "volume already right,
  redistribute at constant power" **must not be carried forward** (today's volume is **2.27x**).
- **P2-4 — the floor taper's four radii, and the mesh gap, which are ONE question.** ⚠ Provenance
  clean: 40 mm is `foot_r_outer` (Oe80 skirt, VT5 drawings), 60 mm the measured mouth ring, 56/76 are
  40+16 and 60+16 where **16 mm is one numerical choice used twice** (8 cells at 2 mm). **The bore is
  quotable.** ⚠⚠ But the shallow bore has a **mesh-tied structural floor**: the blend hands over at a
  centre of 48 mm, so this configuration **cannot cut narrower than Oe ~96-100 mm at 10-50 mm depth**
  — D2v pre-registered it, unread for two weeks. We read **118.0**; Meier wants **91.5**; at 1 mm the
  floor drops to **~88-92 mm**. ⚠ **Write the 1 mm prediction down first.** (1 mm parked, rule 12.)
- **P1g-5 — Leg A, the shape leg.** Roundness by sector plus the Meier Oe(z) profile. ⚠ Today's
  1.003-1.019 was measured under the 8 um flaw and the 1.5 mm cap. **P1h — overhang collapse
  (regime 7)**: breaks `n_voided == 1` by design, 100-150 cells/event.
- ⚠ **promote `lateral_thermal` to default, retire the whole-cell lateral path** — re-pins every anchor.
- ⚠ **P2-1's four named risks**, carried: the **uncounted fallbacks** to the typed 169.3
  (`MMWSpalling.H:5012-5022`, p50 171.89 indistinguishable from 169.3 in a 5 W bin); gates set at a
  **representative bore** not bounding the run's live 37.5-48.7 mm; the scorer's **unauthorised
  narrowing to three depth bands**; the **march's missing `side_face_h_max` cap**. ⚠ Plus the
  **8.1 % kernel-vs-march side-power split**, diagnosed and unfixed.
- ⚠⚠ **MMW — THE TARGET APPLICATION, AND SPALLATION HAS NEVER RUN WITH THE BEAM ON IN A VALIDATED
  CASE.** Leg B and the Kant anchors set `beam.P0 = 0`; Zhang has the beam on with thermal stress
  **inert**. ⚠ P2-2's flux exponent decides whether the model extrapolates there at all: under a beam
  the heat arrives volumetrically with its peak **below** the surface, where none of the five
  absorption paths is guaranteed to hold.
- ⚠ **sources still absent as primary texts** — Hoeser 2016, Rossi 2020, Shen 2018, Wilkinson 1993;
  the Grimsel reference stands in, explicit about its `[estimate]` and `[recalled]` rows.

## Validation Story

- **Unit/analytical** — individual numerical components. **Zhang/Oglesby** — the MMW source and
  thermal model. **Hu** — the prescribed-T chain on `damage_law`.
- **Kant 2017 Central Aare** — the LEFM Sp onset criterion on `sp_weibull`. ⚠⚠ **From P2-2 it is
  also the CALIBRATION SET**: `R` and the flux exponent are fitted here, on
  `sp_kant_onset_meiercfg`, so that Meier stays a pure test set.
- **Meier 2017 pilot** — ROP, hole diameter and volume. ⚠ Live error: **+16.2 mm** depth-matched bore
  and **2.27x** volume (`ARCHIVE_DONE.md` *"P1g-2c"*). ⚠ Quotable only under standing rule 10.
- ⚠ **Rossi 2018** is a prescribed-T mechanics regression of the Sp depth scan, **not** a flake-size
  validation for self-heated drilling (retired 2026-09-15, S1).
- ⚠⚠ **MMW has no validated case with the beam on** — see the last milestone.

## Known Stale/Important Notes

⚠ **Cap: 40 bullets, 3 lines each.** Only caveats that constrain the active or next few steps.
Full prior text: `ARCHIVE_DONE.md`, heading *"ROADMAP `Known Stale/Important Notes` as it stood
before P1g-1"*. ⚠⚠ A compressed bullet is a **pointer**; read the archive before acting on one.

- ⚠⚠ **THE LIVE MEIER ERROR, CORRECTED 2026-10-05 — "+55–65 % volume" IS STALE** and predates the
  D2v funnel fix. Current: **+10–12 mm width at every depth** (Ø 92–100 vs Meier's 84–96) and
  **2× rate** (≈3.2 vs 1.5 m/h); volume per metre **6.3 vs 6.84 L/m, within ~10 %**. ⇒ the rate is
  the larger error and it is **gas-side** (the floor is cap-free); the width is cap-governed and is
  not a prediction. `2026-10-05.md` §2.
- ⚠⚠ **THE LATERAL MODE IS SCOPED OUT OF 44 % OF P1's FREEZE POPULATION.**
  `side_face_h_r_min = 40 mm` splits 7a hole wall (in scope) from **7b floor riser (out of scope**,
  `MMWSpalling.H:2086`/`:2095`). Measured on P1's archive: **152 of 347** frozen columns are inside
  it, including the collider at **r = 26.4 mm**. Report the two populations separately, never pooled.
- ⚠⚠ **NEITHER BREAKING RULE IS FORMULATED FOR VOLUMETRIC HEATING.** `absolute` is local in
  temperature at one point; `differential` needs a thickness. With MMW the peak is interior
  (planner estimate ≈0.15 mm at 300 m⁻¹, ≈0.46 mm at 100 m⁻¹ — both at or beyond the flake), so the
  surface value is not the most-stressed value. The flame-jet literature cannot settle this: those
  experiments are surface-heated.
- ⚠ **REDUCE WHAT A PER-BOX LOOP ACCUMULATED, NEVER WHAT A PER-COLUMN LOOP DID.** Two steps in a row
  booked a per-column sum `rank_count` times (`lat_bank`, then `E_shift`); the second read
  `ledger_err_solid` = **0.0898** against a 5e-3 bar. The energy invariant is what caught both.

**The sub-cell removal group — these four decide P1g-1.**

- ⚠⚠ **`phi` IS NEVER READ BY THE CONDUCTION KERNEL** (`MMWSpalling.H:1699-1700`, absent from the
  `ParallelFor` at `:1701-1835`) ⇒ a receded cell is still a full `dz` of rock thermally, and a naive
  enthalpy subtraction reads `ell` back **1.32× / 3.81× too long** at `dz/ell` = 2 / 1, and **flat** at 0.5.
- ⚠⚠ **PARTIAL RECESSION IS THE CONSISTENT SCHEME, and the machinery exists** —
  `phi_new = phi_old − h_c` (`Removal.H:1475`), with `h_c` forced to `dz` at `:1387`. Voiding a whole
  cell and pushing the remainder's enthalpy to the neighbour **cascades one cell per time step**.
- ⚠⚠ **THE STATE CANNOT EXPRESS LATERAL RECESSION.** Per-column floor from one scalar
  (`phi = base − (k−klo)·dz`, `Removal.H:2818-2821`); a lagging column can only be removed from its
  top, so no sharp floor/wall edge can be rounded. A lateral accumulator must be a **new field**.
- ⚠⚠ **AN UNDERCUT MAKES `surface_mf` WRONG, NOT STALE** — it is still phi-derived, and
  `UpdateSurfaceMaskFromPhi` could not be converted identity-preservingly. Any packet creating one
  runs `surface.follow_mask = 1` or fixes the mask under its own gate; restate beam-closure contiguity.

**The per-cell path.**

- ⚠⚠ **THE FLUX-ANCHORED ESTIMATOR CANCELS THE FLAKE-SCALE CRITERION.** The map **peaks at the pad
  ring** (1.19 MW/m² at 38–42 mm against 0.80 on axis), so anchoring `ell` to flux gives the wall band
  `ell` ≈ 1.1 mm, the cap stops binding and every wall column fires again. **Land one at a time.**
- ⚠⚠ **P1's ROUGHNESS INSTABILITY — the P1g-3 witness.** Aborted t = **177.616 s** at **r = 26.4 mm,
  inside the floor**; **152 of 316** columns stopped firing before 120 s (rev 7: zero), median depth
  **48 vs 180 mm**, **78 %** with ≥ 2 exposed sides. The gradient physics is right, its geometry wrong.
- ⚠ **A NEW REMOVAL REGIME OWES A NEW LEDGER TERM.** Regime 5 once voided cells unrecorded:
  `ledger_err_solid` ran to **−0.258**, pinned at **−8.95 J per cell** against the 8.69 J a 2 mm cell
  holds at 500 K. Any regime added to the winner selection owes the term, **unconditionally**.
- ⚠ **THE THREE `robin_form`s EACH DEFINE "THE SURFACE" DIFFERENTLY, AND NONE IS THE CRITERION'S.**
  `cell` (default) is the raw average; `face` Newton-solves a linear half-cell; `pinned` clamps it.
  Score the form the target configuration actually runs.
- ⚠⚠ **THE KANT ANCHOR RUNS THE WRONG BC.** `sp_kant_onset_meiercfg/test:113-120` never passes
  `surface_patch.robin_form`, so it defaults to `Cell` while Meier runs `surface`. ≈19 of the 21.6 %
  firing-time spread is that BC artefact, **not** the estimator.
- ⚠ **THE WEIBULL ONSET MOVES ~3.7 % WITH MESH BY CONSTRUCTION** — flaws drawn **per cell**
  (`MMWSpalling.H:6589`), `weibull.V0` is `V_cell` at 2 mm only. ⇒ a Weibull leg **cannot** measure
  whether the criterion is mesh-independent; use a deterministic flaw.
- ⚠ **`weibull.a0_*` IS THE WEIBULL SCALE, NOT THE MEAN** — sample mean `a0·Γ(1+1/m)` ≈ `0.965·a0`
  at `m = 15`. Per-cell RNG is a splitmix64 hash of the **global** cell index, bit-identical across
  rank counts and deterministic across AMR regrid.

**Thermal and flux plumbing.**

- ⚠⚠ **LATERAL CONDUCTION IS `k_face`** (`MMWSpalling.H:1537-1571`, harmonic `inv = 0.5/kc + 0.5/k_n`),
  lateral axes only, and **must be symmetric** from either cell or energy is not conserved.
  ⚠ `WallStepDiagnostic` (`:1861-1905`) uses **raw** conductivity and cannot see such a factor.
- ⚠⚠ **CALL ORDER: THE CRITERION RUNS BEFORE THE FLAME IS REBUILT.** `UpdateRemovalAfterCohesive`
  (`:1127`) → `PerCellRemoval` (`:1131`) → swaps → `BuildFlameColumns` (`:1618`) ⇒ anything in
  `Removal.H` reading `flame_h`/`flame_Tg` gets the **previous** step, and `flame_h` is all-zero on step 1.
- ⚠⚠ **`h_expr` IS THE FLOOR `h` UNDER `jet_closure=enthalpy` — NOT inert.** `MMWSpalling.H:2062`/
  `:2067`; floor faces have `s > 0`. Distinct from the D2u blocker, which is **above**-plane only.
- ⚠ **A TAPER IN `h_expr` DOES NOT REMOVE THE HARD MASK** — `r <= surface_patch.radius` is hard-coded
  at `:1579`/`:1626`, so `radius` must sit beyond the taper's outer edge. `h_expr` must land on a
  **positive** floor; `:1788` aborts on `h <= 0`.
- **A CONFINED `h` MAP RAISES THE EXHAUST TEMPERATURE** — less power leaves the block top, so a
  declared above-plane `T_gas` anchored to the frozen **1690 K** sits on the **low** side. State the
  direction; never adjust to match.

**Process and harness.**

- ⚠⚠ **NEVER PUT A STOP RULE FROM AN OPEN-LOOP PREDICTION ON A CLOSED-LOOP LEG.** Meier absorbs a
  rise through two feedbacks the 1-D column lacks (`jet_closure = enthalpy` depletes `T_gas`; `h_expr`
  falls off with a standoff a faster descent opens). **Predict direction on Meier, magnitude on 1-D.**
- ⚠ **THE TREE SWEEP HAS NO RUNNER** — no `run.py --check`/`--repin`. It swaps the binary into place
  and diffs an md5 list against a **declared reference**; the re-pinning procedure must be written
  down before anyone runs it. `witness/{hash.sh,wit.sh,post_all.sh,guarded_post_all.sh}`.
- ⚠ **A PRESCRIBED-FEED LEG CAN RUN INTO THE BOTTOM OF THE BLOCK, AND THE BASE IS ADIABATIC.** Check
  floor-height-above-base at every scored time; the nozzle plane reaches only 283 mm at 800 s, and
  under `prescribed` pads at 28–40 mm **pass through rock** ⇒ Ø and volume void from first contact.
- ⚠⚠ **`q_ratio` DOES NOT EXIST ON THE WALL.** Zero or non-finite on **1521 of 1584** regime-6
  events (96.0 %) in `studies/p1f_lateral/output_legB_fix/FH_removal_events.csv`, against 0 of
  193,322 regime-1. Any wall quantity derived from it is built on 4 % of the population (2026-10-05).
