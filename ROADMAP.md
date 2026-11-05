# ROADMAP.md

This is the short project map. Claude should read this at startup, then read
`ACTIVE_STEP.md` for the detailed coding packet. The full historical plan remains
in `in-main-tex-you-will-quizzical-treasure.md`.

## ⚠ Cap and debt

`CLAUDE.md` budgets this file by **section, not total** (rule set 2026-10-07, replacing the
300 → 400 total cap): the **generated** `## Current Status` list (112 lines at 97 archived
steps) is excluded; **hand-written content ≤ 350 lines** (309 on 2026-10-28 — check:
`awk '/^## Current Status/{g=1} /^Active step:/{g=0} !g{n++} END{print n}' ROADMAP.md`);
`Active step:` + `Next milestones` ≤ 100 lines (58 on 2026-10-28); standing rules 10–20 bullets (raised 2026-10-12); notes
≤ 40 × 3 lines; total alarm **600**. ⚠ Release rule exercised 2026-10-07 — 11 notes bullets
plus `## Context Files` (which duplicated `CLAUDE.md`'s own list) went to `ARCHIVE_DONE.md`
**verbatim** first, heading *"ROADMAP notes released 2026-10-07"*, and the index and status
list were regenerated in the same pass.

## Standing rules

⚠⚠ **Programme invariants, NEVER archived. Hard cap: 20 bullets** (raised 15 -> 20 by the user
2026-10-12). These are not step
caveats — they are traps that will fire again on steps not yet named, so the release rule
governing `## Known Stale/Important Notes` does not apply to them. ⚠ If this block needs a
**21st** bullet, one of these has stopped being an invariant; retire it explicitly.

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
   ⚠⚠ **But a RATIO of two quadrant quantities needs BOTH sides scaled or NEITHER.** P2-2f rev 1
   scaled the volume and not the energy and read the price **0.309** instead of **1.237 GJ/m³** —
   a 4x error produced by obeying the first sentence alone.
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
   inferred second bug that did not exist (10-06 §11). ⚠⚠ **And establish the column's KIND before
   reading it at all — the two errors are OPPOSITES.** `lat_bank_hf` is **per-step** and was read
   at its FINAL value, making a wall charge identically **zero**; `foot_mech_vol` is **cumulative**
   and was **summed**. Obeying "never sum" alone produces the first (P2-2f rev 1, instances 3 and 4).
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
16. ⚠⚠ **A PER-LOCATION VALUE THAT IS ZERO IN A PATTERN IS A DECOMPOSITION SIGNATURE BEFORE
    IT IS PHYSICS.** `ell_q` read zero on every scored onset and both a reviewer and the
    planner diagnosed a mechanism for it — one inventing a "needs recession first" rule not
    in the source. It was positive **iff `col_j <= 11` AND `k_top <= 11`**, 100.00 % of rows
    over four legs: rank 0's box, because the log vector was the one of seven not in the
    `ReduceRealSum` block. Physics cannot depend on the **y index**. Check whether the
    pattern matches the domain split before reasoning about cause. *(Two findings in three
    days; log 2026-10-12.)*
17. ⚠⚠ **BEFORE A NUMBER BECOMES A TARGET, STATE WHAT ITS COLUMN PHYSICALLY MEASURES.** Two
    misattributions in five revisions of one packet: `0.229/0.306` was quoted as *"Kant
    MEASURED"* but is **Wilkinson's**, on **Barre/Westerly**, already the output of a Weibull
    fit — and its citation `[258]` resolved **past the end** of a 151-line file. The planner
    then replaced it with a fit to `theta_s,sens`, which is a **calibration sensor plate that
    replaced the rock** — in the very session that caught the first. ⚠ A citation that cannot
    be opened is not a source, and a suffix (`_sens`) is often the only warning.
    ⚠ Third instance (P2-3, 2026-10-13): the SAME table's `h_fl` was then used as "the measured
    coefficient" of Meier's flame — it is a 3–10 kW threshold burner's ("where at least a single
    spall could be detected", Kant 2016 §1.1); the drilling flame was measured by Höser (7.1 kW/m²K).

18. ⚠⚠ **THE GRIMSEL PROPERTY TABLES STOP AT 873.15 K** and repeat their last value above it, and **about a fifth of
    every calibrating population fires above it** (P2-2f Meier 22.6 %, G1 20.1 %, P2-13 CAL 20.7 %). Report the share
    beside any level a step sets; an offline port moved the mean firing T only ~1 K (parked by the user 2026-10-13;
    the full note: `ARCHIVE_DONE.md` *"ROADMAP notes released 2026-11-02"*).
19. ⚠⚠ **"FIRST FIRING" MUST NAME ITS DIRECTION FILTER.** First top-face (dir 5) and first in any direction differ by
    **7.5 K** on a drilling deck (84 % of G1's columns are first hit from the side, a weakest-of-several selection)
    and not at all on a flat deck (P2-13 DR-2). A packet's refuting-evidence line once compared the two.

20. ⚠⚠ **NAME THE μ: ALONG-LINE OR NORMAL.** On the face list a line runs along a grid axis; its lengths are
    along-line, and a material length (d, a flake thickness, a span) is NORMAL, so the two differ by μ = n·e\*. Twice
    in a row a packet compared them without the μ (P2-14a d/μ, P2-15 x = μ·c_f·ℓ_eff). Every packet row that compares a
    length with a material property states which kind each side is.

⚠ **Seeded 2026-10-02 during the mechanical document pass**, from invariants already in
daily use plus two (5 and 6) that lost prominence in the `Current Status` collapse. The
full triage of the 339-bullet notes section is the **two-reader pass** and is NOT done.

## Current Status

Completed and passing — **one bullet per archive record**, generated so two different
steps can never be merged under a shared label.
⚠ Regenerate with `python3 docs/project/tools/roadmap_status.py` after any archive
  append, and run `archive_index.py` first so the line numbers agree.

- Step 1 - Skeleton Integrator And Build Registration — `ARCHIVE_DONE.md` L196
- Step 2 - Enthalpy Heat Formulation — `ARCHIVE_DONE.md` L208
- Step 3 - MMW Beam Source — `ARCHIVE_DONE.md` L222
- Step 4 - Surface Losses — `ARCHIVE_DONE.md` L233
- Step 4b - Zhang/Oglesby Granite Heating Validation — `ARCHIVE_DONE.md` L244
- Step 5 - Voronoi Microstructure With Mineral Phases — `ARCHIVE_DONE.md` L269
- Step 6 - Heterogeneous And Damage-Modified Conductivity — `ARCHIVE_DONE.md` L314
- Step 6b - H/T Consistency For Heterogeneous Microstructure — `ARCHIVE_DONE.md` L358
- Step 6c - Grain Topology Correction — `ARCHIVE_DONE.md` L398
- Step 6d - Hu Prescribed-Surface-Temperature Conduction — `ARCHIVE_DONE.md` L441
- Step 7 - Heterogeneous Thermoelastic Mechanics — `ARCHIVE_DONE.md` L497
- Step 7b - Hu Thermoelastic Stress Validation, Single-Level — `ARCHIVE_DONE.md` L555
- Step 7c - Early AMR Correctness For Microstructure Regrid — `ARCHIVE_DONE.md` L630
- Step 7d - Hu Thermoelastic Stress Validation With AMR — `ARCHIVE_DONE.md` L679
- Step 8 - Drucker-Prager Failure Criterion And Continuous Damage Evolution — `ARCHIVE_DONE.md` L748
- Step 8b - Hu Breakage-Probability Indicator Without Damage Evolution — `ARCHIVE_DONE.md` L819
- Step 9 - Grain-Boundary Cohesive-Zone Unit Test — `ARCHIVE_DONE.md` L893
- Step 10 - Surface Advancement And Spall Detachment — `ARCHIVE_DONE.md` L947
- Step 11 - Vaporisation Removal And Unified RoP — `ARCHIVE_DONE.md` L1012
- Step 11b - Hu LRST And Onset-Spallation With DP Damage — `ARCHIVE_DONE.md` L1080
- Step 12 - Spallation-Model Toggle + Sp Onset Criterion — `ARCHIVE_DONE.md` L1155
- Step 14 - Hu Flame-Jet End-to-End Validation (damage_law) — `ARCHIVE_DONE.md` L1336
- Step 14b - Damage-Law Spall Sampling Knob — `ARCHIVE_DONE.md` L1417
- Step 15a - Convective Robin Surface-Patch BC — `ARCHIVE_DONE.md` L1491
- Step 15b - Sp Mechanics: Weibull Flaw Distribution + Sign Fix + Depth-Resolved Sampling — `ARCHIVE_DONE.md` L1580
- Step 15b - Kant Onset Validation (sp_kant_onset) — `ARCHIVE_DONE.md` L1695
- Step 15c - Surface labeller, detach_mode A/B, AMR regrid-repair, sensitivity sweep (v2 dual-scoring) — `ARCHIVE_DONE.md` L1784
- Step 16 - Sp-vs-depth h_spall mechanics (Rossi validation deferred to 16b) — `ARCHIVE_DONE.md` L1908
- Step 16b - Rossi crack-event count + sub-cell bisection (depth bands not met at working mesh) — `ARCHIVE_DONE.md` L2010
- Step 16-series Rossi chain — closed, R1 deferred to Step 18+ scope — `ARCHIVE_DONE.md` L2117
- Step 16d - T-derived σ_xx with prescribed-T face BC (option (a)) — `ARCHIVE_DONE.md` L2157
- Step 17 - per-cell v_n with simultaneous spall + vap (mechanics half) — `ARCHIVE_DONE.md` L2267
- Step 16c - h_col CSV log + linear cell-to-cell σ_xx (R1 FAILs on (L) artefact) — `ARCHIVE_DONE.md` L2390
- Step 18 - Kant 2017 confining-pressure sweep validation (revised-approach Stage 5) — `ARCHIVE_DONE.md` L2545
- Step 19 - sp_weibull spall under MMW beam (relax Step 16d surface_patch gate) — `ARCHIVE_DONE.md` L2650
- Refactor R1 — extract UpdateRemovalAfterCohesive into a partial header (completed) — `ARCHIVE_DONE.md` L2718
- P1 — Melt-aware spallation (completed, review-accepted) — `ARCHIVE_DONE.md` L2771
- P4 — Cuttings/debris shielding (WITHDRAWN before implementation, 2026-09-15) — `ARCHIVE_DONE.md` L2838
- E1 — Beam/void gap energy leak fix + deposition closure counter (completed, review-accepted) — `ARCHIVE_DONE.md` L2861
- S1 — 1-D surface-resolution study (fixed-flux and Robin sources) (completed, review-accepted) — `ARCHIVE_DONE.md` L3292
- A1 — Removal fixes (scalar flaw, removal-event log) + S1b firing temperature under production Weibull (completed, review-accepted) — `ARCHIVE_DONE.md` L3859
- C1 — Pinned-surface Robin closure (`robin_form = pinned`) (completed, review-accepted) — `ARCHIVE_DONE.md` L4422
- D2a — Impinging-jet face source: spatial Robin h(r, s), T_gas(s) with a prescribed nozzle path (completed, review-accepted) — `ARCHIVE_DONE.md` L5011
- D2a2 — Energy-conserving wall-jet gas temperature T_gas(r, s); no jet flux above the nozzle plane (completed, review-accepted) — `ARCHIVE_DONE.md` L5650
- D2b — Burner-on-feet descent and the scored Meier ROP / hole-diameter test (completed, review-accepted) — `ARCHIVE_DONE.md` L6390
- D2c — Steady-state and no-funnel closures + pre-registration (code only; completed, review-accepted) — `ARCHIVE_DONE.md` L7335
- D2d — Meier campaign: scored run, far-field control, T_nozzle sweep, funnel isolation, full-domain mesh check (completed, review-accepted) — `ARCHIVE_DONE.md` L8300
- D2e — Speed-up acceptance and the mesh-seed diagnosis (completed; NOT /verify-ed) — `ARCHIVE_DONE.md` L8974
- D2f — Contact support rule as Meier ran it: q = 1.0 mesh gate and quantile sensitivity (completed, review-accepted) — `ARCHIVE_DONE.md` L9615
- D2e — retroactive review (2026-09-19; attaches to the D2e entry above) — `ARCHIVE_DONE.md` L10111
- D2g-1 — Long mesh check of the D2f support rule (q = 0.9) to the 0.40 m bottom stop (completed, review-accepted; GATE FAILED) — `ARCHIVE_DONE.md` L10216
- D2h — Key-only diagnostic pairs on the D2f support rule (2026-09-20; reviewer diagnostic, no packet) — `ARCHIVE_DONE.md` L10529
- D2i — Support-rule mesh gate with the idle clock off (COMPLETED 2026-09-21, review-accepted; GATE FAILED) — `ARCHIVE_DONE.md` L10566
- D2j-0 — The hot-disk edge test with today's model (COMPLETED 2026-09-21; its verdict was WITHDRAWN the same evening by D2j-0b) — `ARCHIVE_DONE.md` L10646
- D2j-0b — The hot disk under a descending exclusion plane (reviewer diagnostic, 2026-09-21 evening; no packet) — `ARCHIVE_DONE.md` L10698
- D2k — Side-face flux on exposed vertical faces (COMPLETED 2026-09-22, review-accepted; arbiter PASS, Meier gate FAIL and inconclusive) — `ARCHIVE_DONE.md` L10766
- D2l — Gate 0 scan: why the model drills 2.7× too fast (COMPLETED 2026-09-22; window EMPTY ⇒ pre-registered STOP, no code written) — `ARCHIVE_DONE.md` L10859
- D2m — Is the cut diameter set by the feet or by the flame? (COMPLETED 2026-09-22; MIXED verdict, reported as mixed) — `ARCHIVE_DONE.md` L10934
- D2o-0 — Prescribed floor/wall heating: does a cylindrical hole appear? (COMPLETED 2026-09-22; review-accepted; NEGATIVE, and the field had a defect found afterwards) — `ARCHIVE_DONE.md` L11012
- D2o-0b — Prescribed floor/wall heating with the field bounded in radius (COMPLETED 2026-09-22; REJECTED by review as incomplete; Q1 FAILS; **closes the prescribed-field line**) — `ARCHIVE_DONE.md` L11229
- D2p-0 — Lateral recession at exposed step faces (CLOSED 2026-09-22 at DESIGN TIME, zero run cost; premise refuted, pre-registered fallback triggered) — `ARCHIVE_DONE.md` L11387
- D2p-1 — Corner regularisation (CLOSED 2026-09-22 at DESIGN TIME, never hashed, never launched; refuted by arithmetic) — `ARCHIVE_DONE.md` L11498
- D2q-1 — Make `removed_mf` authoritative (per-cell void, stage 1: pure refactor) (COMPLETED 2026-09-22, review-accepted; GATE PASSED 2639/2639) — `ARCHIVE_DONE.md` L11609
- D2q-2b (checkpoint stage) — Per-cell removal: the mechanism, built and gated (2026-09-22; arbiter deliberately not launched) — `ARCHIVE_DONE.md` L11707
- D2q-2c — The arbiter: does the corner drain spall the wall? (2026-09-22/23, review-accepted) — `ARCHIVE_DONE.md` L11794
- D2r-0 (pre-flight stage) — closed at item 4, PLANNER ERROR, no code written (2026-09-23) — `ARCHIVE_DONE.md` L11939
- D2r-0b — SUPERSEDED BEFORE LAUNCH by a code fact found in review (2026-09-23) — `ARCHIVE_DONE.md` L12061
- D2r-0c — SUPERSEDED BEFORE LAUNCH, inert by construction (2026-09-23) — `ARCHIVE_DONE.md` L12145
- D2r-1 — Heating the SIDE faces of above-plane columns (COMPLETED 2026-09-23, review-accepted; (a) PASSES) — `ARCHIVE_DONE.md` L12218
- D2t — Retire the mask edge with a smooth radial taper, and score SHAPE (COMPLETED 2026-09-23; (a) PASSES) — `ARCHIVE_DONE.md` L12401
- D2u — Meier shape packet: WITHDRAWN BEFORE LAUNCH 2026-09-23 (fatal design flaw + a prediction refuted by data on disk) — `ARCHIVE_DONE.md` L12503
- D2v — the body-in-the-hole h map (COMPLETED 2026-09-24; (a) PASSES, but the ATTRIBUTION IS CONFOUNDED — planner error) — `ARCHIVE_DONE.md` L12600
- D2w — the 2x2 attribution: the MAP owns the funnel, the SIDE PATH owns the rate (COMPLETED 2026-09-24, review-accepted) — `ARCHIVE_DONE.md` L12684
- D2x — split `f`: the derived area is measured, and it does NOT reach Meier's hole (COMPLETED 2026-09-24) — `ARCHIVE_DONE.md` L12820
- D2y — conservative surface area: PARKED BEFORE LAUNCH 2026-09-24 (never ran, no code) — `ARCHIVE_DONE.md` L12909
- D2z — the step-height closure: PLANNED AND PARKED 2026-09-25, pending D2z-0 (never ran, no code) — `ARCHIVE_DONE.md` L12973
- D2z-0 — retire the half-cell removal increment (COMPLETED 2026-09-25; P1/P2 refuted by the METRIC, P3/P4 held) — `ARCHIVE_DONE.md` L13042
- D2z-1 — whole-cell removal unconditional, the key deleted (COMPLETED 2026-09-27, review-accepted) — `ARCHIVE_DONE.md` L13103
- P4 REFUTED, in the informative direction — `ARCHIVE_DONE.md` L13269
- P3a REFUTED — the planner's crater-wall assumption — `ARCHIVE_DONE.md` L13469
- P3 REFUTED — and the mechanism is the finding: THE SIDE TERM DEFENDS ITSELF — `ARCHIVE_DONE.md` L13634
- P7 CONFIRMED — `PD` jammed, and it was pre-registered — `ARCHIVE_DONE.md` L13667
- P0b — ONE surface temperature for the BC, the criterion and the removal (SOURCE, revision 7) — review-accepted 2026-09-30 — `ARCHIVE_DONE.md` L14061
- P1 — the gradient-aware firing criterion (SOURCE, revision 5) — review-accepted 2026-10-01 — `ARCHIVE_DONE.md` L14191
- P1b — unsupported hot rock comes off (SOURCE, revision 2) — implemented 2026-10-01, **(c) FAILED on a grid artefact** — `ARCHIVE_DONE.md` L14335
- P1c — a grid-isotropic support rule (revision 2) — implemented 2026-10-01, **FAILED both primary gates; the cause was found and P1d is WITHDRAWN** — `ARCHIVE_DONE.md` L14450
- P1e — the staircase side-face correction — implemented 2026-10-01, **correction RIGHT, stall CURED, both shape gates FAILED; the attribution then CONFIRMED by field diagnostic** — `ARCHIVE_DONE.md` L14538
- P1f — lateral flake drive + lateral removal (revision 11) — implemented 2026-10-02, **ACCEPTED, every gate green; the mode is CORRECT and nearly INCONSEQUENTIAL in scope** — `ARCHIVE_DONE.md` L15924
- P1g-0 — the corner-face temperature fix and the azimuthal witness — implemented 2026-10-02, **ACCEPTED; mechanism CONFIRMED, deep deficit only PARTLY closed, and a SECOND mismatch found** — `ARCHIVE_DONE.md` L16157
- P1g-0b — one criterion, one kernel: the coefficient cap and the four-site argument lock — implemented 2026-10-02, **ACCEPTED (revision 2); the third instance of one bug family closed, and (d1) REVERSED** — `ARCHIVE_DONE.md` L16456
- P1g-1 — the lateral flake sink (sub-cell lateral recession) — implemented 2026-10-03, **ACCEPTED AS BOOKKEEPING; PHYSICALLY INERT — the recession RATE is identical to whole-cell** — `ARCHIVE_DONE.md` L17968
- P1g-2 — the lateral flake sink, THERMAL half — implemented 2026-10-03, **ACCEPTED; the thermal half IS real and attributable — it delays the ONSET, and the long-run recession converges back toward whole-cell** — `ARCHIVE_DONE.md` L18239
- P1g-3 — the flake thickness from the rock — implemented 2026-10-05, **REJECTED ON PROVENANCE; the physics findings are reviewer-verified and load-bearing. Closure carried as ord 0 of P2-0** — `ARCHIVE_DONE.md` L18497
- P2-0 — the wall heat-transfer coefficient from the annulus flow — ord 0/1 run 2026-10-06, **SUPERSEDED: the premise was a duct law applied where the flow is a jet. ord 0 and ord 1 are KEPT; the source work was never built** — `ARCHIVE_DONE.md` L18558
- P1g-2b — the lateral sub-cell mode at Meier scale — run 2026-10-06, **ACCEPTED with `(b3)` FAILED on the thermal leg; three measurements nobody had made, and the widening attributed to the right cause at last** — `ARCHIVE_DONE.md` L18639
- P1g-2c — the thermal shift's energy leak, CLOSED; P3 refuted (review-accepted 10-06) — `ARCHIVE_DONE.md` L18706
- P2-1 — the wall heat-transfer coefficient, DERIVED; the typed constant was too LOW, and the wall self-limits (review-accepted 10-07) — `ARCHIVE_DONE.md` L18851
- P2-2a — the stress-free reference, and the retirement of `differential` and the support rule (2026-10-08, rev 7, ACCEPTED — all 11 criteria PASS) — `ARCHIVE_DONE.md` L19003
- P2-2b — the measured Grimsel properties and the relief factor (2026-10-08, rev 2, PARTIAL — ord 0-3 complete, `(c2)` failed as a gate, ord 4 NEVER RUN and deferred) — `ARCHIVE_DONE.md` L19107
- P2-2d — the at-temperature toughness: the LEVEL banked parameter-free, the exponent CHARACTERISED not validated, and four misattributed targets retired (2026-10-12, rev 6, ACCEPTED) — `ARCHIVE_DONE.md` L19170
- P2-2f — the measured Grimsel properties on Meier: the firing temperature TRANSFERS, the price closes the rock half, the residual is DELIVERY (2026-10-12, rev 2, ACCEPTED; archived 2026-10-13) — `ARCHIVE_DONE.md` L19275
- P2-3 — the floor coefficient to Kant, the wall on, the gas marched up the bore (2026-10-13, REJECTED by review on the `b3` gate — a physics negative with the mechanism measured; must-fixes carried into P2-4; archived 2026-10-13) — `ARCHIVE_DONE.md` L19663
- P2-4 — the cold body in the annulus: a second sink in the wall march (2026-10-13, ACCEPTED by review; archived 2026-10-14 with two advisor corrections carried) — `ARCHIVE_DONE.md` L20076
- P2-5 rev 1 — the Höser lab leg: REJECTED at `a1` before any leg ran (2026-10-14; four planner defects found offline at zero compute; revised as rev 2 the same day) — `ARCHIVE_DONE.md` L20482
- P2-5 rev 2 — the Höser lab leg: volume right, shape wrong, no overheating (2026-10-14, ACCEPTED by review; archived 2026-10-19) — `ARCHIVE_DONE.md` L20982
- P2-2e — the constraint depth: d_ref = 2.0 mm reproduces Kant's minimum flux, and the Meier burner jams at the pads (2026-10-24, COMPLETED, NOT reviewed; archived 2026-10-24) — `ARCHIVE_DONE.md` L21585
- P2-7 — Höser at measured heat: her data show a 4.7–6.5× uptake gap and the model takes in 2.7–3.0× her maximum (2026-10-24, ACCEPTED by review 2026-10-25; archived 2026-10-26) — `ARCHIVE_DONE.md` L21987
- P2-8 — the surface face list: one heat rule and one area rule for every face; six dead settings deleted (2026-10-26, rev 3, ACCEPTED by review 2026-10-26 — d2's negative PROVISIONAL; archived 2026-10-26) — `ARCHIVE_DONE.md` L22382
- P2-9 rev 1 — one line per surface cell: BLOCKED at a2 (vi) on P2-8's area rule (2026-10-27, line rule built, a2 i–v PASS, no long leg; archived 2026-10-27) — `ARCHIVE_DONE.md` L22833
- P2-9 rev 2 — one breaking rule per surface cell: normal, edge temperatures and fallbacks fixed; Höser run (2026-10-27, rev 2.1, ACCEPTED by review 2026-10-27 — a5 a recorded FINDING; archived 2026-10-28) — `ARCHIVE_DONE.md` L23220
- P2-10 — the old flame wall machinery retired (5,006 lines); the 45° residue diagnosed, unresolved (2026-10-28, ACCEPTED by review 2026-10-28 — a1 failed as written, exact on substance; archived 2026-10-28) — `ARCHIVE_DONE.md` L23783
- P2-11 — G1 probe: Goldstein's average profile used as a LOCAL value on Höser — GO (2026-10-28, ACCEPTED by review 2026-10-29; archived 2026-10-29) — `ARCHIVE_DONE.md` L24205
- P2-11b — X1, the surface-path heat rule on Höser: depth matches her, Ø freezes at the flange opening; the bounded check closes (2026-10-29, ACCEPTED by review 2026-10-29; archived 2026-10-30) — `ARCHIVE_DONE.md` L24468
- P2-12 — X1's path mode deleted, then restored the same day: net source change zero (2026-10-30, advisor at the user's instruction, no /implement; archived 2026-10-30) — `ARCHIVE_DONE.md` L24861
- P2-13 — One Grimsel rock, one model form: the shared input files and the Kant calibration deck in the production form; the constants stand (2026-10-31, ACCEPTED by review 2026-11-02; archived 2026-11-02) — `ARCHIVE_DONE.md` L24892
- P2-14a — The constraint depth on the face list (the merge's first slice), and one Poisson's ratio: the port is faithful; d = 2 mm puts a −122 K onset-vs-flux trend into the production form (2026-11-02, ACCEPTED by review 2026-11-02; archived 2026-11-02) — `ARCHIVE_DONE.md` L25306
- P2-15 — Two conditions for a spall: the buckling thickness condition, faithful; one span cannot serve both of Kant's rocks under a conductivity-only transfer (2026-11-02, ACCEPTED by review 2026-11-03; archived 2026-11-03) — `ARCHIVE_DONE.md` L25800
- P2-16 — One minimum-flux mechanism: the constraint depth retired; span ∝ 1/UCS survives at equal onsets, but Kant's two granites cannot separate it from one span (2026-11-03, ACCEPTED by review 2026-11-03; archived 2026-11-03) — `ARCHIVE_DONE.md` L26286
- D2z-2 — witness repair, derived conductance, Meier (COMPLETED 2026-09-28) — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)
- D2aa — heat delivery under the burner body (COMPLETED 2026-09-28) — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)
- D2ab — the face assignment and the gap-derived coefficient (COMPLETED — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)
- D2ac — Meier at Meier's kinematics: THE FLOOR-FLUX EXCESS (COMPLETED — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)
- D2ad — calibrate the floor flux on Fig. 8.8's standoff: THE WINDOW IS EMPTY — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)
- P0 — a real surface temperature, and a mesh-independent anchor (COMPLETED — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)
- D2s — WITHDRAWN BEFORE LAUNCH (2026-09-23); premise and arithmetic both — `ARCHIVE_DONE.md` L14624 (no record of its own; in the verbatim Current Status)

Active step:
- **P2-17 — Buckling on in production, the flake cap shown dead as a bound, and the Höser uptake ladder**
  (`ACTIVE_STEP.md`; 2026-11-03; the form is settled before the merge slices, advisor's recommendation, user's /plan).
  - Archived: P2-16 ACCEPTED (`ARCHIVE_DONE.md` *"P2-16"*): d retired; Kant's two granites cannot separate one span
    from span ∝ 1/UCS (DR-1, onsets undetermined); on a flat deck B is a net-flux condition (L_req ≈ k/q_net).
  - `spall.buckle_span = 0.012` goes into `rock_grimsel.inputs` (one span for both Kant granites; 13.2 mm, ∝ 1/UCS,
    runs as the alternative). The onset deck keeps B off by a declared override (a `check_shared` directive).
  - The 1.5 mm cap is dead as a bound under B (both act on μ·c_f·ℓ; crossover T_s* 1273 K at 12 mm), but it also
    sets each line's gather/remap window, so it stays at 0.0015, relabelled; decoupling is a later src change.
  - Höser face-list ladder (form × S, 19 legs of 30 s): does B move the uptake factor? Predicted S*(N)/S*(O) ≈ 1.05,
    a trim of the gap, not a closure.
Next milestones after the active step:
⚠⚠ **FIVE ABSORPTION PATHS STAND MEASURED** — firing rate, energy budget, wall surface temperature,
wall coefficient, standoff (the last only under `nozzle_descent = prescribed`, a modelling choice).
A 15.3 % leak was worth **0.16 mm** of bore; a **x1.430** rise **+9.2 %** of power. ⚠ **An open
research thread, deliberately NOT a standing rule** (§9.5). ⚠⚠ **The wall is closed as the width
lever, five ways**; P2-0's premise is **inverted** — typed 169.3 too low near the floor (`ARCHIVE_DONE.md` *"P1g-2b"*/*"P2-0"*/*"P2-1"*).

- ⚠⚠ **THE DIRECTION (user, 2026-10-14): ONE skin-failure law + one enthalpy field for flame AND beam** —
  `docs/project/DIRECTION_2026-10-14_single_skin_failure_law.md` (⚠ **addendum 2026-11-03: A superseded**; the
  minimum flux is now the buckling condition, P2-15/P2-16). Order: **A** constraint depth (P2-2e, RETIRED in P2-16 on
  form) → **B**
  softening s(T) → 0 toward the solidus, replacing the relief's physics (R0 kept) → **C** α(T) deposition with
  energy closure (`MMWSpalling.H` ~4423; the old `:6166` pointer was stale) → **D** the regime map (power density × absorption depth × spot).
  ⚠ L_span (P2-15) is a rock property: set once, never per case, metres not cells; one calibration with c_f.
- ⚠⚠ **PLAN OF RECORD (user, 2026-10-30; refined 2026-11-02).** Research principles: realistic over matching; predict
  unseen MMW tests; same rock ⇒ same properties, calibrated once (`docs/project/MODEL_CARD.md` §1, §12). In order:
  1. P2-13: one rock, plus the calibration deck — DONE (accepted 2026-11-02).
  2. **The law form first (inserted 2026-11-02):** P2-15 buckling condition B — DONE (accepted 2026-11-03) → P2-16
     retire `constraint_depth` and the span-rule test — DONE (accepted 2026-11-03) → P2-17 B on (active; items 3–4).
     **Then merge the two
     removal paths onto the face list, in slices:** P2-14b Meier's feet (`foot_body_clearance`) → P2-14c MMW (beam,
     vapour, melt, with the four MMW checks below, and the layer depth under an interior peak, which B needs) → retire
     the column path. The user decides port or retire for: `surface_normal`, coherence, `flake_length = flux`, cluster
     detach, the bit, the mechanics stress source (⚠ advisor: keep it, since finite-body ~10 GHz heating needs it);
     AMR when MMW needs it.
  3. The span L_span is a rock property. P2-16: Kant cannot separate one span (≈ 12.0 mm) from ∝ 1/UCS (Grimsel
     13.2 mm); P2-17 declares 12.0 mm in `rock_grimsel.inputs`, 13.2 the alternative. ⚠ span and c_f are ONE
     calibration (B sees L/c_f).
  4. The flake rule BEFORE 5, or the uptake factor absorbs its error (advisor 2026-11-02 §4) = **P2-17**: B ON before
     the uptake fit (user, 2026-11-03; in the rock file, not the form file); the 1.5 mm cap dead as a bound under B;
     c_f needs no separate fit (B sees L/c_f; rate lever ±5 %). The onset deck stays B-off by a declared override.
  5. Flame heat source from rig inputs by one rule (X1's structure; confined-hole standoff and decay; fix the Höser
     double count), plus ONE uptake factor calibrated on Höser 0 MPa.
  6. Check Höser and Meier with B on after 5 (Meier needs the feet slice; d's stall and −60 % sat where delivery was
     weakest).
  7. Blind predictions: Meier, and Höser's field-tool test (thesis Ch. 6); Höser 47 MPa as a pre-registered test of
     the confinement term (the model gives ~half of Kant's measured drop, P2-13 c5).
  - Before any MMW removal run: check the vapour front's energy use, make α(T) conservative, give the criterion a
    meaning for an interior temperature peak, and put Zhang back in the regression set.
- ⚠⚠ **OPEN PROBLEM, kept on the list (user, 2026-10-26): heat uptake at a spalling face.** Slow drilling takes in
  3–6.5× less than a calorimeter implies (Höser 0 MPa 4.7–6.5×, Meier ~3–5×) and nothing measured explains it; fast
  drilling is explained by Wilkinson's hot surface (`studies/heat_budget/`). At the right level the flux law, minimum
  flux and bowl follow unaided (`studies/whatif_uptake/`). To settle: Grimsel's surface T during slow drilling, a
  flux reading on a spalling face, or a hot-jet flow simulation of Höser's ring + cavity. Literature:
  `validation/hoeser/literature_uptake_gap_notes.md`. ETH (Höser/Kant) did not answer the user's emails.
- **PARKED (user, 2026-10-14):** P2-5b (Meier floor re-levelling), P2-6 (wall-only coefficient), the cold body
  and wall march (DELETED in P2-10), and the Höser lab delivery (side-face treatment; P2-5 carry-forwards DR-3/4/6/7/8,
  `ARCHIVE_DONE.md` *"P2-5 rev 2"*). Meier is a regression, not a target, until D.
- ⚠⚠ **MMW — THE TARGET APPLICATION, AND SPALLATION HAS NEVER RUN WITH THE BEAM ON IN A VALIDATED
  CASE.** Leg B and the Kant anchors set `beam.P0 = 0`; Zhang has the beam on with thermal stress
  **inert**. ⚠ P2-2's flux exponent decides whether the model extrapolates there at all: under a beam
  the heat arrives volumetrically with its peak **below** the surface, where none of the five
  absorption paths is guaranteed to hold.
- **Sources ON DISK** (`franck_material/`, notes in `validation/hoeser/`): Höser 2016/2018 (×2), Kant 2016,
  Kant 2017 JGR, Kant SGW 2017, Wilkinson & Tester 1993, Rauenzahn & Tester 1989, Meier/Brkic WGC 2015, LLNL
  Walsh/Lomov 2011–2013 (×4), Saksala; literature maps `validation/hoeser/literature_minimum_flux_notes.md` and
  `literature_uptake_gap_notes.md` (10-26).
  ⚠ Absent: Seagraves ARMA-2019-1799 (closest analogue, paywalled), Rossi 2020, Hu et al. 2018 (not "Shen"), MMW ablation-rate
  data (needed for D); Meier's thesis PDF is outside the repo (`~/Downloads/ThierryMeier_PhD.pdf`).

## Validation Story

- **Unit/analytical** — individual numerical components. **Zhang/Oglesby** — the MMW source and
  thermal model. **Hu** — the prescribed-T chain on `damage_law`.
- **Kant 2017 Central Aare** — the LEFM Sp onset criterion on `sp_weibull`. ⚠⚠ **From P2-2 it is
  also the CALIBRATION SET**: `R` and the flux exponent are fitted here, on
  `sp_kant_onset_meiercfg`, so that Meier stays a pure test set.
- **Höser 2016** — the floor criterion at a MEASURED flux (7.4 MW/m² on the axis, lab scale, no body
  in the hole; `validation/hoeser/`). From P2-5 the Meier floor's level derives from it.
- **Meier 2017 pilot** — ROP, hole diameter and volume. ⚠ Live error, **re-measured under FEET by P2-2f
  (2026-10-13)**: bore **-5.1 mm** (the model is now **narrower**; the +16.2 mm was prescribed
  descent), volume **1.85x**, rate **3.47 vs 1.5 m/h**. ⚠ Quotable only under standing rule 10.
  ⚠ P2-3's legs (10-13) are **rejected on `b3`** and not quotable; their budget read says power
  into rock ROSE 35 % with the wall on (`2026-10-13.md` §13).
- ⚠ **Rossi 2018** is a prescribed-T mechanics regression of the Sp depth scan, **not** a flake-size
  validation for self-heated drilling (retired 2026-09-15, S1).
- ⚠⚠ **MMW has no validated case with the beam on** — see the last milestone.

## Known Stale/Important Notes

⚠ **Cap: 40 bullets, 3 lines each.** Only caveats that constrain the active or next few steps.
Full prior text: `ARCHIVE_DONE.md`, heading *"ROADMAP `Known Stale/Important Notes` as it stood
before P1g-1"*. ⚠⚠ A compressed bullet is a **pointer**; read the archive before acting on one.

- ⚠⚠ **THE MODEL CARD (2026-10-30), `docs/project/MODEL_CARD.md`:** there are two removal paths (the face list refuses
  beam, vapour, melt, d_ref and the feet clearance); the Kant calibration decks are not the production form; P2-2f
  Meier and P2-7 Höser 0 MPa no longer run (§1, §5, §9). Read §12 before planning a validation step.
- ⚠⚠ **THE LITERATURE STUDY (2026-11-02), `Claude_markdowns/reports/Granite spallation and melt physics.md`:** H1 quartz-pinned
  onset, H2 fixed reference depth, H4 solid-state dielectric runaway, H5 surface tension under an interior peak are
  EXCLUDED; H3 two conditions and H6 absorptive melt SUPPORTED. α = 300 m⁻¹ has no measurement at 90–140 GHz; the user's
  rig is ~10 GHz (penetration ~100 mm cold ⇒ a lab block is a FINITE body: the stress references the mean temperature).
- ⚠⚠ **THE FLAT FLOOR ITSELF DRIFTS 10.4 % OVER 2 / 1 / 0.5 mm, on both keys** (P2-9 review DR-21; ℓ ≈ 1.6–1.9 mm
  under-resolved at 2 mm). No cross-mesh claim is cleaner than that; compare a tilt or a shape at the SAME mesh.
- ⚠ **ONE VOLUME ESTIMATOR PER COMPARISON** (P2-9 DR-15/16): face-list legs score V(f); P2-8 and key-off legs scored
  V(voids). S\* 0.1431 on V(f) is 0.1461 on voids; c3's 0.5 mm gap 5.4 vs 8.4 %. Name the estimator in every row.
- ⚠ **`RefillTempGhosts` IS `bc->FillBoundary`, NOT `FillPatch`** (P2-9 DR-3): equal at `max_level = 0` only; it
  omits interpolation from coarse. Fix before any AMR leg on the face list.
- ⚠⚠ **ONSET vs STEADY-STATE spall temperature are different observables** (Wilkinson 1993): the
  ONSET is flat with flux (0.2–7.4 MW/m², five rocks); the steady-state surface MEAN rises (Barre 500 →
  930 °C). Score a flux exponent on the model's time-mean floor T, never on the firing T.
- ⚠⚠ **THE EXHAUST RECIRCULATION HOLDS THE FLOOR POWER UP.** Under `jet_T_ent_mode = exhaust` the
  floor map ×0.743 cut `jet_P_face` only 7 % (1,870 → 1,743 W): entrained mass 2.07 → 3.23 g/s, cap
  +17 %. A floor-coefficient change is not a delivered-power change — read the budget (`2026-10-13.md` §13).
- ⚠⚠ **THE MEIER FLOOR GAS NEVER CAME FROM THE JET CLOSURE'S DECAY.** Under `exhaust` T_mix 1,504 K ⇒
  D_e 8.05 mm, core 8.6 D, T_stag 1,742–1,804 K at 85–107 mm; the SAME deck in `fixed` mode reads 720–831 K.
  P2-5 validates the decay first (bracketed `jet_core_length`); P2-5b must choose Meier's mode and core explicitly (`2026-10-14.md` §6).
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
