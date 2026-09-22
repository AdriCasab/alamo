# Active Step: D2m - Is the cut diameter set by the feet or by the flame?

Status: completed (MIXED verdict — H-feet on min Ø, H-thermal on the erosion
edge; key-only, binary unchanged, packet amended mid-flight and runs redone)

Planned 2026-09-22. **Key-only, no source change, no build, ~1 h.** One
question, one falsifiable prediction, five short 2 mm runs. It decides whether
half the shape error is a support-rule artefact before any heat-transfer work
is written.

## Context

- **D2l's finding that outranks its own verdict:** with the side-face term
  entirely off — the D2i/D2d configuration — the model already removes
  **4.58 cm³/s against Meier's 2.47 (≈ 1.9×)** while the burner rate looks right
  (1.49 vs 1.30–1.5 m/h). **The rate agreement was partly a coincidence of two
  errors.** Shape is the honest signal: the model cuts at **78 mm** (Meier
  85–93) and flares to 145 mm at the mouth.
- **The observation this packet tests.** The cut diameter across D2l's whole
  ladder was **78.0 / 77.1 / 77.7 / 77.5 / 78.0 mm** — flat to 0.9 mm while
  absorbed power **doubled** (1639 → 3204 W) and the rate changed 2.7×. **A
  thermal diameter would move with power. This one does not.** Something
  geometric is pinning it, and the burner's feet span **28–40 mm** in radius,
  i.e. Ø 80 — within 2–3 mm of the observed cut.
- **The mechanism, from the code.** `FeetDescent`
  ([MMWSpalling.H:1716](src/Integrator/MMWSpalling.H#L1716)) reads only live
  columns with r in [`foot_r_inner`, `foot_r_outer`], takes the 0.9 quantile per
  pad, and sets `z_target = z_foot + foot_standoff` (monotone descent).
  `BuildFlameColumns` then gives every column `s = z_nozzle − z_face` and
  `h = h_expr(r, s)`, with `h = 0` where `s ≤ 0`. **So the feet touch the heat
  distribution through exactly one scalar — the nozzle height.** The
  consequence: the 28–40 mm ring is the *only* rock held at the design stand-off
  for the whole run, so pad recession and ROP are locked together, while rock
  outside 40 mm is left behind, its stand-off grows and its heating falls. That
  predicts a sharp erosion-rate transition at r ≈ `foot_r_outer`.
- **Why it matters either way.** Meier's hole was 85–93 mm with the *same* Ø 80
  feet — his rock was removed **beyond** the tool's stance. If our diameter
  tracks the feet, that mechanism is missing and no wall-`h` work will supply
  it. If it does not track, the width is thermal and the gas-path packet is
  properly motivated.

**Caveats carried forward:**
- **Nothing from D2b–D2d may be quoted as a Meier result**; no D2g Stage 2.
- **This is a relative comparison at one mesh.** The support rule is
  mesh-checked only to ≈ 350 s, and D2l showed the cascade is invisible at a
  single mesh, so **no number here may be quoted as a Meier score** — only
  case-to-case differences at fixed mesh and fixed geometry.
- **Report the full Ø(z) profile and the minimum, never one depth** (D2l: the
  profiles cross, and the matched depth landed near a crossing). **Compare at
  matched feet height, never at matched time** — Ø at matched time tracks dwell.
- **The 4–6 J/mm³ "to-rock band" is not a measurement** (it is an assumed 30 %
  delivery) and must not be used as a target. The only independent number on
  that axis is the thermodynamic floor, 1.147 J/mm³.
- `side_face_flux` stays **off**: D2l closed that line (`f_min` 0.5 > `f_rate` 0).
- **Do not rebuild.** The working tree carries unreviewed AMR changes and the
  2639-file sweep is owed before anything compiles. This packet uses the
  existing binary and must not trigger `make`.
- `run.py --kill NAME`, never `pkill -f`; `caffeinate -dis` on every launch;
  load sibling harnesses by path; copy scripts into this study rather than
  running them inside a frozen folder.

## Sources

- **Binary:** `bin/mmwspalling-3d-g++` = `62dea451…`. Verify at start and end.
- **Base configuration:** D2l's `M_f0` = the D2i `LI0` Meier keys at 2 mm with
  `side_face_flux = 0`, stop 250 s. Reproduce it exactly as the control.
- **Frozen, import by path, never edit:** `d2l_scan/` (`run.py`, `analyze.py`,
  and its Ø(z), volume and matched-depth estimators), `d2i_gate/`,
  `d2j0b_plane/`, `d2k_sideface/`.
- **Keys varied (existing, no code):** `surface_patch.foot_r_outer`,
  `surface_patch.foot_r_inner`. Everything else — `foot_rule = pads`,
  `foot_pad_quantile = 0.9`, `foot_body_clearance = 0`, `foot_standoff`,
  `pinned_idle_cycles`, the gas closure — unchanged.
- **Reference numbers:** the control (r_outer 40 mm) expects ROP 1.491 m/h,
  V̇ 4.58 cm³/s, min Ø 78.0 mm, depth-mean Ø 107.6 mm, absorbed 1639 W (quarter).
  Meier: 1.30 m/h, 2.47 cm³/s, Ø 85–93 mm, feet Ø 80, stand-off 50 mm.
- **Cost:** a 2 mm Meier leg to 250 s ≈ 11 min; five legs in parallel ≈ 20 min
  wall, ~1 h with analysis.

## Goal

1. **New study `tests/MMWSpalling/studies/d2m_feet/`** with `run.py`,
   `PREDICTIONS.md` (hashed, read-only, **written before `output/` exists**),
   `analyze.py`, `RESULTS.md`, one figure.
2. **Cases**, all 2 mm, D2i `LI0` keys, `side_face_flux = 0`, stop 250 s:
   - **`R40`** — control, `foot_r_outer = 0.040` (unchanged).
   - **`R32`** — `foot_r_outer = 0.032`.
   - **`R44`** — `foot_r_outer = 0.044`.
   - **`R48`** — `foot_r_outer = 0.048`.
   - **`S3648`** — `foot_r_inner = 0.036`, `foot_r_outer = 0.048`: the ring
     **shifted outward at constant width**, to separate "the outer edge sets it"
     from "the ring as a whole sets it".
3. **Identity guard:** `R40` must be **byte-identical** to D2l's `M_f0` over
   every shared thermo row. If it is not, **stop** — the harness is reaching
   through the frozen studies incorrectly and nothing else here is meaningful.
4. **`PREDICTIONS.md`, hashed before launch**, stating two falsifiable
   alternatives and the reasoning's own sign test:
   - **H-feet:** min Ø tracks the stance, i.e. **min Ø ≈ 2 × `foot_r_outer`
     within 5 mm** (≈ 64 / 78 / 88 / 96 mm for r_outer 32 / 40 / 44 / 48).
   - **H-thermal:** min Ø stays **78 ± 3 mm** at every r_outer.
   - **ROP direction:** a larger ring reaches into cooler, less-eroded rock, so
     its 0.9 quantile sits higher, the burner rides higher, stand-off grows and
     heating falls — **predict ROP falls as `foot_r_outer` rises.** State the
     expected size. **A rising ROP refutes the reasoning even if the diameter
     tracks**, and must be reported as such.
   - **`S3648` vs `R48`:** equal if the outer edge alone sets the diameter;
     different if the whole ring matters.
   - The control's numbers (ROP 1.491, V̇ 4.58, min Ø 78.0) are an identity
     check, not a prediction.
   - **The transition radius (see 5a) is the sharpest form of H-feet and is
     pre-registered as deciding alongside min Ø:** the radius at which the
     recession rate falls to half its central value **tracks `foot_r_outer`
     within one cell (2 mm)** under H-feet, and does not move under H-thermal.
     Min Ø is an integrated outcome; the transition radius is the direct
     signature of the mechanism.
5. **Report per case**, all at **matched feet depth** (the deepest depth every
   case reaches): ROP 150–250 s; absorbed power; V̇; **the full Ø(z) profile at
   20 mm intervals**; **min Ø and the radius where it occurs**; depth-mean Ø;
   J/mm³ to rock against the floor 1.147 (**not** against 4–6); ledger closure;
   `jet_s_c` and `patch_min_standoff`.
5a. **The recession-rate profile — the measurement everything else hinges on.**
   For each case, **rock recession rate against radius, from the axis out to
   70 mm**, in 2 mm annular bins, over the scoring window, plotted as
   `d2m_vprofile.png` (one line per case, each case's `foot_r_outer` marked as a
   vertical line, and Meier's implied floor half-width of 42.5–46.5 mm shaded).
   - **What it answers that nothing else does:** whether the model's floor
     recedes **flat out to ~45 mm**, as Meier's cylindrical, calibrated hole
     requires, or is **already falling at 40 mm**. Without it, D2m can confirm
     that the feet set the minimum diameter and still leave us guessing *why*
     Meier's hole is wider than his tool.
   - Report the central rate, the half-value radius, and the 10–90 % fall-off
     width — the same estimators D2j-0 used for the hot disk, so the two are
     comparable.
   - Take it from the removal log at matched feet depth, not from plotfile
     differencing, so it is not tied to plotfile cadence.
6. **The decision, pre-registered:**
   - **H-feet holds** → the cut diameter is a support-rule artefact: the model's
     hole is the width of its own tool stance while Meier's was wider than his.
     The next packet targets what removes rock **beyond** the feet (the gas path
     up the wall), and **no Meier diameter may be quoted until that exists**.
   - **H-thermal holds** → the feet hypothesis is dead, the width is thermal, and
     the gas-path packet proceeds with the wall as its target.
   - **Mixed** (partial tracking, or `S3648` ≠ `R48`) → report the split; do not
     force a single reading.
7. **`RESULTS.md`:** §0 hashes and binary; §1 identity; §2 the ladder table;
   §3 the Ø(z) profiles; §4 the recession-rate profiles and the transition
   radii; §5 the decision and what the next packet inherits. Figures
   `d2m_profiles.png` (Ø against depth, one line per case, Meier's 85–93 band
   shaded) and `d2m_vprofile.png` (recession rate against radius, per 5a).
8. **Completion notes and takeaways.**

## Guardrails

- **Key-only: no source edit, no build, no `make`.** The binary sha256 must be
  unchanged at the end. If anything appears to need a rebuild, **stop and say
  so** — the 2639-file sweep on the working tree is owed first.
- **`side_face_flux = 0` in every case.**
- **Change only `foot_r_inner` / `foot_r_outer`.** Not the quantile, not the
  stand-off, not the gas closure, not the idle key.
- **Fix `PREDICTIONS.md` before launch**; do not re-window, re-bin or re-choose
  the deciding metric after seeing results. If amended, re-hash, kill and delete
  every affected run, and record it.
- **No Meier score from this packet** — differences between cases only.
- **Frozen:** every earlier `PREDICTIONS.md`/`CRITERION.md`, harness, analysis
  script and output. Import by path; copy, never run in place.
- Do not touch `ext/`, `bin/`, `obj/`, `build/`, `compile_commands.json`,
  `configure`, `LICENSE`, or other integrators. **Do not commit.**

## Commands

```bash
cd /Users/tzetze20/amr_tools/alamo
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
S=tests/MMWSpalling/studies/d2m_feet
shasum -a 256 bin/mmwspalling-3d-g++        # expect 62dea451... (start and end)
# write PREDICTIONS.md, chmod a-w, record its sha256 in RESULTS §0
$PY $S/run.py --list
caffeinate -dis $PY $S/run.py --cases R40 R32 R44 R48 S3648 --jobs 5   # ~20 min
$PY $S/run.py --identity                    # R40 vs d2l_scan M_f0
$PY $S/analyze.py
$PY tests/MMWSpalling/validation/meier/sp_meier_pilot/test_feet
$PY tests/MMWSpalling/validation/meier/sp_meier_pilot/test_feet_d2c
```

## Claude completion notes

**Outcome: MIXED, as the packet's third decision row allows — reported as mixed,
not forced into one reading.** The cut diameter *is* a support-rule artefact
(min Ø tracks 2×`foot_r_outer` to 0.5–2.1 mm in both directions), but the
flame's far-field erosion is stance-invariant (absolute rate at 45–50 mm the
same to 5–7 % across a 16 mm change of stance). **Meier's hole was wider than
his tool and ours cannot be** — the mechanism that removes rock beyond the
stance is missing, and no wall-`h` work supplies it.

**Key-only throughout.** No source edit, no build, no `make`; binary unchanged
at `62dea451…`; nothing committed.

**Files created** (all new, all under `tests/MMWSpalling/studies/d2m_feet/`):
`run.py`, `PREDICTIONS.md` (read-only, sha256 `d735f95f…`, 12:22:10),
`analyze.py`, `analysis.md`, `RESULTS.md`, `d2m_profiles.png`,
`d2m_vprofile.png`, `campaign.log`, `campaign_killed_amendment.log`, `output/`
(5 runs). **No file outside this directory was modified except this one**; no
frozen study was edited (`d2l_scan`, `d2i_gate`, `d2j0b_plane`, `d2k_sideface`,
`d2b_feet_rop`, `d2c_steady`, `d2j0_hotdisk` all untouched, verified by mtime).

**The packet was amended mid-flight and the runs were redone.** `ACTIVE_STEP.md`
changed at 12:16:08 — after `PREDICTIONS.md` was hashed at 12:15:53 and while
the first launch was starting — promoting the transition radius to a **deciding**
metric and adding item 5a. The five in-flight runs were killed at t ≈ 22–28 s of
250 and deleted (**no `.done`, no analysis, no output inspected**),
`PREDICTIONS.md` was amended and re-hashed (`6d1daf69…` → `d735f95f…`), and
everything was relaunched from t = 0 at 12:23. All five completed 12:40–12:42,
after the hash. The amendment changed no input key, so the discarded runs were
physically identical; they were discarded anyway, per the standing rule. Record
at the top of `PREDICTIONS.md`; killed log kept as `campaign_killed_amendment.log`.

**Runs: 5, all rc = 0**, 17.0–18.9 min each, 5 in parallel under `caffeinate -dis`
(harness self-asserts via `keep_awake()`).

**Tests run:**

| check | result |
|---|---|
| `run.py --identity` — `R40` vs `d2i_gate/LI0_2mm`, t ≤ 250 s | **PASS, byte-identical** (2605 rows) |
| `validation/meier/sp_meier_pilot/test_feet` | **PASS** (all outcomes as recorded) |
| `validation/meier/sp_meier_pilot/test_feet_d2c` | **PASS** (all outcomes as recorded) |
| binary sha256 unchanged start → end | **PASS** |
| energy ledger, every run | max abs `ledger_err` ≤ 4.8e-15 |
| `jet_P_face` ≤ `jet_P_cap` | max ratio 0.21 |

**Not run, with reasons:**
- No unit test: the packet changed no source.
- **The 2639-file sweep on the working tree is still owed** before any rebuild
  (D2l Gate 1, never reached). Untouched here because D2m needed no build.
- A `jet_free_surface = 0` stance ladder, which would clean up the ROP numbers
  (see takeaway 6), was **not** run — out of scope ("not the gas closure").

Full results: `studies/d2m_feet/RESULTS.md` (§0 hashes, the amendment record and
the two source couplings; §1 identity; §2 the ladder; §3 `S3648`; §4 the erosion
edge **and a retracted reading**; §5 the cone finding; §6 the free-surface
confound; §7 what the next packet inherits). Figures `d2m_profiles.png`,
`d2m_vprofile.png`.

## Implementation takeaways

1. **The cut diameter is the tool's stance.** min Ø = 63.5 / 79.3 / 87.0 / 93.9 mm
   against stances 64 / 80 / 88 / 96 — tracking to **0.5–2.1 mm over a 32 mm
   swing, in both directions**. `S3648` confirms the **outer edge alone** sets it
   (min Ø +0.9 mm, ROP −2.1 % against `R48`).

2. **But this was pre-registered as weak evidence, and still is.** `z_foot` is the
   **0.9 nearest-rank quantile** of annulus heights, so ~90 % of rock out to
   `foot_r_outer` sits at or below the feet: the hole is ≥ 2·r_outer at the feet
   depth *almost by construction* (`PREDICTIONS.md` §2a, written before the runs).
   **Do not let a future packet quote the tracking as a discovery.** The
   informative results are 3, 4 and 5 below.

3. **My H-asym prediction was refuted on its inward branch.** I expected `R32` to
   stay near 78 mm because `h = h_expr(r, s)` does not depend on the foot radii,
   so rock at 32–40 mm is heated identically. It collapsed to 63.5. **Being
   heated is not sufficient:** rock outside the stance is never brought under the
   nozzle, its stand-off grows, and it lags permanently. The feedback through
   `z_foot` is stronger than the direct heating term.

4. **A retracted reading, kept in `RESULTS.md` §4 for the reviewer.** `r_half`
   spans 16.0 mm for 16 mm of `r_outer` (slope **1.02**), which I first wrote up
   as "the erosion edge moves one-for-one with the stance". **That is an
   artefact.** `r_half` is a crossing of `0.5·v_c`, and `v_c` falls 3.32 → 1.97
   m/h across the ladder, so the *level* moves and the crossing moves with it.
   **Absolute rates at fixed radii settle it: 45–50 mm is invariant to 5–7 %
   while the centre moves 39 %.** — **A normalised metric is not safe to compare
   across cases whose normaliser changes.** This is the third metric trap in this
   line of work, after "Ø at matched time tracks dwell" (D2l) and "Ø at a single
   depth lands near a profile crossing" (D2l). Always check the absolute
   quantity before reading a normalised one.

5. **The floor is a cone with no edge — the answer item 5a was added to get.**
   At the control, v(r)/v_c = 0.83 at 29 mm and **0.62 at 39 mm**, still ≈ 0.36
   at 69 mm, and it never reaches 0.1·v_c inside 70 mm, so **`w_edge` is
   undefined for four of five cases**. Meier's floor half-width is 42.5–46.5 mm;
   ours is removing rock there at roughly half the centre rate. **Meier drilled a
   cylinder, the model drills a cone**, and the support rule does not explain it.
   This is the next packet's target and it lives in the **radial** heat
   distribution — which D2l showed no existing key can reach.

6. **The free-surface confound is bigger than the pre-launch bound.**
   `foot_r_outer` also sets the free-surface boundary (`MMWSpalling.H:2062`), and
   `jet_fs_cols` runs **0 / 18 / 87 / 203 / 247** across the ladder. The bound I
   declared before launch came from D2l's `P_fsoff0`, measured when there were
   ~18 such columns, so **it does not cover the outer end.** min Ø survives by
   ~100×; **the ROP percentages carry the uncertainty and are direction-plus-size,
   not calibrated.** A clean ladder needs `jet_free_surface = 0` throughout.

7. **`foot_r_inner` silently sets `nozzle_collision_radius`**
   (`MMWSpalling.H:5695`, `feet ? foot_r_inner : radius`) whenever the key is not
   given, which the D2i chain never gives. Pinned at 0.028 in every case;
   **identity PASS proves the pin inert.** Any future packet touching
   `foot_r_inner` must pin it too.

8. **ROP falls monotonically with stance** (1.925 → 1.179 m/h) while **absorbed
   power rises** (1258 → 2021 W): a wider stance means a wider hole, more
   absorbed power, and slower descent. The pre-registered direction holds; sizes
   were +29 % (predicted +10…+25) and −21 % (predicted −25…−45).

9. **J/mm³ to rock is 1.40–1.44 (1.22–1.26 × the floor) and flat across the
   ladder** — the same fixed removal cost D2l found, unmoved by stance. Combined
   with D2l: volume rate is set by absorbed power at a fixed conversion, and
   **nothing in the support rule changes that**.

10. **Orientation only, explicitly not a Meier score** (single mesh, 250 s,
    relative study; the packet forbids scoring): `R44` lands at ROP 1.301 m/h and
    min Ø 87.0 mm, both inside Meier's ranges, while removing **5.11 cm³/s
    against 2.47** — because depth-mean Ø is 120 mm. **Matching rate and cut
    diameter together still leaves volume 2× high.** The funnel above the cut is
    where the volume goes, and it is not the cut that needs fixing.

11. **"D2l's `M_f0`" does not exist.** D2l's f = 0 row was read from
    `d2i_gate/output/LI0_2mm`; `R40` reproduces that, named explicitly in
    `run.py` and `RESULTS.md` rather than substituted silently.

## Review findings

Reviewed 2026-09-22, in a separate session from the implementation. Nothing
was rerun: all five runs are present and every number I recomputed matches.

Verdict: accepted

Accepted as work, with the MIXED verdict as reported. The packet's third
decision row allows it, and the split is reported rather than forced.

### What was checked
- **Key-only.** No file under `src/` is newer than the AMR session's 00:01
  edits, the binary still hashes to `62dea451…`, and no frozen study
  (`d2l_scan`, `d2i_gate`, `d2j0b_plane`, `d2k_sideface`, `d2b_feet_rop`,
  `d2c_steady`, `d2j0_hotdisk`) has a file newer than this step's start.
- **The mid-flight amendment, checked rather than taken.** This is the
  claim most worth auditing, and it holds:
  - I recovered the first `PREDICTIONS.md` from the implementation
    transcript (11:15:52Z) and hashed it: exactly
    `6d1daf69895a0ba5…`, the hash the amendment record quotes.
  - The rewrite (11:22:02Z) carried a `12:xx` placeholder, which a `sed`
    replaced with the real time before `chmod 444` — which is why the
    on-disk hash is `d735f95f…` at 12:22:10 rather than the write's own
    hash. Benign, and the sequence is visible in the transcript.
  - The relaunch came at 12:23:29, **after** the re-hash, and the
    implementer verified `output/` was absent first. All five `.done`
    files are 12:40–12:42.
  - `output/` holds only the five scored runs; the killed launch survives
    as `campaign_killed_amendment.log`.
  - The amendment changed no input key, so nothing physical was at stake —
    and the runs were discarded anyway, which is the standing rule.
- **Identity guard.** `R40` is byte-identical to `d2i_gate/LI0_2mm` over
  all 2605 thermo rows to 250 s, so the harness reproduces the control.
- **The ladder, recomputed from thermo:** ROP 1.925 / 1.491 / 1.301 /
  1.179 / 1.154 m/h and absorbed power 1258 / 1639 / 1841 / 2021 / 2052 W
  across R32 / R40 / R44 / R48 / S3648 — matching the notes, and the
  pre-registered ROP direction holds.
- **Minimum Ø, recomputed independently** (cumulative recession from the
  removal log, 9 azimuth sectors, at the matched feet depth): 62.6 / 77.5 /
  84.9 / 91.8 / 93.6 mm against stances of 64 / 80 / 88 / 96 / 96. My
  estimator differs from theirs (events versus the removed mask) and lands
  1–4 mm low, but the tracking is the same and `S3648` matches `R48`
  within 2 mm, confirming the outer edge alone sets it.
- **The retracted reading is correct to retract, and I verified both
  halves:**
  - `r_half` really does track: 37.4 / 46.6 / 50.6 / 53.5 / 54.3 mm, a
    slope of 1.01 over the 16 mm swing, and offsets +5.4 / +6.6 / +6.6 /
    +5.5 mm — so the pre-registered ±2 mm test does fail as written;
  - but `v_c` falls 3.317 → 1.968 m/h across the ladder, so the level the
    crossing is taken at moves with it;
  - the absolute rates settle it, and my numbers reproduce their §4 table
    exactly: at 45 mm the spread is **7 %** and at 50 mm **4 %**, while at
    25 mm it is **37 %**.
  - So H-thermal on the erosion edge and H-feet on min Ø is the right
    reading.
- **The matched-geometry window is the pre-registered one** and is what I
  had to use: on a fixed 150–250 s window the same profiles give a
  completely different picture (R32's 45–50 mm rate reads 2.7× lower),
  because the faster cases have drilled past the comparison depth. The
  cross-check and the disagreement are already reported in RESULTS §4.
- **Other checks:** ledger ≤ 6.6e-15 on every run; `jet_fs_cols` runs
  0 / 18 / 87 / 203 / 247, which is the confound takeaway 6 declares;
  `test_feet` and `test_feet_d2c` PASS.
- **Takeaways are durable.** 2 (the tracking is near-structural), 4 (the
  normalised-metric trap), 5 (the cone), 6 (the free-surface confound) and
  7 (`foot_r_inner` silently pins the collision radius) are all things the
  next planner needs.

### Findings for the next planner (non-blocking)
- **The headline is takeaway 5, not the tracking.** The floor has no edge:
  v(r) is still ≈ 0.36·v_c at 69 mm and never reaches 0.1·v_c inside
  70 mm, so `w_edge` is undefined for four of five cases. Meier drilled a
  cylinder to a half-width of 42.5–46.5 mm; the model drills a cone. That
  is a **radial** distribution problem, and D2l already showed no existing
  key reaches it.
- **Takeaway 2 is the honest caveat and should survive archiving:** with a
  0.9 quantile stance, "hole ≥ 2·r_outer at the feet depth" is close to
  structural. The tracking confirms the mechanism but is weak evidence, and
  a later packet must not quote it as a discovery.
- **The ROP percentages carry the free-surface confound** (`jet_fs_cols`
  0 → 247 across the ladder), as takeaway 6 states. A clean stance ladder
  needs `jet_free_surface = 0` throughout; it was correctly left out of
  scope here.
- **Takeaway 10 is the number to carry into any future scoring:** at R44
  the rate (1.301 m/h) and the cut diameter (87 mm) are both inside
  Meier's ranges while volume is 5.11 cm³/s against 2.47, because the
  depth-mean is 120 mm. Matching rate and cut diameter together still
  leaves volume 2× high, so the funnel above the cut is where the error
  lives.
- **Still owed before anything compiles:** the 2639-file reference sweep on
  the working tree with its unreviewed AMR delta. Three key-only packets in
  a row have deferred it; the next packet that builds cannot.
