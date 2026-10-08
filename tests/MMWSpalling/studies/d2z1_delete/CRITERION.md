# D2z-1 CRITERION — hashed before the first scored run

Whole-cell removal made unconditional; `spall.removal_increment` deleted. Every
bar below is transcribed from `ACTIVE_STEP.md`; none was chosen by me.

Metric for (c)/(d): `col_rate` and the metadata reader are **imported unchanged**
from `studies/d2j0b_plane/analyze.py`, the same estimator the frozen `A_f02`
baseline was scored with. Frozen studies are read by path and never re-run;
re-baselined outputs land under this study's own `output/`.

Legs: `W2`/`W1` (arbiter at 2/1 mm, frozen argv, no key), `RB_C1_2mm`, `RB_D`.

---

## (a) SCORED — every test passes

All `tests/MMWSpalling/*` regressions plus the **Hu** and **Kant/von Rohr**
validations. **Zero tolerance failures.** No test in the tree carries byte-exact
reference data (verified: **0** `reference/` directories), so a +0.4 % timing
shift should move none of them.

## (b) SCORED — the tree re-pin audit

Every output that changed must satisfy **all** of:

* **no column appears or disappears** in any `thermo.dat` or CSV header, except
  the two Goal 3/4 additions, which must be listed by file:
  * `firings_per_void` — added **unconditionally** (Goal 4), so it is expected in
    **every** `thermo.dat`;
  * `sn_recession_superseded` — added **only** under `spall.surface_normal = 1`
    (Goal 3ii), so it is expected in exactly the sweep members that set it.
* **no test crosses its own tolerance** (that is (a));
* **the sign is consistent** — aggregate cells voided **up or equal**; ⚠ in a
  tiny unit domain a one-event reshuffle can read **−1**, so a decrease of **at
  most one event** is acceptable **with an explanation**, and anything larger is
  a failure;
* **`firings_per_void == 1`** wherever reported, and `n_voided == 1` asserted.

⚠⚠ **DO NOT AUDIT ON "aggregate rate up 0.2–0.6 %".** D2z-0 established that band
is **below the metric's resolution**: block spread was ±2 pts at 2 mm and two
reasonable estimators disagreed by ~0.4 pts on the aggregate — the size of the
effect. **Direction and structure are resolvable at this effect size; magnitude
is not.**

### Re-pin procedure (written before the sweep ran; Goal 6)

⚠ There is **no `run.py`** and no `--check`/`--repin` in `studies/tree_sweep_0922/`.
The folder holds `RESULTS.md` and `witness/{hash.sh, wit.sh, post_all.sh,
guarded_post_all.sh}` plus md5 lists. The sweep runs the tests against the binary
at its hardcoded path and writes an md5 list; a re-pin is a **declaration** about
which list is the reference, not a script flag.

1. `bash witness/guarded_post_all.sh d2z1` → writes `witness/d2z1.md5` and
   `witness/d2z1_wit/`. The guard fingerprints the binary before and after and
   fails loudly if it moved mid-sweep.
2. **Audit** `diff witness/tree0922.md5 witness/d2z1.md5`, classifying every
   differing path against the four bullets above. The diff is kept at
   `studies/d2z1_delete/tree_repin_audit.txt`.
3. **The new reference is `witness/d2z1.md5`**, declared in
   `witness/REFERENCE.md`.
4. **`witness/tree0922.md5` is preserved untouched** as the superseded
   pre-D2z-1 reference. Nothing is deleted or overwritten.

## (c0) SCORED — THE CHEAPEST EXACT CHECK: byte-identity against D2z-0

With the key deleted, `W2`/`W1` must be **byte-identical to D2z-0's `F2`/`F1`**
in the removal log and **every cell plotfile**, with `thermo.dat` differing
**only** by the header columns Goal 3/4 add. Same physics, key removed, identical
outputs — **anything else means the deletion changed behaviour.**

## (c) SCORED — cells voided, not rate

Over the **single 50–200 s window** (pre-registered length, per D2z-0
takeaway 5), total cells voided inside r < 25 mm on the arbiter at both meshes,
against the frozen `A_f02` pair. **Bar: +0.2 % to +1.0 %.**

Counting cells is the estimator that is not window-sensitive: raw cells per 50 s
block read +0.08 / +2.03 / +0.07 % where `col_rate` read −1.73 / +4.35 / −1.38 %.
**If a criterion is about how much rock left, count cells.**

## (d) REPORTED, not gated — aggregate rate and mesh gap

50–200 s single window, r < 25 mm, against the frozen pair. Reference:
**+0.352 % / +0.438 %** and an aggregate gap shift of **+0.085 pts**.
⚠ **Reported only** — a bar here is what failed D2z-0.

---

## Pre-registered predictions — contiguous, no dead zones

* **P1** all regressions and both validations pass: **zero** failures held; any
  failure refuted.
* **P2** cells voided **+0.2…+1.0 %** at both meshes held; outside refuted.
* **P3** `firings_per_void == 1` and `n_voided == 1` everywhere held; any other
  value refuted.
* **P4** ⚠ **RESTATED.** `ACTIVE_STEP.md` is internally inconsistent here: Goal 3
  specifies "a loud parse-time message plus an on-disk witness, **not an
  abort**", with a verified reason (two sweep members set the key, so an abort
  makes (a) and Goal 6 unmeetable), while the Guardrails, Goal 5 and P4 still say
  "abort" from an earlier draft. **Goal 3's reasoning is verified on disk and is
  taken as authoritative.** P4 is therefore: with `spall.surface_normal = 1`, the
  parse-time message fires naming the superseded scaling **and** the
  `sn_recession_superseded` column appears, **and the run proceeds** — held;
  running with no message or no column, or aborting, refuted.
* **P5** the tree re-pin audit finds **no** column change beyond the two
  deliberate additions above held; any other header change refuted.

**Refuting evidence on disk:** D2z-0's `F2`/`F1` already ran whole-cell at both
meshes with `surface_normal` **off**. If any regression moves by more than its
tolerance, compare against those two before looking for a physical cause.
