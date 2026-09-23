# D2q-2b results — per-cell removal

**Status: IN PROGRESS.** Item 0 (verification scaffolding) coded; items 1–3
(mechanism) and item 4 (arbiter) pending.

## 0. The two questions, answered before coding the mechanism

The packet requires both answered explicitly, with the line that implements the
answer. Both are about removing a cell that is **not** the column top — which
nothing in this codebase has ever done.

### Q1 — What happens to `phi` when a non-top cell is removed?

**Answer: nothing is written to `phi`, and that is correct — but only because of
two properties that must not silently change.**

`phi` is the per-column sub-cell remainder: `phi(k) = base − (k−klo)·dz`,
decremented for the whole column by `h_c` each step, flipping a cell to void
when it crosses zero ([Removal.H:1418](../../../../src/Integrator/MMWSpalling/Removal.H#L1418)).
An interior removal does not move `k_top`, so there is no remainder to carry and
no correct value to write. Leaving `phi` alone is the answer, and it is safe
because:

1. **Every `k_top` decision reads `rem` alone** after D2q-1, so the scans still
   find the true top above an interior void.
2. **The already-void guard is on `rem`, not `phi`**
   ([Removal.H:1418](../../../../src/Integrator/MMWSpalling/Removal.H#L1418),
   `if (rem(i,j,k) > 0.5) continue;`), so a cell removed laterally is not flipped
   or counted a second time when the column later recedes through it.

**The hazard this exposes, and it is real.** `phi` is still the sole input to
`UpdateSurfaceMaskFromPhi`
([MMWSpalling.H:4655](../../../../src/Integrator/MMWSpalling.H#L4655)) — D2q-1's
unconverted site 7. For an *interior* removal the mask is unaffected, because the
top cell is unchanged and still owns the window `0 ≤ phi < dz`. But if per-cell
removal ever voids a **top** cell without decrementing `phi`, the window still
points at that now-void cell and the surface flux would be applied to a void.

**Why that cannot bite here:** `surface.follow_mask = 1` is set in the arbiter's
own input (`d2j0_hotdisk/input_hotdisk:135`), in both Meier inputs
(`input_feet:181`, `input_feet_d2c:187`) — verified by grep — and now in the
item-0 probe. Under it the surface cell is `k == k_top_c` from the rem-based scan
and `surface_mf` is never read
([MMWSpalling.H:1501](../../../../src/Integrator/MMWSpalling.H#L1501)).

**Consequence to carry forward: per-cell removal is only safe under
`surface.follow_mask = 1`.** That is a real coupling, not a stylistic
preference, and the mechanism's parse must enforce it rather than trust the
input.

### Q2 — Does the energy ledger account correctly?

**Answer: yes, cell-wise and by construction, because removal writes no `H`.**

The ledger is `E_H = Σ H·V` over **all** cells, void included, with void `H`
frozen ([MMWSpalling.H:2554](../../../../src/Integrator/MMWSpalling.H#L2554)) —
a sum over cells, never over columns, so it has no top-only path to rely on.

The load-bearing fact is that removal deliberately does **not** reset the voided
cell's thermal state. [Removal.H:1447–1464](../../../../src/Integrator/MMWSpalling/Removal.H#L1447)
records the reasoning: the cell holds no material and is insulating (the
void-aware stencil carries zero flux across any void face and zeroes the void's
own conductivity), so its stored `T`/`H` is physically inert; an earlier reset to
`T_amb` produced a spurious 280 K checkerboard. Only `D` is zeroed, which the
ledger does not read.

So removing a cell changes no term in `E_H`, wherever that cell sits. **Removal
is exactly ledger-neutral at any topology** — the top-only assumption was never
in the ledger to begin with.

**Still owed, and not claimed here:** `ledger_err` at round-off on a run where
lateral removal *actually fires*. That is item 4. The item-0 probe supplies the
weaker interim datum — **max |ledger_err| = 1.788e-13 with an interior void
present**, which is round-off but is *not* the same statement as "closes when
lateral removal fires".

## 1. Item 0 — the verification scaffolding

| # | Deliverable | Where |
|---|---|---|
| 0.1 | `spall.seed_void`, default off, no output column | parse + `ApplySeedVoid` ([MMWSpalling.H:4678](../../../../src/Integrator/MMWSpalling.H#L4678)) |
| 0.2 | seeded-void probe at real call sites | `tests/MMWSpalling/unit/seed_void/` |
| 0.3 | beam invariant (a) restated as `n_above == 0` iff `k == k_top` | [MMWSpalling.H:3127](../../../../src/Integrator/MMWSpalling.H#L3127) |
| 0.4 | 2639/2639 with every new key off | _pending_ |

### The probe found two real defects. Item 0 paid for itself.

**Defect 1 — in my own code, caught by the guard rather than the clever checks.**
The first version parsed `queryarr("seed_void", …)`, but `ParseSpallSettings`
receives the **top-level** `ParmParse` and every key in it carries its own
`spall.` prefix. AMReX only *warns* on an unconsumed key ("Unused ParmParse
Variables"), so the run completed with the feature silently off — **and S1 and S2
both passed against a state that did not exist.** Only S3, which asserts the
seeded box really is void, caught it.

That is D2q-1 takeaway 4 again, from the other side: a test that asserts only the
interesting property will happily certify an empty experiment. **Every probe of
this kind needs a check that the state under test was actually constructed.**

**Defect 2 — a genuine accounting bug in `BeamEnergyBalance`, found by the probe
doing its job.** With the void seeded, `beam_closure_err` came out at
**1.429e-01**. Cause: `P_inc` was accumulated wherever `is_top_solid` held —
"the cell directly above is void, or this is the domain top" — which is a
**contiguity** test. With an interior void, **two** cells per column satisfy it:
the true top at k = 29, and the cell sitting just under the seeded void at
k = 19. The incident beam power was counted twice for every seeded column.

The arithmetic confirms the mechanism rather than merely fitting it: 20 of the
120 lit columns are seeded, so `P_inc` is inflated by 140/120, predicting a
closure error of (1 − 1.1667)/1.1667 = **−14.29 %** against the measured
**1.429e-01**.

Fixed by the same restatement invariant (a) needed — count at `k == k_top`, not
at "void above" ([MMWSpalling.H:3116](../../../../src/Integrator/MMWSpalling.H#L3116),
[:3158](../../../../src/Integrator/MMWSpalling.H#L3158)). Identical while void is
contiguous, so identity is preserved. After the fix
**`beam_closure_err` = 1.957e-15**.

**This is the finding that justifies item 0 existing.** The bug was latent in
D2q-1: shipped, byte-identical, and wrong the instant anything created the state
D2q-1 was built to allow. The mechanism in items 1–3 would have hit it
immediately, in a run where the closure error would have been one symptom among
many.

| probe check | before fix | after fix |
|---|---|---|
| S1 invariant (a) restated | pass | pass — 0 violations over 200 rows |
| S2 beam path counts solid cells | **FAIL 1.429e-01** | **pass 1.957e-15** |
| S3 seeded box is void | **FAIL (key never parsed)** | pass — 40 cells |
| S3b solid survives above | pass | pass — 8 layers, k = 22…29 |
| S4 `k_top` above the void | pass | pass — k_top = 29 |
| S5 undercut survives the run | fail | pass — 21 plotfiles |

**The seed writes `removed` only — `phi` is left saying "solid".** That
disagreement is deliberate: it is the state Q1 describes, and any consumer that
quietly fell back to the level set fails the probe immediately.

**Why invariant (a) had to be restated.** The old form tested *contiguity* —
"the cell directly above is void" iff `n_above == 0`. Each of the eight solid
cells stranded above the seeded void satisfies the left side and not the right,
so the old form flags them as violations the moment an undercut exists. The
restatement is exact at any topology, because `n_above` counts solid cells in
`(k, k_top]` and `k_top` is itself solid, so the count vanishes exactly at the
top cell. With void contiguous the two agree — which is what preserves identity.

## 2. Identity — the gate, with every key off

Batched bisectably (D2q-1 takeaway 1). Each row is a full `witness/post_all.sh`
sweep against `ref_post_d2e.md5`.

| batch | contents | runs | identity |
|---|---|---|---|
| item 0 | `seed_void`, invariant (a) restated, the `P_inc` fix | 11/11 rc=0 | **2639/2639 PASS** |
| items 5 + 7 | `per_cell_removal` and `side_face_h_max` keys | 11/11 rc=0 | **2639/2639 PASS** |
| **CHECKPOINT** (items 6 + 8) | criterion on face normals, `PerCellRemoval`, counters | 11/11 rc=0 | **2639/2639 PASS** |

**Process failure to record: one sweep was invalidated by my own error.** I
rebuilt while the items-6+8 sweep was still executing, swapping the binary
mid-sweep — the exact hazard written into this packet's guardrails. That run was
killed, its `mech_wit`/`mech.md5` artefacts deleted, orphaned processes checked,
and the gate re-run from scratch on the settled binary. **No result was reported
from the corrupted run.** The guardrail is correct and I broke it anyway; the
only defence that worked was noticing.

## 2b. The feature actually exercised, not merely proved inert

Identity proves the keys are inert when **off**. These checks exercise them
**on**. None is the arbiter; the campaign is deliberately not launched.

**Guards fire, with their reasons:**

| configuration | result |
|---|---|
| `per_cell_removal=1` without `side_face_flux=1` | **aborts** ([:693](../../../../src/Integrator/MMWSpalling.H#L693)) — "the criterion evaluated there reads bulk rock and the mechanism is silently inert" |
| `per_cell_removal=1` + `connected_cluster` + `A_crit>0` | **aborts** ([:715](../../../../src/Integrator/MMWSpalling.H#L715)) — "the 4-connected cluster filter is defined over the top-z plane and has no lateral analogue" |

**Key on, arbiter C1 keys, 90 s (2 mm):** runs clean; `pcr_overhang_cells` and
`pcr_lateral_cells` are registered **after** `side_*`, so every pre-existing
`thermo.dat` column keeps its position.

| quantity | value |
|---|---|
| `side_faces` | **323** — real exposed lateral faces exist |
| `pcr_lateral_cells` | 0 |
| `pcr_overhang_cells` | 0 |
| `ledger_err` | max 2.009e-15 |
| `beam_invariant_viol` | 0 |
| `beam_closure_err` | 0 |

**A zero firing count is ambiguous, so it was resolved rather than reported.**
Reading the last plotfile directly:

| | |
|---|---|
| candidate cells (solid, lateral face exposed, **not** column top) | **300** |
| of those carrying a written `Sp` | **299** |
| `Sp` on candidates | max **0.4533**, mean 0.1536 |
| `T` on candidates | max **656.6 K** (firing point ≈ 821 K) |
| `Sp` at column tops | max 0.9964 |

So **the criterion is live on lateral faces and correctly declines to fire** —
the wall is 165 K below threshold, which is precisely what the 80 W/m²K cap
exists to enforce (it sits below `h*`, so the wall provably cannot heat itself
into spalling). Whether the **corner drain** carries it over is the arbiter's
question at 250 s and is **not claimed here**.

*One candidate of 300 carries `Sp = 0` exactly. Not chased at the checkpoint;
flagged for the arbiter's diagnostics.*

## 3. The mechanism (items 5–8)

| # | Deliverable | Where |
|---|---|---|
| 5 | `spall.per_cell_removal`, default 0, preconditions enforced | parse [:3901](../../../../src/Integrator/MMWSpalling.H#L3901), cross-block checks [:685](../../../../src/Integrator/MMWSpalling.H#L685) |
| 6 | `Sp >= 1` on the face normal, at every exposed cell | sampler [:4268](../../../../src/Integrator/MMWSpalling.H#L4268), exposure [:4218](../../../../src/Integrator/MMWSpalling.H#L4218) |
| 6 | the lateral removal pass | `PerCellRemoval` [:4820](../../../../src/Integrator/MMWSpalling.H#L4820), called [:998](../../../../src/Integrator/MMWSpalling.H#L998) |
| 7 | `surface_patch.side_face_h_max`, default 80 W/m²K | parse [:5663](../../../../src/Integrator/MMWSpalling.H#L5663), capped in kernel **and** ledger |
| 8 | `pcr_overhang_cells`, `pcr_lateral_cells` | registered only with the key on [:709](../../../../src/Integrator/MMWSpalling.H#L709) |

### The gate is the same gate — verified, not assumed

The packet's claim was re-checked with my own greps before relying on it:
`weibull_detach_mode` defaults to **`PerFace`**
([:6420](../../../../src/Integrator/MMWSpalling.H#L6420)); `weibull_A_crit`
defaults to 0 and is `query_required` only under `ConnectedCluster`
([:6618](../../../../src/Integrator/MMWSpalling.H#L6618)); the cluster area
filter needs **both** ([:4462](../../../../src/Integrator/MMWSpalling.H#L4462));
and **none** of `d2j0_hotdisk/input_hotdisk`, `input_feet` or `input_feet_d2c`
sets either. So `clu(i,j,k_top) > 0.0` reduces to `Sp_top >= 1` in every run on
this line, and `Sp >= 1` along the face normal is the *same* rule, not an
approximation of it.

**The divergent case aborts rather than falling back.** `per_cell_removal = 1`
with `connected_cluster` *and* `A_crit > 0` is refused with the reason: a
4-connected labeller over the top-z plane has no lateral analogue. Generalising
it is a separate packet.

### Scope choice that makes criterion (f) exact

`PerCellRemoval` considers **only** cells that are not their column's top and
that have an exposed **lateral** face. The column top keeps its existing path
entirely — depth scan, flake coherence, `h_col` remainder. A **flat interior
therefore has no candidate cells at all**, so with the key on its rate is not
merely *close* to the key-off rate, it is produced by the identical code.
Criterion (f) then measures a real difference rather than numerical drift.

### How a non-top cell is voided MPI-consistently

This has no precedent in the codebase, so it is set out in full.

**The criterion is local.** `Sp` is computed in `UpdateSpAfterMechanics` from
the cell's own `T` and `a_f`. The `K_I` integral runs over the flaw length
`a0 = 20 µm` against `dx = 2 mm` — a ratio of 100 — so `floor(dv/dz) = 0` for
every sample and each one lands in the scored cell. (That is the same fact S1
recorded as "the criterion is a cell-average `T` threshold".) No column gather
is required, unlike the top-cell path, where a column spans ranks and `k_top`
must be `ReduceIntMax`-ed.

**The only non-local read is the six-neighbour void test**, taken from
`removed_mf`'s one-ghost layer, refilled by `Util::RealFillBoundary` inside
`UpdateSpAfterMechanics` before the scoring loop.

**Decide, then flip.** `Sp` is written from the frozen pre-pass state; the
removal loop reads only `Sp` at an owned cell when deciding. So no cell's fate
depends on whether a neighbour was flipped first, and the outcome is independent
of box decomposition and of iteration order within a box. Flipping in place
while re-testing exposure would have made a rank's sweep order observable — the
same class of defect as the 0-ghost `removed_mf` seam bug recorded earlier in
this project. After the flip the ghost layer is refilled and both counters are
`ReduceRealSum`-ed, so every rank reports the same census.

**The voided cell's thermal state is left alone**, matching the column path
([Removal.H:1447](../../../../src/Integrator/MMWSpalling/Removal.H#L1447)) — which
is also what keeps the ledger neutral (Q2).

## 4. Regressions — 8/8 PASS

| test | result |
|---|---|
| `unit/spall_event` | PASS |
| `unit/sp_lefm` | PASS |
| `unit/regime_low_high_power` | PASS |
| `unit/amr_microstructure_regrid` | PASS |
| `validation/hu/hu_end_to_end` | PASS |
| `sp_meier_pilot/test_feet` | PASS (check-only) |
| `sp_meier_pilot/test_feet_d2c` | PASS (check-only) |
| `unit/seed_void` (new) | PASS 6/6 |

The two `test_feet*` scripts re-assert recorded outcomes — including the
standing D2g/D2i "expected-fail: mesh 2 → 1 mm" verdicts — rather than running a
simulation. Their PASS means "nothing recorded has moved", **not** "the Meier
gate passes".

## 5. Arbiter — NOT LAUNCHED

Stopped at the packet's checkpoint. Item 4 in full remains: `CRITERION.md` with
hashed predictions, `run.py` / `analyze.py`, `C1_2mm` / `C1_1mm` / `C0_2mm` /
`C0_1mm` at 250 s, the six criteria, and `C1_0p5mm` only on a pass.

**What the checkpoint evidence does and does not support.** It supports: the
plumbing is inert with the keys off (three sweeps, 2639/2639); the guards fire;
the criterion is evaluated on lateral faces and correctly declines to fire at
656 K. It does **not** support any claim that the corner drain will spall the
wall — `pcr_lateral_cells` is 0 at 90 s and nothing here shows what happens at
250 s. That is exactly what the arbiter is for.

**Carried into the arbiter's diagnostics:**
- the one candidate cell of 300 carrying `Sp = 0` exactly;
- `ledger_err` on a run where lateral removal *actually fires* — still owed,
  and the weaker "void merely present" datum must not be quoted for it;
- the overhang census, where a deep undercut is a **finding**, not a bug.

## 5. Recorded gaps

- **The regrid-repair site is exercised by no run in the identity set** (it is
  single-level). `unit/amr_microstructure_regrid` is **not** a check on it — that
  test has no removal. The site remains unverified. AMR is parked, so this is a
  recorded gap rather than a task, and it must not be described as covered.
