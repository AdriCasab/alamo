# D2q-1 results — make `removed_mf` authoritative (stage 1, pure refactor)

**Status: IN PROGRESS.** Gate = 2639/2639 byte-identity and nothing else.

## 0. Provenance

| item | value |
|---|---|
| binary **before** (reverted D2p-1 tree, rebuilt) | `83603563b1e138c489aa8a94c4db04056fd51af72aa9f11642913bca6fa288b5` |
| baseline sweep on that binary | **2639/2639 byte-identical to `ref_post_d2e.md5`** |
| binary **after** batch 1 | `dc9c516b1cd18fa998ba7fead23b49175b2a5dae5b55abe373aa58c4b019f4d2` |
| binary **after** (final) | _pending_ |

**The binary sha256 is provenance, not an oracle.** Rebuilding the *unchanged,
git-clean* source produced `83603563…` where the D2o-0b campaign recorded
`62dea451…`. Same source, different hash ⇒ the link step is not reproducible on
this machine (LTO + parallel link). The acceptance criterion is therefore the
**2639 output files**, never the binary hash. Recorded here because this packet's
whole gate is identity, so it matters that one plausible identity check is
useless.

## 1. The four sites, and their treatment

`removed_mf` becomes the sole authority on "is this cell void". `phi` is demoted
from void oracle to **sub-cell recession accumulator** — it still carries the
fractional remainder that decides *when* the next cell flips, but no decision
path asks it *whether* a cell is void.

| # | Site | Old assumption | New handling | Why byte-identical while void is contiguous |
|---|---|---|---|---|
| 1 | `ColumnTopSolid` (`MMWSpalling.H:2878`) | none — already scans `rem` alone | **unchanged** | — |
| 2 | `RepairRemovalStateAfterRegrid` (`:2769–2830`) | rebuilds a column-linear `phi`, then sets `rem = (expect < 0)` | `rem` is preserved and only *cleaned* to an exact 0/1 by the same `> 0.5` test every reader uses; `phi` is re-imposed as the linear accumulator anchored on `rem`'s own `k_top` | `expect < 0` and `rem > 0.5` agree cell-for-cell while void is contiguous |
| 3 | Beam path integral (`:1589`, `:2960`, `SolidFacePath` `:2737`) | `n_above = k_top − k`, i.e. geometric distance | count **solid cells** above `k` | equal when every cell in `(k, k_top]` is solid |
| 4 | Beam-closure invariant (b) (`:2971`) | asserts `phi − n·dz ∈ [0, dz)` — a phi/mask *sync* check | retired; contiguity check (a) restated on `rem` | `n_viol = 0` under both forms today, and `n_viol` feeds `beam_invariant_viol` in `thermo.dat` |
| 5 | Depth scan (`Removal.H` ~`:640–976`) | strides down from `k_top` by cell index, every cell solid | stops at the first void below `k_top` | no interior void exists |
| 6 | `k_top` decision scans (`:2780`, `:3834`, `Removal.H:240`) | `rem <= 0.5 && phi >= 0.0` | `rem <= 0.5` | the `phi` conjunct is redundant |
| 7 | `UpdateSurfaceMaskFromPhi` (`:4452`) | mask = `0 <= phi < dz` | **NOT CONVERTED — ESCALATED** | it is not byte-identical; see below |
| 8 | `EnergyLedger` (`:2537+`) | `E_H = Σ H·V` over **all** cells, void included, `H` frozen | **no change needed** — already cell-wise, never column-wise | — |

### Site 7 — THE ESCALATION. Not converted; left on `phi`.

**This is the one site that cannot be made `rem`-authoritative without changing
behaviour, so per the packet's guardrail it was reverted and escalated rather
than approximated.**

The pre-check said it was safe and the pre-check was wrong, which is worth
recording. The code's own comment documents a known defect — `surface_missing_cols`
counts columns with a mask-top solid cell but **no** `surface_mf` cell, "the
phi-window dropout at `phi_top ~ dz` ties". I swept every `thermo.dat` in the
tree: that counter is **0 everywhere**, so I predicted the conversion was
identity-preserving. **The gate refuted it.**

**Measured:** replacing `0 <= phi < dz` with "solid here, void (or domain top)
directly above" changed **24 files, all in `validation/meier/sp_meier_pilot/output/dev2d`**.

**It is not a cosmetic difference in a written field — the trajectory diverged.**
The differing plotfiles *migrate between* the `40000cell`, `48000cell` and
`56000cell` output directories, i.e. the number of removed cells at a given plot
time changed. `input_2d_dev` sets no `surface.follow_mask`, so `surface_mf` is a
**live decision path** feeding `SurfaceCellFlux`, not a diagnostic: change which
cell is "the surface" and you change where the heat goes.

**Why the pre-check missed it:** `surface_missing_cols` only counts the
*dropout* direction (a top solid cell with no surface cell). It does not count
the opposite error — a surface cell that is **not** the top solid cell. At a
`phi_top ≈ dz` tie the window `[0, dz)` can select the **void** cell above the
top solid cell, which the counter is blind to. That is consistent with the E1
finding already in the record ("the `floor(phi/dz)` form failed at round-off
ties, `phi = dz` exactly, which is common"). **This mechanism is a hypothesis
fitted to the evidence, not something I instrumented and confirmed** — the
established facts are the 24 files and the directory migration.

**Consequence for stage 2.** `surface_mf` remains a phi-derived decision path.
Stage 2 must either run with `surface.follow_mask = 1` (where `surf_c` is
computed from `k_top` and this mask is bypassed entirely), or carry a deliberate,
separately-gated change to this mask — it cannot be folded into an
identity-gated refactor, because the correct behaviour and the current behaviour
differ. **A future planner should treat this as a latent bug with a known
reproducer**, not as a limitation of the refactor.

## 2. `phi_mf`: retained, demoted

Retained as the sub-cell accumulator and as a plotfile diagnostic. It is
**written** by the removal producers and by the regrid repair; it is **not read
by any decision path** after this step. The two remaining `phi` tests
(`Removal.H:200`, `:204`) are inside the `foot_body_clearance` producer, where
they *maintain* the phi/rem sync rather than consult it.

## 3. Identity — THE GATE

Run in three batches so a failure bisects to one group rather than to "the
refactor". Each row is a full `witness/post_all.sh` sweep against
`ref_post_d2e.md5`.

| batch | sites | binary | runs | identity |
|---|---|---|---|---|
| baseline (no edits) | — | `83603563…` | 11/11 rc=0 | **2639/2639 PASS** |
| 1 | 2, 6 — regrid repair + the three `k_top` scans | `dc9c516b…` | 11/11 rc=0 | **2639/2639 PASS** |
| 2a | 4, 7 | `9f0cb3cf…` | 11/11 rc=0 | **FAIL — 24 files, all `dev2d`** |
| 2b | 1–6, 8 (site 7 reverted; adds the beam path + depth-scan gathers) | `7fe8d434…` | 11/11 rc=0 | **2639/2639 PASS** |

**Final state: sites 1–6 and 8 shipped, byte-identical. Site 7 reverted and
escalated** (above). Batch 2a is the only failure and it isolates cleanly: site 4
can only move `n_viol` into `thermo.dat`, while every differing file was a
`dev2d` plotfile, so site 7 is the cause — confirmed by 2b passing with site 4
still in.

**Note on the binary hash as an oracle: it is not one.** Three passing sweeps
produced three different binaries, and rebuilding *unchanged* source changed the
hash too (§0). Only the 2639 output files decide.

## 4. Regressions — 7/7 PASS

`unit/spall_event`, `unit/sp_lefm`, `unit/regime_low_high_power`,
`unit/amr_microstructure_regrid`, `validation/hu/hu_end_to_end`,
`sp_meier_pilot/test_feet`, `test_feet_d2c` — all rc=0 and all reporting PASS.

`amr_microstructure_regrid` is the one that matters most here: it is the direct
check on site 2, the only site where `removed_mf` has to survive interpolation
rather than merely be read in place of `phi`.

The two `test_feet*` scripts are **check-only** — they re-assert recorded
outcomes (including the standing D2g/D2i "expected-fail: mesh 2 → 1 mm"
verdicts) rather than running a simulation, so their PASS means "nothing
recorded has moved", not "the Meier gate passes".

**Energy ledger:** `max |ledger_err|` = 1.6e-15 … 5.3e-15 across the
`robin_feet` runs (the `sp_meier_pilot/input_feet` keys at short stop times).
Round-off, and these `thermo.dat` files are inside the byte-identical set, so
the value is unchanged from the pre-change binary — which is what item 4 asked
for. Note `dev2d` contributes **no** `thermo.dat` to the hashed set (its 24
witness files are all `Cell_D_*`), so the ledger evidence comes from
`robin_feet`, not from `dev2d`.

## 5. What stage 2 can now do that it could not before

**It can leave solid above void in one column — an undercut.**

Before this step, `rem` was reconstructed from a column-linear `phi` on every
regrid (`rem = (expect < 0)`), so an interior void was erased at the next
regrid even if something had created one. Every column-top scan additionally
required `phi >= 0`, so a cell below a void would not have been found as a
surface. The beam measured its optical path as geometric distance `k_top - k`,
so it would have attenuated through a void as if it were rock; and the K_I depth
scan strode downward reading every cell as material.

After it: `removed_mf` is authoritative and survives regrid (it is only cleaned
to an exact 0/1, never re-derived); every `k_top` decision consults `rem` alone;
the beam counts **solid cells**; and the depth scan is clamped to the contiguous
solid run below the top. So when stage 2 flips a cell to void at an exposed
lateral face — the wall column heated **at its base** by the corner drain, which
is the configuration D2p-0 proved the height function could not express — the
state survives, the surface is found beneath it, the beam does not heat through
it, and the crack scan stops at it.

**What stage 2 still has to do**, and must not assume is done:
1. **`surface_mf` is still phi-derived** (site 7). Either run with
   `surface.follow_mask = 1`, which bypasses that mask entirely, or change it
   under its own gate — it cannot ride along inside an identity-gated refactor.
2. **Beam invariant (a) stops being an invariant.** It asserts contiguity, which
   an undercut violates by construction. It is untouched here because nothing
   yet creates one; stage 2 must restate it or it will count violations against
   `beam_invariant_viol` forever.
3. **Nothing creates an interior void yet**, so none of the new handling has been
   exercised against a real one — only proved inert. See the note on the unit
   test in `ACTIVE_STEP.md`.
