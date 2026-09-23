# D2q-2c — criteria and predictions
# Written and hashed BEFORE any run, before `output/` existed. NEVER edit.

Binary `bin/mmwspalling-3d-g++` =
`20ab657331d72dd25b4b34eddb8698295c90306c63e9400b4b2569adebab18fc`.
**Provenance only, not an oracle** — unchanged source rebuilds to different
hashes on this machine (D2q-1 takeaway 3).

Identity with every key off: **2639/2639 PASS**, run through
`guarded_post_all.sh`, which confirmed the binary did not move mid-sweep.

## 0. Every reference value reproduced with my own estimator first

Per the guardrail, before trusting anything the packet states:

| packet states | my estimator | verdict |
|---|---|---|
| ρc = 2.1725e6, α = 6.9045e-7 m²/s | 2.1725e6, 6.90449e-7 | ✓ |
| κ/v = 1.604 mm on the prescribed feed | 1.6036 mm | ✓ |
| `h*` = 88.5 W/m²K reaches 821 K in the 121 s dwell at T_gas = 1600 | **h = 88.5 → 822 K at 121.9 s** | ✓ |
| cap 80 W/m²K gives 786.9 K in 121 s | **786.9 K** | ✓ |
| probe: Sp ≤ 0.4533 at T ≤ 656.6 K on lateral faces | 0.4533 / 656.6 K | ✓ |

The `h*` row is the strongest check: I solved the semi-infinite Robin relation
`T_s = T_0 + (T_g − T_0)[1 − erfcx(β)]`, `β = h√(αt)/k`, for the dwell at which
h = 88.5 reaches firing, and got **121.9 s** against the packet's stated 121 s
dwell — reproducing the packet's constant without being given it.

## 1. THE STOP CHECK (pre-flight 1) — result: PROCEED

**At what temperature does the lateral `Sp` reach 1?**

Measured from the 90 s probe, not extrapolated: the 12 cells whose `Sp` is
nearest 1 sit at **T = 822.2 ± 9.5 K**. That matches the ~821 K at which the
*vertical* criterion fires, which it must — at `a0 = 20 µm` against `dx = 2 mm`
every `K_I` sample lands in the scored cell, so both orientations reduce to the
same local function of `T`.

*(A quadratic fit over all 1899 scored cells extrapolates to 858.9 K. The
822.2 K figure is interpolation at Sp ≈ 1 and is the one used. The discrepancy
is recorded, not hidden.)*

**Is that reachable under the cap?**

| h [W/m²K] | exposure to reach 822.2 K |
|---|---|
| 80 (the cap) | **149.1 s** |
| 88.5 (`h*`) | 121.9 s |

The arbiter runs **250 s**, so firing is **not excluded by construction** and
the campaign is worth running. **But 149.1 s is *continuous* exposure at the
full cap coefficient, and a lateral face only becomes exposed partway through
the run** — so the available exposure is well under 250 s, and the corner drain
must make up the difference. That is precisely the open question, and it is why
P4 below is a live prediction rather than a formality.

## 2. Pre-flight 2 — the contiguity grep

Searched `MMWSpalling.H` and `Removal.H` for: every `rem(...k+1)` read; every
"cell above" / "above is void" / `is_top_solid` phrasing in code and comments;
and counted the correct `k == k_top` idiom.

**Found: no third contiguity bug.** Two live `rem(k+1)` reads remain, both
audited:

| site | use | verdict |
|---|---|---|
| `MMWSpalling.H:4327` | "does this cell have an exposed **top face**" — picks the scoring normal | correct question, not a column-top proxy |
| `MMWSpalling.H:4934` | "void here, solid above" | this **is** the undercut definition |

One **stale comment** was found and corrected: `BeamEnergyBalance`'s header
still documented invariant (a) in its old contiguity form and invariant (b) as
live after D2q-1 retired it. That is exactly the note that would have led the
next reader to reintroduce the bug.

## 3. Cases

`C1_2mm`, `C1_1mm`, `C0_2mm`, `C0_1mm`, all with
`side_face_flux = 1`, `side_face_h_max = 80.0`, `per_cell_removal = 1`.
**The pair decides; a failure there is decisive and `C1_0p5mm` is not run.**

**Key-off anchors, already on disk, never re-run:** `d2j0b_plane/output/`
(C0/C1 at 2 and 1 mm) and `d2l_scan/output/` (`A_f02`, `A_f02_1mm`, `A_f01`,
`A_f01_1mm`) as the attribution baseline, so a pass belongs to per-cell removal
and not to side heating.

## 4. Metrics — `d2j0b_plane/CRITERION.md` §Metrics, reused UNCHANGED

Disk-mean rate over columns with r < 30 mm from `<pf>_removal_events.csv` per
50 s block; `t_stop`, `t_cross`, `Δ_outer`, `r_front`; 2 mm bands 18–30 mm;
`v_c` over r < 6 mm. Same estimator as the Meier gate, so a pass transfers.

## 5. PASS requires all six

- **(a)** disk-mean gap within **5 %** in every block 50–100, 100–150,
  150–200 s, for both mesh pairs.
- **(b)** `Δ_outer` mesh-independent within 20 %, **or no band stalls**.
- **(c)** no inward-walking `r_front` at the finer meshes.
- **(d) INTEGRATED, never pointwise.** The corner flux is *physically* singular
  (reflex corner, `r^(-1/3)` for an isotherm boundary) and will not converge
  however correct the model is. Gate on **removed volume per unit perimeter
  across the corner** and on the **disk-mean rate**, ratio ≤ 1.25 per mesh
  halving. `k·ΔT/dx` is reported as a **diagnostic only**, against D2o-0b's
  302 kW/m² / 58 %.
- **(e) REMOVAL CEILING, anchored to the ARTEFACT-FREE control.** Disk-mean rate
  must not exceed **`C0` KEY-OFF at the same mesh by more than 15 %**.
  Pre-registered anchor, with its numbers fixed here:

  | block | C0_2mm key-off | ceiling (×1.15) | C0_1mm key-off | ceiling (×1.15) |
  |---|---|---|---|---|
  | 50–100 | 1.211 | **1.393** | 1.085 | **1.248** |
  | 100–150 | 1.119 | **1.287** | 0.948 | **1.090** |
  | 150–200 | 1.048 | **1.205** | 0.838 | **0.964** |

  **Anchored to C0 key-off, never to a cascade-suppressed C1 run** — that made
  (a) and (e) jointly unsatisfiable and killed D2p-1. This is also the packet's
  own worked example: D2k's 1 mm block-1 rate of 1.651 is +52 % over C0's 1.085
  and is rejected. Any block where C1 exceeds C0 at all must be explained even
  if it passes: the exclusion plane can only *remove* heating.
- **(f) EXACT where it can be.** `PerCellRemoval` touches only non-top,
  laterally-exposed cells, so **until `pcr_lateral_cells` first goes non-zero
  the key-on run must be BYTE-IDENTICAL to the key-off run** — asserted
  directly on the `removal_events.csv`, not to a tolerance. After the first
  lateral firing, the flat-interior disk-mean rate must still match key-off
  within ≈ 1–2 %. **If (f) fails, nothing else here is scoreable.**

## 6. Predictions

- **P1 — lateral over-width 2–8 mm per side** (Meier: skirt + 2.5–7.5).
- **P2 — wall 600–800 K** where the side coefficient sits at the cap: altered,
  not spalled.
- **P3 — the drain becomes removal rather than loss.** Report the fraction of
  the drain appearing as removed volume against stored enthalpy, plus
  `pcr_lateral_cells` and `pcr_overhang_cells` per block and `ledger_err`.
- **P4 — THE INERT BRANCH IS A REAL OUTCOME, and I expect it.** My own reading:
  `Sp` maxes at 0.4533 at 90 s; the cap needs 149.1 s of *continuous* exposure
  to reach firing; a lateral face is exposed for well under that within a 250 s
  run. **So I predict `pcr_lateral_cells` stays 0 or near 0 for the whole
  campaign**, and that (a)–(e) therefore read close to the key-off C1 cascade.
  If so the result is **negative and honest** — the drain does not carry the
  wall to firing under a physically-derived cap — and **the cap must NOT be
  raised to rescue it.**

**What would refute me:** any block with `pcr_lateral_cells > 0` and a
disk-mean gap inside 5 %.

## 7. Decision, pre-registered

| outcome | reading |
|---|---|
| **(a)–(f) all hold** | the representation change worked; next is the Meier gas closure, scored on the **finished** wall only |
| **(f) fails** | the per-cell path is not the existing physics → **stop** |
| **(e) fails** | it converged by drilling faster → **stop**, do not touch the cap |
| **inert (P4)** | report it; the next move is a gas closure that raises the wall temperature **honestly**, not a bigger cap |

Nothing re-windowed, re-binned or re-chosen after results are seen. Neither
`side_face_h_max` nor `side_face_factor` is ever adjusted.
