> ⚠ **SUPERSEDED IN PART — this is revision 2's record, which review
> rejected.** The physics it reports stands; the numbers were re-measured on
> the revision-5 binary and (e), (d) and the estimator were fixed. Read
> **`RESULTS_REV5.md`** first. Two claims here are now known to be wrong:
> §5b's `surf_side_share` = 0.618 is the END-OF-RUN figure (0.470 over the
> rate window), and (a)'s 2.57 % came from an estimator whose answer moves
> 0.29–7.77 % with the window edges.

# P0b — ONE surface temperature: the rate is mesh-independent, with two open faults

**The 1-D steady rate's mesh spread is 2.57 %**, against a 10 % bar, **86.2 %
after P0** and **7.1 % pre-P0**. All four scored quantities that depend on the
new physics pass: the spread (a), the scheme's own energy balance (c′), the
domain ledger with the sink (c″), and byte-identity with the keys absent (c).

⚠⚠ **Two criteria do not pass, and neither is cosmetic:**
* **(e) decomposition FAILS** — the `flake`/remap path is not layout-invariant.
  Isolated: keys off and reconstruction-only are bit-identical; the dependence
  needs `removal_energy = flake`. See `DECOMPOSITION.txt`.
* **(d) is split** — the Kant anchor holds at `h_conv = 150` (0.169 %, bar 0.2 %)
  but moves to **2.089 %** at `h_conv = 500`, because deliverable 5.1's τ = 1 s
  smoothing is *longer than the onset time* there (0.3–0.4 s).

## 0. Provenance

* Binary **`cd764eec…`** (P0 was `f65d10eb…`). Tree gate vs `witness/d2z1.md5`
  **2639/2639, zero differing** on all FOUR binaries built during this packet,
  `p0b` … `p0b4`, the last being the final `cd764eec…`. Each gate used `--force`, and each was
  re-run from scratch after every rebuild — never reused.
* Nothing committed.

## 1. (a) — the rate, and the controls that make it mean something

`validation/p0_steady_column`, 8 × 8 × 120 mm, `h_conv` 600 / `T_flame` 1673.15
(0.55 MW/m² at the firing surface), steady window 20–55 s.

| leg | 2.0 mm | 1.0 mm | 0.5 mm | **spread** |
|---|---|---|---|---|
| **`surface`, c_f = 0.4** | 2.544 | 2.555 | 2.490 | **2.57 %  PASS** |
| `cell` | 3.139 | 2.883 | 2.729 | 14.06 % |
| `face` | 2.189 | 2.379 | 2.482 | 12.47 % |
| `pinned` | 3.039 | 2.791 | 2.602 | 15.56 % |
| `surface`, c_f = 0.3 | 2.414 | 2.405 | 2.392 | 0.91 % |
| `surface`, c_f = 0.5 | 2.761 | 2.670 | 2.661 | 3.71 % |

✅ **The controls still show the spread** (12–16 %), so `surface` is not falling
through to an existing form — the failure mode §14.4 of the brief warned about.
Making it a **fourth** form rather than a branch selected by `sp_surface_T_apply`
is what preserved that comparison.

⚠ **Score the spread, never the level.** The ~1.2× offset above the continuum
ablation limit is physics: a thick flake leaves cooler than the surface.

## 2. (c′) and (c″) — the scheme is energetically right, not merely flat

**(c′)** `v·ρc·(T̄_f − T0)` against the ledger's net flux: **−0.15 / +1.37 /
−1.09 %**, bar 5 %. **PASS.** Spread alone would have passed a mesh-independent
but energetically wrong scheme; this is what rules that out.

**(c″)** ⚠⚠ **and it exposed a defect in the ledger itself.** `ledger_E_H` sums
H over **every** cell including void ones — whose H is deliberately frozen, not
reset — so **removal has never been a sink in the ledger's eyes**. Measured:
`dE_H` = 1363.9 J against `E_src` = 2264.3 J and `E_flake` = 2193.0 J, which
cannot be reconciled in any combination. Worse, on a *voiding* step the remap's
index shift slides the profile down one cell while the voided cell keeps its own
H, so the all-cell sum gains a spurious cell of enthalpy per void.

A solid-only sum (`ledger_E_H_solid`, new) removes both effects, and then the
flake energy is the sink **exactly**:

| mesh | `ledger_err_solid` | all-cell `ledger_err` |
|---|---|---|
| 2.0 mm | **5.99e-4** | −0.397 |
| 1.0 mm | **6.15e-4** | −0.192 |
| 0.5 mm | **3.65e-4** | −0.043 |

Bar 5e-3. **PASS 3/3.** The all-cell column shrinking with `dz` is the same
defect seen from the other side: less is removed per void at finer mesh.

⚠ **(b), the power identity**, is covered for the kernel↔ledger pair *by* (c″):
if deliverable 1 had reached the kernel and not the ledger, the closure could
not be 6e-4. The **jet-march** callers (`:2942`, `:2976`) carry no flux in a
1-D column and are exercised only by the Meier rerun.

## 3. The pre-registered predictions

| | prediction | measured |
|---|---|---|
| **P6** | `T̄_f/(T_surf−T0)` = `(1−e^-cf)/cf` = 0.8639 / 0.8242 / 0.7869 | **0.8639 / 0.8242 / 0.7869**, mesh-independent to 4 dp. **EXACT** |
| **P2** | level vs continuum = 1.157 / 1.213 / 1.271 | 1.166 / 1.212 / 1.299 ✅ (revision 1's 1.13/1.18/1.23 would have been wrong) |
| **P3** | `ell_eff` continuous across firings | median \|Δℓ/ℓ\| **3.0e-4 / 8.2e-3 / 6.8e-3** ✅ — no jump of order e^{c_f} |
| **P1** | (a) ≤ 1 % | **2.57 %** — above, but matching none of the failure ladder (20 / 70 / 86 %). At c_f = 0.3 it is **0.91 %**, so the residual is c_f-dependent |
| **P7** | firing window narrows | pending the Meier rerun |

## 4. ⚠⚠ Two faults found by running it

### 4.1 `h_f = c_f·ell_eff` is UNBOUNDED ABOVE — caught by the assertion

`ell_eff` is clamped from below (≥ dz/4) but **not from above**, so a
nearly-flat column — `x` small but not zero — gives an arbitrarily thick flake.
The Meier rerun aborted at **t = 129.4 s** with column (25, 13) trying to void
**287 cells (574 mm, more than the whole 500 mm block)** off a two-cell ratio of
**1.0014**. The restated "voids == crossings" assertion caught it rather than
letting half the domain disappear.

**Fixed** by capping `h_f` at **`spall.h_spall_max`** (0.015 m), the key that
already defines how deep this model looks for a flake at all. ⚠ **That is not
the `dz` clamp the packet forbids**: `dz` is a *mesh* quantity and clamping
there is exactly what made `T̄` mesh-dependent in P0; `h_spall_max` is
*physical* and mesh-independent, so the cap cannot reintroduce mesh dependence.
Verified inert where it should be: **0 of 88 firings capped** on the 2 mm
column, rate unchanged (2.543 vs 2.544 m/h). Bind fraction is reported as
`recon_hf_capped` and flag bit 8.

### 4.2 `pinned_idle_cycles` aborts the run under the new form

Deliverable 1.4 asks for it to be logged as inert. ⚠ **Not reading it is not an
option**: `ParmParse::queryclass` strict-checks unused keys, and **every frozen
leg from `d2v` to `d2ad` sets `pinned_idle_cycles = 1.0e9`** — so the first
Meier rerun of this packet aborted at parse before a single step ran. The key is
now read under `surface` too, and the "inert" message is what deliverable 1.4
actually asked for.

## 5. ⚠⚠ (e) FAILS — and it needs `removal_energy = flake`

Full isolation in `DECOMPOSITION.txt`. In one line: keys off ✅ bit-identical,
reconstruction-only ✅ bit-identical, **+ flake/remap ✗**. The first firing is
bit-identical between layouts; the field then diverges at **3.9e-11 relative**
in `T_top`, and the `Sp ≥ 1` threshold amplifies that into a 3.6 % difference in
firing count. The layouts that agree keep x and y whole; the one that differs is
the only one that splits them.

⚠ **Do not assume it is inherited.** P0's decomposition witness ran the Kant
case, which uses `removal_energy = cell` — so **P0 never tested the flake path
for layout invariance**, and this may predate P0b's remap. Telling them apart
needs a build with the remap off and P0's sub-cell `h_c` on, which was not run.

## 5b. (f) — the Meier `FH` rerun. REPORTED, NOT SCORED.

D2ac's frozen `FH` argv verbatim plus `robin_form = surface`,
`surface_T = reconstructed`, `removal_energy = flake`, `flake_c_f = 0.4` and the
P0 anchor (`a0_ig = 8.0e-6`, `a0_gb = 40.0e-6`). Ran the full 500 s, `rc = 0`.

| | frozen (pre-P0) | + P0 | **+ P0b** |
|---|---|---|---|
| axis ROP, 60–160 s | 2.766 m/h | 3.690 (+33.4 %) | **3.233 m/h (+16.9 %)** |
| axis depth at ≈499 s | 332.0 mm | 494.0 (at the base) | **392.5 mm** |
| firing window p10–p90 | 19.4 K (`T_top`) | 86.9 K | **15.8 K** |
| `flat` fallback at firings | — | 12.8 % | **0.0 %** |
| `T_surf` − `T_top` at firing | — | +73.0 K | **+188.9 K** |

✅ **P4 HELD** — the rate comes *down* from P0's +33.4 %, as predicted, because
the BC now drives on ~762 K instead of the cell average and removal costs
enthalpy. ✅ **P7 HELD** — the window narrows from 86.9 K to **15.8 K**, tighter
even than the frozen pre-P0 19.4 K. ✅ **P5 HELD** — deliverable 1.3's hold
takes the `flat` fraction at firings to **0.0 %**.

⚠⚠ **SIDE FACES CARRY 62 % OF `P_robin` ON MEIER**, and deliverable 1.2 leaves
them on the **cell average** (`surf_side_share` = 0.618). This is the single
biggest caveat on the Meier number: the top face has one surface temperature and
the walls have another, so *most* of the delivered power on this configuration
has not been converted. The 1-D column has no side faces, which is why this had
to be measured rather than assumed — it is exactly what deliverable 1.2 asked
for and it is not a small correction.

⚠ **The `h_f` cap binds on 85 % of COLUMNS but 0 % of FIRINGS**
(`recon_hf_capped` ≈ 3060 of 3600). It catches the cold, nearly-flat columns
outside the heated patch — the pathological population — and never a column that
actually fires. That is the cap behaving exactly as intended.

⚠ **The ledger does NOT close on Meier**: `ledger_err_solid` = −0.160,
`ledger_err` = −0.083, against 6e-4 for the solid form on the 1-D column. Meier
adds the jet-enthalpy march, side-face flux and a beam term to `E_src`, and the
flake sink is only wired for the column path. **(c″) is scored on the 1-D
column, where the configuration is clean**; the Meier residual is an open item
for whoever wires the sink through the march.

## 6. Scorecard

| | result |
|---|---|
| **(a)** rate mesh spread ≤ 10 % | **PASS — 2.57 %** (86.2 % after P0, 7.1 % pre-P0) |
| **(c′)** scheme's own balance ±5 % | **PASS 3/3** (−0.15 / +1.37 / −1.09 %) |
| **(c″)** ledger closes with the sink < 0.5 % | **PASS 3/3** (6.0e-4 / 6.1e-4 / 3.7e-4) |
| **(b)** power identity | **covered for kernel↔ledger by (c″)**; march path only via (f) |
| **(c)** identity, keys absent | **PASS 2639/2639, 0 differing** on all four binaries built |
| **(d)** no criterion regression | **SPLIT** — Kant `h150` 0.169 % PASS, `h500` **2.089 % FAIL**; `hu_spall_onset` and `sp_kant_pressure_sweep` both **PASS unchanged** on the final binary |
| **(e)** decomposition | **FAIL**, isolated to the flake/remap path |
| **(f)** Meier `FH` rerun | **reported, §5b** — ROP +16.9 % (down from P0's +33.4 %), window 86.9 → 15.8 K, `flat` 12.8 → 0.0 %. ⚠ side faces carry 62 % of `P_robin` on the cell average |
| **(x)** no cascade | **PASS 18/18** |

## 7. Named debt

1. ⚠⚠ **(e) must be resolved before acceptance.** This repo has shipped
   decomposition bugs twice. The remap writes only from reduced per-column
   buffers, so the candidate is upstream of the write.
2. ⚠ **Smoothing (5.1) costs accuracy on fast transients.** τ = 1 s against a
   0.3–0.4 s onset at `h_conv` = 500 is why (d) splits. The packet's own words:
   *"if 5.1/5.2 change the answer materially, the remap is not conservative —
   that is a witness, not a tuning knob."* Here it is neither: the remap is
   conservative (P3, c′, c″ all hold) and the smoothing is simply mis-scaled for
   a sub-second transient. Smoothed vs raw `ell_eff` on the scored legs:
   **1.237 vs 1.257 / 1.129 vs 1.180 / 1.043 vs 1.143 mm** — a 1.6 / 4.3 / 9.6 %
   gap that grows as firings get more frequent relative to τ. A successor should
   scale τ to the firing interval, or drop 5.1 now that the remap makes it
   unnecessary.
3. ⚠ **`c_f` remains the one empirical constant.** Under `flake`,
   `Removal.H` bypasses the depth scan's `h_col`, so the thickness is **imposed,
   not derived from fracture**. Its lever on the level is weak (P2: 1.17 → 1.30
   across 0.3–0.5), which is the real argument that it is not a tuning knob —
   but it *does* move the residual mesh spread (0.91 % at 0.3, 2.57 % at 0.4,
   3.71 % at 0.5), so it is not inert either. The successor is P1's `K_I(z)`
   scan on the reconstructed exponential, which makes `c_f` an **output**.
4. ⚠ **Side faces keep the cell average** (deliverable 1.2). The 1-D column has
   none, so `surf_side_share` reads 0.000 there and the witness is vacuous until
   the Meier rerun.
5. ⚠⚠ **Side faces carry 62 % of `P_robin` on Meier and still use the cell
   average.** Deliverable 1.2 accepted this because a lateral face has no column
   profile — but at 62 % it is the dominant remaining term on the configuration
   the programme is aimed at. A successor needs a surface temperature for
   lateral faces before any Meier rate is quotable.
6. ⚠ **The flake sink is not wired through the jet march**, so the ledger closes
   to 6e-4 on the 1-D column and −0.16 on Meier.
7. **Feet mode is still owed to P4**, unchanged.
