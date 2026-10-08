# D2v pre-flight — written BEFORE any campaign leg launched

Three obligations: **item 0** (the standing rule, six packets bought it), the
**reference reproductions**, and the **design readings** I had to make where the
packet was ambiguous. The compiled-map trace is a **separately hashed**
pre-flight: `HMAP_TRACE.md` / `.sha256` (`90d84f2f8b2172c7…`), **PASSED T1–T10**.

---

## 0. Reachability ledger — every leg, every call site

| # | behaviour | call site | reached because | legs |
|---|---|---|---|---|
| 1 | confined `h` below the plane | `:1778` true branch of `(s > 0.0 \|\| jet_above_h)` | every leg sets `h_expr`; `jet_closure = enthalpy` so this is the branch taken | all |
| 2 | confined `h` **above** the plane | `:1778`, via `jet_above_h` | `M2`/`TP` set it; `M0`/`M1`/`TF` do not | `M2`, `TP` |
| 3 | above-plane `T_gas` | `:1792-1793`, gated `jet_T_above >= 0` | `M2`/`TP`; parser aborts if set without `jet_above_h` | `M2`, `TP` |
| 4 | the `hj >= 0` abort | `:1781-1786` | evaluated for every in-patch live column; floor is `1e-3 > 0` | all |
| 5 | side-face term | `:1655-1656` inside `col_in_patch` → `nfx+nfy>0` | `M1`/`M2` set `side_face_flux=1`; `radius = 0.2` covers the whole domain | `M1`, `M2` |
| 6 | per-region side factor | `:1653-1654` on `fs[col_c] > 0.0` | `M1` sets 0.0, `M2` sets 1.0 | `M1`, `M2` |
| 7 | ledger twin | `:2754-2885` | `energy_ledger.enabled = 1`; witnessed by `ledger_err` | all |

### Empirical demonstration — not a grep

**The five new guards were each made to fire**, on the real Meier key set:

| set | result |
|---|---|
| `jet_above_h=1` alone | **ABORT** `:6280` "requires surface_patch.jet_T_above" |
| `jet_T_above=1690` alone | **ABORT** `:6284` "requires surface_patch.jet_above_h = 1" |
| `jet_above_h=1`, `jet_T_above=200` (≤ `jet_T_ent`) | **ABORT** `:6288` "must be finite and > jet_T_ent" |
| both, correctly paired | **runs, rc = 0** |
| either, under `jet_closure = none` | **ABORT** `:6313` "requires jet_closure = enthalpy" |

**And the map itself was measured, not asserted** (item 3): the floor branch
inverts to **1407.5 W/m²K against a known 1449** (+2.86 %), the wall branch to
**142.00 flat** where a blend would read 559–695, and `max dT` beyond
`r = 76 mm` is **0.013 K** over 2463 columns.

**No leg fails to reach any site, so the rule does not stop this packet.**

### ⚠ A consequence of the change that the packet does not state

`JetEnthalpyMarch` **skips above-plane columns** at both its column pass
(`:2126`) and its bin pass (`:2212`), testing `!(flame_s[c] > 0.0)`. Two
consequences, both verified:

1. **`jet_T_above` survives the march** — nothing overwrites `flame_Tg` for
   those columns. This is why the key works at all.
2. **The heat delivered above the plane is NOT deducted from the jet's enthalpy
   budget.** The march's invariants (`P_face <= P_cap`, `P_face + P_exhaust ==
   P_cap`) cover below-plane columns only. So on `M2` the above-plane source is
   **unbounded by the jet**: it cannot be capped by the gas's own enthalpy.
   **That is exactly the mechanism criterion (d) exists to catch**, and it is
   why (d) is scored on every leg. The global energy ledger still balances
   (the flux is counted in `P_robin`), so `ledger_err` will not see it.

---

## 1. Reference reproductions

**Reproduced exactly:**

* The stale `h_expr` evaluated at `s = 50 mm`: **1449 / 1449 / 1092 / 912 / 686
  / 549 / 361 / 229 / 164** at r = 5 / 18 / 25 / 30 / 40 / 50 / 76 / 120 / 168 mm
  (packet: 362 at 76 mm, rounding). It decays as ~1/r and never floors ✓.
* The frozen key set: `prob_hi = 0.12 0.12 0.4`, `n_cell = 60 60 200` ⇒ dz = 2 mm,
  `timestep = 0.016`, `radius = 0.2`, `foot_r_outer = 0.040`,
  `jet_closure = enthalpy`, `jet_T_ent = 293.15`, `jet_T_ent_mode = exhaust`,
  `nozzle_descent = feet`, `pinned_idle_cycles = 1.0e9`,
  `foot_body_clearance = 0`, `jet_free_surface = 1` ✓. `side_face` appears **zero
  times** ✓.
* The Meier csv: 96 / 92 / 88 / 87 / 87 / 86 / 85 at 10–200 mm ✓, with its own
  header stating the visible width is a **lower bound** at **±5 mm**.
* `jet_T_exhaust` = **1689.5 K** at 572 s (packet: 1690), `jet_T_rec` = 1689.5
  (packet: 1689) ✓ — uncooled.

### ⚠ Three disagreements, reported not replaced

1. **The model Ø(z) table.** With my own interpolated metric I measure
   **144.3 / 129.1 / 114.3 / 105.3 / 99.5 / 92.6 / 86.6 mm** at
   10 / 25 / 50 / 75 / 100 / 150 / 200 mm, against the packet's
   **141 / 128 / 113 / 105 / 99 / 92 / 86** — a **+0.3 to +3.3 mm** offset.
   I could not reproduce the packet's values with **five** conventions (mean /
   median / min / max per bin, bin centres vs lower edges, corner-referenced
   radius, coverage-50 % crossing). The face reconstruction itself is
   unambiguous: cumulative `n_voided` and `last k_top − n_voided` agree on
   **all 3600 columns**.
   **Consequence:** criterion (a)'s `M0` reads **37.2 mm**, not the packet's
   35.3. This changes no decision — the bar is ≤ 20 mm against an M0 of 35–37,
   and **the same metric is applied to all three legs**, so `M1 − M0` and
   `M2 − M1` (the packet's actual purpose) are untouched.
2. **(d)'s baseline.** Beyond r = 76 mm the frozen run gives **mean 0.216 mm,
   max 6.00 mm** (three single-cell events in one column), not the packet's
   0.08 / 4.00. **Exactly 0 beyond 90 mm** ✓ as stated. Both are far below the
   2 mm bar, so the bar is safe either way.
3. **The 50–76 mm ring** strips **32.5 mm**, not 33.6.

Minor: the packet's `w1 = 0.625 at r = 46 mm` is a **linear** ramp value; the
raised cosine this packet requires gives **0.691**, so the rejected blended form
would have been ~461 W/m²K there rather than ~430 — the packet's argument
against it is **stronger**, not weaker, than stated.

---

## 2. Design readings recorded before launch

**`input_drilling`, not `input_feet_d2c`.** The packet's configuration table
says to use `input_feet_d2c`, but its own item 5 says `M0` **must** come from
the frozen run's exact keys. The frozen `d2i_gate/LI0_2mm` run uses
**`input_drilling`** with all 60 keys supplied on the command line. Item 5 wins;
the harness imports `d2i_gate/run.py` by path and takes its argv verbatim.

**`jet_T_above = 1690.0 K`, declared and never fitted.** The frozen run's own
uncooled exhaust is 1689.5 K (1676.5–1689.7 over t ≥ 500 s). Meier's water
jacket is not in the model, and **its "186 W/K" is the WATER stream's `mcp`**
(160 L/h × 4186), not a conductance — the water ΔT is **unreported**, so the gas
drop is **19 K** for ΔT_water = 5 K and **194 K** for 50 K, against the run's own
marched `mcp` of 11.99 W/K per quarter sector (**47.96 W/K** full). Declaring a
drop would mean **inventing ΔT_water**, which the guardrails forbid. So the
measured uncooled value is used, and it is:
* an **upper** bound on the physical gas temperature (no jacket), and
* a **lower** bound relative to `M1`/`M2`'s own exhaust, which will run hotter
  because the confined map extracts less power from the block top.
Both directions are stated. If (d) fails at 1690 K **that is the finding**, not
a reason to lower it.

**`min`/`max` in the blends.** The packet says "No `min`/`max` in the blends"
but also requires `w1 = 1` exactly inside `r_imp` and `0` exactly by `r_ann_r` —
which cannot be done without clamping. I use the **numeric-first** clamped
raised cosine, the form already proven against this parser by D2t's compiled
trace and used by the frozen `h_expr` itself (`min(0.09,max(0.015,s))`,
`max(0.01875,r)`). No denominator is a `max()`. **Item 3's trace is the
empirical check and it passed**, so the tension is resolved by measurement.

**`h_out = 1.0e-3`.** Under `jet_closure = enthalpy` the abort at `:1781` admits
`hj >= 0`, so zero would in fact be legal here — the `> 0` rule the packet cites
belongs to the **non-jet** branch at `:1787`. `1e-3` is kept anyway: it is
C1's value and it costs nothing (measured `max dT` beyond the mouth: 0.013 K).

**The metric.** One definition, fixed here and applied identically to every leg:
final `k_top` from `_removal_events.csv`, **fixed physical 2 mm** radial bins,
**mean** face height per bin, and the crossing radius **linearly interpolated**
between bracketing bins. D2t's `D = 2·r_centre` snapped to a 4 mm lattice —
4 mm is 4.7 % of an 85 mm hole and equals the reference's own ±5 mm.
