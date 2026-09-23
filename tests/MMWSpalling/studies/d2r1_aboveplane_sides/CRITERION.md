# D2r-1 — criteria and pre-registered predictions
# Hashed BEFORE either 1 mm leg finished. NEVER edit.

## 0. Disclosure of what was already known when this was written

Honesty about the pre-registration rule matters more than the appearance of it.
When this file was written and hashed, **the three 2 mm legs had already run**
(`B_2mm`, `T_2mm`, `S_2mm`) because criterion (c) requires the plumbing check
*before* spending ~28 min on the 1 mm legs. **Neither 1 mm leg had finished**,
and **criterion (a) — the only scored criterion — depends entirely on them.**

Known at write time:

| block | `B_2mm` | `T_2mm` | `S_2mm` | T/B − 1 | S/B − 1 |
|---|---|---|---|---|---|
| 50–100 | 1.5985 | 1.6001 | 1.5986 | **+0.1 %** | +0.0 % |
| 100–150 | 1.5563 | 1.5571 | 1.5560 | **+0.1 %** | −0.0 % |
| 150–200 | 1.6178 | 1.6185 | 1.6179 | **+0.0 %** | +0.0 % |

Not known: `B_1mm`, `T_1mm`, and therefore every gap.

## 1. Goal item 1 — the baseline-validity proof: **PASSED**

`B_2mm` reproduces D2l's frozen `A_f02` **byte for byte** —
`_removal_events.csv`, `thermo.dat`, and every compared `Level_0/Cell_D_00000`.
So D2q-2b/2c's changes to the side-face path did **not** perturb it, and the
frozen **6 / 8 / 14 %** remains a valid baseline. The packet's stop condition
does not fire.

## 2. Pre-flight item 0 — call sites, and proof each leg reaches them

| changed behaviour | call site | reached because |
|---|---|---|
| per-region side factor, kernel | `MMWSpalling.H` side-face branch, inside `if (col_in_patch)` → `if (nfx + nfy > 0)` | every leg sets `side_face_flux = 1`; arbiter is `convective_flame` (not prescribed); scoring disk r < 30 mm is inside `surface_patch.radius = 0.030`; spall removal on ⇒ void cells exist |
| per-region side factor, **ledger twin** | the mirrored side-face block | same predicate; **verified by `ledger_err` ≤ 5e-15 on both new legs** |
| `side_P_face_above`, `side_cols_above` | ledger twin, gated `side_face_factor_above >= 0` | `T`/`S` set it; `B` does not (deliberate — see §3) |
| `WallStepDiagnostic` | `BuildFlameColumns` tail, same gate | `T`/`S` set it |

**Demonstrated empirically, not asserted:** `side_P_face_above` = **43.57 W** on
`T_2mm` and **exactly 0.00 W** on `S_2mm`, with ~103 above-plane faces present in
both. The branch is reached, and the separator zeroes exactly the intended term.

## 3. Why `B` carries no witness column, and what pays for it

The new thermo variables are registered **only when `side_face_factor_above` is
set**. A new column changes `thermo.dat`, which would have broken the byte
comparison in §1 — the packet's own stop condition. So `B` cannot carry
`side_P_face_above`, and criterion (c)'s "= 0 on `B`" is discharged by **`S`**
instead, which carries the column and reads exactly 0. This is a deviation from
the packet's wording, forced by item 1, and it is recorded rather than hidden.

## 4. Criteria

- **(a) THE SCORED ONE.** `T_2mm` vs `T_1mm`, disk-mean rate over r < 30 mm,
  blocks 50–100 / 100–150 / 150–200, against the baseline's **6 / 8 / 14 %**.
  **Pass: all three ≤ 5 %.**
- **(b) THE SEPARATOR.** `S_2mm` must reproduce `B_2mm` within noise.
  **Already satisfied: +0.0 / −0.0 / +0.0 %.** If `S` had closed the gap the
  decisive term would be the tops and D2j-0b's `C2` would be contradicted.
- **(c) THE WITNESS.** `side_P_face_above` > 0 on `T`, = 0 on `S`.
  **Already satisfied (43.57 / 0.00 W).**
- **(d) THE DRAIN.** `wall_step_q`, `T` vs `B` at the same mesh, against
  D2o-0b's 302 kW/m² / 58 % of q_pin.
- **(e) THE CEILING.** Absorbed power and disk-mean rate, `T` vs `B` at the same
  mesh (one key apart). **Pre-registered bound: `T` must not exceed `B` by more
  than 15 % in any block.** At 2 mm the measured excess is +0.1 %.

## 5. Predictions

- **P1 — I predict the test does NOT close the gap: the null branch P4.**
  The 2 mm data already in hand shows the mechanism moves the rate by **+0.1 %**
  while delivering only **43.6 W** to above-plane side faces. For the gaps to
  close, `T_1mm` would have to exceed `B_1mm` by ~6 / 8 / 14 %. **I predict
  `T_1mm`/`B_1mm` − 1 lands in +0.1 … +1.5 %**, leaving gaps of roughly
  **5–6 / 7–8 / 12–14 %** — improved by at most ~1–2 points, not cleared.
  **Pass on (a) is therefore predicted to FAIL.**
- **P2 — the separator reproduces the baseline.** Already confirmed to ±0.05 %.
  Consistent with D2j-0b's `C2`.
- **P3 — the ceiling holds easily.** Predicted excess < 2 % at both meshes, far
  inside the 15 % bound. Above-plane side heating adds little absorbed power
  because the faces are few (~103) and `h` = 142 is modest.
- **P4 — if P1 is right, the reading is:** above-plane **side** heating is not
  the decisive term either. Combined with D2j-0b `C2` (tops inert) and `C0` /
  D2q-2c no-plane (everything converges), that would mean the convergence in the
  `C0` cases comes from something **other** than above-plane heating per se —
  most likely the far larger **below-plane** coefficient those runs apply to the
  same faces (700 vs 142, and `f` = 1 vs 0.2), i.e. a magnitude effect, not a
  region effect. **Do NOT raise `h` above 142 to rescue it.**
- **What would refute P1:** any block where `T_1mm`/`B_1mm` − 1 exceeds +4 %.

## 6. Guardrails restated

`side_face_factor = 0.2` below the plane is a frozen stand-in, never scanned.
142 is the annulus derivation's value, never varied. `side_face_h_max` is never
set. `per_cell_removal` stays off. **No absolute number here is physics — only
`T − B` and `T − S`.** The arbiter is a flat plate and carries **no shape
result**.
