# D2z-0 PRE-FLIGHT — step zero, before any scored leg

Everything below was checked before `CRITERION.md` was hashed
(`ce81f371962f6588…`) and before the first scored leg was launched.

## Item 0 (mandatory) — name the call site, show each leg reaches it

`full_cell` acts at exactly one place, the per-column winner loop in
`src/Integrator/MMWSpalling/Removal.H`, immediately before
`h_chosen_col[idx] = h_c`:

```cpp
if (inc_full && reg_choice == 1 && !sp_weibull_presc_T_face)
    h_c = dz;
```

Three conditions, demonstrated from each leg's **literal key set**:

| condition | how each scored leg satisfies it |
|---|---|
| `inc_full` | `spall.removal_increment=full_cell` is on `F2`/`F1`'s argv (`run.py --list`), absent on `H2`'s |
| `reg_choice == 1` | `spall.enabled = 1`, `vapor.enabled = 0` in `d2j0_hotdisk/input_hotdisk:83` ⇒ the spall regime is the only one that can win |
| `!sp_weibull_presc_T_face` | `surface_patch.mode = convective_flame` (`input_hotdisk:127`), so the flag is false by its definition at `Removal.H:410-412` |

and the pair that makes the artefact exist at all —
`surface_patch.mode = convective_flame` (`:127`) with
`spall.stress_source = local_thermoelastic` (`:70`) — is what collapses the
depth scan's surface lerp to `T_face = T_top`.

`spall.removal_events_csv = 1` (`:88`) supplies the event log criterion (c)
reads. **No leg fails a row.**

## The parse sits BEFORE the early return — a defect found and fixed here

The key was first added next to `spall.per_cell_removal`, which is **after**
`if (!value.spall_enabled) return;` (`MMWSpalling.H:4145`). With spall off,
`full_cell` was then never read, the `removal_enabled()` cross-check could never
fire, and AMReX only **warns** about an unused ParmParse variable — so such a
leg would run silently in `half_cell` and look like a null result. Caught by the
regression's abort variant, not by inspection. The parse now sits before the
early return, next to the CSV keys that are there for the same reason.

## Reference values — reproduced independently before trusting them

Reproduced from the frozen event logs only (read-only; nothing re-run):

| quantity | `ACTIVE_STEP.md` | reproduced here | verdict |
|---|---|---|---|
| `a_f` range / median (Meier `A`) | 2.20–4.49 µm, med 3.93 | 2.20–4.49, med **3.93** | **exact** |
| events with `a_f == 0` | zero of 91,473 | **0 of 91,473** | **exact** |
| `firings_per_void` (Meier `A`) | 1.977 | **1.977** | **exact** |
| gate gaps, r < 25 mm | −0.2 / +3.0 / −1.8 % | **−0.2 / +3.0 / −1.8 %** | **exact** |
| gate gaps, r < 30 mm | −6.4 / −8.0 / −14.0 % | **−6.4 / −8.0 / −14.0 %** | **exact** |
| cycle period | 4.5 s (2 mm) / 2.2 s (1 mm) | **4.48 / 2.24 s** | **exact** |
| one-dt saving | +0.36 % at both meshes | **+0.357 % / +0.357 %** | **exact** |
| two half-firings one dt apart | 100 % of 45,016 cases | **45,177 gaps, and the gap histogram is strictly bimodal — exactly 1 dt or ≥ 6 dt, with 2/3/4/5 dt all empty** | confirms |

### Disagreements — reported, not substituted

1. **`firings_per_void` baseline.** The packet's **1.977** is **Meier `A`'s**.
   This packet runs the **arbiter**, where the same estimator reads **1.957**
   (`A_f02`) and **1.962** (`A_f02_1mm`). Criterion (c)'s destination (exactly
   1.000) is unaffected; the origin is restated in `CRITERION.md`.
2. **Median residual at `Sp − 1 → 0`.** Packet: **+4.1 µm** at all three.
   Reproduced: **+11.8 µm** (arbiter 2 mm), **+8.2** (arbiter 1 mm),
   **+11.7** (Meier `A`) — 2–3× larger. Narrative only; no criterion uses it.
3. **`h_applied` p5 at 2 mm.** Packet: observed **0.992 mm**. Reproduced, the
   **minimum** is 0.99219 mm (so the packet's figure is the min, not the p5);
   the actual p5 is **1.0078 mm**, i.e. above `dz/2`, not below it.
4. **`Sp − 1` median.** Packet: 0.003–0.006. Reproduced: 0.0048 (Meier `A`,
   in range) but **0.0107 / 0.0077** on the arbiter at 2 / 1 mm.
5. **Command paths.** Four of the five regressions in the packet's command
   block are written `tests/MMWSpalling/<name>/test`; they live under
   `tests/MMWSpalling/unit/<name>/test`. And the tree gate is
   `witness/guarded_post_all.sh <label>`, not
   `studies/tree_sweep_0922/run.py --check`, which does not exist.

## An observation the packet's mechanism story does not cover

The packet scopes `full_cell` away from `prescribed_T` because "16d's fix
already applies and the flat half-cell does not arise". In the new regression's
`prescribed_T` control the surface is genuinely hot (≈1100 K) against a 655 K
top cell, so the lerp is a real gradient — and the scan **still** returns
`h_scan = 0.98437 mm ≈ dz/2 − a_f`, the same landing point as the computed-T
case.

That suggests the deeper cause is not the flat lerp but the **gate location**:
`Sp_top >= 1` is evaluated at the top cell's **centre** (`z = dz/2`) while the
removal depth is measured from the **face**, so a marginal firing lands at
`≈ dz/2` whichever mode supplies `T_face`. This is reported, not acted on —
the packet's guardrail forbids touching the scan, and `prescribed_T` is
excluded as specified. It belongs to whoever revisits Step 16d.

## Budget arithmetic behind criterion (a)

Mean `h_applied / dz` per firing on the frozen legs: **0.5085** (arbiter 2 mm),
**0.5085** (1 mm), **0.5049** (Meier `A`) — i.e. each firing removes slightly
**more** than half a cell, and `Σh / (cells × dz)` is 0.995–0.998. Under
`full_cell` a firing removes exactly one cell. The cadence is set by reheating
(the ≥ 6 dt gap), so the only change is the one dt saved per cell, and
`dt / T_cycle` is **0.016/4.48 = 0.357 %** at 2 mm and **0.008/2.24 = 0.357 %**
at 1 mm — the same number at both meshes, which is why the packet's band is
mesh-independent.

## Smoke (run, then deleted before hashing)

`H2` and `F2` to `stop_time = 25 s` at 2 mm:

* `H2`'s first `removal_events.csv` row is **byte-identical** to frozen
  `A_f02`'s — the lifted argv is right, including the quoted `h_expr`.
* `H2` registers **neither** new thermo column; `F2` registers **both**.
* `F2`: 655 regime-1 firings, 655 voids, `firings_per_void` **1.0000**,
  `removal_increment_cells` **1.0000**, `n_voided` ∈ {1} on every row.
* `H2`: 1137 firings / 654 voids = 1.7385 (the early transient; the full-run
  figure is 1.957).

Witness non-zero on the test leg and absent on the control, as required.

## Gates already passed

* **Tree identity 2639/2639** with the key unset
  (`witness/d2z0.md5` vs `witness/tree0922.md5`), run after the last source
  edit, under `guarded_post_all.sh` — the guard confirms the binary did not
  move mid-sweep.
* **New regression `unit/removal_increment`: 14/14 PASS.**
* `unit/spall_event`, `unit/side_face_flux`, `unit/convective_patch`,
  `unit/jet_enthalpy`, `unit/side_face_area`: **PASS**.
