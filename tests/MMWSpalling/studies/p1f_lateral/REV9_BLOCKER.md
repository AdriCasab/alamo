# P1f revision 9 — deliverables 9/10/11 LAND; deliverable 12 and both order-1 gates are BLOCKED

## Step zero: the packet's ledger is sound

Every row resolves, and the reference values reproduce with my own estimator:

| row | packet | verified |
|---|---|---|
| the bug site | `Removal.H:624`, `K_Ic` MPa vs `K_I` Pa, no conversion | ✅ confirmed verbatim |
| vertical converts | `:1333`, `:1357` | ✅ both via `MPa_sqrtm_to_Pa_sqrtm` |
| **only site in `src/`** | `:624` | ✅ swept `KIcNasseri\|KIcLinear\|MPa_sqrtm_to_Pa_sqrtm`: 3 consumer groups, only `:624` lacks it |
| both lateral branches MPa | table **and** `kic.K_Ic0` | ✅ `:622` fallback and `:624` table |
| `C_th` | `:488`, Pa K⁻¹ | ✅ `phase_E_loc` in Pa |
| `sp_driving_diff` | `MMWSpalling.H:5632`, default `absolute`, read `Removal.H:1416` | ✅ |
| **(b2) thresholds** | absolute+8 µm **756 K** any `ell`; differential+8 µm **942/998/1085** at 3/5/8 mm | ✅ **756.4**; **942.0 / 998.4 / 1084.9** — to the kelvin |

⚠ My first threshold pass gave 689 K because I built the Nasseri table from the packet's prose
("1.43 … 0.22") instead of reading it. The real anchors are **{1.43, 1.35, 0.98, 0.43, 0.22}**
(`SpCriterion.H:102`). A packet that quotes a table's endpoints is not quoting the table.

## Landed and verified

* **9 — the unit fix.** `K_Ic_MPa` on both branches, then one
  `MPa_sqrtm_to_Pa_sqrtm`. Build clean.
* **10 — the four diagnostics**, `K_I_lat` / `K_Ic_lat` / `q_lat` / `sigma_lat` (Pa·√m, Pa·√m,
  dimensionless, Pa at ξ = 0), appended **inside the existing `lat_csv` guard**.
  ⚠ **Identity-safe, no re-pin owed**: with the key off the header is character-for-character
  what it was, and 5 pinned `*_removal_events.csv` reproduce **byte-for-byte** against `p1f.md5`
  on the rebuilt binary. The gate hashes 47 of these files, so this mattered.
* **11 — the form follows the run** via `lat_diff_loc = sp_driving_diff`; `absolute` drops the
  `(1 − e^{−h_f/ell})` factor. No abort, no key on any control's argv.
* Targeted cover: `robin_face` **PASS**, `robin_pinned` **PASS**.

## ⚠⚠ THE BLOCKER — a conjunction, exactly as in revision 3

**(w)** is a **gate** requiring `0 < lat_fires < lat_cand` on `RISER`. With the units corrected,
`RISER` fires **0 of 80 candidates**. Measured, not inferred:

| | |
|---|---|
| firing threshold, `a0 = 20 µm`, `absolute` (the fixture's own config) | **657.3 K** |
| hottest open lateral face the fixture reaches at 60 s (80 faces) | **491.4 K** |
| p50 / p95 | 480.8 / 491.4 K |
| **shortfall** | **165.9 K** |

The fixture cannot reach its own criterion. `surface_patch.side_face_T = 1200` and
`side_face_h = 169.3`; the rock's own conductance at 60 s, inferred from the measured face
temperature, is **G ≈ 605 W m⁻² K⁻¹**, so the face sits at only `T0 + 0.219·(T_gas − T0)`.

**Three consequences, all order-1:**

1. **(w) is unreachable**, not failed — there is no firing to bound.
2. ⚠⚠ **(v) becomes VACUOUS.** It was specified as "the negative case the suite never had", but
   post-fix **every** leg is a cold-face leg, so a cold leg scoring 0 firings proves nothing. This
   is the non-vacuity trap: the probe must measure something the mechanism moves.
3. **(b2) is unobservable.** Distinguishing 756 K from 942 K requires faces that reach them;
   these reach 491 K, so both forms give 0 firings and the legs are trivially equal.

**What the guardrails forbid, and why I stopped rather than proceeded.** The levers are
`surface_patch.side_face_h` / `side_face_T` — and the packet says *"Do not touch … any
`surface_patch.*` key beyond reading `side_face_staircase`."* Lowering `a0` is also forbidden
("reuse the shipped flaw"), and would be tuning the calibration to make a test fire.

## The arithmetic for whoever authorises the fix

To put the hottest face at **700 K** (comfortably past 657.3), holding `T_gas = 1200`:

* `side_face_h ≈ 500` (from `h/(h+G) = 406.85/906.85` at `G = 605`) — a **×2.95** increase, and
  the cheapest lever.
* or `side_face_T ≈ 2150` at the shipped `h` — synthetic, and further from anything physical.
* or a **×16 longer soak**: `G ∝ t^{−1/2}`, so 240 s reaches only **619 K** and ~960 s is needed,
  which costs compute and risks the containment the (c2) sweep already has to filter for.

⚠ **Cleanest route, touching no shipped input:** add the (w)/(v)/(b2) legs with a **per-leg argv
override** (`run_sweep` already overrides `amr.n_cell`, `geometry.prob_hi`, `stop_time`,
`spall.seed_void` this way), leaving `unit/lateral_drive/input` and every production input alone.
That changes the **drive** for one synthetic leg, never the flaw, toughness, `c_f` or `flake_h_max`.
I have not done it: it needs the guardrail lifted in one line, and this is a gate.

## Why Leg B was NOT launched

(w) and (b2) are **order 1**; Leg B is order 2. The run-order rule stops long legs when an order-1
gate does not pass, and these two are the only cover against a second over-firing bug — the exact
cover whose absence let 224,359 spurious events through. ⚠ Leg B itself is **not** blocked on
physics: the packet's P1 expects sparse firing (441 faces above 756 K), so the Meier geometry does
reach the threshold. It is blocked on having no two-sided assertion in front of it.
