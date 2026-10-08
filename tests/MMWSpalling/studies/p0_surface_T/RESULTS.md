# P0 — a real surface temperature: THE CRITERION IS FIXED AND IT EXPOSED A BIGGER ONE

**The firing criterion is now mesh-independent to 0.16–0.35 %, measured in two
independent configurations.** The pre-P0 criterion fired the same rock at a
true surface temperature that moved **14 %** across 2 / 1 / 0.5 mm in the Kant
anchor and **38 %** in the 1-D column; under `spall.surface_T = reconstructed`
it fires at the same surface temperature on every mesh.

**And that is not the headline.** ⚠⚠ **The Robin boundary condition reads the
CELL AVERAGE, not the surface**, so with the criterion fixed the two no longer
cancel: the 1-D steady rate's mesh spread goes from **7 % to 86 %**. P0 on its
own makes the rate *worse*. It must be paired with a surface-temperature
boundary condition before any rate is quoted. Details in §4 — **this is the
single most important thing in this packet for the next planner.**

## 0. Provenance

* Binary **`f65d10eb…`** (pre-P0 `ea833f2f…`). `BINARY.sha256`, `GATE.txt`.
* Tree gate vs `witness/d2z1.md5`: **2639/2639, zero differing, zero added,
  zero removed**, on both binaries built during this packet, both with
  `--force`. The first gate was **discarded**, not reused, when the
  `diagnostic` mode was added.
* Decomposition witness: `removal_events.csv` **bit-identical** across four
  (`max_grid_size`, ranks) layouts — see `GATE.txt`.
* Nothing committed.

## 1. ⚠ THE PACKET NAMED ONE CRITERION SITE. THERE ARE TWO.

`ACTIVE_STEP.md` locates the defect at `Removal.H:846-858` — correctly, but
that is the **depth scan**, which only runs once a column is already spall-ready
and only governs the flake thickness. The **onset** criterion, and the quantity
the Kant anchor is scored on, is `Sp_field`, computed in
`UpdateSpAfterMechanics` (`MMWSpalling.H:5186`). That site has **no surface
lerp at all**: with `a_f` = 4–20 µm against `dz` = 0.5–2 mm every `K_I` sample
lands in the top cell, so σ_xx is flat at the cell average and `K_Ic` is
evaluated at the same value.

Implementing the packet as literally written would have re-anchored `a0_ig` on
a criterion the Meier run does not use — *"anchoring on a different formula
from the one in use, the exact disease this packet exists to cure"*, in the
packet's own words about revision 1. **Both sites are changed**, and they share
one function (`ReconstructSurfaceT`) so they cannot drift.

⚠ The two sites needed different plumbing for the same reason: neither may read
`k-1` directly. The pre-P0 code deliberately clamps its sample index into the
validbox, so a direct read would make the criterion depend on
`amr.max_grid_size` at every box seam. `Removal.H` already had a gathered
column profile (PASS 1b); `UpdateSpAfterMechanics` did not, and one was added.

## 2. The anchor — `validation/kant/sp_kant_onset_meiercfg/`

24 mm Kant cube, lateral rollers, Meier material / flaw / stress settings, at
2 / 1 / 0.5 mm and at `h_conv` = 150 and 500 (Kant's `h_fl,min`), **one
timestep (0.0125 s) on all three meshes** so `dz` is the only variable.

| leg | t [s] | `T_top` | `T_surf` | **ΔT_fire** | ℓ_eff [mm] | clamps |
|---|---|---|---|---|---|---|
| R150 2.0 | 19.000 | 649.85 | 756.46 | **463.31** | 3.650 | none |
| R150 1.0 | 21.775 | 706.33 | 755.88 | **462.73** | 4.330 | none |
| R150 0.5 | 23.100 | 731.89 | 755.72 | **462.57** | 4.685 | none |
| R500 2.0 | 1.437 | 482.05 | 758.95 | **465.80** | 0.913 | none |
| R500 1.0 | 1.450 | 575.26 | 758.79 | **465.64** | 0.907 | none |
| R500 0.5 | 1.812 | 671.68 | 758.40 | **465.25** | 1.169 | none |

**(a) PASS** — spread **0.161 %** at h 150 and **0.116 %** at h 500, bar 3 %.
**(b) PASS** — all six inside Kant's 390–560 K band.

⚠ The cell average `T_top` moves **12.6 %** across the same three meshes
(649.85 → 731.89). That is the defect, measured.

### 2.1 ⚠ THE OBVIOUS WITNESS IS VACUOUS, AND THE REAL ONE IS NOT

A threshold criterion fires when **its own** temperature variable reaches the
threshold. So the pre-P0 criterion is *also* flat — it fires at a cell-average
rise of 461.70 / 461.73 / 461.67 K, spread **0.01 %**. Scoring (a) alone would
have passed vacuously, and the first version of this harness did exactly that.

The quantity that actually differs is the **physical state at firing**. A third
mode, **`spall.surface_T = diagnostic`**, computes and logs the reconstruction
while leaving the criterion bit-for-bit the pre-P0 criterion, so the control leg
reports its own true surface temperature:

| | 2.0 mm | 1.0 mm | 0.5 mm | spread |
|---|---|---|---|---|
| control fires at cell-average ΔT | 461.70 | 461.73 | 461.67 | 0.01 % |
| **control's TRUE surface ΔT** | **557.62** | **508.75** | **484.92** | **14.06 %** |
| reconstruction's surface ΔT | 463.31 | 462.73 | 462.57 | **0.16 %** |

The control's true surface temperature converges **downward toward the
reconstruction's value** as the mesh refines — which is what "the coarse mesh
was firing hotter rock than it thought" looks like. **(w) PASS.**

### 2.2 Deliverable 7 — the re-anchor

Bisection at 1 mm, deterministic flaw, target ΔT_fire = 461 ± 20 K:

| `a0` | ΔT_fire [K] |
|---|---|
| 6.0 µm | 489.63 |
| **8.0 µm** | **462.73** ✅ |
| 11.0 µm | 432.34 |

**`a0_ig` = 8.0 µm** (`a0_gb` = 40 µm = 5×). The same value lands the Weibull
leg at 463.39 K, so one number anchors both flaw models.

✅ **P1 HELD, and close**: predicted 6–20 µm with a point estimate of ≈8.5 µm;
measured **8.0 µm**. Well inside the 1–100 µm escalation window, so the 2 mm
mesh can represent the layer and there is no AMR escalation.

⚠ **The Weibull legs are REPORTED, not scored**: 455.67 / 463.39 / 472.65 K,
spread **3.66 %**. `V0 = 8e-9` is `V_cell` at 2 mm **only**, and the first
firing is the luckiest of 4× as many surface cells at each refinement. Both
move the onset *flaw* with `dz` by construction, so a Weibull leg cannot
measure whether the *criterion* is mesh-independent. `V0` is held fixed across
the three meshes as the packet directs, and the size effect is quantified here
rather than hidden.

## 3. The 1-D gate — `validation/p0_steady_column/`

8 × 8 × 120 mm laterally uniform column, convective BC at `h_conv` = 600 /
`T_flame` = 1673.15 (**0.55 MW/m² at the firing surface** — the packet asks for
a uniform-flux column and this model has no prescribed-flux BC; revision 1 was
corrected for the same reason).

| leg | v [m/h] | `T_surf` | `T_top` | `T̄_flake` | ℓ [mm] | firings/void |
|---|---|---|---|---|---|---|
| F 2.0 | 7.142 | 761.21 | 471.01 | 589.02 | 0.837 | **2.371** |
| F 1.0 | 4.579 | 760.37 | 559.01 | 588.49 | 0.795 | **1.273** |
| F 0.5 | 2.931 | 759.58 | 659.48 | 659.48 | 0.992 | 1.000 |
| C 2.0 | 7.200 | 760.61 | 470.61 | 588.64 | 0.836 | 1.000 |
| C 1.0 | 4.526 | 759.89 | 558.53 | 588.19 | 0.794 | 1.000 |
| C 0.5 | 2.931 | 759.58 | 659.48 | 659.48 | 0.992 | 1.000 |
| **D 2.0** (pre-P0) | **2.263** | **1103.30** | 755.61 | 805.26 | 1.600 | 1.000 |
| **D 1.0** (pre-P0) | **2.160** | **938.98** | 755.93 | 755.93 | 1.412 | 1.000 |
| **D 0.5** (pre-P0) | **2.109** | **847.01** | 755.67 | 755.67 | 1.344 | 1.000 |

**(a) PASS in a second, independent configuration** — active removal, a moving
front, Meier-like flux: ΔT_fire = 468.1 / 467.2 / 466.4 K, spread **0.35 %**.

⚠ And the `D` column is the packet's opening table, **measured**: the pre-P0
criterion fires the rock at a true surface temperature of **1103 / 939 / 847 K**
— a **38 % spread in the rise**. The packet predicted 1215 K at 2 mm and 1003 K
at 1 mm for a 823 K cell average; at this case's 756 K cell average the same
formula gives 1103 and 939. **The derivation is confirmed by measurement.**

✅ **P6 HELD**: `firings_per_void` rises above 1 (2.371 / 1.273 / 1.000) and
`h/dz` falls below 1 **under `flake` only**; the `cell` path stays at exactly
1.000 with `h_applied == dz` on every row.

## 4. ⚠⚠ (c) FAILS, AND THE REASON IS THE POINT OF THIS PACKET

| mesh | v measured | packet's balance | error |
|---|---|---|---|
| 2.0 mm | 1.984e-3 | 8.512e-4 | **+133.1 %** |
| 1.0 mm | 1.272e-3 | 8.536e-4 | **+49.0 %** |
| 0.5 mm | 8.143e-4 | 6.887e-4 | **+18.2 %** |

The packet's target `v = h(T_g − T_surf,fire) / (ρc(T̄_flake − T0))` makes two
assumptions about the code, **and neither holds**:

**(i) ⚠⚠ THE BC'S SURFACE TEMPERATURE IS NOT THE CRITERION'S.** ⚠ *Corrected
after the first draft of this file, which said flatly that "the BC reads the
cell average". That is true only of `robin_form = cell`, which is the DEFAULT
and is what these new cases were built with — it is **not** what the Meier legs
use. The three forms differ and the correction matters:*

| `surface_patch.robin_form` | what it calls "the surface" |
|---|---|
| **`cell`** (default, and what §3's legs use) | the cell average. `q = h(T_g − T_c)`, no half-cell resistance at all |
| **`face`** | a real face temperature from a Newton solve, but on a **LINEAR** half-cell profile: `(2k/dz)(T_f − T_c) = h(T_g − T_f) − losses` |
| **`pinned`** (the Meier legs, D2a C1) | that face Newton, clamped to `T_pin` — which `Removal.H:1544` records as `rev_T_top`, **the CELL AVERAGE at firing** |

None of the three is the criterion's surface temperature, and P0 widened the
disagreement in all three. At 2 mm under `cell` the criterion says 761 K while
the BC drives on 471 K, so it delivers `600·(1673 − 471)` = 721 kW/m² where a
surface BC would deliver 547 — **32 % too much heat, by a mesh-dependent
factor**. Under `face` the gap is smaller but the same in kind: the half-cell
model is **linear** while the reconstruction is **exponential**, and at these
rates ℓ_eff (0.8–1.0 mm) is below dz/2, which is exactly where a linear
half-cell profile understates the surface most.

**(ii) `flake` CANNOT CHANGE THE ENERGETICS.** The code has no enthalpy
subtraction on removal at all: a voided cell keeps its thermal state and simply
becomes inert, so heat leaves exactly when geometry does. Partial recession
moves `phi` **without changing T**, so the column fires again on the very next
step and the cell still voids carrying its whole enthalpy. Measured: F and C
give the same rate to 0.8 % (7.142 vs 7.200, 4.579 vs 4.526, 2.931 vs 2.931).
⚠ **Deliverable 5's "subtract only ρc·h_f·(T̄_flake − T0)·A_cell" describes a
subtraction that does not exist anywhere in this code.** That is consistent
with D2z-0's measured +0.35/+0.44 % rate effect for the same class of change.

The balance that **does** close is whole-cell removal driven by a cell-average
flux: **+6.3 / +9.9 / +6.6 %** across the three meshes. That identifies (i) and
(ii) as the whole of the discrepancy.

### 4.1 ⚠⚠ THE RATE'S MESH DEPENDENCE GOT WORSE, NOT BETTER

| | 2.0 mm | 1.0 mm | 0.5 mm | **spread** |
|---|---|---|---|---|
| `F` reconstructed, `robin_form = cell` | 7.142 | 4.579 | 2.931 | **86.2 %** |
| `C` reconstructed + whole cell, `cell` | 7.200 | 4.526 | 2.931 | **87.4 %** |
| **reconstructed, `robin_form = face`** | **4.141** | **3.291** | **2.520** | **48.9 %** |
| `D` **pre-P0**, `robin_form = cell` | 2.263 | 2.160 | 2.109 | **7.1 %** |

⚠ **The `face` row is the honest comparison and it was added after the first
draft.** Scoring the 1-D gate on the *default* `cell` form overstated the
effect: the correct finite-volume BC halves the mesh spread, 86 % → **49 %**.
It does not remove it, and 49 % is still seven times the pre-P0 7.1 %, so the
conclusion stands — but the number that was quoted first was the crudest of
the three forms, not the model as it is used on Meier.

Pre-P0 the criterion and the BC used the **same** cell average, so their errors
were coupled and largely cancelled in the rate. P0 fixes one of them. The
column now fires when the *surface* reaches 760 K, which happens while the cell
average is only 471 K at 2 mm and 659 K at 0.5 mm — so the BC's driving ΔT, and
hence the delivered power, becomes strongly mesh-dependent.

> **P0 IS NOT DELIVERABLE AS A RATE IMPROVEMENT ON ITS OWN.** The surface
> temperature must reach the Robin BC before any rate from this branch is
> quotable. This is a *criterion* packet and it does what it claims for the
> criterion; the rate consequence is an escalation, not a regression to hide.

✅ The **direction** is what P2 predicted (the rate goes sharply up: +216 % at
2 mm in this column), and the packet's stop rule — *"a rise below +50 %, or any
fall, means the reconstruction did not reach the criterion"* — is passed by a
wide margin. But the **cause** is not the one P2 gives.

**(x) no cascade** on any leg: `firings_per_void` ≥ 1 everywhere. Revision 1's
withdrawn enthalpy-transfer scheme would read as *fewer* than one firing per
void; the partial-recession signature is the opposite.

## 5. The witnesses the packet asked for

| witness | reading |
|---|---|
| the endpoint edit missed | `T_surf` − `T_top` = **+106.6 / +49.6 / +23.8 K** at 2 / 1 / 0.5 mm — non-zero, and **shrinking with `dz` exactly as it must** |
| second cell void/garbage | `recon_fallback` = **0** on every step of the Meier rerun; guarded by the same contiguous-solid-run test that bounds the scan |
| ℓ_eff clamp always binds | `recon_ell_clamped` = **0.0000** over 3600 columns × 838 thermo steps. **P4 HELD** (predicted < 5 %); the gas ceiling also never bound |
| ℓ_eff is noise, not a profile | `d_ell_rel` = **3.0e-4 … 5.1e-3** per step — smooth, not noise |
| box-local two-cell form | **bit-identical** over four (`max_grid_size`, ranks) layouts, `GATE.txt` |
| partial recession never recedes partially | `h_applied` < `dz` on 176 rows, `firings_per_void` 2.371 |
| the revision-1 cascade | absent — §4, **(x)** |

⚠ `recon_flat` reads ≈100 of 3600 columns (**2.8 %**) on the Meier rerun: the
cold columns outside the heated patch, which have no usable gradient and
correctly degrade to the cell average with a factor of exactly 1.0.

## 6. Scorecard

| | result |
|---|---|
| **(a)** firing `T_surf` ±3 % over 2/1/0.5 mm | **PASS** — 0.161 / 0.116 % (Kant), 0.35 % (column) |
| **(b)** Kant onset in 390–560 K, both `h_conv` | **PASS 6/6** |
| **(c)** 1-D steady state ±5 % | **FAIL** +133 / +49 / +18 % — diagnosed, §4 |
| **(d)** Hu LRST + Kant pressure sweep | **PASS both**, §6.1 |
| **(e)** identity, both keys | **PASS 2639/2639, 0 differing** |
| **(f)** Meier `FH` rerun | reported, §6.2 |
| unit `unit/surface_T_recon` | **PASS 6/6** |

### 6.1 (d) — the regressions the keys must not move

Both run on the P0 binary with **neither key set**, so both exercise the
default path end to end rather than trusting the tree gate alone.

* **`validation/hu/hu_spall_onset` PASS.** Granite onset **37.6 s**, LRST
  **778.84 K** (Hu primary 37.0 s / 773.15 K); sandstone **93.0 s**,
  **913.42 K** (89.0 s / 893.15 K); material ordering preserved.
* **`validation/kant/sp_kant_pressure_sweep` PASS.** FEM vs Kant Eq. 15 at
  p = 0 / 27 / 48 MPa: 446.1 / 421.1 / 401.6 C, |err| **0.1 %** on all three;
  strict monotone decrease held; slope −0.9265 vs −0.9286 K/MPa (0.2 %).
* **`unit/removal_increment` PASS 15/15** — the assertion restated in §3 keeps
  its strict form (`n_voided == 1`, `h_applied/dz == 1.000000000000`) on the
  default path.
* **`validation/hu/hu_end_to_end` NOT RUN** — the record says it never has
  been, and it is the D2z-1 over-broad-override blast radius. Neither key can
  reach it: `removal_energy = flake` **aborts at parse time** on
  `spallation.model = damage_law`, and `unit/spall_event` (the damage-law
  multi-cell case) is inside the 2639-file gate and byte-identical.

### 6.2 (f) — the Meier `FH` rerun. REPORTED, NOT SCORED.

D2ac's frozen `FH` argv verbatim (lifted through `d2ac_feed/run.py`'s own
`frozen_argv`), plus three keys: `spall.surface_T = reconstructed`,
`weibull.a0_ig = 8.0e-6`, `weibull.a0_gb = 40.0e-6`.

Both legs ran the full 500 s and exited cleanly (`rc = 0`).
✅ **P2's central prediction is confirmed: the axis runs to the block base** —
494.0 mm of 500 mm, against the frozen leg's 332.0 mm.

⚠ **Rates are compared over 60–160 s, not D2ad's 150–400 s**, because the P0
axis is within 6 mm of the base by the end and its late rate is floor-limited,
not physics-limited. 60–160 s is the longest window both legs share with the
axis well clear of the base.

| | frozen `FH` (pre-P0) | `FH` + P0 |
|---|---|---|
| **axis depth at t = 499 s** | 332.0 mm | **494.0 mm of a 500 mm block** |
| axis depth at ≈163 s | 148.0 mm | **276.0 mm** |
| axis ROP, 60–160 s | 2.766 m/h | **3.690 m/h (+33.4 %)** |
| firing temperature, mean | `T_top` 826.5 K | `T_surf` **775.2 K** |
| firing window p10–p90 | **19.4 K** | **86.9 K** |
| `T_surf` − `T_top` at firing | — | **+73.0 K** mean |

⚠⚠ **+33 % IS BELOW THE PACKET'S +50 % STOP RULE, AND THE STOP RULE'S
INFERENCE IS FALSE HERE.** It reads *"a rise below +50 % … means the
reconstruction did not reach the criterion"*. The reconstruction demonstrably
reached the criterion: the 1-D column rises **+216 %** at the same mesh (§4.1),
the control's true surface temperature moves 14 % where the reconstruction's
moves 0.16 % (§2.1), and `T_surf` − `T_top` is +73 K on the Meier firings
themselves. **Meier absorbs the rise through two feedbacks the 1-D column does
not have**: `jet_closure = enthalpy` depletes the gas temperature in proportion
to the heat it delivers, and the `h_expr` map falls off with standoff, which a
faster descent opens — D2ac's *"the model self-corrects by retreating"*. A
stop rule written on an open-loop prediction cannot be applied to a closed-loop
leg.

⚠ **The firing window widened 19.4 → 86.9 K (4.5×).** The packet flags this as
*"a signal worth reporting"* and it is: the pre-P0 criterion pinned the cell
average hard (that is what made the model flux-blind), and the reconstruction
does not pin it — `T_top` at firing now spans **230 K** p10–p90.

⚠ **12.8 % of Meier firings (3113 of 24296) fall back to the cell average** on
the `flat` branch — no usable two-cell gradient. The `ell_eff ≥ dz/4` floor and
the gas ceiling never bind (0 of 3600 columns × 838 thermo steps), so P4 holds
as stated, but the reconstruction is **inert on one firing in eight** on this
configuration, against 0 % on both validation cases.

⚠ **On Meier, `ell_eff` IS partly noise.** Median |`d_ell_rel`| is 3.0e-3 but
the p90 is **0.836** — an 84 % step-to-step jump at the tail, against 3e-4 on
the clean 1-D cases. Firing events disturb the two-cell ratio exactly as the
packet warned, and on a cratered, laterally heated floor they do so often.

**None of this is a physics result**, for the reason in §4: the boundary
condition still drives on the cell average, so part of the +33 % is delivered
power that a surface BC would not have delivered.

## 7. Named debt and what the next packet must own

1. ⚠⚠ **THE BC AND THE CRITERION MUST SHARE ONE SURFACE TEMPERATURE.** They do
   not, in any of the three `robin_form`s (§4(i)), and P0 widened the gap by
   fixing one side of it. Measured rate mesh spread: **7.1 %** pre-P0,
   **48.9 %** with the correct `face` BC, **86.2 %** with the default `cell`
   BC. Until this closes, **no rate from this branch is quotable**, including
   the Meier rerun in §6.2.
   ⚠ The cheapest concrete step is `robin_form = pinned`, which the Meier legs
   already use: `Removal.H:1544` pins the surface at `rev_T_top`, the CELL
   AVERAGE at firing. Under P0 the criterion fires on `T_surf`, so the pin now
   records a temperature the criterion never used. Pinning at `rc_Tsurf`
   instead is a one-line change with a real physical meaning — but it is a
   behaviour change on a frozen-study path, so it belongs in its own packet
   with its own gate, not here.
1b. ⚠ **The `face` form's half-cell profile is LINEAR while the criterion's is
   EXPONENTIAL.** At these rates ℓ_eff = 0.8–1.0 mm sits below dz/2 at 2 mm,
   which is where the two disagree most. Giving the face Newton the same
   exponential profile is the principled fix and is strictly larger in scope
   than swapping one variable.
2. ⚠ **`removal_energy = flake` is mechanically correct and energetically
   inert.** It splits one whole-cell void into several `phi` steps; it does not
   change what enthalpy leaves, because nothing in this code subtracts enthalpy
   on removal. Making it matter needs an explicit energy sink at the firing —
   which is a different packet, and the reason revision 1's version cascaded.
3. ⚠ **`h_f = min(ℓ_eff, dz)` is thick.** Once ℓ_eff < dz the flake average is
   pinned at `(1 − e⁻¹) = 0.632` of the surface rise, mesh-independently — so
   the scheme's steady state sits a fixed **1.58×** above the continuum ablation
   limit `q = ρcvΔT_s` and does **not** converge to it. Whole-cell removal at
   least converges (slowly). This is a definition question for P1/P2.
4. **Feet mode is still unresolved and still owed to P4**, exactly as the packet
   states: calibration here runs `prescribed`/fixed geometry, the rate
   prediction needs `nozzle_descent = feet`.
5. ⚠ `surface.follow_mask = 1` is **required** on any case with active removal
   and a `surface_patch`: without it the patch is applied at the fixed domain
   top, so heating stops dead at the first void. Measured while building the
   1-D case — the new top cell sat at **295.9 K after 59 s** of nominal
   0.55 MW/m². Not a P0 bug, but it silently invalidated the first build of
   this case and will do it again.
6. `beam_void_closure` **T3** — unchanged, still the only red.
