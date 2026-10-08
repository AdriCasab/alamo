# D2x criteria — fixed and hashed BEFORE the first leg launched

## 0. Disclosure

**No leg (`C0`, `C`, `D`, `E`) has been launched when this file is hashed.**
What has run: the new `side_face_area` regression and the 2639-file tree sweep,
both recorded in `PREFLIGHT.md`. `A`, `B`, `M0`, `M1` are **imported by path and
never re-run**; every number this packet rests on was reproduced exactly.

## 1. Scoring is at EQUAL BURNER DEPTH, never at equal time

Each leg is sampled at the first time its burner descent
(`nozzle_z[0] − nozzle_z(t)`) reaches the target, and **only depths at least
50 mm above the pads are scored** — the pads start at the rock surface, so the
cap is `descent − 50 mm`. Rock the nozzle plane has not passed is unfinished
pit, not wall (standing rule, D2o-0). **Equal-time scoring across legs with
different rates compares unfinished pit against finished wall — the error D2w
made.**

* **Primary: descent 150 mm**, depths **10/25/50/75/100 mm**.
* **Secondary, reported where reached: descent 180 mm**, depths
  **10/25/50/75/100/125 mm**. ⚠ `A` reaches 180 mm only at 491.9 s and tops out
  at 183.3 mm, so any leg slower than `A` cannot be scored there. **Not gated.**
* **No "common end time" anywhere.** `D` is expected to stop early on the bottom
  watchdog; equal-depth scoring is time-free.

Metric imported unchanged from `d2v_hmap/analyze.py`: fixed physical 2 mm radial
bins, diameter = 2 × interpolated crossing radius. **Pit lead and volume are
reported at the PRIMARY target**, because both move with the sampling point
(`A`: lead 128 mm and 1.24× at 150 mm, but 136 mm and 1.19× at 180 mm).

## 2. On-disk references, reproduced this session

| | `M0` | **`A`** | `B` | `M1` |
|---|---|---|---|---|
| (a) at descent 150 mm | 28.40 | **2.27** | 19.93 | 2.73 |
| (a) at descent 180 mm | 25.12 | **2.51** | 17.36 | 2.55 |
| burner descent [m/h] | 1.502 | **1.296** | 3.087 | 3.123 |

`A`'s four defects: Ø **84.4** at 100 mm and **82.7** at 125 mm against an 87 mm
lower bound; rings 40–48 mm at **−29 / −59 / −87 / −96 %** of `M0`; pit lead
**128 mm**; volume **1.24×** on the descent basis.

## 3. Criteria

* **(a) SCORED — shape at equal burner depth.** ⚠ **`A` IS THE CONTROL THAT
  ALREADY PASSES; the bar is not "beat `M0`".** A leg succeeds only by closing
  `A`'s wall-band narrowing **without losing what `A` has**, i.e. ALL FOUR of:
  1. **Ø at 100 and 125 mm depth within 2 mm of 87** (`A`: 84.4, 82.7);
  2. **(a) no worse than `A`'s** at the same descent (≤ 2.27 at 150 mm);
  3. **burner descent in 1.3–1.6 m/h** (`A`: 1.296, just below);
  4. **pit lead ≤ `A`'s 128 mm** and **volume closer to 1.0× than `A`'s 1.24×**
     on the descent basis.
* **(b) REPORTED — volume on BOTH bases**, pit and descent, against Meier
  3.42 L / 0.5 m pro-rated. A leg is a **candidate** only if it lands within
  0.75–1.25× on at least one basis.
* **(c) REPORTED — burner descent**, least-squares slope of `nozzle_z` over the
  **100–150 mm descent interval**, so every leg is measured over the same
  physical stretch. **Never pit depth.** Report **pit lead** at both targets.
* **(d) REQUIRED WITNESS — ring-resolved recession, 2 mm rings, r = 30–50 mm**,
  each leg against `A` and `M0`. ⚠ **D2w's three required witnesses were all
  blind** because the pad-ring witness averaged healthy inner rings against the
  collapsing outer one. **Per ring, printed, never averaged.**
* **(e) REQUIRED WITNESS — power split and the realised weight.**
  `side_P_face`, `side_P_loss`, `patch_P_robin`, `side_faces`, `side_w_mean`,
  and the identity `jet_P_face = patch_P_robin + side_P_face + side_P_loss` to
  **±0.01 W on every leg**.

## 4. Pre-registration

* **P1** `D` burner descent **> 4.0 m/h** (refuted below 3.5).
* **P2** `side_w_mean` on `D` in **0.55–0.75** (refuted outside 0.45–0.85). This
  is the derived number that replaces `f = 0.2`.
* **P3** `C` restores the collapsed rings: recession at 42–46 mm within **20 %**
  of `M0`'s (refuted if still below −50 %).
* **P4** `C` closes the wall band — Ø at 100 and 125 mm (descent 180) **within
  2 mm of 87**, against `A`'s 84.4 and 82.7 (refuted if it gains < 1 mm at either).
* **P5** `C` burner descent in **1.3–1.9 m/h** (refuted above 2.4).
* **P6** `E` lies between `C` and `D` on descent.

**⚠ P4 and P5 together are the packet.** If both hold, the wall band was a
drain, `f` is retired, and the model matches Meier's shape and rate from derived
quantities alone. If P4 holds and P5 fails, wall support and rate are
inseparable through the same heat and the next question is the **skirt**.

## 5. Inertness witnesses

| failure | witness | bar |
|---|---|---|
| `r_min` never reaches the block | `side_faces` on `C` vs `M1` | `C` ≥ **30 % fewer** |
| `area_exact` never reaches | `side_P_face` on `D` vs `M1` | `D` ≥ **2×** `M1` |
| weight formula wrong | `side_w_mean` | in **[0.4142, 1.0]**; 0.2 = inert |
| above-plane leakage | `side_P_face_above` | **exactly 0.0 W** on every leg |
| ledger twin not mirrored | `ledger_err` | **≤ 1e-13** on every leg |
| new keys not default-off | `C0` vs `M1` | `thermo.dat` + `removal_events.csv` **byte-identical to 200 s** |

⚠ **`C0` is mandatory and the tree gate cannot replace it**: no input in the
2639-file set turns `side_face_flux` on, so 2639/2639 says nothing about whether
these keys are inert on a side-path run.

## 6. Guardrails binding on the scoring

* **Do not touch the map.** `h_expr` copied verbatim from `A.log` and asserted
  byte-identical to `M1.log`'s 751-character token before launch.
* **`jet_above_h` / `jet_T_above` stay unset** — above-plane heating is `M2`'s
  territory and the march does not deduct it from the jet.
* **`side_face_factor_above = 0.0` on every leg, explicitly.**
* **`f = 0.2` frozen, never scanned.** `D`/`E` must not set it at all.
* **No 1 mm run, no convergence claim. No Meier rate, diameter, volume or
  energy claim** — Fig. 8.8 is a lower bound at ±5 mm, fuel is methane/air at
  λ = 1.2.
