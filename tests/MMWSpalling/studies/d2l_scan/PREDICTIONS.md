# D2l Gate 0 — predictions and decision rule
# (written and hashed BEFORE any run; never edit)

Key-only. No source change, no build. Binary `bin/mmwspalling-3d-g++` =
`62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa` (the D2k campaign build).

**`f` (`surface_patch.side_face_factor`) is a diagnostic dial in this packet, not a calibration.**
It scales `q_side` linearly at exactly the place an `h` scale would, which is why it maps the
window. **Its physical meaning is still AREA. No value found here may be kept as a calibrated
setting** — a closure value must be derived from a named heat-transfer source.

## 0. What is already known (not predicted)

| quantity | f = 0 (D2i, side off) | f = 1 (D2k) |
|---|---|---|
| burner ROP, 150–250 s | 1.491 m/h | 4.066 m/h |
| absorbed power (quarter) | 1638.6 W | 3204.2 W |
| `jet_s_c` (centre standoff) | 88 → 166 mm over the run | pinned 67–69 mm |
| `jet_fs_cols`, 150–250 s | 17.9 | **0.0** |
| arbiter mesh gap by block | −13.5 / −28.7 / −45.5 % | +3.4 / +1.5 / −3.2 % |

## 1. Deviation from ACTIVE_STEP, declared before launch: 0c is run at BOTH ends

ACTIVE_STEP puts 0c at f = 1. **Both probes are inert there**, and the table above is why — these
are properties of the geometry the side flux produces, not of the keys:

- the far law is `pow(12D/max(12D, s), n)` with 12D = **90 mm**. At f = 1 the hole is self-similar
  with `jet_s_c` 67–69 mm and `patch_min_standoff` 63 mm, so **every column sits in the clamp
  region where the factor is exactly 1 for any `n`**;
- `jet_fs_cols` = 0.0 and `m_ratio` = 1.000 through the entire scoring window, so the free-surface
  dilution is **already off**.

With the side key off both are live (`s_c` past the clamp, `fs_cols` = 17.9), and that is also the
configuration the competing explanation describes — the runaway centre pit with a skirt. So each
probe is run **at f = 1 as specified** (to settle inertness by byte-identity rather than by
argument) **and at f = 0** (`P_fsoff0`, `P_flat0`), where it can move something. No new key is
used. Nothing from the packet is dropped.

## 2. Pre-registered predictions

### 0a — `f_min`, the smallest f at which the cascade stays suppressed

The step-face sink is `k(T_top − T_wall)/dx` ≈ **240 kW/m² at 2 mm** (480 at 1 mm) against a pinned
top flux `q_pin` ≈ 520 kW/m². A side face is given the same `h` and the same gas temperature as the
top face, so it supplies ≈ `f · 520` kW/m². Suppression needs supply ≳ sink:

- **P1. `f_min` = 0.5 at 2 mm.** f = 1.0 and f = 0.5 suppress (all deciding blocks within 5 %, no
  inward `r_front`, no band stalling); **f = 0.2 and f = 0.1 do not** — band stalling reappears at
  2 mm as the single-mesh signature, and the boundary pair confirms it on the mesh gap.
- **P2. `f_min` is itself mesh-dependent, and rises with refinement.** The sink scales as 1/dx and
  the supply does not: 240 vs 520 at 2 mm, 480 vs 520 at 1 mm. So the 1 mm leg at f = 0.5 is
  predicted to be **worse than its 2 mm partner**, and f = 1 is predicted to be marginal at 1 mm and
  **to fail at 0.5 mm**. *(The 0.5 mm leg is 5–7 h and is NOT run in Gate 0; it is named here so
  the prediction is on record before the queued overnight run.)* **If P2 holds, D2k's arbiter PASS
  is a two-mesh statement about a remedy whose sufficiency is mesh-dependent** — the week's
  recurring failure mode — and that is a finding independent of the rest of this gate.

### 0b — `f_rate`, the f at which the burner rate re-enters 1.3–1.6 m/h

Absorbed power runs 1639 → 3204 W (2.0×) and the rate 1.49 → 4.07 m/h (2.7×) across f = 0 → 1.

- **P3. `f_rate` ≈ 0.1, and certainly `f_rate` < 0.2.** f = 0.1 lands in or just above the band;
  f = 0.2 is above it; f = 0.5 and f = 1.0 are far above it (> 2.5 m/h).
- **P4. Ø at matched feet height is nearly flat in f**, changing by **< 10 mm** across the whole
  ladder. D2k already showed the term is a feet-descent accelerator, not a hole-widener (+1.4…+4.0
  mm at matched feet height at f = 1). So **cutting f buys rate without buying shape** — the D2b
  failure mode restated. Ø is predicted to stay **below 85 mm**, i.e. below Meier's band, at every f.
- **P5. J/mm³.** Two conventions, and they must not be mixed. Meier's **15.4 J/mm³ is
  combustion-enthalpy** (38 kW HHV ÷ cavity volume); the model absorbs ~3.2 kW, so scoring absorbed
  power against 15.4 would be wrong by about an order of magnitude. Reported as both:
  **TSE(HHV) = 38 kW / V̇** against 15.4, and **TSE(to-rock) = 4·`jet_P_face` / V̇** against the
  reference's 4–6 J/mm³ band. Predicted: TSE(HHV) ≈ **4–5 J/mm³ at f = 1** (far below 15.4, i.e. the
  model removes far too much rock per unit fuel) rising toward ≈ **14–15 at f = 0.1**;
  TSE(to-rock) ≈ **1.4–2.5 J/mm³ at every f**, i.e. **below the 4–6 band throughout**, never more
  than ≈ 2× the thermodynamic floor ρc_p ΔT_fire = 1.147 J/mm³. **So the to-rock band fails at every
  f, including the f that fixes the rate** — the rate is matched by delivering less power, not by
  removing rock at the right cost.

### 0c — the competing explanation

- **P6. `P_flat` is byte-identical to `M_f10`** (thermo.dat identical to the last row). Falsified by
  any difference at all.
- **P7. `P_fsoff` differs from `M_f10` only through the early transient**, with the 150–250 s rate
  within **5 %** and Ø at matched feet height within **5 mm**.
- **P8. At f = 0 the probes are live but do NOT show the discriminating signature.** `P_flat0`
  (weaker s-decay ⇒ *more* flux at the deep runaway centre, which is where s is largest) raises the
  rate and **narrows** the hole; `P_fsoff0` (no ambient dilution ⇒ hotter gas at the skirt) raises
  both the rate and Ø. **Neither gives rate DOWN with Ø UP.**
- **P9. Therefore the distribution knobs that exist today cannot produce the signature**, in either
  direction of the dial. If any run does show rate down *and* Ø up together, P8 is refuted and the
  distribution becomes the stronger lead, per the packet.

### The decision

- **P10. The window is EMPTY: `f_rate` ≈ 0.1 < `f_min` = 0.5.** The same parameter controls
  cascade suppression and the power excess and they pull opposite ways, so **no single `f` satisfies
  both gates**. Predicted outcome: **STOP — write no code, do not run Gate 1, do not write the
  `side_face_h_scale` closure.**

**The decision rule, from ACTIVE_STEP, applied to whatever comes out:**

| outcome | action |
|---|---|
| `f_rate ≥ f_min` and 0c shows nothing | explanation 1 — proceed to Gate 1 and the closure |
| `f_rate < f_min` | **STOP. Write no code.** Report `f_min`, `f_rate`, 0c, and name the next packet |
| 0c shows the signature (rate down *and* Ø up) | report it as the stronger lead and stop before the closure, even if the window exists |

## 3. Scoring definitions (fixed here, before any result is seen)

- **Burner ROP** = least-squares slope of `nozzle_z` over thermo rows in [150, 250] s, m/h — the
  D2f/D2g/D2h/D2i/D2k estimator, unchanged.
- **Cascade suppression (0a)**, the D2k `CRITERION.md` conditions, unchanged: disk-mean rate
  (r < 30 mm) gap within **5 %** in every block 50–100, 100–150, 150–200 s; no band inside 28 mm
  stalling with Δ_own > 10 s (a band with `t_stop ≥ t_end − 5` is still firing, not stalled); no
  inward-walking `r_front` at 1 mm.
- **`f_min`** = the smallest f on the ladder {1.0, 0.5, 0.2, 0.1} meeting all three.
- **Ø at matched feet height** = 2 × azimuth-mean `r_wall` (9 sectors, the D2c `score.py`
  definition) at **one common feet depth reached by every run on the ladder**, chosen as the
  shallowest of their end-of-run feet depths. Never at matched time.
- **V̇** = 4 × Σ `h_applied`·dx·dy over the removal rows in [150, 250] s, ÷ 100 s.
- **`f_rate`** = the f whose 150–250 s ROP falls in [1.3, 1.6] m/h, by linear interpolation in f
  between the bracketing runs if none lands inside.
- No other window, block length or offset decides anything. **Nothing is re-windowed after results
  are seen, and `f` is not tuned.**
