# D2k conditional Meier pair — criterion and pre-registered direction
# (written and hashed BEFORE launch; never edit)

Run only because the arbiter passed (`CRITERION.md` sha256 `6322064a…`; analysis in
`analysis.md`: blocks +3.4 / +1.5 / −3.2 % against the key-off C1 −13.5 / −28.7 / −45.5 %,
no band stalls, no inward front).

**Runs.** `LI0_2mm` / `LI0_1mm`: the **D2i `LI0` keys verbatim** (D2f Q09 support rule — feet,
`foot_rule = pads`, `foot_pad_quantile = 0.9`, `foot_body_clearance = 0` — with
`pinned_idle_cycles = 1e9`), loaded through `d2i_gate/run.py`, plus
**`surface_patch.side_face_flux = 1`, `side_face_factor = 1.0`** and a stop of **450 s**.
Nothing else changes. Binary: the D2k build (sha in `RESULTS.md` §0).

## The gate — the D2i criterion, unchanged

- **Burner ROP** = least-squares slope of `nozzle_z` over the thermo rows with a ≤ t ≤ b,
  t > 0, in m/h — the same estimator as D2f, D2g-1, D2h and D2i.
- **t_m** = min(t_end of the pair). **Gap** = ROP(1 mm)/ROP(2 mm) − 1.
- **PASS requires both:**
  - **(a)** the matched-window gap over [150 s, t_m] within **5 %**;
  - **(b)** every **100 s** block from 150 s within **5 %**; the trailing block is dropped if
    shorter than 50 s.
- **Otherwise FAIL.** No other window, block length or offset decides.

D2i's numbers to beat, same keys with the side key off: matched window 150–572 s **−13.9 %**;
100 s blocks **−2.1, −4.0, −16.3, −36.1 %**. To 450 s the comparable blocks are −2.1, −4.0,
−16.3 %.

## Pre-registered direction (written before the run)

With the budget binding — side power is drawn from the enthalpy march, so `jet_P_face` rises
toward the same `jet_P_cap` — the side term **moves energy from the pit to the wall**. Therefore:

1. **The centre rate goes DOWN somewhat.** The cap is the same, so what the wall absorbs the
   face cannot. Expect the 2 mm burner ROP over 150–450 s to fall by **a few per cent**
   relative to D2i's 1.49–1.51 m/h, not to rise. **A rise would mean power is entering from
   outside the budget** and must be investigated, not reported as an improvement.
2. **The hole WIDENS at the feet depth.** This is the quantity D2d failed
   (min Ø 77.6 mm to the feet, against Meier's 85–93 mm).
3. **`jet_P_face` ≤ `jet_P_cap` on every block** (the in-code invariant), with the margin
   reported per block.

### Hand estimate of the widening, before the run

The arbiter measured side power at **7 % of top-face power** (6.9 % at 2 mm, 7.6 % at 1 mm).
Taking the same fraction on the Meier whole-excavation rate of ≈ 2.5 cm³/s gives a lateral
removal of ≈ **0.18 cm³/s**. The exposed wall below the nozzle plane in the quarter domain is
roughly 2πr·h_active/4 with r ≈ 42 mm and h_active ≈ 50 mm (nozzle plane to floor), i.e.
≈ 3.3e-3 m². That is a radial recession of

    v_r = 0.18e-6 / 3.3e-3 = 5.5e-5 m/s = 0.055 mm/s,

and a given depth stays in the active zone for h_active/ROP ≈ 50 mm / 0.42 mm/s ≈ 120 s, so

    **Δradius ≈ 6.5 mm, i.e. Δdiameter ≈ 13 mm at the feet depth.**

**Pre-registered as an order of magnitude, 5–15 mm in diameter**, on the explicit caveat that
the 7 % fraction was measured on the arbiter's small static disk and the Meier hole has more
wall below the plane, so the true fraction may be larger. On D2d's 77.6 mm this would give
≈ 85–93 mm, which is Meier's band — **that coincidence is a prediction, not a target, and the
hole is reported, never scored, in this packet.**

## Outcomes

- **PASS →** the support rule is mesh-converged to t_m at 0.40 m with the side-face term, and
  D2g Stage 2 can be planned on it. Report the widening against the estimate above.
- **FAIL →** stop. Report the residual gap, its block structure, `jet_P_face` against the cap,
  `side_P_face` as a fraction of top, and where the remaining mesh dependence sits. **Do not
  tune f.** The named fallback is the isotherm front plus side-face flux, for a later packet.

Not scored against Meier's ROP here. Nothing from D2b–D2d may be quoted as a Meier result.
