# D2k arbiter criterion — side-face flux against the D2j-0b C1 cascade
# (fixed before launch; never edit)

**Runs.** `C1_2mm` / `C1_1mm`: the D2j-0b **C1** keys verbatim (static disk under a descending
exclusion plane, h = 700 below it and 1e-3 above, T_gas = 1600, quarter domain 0.08 × 0.08 ×
0.15 m, stop 250 s, dt 16 ms at 2 mm and 8 ms at 1 mm, `blocking_factor = 1`, V0 = dz³) plus
**`surface_patch.side_face_flux = 1`, `side_face_factor = 1.0`** and nothing else.
`C1f07_2mm` is the same at f = 0.7. Binary: the D2k build (sha in `RESULTS.md` §0).

**f is pre-registered at 1.0, with 0.7 as a sensitivity. It is never fitted, and the arbiter
verdict is read off f = 1.0 alone.**

**Metrics — the D2j-0b definitions, reused unchanged** (`d2j0b_plane/CRITERION.md` §Metrics,
`analyze.py`):
- **disk-mean rate**: mean over columns with r < 30 mm of the per-column recession rate from
  `<pf>_removal_events.csv`, per block. **Gap = 1 mm / 2 mm − 1.** This is deliberately the
  **same estimator as the Meier gate**, so a pass transfers.
- **t_stop(band)**: median over a 2 mm band's columns of each column's last removal time
  (columns still firing in the last 5 s count as t_end). **t_cross(band)**: median over the
  band's columns of the first time s(t) ≤ 0. **Δ_outer = t_stop(band) − t_cross(next band
  outward)**.
- **r_front(t)**: inner edge of the contiguous run of 2 mm bands reaching 30 mm whose mean rate
  over a trailing 25 s is < 0.1 v_c.

**Numbers to beat (D2j-0b, key off).** Disk-mean gap by block 50–100 / 100–150 / 150–200 s:
C0 −10.4 / −15.3 / −20.1 %; **C1 −13.5 / −28.7 / −45.5 %**. C1 `r_front` at 1 mm walks
26 → 22 → 18 mm.

**PASS requires all three:**
- **(a)** the disk-mean gap is within **5 %** in **every** block 50–100, 100–150, 150–200 s;
- **(b)** Δ_outer is mesh-independent within **20 %**, **or no band stalls at all** (the better
  outcome: no band inside 28 mm stops with Δ_own > 10 s at either mesh);
- **(c)** **no inward-walking `r_front` at 1 mm** — r_front is either absent or does not
  decrease across the three sample times.

**Otherwise FAIL.** No other window, block length or band width is evaluated for the verdict.

**Reported, not deciding:** total side power as a fraction of top-face power per block
(`side_P_face` against `patch_P_robin`); the f = 0.7 sensitivity; the band table in the D2j-0b
format, so it sits beside C0/C1/C2; `side_faces`; `ledger_err`.

**Known bias, stated in advance and not to be tuned away:** a staircase of one-cell steps has
lateral area dx + Δz per column against a true sloped √(dx² + Δz²), so f = 1 over-heats the
wall by up to 41 % at 45° and ≈ 16 % at the observed ~80° cone.

**Outcomes:**
- **PASS →** run the conditional Meier pair `LI0_2mm` / `LI0_1mm` (the D2i keys with the key on,
  f = 1, to 450 s) under the **D2i criterion** — matched window from 150 s and every 100 s block
  within 5 % — hashed separately before that launch, with the direction pre-registered in
  writing first.
- **FAIL → stop.** Do not run the Meier pair, do not tune f, and name what the residual looks
  like. The named fallback is the isotherm front **plus** side-face flux, for a later packet.
