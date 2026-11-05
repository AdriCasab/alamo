# What-if: is the heat LEVEL the only thing wrong on Höser? (run 2026-10-26, user's request, OUTSIDE /plan)

Not a packet step and not a calibration. The heat scale S below is a **what-if lever, never a model setting**.
Base = P2-7 leg HG (`p2_7_hoeser_heat/legs/HG/argv.sh`: Goldstein profile at her measured h, 150 mm block),
with ONLY `surface_patch.h_expr`'s leading constant multiplied by S, `amr.plot_int=3750` (end plotfile only;
output-only key) and `plot_file`. Same binary 9f5aa6de7cf36f52855cf625, no src change. HG itself is S = 1.

Legs (S fixed before any run; 0.581 = 4.3/7.4 = her flux-series range):
  s0581 0.581 | s150 0.15 | s087 0.15x0.581 | s100 0.10 | s058 0.0581 | s034 0.0581x0.581
Pairs scored for the flux exponent n = ln(V_hi/V_lo)/ln(1/0.581): (HG, s0581), (s150, s087), (s100, s058), (s058, s034).

## Pre-registered predictions (written before launch)
1-D receding face, const Höser 2018 props, fires at 850 K, no radiation, mean q = S x 6261.6 x (1866-800)
(scratchpad recede1d_pred.py):
  S 1.0 -> 120 mm | 0.581 -> 75 | 0.15 -> 17.6 | 0.0872 -> 8.3 | 0.10 -> 10.3 | 0.0581 -> 3.4 | 0.0349 -> no fire
  1-D exponents: full-heat pair 0.87 | (0.15, 0.087) 1.39 | (0.10, 0.058) 2.07
3-D expectations (sideways heat loss out of the Ø40 spot makes low-heat legs shallower than 1-D):
  P1  full-heat pair (HG, s0581): n 0.9-1.2 (model-like, ~linear)
  P2  (s150, s087): n 1.3-1.8 ; (s100, s058): n >= 1.5        [her n_q 1.44-1.83, but measured at 47 MPa]
  P3  s100 volume 6-13 cm3 (her 0 MPa 12.9) ; absorbed 20-30 kJ, inside her bound <= 42-47 kJ
      ⚠ weak evidence: S ~ 0.10 was chosen from her depth via the 1-D model, so volume is half-circular.
  P4  depth ratio at r 17.5 mm / axis at s100: <= 0.75 if the heat level alone rounds the bowl (HG 0.886, her 0.61)
  P5  removal rate RISES with time at S <= 0.15 (HG falls 294 -> 108 cm3/min; her 47 MPa time series rises 64 -> 90)
## Decision rule
  P2 + P4 + P5 hold  -> the level is the main missing piece; next = physics that sets the level (nothing fitted).
  P2 holds, P4 fails -> the stepped side faces flatten the floor (c2 finding, G8); fix that before any shape claim.
  P2 fails (n <= 1.15 at low heat) -> uptake itself must depend on flux.

## Result (scored 2026-10-26 against the predictions above; numbers in RESULTS.txt, reducer = reduce_p2_7.one)
| # | measured | verdict |
|---|---|---|
| P1 | full-heat pair n = **0.58** on S (1.08 on absorbed energy): the deeper hole carries the floor away from the nozzle | MISS as written (bar 0.9-1.2 on S) |
| P2 | (s150,s087) **1.31**; (s100,s058) **1.65**; (s058,s034) 4.41 = the minimum flux (s034 onset 21 s, 0.48 cm3) | PASS (her 1.44-1.83) |
| P3 | s100 **13.05 cm3** = 1.01x her 12.9 (bar 6-13: over by 0.05); absorbed **26.3 kJ**, inside her bound <= 41.5-47.4 | volume over the bar by a hair; energy PASS. Half-circular |
| P4 | s100 depth ratio r 17.5 / axis **0.643** (14 -> 9 mm; ±1 mm cells ⇒ ±~0.07); s087 0.54, s150 0.76; HG 0.886; hers 0.61 | PASS |
| P5 | s100 per 5 s window 15 -> 28.5 -> 27.3 (rises, then level); cumulative average 23.9 (15 s) -> 26.1 (30 s), +9 % (hers 64 -> 66, +3 %); HG falls 294 -> 108 | PASS on the cumulative average (her convention); per window rises then levels |
Witnesses: ledger_err_solid max <= 8.4e-4 on every low leg; Ø at 10 mm by sector spread 2.6-2.9 mm on s087/s100 (a 10 mm cut through a
12-14 mm bowl's sloping wall, so 1 mm depth steps move it several mm); flake share falls 92 % (S 1) -> 61 % (s100) -> 43 % (s058).
**Decision rule outcome: P2 + P4 + P5 hold -> the heat LEVEL is the main missing piece.** At her conditions the rock takes in
~26 kJ of the ~240 kJ her calorimeter implies at a fixed face (~1/9). The flux law, the minimum flux, the bowl and the
rate history then come out of the heat equation the model already solves.
⚠ Kant onset passes at FULL measured heat, so whatever cuts uptake must switch on only once spalling has started
(points at flakes/debris shielding the face, Rauenzahn §7.4), not at a wrong calorimeter h.
⚠ Targets n_q and the bowl are her 47 MPa data; the model has no load. Width untested (the r 20 mm lid still holds).
