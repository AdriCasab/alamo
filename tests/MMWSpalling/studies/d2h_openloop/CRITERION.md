# D2h pre-registered reading (written before launch, 2026-09-20)

Question: in D2g-1 the 1 mm ring goes silent from the wall inward after ~300 s while the 2 mm
ring does not (100 s block gap 250–350 s: −12.5 %; silent fraction of the 36–40 mm band at
350 s: 1 mm 9 %, 2 mm 0 %; at 450 s 48 % vs 2 %). Is that front (a) local and lateral (cold
rock beside the rim + per-cell criterion), (b) driven by the gas feedback (nozzle path, T_rec,
s_c), or (c) driven by the idle clock (t_cell ∝ dz)?

Pairs (2 mm / 1 mm, 350 s, same keys as D2g L09 except one switch):
- I0: idle clock off. Tests (c).
- RF: recirculation temperature fixed at 1600 K. Tests the T_rec half of (b).
- OL: RF + prescribed nozzle path at L09_2mm's mean feed; free surface off (code constraint).
  Tests (b) as a whole except the s_c → φ path, which remains.

Metrics, computed the same way for every pair (rim.py-style reconstruction from
removal_events.csv; thermo.dat):
1. Pad-band (28–40 mm) spall recession per 100 s block, 150–250 and 250–350 s, 1 mm vs 2 mm.
   For I0 and RF also the burner ROP blocks (nozzle_z fit). For OL the nozzle is prescribed,
   so only the recession metric applies.
2. Silent fraction (time since last firing > 2·t_cell, t_cell = ρc_p dz·528/3e5) by 4 mm band
   at 350 s, and firings per pad column per 100 s.
3. T_rec, s_c, jet_r_reach, pinned share at 350 s.
4. Fields at 300 and 350 s (both meshes): T_top by band, T at the rim-floor level in the
   40–46 mm bands.

Reference: L09 pair (D2g-1) at the same times.

Reading rule (fixed now):
- A switch "removes the front" if the 1 mm 250–350 s pad recession is within 5 % of 2 mm AND the
  36–40 mm silent fraction at 350 s is ≤ 5 % at 1 mm.
- If OL removes the front → the feedback is doing the work; (a) is not sufficient.
- If OL does not remove the front → the front is local/lateral (a); the wall/edge physics is
  the target, and the gas loop is downstream.
- If I0 removes the front and OL does not → the idle clock is the discrete mechanism at the
  edge; make it physical first.
- If RF alone removes the front → the recirculation ratchet is the amplifier; check whether
  the fixed value (1600 K) merely lowered the gas enough to slow both meshes (compare 2 mm RF
  with L09_2mm).
Prior (reviewer, from the 09-20 fields): OL does not remove the front; I0 does not; RF does not.

Not scored against Meier. Not a candidate configuration. Binary 6d047b50… (no build).
