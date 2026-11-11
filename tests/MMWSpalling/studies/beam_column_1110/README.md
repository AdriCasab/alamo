# Beam column resolution study (run 2026-11-10/11, directly at the user's request; not a packet)

Binary `e705b8591d083230` (P2-31's, implemented, not yet reviewed, uncommitted). Fixture, legs and keys: `run.py`.
Tables: `RESULTS.txt` (`reduce.py`). Six points (absorption depth D 0.5, 1, 2, 4, 8 mm at 2 MW/m2 absorbed; D 2 mm at
4 MW/m2), six cells (3.125 ... 0.098 mm), two legs each (A the model as it is, H heat only). 68 legs rc 0; four stopped
on purpose (`SKIPPED.txt`). `diag/`: D = 1 mm, 2.5 s, six cells, the buckling log on.

## Reading

1. THE HEAT SOLUTION CONVERGES (H). Error against the extrapolated limit, by cells per absorption depth dz/D:
   top-cell temperature at 1 s: 8-20 % low at dz/D 1-2, 2-8 % at 0.5-1, <= 2 % at 0.25-0.5, <= 0.4 % below 0.25;
   melt onset: 35-71 % LATE at 1-2, 8-23 % at 0.5-1, 1-12 % at 0.25-0.5, <= 4.4 % at 0.125-0.25, <= 1.7 % below.
   At the regime map's 1.5625 mm cell the melt onset is late by 131 / 54 / 17 / 7 / 4 % at D 0.5 / 1 / 2 / 4 / 8 mm.
   Working rule: dz <= D/4 for the temperature (2 %), dz <= D/8 for the melt onset (2-4 %).
2. THE ONSET CRITERION CONVERGES, FIRST ORDER, FROM THE EARLY SIDE. The reconstruction (State.H, RawSurfaceX /
   ReconFromX) extrapolates an exponential from the top two cells; under a beam the profile is flat at the surface, so
   on a coarse cell it overshoots (D 0.5 mm at 3.125 mm: top cell 477 K, "surface" 849 K). Criterion met 22-33 % early
   at dz/D >= 1, 11-19 % at 0.5-1, 4-11 % at 0.25-0.5, 1-5 % at 0.125-0.25.
3. THE THERMAL LENGTH HAS NO CONVERGED VALUE. ell_eff = dz / ln(rise_top / rise_1) on a profile with zero slope at the
   surface goes as 1/dz (D 1 mm: 1.45, 1.62, 2.15, 3.5, 7.0, 19 mm at the six cells; `diag/`). The flake is
   c_f ell_eff and the buckling condition reads it, so every firing at onset is withdrawn once the cell is small.
4. THE CLASS IS NOT MONOTONE IN THE CELL. Spall / melt by cell, coarse to fine:
   D0.5 Q2  S S S S S M      D1 Q2  S S S M M S      D2 Q2  S M M M M M
   D4 Q2    M M M M M -      D8 Q2  M M M M M -      D2 Q4  S S M M S S
   Fine-cell firings are late and hot (1,100-1,620 K, some inside the fusion band 1,410-1,610 K), and at the finest
   cells they come through the FALLBACK branch (surface cell cooler than the one below, the face losing heat), which
   takes the profile's own e-fold depth ell_p instead -- a mesh-independent length (D 1 mm: 2.14-2.17 mm at every
   cell from 0.78 mm down).
5. ONLY "D >= 4 mm melts" is the same at every cell. The spall side of the map is not established at any cell.

Not run, a reading of logged numbers only: with ell_p in the buckling condition (bk = 0.428 x ell[mm] at the onset
temperature, from the logged pairs), D 0.5 and 1 mm pass at onset (0.59, 0.92), D 2 mm is withdrawn (1.3-1.5), D 4
and 8 mm are withdrawn (2.7, 5.0).
