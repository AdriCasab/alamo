# P1f — lateral flake drive + lateral removal

**Verdict: the lateral mode is implemented, fires, voids exactly one cell and
closes the ledger to machine precision — but criterion (c), its isotropy gate,
is NOT RESOLVABLE at the stated 5 % bar, so no long leg was launched.**

⚠⚠ **The packet's stopping rule is NOT triggered.** It says to build the tracked
front if **(f)** fails on Leg A. (f) was never scored: (q) ruled Leg A inert by
scope and (c) stopped Leg B. Nothing here is evidence about whether a per-face
criterion can remove the √2 bias.

## What the offline gates decided, at zero compute

**(q) — Leg A is inert BY SCOPE, not by soak.** Viability is `ell_lat <=
flake_h_max/c_f = 3.75 mm` (the full-drive condition; the earlier 8.0 mm was the
laxer "capped drive still reaches 77 K" figure — both now reported).

| leg | viable faces | inside `r < 40 mm` | **in scope** | viable radius p50 |
|---|---|---|---|---|
| P1e (Leg A's launch state) | 180 | 160 (**88.9 %**) | 20 (11.1 %) | **25.8 mm** |
| rev7sc (Leg B's launch state) | 402 | 299 (74.4 %) | **103 (25.6 %)** | 38.1 mm |

P1e's viable faces are **floor-corner risers inside the cut** — the 7b branch the
packet excludes for a measured reason. So Leg A would fire ~20 faces for a reason
unrelated to the thesis. ⇒ Leg A not launched.

Per-band, the mechanism is plain: P1e's wall band reads `ell_lat` p50 **55–77 mm**
against rev 7's **3.3–5.7 mm**, i.e. P1e's wall has soaked and rev 7's is an
actively cutting corner. Full tables in `ELL_LAT_PRECHECK.txt`.

## ⚠⚠ (c) — the isotropy gate is below this test's resolution

Measured on one prescribed geometry, varying only the soak time, with removal OFF
so the geometry is static:

| `ell_lat / dz` | axis | 45° | gap | note |
|---|---|---|---|---|
| 4.02 (60 s, shipped) | 8.03 mm | 8.71 mm | **8.4 %** | the scored leg |
| 7.83 (240 s) | — | — | **3.4 %** | |
| 8.66 (960 s) | — | — | 6.3 % | ⚠ the 17 mm profile exceeds the 16 mm domain; the Neumann far side contaminates it |

⇒ **±3 % scatter straddling a 5 % bar, not a systematic anisotropy.** The gap does
not shrink monotonically with resolution, so it is not a pure resolution artefact
either. ⚠ Reporting the 240 s point as a pass would be choosing the soak time that
clears the bar — that is tuning, and the shipped 60 s leg is what is scored.

⚠ Two measurement conditions had to be fixed first, and both were mine:

1. **With removal ON the profile never resolves.** Every firing re-exposes a fresh
   face, pinning `ell_lat` near `dz/4` (0.46–0.53 mm), where a two-cell log ratio
   reads raw ratios of 22–44 and cannot be isotropic. Isotropy must be measured on
   a **static** geometry. ⚠ This is also the real wall's condition: `ell_lat` there
   is 3.75–8 mm at `dz` = 2 mm, i.e. `ell/dz` ≈ 2–4 — the same marginal regime, so
   8.4 % is representative rather than an artefact of a badly chosen test.
2. **The first attempt soaked a cold slab for 2 s**, giving `ell_lat` = 0.28 mm =
   `dz`/7 and an 18.7 % gap. The sequence 18.7 → 13.6 → 8.4 % as resolution
   improved is what identified the cause.

## What does work, measured

| | measured |
|---|---|
| the mode fires (inertness case 1) | **1312 candidates, 240 firings** on a 10-cell riser |
| one cell per firing | **240 voided of 240 fired**; `n_voided` histogram `{1: 240}` |
| `n_voided == 1` invariant | **intact, unrelaxed** — which is why collapse is P1g |
| the face temperature arrives | `lat_no_faceT` = **0** on every face |
| ledger, `removal_energy = cell` | worst **1.47e-15** on all three removal legs |
| a flat slab | **0** firings, **0** regime-6 events |
| firing depths | read from the log: 0, 2, 4, … 14 mm below the column top |

⚠⚠ **The ledger found a latent bug that both Meier legs would have hidden.** The
regime-6 term was first added inside the existing `if (remap_on)` block. The unit
test runs `removal_energy = cell`, so `remap_on` is false, the term was never
added, and `ledger_err_solid` ran to **−1.0**. Both Meier legs run
`removal_energy = flake`, where `remap_on` is true — so the bug would have **passed
(e)** and been silently wrong on every cell-energy run. The term is now
unconditional.

## Carried forward

* ⚠⚠ **(c)'s bar needs a test that can resolve it, or a different statistic.** A
  two-cell log ratio on a 2 mm grid is a ±3 % instrument; a 5 % bar is inside its
  noise. Options: a finer unit mesh (the stencil claim is mesh-convergent, so
  1 mm would halve the scatter), or score the *drive* integrated over a face
  population rather than a median `ell_lat`.
* ⚠ **(q) may overturn the floor-riser exclusion.** P1e's entire viable population
  is inside 40 mm. If the lateral mode is to act where the rock is actually fresh,
  the 7b branch is where it has to act — and the 25 %-per-event floor-rate cost of
  whole-cell removal is then the thing to fix, not to avoid.
* The cell-form branch of `T_face_out` has **no runtime coverage** from the lateral
  path, which only ever calls the face form. It is asserted at source level (three
  write sites, none after the cell-form marker); a direct C++ call is what would
  close it.
* ⚠⚠ No Meier rate, volume, width or standoff is quotable from P1 … P1f.
