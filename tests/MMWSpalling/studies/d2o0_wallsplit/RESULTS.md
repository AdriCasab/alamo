# D2o-0 results — prescribed floor/wall heating (2026-09-22)

**No cylinder appears. P1 and P2 are unscoreable, P3 and P4 fail low, and the prescribed split
does something nobody predicted: it jams the burner.**

- **The tool stalls.** ROP falls from the control's 1.491 m/h to **0.108 / 0.053 / 0.025 m/h**
  (grazing / zero-skirt / h = 20). The nozzle descends ~27 mm and stops; `foot_stall_time` reaches
  219 s. The hole never reaches 10 cm depth, so **P1 and P2 have no shaft to be scored on.**
- **The cause is where the transition was placed, not how strong it is.** The feet rest on the
  **0.9 quantile** of annulus heights, so ~10 % of annulus rock sits above the pad plane — at
  s < the 50 mm stand-off. **The floor/wall switch at s = 45 mm cuts straight through that
  population, 5 mm below the plane the feet themselves define.** Those columns get h = 40, sit at
  **444–482 K**, can never reach the 821 K firing temperature, and hold the burner up for the rest
  of the run. In the control the same columns sit at **680–814 K** — at threshold — and keep firing.
- **The wall is inert, exactly as predicted.** `v(r)` beyond 45 mm is **0.000 m/h** in every
  variant against the control's 1.318, and the band wall ends at **565 K** against a predicted
  578 K. **P5 does not fire** (209 K of margin), so no knife-edge and no mesh-divergence risk from
  that mechanism.
- **Per the packet's decision rule this is the "P1 or P2 fails" branch: report what a prescribed
  field could not produce. The field was not tuned to rescue it**, and must not be.

## 0. Hashes, binary, provenance

- **Binary** `bin/mmwspalling-3d-g++` = `62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`,
  **unchanged at the end**. No source edit, no build, no `make`, nothing committed.
- **`PREDICTIONS.md`** sha256 `855c70e90a491509151e4a630d618423e390b959410a36415bbadc230c5ecc97`,
  read-only, **14:53:01**. All runs completed 15:07–15:11, after the hash.
- **Frozen, imported by path, none edited:** `d2i_gate/` (the control and the LI0 keys),
  `d2j0b_plane/` (the runner), `d2l_scan/analyze.py` (Ø(z), V̇, matched depth),
  `d2m_feet/analyze.py` (D2j-0's v(r) estimators), `d2b_feet_rop/analyze.py`. **The control was not
  re-run.**

### Amendment before any result (full record at the top of `PREDICTIONS.md`)

First hash `9bbf85bc…` at 14:47:16. Amended because **my P2 definition did not reproduce the
packet's own control reference** — it gave −0.5 mm/side where ACTIVE_STEP states ≈33. The packet's
number is the **mouth** gain: the control's Ø at z = 10 mm is 144.9, so (144.9 − 80)/2 = **+32.4**,
an exact match. Found by checking the definition against that reference, **not by looking at any
result**: the runs were killed ~2 min in, `output/` was deleted, no `.done` was written and no
D2o-0 output was inspected. Killed log kept as `campaign_killed_defn.log`.

Two things that check exposed, both on record before results:
- **P1 on the shaft mean is already satisfied by the control** (control shaft mean 90.0 mm, inside
  85–95), so it only discriminates if required of the whole shaft profile.
- **P2 at the mouth is largely a start-up comparison.** The control's *steady* band gain is
  **+5.0 mm/side — already inside P2's 2.5–7.5 target.**

### Deviation: the 1 mm leg was deferred, on the user's prompting and on a stated criterion

ACTIVE_STEP item 4 requires both meshes for the primary. The 1 mm leg was killed ~3 min in (no
`.done`, output removed) because it was projecting to ~6 h against the packet's 2.5 h estimate.
**The criterion for running it was fixed before the 2 mm results were read:** run it if P5 fires
(wall within 50 K of 821 K), or if the result is *positive* and therefore unquotable on one mesh.

**Neither holds.** The wall ends **209 K** clear of the threshold, so P5 does not fire; and the
result is negative. The negative conclusion is also mesh-robust in direction: by the D2j-0b
step-face mechanism, refinement makes lagging columns **colder**, which pushes the frozen annulus
columns *further* from the firing temperature and makes the jam **more** likely, not less.
**Nonetheless this is a single-mesh study and no number in it may be quoted as converged.**

## 1. The result: the prescribed split jams the burner

| case | h_wall | skirt | ROP 150–250 s [m/h] | vs control | nozzle descent in 450 s | `foot_stall_time` at end [s] | `patch_P_robin` [W] |
|---|---|---|---|---|---|---|---|
| control | — | — | **1.491** | — | 236.0 mm | 1 | 1639 |
| `W_2mm` | 40 | grazing | **0.108** | 0.07× | 55.3 mm | 23 | 1137 |
| `Wz_2mm` | 40 | zero | **0.053** | 0.04× | 48.0 mm | 219 | 775 |
| `Wsens_2mm` | 20 | grazing | **0.025** | 0.02× | 50.7 mm | 147 | 1067 |

`patch_min_standoff` grows 60 → 181 mm over the run: **the centre pit keeps deepening** (it stays
in the floor zone and stays impingement-heated) **while the feet are held up**. The tool and the
hole decouple.

### The locking ratchet

| run | annulus columns in the wall zone (s ≤ 45 mm) | their surface T | fate |
|---|---|---|---|
| control | 5 of 162 | **680–814 K** | at the 821 K threshold — keep firing, tool descends |
| `W_2mm` | 15 of 162 | **444–482 K** | 340 K below firing — **frozen permanently, tool jams** |

`foot_carry_cols` = 18 in `W_2mm` against 20–26 in the control, and those carriers cannot spall at
h = 40. **Any annulus column that drifts below the switch freezes and holds the burner up for the
rest of the run.**

**This is a property of the transition's position, not of the coefficient.** All three variants
stall — h = 40 grazing, h = 40 with a zero skirt, and h = 20 — and their hole shapes are identical
to 0.1 mm. Halving the coefficient changes the stall rate, not the outcome.

## 2. Shape — unscoreable, and why

**P1 and P2 are UNSCOREABLE** for every case: the hole never reached 10 cm depth (feet at 40 mm at
end of run), so there is no shaft. The matched comparison depth collapsed to **48 mm**, which puts
every case inside its own start-up funnel — including the control, whose gains read +28.4 / +24.7
there instead of its own +32.4 / +5.0. **The shape table decides nothing in this packet** and is
reported only for completeness.

Ø(z) at the 48 mm matched depth: control 118 / 87 mm at z = 20 / 40; all three variants
**85 / 78 mm**, identical to each other.

## 3. The wall — P3, P4, P5

| case | band cols | T_wall min/mean/max [K] | 821 − T_max | wall heat [W] | P3 (600–800 K) | P4 (1.5–3.5 kW) | P5 |
|---|---|---|---|---|---|---|---|
| `W_2mm` | 239 | 444 / **565** / 612 | **+209 K** | **160** | FAIL low | FAIL low | no |
| `Wz_2mm` | 3298 | 320 / **333** / 441 | +380 K | 11 | FAIL low | FAIL low | no |
| `Wsens_2mm` | 400 | 417 / **496** / 526 | +295 K | 143 | FAIL low | FAIL low | no |
| control | 84 | 589 / **749** / 840 | **−19 K** | — (impinging field) | PASS | — | **fires** |

- **P3 FAILS low at 565 K** against 600–800. But **565 K is within 5 K of the ≈570 K alteration
  threshold**, so the *qualitative* constraint Meier's photograph supports — altered, not spalled —
  is essentially reproduced even though the number misses the band's lower edge. These are
  different claims and are scored separately.
- **P4 FAILS low at 160 W** against 1.5–3.5 kW — an order of magnitude short, and 2.5–5× below even
  my own pre-registered 0.4–0.9 kW.
- **P5 does not fire.** 209 K of margin, so the knife-edge the packet anticipated is absent and the
  step-face cascade is not the operative risk here.
- **The control's own wall sits at 749 K mean, 840 K max — above the firing threshold.** That is
  precisely why the control's hole widens and this one does not.

## 4. v(r) — the cone gets narrower, not flatter

| case | v_c [m/h] | 25 mm | 35 mm | **45 mm** | 50 mm | 60 mm | v(45)/v_c |
|---|---|---|---|---|---|---|---|
| control | 3.595 | 2.696 | 1.848 | **1.318** | 1.118 | 0.658 | 0.367 |
| `W_2mm` | 2.813 | 2.357 | 1.700 | **0.000** | 0.000 | 0.000 | **0.000** |
| `Wz_2mm` | 1.970 | 1.666 | 1.395 | **0.000** | 0.000 | 0.000 | **0.000** |
| `Wsens_2mm` | 2.669 | 2.301 | 1.680 | **0.000** | 0.000 | 0.000 | **0.000** |

**Q6 confirmed.** Erosion beyond 45 mm is exactly zero. Meier's implied floor half-width is
42.5–46.5 mm — precisely where this field stops removing rock. **A single grazing coefficient on
the wall does not turn the cone into a cylinder; it truncates the cone.**

## 5. Pre-registration scorecard

| prediction | outcome |
|---|---|
| **Q1** wall does not spall, P5 does not fire | **HELD** — 565 K, 209 K of margin, v(r>45) ≡ 0 |
| **Q3** shaft Ø ≈ 80 mm, not 85–95 | **unscoreable** (no shaft); Ø at the depths reached is 78–85 |
| **Q4** P3 marginal, 550–650 K | **HELD** — 565 K, and the ≈570 K alteration threshold is straddled |
| **Q5** P4 fails low, 0.4–0.9 kW | **direction HELD, magnitude wrong** — 160 W, 2.5–5× below my own estimate |
| **Q6** v(r) steeper, not flatter | **HELD** — v(45)/v_c 0.367 → 0.000 |
| **Q2** mouth gain +5…+15, steady 0–2 mm/side | **unscoreable** — the comparison depth collapsed to 48 mm |
| **the stall** | **NOT PREDICTED.** I predicted an inert wall; I did not predict that an inert wall would jam the descent |

The semi-infinite Robin estimate in `PREDICTIONS.md` §4 (578 K) came out within **13 K** of the
measured 565 K, which is why the wall-side predictions held and is a useful check on the estimate
itself.

## 6. Ledger

Max abs `ledger_err` ≤ 5e-15 on every run; energy closes. J/mm³ to rock is reported against the
thermodynamic floor 1.147 only — **the 4–6 J/mm³ band is an assumed 30 % delivery, not a
measurement, and was not used as a target.**

## 7. What the next packet inherits

1. **A hard design constraint on the D2o closure, which is the real product of this packet: the
   floor/wall transition must not fall inside the height spread of the support annulus.** The feet
   rest on the highest rock under the pads; if the transition sits within that spread, the carrying
   columns land on the cold side and the tool locks. Here the switch (45 mm) sat 5 mm below the
   stand-off (50 mm) and the 0.9 quantile spans exactly that gap. A gas-path closure must either
   place the transition well below the pad plane, or make it smooth enough that no carrying column
   can freeze.
2. **A single grazing coefficient cannot produce Meier's cylinder** — it zeroes erosion beyond
   45 mm and truncates the cone. The widening mechanism Meier's hole requires is still missing, and
   it is not "less heat on the wall".
3. **The wall temperature target is nearly reachable**: 565 K against an alteration threshold of
   ≈570 K and a spall threshold of 821 K. A closure that puts the wall in 600–800 K would reproduce
   "altered, not spalled" without spalling it — but on this evidence that alone will not widen the
   hole, because at those temperatures the wall still does not fire.
4. **Owed:** the 1 mm leg, if any positive result is ever claimed from this configuration; and the
   2639-file sweep before any rebuild (the tree passed it on 09-22, so a build is permitted, but
   this packet needed none).
