# D2z-2 RESULTS — the three measurement defects, the conductance, and Meier

`CRITERION.md` sha `46571c3c884b976a4bbb56a45f571a1ff4217c085a5265d76bb985a096b4c118`,
hashed with **zero outputs on disk**. Full scoring in `analysis.txt`.

**Headline: the repairs are inert on physics, the arithmetic fix does NOT
reopen the rim D2z cured, and `wetted` wins on Meier — but `SW` over-delivers
rate by 2.4x the Meier band, and `SW0` proves that is not the conductance form.**

## Scorecard

| criterion | result |
|---|---|
| (a) repairs are diagnostics-only | **PASS** — `R2` byte-identical to frozen `S2` |
| (b) the estimator | **FAIL as pre-registered**; **−6.62 % (inside the bar)** under the measured live population |
| (c) THE RIM TEST | **PASS both discs** — r<25 `+1.39 %` (0.00 pts from frozen), r<30 `−0.40 %` |
| (d) drain meter self-consistency | **FAIL as written** (flagged in `CRITERION.md` before launch); exact identity PASSES |
| (e) power identity | **PASS on all four Meier legs** — 0.0134–0.0141 W against 0.02 |
| (f) identity | **PASS** — header delta EMPTY; tree **2639/2639** vs `d2z1.md5` |
| (g) Meier shape | reported below |
| (h) law A/B | **DISTINGUISHED: `wetted` wins**, 2.80 vs 4.90 mm (bar 0.25) |
| (i) conductance regression | reported — worth **0.08 mm and 0.6 %** |
| (j) multi-cell witness | **f_multi 0.73 on `SW`** — the branch is finally exercised |

| prediction | verdict |
|---|---|
| **P1** `R2` ≡ `S2` | **HELD** |
| **P2** live-vs-offline ≤ 10 % | **REFUTED** — population, not physics (below) |
| **P3a** derived inert at 2 mm | **REFUTED**, with a measured mechanism |
| **P3b** `R1d` ≠ `S1` | **HELD** — 76,797 vs 76,558 rows |
| **P4** the rim holds | **HELD** — under both readings of the bar |
| **P5** `f_multi` ≥ 20 % on `SW` | **HELD** — 73.2 % |
| **P6** `SW` Ø@100 < 87 mm | **REFUTED** — 87.6 mm at 150, 88.2 at 180 |
| **P7** `SW` descent 1.3–1.6 m/h | **REFUTED** — 3.694 m/h |

---

## 1. The repairs changed no physics — three independent proofs

* **(a)** `R2_removal_events.csv` is **byte-identical** to the frozen
  `d2z_stepheight/output/S2_removal_events.csv`. Deliverables 1–3 are
  diagnostics-only, which is what licenses scoring `R2d`/`R1d` against the
  frozen legs.
* **(f)** `S0` is byte-identical to `d2z1_delete/output/RB_D` — 89,775 removal
  rows, 1,732 jet-profile rows, **thermo header delta EMPTY**. None of the eight
  new columns appears with the keys off.
* **Tree gate: 2639/2639 byte-identical to `witness/d2z1.md5`.** Zero differing
  hashes, no path added or removed. `guarded_post_all.sh` confirmed the binary
  never moved (`GUARD: binary unchanged across the sweep`).

## 2. Deliverable 2 is vindicated — the drain meter was reading the patch edge

`R2`, ungated **0.8427** against in-gate **0.1967**. D2z's 0.843 was *exactly*
the patch-edge contamination the packet diagnosed: cliffs whose outer column has
never fired enter at `g = 1` by construction and outnumber the real faces
(`n` 70.2 against `n_in` 13.0). The real in-patch reduction is **~5x**, not the
1.19x D2z reported.

## 3. THE RIM TEST — P4 HELD, and this is the packet's point

| disc | `R1d` | `R2d` | gap | frozen `S1`/`S2` |
|---|---|---|---|---|
| r < 25 mm | 1.6225 | 1.6003 | **+1.39 %** | +1.39 % |
| r < 30 mm | 1.5932 | 1.5996 | **−0.40 %** | −0.71 % |

The interior reproduces the frozen gap to **0.00 pts**, exactly as predicted from
`n = 1` being exact there. At the rim the gap is **−0.40 %**, against a bar of
5 % under either reading (own gap 0.40 %, shift from frozen 0.31 pts).

**⇒ D2z's rim cure was the population and drain correction. The double count was
not doing the work.**

## 4. P3a REFUTED — the derived form is NOT inert at 2 mm

`R2d` differs from `R2`/`S2` from row 343. The packet's basis was that at the
crater wall a never-fired column gives `a = 0` and `delta = n`, where both forms
return 1. Measured on `R2`'s own log over 50–200 s:

* **42.9 %** of stepped lateral faces at 2 mm are crater wall (taller column
  never fired ⇒ `col_dt = 0` ⇒ `a = 0`);
* the two forms **disagree on 71.7 % of those**.

The reason: `phi_T = 1` only because the **taller** column never fired. The
**shorter** column is active, so `phi_S < 1` and `delta = n + phi_S − 1 < n`.
D2z's `a → 0` branch returns `min(delta, 1) = 1`; the derived branch returns
`delta/n ≈ 0.93`. Small per face, but on ~31 % of all stepped faces.

⚠ **This means the derived form also acts on the crater wall**, which neither
packet anticipated. It is a reduction, not an amplification, and the rim test
(§3) shows it costs nothing.

## 5. (b) FAILS as pre-registered — two population errors, both identified

Live `delta_mean` on `R2w` = **0.5043**.

| radius | pre-registered offline | live-population offline | rel to live |
|---|---|---|---|
| r < 25 | 0.2145 | 0.2145 (identical) | **+135.1 %** |
| r < 30 | **3.1231** | **0.5401** | −83.9 % / **−6.62 %** |

Two distinct causes, neither a code defect:

1. **The radius.** `delta_mean` has **no 25 mm gate** — `step_diag_radius` gates
   only `WallStepDiagnostic`. `delta_mean` inherits `surface_patch.radius`,
   which is **30 mm** on the arbiter. Scoring at r<25 compares different areas.
2. **Never-fired columns.** The packet's instruction to include them at `phi = 1`
   on both sides is right about the *phi convention* and wrong about *inclusion*:
   a never-fired column at the patch edge sits beside one that has receded 40
   cells, so the pair carries **delta = 40** — the crater wall, not a phase lag.
   That drives the offline mean from 0.540 to **3.123**. The offline fallback
   fraction is **11.47 %** there, while live `delta_fallback_cols` is **exactly
   0** on every closure leg — independent confirmation that the live population
   contains no such column.

**With both corrected the agreement is −6.62 %, inside the 10 % bar.** Scored as
pre-registered (FAIL); the corrected reading is reported, not substituted.

## 6. (d) FAILS as written — recorded before launch, not after

| leg | flux-weighted `g` | `step_g_mean` (face-weighted) | rel |
|---|---|---|---|
| `R2` | 0.1967 | 0.1173 | 67.7 % |
| `R2w` | 0.3789 | 0.3130 | 21.1 % |
| `R2d` | 0.1960 | 0.1169 | 67.7 % |
| `R1d` | 0.1533 | 0.1217 | 31.6 % |

These are two *different* means of the same `g`, and `q_raw` changes sign where
the neighbour is hotter, so the flux-weighted mean is not even bracketed by the
range of `g`. The bar measures their **correlation**. The check with teeth — the
exact identity `wall_step_q_in == wall_step_q` and
`wall_step_q_eff_in == wall_step_q_eff` with the gate off — **passes bit for
bit** (`unit/side_face_step_height` k5).

## 7. Meier — (g), (h), (i)

At equal burner descent (150 mm), depths capped 50 mm above the burner.
**Rate definition: burner DESCENT slope over the 100–150 mm stretch.**

| leg | mean \|ΔØ\| | Ø@100 | Ø@125 | ROP | pit lead | vol | x_pit | x_desc |
|---|---|---|---|---|---|---|---|---|
| `SW` wetted, derived | **2.80** | 87.6 | 85.8 | **3.694** | 38.0 | 1.1 | 0.84 | 1.05 |
| `SWs` step, derived | 4.90 | 90.6 | 87.7 | 3.836 | 28.0 | 1.1 | 0.90 | 1.07 |
| `SW0` wetted, D2z `g` | 2.72 | 87.5 | 85.4 | 3.716 | 38.0 | 1.1 | 0.83 | 1.04 |
| `A` map, no side | 2.27 | 83.8 | 81.3 | 1.251 | 128.0 | 1.3 | 0.67 | 1.24 |
| `D` map, exact | 7.56 | 94.2 | 90.6 | 3.799 | 32.0 | 1.2 | 0.93 | 1.13 |

Ring table, 2 mm rings, **per ring, never averaged** [mm recession]:

| leg | 30 | 32 | 34 | 36 | 38 | 40 | 42 | 44 | 46 | 48 |
|---|---|---|---|---|---|---|---|---|---|---|
| `SW` | 256.3 | 259.4 | 258.4 | 258.1 | 258.3 | 247.8 | 62.1 | 1.7 | 0.0 | 0.0 |
| `SWs` | 246.9 | 250.6 | 254.6 | 258.2 | 262.6 | 267.7 | 187.5 | 21.2 | 0.0 | 0.0 |
| `SW0` | 256.2 | 259.0 | 258.4 | 258.3 | 258.1 | 250.8 | 72.8 | 1.2 | 0.0 | 0.0 |
| `A` | 124.9 | 119.2 | 112.7 | 104.3 | 90.1 | 69.1 | 37.9 | 11.7 | 3.4 | 0.2 |
| `D` | 199.2 | 199.3 | 199.1 | 199.1 | 199.0 | 199.7 | 199.2 | 199.5 | 142.2 | 28.0 |

### (h) `wetted` wins, and it is a real separation
2.80 against 4.90 mm, **2.10 mm apart** against a 0.25 mm scoring-depth floor.
D2z's UNDECIDED tie is broken on Meier in favour of the pre-registered
preference. The mechanism is visible in the ring table: `SWs` (`step`) keeps
cutting out to ring 44 where `SW` stops at 42.

### (i) the amplification error was worth almost nothing
`SW` against `SW0` — the only difference is the conductance form:
mean |ΔØ| **2.80 vs 2.72 mm**, ROP **3.694 vs 3.716 m/h**, volume identical.
`step_g_max` **1.2831 vs 2.7203** confirms the two forms really did differ (and
that `n` reaches `StepConductance`), and `f_multi = 0.73` confirms the branch was
exercised — yet the outcome moves by 0.08 mm and 0.6 %.

**⇒ The normalisation was an arithmetic error and is now correct, but it was not
the cause of anything.** Stated plainly so it is not credited with more.

### P6 REFUTED — `SW` DOES recover the wall band
Ø@100 = **87.6 mm** (88.2 at 180 mm descent) against Meier's 87. This is the
first Meier leg with a side-face path to reach the band without `D`'s overshoot
(94.2). ⚠ **No agreement claim tighter than ±5 mm** — Fig. 8.8 is a visible
width, i.e. a LOWER BOUND on diameter; band 85–95 ± 5 mm.

### P7 REFUTED, and it is the result that matters
`SW` descends at **3.694 m/h** against the 1.3–1.6 m/h band — **2.4x the top of
the band**, the same over-delivery `f = 0.2` produced (1.50 → 3.09). The ring
table shows why: `SW` removes ~258 mm across rings 30–40 where `A` removes
125 → 69.

**`SW0` reads 3.716 m/h**, so **the rate error is not the conductance form**. It
is the side-face heating path itself: `SW` draws **1022.6 W** of side power
(one-cell 127.6 W + multi-cell 895.0 W) against `A`'s zero. The packet's
"refuting evidence" check passes — `SW`'s one-cell draw is far below `D`'s
1234.7 W, so `delta` is not inverted — the problem is the *amount*, not the sign.

## 8. ⚠ A reference value that does not reproduce

Per the standing rule, reported rather than replaced.

| quantity | packet | reproduced |
|---|---|---|
| `A` mean \|ΔØ\| at 150 / 180 mm | 2.27 / 2.51 | **2.27 / 2.51** ✓ |
| `A` Ø@100 | 84.4 | **84.4 at 180 mm** ✓ (83.8 at 150) |
| `A` ROP over 100–150 mm | 1.251 | **1.251** ✓ |
| D2x side power `C`/`D`/`E` | 223.2 / 1234.7 / 453.8 | **223.2 / 1234.7 / 453.8** ✓ |
| **`D` Ø@100** | **97.1** | **94.2 at 150 mm, 95.1 at 180 mm** ✗ |

Everything reproduces except `D`'s Ø@100, which is **2.0–2.9 mm low at either
scoring depth**. `A`'s 84.4 is the **180 mm** reading, so the packet's P6 line
mixes scoring depths. **No scored outcome turns on this**: P6 is refuted at both
depths (87.6 / 88.2), and `D` is a reference leg only.

## 9. What this leaves open

* **The rate.** `SW` is 2.4x the band and `SW0` rules out the conductance form.
  The side-face heating path delivers ~1 kW that `A` does not, and the multi-cell
  faces carry **87.5 %** of it. That is the next question, and it is about **how
  much heat a lateral face should receive**, not about `delta`.
* **`step` vs `wetted` is now settled on shape** (2.10 mm), but both laws
  over-deliver rate almost equally (3.694 vs 3.836), so the law is not the lever.
* **The arbiter cannot test any of this.** `wall_step_f_multi` is **0.0000 on
  every arbiter leg at both meshes**, including `R1d`. Only Meier has the
  population.
