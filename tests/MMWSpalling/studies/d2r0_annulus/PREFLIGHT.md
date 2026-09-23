# D2r-0 pre-flight 0 — arithmetic and existing output only

**Result: items 1, 2, 3 PASS. Item 4 STOPS. No code was written.**

Per the packet: *"Any STOP means close the packet and report, as D2p-0 and D2p-1
were closed."* Reproduce with `python3 preflight.py` (reads frozen study output
by path; runs nothing).

## Scorecard

| item | test | result |
|---|---|---|
| 1 | `jet_r_reach` collapses, `jet_T_exhaust` high | **PASS** — 168 → 45.8 mm, 1689 K |
| 2 | honest lower bound on `h_wall` ≥ ~20 W/m²K | **PASS** — 81 W/m²K |
| 3 | stream must not deplete before the mouth | **PASS** — 8.5–33 % drawn |
| 4 | predicted drain reduction ≥ 2× at the band's **lower** end | **STOP — 1.46×** |
| 5 | predicted wall temperature | 545 K at `A_low`, **misses (c)'s 600 K bar** |

## 1. The exhaust budget — reproduced exactly

| t [s] | `jet_P_face` | `jet_P_exhaust` | `jet_T_exhaust` | `jet_r_reach` |
|---|---|---|---|---|
| 95 | 2756.7 W | 8892.4 W | 672.2 K | 168.3 mm |
| 285 | 1359.2 W | 14151.2 W | 1612.5 K | 53.8 mm |
| 571 | **997.0 W** | **17014.1 W** | **1689.1 K** | **45.8 mm** |

The packet's table is reproduced to the digit. From the march's own invariant
`jet_P_cap = m·cp·(T_stag − T_ent)`, the stream is **mcp = 12.188 W/K**
(m = 9.751 g/s, m_ratio 2.58), and `mcp·(T_exh − T_ent) = 17014.0 W` against the
logged 17014.1 — so `jet_P_exhaust` is referenced to `T_ent` and the 17 kW is
genuinely available enthalpy, not an absolute-zero artefact. **NO STOP.**

## 2. The wall `h` — both derivations, and where the band's floor really sits

The floor anchor first, because the packet quotes a **simplified stand-in**:

| | `h` at r = 45.8 mm |
|---|---|
| the `h_expr` `LI0_2mm` **actually ran** (D2a2 wall-jet form) | **605.4 W/m²K** |
| the packet's stated `1500·0.01875/max(0.01875,r)` | 614.1 |
| the packet's stated value | 611 |

Agreement within 1.5 %, so the stand-in is faithful **at this radius** and the
anchor holds.

- **(i) area expansion.** 2 mm wall jet → 12.3 mm annulus = **6.15×**; `u/6.15`,
  `h ~ u^0.8` ⇒ **142 W/m²K**.
- **(ii) Nusselt.** Quarter flow area 7.66e-4 m², `D_h` = 24.6 mm,
  `mu` = 5.62e-5 (Sutherland), `k_gas` = `mu·cp/Pr` = 0.100 W/mK — taken from the
  **run's own cp**, not a textbook conductivity. **Re = 5568** ⇒ Dittus-Boelter
  Nu 19.8 ⇒ **81 W/m²K**, ×1.24 for entry length (L/D_h = 7.5) ⇒ **100 W/m²K**.

Both land inside the band. **But the band's lower end does not.** A
fully-developed *laminar* annulus (Nu = 5.7) gives 23 W/m²K, which clears the
20 W/m²K threshold — **yet Re = 5568 is not a laminar flow**, so 23 is an
arithmetic floor, not a regime this annulus is in. **The lowest value either
derivation actually supports is 81 W/m²K.** The packet's 50 sits below both.
**NO STOP**, but this is the root of the item-4 STOP.

## 3. Affordability — the stream is not the constraint, and the bound does not bind

Wall area measured from the frozen removal-event log, not assumed: the final
surface is a **cone 360 mm deep** with a mouth near r = 78 mm, giving an exposed
vertical (staircase) face area of **0.0310 m² in the quarter**. Counter-current
march `T_out = T_w + (T_gas − T_w)·exp(−hA/mcp)`:

| h | draw | % of 17.0 kW | gas cooling | `T_gas` at the mouth |
|---|---|---|---|---|
| 50 | 1442 W | 8.5 % | 118 K | 1571 K |
| 150 | 3829 W | 22.5 % | 314 K | 1375 K |
| 250 | 5679 W | 33.4 % | 465 K | 1223 K |

**NO STOP** — the stream does not deplete.

**⚠ But this falsifies prediction P3 and changes what criterion (g) is testing.**
The packet's item 3 says to prove the item-1.1 restriction inert by showing
`T_gas` at the last wall bin is *"within a few K of `T_ent`"*. It will be
**1223–1571 K**. The enthalpy bound is real in total (≤17 kW) but **it does not
bind inside the hole**, so the "recessed below the original surface" restriction
is **doing physics, not bounding cost**: the excluded flat rock outside the mouth
would receive a large flux. `A_wide` should therefore differ **materially** from
`A_mid`, and criterion (g)'s "must be justified" branch is the live one.

A justification does exist and should be stated rather than discovered: gas past
the mouth is no longer confined in an annulus, and Meier's sealed water-cooled
wellhead ducts it away. That is a physical modelling choice, not a cost bound.

**Planner arithmetic that does not reproduce:** the packet estimates ~1.15 kW and
~100 K of cooling. On the measured wall area the draw is **1.4–5.7 kW** and the
cooling **118–465 K**. The verdict (no depletion) is unchanged.

## 4 + 5. The wall temperature and the drain — **STOP**

My first two attempts at a closed-form wall temperature were wrong and are
recorded as such: a semi-infinite Robin without losses returns 915–1270 K (above
the firing point), and a lumped steady balance with a conduction source fails its
own sanity check by 372 K (predicts 790 K where D2o-0b measured 418 K).

**The project already measured this exact quantity.** D2q-2c ran a vertical wall
under a side-face Robin at **h = 80 W/m²K against `T_gas` = 1600 K** (d2j0b's
`H_IN`/`T_IN`) in this code, with these losses, at 2 mm, and the wall equilibrated
**steady at 630–650 K**. Calibrating the quasi-steady rock sink on that
measurement — `h(T_g − T_s) = loss(T_s) + S(T_s − T_amb)` ⇒ **S = 190 W/m²K**
(181–200 across the 630–650 K anchor):

| h | `T_gas` | `T_wall` | step vs 820 K | drain | % `q_pin` | reduction |
|---|---|---|---|---|---|---|
| **50** | 1630 | **545 K** | 275 K | **206 kW/m²** | 40 % | **1.46×** |
| 80 | 1598 | 639 | 181 | 136 | 26 % | 2.23× |
| 150 | 1532 | 777 | 43 | 32 | 6 % | 9.39× |
| 250 | 1456 | 879 | −59 | −44 | −9 % | flux reverses |

The h = 80 row returns the anchor, so the model is self-consistent where it was
fitted.

- **STOP (item 4): reduction at the band's lower end = 1.46×, below the 2× bar.**
- **Criterion (c) also misses at `A_low`: 545 K against the > 600 K bar.**
- The 2× bar is met from **h ≥ 73 W/m²K**.
- The wall reaches `T_fire` = 822.2 K only at **h = 186 W/m²K**. So across the
  whole band the wall is **warmed, not spalled** — which is what criterion (c)
  intends, and it means `per_cell_removal` staying off costs nothing.

### Why this is a criterion problem, not a physics problem

**The mechanism works at both honest derivations of `h`** (81 → 2.3×, 100 → 3.4×,
142 → 8×). It fails only at 50, a value **below the lowest either derivation
supports**. The packet chose to score at 50 on the sound principle that the drain
falls monotonically with `h`; the pre-flight's finding is that the chosen floor
is not defensible as the conservative end of *this* band.

### ⚠ A second, larger problem with criterion (b) as written

**The ≥ 2× bar is already met by code that exists today.** D2o-0b's 402 K step is
an **h = 0** wall. D2q-2c *measured* the wall at 630–650 K with the existing D2k
side-face path at its default `side_face_h_max` = 80 — a 181 K step, i.e.
**2.23× — a pass on criterion (b) with the annulus march switched off entirely.**

Criterion (b) is anchored to `C0` (`jet_annulus = 0`), and the packet does not
say whether `C0` carries `side_face_flux`. If it does not, **(b) measures the D2k
side-face path, not the annulus march**, and a pass would be unattributable in
exactly the way D2q-2c's §2 attribution existed to prevent.

## 6. Design notes the planner should have before rewriting

- **A single `h_ann` at the burner z-plane is not valid for the whole wall.** The
  hole flares from r = 45.8 mm at the nozzle plane to ~52 mm at mid-wall and
  ~78 mm at the mouth. Flow area grows 1.7× to 5×, so `u` and `h ~ u^0.8` fall by
  **1.5× over the confined part and up to 3.6× at the mouth**. Item 1's "computed,
  not set" should be computed **per bin**, which the ascending-z binning already
  makes natural.
- **Constants that do not match the run.** The packet states `jet_mdot = 1e-3
  (quarter)`; `input_feet` and the `LI0_2mm` command line both carry
  **3.7826e-3**, and the march's entrained `m` is 9.751e-3. Nothing above depends
  on it — I used the measured `jet_P_exhaust` throughout — but the stated constant
  is wrong.
- **The reference geometry is the D2i runaway**, a 360 mm cone, not a 95 mm
  cylinder. Every area above is measured on it and inherits that.
