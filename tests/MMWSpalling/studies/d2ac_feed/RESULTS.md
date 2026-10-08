# D2ac — Meier at Meier's kinematics: the floor-flux excess

**The deliverable is a number: the model's floor flux exceeds the one that would
give Meier's rate at Meier's standoff by 1.98 – 2.23×.** Pre-registered
2.15 – 2.60×; measured inside that, on four legs, with the rate as an *input*.

Driven at Meier's own feed, the model does not hold Meier's standoff: it opens to
170 – 181 mm by 500 s and is **still growing**, against Fig. 8.8's ≈100 mm
ceiling. Each leg's own decay-corrected rate implies an equilibrium of
**198 – 223 mm**, roughly twice the ceiling.

## 0. Provenance

* **`PREFLIGHT.md`** `6b4fc3cd2173a038446c4991b0256d9c1c52ec481e180eda91148d459815f593`,
  **`preflight.py`** `3be11f33c589ae4f4f30912cd5314c6471de298333b1c0ba941801b066d3a30d`
  — frozen chmod 444 before `CRITERION.md`, zero outputs on disk.
* **`CRITERION.md`** `4907c247e46a64c1db86f242358e3a6cf06297af62f2c3f9bb353432a8c99e2a`,
  chmod 444, zero outputs on disk when hashed.
* **`PREFLIGHT_ADDENDUM.md`** — records an **error in the frozen pre-flight**,
  found by running. The frozen files are **not modified**.
* **Binary** `ea833f2f03855145a93eb98508e04a92cbc1aff80bfda40b1cc32bf96776c471`
  — **no source change in this packet**; D2ab's post-guard build, unchanged.
* Nothing committed. Frozen studies imported by path, never re-run.

## 1. Two errors found before any result was read

### (i) The packet's leg spec was incomplete — `prescribed` refuses ten keys

The packet says the shared block "replaces `nozzle_descent=feet`". It is not
enough: `A`'s argv carries a feet-only block that the code refuses. Measured with
a parse probe, not reasoned:

| attempt | result |
|---|---|
| swap `nozzle_descent` only | rc=6, `:7593` *"`jet_free_surface = 1` requires `nozzle_descent = feet`"* |
| + drop `jet_fs_*` | rc=6 (`foot_body_clearance`, `:7596`) |
| + drop the whole block | **rc=0, unused-key list EMPTY** |

Ten keys dropped. **Forced by the code, not chosen** — but it means **`FM`
differs from `A` by more than the feed**, stated in every comparison below.
Bounded: free-surface dilution barely fired on the frozen legs (`A`
`jet_fs_cols` mean 0.019, `m_ratio` max 1.040; `SWc` exactly 0 and 1.000).

⚠ Consequence: `foot_r_inner/outer` fall to 0.0, so the D2aa region split bins
**every** charged face as `wall`. Those columns are **meaningless here** and are
used only as a reached-a-face witness, never as a region split.

### (ii) The collision guard is the OPPOSITE of inert — both the packet and my own pre-flight were wrong

`ACTIVE_STEP.md` item 0(b) says the guard "cannot fire at any feed"; `PREFLIGHT.md`
(b) said it "defaults to 0.0 and cannot fire". **Both wrong.** At `:7676`:

```cpp
pp.query_default("nozzle_collision_radius", v.nozzle_collision_radius,
                 v.feet ? v.foot_r_inner : v.radius);
```

Under **`prescribed`** the default is `surface_patch.radius` = **0.2 m**, so
every column in the domain is a collision column. **All four legs aborted at
t = 120 s**, at r ≈ 59 mm — unrecessed rock beside the hole, which the code
itself calls normal (`:7670`). The packet's own "Refuting evidence" had flagged
D2a2 hitting this guard; both the planner and I read it as a drain-era artefact.

Fixed with `nozzle_collision_radius = 0.028` — the code's **own feet-mode
default** (`foot_r_inner`, the Ø56 bore radius). `rc2` appears at exactly two
lines, its definition and the abort, so it **moves no trajectory**. `S0` ran
before the fix and is unaffected; no rebuild, so no identity re-proof is owed.

## 2. Scorecard

| | result | |
|---|---|---|
| **(a)** floor-flux excess | **the product** | 1.98–2.23× vs 2.15–2.60× pre-registered |
| **(b)** pad contact, offline | **4/4 contact** | transient: 4.1 s, 5–17 % of run |
| **(c)** Ø 150–300 mm vs 85±5 | **2/4** | FL, FH pass; FM, FM13 fail high |
| **(d)** volume ≤ 3.42 L/0.5 m | **2/4** | FL, FH pass; FM, FM13 exceed |
| **(e)** standoff ≤ 100 mm | **0/4 — PRE-REGISTERED FAILURE** | 170–181 mm |
| **(f)** power identity | **5/5 PASS** | ≤1.2e-2 W, `ledger_err` ≤8.8e-15 |
| **(g)** identity | **PASS** | `S0` byte-identical; tree 2639/2639 |

## 3. (a) THE DELIVERABLE — the floor-flux excess

| leg | feed | s @500 s | floor R | excess at s | **R @100 mm** | **excess @100 mm** |
|---|---|---|---|---|---|---|
| `FM` | 1.50 | 169.7 | 2.049 | 1.37× | 2.969 | **1.98×** |
| `FL` | 1.50 | 177.7 | 2.096 | 1.40× | 3.176 | **2.12×** |
| `FH` | 1.50 | 173.7 | 2.048 | 1.37× | 3.030 | **2.02×** |
| `FM13` | 1.30 | 181.5 | 1.905 | 1.47× | 2.894 | **2.23×** |

Standoffs in mm, rates in m/h. `R @100 mm` de-scales each leg's measured floor
rate to a common 100 mm standoff through **the map's own decay term**
`pow(0.09/max(0.09,s),1.0)`, so the number is on item 0's basis.

**Measured 1.98 – 2.23× against a pre-registered 2.15 – 2.60×.** The measurement
sits at the low end of, and overlaps, the prediction.

⚠ **The raw excess at the realised standoff is only 1.37 – 1.47×**, because the
hole has already opened to ~170 mm and the map's 1/s decay has throttled the
flux. That is the model *partly self-correcting by retreating* — which is
precisely why the comparison must be made at Meier's standoff, not at the
model's.

## 4. (e) The standoff, and why P1 needs care

| leg | 50 s | 100 s | 200 s | 300 s | 400 s | 500 s |
|---|---|---|---|---|---|---|
| `FM` | 87.2 | 108.3 | 132.7 | 149.0 | 161.3 | **169.7** |
| `FL` | 87.2 | 110.3 | 138.7 | 157.0 | 169.3 | **177.7** |
| `FH` | 87.2 | 110.3 | 136.7 | 153.0 | 165.3 | **173.7** |
| `FM13` | 88.0 | 111.9 | 139.8 | 157.7 | 169.6 | **181.5** |

**Every leg fails the ≤100 mm ceiling — the pre-registered outcome.** It crosses
100 mm before t = 100 s and never returns.

⚠ **The run does not reach equilibrium.** Growth is still 8–12 mm per 100 s at
t = 500, and the measured floor rate (2.05 m/h) is still 1.37× the feed. So
170–181 mm is a **lower bound**. Each leg's own implied `R90` gives the
equilibrium it is heading for:

| leg | feed | R measured | R90 implied | **s_eq implied** | pre-registered |
|---|---|---|---|---|---|
| `FM` | 1.50 | 2.049 | 3.299 | **198 mm** | 215–222 |
| `FL` | 1.50 | 2.096 | 3.528 | **212 mm** | 215–222 |
| `FH` | 1.50 | 2.048 | 3.366 | **202 mm** | 215–222 |
| `FM13` | 1.30 | 1.905 | 3.216 | **223 mm** | 248–256 |

So **P1 is held on the quantity it predicts** (equilibrium, 198–223 inside
180–260) and refuted on the direct 500 s reading for three legs. Scored as
measured; the equilibrium reading is reported alongside, not substituted.

## 5. (b) Pad contact — held, but transient

| leg | first contact | contact fraction | worst interference | clearance @500 s |
|---|---|---|---|---|
| `FM` | 4.1 s | 8.3 % | +2.6 mm | −93.9 mm |
| `FL` | 4.1 s | 17.4 % | +5.8 mm | −77.9 mm |
| `FH` | 4.1 s | 11.6 % | +2.9 mm | −95.9 mm |
| `FM13` | 4.1 s | 5.0 % | +1.5 mm | −103.6 mm |

**P3 HELD** — the pads do touch, within the first seconds and intermittently for
5–17 % of the run, by up to 5.8 mm. But the ring then **clears permanently**: by
500 s the feet are 78–104 mm above the rock. The burner ends up hanging, which is
the same statement as (e) — the floor outruns it.

⚠ **Contact at t = 0 is by construction** (`nozzle_z0 = lz + foot_standoff` puts
the feet exactly on a flat surface) and is excluded. An earlier revision of this
scorer reported "never" for all four legs through a falsy-zero bug — `if first`
rather than `if first is not None` — which would have inverted P3.

## 6. (c) and (d) — shape and volume

| leg | calibrated Ø 150–300 mm | tip-excluded | vs 85±5 | volume L/0.5 m | vs 3.42 |
|---|---|---|---|---|---|
| `FM` | 93.72 | 95.44 | FAIL high | 3.540 | FAIL |
| `FL` | 80.24 | 80.57 | PASS | 2.710 | PASS |
| `FH` | 87.14 | 88.48 | PASS | 3.151 | PASS |
| `FM13` | 90.14 | 96.78 | FAIL high | 3.632 | FAIL |

⚠ The hole only reaches 312–336 mm in 500 s, so the top of the 150–300 mm window
can sit inside the **rounded tip**, where Ø collapses (`FM13` reads 58.1 mm at
300 mm). The tip-excluded column drops bins within 40 mm of the floor and is
**reported alongside** the pre-registered score, never instead of it.

* **P4 HELD** — `FM` exceeds 85 mm (93.7; exposure scaling predicted ≈99).
* **P5**: `FL` 80.2 and `FH` 87.1 both land in the open 80–95 mm band. ⚠ **This
  is not a derived agreement.** P5 was open by construction and 85 is Meier's own
  number; `FH` landing 2 mm from it is not evidence.
* **P6 REFUTED on `FM`/`FM13`** — and this is a result, not a miss.
  **The model removes MORE rock per metre than Meier's water-filled 3.42 L**,
  even though Meier's figure *includes* a mouth and a pause underream the model
  never makes. The excess is the flux excess showing up in volume.

Only `FL` — the annulus wall throttled to 38.4 W/m²K — is comfortably under on
both shape and volume, and it is the narrowest hole.

## 7. Identity and witnesses

* **`S0` vs `d2z1_delete/RB_D`**: `_removal_events.csv` and `_jet_profile.csv`
  byte-identical, thermo **header delta EMPTY**, 0 columns differ.
* **Tree gate 2639/2639** vs `witness/d2z1.md5` — see `GATE.txt`.
* **(f) power identity** ≤ 1.2e-2 W on all five legs, `ledger_err` ≤ 8.8e-15.
* **Feed applied**: exact by construction — `nozzle_z` is registered only under
  the feet rule (`:629`), so under `prescribed` z_n is analytic and the standoff
  is read from `jet_s_c`.
* **Pit vs domain**: 312–336 mm of 500 mm. The floor never reached the boundary,
  so no leg is truncated.

## 8. What this establishes

**The number.** At Meier's standoff the model's floor flux is **≈2× too high**
(1.98–2.23×). Driven at Meier's feed the model cannot hold Meier's geometry: it
retreats to a standoff roughly twice the photograph's ceiling, and even there it
removes more rock per metre than the water-filling measurement allows.

**The follow-on is unconditional**, and the candidates must now be sized against
a factor of ≈2:

* nozzle gas temperature over Meier's 1450–1850 K — ≈30 % of the flux;
* gas-to-body heat loss to the water-cooled tool, 160 L/h = 186 W/K; 10–30 K
  gives 1.9–5.6 kW — real, unmeasured, therefore a **bound**;
* under-expanded supersonic impingement from the Laval nozzle (6 mm throat);
* spall temperature within Rauenzahn's measured band — the upper end raises the
  specific energy by up to ≈1.7×;
* gas radiation — absent, and worth ≈10 % **in the wrong direction**.

No single candidate covers 2×; the temperature and spall-temperature terms
together plausibly do.

## 9. Named debt

* **The D2l factor-of-four audit is NOT owed.** The packet lists it as owed by
  five packets; **D2ab discharged it** (reproduced: 3204.2 W quarter × 4 =
  12.82 kW, matching `d2l_scan/RESULTS.md:190`). Pointer in `PREFLIGHT.md` (c).
* **Lateral recession** — still the critical path (blocks the pad ring, D2j's
  continuous front, and D2ab's narrower-than-the-burner hole).
* `beam_void_closure` **T3** — the only red in the tree, unchanged.
* `A` at 1 mm — user-parked.
* ⚠ **The guard cache is still armed for a false pass** — it is keyed on the
  binary alone and `exit 0`s with no `<label>.md5` on a hit, so an input-only
  change would pass silently. Run with `--force` when inputs change, as here.
* ⚠ **`nozzle_collision_radius`'s prescribed default (`radius`, 0.2 m) is a
  trap** for any future prescribed leg, and nothing in the source says so at the
  call site.
