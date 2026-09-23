# D2t pre-flight — written BEFORE any leg launched

Two separate obligations live here:

* **Item 0**, the standing rule: name the call site of every new or changed
  behaviour and demonstrate, from the literal key set of each leg, that it is
  reached. Three packets (D2r-0 / 0b / 0c) and D2s died on this.
* **The reference reproductions.** Every number this packet states was
  recomputed here with my own estimator before it was trusted. Disagreements
  are reported, never silently replaced.

The compiled-taper thermal trace (Goal item 2) is a **separate pre-flight with
its own hash**: `TAPER_TRACE.md` / `TAPER_TRACE.sha256`.

---

## 0. Reachability ledger — every leg, every call site

No `src/` change in this packet, so "new behaviour" means: the taper inside
`h_expr`, and the wider `surface_patch.radius`. Both are key-only.

The arbiter never sets `surface_patch.jet_closure`, so it defaults to `none`
and `jet_enthalpy = false` (`:6123-6126`). **That is what makes the
`else` branch at `:1785-1786` the reached one and `JetEnthalpyMarch` (`:1797`)
dead** — the same fact that killed D2r-0c, working in this packet's favour.

| # | behaviour | call site | reached because | legs |
|---|---|---|---|---|
| 1 | radial `h` | `h = surface_patch.h_f(x, y, r, s, time)`, **`:1785`** | every leg sets `surface_patch.h_expr`; `jet_enthalpy` false so `:1773`'s jet branch is not taken | all 6 |
| 2 | `T_gas` | `Tg = surface_patch.Tg_f(...)`, **`:1786`** | every leg sets `T_flame_expr="1600.0"` | all 6 |
| 3 | the `h > 0` abort | **`:1787-1792`** | evaluated for every in-patch live column; taper floor is `1e-3` (below plane) / `2.03e-4` (above), both `> 0` | all 6 |
| 4 | top-face heating across the taper | `in_patch`, `r2 <= sp_r2`, **`:1579`** | `surface_patch.radius = 0.050` > `r_taper = 0.034` | all 6 |
| 5 | side-face term across the taper | `col_in_patch`, same cutoff, **`:1626`** | same radius; every leg sets `side_face_flux=1`, `side_face_factor=0.2` | all 6 |
| 6 | above-plane side factor | selector on `fs[col_c] > 0.0`, **`:1651-1656`** | `PA_*` set `side_face_factor_above=1.0`; `P_*` leave it unset (sentinel `-1` ⇒ uniform `0.2`) | `PA_2mm`, `PA_1mm` |
| 7 | the `s` switch, above branch | `if(s>0.0,700.0,142.0)` inside `h_expr`, evaluated at `:1785` | `flame_s[col] = z_n - z_face` (`:1759`) is set for **every** column, before the `r2 > r2max` cut | `TR_above` by construction; `P/PA` once the band falls behind |
| 8 | ledger twin of the side term | **`:2762+`** | `energy_ledger.enabled = 1` in `input_hotdisk`; mirrors #5/#6 or `ledger_err` breaks | all 6 |

### The empirical demonstration (literal keys, from `run.py --list` / `job()`)

Every leg, verbatim:

```
surface_patch.radius=0.05
surface_patch.h_expr="(1.4285714285714286e-06+0.9999985714285714*0.5*(1.0+cos(
   3.141592653589793*(min(0.034,max(0.018,r))-0.018)/0.016)))*if(s>0.0,700.0,142.0)"
surface_patch.T_flame_expr="1600.0"
surface_patch.nozzle_descent=prescribed
surface_patch.nozzle_z0=0.17  nozzle_feed=0.00043056  nozzle_collision_radius=0.002
surface_patch.side_face_flux=1
surface_patch.side_face_factor=0.2
```

plus `surface_patch.side_face_factor_above=1.0` on `PA_2mm` / `PA_1mm` only, and
for the two trace legs `stop_time=0.16`, `amr.plot_int=10`, with `TR_above`
additionally at `nozzle_z0=0.14`, `nozzle_collision_radius=0.0`.

`jet_closure` appears in no leg. **No leg fails to reach any site, so the rule
does not stop this packet.**

### Trace overrides replace, they do not duplicate

`TR_above` needs `nozzle_z0` and `nozzle_collision_radius` different from the
frozen base. The harness **replaces the key in place** and raises if the base
ever carried it twice, so nothing here depends on how ParmParse resolves a
duplicated name. Verified: each of those keys appears exactly once in each
trace leg's argv.

---

## 1. Reference reproductions — my own estimator, from the frozen event logs

Metric: per-column recession summed from `_removal_events.csv` column 9 over the
block, divided by the block length, disk-mean over the stated radius mask; gap =
`1mm/2mm - 1`. Blocks 50-100 / 100-150 / 150-200 s.

| leg | mask | packet | measured | |
|---|---|---|---|---|
| `B` (= D2l `A_f02`) | `r < 30` | −6.4 / −8.0 / −14.0 | **−6.4 / −8.0 / −14.0** | ✓ exact |
| `B` | `r < 25` core | −0.2 / +3.0 / −1.8 | **−0.2 / +3.0 / −1.8** | ✓ exact |
| `B` | 25–30 rim | −21.4 / −33.3 / −44.1 | **−21.4 / −33.3 / −44.1** | ✓ exact |
| `T` (D2r-1) | `r < 30` | −2.3 / +0.7 / −2.0 | **−2.3 / +0.7 / −2.0** | ✓ exact |
| `T` | `r < 25` core | −0.0 / +3.3 / −1.5 | **−0.0 / +3.3 / −1.5** | ✓ exact |
| `T` | 25–30 rim | −7.8 / −5.2 / −3.7 | **−7.8 / −5.2 / −3.7** | ✓ exact |

**Criterion (d)'s comparator**, recomputed on the `r < 18 mm` mask, block
150–200 s, as the packet demands: `T_1mm` **1.6349**, `T_2mm` **1.6667** m/h —
both exact to the quoted 4 figures. Against the feed `4.3056e-4 m/s` =
**1.5500 m/h**, the margins are **+5.48 %** and **+7.53 %** (packet: +5.5 / +7.5).
**The manufactured lagging band is therefore real and predicted correctly.**

`z_nozzle(250 s) = 0.170 − 4.3056e-4 × 250 =` **62.36 mm** (packet: 62.4) ✓.

`T_2mm`: **1** in-patch above-plane column at r = 29.83 mm; `T_1mm`: **13**,
r = 29.71–29.91 mm, **11 never spalled** — all three exact ✓.

P1's population estimate also reproduces: columns with 19 ≤ r < 34 mm number
**156 at 2 mm / 624 at 1 mm** by direct grid count (packet: 159 / 629; the
difference is edge rounding). **P1 is sound.**

### Three disagreements, reported not replaced

1. **"111 in-patch columns are above-plane in `B_1mm`, spanning r = 27.3–29.9 mm."**
   I measure **110**, spanning **27.50–29.91 mm**. One column and 0.2 mm on the
   inner edge — a `s <= 0` vs `s < 0` boundary convention. Changes nothing.
2. **"the front sits 18–21 mm BELOW [the plane]."** Measured mean `z_face` over
   `r < 18 mm` at t = 250 s sits **21.9 mm** (2 mm) / **23.6 mm** (1 mm) below.
   Slightly deeper than stated.
3. **"the front reaches ~104 mm at 250 s."** Measured recession is
   **109.5–111.2 mm** mean, **110–112 mm** at the deepest column. **This is the
   one that matters**, because the packet uses it for the domain watch. Floor
   clearance is therefore **38–40 mm**, not ~46 — still far above the 10 mm
   assert, so the campaign is safe as specified. Asserted again after the run.

---

## 2. Design readings recorded before launch

**The floor `1.0e-3` is read as a floor on `h`, not on the multiplier.**
`:1787` aborts unless `h` is finite and `> 0`, so the taper cannot reach zero;
the packet fixes the floor at `1.0e-3` and justifies it by citing `C1`, whose
`1.0e-3` is a **value of h**. Taken instead as a multiplier floor, `h` would
land at `0.7 W/m²K`; against `T_gas = 1600 K` that is ~915 W/m² into ambient
rock, which balances the `ε = 0.8` radiative and `h_c = 10` losses at
**≈ 350 K**, warming the whole unheated far field by ~57 K. That is not
"emulating `h = 0`". So `F0 = 1e-3/700 = 1.4286e-6`, giving `h = 1.000e-3`
below the plane and `2.03e-4` above — both positive, both inert.

**`surface_patch.radius`: the packet says 0.050 in Goal item 1 and 0.040 in the
item-5 ledger.** Both satisfy the stated constraint (`>= r_taper + 2` coarse
cells = 0.038). **Taking 0.050**, which is the Goal section's own recommended
set and carries its own justification ("0.050 costs nothing"); the ledger row
is the stale one. It puts the hard cutoff **8 coarse cells** into the floor.

**Taper form: raised cosine**, C1 at both ends as the packet prefers over a
ramp. It is also the exact form already proven against this parser by the
frozen `d2j0_hotdisk` tapered legs (`H_EXPR`), including the numeric-first
`min`/`max` and a literal — never a `max()` — in the denominator.

**Geometry check.** `r_flat = 18 mm`, `r_taper = 34 mm`, `W = 16 mm` = **8
coarse cells / 16 fine cells**, a fixed *physical* width at both meshes.
`radius = 50 mm` fits inside the 80 mm quarter domain with 30 mm to spare.
`h` falls 5 % below peak at **r = 20.3 mm**, so with the measured +5.5 % core
margin the lagging band starts there — consistent with the packet's "19–34 mm".

**What the trace can and cannot see.** At t ≤ 0.16 s every column is at
`z_face = 0.15`, so `TR_below` (`nozzle_z0 = 0.17`) reads the `s > 0` branch,
`h = 700·taper`. `TR_above` (`nozzle_z0 = 0.14`) puts every column above the
plane and reads `h = 142·taper` from the **same compiled expression**. Their
ratio must be **142/700 = 0.202857** at every radius: that is the check that
the taper multiplies the whole `if()` and is not mis-associated with one branch.
