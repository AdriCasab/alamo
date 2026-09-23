# D2r-1 results — heating the SIDE faces of above-plane columns

## Verdict: **(a) PASSES. The mesh gap closes, and it is the SIDES that do it.**

| block | frozen baseline | test `T` | tops alone `S` |
|---|---|---|---|
| 50–100 | −6.4 % | **−2.3 %** | −6.4 % |
| 100–150 | −8.0 % | **+0.7 %** | −7.9 % |
| 150–200 | −14.0 % | **−2.0 %** | −13.9 % |

All three test gaps are inside 5 %. **`S` — which heats the above-plane tops at
the same h = 142 but forces the above-plane side factor to zero — reproduces the
baseline to 0.08 %, at BOTH meshes.** So the tops are inert and the side faces
carry the entire effect. This is the packet's title question, answered.

## 0. Provenance

| item | value |
|---|---|
| `CRITERION.md` sha256 | `33b27e479f94844c57146a4709504fdd6d58342caac33bcf113641beef5cbaef`, read-only, hashed **14:23:55**, before either 1 mm leg finished (`shasum -c` verified at scoring) |
| identity, key unset | **2639/2639 byte-identical** to `d2q2c.md5`, via `guarded_post_all.sh` (binary verified unmoved mid-sweep) |
| regressions | **8/8 PASS** |
| `ledger_err` | ≤ **1.0e-14** on all six legs |
| baseline validity | `B_2mm`/`B_1mm` reproduce frozen `A_f02`/`A_f02_1mm` **exactly** (0.000e+00) |

**Frozen, imported by path, never re-run:** `d2l_scan/output/A_f02*`,
`d2j0b_plane/` (keys, metric, `analyze.py`).

## 1. Goal item 1 — the stop condition that did NOT fire

**No input in the 2639-file identity set turns `side_face_flux` on**, and
D2q-2b/2c both modified that path — so tree byte-identity said *nothing* about
the frozen `A_f02` baseline. Re-running its keys on the current binary
reproduces it **byte for byte** (`_removal_events.csv`, `thermo.dat`, every
compared `Level_0/Cell_D_00000`). The 6 / 8 / 14 % baseline is valid.

## 2. Attribution — tops vs sides, measured at both meshes

`S` isolates the tops (`S − B`); `T − S` isolates the sides.

| | 2 mm tops | 2 mm SIDES | 1 mm tops | 1 mm SIDES |
|---|---|---|---|---|
| 50–100 | +0.00 % | +0.09 % | +0.02 % | **+4.50 %** |
| 100–150 | −0.02 % | +0.07 % | +0.07 % | **+9.49 %** |
| 150–200 | +0.00 % | +0.03 % | +0.08 % | **+13.86 %** |

**The tops contribute nothing at either mesh**, independently confirming
D2j-0b's `C2` (which held above-plane tops warm and was inert). The sides
contribute everything, and only at the fine mesh.

## 3. Why it works — the remedy self-scales with the artefact

| | 2 mm | 1 mm |
|---|---|---|
| above-plane exposed side faces (`side_cols_above`) | 103 | **2441** |
| power into them (`side_P_face_above`) | 43.6 W | 246 W |
| rate change vs baseline | +0.1 % | +4.5 / +9.6 / +13.9 % |

Refinement multiplies the staircase faces ~24×. The drain is the only 1/dx term
in the model, and this remedy acts on exactly the faces the drain flows into — so
it is negligible where the drain is negligible and large where the drain is
large. **The coarse mesh barely moves while the fine mesh rises to meet it**,
which is the signature of removing a mesh artefact rather than of adding heat
(that would have lifted both meshes together).

## 4. The drain, measured (owed since D2q-2c)

`wall_step_q`, the mean `k·dT/dx` across the wall step:

| | 2 mm | 1 mm | ratio |
|---|---|---|---|
| `S` (sides off) | 116.6 kW/m² | 131.5 | ×1.13 |
| `T` (sides on) | 115.4 kW/m² | 160.8 | ×1.39 |

**The drain is compensated, not eliminated** — it still grows with refinement,
though sub-linearly rather than the ×2 of a pure 1/dx term. What closes the gap
is that the heat now arrives at the receiving rock, not that the lateral flux
stopped. (Not directly comparable to D2o-0b's 302 kW/m² / 58 %: different
configuration. Reported because the packet requires it whatever the verdict.)

## 5. ⚠ Where I was WRONG, caught by my own pre-registration

**P1 predicted the null branch.** I wrote that `T_1mm`/`B_1mm` − 1 would land in
**+0.1 … +1.5 %**, that (a) would **FAIL**, and that **"any block above +4 %
refutes me"**. Measured: **+4.52 / +9.56 / +13.94 %** — all three blocks.
**Refuted on my own stated test.**

The error: I extrapolated linearly from the 2 mm legs, which were already in hand
when `CRITERION.md` was written (§0 of that file discloses this). The 2 mm run
has 24× too few of the faces the mechanism acts on, so it understates it
systematically. **P3 ("excess < 2 %") was wrong for the same reason.**

## 6. Caveats that survive the pass

- **(e) passes narrowly and monotonically.** +4.5 / +9.6 / **+13.94 %** against
  the pre-registered 15 % bound, rising block over block. The bound was fixed
  before the result and holds, but "bought some convergence by drilling faster"
  is **not fully excluded at 1 mm**, and a longer run would likely breach it.
- **No absolute number here is physics.** `side_face_factor = 0.2` below the
  plane is a frozen stand-in and the area treatment differs between regions.
  Only `T − B` and `T − S` are meaningful.
- **`B` carries no `side_P_face_above` / `wall_step_q` column**, because
  registering it would change `thermo.dat` and break item 1's byte comparison.
  `S` discharges criterion (c)'s "= 0" instead. Recorded in `CRITERION.md` §3.
- **THE ARBITER IS A FLAT PLATE — no shape result.** Cylindricity is the Meier
  follow-on.
- **`S_1mm` was not in the packet's five-leg table.** It was added because
  criterion (b) is the falsification guard and `S_2mm` had no discriminating
  power (the whole 2 mm effect is +0.1 %). Without it the result would have read
  "above-plane heating closes the gap" rather than "above-plane **side** heating
  closes the gap".

## 7. What this does and does not establish

**Does:** under an arbiter that isolates the cascade, heating the side faces of
above-plane columns — the faces the drain flows into, which have received zero
heat in every run ever done — closes the 6 / 8 / 14 % mesh gap to ≤ 5 %, and the
effect is attributable to those faces specifically and not to the column tops.

**Does not:** establish any drilling rate, hole shape or Meier comparison; the
arbiter is a flat plate. Nor does it establish that 142 W/m²K is the right
coefficient — it is the annulus derivation's value, used once, never varied.
