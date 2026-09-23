# D2r-1 pre-flight item 0 — call sites, and proof every leg reaches them

The packet makes this mandatory: *"name the call site of every new or changed
behaviour and demonstrate, from the key set of each leg, that it is reached. If
any leg cannot reach it, STOP and report."* Three consecutive packets died here
(D2r-0 anchored to unchecked C0 keys; D2r-0b specified a cap the code ignores in
the mandated state; D2r-0c specified a pass inside `JetEnthalpyMarch`, which the
arbiter never calls). **This one is demonstrated by measurement, not argument.**

## The arbiter's key set, as the harness actually emits it

```
B_2mm/B_1mm  surface_patch.h_expr="if(s>0.0,700.0,0.001)"
             surface_patch.side_face_flux=1  surface_patch.side_face_factor=0.2
T_2mm/T_1mm  surface_patch.h_expr="if(s>0.0,700.0,142.0)"
             ... side_face_factor=0.2  side_face_factor_above=1.0
S_2mm        ... same as T, side_face_factor_above=0.0
```

## Call sites

| # | changed behaviour | call site | predicate that must hold |
|---|---|---|---|
| 1 | per-region side factor `sf_c` | conduction kernel, side-face branch: `if (col_in_patch)` → `if (nfx + nfy > 0)` | `side_face_flux = 1`; `sp_on`; **not** `PrescribedTemperature`; column inside `surface_patch.radius`; cell has an exposed lateral face (needs void ⇒ removal on) |
| 2 | the same factor in the **ledger twin** | mirrored side-face block | identical predicate |
| 3 | `side_P_face_above`, `side_cols_above` | ledger twin, gated `side_face_factor_above >= 0` | only `T`/`S` set the key |
| 4 | `WallStepDiagnostic` | `BuildFlameColumns` tail | same gate |

**Why each predicate holds on every leg:** the arbiter is
`mode = convective_flame` (so not the prescribed branch), `surface_patch.radius
= 0.030` contains the whole r < 30 mm scoring disk, spall removal is on so void
cells and therefore exposed lateral faces exist, and every leg sets
`side_face_flux = 1`.

⚠ **Site 4 is the one that killed D2r-0c.** It sits at the tail of
`BuildFlameColumns`, which the arbiter **does** execute, rather than inside
`JetEnthalpyMarch`, which is guarded by `if (surface_patch.jet_enthalpy)` and
which the arbiter — carrying no `jet_closure` — never calls.

## Demonstration (measured, 2 mm, 250 s)

| leg | `side_P_face_above` max | `side_cols_above` max | `ledger_err` max |
|---|---|---|---|
| `T_2mm` | **43.57 W** | 103 | 4.3e-15 |
| `S_2mm` | **0.00 W** | 102 | 5.0e-15 |

Both legs see the same ~103 above-plane faces; only the factor differs. So the
branch is reached, the key controls exactly the intended term, and forcing it to
zero removes that term and nothing else. **No leg fails to reach its call site;
the packet does not stop here.**

`ledger_err` at 1e-15 is the independent proof that site 2 mirrors site 1: if the
twin had not been updated in step with the kernel, the ledger would stop
reproducing the flux and this number would blow up.

## Goal item 1 — the baseline-validity proof: **PASSED**

`B_2mm` re-run with the current binary reproduces D2l's frozen `A_f02`
**byte for byte**: `_removal_events.csv`, `thermo.dat`, and every compared
`Level_0/Cell_D_00000`. This mattered because **no input in the 2639-file
identity set turns `side_face_flux` on**, and D2q-2b/2c both modified that path —
so byte-identity of the tree said nothing about this baseline. It does now.
