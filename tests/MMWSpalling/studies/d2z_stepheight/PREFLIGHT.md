# D2z PRE-FLIGHT — item 0, read-only, no build, no runs

Everything here was produced from frozen/re-baselined event logs only. Nothing
was rebuilt and no leg was launched.

## Item 0a — `delta` re-measured on the RE-BASELINED anchors

Every number in the packet was measured on **half-cell** logs. Re-measured on
`studies/d2z1_delete/output/{W2,W1}` (whole-cell) with the both-column estimator
and the **trailing**-interval phase the live code must use, filtering
`n_voided > 0.5`. Estimator: `studies/d2z_stepheight/delta_offline.py`.

| leg | window | step % | ≥2-cell % | birth med | **t-wtd MEAN** | fallback % |
|---|---|---|---|---|---|---|
| `W2` 2 mm | 50–100 | 6.5 | **0.0** | 0.049 | **0.137** | 0.00 |
| | 100–150 | 8.2 | **0.0** | 0.060 | **0.159** | 0.00 |
| | 150–200 | 9.1 | **0.0** | 0.068 | **0.171** | 0.00 |
| | **50–200** | 7.9 | **0.0** | 0.059 | **0.157** | 0.00 |
| `W1` 1 mm | 50–100 | 8.2 | **0.0** | 0.062 | **0.197** | 0.00 |
| | 100–150 | 11.4 | **0.0** | 0.082 | **0.262** | 0.00 |
| | 150–200 | 14.1 | **0.0** | 0.097 | **0.319** | 0.00 |
| | **50–200** | 11.2 | **0.0** | 0.079 | **0.270** | 0.00 |

**The band survives the re-baseline**, as the packet expected. The packet quotes
0.16–0.20 (2 mm) and 0.19–0.31 (1 mm); re-measured, 2 mm reads **0.137–0.171**
(marginally below the quoted band) and 1 mm **0.197–0.319** (matching). The
values carried into `CRITERION.md` are the re-measured ones.

Two facts worth stating separately:

* **`≥2-cell` is 0.0 % at BOTH meshes** on the arbiter's inner disc. The arbiter's
  floor is fully in the phase-lag state; the multi-cell population the packet
  worries about (53.2 % on Meier `A`) does not exist here. **(f2) must therefore
  be read on `SW` before any claim about `A`'s floor.**
* **`delta_fallback_cols` will be ≈ 0**: the fallback fraction is **0.00 %** in
  every block at both meshes, far under the 5 % bar. A column with fewer than two
  voids is essentially absent once the run is past its transient.

Mesh trend, which is the packet's whole argument for the closure: the
time-weighted `delta` rises **0.157 → 0.270** from 2 mm to 1 mm, so **no fixed `f`
can serve both meshes**, and `f = 0.2` sits between them.

## Item 0b — the gate location (REPORT ONLY; not a blocker, not fixed here)

`Sp_top >= 1` is tested on the **top cell** (`Removal.H:650-652`) while the scan
measures depth from the **face**. Measured on the whole-cell legs:

| leg | `h_scan/dz` median | range | `h_applied/dz` | `Sp_top − 1` med | `(dz/2 − a_f)/dz` |
|---|---|---|---|---|---|
| `W2` 2 mm | **0.50391** | 0.496–0.520 | **1.00000** | 0.00336 | 0.49804 |
| `W1` 1 mm | **0.50000** | 0.492–0.523 | **1.00000** | 0.00238 | 0.49607 |
| `RB_C1_2mm` | **0.50000** | 0.496–0.508 | **1.00000** | 0.00160 | 0.49804 |

**The landing is `dz/2` at both meshes and on the cascade reference.** It is
**mesh-invariant in units of `dz`**, i.e. it scales with the cell and is not a
physical length — which is what a gate evaluated at the cell centre against a
depth measured from the face would produce. Together with D2z-1's `prescribed_T`
control (a real 1100 K face against a 655 K top cell, still `dz/2 − a_f`), the
lerp is not the cause.

⚠ **Scope confirmed, and it is not this packet's problem:** `h_applied/dz` is
**exactly 1.0**, so whole-cell discards the scan's depth on the computed-T spall
path. This cannot reach D2z's legs or the Meier line. It affects **`prescribed_T`
spall (Hu, Kant)** and the damage-law path. Recorded for whoever revisits
Step 16d; **not fixed here.**

## Reachability ledger — verified, with citation drift noted

All sites exist. Line numbers drifted from the packet by ~+11 because of D2z-0 /
D2z-1 edits to the same file; the packet's own exact citations (`k_face`,
`WallStepDiagnostic`, the registration) are unchanged.

| site | packet | actual |
|---|---|---|
| kernel deposit | `:1675-1678` | **1686-1689** |
| march weight | `:2218-2221` | **2229-2232** |
| **march factor** | `:2354-2355` | **2365** |
| ledger twin | `:2958-2961` | **2969-2972** |
| `RiserWeight` | `:2720-2725` | **2731** |
| **`k_face` (conductance target)** | `:1537-1571` | **1537 / 1557 / 1570** ✓ exact |
| `WallStepDiagnostic` | `:1861-1905` | **1861** ✓ exact |
| `wall_step_q`/`_n` assignment | `:1872-1873` | **1910-1911** |
| drain registration | `:692-693` | **692** ✓ exact |

⚠ **The packet contradicts itself on the conductance site.** The Sources table
and Step 3 both name **`MMWSpalling.H:1537-1571` (`k_face`)** and Step 3 states
explicitly that an earlier draft's `Removal.H:1441-1452` "is wrong — that line is
`phi_new = phi_old - h_c`, the level-set shift, i.e. recession, not conduction."
But the **reachability ledger table still carries the stale
`conductance Removal.H:1441`**. Step 3's correction is taken as authoritative and
verified: `Removal.H:1444` is indeed `phi_new = phi_old - h_c`.

## Step 3b — the per-cell split, derived before being trusted

The packet's split — `(1 - phi_T)` top, 1 per middle, `phi_S` lowest — **does not
sum to `delta`**, and the packet says to derive it rather than take it on trust.
For an `n`-cell step the taller column's exposed cells run from `k_top(S)+1` to
`k_top(T)`, which is exactly **`n` cells**:

```
  packet's split : (1 - phi_T) + (n-2)*1 + phi_S = n - 1 - phi_T + phi_S
  delta          :                                  n     - phi_T + phi_S
  -> short by exactly ONE CELL
```

The split that closes is **top `(1 - phi_T)`, each middle 1, and the lowest
`1 + phi_S`** (the `phi_S` piece folded into the lowest exposed cell, since
geometrically it belongs to a cell with no void neighbour):

```
  (1 - phi_T) + (n-2)*1 + (1 + phi_S) = n - phi_T + phi_S = delta   ✓
```

This is asserted per pair in code (pieces sum to `delta` to round-off), as the
packet requires.

## Derived constants (no new parameter)

`a = dz * v / alpha`, `v = dz / dt_i`, `alpha = kappa / (rho*c_p)` — the thermal
**diffusivity**, not the codebase's `kappa` (conductivity). Taken from where they
are consumed: `mat->kappa(T)` and `mat->effective_rhocp(T)`, both evaluated at the
**same cell temperature the conduction stencil uses** (`T_old(i,j,k)`), stated
here as the packet requires.

Sanity, at the arbiter's measured cadence (median flip interval 4.48 s at 2 mm):
`v = 2 mm / 4.48 s = 4.46e-4 m/s`; `alpha/v = 6.9045e-7 / 4.46e-4 =` **1.55 mm**,
matching the packet's "1.6 mm, NOT 3.2". Then `a = dz/(alpha/v) = 1.29`, and at
the re-measured `delta = 0.157`:

```
  g(0.157) = (1 - e^{-0.157*1.29}) / (1 - e^{-1.29}) = 0.25   ->  4.0x per-face
```

against the packet's "≈ 0.28, a 3.5x reduction" — consistent, and the small
difference follows from the re-measured `delta` being slightly lower than the
quoted band.
