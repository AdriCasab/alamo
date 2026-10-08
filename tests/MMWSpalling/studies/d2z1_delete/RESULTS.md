# D2z-1 RESULTS — whole-cell removal unconditional, the key deleted

Criteria: `CRITERION.md`, sha256 `12db5609fec1652c…`, hashed read-only with zero
scored output on disk. Metric imported **unchanged** from
`studies/d2j0b_plane/analyze.py`. Frozen studies read by path, never re-run.

Legs (final binary): `W2` 1.8 min, `W1` 18.5 min, `RB_C1_2mm` 1.8 min,
`RB_D` 12.3 min. Full scoring in `analysis.txt`; audit in
`tree_repin_audit.txt`.

---

## Verdict

| criterion | result |
|---|---|
| **(a)** every test passes | **FAIL** — one assertion, `beam_void_closure` T3 |
| **(b)** tree re-pin audit | **PASS** |
| **(c0)** byte-identity vs D2z-0 | **PASS** |
| **(c)** cells voided +0.2…+1.0 % | **PASS** — +0.72 % / +0.32 % |
| **(d)** reported | +0.352 % / +0.438 %, gap shift +0.085 pts |

| prediction | outcome |
|---|---|
| **P1** zero test failures | **REFUTED** — one, and it is structural (below) |
| **P2** cells voided +0.2…+1.0 % | **HELD** |
| **P3** `n_voided == 1` everywhere | **HELD** |
| **P4** message + witness, run proceeds | **HELD** |
| **P5** no unexpected column change | **HELD** |

**The packet landed. One assertion is red, and it is red because the behaviour
it asserts has been deleted on purpose.**

---

## (c0) — the deletion is behaviour-preserving, by measurement

`W2`/`W1` against D2z-0's `F2`/`F1` (same physics, key present there, absent
here): **removal logs byte-identical** (1.4 MB, 10.5 MB), **all 24 cell
plotfiles identical** at both meshes, thermo identical on all 28 shared columns.
The only header delta is two **removed** columns, `firings_per_void` and
`removal_increment_cells`, both deliberately dropped.

## (b) — the re-pin audit

`witness/tree0922.md5` (superseded, preserved) vs `witness/d2z1.md5` (new
reference, declared in `witness/REFERENCE.md`):

* 2639 files both sides; **no path appears or disappears**;
* 837 differing: 720 `Cell_D`, 44 `thermo.dat`, 41 `h_col_events.csv`,
  32 `removal_events.csv`;
* **exactly one file gains a column**: `unit/beam_void_closure/output/thermo.dat`
  gains `sn_recession_superseded`, registered only under
  `spall.surface_normal = 1`. No other header change anywhere.

Goal 4's revision earned this: with `firings_per_void` dropped, the audit needs
**no exception carved out for our own addition** — 63 differing `thermo.dat`
files fell to 44, and the bullet is absolute rather than "clean except for ours".

**Bullet 4 (`n_voided == 1`)** holds tree-wide: every differing removal log now
reads rows == cells voided, ratio 1.000. On the arbiter, 9,606 and 69,702
regime-1 rows, `n_voided ∈ {1}` on every one.

**Bullet 3 (sign consistent) — verified where recoverable, and a tooling gap.**
The sweep overwrites the outputs it hashes and the reference stores hashes only,
so old counts survive only where a witness log printed them:

| test | old (half-cell) | new (whole-cell) |
|---|---|---|
| `scalar_flaw` | 52 rows, `min h` 9.844e-04 m, **26 cells** | 26 rows, `min h` 2.000e-03 m, **26 cells** |
| `beam_void_closure` | 18–18 layers | 18–18 layers |
| arbiter `W2`/`W1` | — | **+0.72 % / +0.32 %** |

All "up or equal". For the remaining unit tests the old counts are **not
recoverable** — a real limit of an md5-only reference, which can answer "did this
file change" but not "did it change in the right direction".

`scalar_flaw` is the packet in one line: **same 26 cells removed, half the rows,
`h` exactly `dz`.**

## (a) — the one failure, and why it is not a defect

`beam_void_closure` **T3**: `phi_top` must be observed in **both** sub-cell
halves; it now reads `[0, dz/2) = 2400, [dz/2, dz) = 0`.

Measured: **`phi_top` is exactly `dz/2` at every plotfile, with zero drift** —
whole-cell conserves it bit-for-bit, because `h = dz` exactly. The histogram is
therefore single-valued and **cannot populate both bins under any parameter
choice of this test**; which bin it lands in is decided by one ulp.

⚠ **This is not a generic invariant. T3's own header says the rejected formula
`floor((phi-dz/2)/dz)` double-deposits in the lower half and "the original leak
lived in the upper half". The bin now reading zero is `[dz/2, dz)` — the half
E1 was in.** So what is lost is the guard on the exact region of a bug that
shipped, and T3 was **not** relaxed.

Two restorations were attempted and both fail on this test:

* **`prescribed_T` variant** — unavailable: the test runs
  `surface_patch.enabled = 0`, so there is no patch to put in that mode. Adding
  one changes the driver of a deliberately beam-only test.
* **vapour-driven variant** — `spall.enabled=0 vapor.enabled=1 vapor.A_beam=1.0
  vapor.T0=300.0` ran to completion but removed **nothing**, so it does not
  produce the state either. Tuning vapour parameters until it did would be
  inventing a physics configuration to satisfy a test.

⇒ **Restoring this coverage needs a designed sub-cell removal configuration and
is a packet of its own.** Left red and documented, deliberately.

### The four assertions that are no longer red

* `spall_event` (5 assertions) — **a real defect in this packet, found and
  fixed**: the override was over-broad. It was catching the damage-law /
  `spall.prescribed_ratio` path, where `h_spall = 2·dz` is a *calibrated* two-cell
  thickness, not the `sp_weibull` depth scan's half-cell artefact. Now scoped
  `sp_weibull_spall && reg_choice == 1 && !sp_weibull_presc_T_face`, with the
  `n_voided == 1` assertion scoped identically or a damage-law run would trip it
  as a hard abort. **All four study legs are byte-identical across that rebuild**,
  so the added term provably touches only the damage-law path.
* `jet_enthalpy` (a) — resolved by **dropping `firings_per_void`**, not by
  editing the test.
* `jet_enthalpy` (d) — **diagnosed, not edited around.** The column fires and
  voids (3 events, 2 in window, `Sp_top` 1.0001): not starvation. `fit_rop`
  returns a **sentinel 0.0 below 3 samples**, a threshold implicitly calibrated
  on *firings*, which whole-cell halves. From those 2 rows the rate is
  **0.7279 → 0.7178 m/h** against the half-cell run's 0.7120 → 0.7022, and the
  bin-width change the check scores is **−1.387 % against the old −1.378 %**.
  Physics unchanged; threshold lowered to 2, which restores the original
  *physical* sampling requirement.
* `scalar_flaw` — `removal_rows > h_col_rows` → `==`. That inequality *was* the
  twice-per-cell artefact this packet deletes.

## Stage C — done

Re-baselined under whole-cell, written under **this study's** `output/`; the
frozen studies are untouched:

| anchor | leg |
|---|---|
| `d2l_scan/A_f02`, `A_f02_1mm` | **`W2` / `W1`** — same argv, so the gate pair *is* the re-baseline |
| `d2j0b_plane/C1_2mm` | `RB_C1_2mm` |
| `d2x_areasplit/D` | `RB_D` |

Total **34 min**, not the estimated 2.5–3 h, because `W2`/`W1` already are the
`d2l_scan` pair. Meier `A`/`M1` left frozen, as the packet directs.

## What D2z gets

A column now fires **once per cell**: no void-event filter, no half-cell warning,
`delta = |n + phi_S − phi_T|` built on firings that are removals. The tree has
**one** removal increment and no key.
