# D2ad — scoring criterion, FIXED BEFORE THE FIRST LADDER LEG

Hashed after item 0 fixed the scales and measured C1, `G1600`'s `T_g` and `R90`,
as the packet directs. `PREFLIGHT.md`
`427d727be0548d38029331d6fcf876d8e365228624ca64fd974d184256a73eeb`,
`preflight.py`
`a684615cc66d30087625a111196fd4936f5eeb0278299ddeca77735f145d2197` — both frozen
chmod 444 first. **No re-windowing and no re-choosing a metric after a result.**

**NO SOURCE CHANGE.** Binary `ea833f2f03855145a93eb98508e04a92cbc1aff80bfda40b1cc32bf96776c471`,
to be re-hashed at the end.

## ⚠⚠ THE LADDER IS RECOMPUTED — mandated by the packet itself

Item 0 fired **both** of the packet's recompute triggers:

* **(c)** worst |measured − C1| = **53.5 K**, bar 40 K ⇒ *"recompute F, s* and
  q100 and re-hash"*. C1's linear `1900 − 235.3·F_q` is replaced by a quadratic
  fitted **exactly through the three measured points** (F_q = 0.35/0.50/1.00 →
  T_g = 1871.2/1834.2/1664.7 K). Every leg is interior to that range.
* **(e)** **`R90` = 4.237 m/h measured**, not 3.3 — **+28.4 %**, bar 10 % ⇒
  *"recompute the ladder and re-hash"*. The packet called 3.3 "the validated
  one"; it is not. ⇒ **`s* = 254.2·F` mm, not `198·F`**, and **contact when
  `F ≤ 0.3540`**, not 0.4545.

**Corrected window: `0.3540 < F ≤ 0.3934`, i.e. `s*` in (90, 100] mm.**

| leg | F_q | scale | T_s | T_g | F | s* mm | q100 | pre-registered |
|---|---|---|---|---|---|---|---|---|
| `C0` | 1.00 | 1.0 | 823.3 | 1664.7 | 1.000 | 254 | 1.02 | fail; above Kant's 0.87 |
| `G1600` | 1.00 | 1.0 | 823.3 | **1465.8** | 0.764 | 194 | 0.78 | fail — degeneracy witness |
| `Q5` | 0.50 | 1.0 | 823.3 | 1834.2 | 0.601 | 153 | 0.61 | fail |
| `Q5T870` | 0.50 | 1.3523 | 870 | 1834.2 | 0.527 | 134 | 0.59 | **the packet's predicted pass — FAILS** |
| `Q35` | 0.35 | 1.0 | 823.3 | 1871.2 | 0.436 | 111 | 0.45 | fail — nearest bracket above |
| **`Q35T870`** | **0.35** | **1.3523** | **870** | 1871.2 | **0.383** | **97** | **0.43** | ***** PASS ***** |
| **`Q30`** | **0.30** | **1.0** | **823.3** | 1882.1 | **0.378** | **96** | **0.39** | ***** PASS, other lever ***** |
| `Q30T870` | 0.30 | 1.3523 | 870 | 1882.1 | 0.332 | 84 | 0.37 | **CONTACT** — bracket below |
| `S0` | — | — | — | — | — | — | — | identity |

`G1600`'s `T_g` is **measured** (1465.8 K), not the packet's estimated 1410.
`kic.scale` 1.3523 → 870 K and 1.9088 → 920 K are **measured** by item 0(b)
(T_s = 821.3/864.6/895.5/919.1/950.1 K at scale 1.0/1.3/1.6/1.9/2.2).

**Dropped** from the packet's table: `Q55T870` (145 mm) and `Q5T920` (117 mm) —
both now fail high and carry less than the legs bracketing the window.
**`Q35T870` is restored**: the packet dropped it as "contact under C4", a verdict
that used the refuted `R90` = 3.3.

⚠ **`Q35T870` and `Q30` reach the same `s*` by different lever pairs** — that is
criterion (i)'s degeneracy, measured rather than asserted.

## ⚠⚠ THE GEOMETRY WINDOW AND THE FLUX GATE NO LONGER COINCIDE

Revision 2 had the predicted pass at `q100` = 0.55, the centre of Meier's
required 0.5–0.6, and called that coincidence "a real result" if it held. **With
the measured `R90` it does not hold.** The legs reaching the standoff window sit
at **`q100` = 0.39–0.43 MW/m²**, *below* Kant's measured spallation band
(0.45–0.87) and below Meier's required 0.5–0.6. **Registered here, before the
ladder runs.** If it survives, holding Fig. 8.8's standoff would require a floor
flux *below the threshold at which Kant measured spallation to occur at all* —
a result, not a defect, and **not to be resolved by moving either bar.**

---

## Criteria

**(a) SCORED — THE WINDOW, both halves.** (i) standoff plateau (slope < 2 mm per
100 s over the last 200 s) **≤ 100 mm**, AND (ii) no pad contact after t = 60 s
(max top face over r ∈ 28–40 mm against `z_n − 0.050`).
⚠ Contact at t = 0 is **by construction** and excluded. ⚠ No `nozzle_z` under
`prescribed`; standoff is **`jet_s_c`**. ⚠ **`C0` is scored 400–600 s** (C3).
⚠⚠ **Report `s*` with ±5 %; 95–105 mm is NOT a decisive verdict.** The
deliverable is **`F` to about ±10 % with an equivalent nozzle temperature**, not
a binary pass.

**(b) SCORED — THE FLUX GATE.** `q100 = F_q·1216·(T_g − T_s)` using each leg's
**measured** floor `T_g`, against **Kant 0.45–0.87** and **Meier 0.5–0.6 MW/m²**.
⚠ At a **common 100 mm standoff**, never at the leg's own `s*`.

**(c) SCORED — Ø(z) over 120–250 mm** within **85 ± 5 mm**; 150–300 mm reported
alongside. ⚠ Bind `LZ` **per leg** from `.done`.

**(d) SCORED — volume ≤ 3.42 L per 0.5 m**, shortfall stated.

**⚠⚠ (c) AND (d) ARE VOID FROM SUSTAINED PAD CONTACT.** Under `prescribed` the
burner descends on a schedule with no force balance and the collision guard only
covers `r ≤ 28 mm`, so pads at 28–40 mm **pass through rock** unphysically. On
any leg meeting (a)(ii)'s contact condition, report Ø and volume **only up to
first sustained contact** and say the rest is void.

**(e) SCORED — firing window** `p90 − p10 ≤ 30 K` on every leg, 150–400 s.
Item 0(b) already measured 10.5–18.8 K across the scale ladder (PASS).

**(f) SCORED — power identity** ≤ **0.02 W** every leg. Bar is above the 0.015 W
six-sig-fig floor; do not tighten.

**(g) SCORED — identity.** `S0` byte-identical to `d2z1_delete/output/RB_D`;
tree gate vs `witness/d2z1.md5` with **`--force`** (inputs-only change).

**(h) REPORTED — the calibrated pair** `(F_q, T_s)`: T_s against 770/870/1190/1273 K;
`q100` against Kant and Meier.

**(i) REPORTED — THE DEGENERACY, with a measured slope.** State the calibrated
`F_q` as an **equivalent nozzle temperature** using `G1600`'s **measured** flux
reduction as the slope, against the thesis's **1450–1850 K** and the
**1.9–5.6 kW cooling loss (−100…−300 K)**. This is a band, not a result.

## Pre-registered predictions

* **P1** (packet) the window closes at `Q5T870`. **Item 0 predicts this is
  REFUTED**: `Q5T870` lands at 134 mm. The window closes at `Q35T870`/`Q30`.
* **P2** `Q30T870` makes contact; `C0`, `G1600`, `Q5`, `Q5T870`, `Q35` plateau
  above 100 mm in descending `s*` order (254/194/153/134/111).
* **P3** `q100` above Kant's 0.87 on `C0` and the passing leg inside 0.5–0.6.
  **Item 0 predicts the third clause is REFUTED** — the passing legs sit at
  0.39–0.43, below both bands.
* **P4** the passing leg's Ø over 120–250 mm lands **82–87 mm**.
* **P5** the firing window stays ≤ 30 K. Already measured 10.5–18.8 (PASS).
* **P6** **a required T_s > 920 K is a REFUTATION of the rock-side explanation.**
  The corrected ladder needs only 870 K, so this is expected to hold.
* **P7** volume on the passing leg **2.6–3.0 L per 0.5 m**.
* **P8** the equivalent nozzle temperature lands **inside 1450–1850 K**; below
  1450 refutes the flux explanation.

## Guardrails in force

* **NO SOURCE CHANGE**; binary must still hash to `ea833f2f…`.
* Only three keys move: the `h_expr` prefix, `spallation.sp.kic.scale`, and
  `jet_T_nozzle` on `G1600` alone. `side_face_h` stays 169.3; `weibull.*`
  untouched.
* `nozzle_collision_radius = 0.028` on every leg (carried from D2ac's `FH`).
* The ten feet-only keys stay dropped.
* Bind `LZ` per leg from `.done` in every imported analyzer.
* Every launch under `caffeinate -dis`.
* AMR parked; 1 mm parked; no third mesh.
* The user runs `git add/commit/push`. Do not.
