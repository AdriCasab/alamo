# D2t criteria — fixed and hashed BEFORE the first campaign leg launched

## 0. Disclosure

**No campaign leg (`P_2mm`, `P_1mm`, `PA_2mm`, `PA_1mm`) has been launched when
this file is hashed.** D2r-1 hashed after its 2 mm legs had finished; the packet
says not to repeat that, and it is not repeated.

What *has* run: the two ~10 s trace legs `TR_below` / `TR_above` of Goal item 2,
whose own criteria were hashed separately in `TAPER_TRACE.md`
(`abf5fdc167783e05...`) before any temperature was read. They passed T1–T8:
the compiled `h(r)` is the intended raised cosine (self-calibration returned
**699.99 W/m²K against a known 700**), it applies to **both** branches of
`if(s>0,700,142)` (ratio **0.202857** = 142/700 exact), and **nothing is heated
beyond `surface_patch.radius`** (max ΔT over 1106 outside columns =
**1.13e-6 K**). N1 and N2 are therefore excluded before launch, not after.

## 1. The taper, as compiled and verified

`h(r, s) = taper(r) · if(s > 0, 700, 142)`, raised cosine, full strength to
`r_flat = 18 mm`, at its floor (`h = 1.0e-3`) by `r_taper = 34 mm`,
`W = 16 mm` = **8 coarse / 16 fine cells**, a fixed *physical* width.
`surface_patch.radius = 50 mm`, so the hard mask sits **8 coarse cells** into
the floor region. Landmarks used below, computed from the taper, not measured:

* `h >= 0.95 ×` peak for **r <= 20.30 mm** — with the measured core margin of
  only +5.5 % over feed, this is where columns start to fall behind the plane.
* `h >= 0.50 ×` peak for **r <= 26.00 mm**.
* Columns in the 19–34 mm band: **156 at 2 mm, 624 at 1 mm** (direct grid count).

## 2. The shape metric — no disk-mean rate is scored in this packet

1. Final `k_top` per column from `_removal_events.csv`: `k_top = k_top(0) −
   Σ n_voided`, columns with no event keep the initial top.
   `z_face = (k_top + 1)·dz`, so `z_face = LZ − Σ n_voided · dz`.
2. Binned by radius into **fixed 2 mm physical bins**, identical at both meshes.
   **Asserted in the analyzer** (N4: a mesh-dependent binning would make the
   metric itself mesh-dependent).
3. `z(r)` = mean `z_face` per bin; **bin population always reported**.
4. `D(z) = 2 · max{ r : z(r) <= z }` on a **fixed depth ladder: 10, 20, …,
   100 mm** below the original surface (`z = LZ − depth`). The ladder is
   absolute and fixed here, so no depth range is chosen after seeing the data;
   110 mm is reported as an extra, **ungated** row.
5. Radial rate decomposition: `r < 18` / `18–34` (taper band) / `>= 34` mm,
   by 50 s block.

## 3. Criteria

* **(a) THE SCORED ONE — SHAPE CONVERGENCE, `P_1mm` vs `P_2mm`.**
  * **Core:** `max |z_1mm(r) − z_2mm(r)|` over bins with `r < 18 mm`.
    **Bar: <= 4 mm** (2 coarse cells). Anchor: the core agrees on rate to ~3 %
    and 3 % of the ~110 mm recession is 3.3 mm.
  * **Band:** scored through **`D(z)` on the fixed ladder**, not per-bin `z(r)`.
    **Bar: <= 1 bin = 4 mm in diameter** at every ladder depth.
* **(b) THE SMOOTHNESS WITNESS.** An above-plane column count cannot separate
  "stalled against a cliff" from "lagging because `h` is lower", and the taper
  guarantees the second, so it is not used here.
  * largest step in `z(r)` between **adjacent 2 mm bins** inside
    `surface_patch.radius`, reported at both meshes;
  * **zero never-spalled columns with `r <= 26.00 mm`** (where `h >= 0.5 ×` peak).
    **This is the gated half.**
* **(c) DOES THE WALL CONVERGE WITH THE LEVER ON.** `PA_1mm` vs `PA_2mm` on the
  (a) metric, **same bars**. Separately, `PA − P` on `D(z)` is **reported as the
  first measurement of what annulus heating does to a wall slope** — it is *not*
  a null check.
* **(d) THE CORE IS NOT BROKEN.** Rate over `r < 18 mm`, block 150–200 s,
  against the reference **recomputed on the same mask**: `T_2mm` **1.6667**,
  `T_1mm` **1.6349** m/h (reproduced exactly in `PREFLIGHT.md` §1).
  **Bar: within 5 % of the recomputed reference.** A regression here is a taper
  bug, not a finding.
* **(e) THE LATCH — reported, not gated.** `h` steps **142 → 700, 4.9×, at the
  nozzle plane**. Below it a column drills 1.63 > 1.55 feed and stays below;
  above it drills far under feed and falls further behind, so the two states
  are separated by the feed rate and each column **latches**. Report per block
  and per mesh: above-plane column count, and number of plane crossings. This
  is D2o-0b's sharp-corner pathology moved from `r` into `s`.

## 4. Pre-registration, with numbers

* **P1 — the band becomes an above-plane wall**, ~**156 at 2 mm / 624 at 1 mm**
  in 19–34 mm (my grid count; the packet says 159/629). **This is expected and
  is NOT the cliff.** Refuted if the 1 mm count is under 100.
* **P2 — `PA − P` is LARGE, not small.** Refuted if `max |PA − P|` on `D(z)`
  is under 1 bin (4 mm).
* **P3 — (d) passes.** The core is untouched inside `r_flat`: `taper = 1` there,
  and the side-face term cannot add heat to a core column, because the
  side-face branch fires on a **solid cell with a VOID lateral neighbour** and
  the core column is the *deeper* one, whose neighbours at its own `k_top` are
  solid. The extra side flux goes to the **band** column, not the core.
* **P4 — (b) passes:** no step in `z(r)` larger than ~2 bins, and no
  never-spalled column inside 26 mm.
* **P5 — (a) IS GENUINELY OPEN. I do not predict it.** No shape metric has ever
  been run at two meshes in this project. Recording "unknown" is the honest
  pre-registration; inventing a number here is what killed D2s.
* **P6 — the latch may prevent band convergence entirely.** If (a)'s band bar
  fails while the core passes and (b) is clean, the reading is that **`h`
  stepping at the nozzle plane cannot represent a wall**, and the successor
  smooths `h` in `s`. **It will NOT be rescued by re-tuning the taper.**

## 5. Guardrails binding on the scoring

* **No disk-mean rate is scored.** (d) uses one as a *regression check on a
  known reference*, never as the result.
* **Bins are a fixed physical 2 mm at both meshes**, asserted in code.
* **The depth ladder is fixed here (10–100 mm)**, so no window is chosen after
  the fact.
* `side_face_factor = 0.2` is frozen and is **not scanned**.
* **Domain watch:** assert the deepest column stays **>= 10 mm above the domain
  floor**. Measured on the frozen sharp-edge legs the core reaches 110–112 mm
  of 150 mm, leaving 38–40 mm, so this is expected to hold comfortably; if it
  fails the depth result is void.
* **No Meier number, no validation claim, no shape-correctness claim.** The
  arbiter is a flat plate under a descending disk: it can say whether shape is
  **numerically convergent**, never whether the shape is **right**.
