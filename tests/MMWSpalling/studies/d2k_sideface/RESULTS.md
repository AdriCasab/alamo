# D2k results — side-face flux on exposed vertical faces (2026-09-21/22)

**Two separable outcomes, and they must not be read as one.**

1. **The arbiter PASSES.** On the D2j-0b C1 cascade the side-face term removes the mesh
   divergence outright: disk-mean block gaps go from **−13.5 / −28.7 / −45.5 %** (key off) to
   **+3.4 / +1.5 / −3.2 %**, no band stalls, no inward-walking front. It does this drawing only
   **7 %** of face power, and the pre-registered f = 0.7 sensitivity moves the rate by **0.8 %**.
2. **The conditional Meier gate FAILS**, and the run **refutes the pre-registered direction**:
   the burner rate rose **2.7×** (4.07 vs D2i's 1.49 m/h) instead of falling a few per cent, and
   both legs **drilled out of the 0.40 m domain** before the D2i divergence window was reached.
   The gate is therefore FAIL by the written rule and **inconclusive on the question it was
   built to ask**.

The cause of (2) is a modelling assumption in the packet's rule, not a code defect: a side face
is given the **same h as the top face**, i.e. full stagnation-impingement heat transfer on a
vertical wall. In the Meier geometry the exposed wall below the nozzle plane is comparable in
area to the footprint, so absorbed power nearly doubles. **`f` corrects the staircase area, not
`h`, and `f` was not tuned.**

## 0. Hashes and binaries

- **Binary before:** `6d047b50bdcde619952aff907d88f89ed0cae8f9d3ed645b43ba0ac2cd647005`.
  **After (this step builds):** `62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`,
  verified unchanged at the end.
- **Arbiter criterion** `CRITERION.md`, sha256
  `6322064a1287e4a00b9f07eca1b8d8ec869a31952542070d06250cee9d067e61`, 2026-09-21 19:21:36,
  read-only, written before the output directory existed.
- **Meier criterion + pre-registered direction** `CRITERION_MEIER.md`, sha256
  `90af5275f0e2028742a48e682673f879e3458a15a26895a69a94d26a0da646bf`, 2026-09-21 22:16:23,
  read-only, written before any `LI0` output existed.
- Frozen studies (d2c … d2j0b) imported by path, none edited. The witness sweep was **copied**
  into `witness/` here rather than run in `d2e_mesh/witness/`, which is frozen.

## 1. Reference-set identity with the key off — the gate on the source change

**PASS: all 2639 files byte-identical** against the new binary
(`witness/d2k.md5` vs `witness/ref_post_d2e.md5`, `diff` empty). Every test in the sweep rc = 0:
`jet_d2c`, `dev2d`, `robin_face`, `robin_pinned`, `robin_jet`, `jet_enthalpy`, `robin_feet`,
`beam_void_closure`, `scalar_flaw`, `spall_event`, `s1_smoke`.

So `side_face_flux = 0` is inert and every result below is attributable to the new term alone.

## 2. Unit tests — `tests/MMWSpalling/unit/side_face_flux/`, 12/12 PASS

| check | result |
|---|---|
| (a) 1-D column, key on ≡ key off | all 25 common thermo columns, 24 `Level_0 Cell_D` files and the removal log byte-identical (key-on adds only `side_P_face`, `side_P_loss`, `side_faces`) |
| (b) ON ≠ OFF; ledger closes | max `ledger_err` **1.9e-14** (ON), 3.2e-14 (HALF), limit 1e-13 |
| (c) `side_faces` vs geometry | matches a staircase reconstruction from `removal_events.csv` on every sampled row |
| (d) exact linearity in f | `side_P_face`(ON) = 2 × (HALF) to **rel 0.00e+00** |
| (e) admissible range | 0 < 2.7 W ≤ 11.3 W bound |
| (f) exclusion | flame gain **identically zero** above the plane (`side_P_face + side_P_loss == 0.0` exactly) while the face still radiates; ON gains on every face row |

(f) is sharper than the packet asked for. An excluded side face takes nothing from the flame but
still loses heat to ambient — exactly the treatment an excluded **top** face gets, because both
go through the same `SurfaceCellFlux` call. My first version of the check asserted
`|side_P_face| < 1e-20` and failed; the code was right and the expectation was wrong.

## 3. The arbiter — PASS on all three conditions

| block | 2 mm | 1 mm | **D2k gap** | D2j-0b C1 (key off) | C0 |
|---|---|---|---|---|---|
| 50–100 s | 1.5971 | 1.6509 | **+3.4 %** | −13.5 % | −10.4 % |
| 100–150 s | 1.6207 | 1.6455 | **+1.5 %** | −28.7 % | −15.3 % |
| 150–200 s | 1.6761 | 1.6233 | **−3.2 %** | −45.5 % | −20.1 % |

- **(a) PASS** — every deciding block within 5 %.
- **(b) PASS**, by the better branch — **no band stalls at all**. Every band is still firing at
  250 s and sits 22–24 mm *below* the plane. In the key-off C1 the bands stopped at 40 / 68 / 91
  / 118 / 139 / 162 s with Δ_own = +12…+15 s, i.e. they stopped while still 5–6 mm below it.
- **(c) PASS** — `r_front` absent at both meshes; C1 walked 26 → 22 → 18 mm.

**Scoring correction, made after seeing results — a reviewer should check this.** My script first
reported (b) FAIL because it computed Δ_own = t_cross − t_stop = ∞ and tripped a `> 10 s` test.
But `t_stop = t_end` is D2j-0b's **sentinel for "still firing"** (it prints `250 (firing)`), and a
band still firing has not stopped. I added `t_stop ≥ t_end − 5 ⇒ not stalled`, which restores the
criterion's stated meaning ("a band stops firing before its own face crosses the plane"). The
criterion text is unedited and still hashed at `6322064a…`.

### Side power and the f sensitivity (reported, not deciding)

| run | block | `side_P_face` | `patch_P_robin` | side/top | `side_faces` | max abs `ledger_err` |
|---|---|---|---|---|---|---|
| C1_2mm | 150–200 s | 29.0 W | 387.8 W | 7.5 % | 14 | 2.9e-15 |
| C1_1mm | 150–200 s | 29.6 W | 387.1 W | 7.6 % | 57 | 5.5e-15 |

Face count 13 → 57 ≈ 4× as 1/dx predicts, while the *fraction* of power is nearly
mesh-independent. **f = 0.7 vs f = 1.0 at 2 mm: −0.8 %, −0.0 %, −0.8 %** — the result is not
sensitive to the one number that could have been fitted.

Figure: `d2k_bands.png`. Full tables: `analysis.md`.

## 4. The conditional Meier pair — FAIL, and inconclusive

`LI0_2mm` / `LI0_1mm`: the D2i `LI0` keys (D2f Q09 support rule, `pinned_idle_cycles = 1e9`) with
`side_face_flux = 1`, f = 1, to 450 s. 2 mm: 449.95 s in 10.7 min. 1 mm: 449.95 s in 154.4 min.

| block | 2 mm | 1 mm | **D2k gap** | D2i (side off) |
|---|---|---|---|---|
| 150–250 s | 4.066 | 4.179 | **+2.8 %** | −2.1 % |
| 250–350 s | 4.072 | 4.212 | **+3.5 %** | −4.0 % |
| 350–450 s | 0.278 | 0.041 | **−85.2 %** | −16.3 % |
| matched window 150–450 s | 3.303 | 3.239 | **−1.9 %** | — |

- **(a) PASSES** (−1.9 %, limit 5 %). **(b) FAILS** on the last block. **Verdict: FAIL.**
- **The failing block measures nothing.** Both legs have drilled out of the box — `nozzle_z` ends
  at 52 mm (2 mm) and 51 mm (1 mm) above a floor at 0 — so 350–450 s is the bottom stop, not the
  mesh. **I have not re-windowed to rescue it.**
- **The test never reached its own question.** D2i passed 150–250 and 250–350 s too (−2.1, −4.0 %);
  its failure appeared at 350–450 s (−16.3 %) and 450–550 s (−36.1 %). That window is exactly what
  drilling out of the domain destroyed here.

### The pre-registered direction is refuted

`CRITERION_MEIER.md` predicted the centre rate would fall a few per cent and the hole would widen
5–15 mm in diameter. Neither held.

| prediction | outcome |
|---|---|
| burner rate **down** a few % | **up 2.7×** (4.07 vs 1.49 m/h) |
| hole **wider by 5–15 mm** Ø at the feet | **+1.4 … +4.0 mm** at matched feet height |
| `jet_P_face` ≤ `jet_P_cap` per block | **holds**, max ratio 0.41 (2 mm) / 0.38 (1 mm) |

**It is not a budget leak** — which is what my own written rule said a rise would imply, and which
I therefore checked first. `jet_P_face` equals `ledger_P_robin` to 0.1 W (3204.2 vs 3204.2 W over
150–250 s), so the march accounts for every watt the rock absorbs, and the cap is never reached.

**What actually happened.** The jet was **area-limited, not enthalpy-limited**. Side faces carry
**39 % (2 mm) / 35 % (1 mm)** of face power and nearly double the absorbed power (3204 W against
D2i's 1639 W). And the energy does not go where I predicted: compared at the **same feet height**
the hole is only ~4 mm wider, because the feet rest on the highest rock under the pads and side
heating attacks exactly those protruding staircase columns. **The term acts as a feet-descent
accelerator, not a hole-widener.** Comparing at matched *time* (108 vs 78 mm Ø) is misleading —
that is the D2k hole having drilled three times further, and must not be quoted as widening.

## 5. Verdict and what the next packet inherits

- **The mechanism identified by D2j-0b is real and the remedy works on it.** The step-face flux
  was the mesh-divergent term, and heating those faces removes the cascade at constant f.
- **The rule as specified is too strong for a real hole.** Applying the impinging-jet `h` to a
  vertical wall doubles delivered power in the Meier geometry. The next packet needs a **wall
  heat-transfer closure** — a grazing-flow `h_wall(r, s)` well below the stagnation value, or an
  explicit `side_face_h_scale` — before any Meier number is scored. That is a *physics* decision
  and must not be made by tuning `f`, which is an area correction.
- **A Meier re-test needs a taller domain.** At 4 m/h the 0.40 m box is consumed in ~400 s. Either
  fix the wall `h` first (which will bring the rate back down) or use the 0.70 m domain.
- **Nothing from D2b–D2d may be quoted as a Meier result**, and D2g Stage 2 remains unplanned:
  the support rule is still only mesh-checked to ≈ 350 s.
- Keys are default-off and the 2639-file set is byte-identical with them off, so the source change
  is safe to carry regardless of the Meier outcome.
