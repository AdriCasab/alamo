# D2ac — scoring criterion, FIXED BEFORE THE FIRST LEG

Hashed with **zero outputs on disk**. `PREFLIGHT.md`
(`6b4fc3cd2173a038446c4991b0256d9c1c52ec481e180eda91148d459815f593`) and
`preflight.py`
(`3be11f33c589ae4f4f30912cd5314c6471de298333b1c0ba941801b066d3a30d`) were frozen
first. **No re-windowing and no re-choosing a metric after a result is seen.**

**No source change in this packet.** `nozzle_descent` defaults to `prescribed`
(`:7584`, verified); `nozzle_z0` is `query_required` (`:7665`) and `nozzle_feed`
defaults to 0 (`:7666`). The pad-contact diagnostic is in `analyze.py`.
Binary `ea833f2f03855145a93eb98508e04a92cbc1aff80bfda40b1cc32bf96776c471`,
unchanged from D2ab's post-guard build.

## Legs

Base `A` (`d2w_attribution/output/A.log`) + D2ab's closure block, `lz` 0.4 → 0.5,
`amr.n_cell` z 200 → 250 (dz stays 2 mm and cubic), `nozzle_z0 = 0.550`.

| leg | feed [m/h] | side-face treatment |
|---|---|---|
| `FM` | 1.5 | inherited map |
| `FL` | 1.5 | `side_face_h=38.4`, `side_face_T=live`, `r_min=0.040` |
| `FH` | 1.5 | `side_face_h=169.3`, same |
| `FM13` | 1.3 | inherited map |
| `S0` | — | identity, all new keys absent |

**⚠ FORCED KEY DROPS** (item 0(d), measured not reasoned). `prescribed` aborts on
A's feet block: `:7593` refuses `jet_free_surface=1`, `:7596` refuses
`foot_body_clearance`. Ten keys are dropped — `jet_free_surface`,
`jet_fs_aspect`, `jet_fs_exponent`, `foot_r_inner`, `foot_r_outer`,
`foot_standoff`, `foot_rule`, `foot_npads`, `foot_pad_quantile`,
`foot_body_clearance`. With the full block dropped the unused-key list is EMPTY.
**`FM` therefore differs from `A` by more than the feed**, and that is stated in
every comparison. Bounded: free-surface dilution barely fired on the frozen legs
(`A` `jet_fs_cols` mean 0.019, `m_ratio` max 1.040; `SWc` exactly 0 and 1.000).

**⚠ The region columns `side_P_{bore,pads,wall}` are MEANINGLESS on these legs** —
`foot_r_inner/outer` fall to 0.0, so every charged face bins as `wall`. They are
not used for scoring. This is why (b) is computed offline.

**Imported by path, never re-run:** `A`, `SW`, `SWc`, `GL`, `GH`, `WLr`, `WLv`,
`WL`, `WL81`, `D`, `RB_D`.

## ⚠ The LZ trap

`d2v_hmap/analyze.py:30` and `d2x_areasplit/analyze.py:47` both hard-code
`LZ = 0.4` at module level, and `zface()` returns `LZ − nv*dz`. Against a 0.5 m
block that misplaces every depth by 100 mm. `analyze.py` binds `V.LZ`/`X.LZ`
**per leg** from that leg's own `.done` metadata before each call and restores
them after. The imported functions are untouched; only the constant they read is
bound per leg.

---

## (a) SCORED — THE FLOOR-FLUX EXCESS. This is the packet's product.

Per leg: realised **floor** recession rate at the realised standoff, against the
feed. Reported as the ratio.

**Pre-registered (item 0(a)): 2.15 – 2.60×** at a 100 mm standoff.
⚠ The packet states 2.22–2.60× from the fitted anchor (`SWc` ROP 3.694 over the
100–150 mm stretch); item 0 carries the endpoint anchor (3.578) as well, widening
the low end to 2.15. Both lie inside the packet's P2 band of 2.0–3.0×.

**This is expected to be a FAILURE against Meier and a MEASUREMENT of the missing
factor.** It is not a pass/fail bar on the model.

## (b) SCORED — pad contact, computed OFFLINE

Foot tips at `z_tip(t) = nozzle_z0 − nozzle_feed·t − 0.050`. Contact when the
**max top-face height over 28 ≤ r ≤ 40 mm** reaches `z_tip`. Computed from the
removal log with the same `zface()` the shape metric uses, at each output time.

Report: first contact time, contact fraction, worst interference [mm].
⚠ **`nozzle_collision_radius` is NOT the test** — it is unset in A's argv so it
defaults to 0.0 and cannot fire at all.

## (c) SCORED — Ø(z) over the calibrated section 150–300 mm

Against **85 ± 5 mm** (Tab. 8.2 **text**). The mouth (0–40 mm) and the ≈100 mm
underream are **reported, not scored** — Meier attributes the underream to a
pause, and a steady-feed run has none.
⚠ **No photo-derived diameter claim.** The figure does not support measurement
better than ±10 mm.

## (d) SCORED — volume ≤ 3.42 L per 0.5 m

Shortfall attributed to the absent mouth and pause.
⚠ **93.3 mm is NOT a second diameter target** — the water-filled volume includes
the mouth, the pause underream and the rounded tip, so a uniform 93 mm bore would
**contradict** Fig. 8.8.

## (e) SCORED — the standoff against the photograph's ceiling

≤ 100 mm passes. **215–222 mm (1.5 m/h) / 248–256 mm (1.3) is the pre-registered
failure.** Report the full trace.

## (f) SCORED — the power identity

`|jet_P_face − (tops + sides + side_loss)| ≤ 0.02 W` on every leg.

## (g) SCORED — identity

`S0` byte-identical to `d2z1_delete/output/RB_D`, header delta empty; tree gate
vs `witness/d2z1.md5`. ⚠ No source change, so the gate is expected to hit the
memoised cache on an unchanged binary hash.

---

## Pre-registered predictions

* **P1** the standoff settles at **180–260 mm**, failing the ≤100 mm ceiling.
  Item 0 derives 215–222 (1.5) and 248–256 (1.3), inside that band.
* **P2** the floor-flux excess is **2.0–3.0×**. Item 0 derives 2.15–2.60×.
* **P3** the pads make contact at 1.5 m/h (D2x put the pad-ring rate at ≈2.1 m/h;
  marginal). Refuted ⇒ the ring clears and the burner runs free, the stronger
  result.
* **P4** `FM` Ø in the calibrated section **exceeds 85 mm** — exposure scaling
  from `SW` (+7.6 mm over Ø80 at 3.694 m/h) gives ≈99 mm at 1.5 m/h.
* **P5** `FL`/`FH` are **OPEN, 80–95 mm, no preferred value**. `GL`/`GH` carried
  zero wall power so their exposure scaling is void. ⚠ **A landing near 85 mm is
  NOT to be presented as derived** — 85 is Meier's own number.
* **P6** volume **below 3.42 L**, the shortfall being mouth plus pause.

## Guardrails in force

* **The rate is an INPUT. Do not score burner descent, and do not tune to it.**
* Do not score the mouth or the ≈100 mm underream.
* Do not use `nozzle_collision_radius` as the contact test.
* Do not claim a photo-derived diameter; 85 mm is Tab. 8.2 text.
* Do not present `FL`/`FH` near 85 mm as derived.
* The pad ring stays the open-floor map (`PD` jammed on lateral recession).
* Every leg states its `side_face_area_law` (default is `step`, `:6819`).
* `spall.per_cell_removal` OFF; `side_face_h_max` never set.
* Build once; a rebuild requires re-running the identity leg and proving byte
  equality, as D2ab did.
