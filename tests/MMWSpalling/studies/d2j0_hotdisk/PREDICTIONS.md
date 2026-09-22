# D2j-0 pre-registered predictions — the hot-disk edge test (written before launch; never edit)

Binary `6d047b50…`, no build. Configuration: quarter domain 0.08 × 0.08 × 0.10 m, axis at the
xlo/ylo corner, fixed Robin disk, **no jet, no feet, no descent**, `pinned_idle_cycles = 1e9`,
`t_end = 150 s`, `dt = 4 ms` at every mesh, `weibull.V0 = V_cell = dz³` (the scored convention,
so it moves with the mesh). Sharp set `S4`/`S2`/`S1` (+ conditional `S05`) at dx = dz = 4, 2, 1
(, 0.5) mm with `h_conv = 700`, `T_flame = 1600` inside `radius = 0.030`. Tapered set
`T4`/`T2`/`T1`, identical but with h ramped by a **raised cosine** from 700 at r = 27.5 mm to 0
at r = 32.5 mm, patch radius 0.0325 so the ramp is not truncated.

## Hand estimate of the operating point (not a prediction under test)

At a pinned surface of 821 K (the A1/S1b firing temperature for this Weibull draw):
q_net = 700·(1600 − 821) − 0.8σ(821⁴ − 293.15⁴) − 10·(821 − 293.15) = **0.520 MW/m²**, so
v = q_net/(ρc_p·ΔT_fire) = **1.632 m/h** and ≈ 68 mm of recession in 150 s, leaving ≈ 32 mm of
rock under the floor of a 100 mm box.

**The thermal length is shorter than the packet assumed.** α/v = 6.9022e-7 / 4.53e-4 =
**1.52 mm**, not 3.2 mm, so the ablation layer spans **0.38 cells at 4 mm, 0.76 at 2 mm, 1.52 at
1 mm and 3.0 at 0.5 mm**. That makes the mesh question sharper, not weaker: on this estimate
**no mesh in the primary sweep resolves the layer**, and `S05` is the first that begins to.
Recorded here before the runs so it cannot be read back into the result.

## Metrics (fixed now; no re-windowing or re-binning after the fact)

All from the per-column cumulative recession in `<pf>_removal_events.csv`, which logs every
removal with its column and its applied increment `h_applied`.

- **Window 50–150 s** (the first 50 s are the start-up transient). Per column,
  `v_col` = Σ `h_applied` over events with 50 ≤ t ≤ 150 s, divided by 100 s, in m/h. Columns
  with no event score 0 and are kept.
- **Bins:** 2 mm annuli in the column-centre radius r, edges 0, 2, …, 40 mm. `v(r)` = the mean
  of `v_col` over every column whose centre falls in the bin.
- **`v_c`** = mean `v_col` over columns with r < 6 mm.
- **`r_half`** = the radius where `v(r)` crosses 0.5·`v_c`, by linear interpolation between bin
  centres, taking the outermost crossing.
- **`w_edge`** = r(0.1·`v_c`) − r(0.9·`v_c`), both by the same interpolation, outermost
  crossing. This is the primary mesh quantity.
- **`f_stall`** = the fraction of columns with r < 30 mm whose `v_col` < 0.1·`v_c`.
- **`r_stall_inner(t)`** for t = 50, 100, 150 s: recompute the binned `v(r)` over the trailing
  window [t − 50, t], then take the inner edge of the **contiguous** run of bins that reaches
  r = 30 mm and has mean v < 0.1·`v_c`. Report "none" if no such run exists. (At t = 50 s the
  window is [0, 50] and therefore includes the transient; it is a baseline, not a trend point.)

The silence heuristic used in D2h/D2i (no event within 2·t_cell) is **not** used here: every
metric above is a recession rate and needs no firing-time threshold.

## Predictions

- **P1 — `v_c` agrees within 5 % across 4 / 2 / 1 mm.** It is the 1-D energy balance, which S1
  showed to be mesh-independent at fixed flux. **A P1 failure means the configuration is wrong,
  not the edge: stop and say so**, and do not score P2–P4.
- **P2 — under the event model a stalled annulus exists at the patch edge at every mesh, and
  `f_stall` grows as dx falls.**
- **P3 — `w_edge` does not converge from 2 mm to 1 mm** (it changes by more than 20 %).
- **P4 — `r_stall_inner` moves inward with time at 1 mm**, and does so less, or not at all, at
  4 mm.

## The sharp-vs-tapered discriminator

If the sharp set shows a dx-scaling edge and **the tapered set does not**, then the D2i front is
driven by the **sharpness of the h edge** — in the Meier loop, by `jet_r_reach` — which points
at the gas closure rather than at the removal law. That is to be stated plainly if it happens,
because it would redirect the next packet away from the front redesign.

## The decision this packet makes

- **P2 holds →** the D2i mechanism reproduces in open loop, the continuous-front redesign is
  aimed at the right target, and **`w_edge`(dx) measured here becomes the pre-registered
  prediction** for the residual mesh gap in the later Meier pair. Next packet: D2j-a.
- **P2 fails** (no stalled annulus, or `f_stall` independent of dx) **→ stop.** The edge alone
  does not explain D2i, and the next packet must isolate the **stand-off feedback** instead —
  the D2h open-loop pair, whose `nozzle_collision_radius = 0.028` fix is already in its
  `run.py` and which was never rerun. Do not start the front on this evidence.

## Conditional run

`S05` (0.5 mm, ≈ 2–3 h) runs **only if** `w_edge` changes by more than 20 % from 2 mm to 1 mm,
i.e. only if there is something to converge. Two meshes give a difference; three give an order.
