#!/usr/bin/env python3
"""Write and hash d2f_support/PREDICTIONS.md (D2f Goal 2), before any run.

Hand model: studies/d2c_steady/handmodel.py with the scored closures (J-M
1900 K, far power n = 1, momentum D_e, core 8, exhaust recirculation,
entrained mass, free surface aspect 1 / exponent 1), the ring moved to the
radius the quantile selects in an axisymmetric hole,
    R_Q(q) = sqrt(ri^2 + q (ro^2 - ri^2)),   ri = 28 mm, ro = 40 mm,
by rebinding feetmodel.R_Q at run time (the frozen files are not edited), and
the clip off (handmodel.transient(clear=False)). Transient: 0-300 s at 0.25 s
on the 2 mm outer grid, stopped if the centre reaches 0.36 m (the 0.40 m
domain's bottom rule). Refuses to overwrite an existing PREDICTIONS.md.
"""
import contextlib
import datetime
import hashlib
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D2C = os.path.join(os.path.dirname(HERE), "d2c_steady")
sys.path.insert(0, D2C)
import handmodel as hm  # noqa: E402

fm = hm.fm
OUTF = os.path.join(HERE, "PREDICTIONS.md")
T_NOZ = 1900.0
QS = (1.0, 0.9, 0.8)
R1_D2D = 1.65          # D2d R1 steady ROP (q 0.9 + clip, 0.70 m), m/h
FS = hm.FS(on=True, aspect=1.0, exponent=1.0)


def r_q(q):
    return math.sqrt(fm.RI ** 2 + q * (fm.RO ** 2 - fm.RI ** 2))


@contextlib.contextmanager
def ring_at(q):
    old = fm.R_Q
    fm.R_Q = r_q(q)
    try:
        yield
    finally:
        fm.R_Q = old


def fit(t, y, a, b):
    s = (t >= a) & (t <= b)
    return float(np.polyfit(t[s], y[s], 1)[0])


def evaluate(q, clear=False):
    cfg = hm.cfg_for(T_NOZ)
    with ring_at(q):
        st = hm.steady(T_NOZ, cfg)
        res = hm.transient(T_NOZ, cfg, fs=FS, t_end=300.0, z_stop=0.36, clear=clear)
        t, d = res["t"], res["d"]
        w25, w50 = hm.wall_profile(res, [0.025, 0.050])
        zg = np.linspace(0.0, d[-1], 401)
        rw = hm.wall_profile(res, zg)
    return dict(q=q, R_Q=r_q(q), v_steady=st["v"], sc_steady=st["s_c"], Trec_steady=float(st["T_rec"]),
                rop=fit(t, d, 150.0, t[-1]) * 3600.0, t_end=float(t[-1]),
                sc300=float(res["s_c"][-1]), Trec300=float(res["T_rec"][-1]),
                d25=2e3 * w25, d50=2e3 * w50, dmin=2e3 * float(rw.min()), feet=d[-1] * 1e3,
                V_mech=float(res["V_mech"][-1]) * 1e6)


def main():
    if os.path.exists(OUTF):
        sys.exit(f"{OUTF} exists; PREDICTIONS.md is written once")
    E = {q: evaluate(q) for q in QS}
    C = evaluate(0.9, clear=True)
    mean_rop = np.mean([E[q]["rop"] for q in QS])
    spread = max(abs(E[q]["rop"] / mean_rop - 1.0) for q in QS)
    per01 = (E[1.0]["rop"] - E[0.8]["rop"]) / 2.0
    per01_st = (E[1.0]["v_steady"] - E[0.8]["v_steady"]) / 2.0
    ii = "hold" if spread <= 0.05 else "fail"
    rows = "\n".join(
        f"| {q:.1f} | {E[q]['R_Q'] * 1e3:.2f} | {E[q]['v_steady']:.3f} | {E[q]['rop']:.3f} | "
        f"{E[q]['sc300'] * 1e3:.1f} / {E[q]['Trec300']:.0f} | {E[q]['sc_steady'] * 1e3:.0f} / "
        f"{E[q]['Trec_steady']:.0f} | {E[q]['d25']:.0f} / {E[q]['d50']:.0f} | {E[q]['dmin']:.1f} |"
        for q in QS)
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    L = f"""# D2f pre-registration: the contact support rule (pads, no clip), q = 1.0 / 0.9 / 0.8

Written {now} by `predict.py`, before any D2f run. **Frozen: never edit.**
Its sha256 and mtime are in RESULTS §0. These are diagnostics. No number
here or in the runs is scored against Meier's band.

## Hand model

- **Model:** `d2c_steady/handmodel.py` with the scored closures:
  - J-M at 1900 K;
  - far law power n = 1;
  - momentum D_e, core 8;
  - exhaust recirculation with entrained mass;
  - free surface with aspect 1, exponent 1.
- **The ring is at the radius the quantile selects in an axisymmetric hole,**
  R_Q(q) = √(r_i² + q(r_o² − r_i²)), set by rebinding `feetmodel.R_Q`.
- **The clip is off** (`transient(clear=False)`).
- **Steady:** `feetmodel.equilibria`.
- **Transient:** 0–300 s on the 2 mm outer grid, with the 0.36 m centre stop.
  The ROP is a fit of the feet depth over 150 s to the end.

| q | R_Q [mm] | steady ring ROP [m/h] | transient ROP 150–300 s [m/h] | s_c / T_rec at 300 s [mm / K] | steady s_c / T_rec | Ø 25 / 50 mm at 300 s | min Ø to the feet depth [mm] |
|---|---|---|---|---|---|---|---|
{rows}

- **q-sensitivity (hand):**
  - transient **{per01:+.3f} m/h per +0.1 of q** ({per01 / mean_rop * 100:+.1f} % of the mean);
  - steady {per01_st:+.3f} m/h per +0.1 of q.
- **Spread** of the three transient ROPs about their mean ({mean_rop:.3f} m/h):
  **±{spread * 100:.1f} %**.
- **Gate (ii) predicted to {ii.upper()}** (limit ±5 %).
- **Direction against D2d R1** (1.65 m/h steady, q 0.9 + the clip, 0.70 m):
  - the hand model with the clip on, at q 0.9, gives a transient
    {C['rop']:.3f} m/h, against {E[0.9]['rop']:.3f} m/h with it off;
  - in the hand model the clip does not touch the ring rate. It only cuts
    the annulus band, cutting {C['V_mech']:.1f} cm³ by 300 s;
  - so the hand model predicts **no change** in ROP from removing the clip.
    D2e's P-a at 2 mm ran +1–2 % faster than B0 (the simulation, not the
    hand model).
- **The minimum Ø to the feet depth is 2 R_Q in the hand model by
  construction** (the wall is floored at the ring). It carries no
  information beyond R_Q.

## Limits (stated)

- **The hand model is axisymmetric with no column scatter.** It cannot
  predict the mesh behaviour of the per-pad maximum, which is set by the
  single slowest column. **Gate (i) has no hand prediction.** The expectation
  below comes from D2e's diagnosis, not from the hand model.
- **Its q-dependence is geometric only:** a larger R_Q gives a lower h and
  cooler gas at the ring. In the simulation, q also picks which columns'
  scatter sets the pad. The hand spread is therefore a lower bound on the
  simulated q-sensitivity.
- **The hole numbers are direction only** (packet: the wellhead, two phases,
  and Fig. 8.8's lower-bound widths).

## Expectations and the decision rule (the packet's, fixed before any run)

| id | expectation | criterion |
|---|---|---|
| E1 (gate i) | Q10 converges | Q10 burner ROP, 1 mm / 2 mm over 150–300 s, within **5 %**, and every 50 s block after 50 s within 5 %. No hand number: the maximum rides on the slowest column. D2e's P-a (q 0.9) converged to +0.1 %, and whether q 1.0 keeps that is the open question |
| E2 (gate ii) | q does not control the rate | the 2 mm burner ROP over 150–300 s for q 0.8, 0.9 and 1.0 is within **±5 %** of their mean. Hand: ±{spread * 100:.1f} %, so predicted to {ii} |
| E3 | the ordering | the ROP falls as q rises (a larger R_Q gives a lower h and cooler gas). Hand: {E[0.8]['rop']:.3f} > {E[0.9]['rop']:.3f} > {E[1.0]['rop']:.3f} m/h |
| E4 | the transient ROP, 2 mm, 150–300 s, per q | the hand value ±27 % (the D2c P6 rate band): q 1.0 {E[1.0]['rop']:.2f}, q 0.9 {E[0.9]['rop']:.2f}, q 0.8 {E[0.8]['rop']:.2f} m/h |
| E5 | s_c at 300 s, 2 mm, per q | the hand value ±20 mm (the D2c P6 length band): {E[1.0]['sc300'] * 1e3:.0f} / {E[0.9]['sc300'] * 1e3:.0f} / {E[0.8]['sc300'] * 1e3:.0f} mm |
| E6 | Q09_2mm is P-a's configuration | thermo and CSV rows to t ≤ 200 s byte-identical to `d2e_mesh/output/Pa_2mm` (harness) |

**Outcomes:**
- **(a) E1 and E2 hold:** D2g scores **q = 1.0**, the literal "highest rock
  under each foot". The jolt allowance is shown to be immaterial.
- **(b) E2 holds, E1 fails:** D2g scores **q = 0.9**, justified as the jolt
  allowance and shown converged by P-a (and by Q09_1mm to 300 s, if run).
  - The report covers the maximum's mesh failure and its mechanism, from
    `rim.py`.
  - If Q09_1mm was not run, the adoption is conditional on running it in
    D2g.
- **(c) E2 fails:** q controls the rate. The report gives the ROP per 0.1 of
  q and stops. D2g then needs a physical basis for the allowance. **q is not
  picked to hit the band.**
"""
    with open(OUTF, "w") as f:
        f.write(L)
    sha = hashlib.sha256(open(OUTF, "rb").read()).hexdigest()
    mt = datetime.datetime.fromtimestamp(os.path.getmtime(OUTF)).strftime("%Y-%m-%d %H:%M:%S")
    print(f"PREDICTIONS.md sha256 {sha}\nmtime {mt}")


if __name__ == "__main__":
    main()
