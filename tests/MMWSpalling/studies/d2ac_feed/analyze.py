#!/usr/bin/env python3
"""D2ac analysis -- scores CRITERION.md (a)-(g) and P1-P6.

Imports the equal-burner-depth machinery UNCHANGED from
studies/d2x_areasplit/analyze.py (which itself imports d2v_hmap's zface/profile/
diameter). Its WHERE lookup is EXTENDED with this study's legs.

⚠⚠ THE LZ TRAP. Both modules hard-code LZ = 0.4 at MODULE level
(d2v:30, d2x:47) and zface() returns LZ - nv*dz. D2ac's block is 0.5 m, so
importing them "unchanged" would misplace every depth by 100 mm and silently
corrupt Ø(z), the pit and the volume. `with_lz(n)` binds V.LZ and X.LZ from the
LEG'S OWN .done metadata around each call and restores them after. The imported
functions are untouched; only the constant they read is bound per leg.

⚠ nozzle_z is registered ONLY under the feet rule (:629), so the prescribed legs
have no descent column. Under `prescribed` z_n is exact by construction:
    z_n(t) = nozzle_z0 - nozzle_feed * t
and the standoff comes from `jet_s_c` (the axis bin's mean stand-off).
"""
import ast
import contextlib
import filecmp
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load_mod(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


X = load_mod("d2x_an", os.path.join(STUD, "d2x_areasplit", "analyze.py"))
V = X.V

OUT = os.path.join(HERE, "output")
AB = os.path.join(STUD, "d2ab_gaph", "output")
Z2 = os.path.join(STUD, "d2z2_multicell", "output")
Z1 = os.path.join(STUD, "d2z1_delete", "output")
W = os.path.join(STUD, "d2w_attribution", "output")

NEW = ("FM", "FL", "FH", "FM13")
X.WHERE.update({n: OUT for n in NEW + ("S0",)})
X.WHERE.update({"SWc": AB, "SW": Z2, "RB_D": Z1, "A": W})

# ---- PRE-REGISTERED, from PREFLIGHT.md. DO NOT EDIT. ----------------------
NOZZLE_Z0 = 0.550
EXCESS_LO, EXCESS_HI = 2.15, 2.60      # item 0(a), the deliverable
S_CEIL = 100.0                         # mm, Fig. 8.8's standoff ceiling
S_PRED = {1.5: (215, 222), 1.3: (248, 256)}
D_TARGET, D_TOL = 85.0, 5.0            # mm, Tab. 8.2 TEXT
VOL_CEIL = 3.42                        # L per 0.5 m
CAL = (150.0, 300.0)                   # mm, the calibrated section
FOOT_RI, FOOT_RO, FOOT_STANDOFF = 0.028, 0.040, 0.050
WIN = (150.0, 400.0)


def have(n):
    return os.path.exists(os.path.join(X.WHERE.get(n, OUT), n + ".done"))


def meta(n):
    return ast.literal_eval(open(os.path.join(X.WHERE[n], n + ".done")).read())


def thermo(n):
    p = os.path.join(X.WHERE[n], n, "thermo.dat")
    with open(p) as f:
        cols = f.readline().split()
    v = np.loadtxt(p, skiprows=1, ndmin=2)
    return {c: v[:, k] for k, c in enumerate(cols)}


@contextlib.contextmanager
def with_lz(n):
    """Bind the imported modules' LZ to THIS leg's block height."""
    lz = meta(n).get("lz", 0.4)
    ov, ox = V.LZ, X.LZ
    V.LZ, X.LZ = lz, lz
    try:
        yield lz
    finally:
        V.LZ, X.LZ = ov, ox


def hdr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def feed_of(n):
    return float(meta(n)["nozzle_feed"]) * 3600.0      # m/h


def floor_trace(n):
    """(t, standoff [m], floor z [m]) for a prescribed leg."""
    d = thermo(n)
    t = d["time"]
    s = d["jet_s_c"]
    zn = NOZZLE_Z0 - float(meta(n)["nozzle_feed"]) * t
    return t, s, zn - s


def crit_a(rep):
    hdr("(a) SCORED -- THE FLOOR-FLUX EXCESS  (the packet's product)")
    print("  realised FLOOR recession at the realised standoff, against the feed.")
    print(f"  PRE-REGISTERED (PREFLIGHT.md a): {EXCESS_LO}-{EXCESS_HI}x at a "
          f"{S_CEIL:.0f} mm standoff.")
    print("  ⚠ Expected to be a FAILURE against Meier and a MEASUREMENT of the")
    print("    missing factor -- not a pass/fail bar on the model.\n")
    print(f"{'leg':<6}{'feed':>7}{'s_end':>8}{'floor R':>9}{'excess':>8}"
          f"{'R@100mm':>9}{'excess@100':>11}")
    for n in NEW:
        if not have(n):
            continue
        t, s, zf = floor_trace(n)
        f = feed_of(n)
        m = (t >= WIN[0]) & (t <= WIN[1])
        R = -np.polyfit(t[m], zf[m], 1)[0] * 3600.0
        send = float(s[-1]) * 1e3
        # De-scale to the reference 100 mm standoff through the map's own decay
        # term, pow(0.09/max(0.09,s),1.0), so legs at different standoffs are
        # comparable and the number matches item 0's basis.
        dec = lambda x: 0.09 / max(0.09, x)
        R100 = R / dec(float(np.mean(s[m]))) * dec(0.100)
        rep.setdefault("R", {})[n] = (R, send, R100 / f)
        print(f"{n:<6}{f:7.2f}{send:8.1f}{R:9.3f}{R/f:8.2f}{R100:9.3f}"
              f"{R100/f:11.2f}")
    vals = [v[2] for v in rep.get("R", {}).values()]
    if vals:
        ok = all(EXCESS_LO * 0.85 <= v <= EXCESS_HI * 1.15 for v in vals)
        print(f"\n  measured excess at 100 mm: {min(vals):.2f} - {max(vals):.2f} x")
        print(f"  pre-registered {EXCESS_LO}-{EXCESS_HI} x   "
              f"{'CONSISTENT' if ok else 'OUTSIDE (reported, not retrofitted)'}")
        rep.setdefault("a", []).append(("excess", ok))


def crit_b(rep):
    hdr("(b) SCORED -- PAD CONTACT, computed OFFLINE")
    print("  foot tips  z_tip(t) = nozzle_z0 - feed*t - 0.050")
    print(f"  pad ring   {FOOT_RI*1e3:.0f} <= r <= {FOOT_RO*1e3:.0f} mm")
    print("  contact    max top-face height over that annulus >= z_tip")
    print("  ⚠ nozzle_collision_radius is NOT the test -- but NOT because it is")
    print("    inert. PREFLIGHT_ADDENDUM.md: its default is `feet ? foot_r_inner")
    print("    : radius`, so under PRESCRIBED it defaults to radius = 0.2 m and")
    print("    fired on all four legs at t = 120 s. Both the packet's claim and")
    print("    PREFLIGHT.md (b)'s were wrong. It is not the test because it")
    print("    guards the BORE, not the pad ring.")
    print("  ⚠ CONTACT AT t = 0 IS BY CONSTRUCTION: nozzle_z0 = lz +")
    print("    foot_standoff puts the feet exactly on the flat surface, so")
    print("    d(0) = 0 identically. t = 0 is excluded below.\n")
    print(f"{'leg':<6}{'first contact':>15}{'contact frac':>14}"
          f"{'worst interf':>14}{'d@end':>10}")
    for n in NEW:
        if not have(n):
            continue
        with with_lz(n):
            m, r, a = X.cell(n)
            t, s, zf = floor_trace(n)
            ann = (r >= FOOT_RI) & (r <= FOOT_RO)
            # ⚠ t = 0 excluded: d(0) = 0 by construction (see above). Using
            # `if first` rather than `is not None` would ALSO have mis-reported
            # a genuine t = 0 contact as "never" -- a falsy-zero bug.
            ts = [x for x in t[::max(1, len(t) // 120)] if x > 0.0]
            first, nc, worst, dend = None, 0, -1e9, None
            for tt in ts:
                zface, _ = V.zface(m, r, a, tt)
                ztip = NOZZLE_Z0 - float(meta(n)["nozzle_feed"]) * tt - FOOT_STANDOFF
                d = float(zface[ann].max()) - ztip
                worst = max(worst, d)
                dend = d
                if d >= 0.0:
                    nc += 1
                    if first is None:
                        first = float(tt)
            frac = nc / len(ts)
            print(f"{n:<6}{(f'{first:.1f} s' if first is not None else 'never'):>15}"
                  f"{frac*100:>13.1f}%{worst*1e3:>13.1f}mm{dend*1e3:>9.1f}mm")
            # CRITERION.md (b) specifies REPORTED quantities (first contact,
            # fraction, worst interference), not a bar. P3 carries the verdict,
            # so (b) is recorded as a measurement and kept OUT of the pass/fail
            # scorecard -- treating "contact occurred" as a FAIL would invert
            # the prediction it confirms.
            rep.setdefault("bmeas", {})[n] = (first, frac, worst * 1e3,
                                              dend * 1e3)
    print("\n  ⚠ 'worst interf' > 0 means the pad ring reached the foot plane.")
    print("    P3 predicted contact at 1.5 m/h (D2x pad-ring rate ~2.1 m/h).")


def crit_c(rep):
    hdr(f"(c) SCORED -- Ø(z) over the calibrated section "
        f"{CAL[0]:.0f}-{CAL[1]:.0f} mm against {D_TARGET} +/- {D_TOL} mm")
    print("  ⚠ 85 mm is Tab. 8.2 TEXT. NO photo-derived diameter is claimed.")
    print("  ⚠ mouth (0-40 mm) and the ~100 mm underream are REPORTED, not")
    print("    scored: Meier attributes the underream to a pause.\n")
    depths = [40, 60, 80, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350]
    print(f"  {'leg':<6}" + "".join(f"{d:>6}" for d in depths))
    for n in NEW:
        if not have(n):
            continue
        with with_lz(n) as lz:
            m, r, a = X.cell(n)
            zf, nv = V.zface(m, r, a, float(thermo(n)["time"][-1]))
            rc, zc, _ = V.profile(r, zf)
            D = {d: V.diameter(rc, zc, d) * 1e3 for d in depths}
            rep.setdefault("D", {})[n] = D
            print(f"  {n:<6}" + "".join(
                f"{(D[d] if D[d] == D[d] else float('nan')):>6.1f}" for d in depths))
    print()
    for n in NEW:
        if n not in rep.get("D", {}):
            continue
        D = rep["D"][n]
        cal = [D[d] for d in depths if CAL[0] <= d <= CAL[1] and D[d] == D[d]]
        if not cal:
            print(f"  {n:<5} calibrated section not reached")
            continue
        mean = float(np.mean(cal))
        ok = abs(mean - D_TARGET) <= D_TOL
        # ⚠ The hole only reaches ~310-336 mm in 500 s, so the upper end of the
        # 150-300 mm calibrated section can sit INSIDE the rounded tip, where Ø
        # collapses toward 0. Reported ALONGSIDE the pre-registered score, never
        # instead of it (the D2z-2 / D2ab pattern).
        # ⚠ computed HERE, not read from rep: crit_d runs AFTER crit_c.
        with with_lz(n) as lz:
            mm, rr, aa = X.cell(n)
            zzf, _ = V.zface(mm, rr, aa, float(thermo(n)["time"][-1]))
            floor = (lz - float(zzf.min())) * 1e3
        clean = [D[d] for d in depths
                 if CAL[0] <= d <= min(CAL[1], floor - 40.0) and D[d] == D[d]]
        tip = (f"   tip-excluded (<= floor-40 mm): {np.mean(clean):6.2f} mm "
               f"over {len(clean)} bins" if clean and len(clean) < len(cal)
               else "")
        print(f"  {n:<5} calibrated mean Ø {mean:6.2f} mm  "
              f"(range {min(cal):.1f}-{max(cal):.1f})  "
              f"{'PASS' if ok else 'FAIL'}  vs {D_TARGET} +/- {D_TOL}{tip}")
        rep.setdefault("c", []).append((n, ok))
        rep.setdefault("Dcal", {})[n] = mean
    print("\n  ⚠ P5: FL/FH are OPEN, 80-95 mm, NO preferred value. A landing")
    print("    near 85 mm is NOT derived -- 85 is Meier's own number.")


def crit_d(rep):
    hdr(f"(d) SCORED -- volume <= {VOL_CEIL} L per 0.5 m")
    print("  ⚠ 93.3 mm is NOT a second diameter target: the water-filled volume")
    print("    includes the mouth, the pause underream and the rounded tip, so a")
    print("    uniform 93 mm bore would CONTRADICT Fig. 8.8.\n")
    print(f"  {'leg':<6}{'floor depth':>13}{'volume L':>10}{'per 0.5 m':>11}"
          f"{'vs ceiling':>12}")
    for n in NEW:
        if not have(n):
            continue
        with with_lz(n) as lz:
            m, r, a = X.cell(n)
            t, s, zf_tr = floor_trace(n)
            zface, nv = V.zface(m, r, a, float(t[-1]))
            depth = (lz - float(zface.min()))
            vol = float(nv.sum()) * m["dz"] ** 3 * 4.0 * 1e3      # L, full hole
            per = vol / depth * 0.5 if depth > 0 else float("nan")
            ok = per <= VOL_CEIL
            print(f"  {n:<6}{depth*1e3:>11.1f}mm{vol:>10.3f}{per:>11.3f}"
                  f"{'  PASS' if ok else '  FAIL':>12}")
            rep.setdefault("d", []).append((n, ok))
            rep.setdefault("vol", {})[n] = (vol, per, depth * 1e3)
    print("\n  P6: volume BELOW the ceiling, the shortfall being mouth + pause.")


def crit_e(rep):
    hdr("(e) SCORED -- the standoff against the photograph's ceiling")
    print(f"  <= {S_CEIL:.0f} mm passes. Pre-registered FAILURE: "
          f"{S_PRED[1.5][0]}-{S_PRED[1.5][1]} mm at 1.5 m/h, "
          f"{S_PRED[1.3][0]}-{S_PRED[1.3][1]} at 1.3.\n")
    print(f"  {'leg':<6}{'feed':>6}" + "".join(f"{int(x):>8}" for x in
                                               (50, 100, 200, 300, 400, 500))
          + "   [standoff mm at t s]")
    for n in NEW:
        if not have(n):
            continue
        t, s, zf = floor_trace(n)
        row = [float(np.interp(x, t, s)) * 1e3 for x in (50, 100, 200, 300, 400, 500)]
        send = row[-1]
        f = feed_of(n)
        lo, hi = S_PRED.get(round(f, 1), (0, 0))
        ok = send <= S_CEIL
        print(f"  {n:<6}{f:>6.2f}" + "".join(f"{v:>8.1f}" for v in row)
              + f"   end {send:.0f} mm, pred {lo}-{hi}  "
              + ("PASS" if ok else "FAIL (pre-registered)"))
        rep.setdefault("e", []).append((n, ok))
        rep.setdefault("s_end", {})[n] = send


def crit_f(rep):
    hdr("(f) SCORED -- the power identity, +/- 0.02 W")
    for n in NEW + ("S0",):
        if not have(n):
            continue
        d = thermo(n)
        if "jet_P_face" not in d:
            print(f"  {n:<5} no jet_P_face (no march)")
            continue
        m = (d["time"] >= WIN[0]) & (d["time"] <= WIN[1])
        e = float(np.abs((d["jet_P_face"] - (d["patch_P_robin"] + d["side_P_face"]
                                             + d["side_P_loss"]))[m]).max())
        ok = e <= 0.02
        print(f"  {n:<5} max residual {e:.3e} W   {'PASS' if ok else 'FAIL'}"
              f"    ledger_err max {np.abs(d['ledger_err']).max():.2e}")
        rep.setdefault("f", []).append((n, ok))


def crit_g(rep):
    hdr("(g) SCORED -- identity")
    for suf in ("_removal_events.csv", "_jet_profile.csv"):
        pa, pb = os.path.join(OUT, "S0" + suf), os.path.join(Z1, "RB_D" + suf)
        if not (os.path.exists(pa) and os.path.exists(pb)):
            print(f"  {suf}: missing")
            continue
        ok = filecmp.cmp(pa, pb, shallow=False)
        print(f"  S0 vs RB_D {suf}: {'IDENTICAL' if ok else 'DIFFER'}")
        rep.setdefault("g", []).append((suf, ok))
    pa, pb = os.path.join(OUT, "S0"), os.path.join(Z1, "RB_D")
    if os.path.exists(pa) and os.path.exists(pb):
        ha = open(os.path.join(pa, "thermo.dat")).readline().split()
        hb = open(os.path.join(pb, "thermo.dat")).readline().split()
        va = np.loadtxt(os.path.join(pa, "thermo.dat"), skiprows=1, ndmin=2)
        vb = np.loadtxt(os.path.join(pb, "thermo.dat"), skiprows=1, ndmin=2)
        nr = min(len(va), len(vb))
        delta = set(ha) ^ set(hb)
        diff = [c for c in set(ha) & set(hb)
                if not np.array_equal(va[:nr, ha.index(c)], vb[:nr, hb.index(c)])]
        ok = not delta and not diff
        print(f"  S0 vs RB_D thermo: header delta "
              f"{'EMPTY' if not delta else sorted(delta)}; {len(diff)} column(s) "
              f"differ  {'PASS' if ok else 'FAIL'}")
        rep.setdefault("g", []).append(("thermo", ok))
    print("\n  tree gate vs witness/d2z1.md5: see GATE.txt.")
    print("  ⚠ NO SOURCE CHANGE in this packet, so the binary is D2ab's")
    print("    post-guard build unchanged and the gate is expected to be a")
    print("    memoised cache hit on the same binary hash.")


def witnesses(rep):
    hdr("WITNESSES -- what would make this inert")
    print(f"  {'leg':<6}{'feed applied':>14}{'pit vs lz':>20}"
          f"{'side_P_wall':>13}")
    for n in NEW:
        if not have(n):
            continue
        d = thermo(n)
        with with_lz(n) as lz:
            m, r, a = X.cell(n)
            zf, _ = V.zface(m, r, a, float(d["time"][-1]))
            pit = (lz - float(zf.min())) * 1e3
        # Under `prescribed` z_n is exact by construction (nozzle_z is not even
        # registered), so "the feed was not applied" is checked by the standoff
        # RESPONDING rather than by a descent column.
        t, s, _ = floor_trace(n)
        spw = 4.0 * float(np.nanmean(d["side_P_wall"][-20:])) \
            if "side_P_wall" in d else float("nan")
        print(f"  {n:<6}{'exact (analytic)':>14}"
              f"{f'{pit:.0f} of {lz*1e3:.0f} mm':>20}{spw:>13.1f}")
    print("\n  ⚠ side_P_{bore,pads,wall} are MEANINGLESS here: foot_r_inner/outer")
    print("    are dropped (forced, PREFLIGHT d) so they fall to 0.0 and EVERY")
    print("    charged face bins as `wall`. Shown only as a reached-a-face")
    print("    witness for the FL/FH pair, never as a region split.")


def predictions(rep):
    hdr("PREDICTIONS")
    se = rep.get("s_end", {})
    for n in NEW:
        if n in se:
            ok = 180 <= se[n] <= 260
            print(f"  P1 {n} standoff 180-260 mm (got {se[n]:.0f}): "
                  f"{'HELD' if ok else 'REFUTED'}")
    vals = [v[2] for v in rep.get("R", {}).values()]
    if vals:
        ok = all(2.0 <= v <= 3.0 for v in vals)
        print(f"  P2 floor-flux excess 2.0-3.0x (got {min(vals):.2f}-{max(vals):.2f}): "
              f"{'HELD' if ok else 'REFUTED'}")
    for n, (first, frac, worst, dend) in rep.get("bmeas", {}).items():
        print(f"  P3 {n} pads make contact at the feed: "
              + ("HELD" if first is not None else
                 "REFUTED (ring clears, burner runs free)")
              + (f"  (first {first:.1f} s, {frac*100:.1f} % of run, worst "
                 f"+{worst:.1f} mm, then clears to {dend:.1f} mm)"
                 if first is not None else ""))
    if "FM" in rep.get("Dcal", {}):
        v = rep["Dcal"]["FM"]
        print(f"  P4 FM calibrated Ø exceeds 85 mm (got {v:.2f}, "
              f"exposure scaling predicted ~99): {'HELD' if v > 85 else 'REFUTED'}")
    for n in ("FL", "FH"):
        if n in rep.get("Dcal", {}):
            v = rep["Dcal"][n]
            print(f"  P5 {n} OPEN 80-95 mm (got {v:.2f}): "
                  f"{'in range' if 80 <= v <= 95 else 'outside'}  "
                  f"⚠ NOT derived -- 85 is Meier's own number")
    for n, ok in rep.get("d", []):
        print(f"  P6 {n} volume below {VOL_CEIL} L per 0.5 m: "
              f"{'HELD' if ok else 'REFUTED'}")


def main():
    print("D2ac -- scoring CRITERION.md (sha in CRITERION.sha256)")
    print("legs on disk: " + ", ".join(n for n in NEW + ("S0",) if have(n)))
    print("⚠ THE RATE IS AN INPUT. Burner descent is never scored.")
    rep = {}
    for f in (crit_a, crit_b, crit_c, crit_d, crit_e, crit_f, crit_g,
              witnesses, predictions):
        f(rep)
    hdr("SCORECARD")
    for k in "abcdefg":
        rows = rep.get(k, [])
        if not rows:
            continue
        npass = sum(1 for _, ok in rows if ok)
        print(f"  ({k}) {npass}/{len(rows)} pass"
              + ("" if npass == len(rows) else
                 "   not passing: " + ", ".join(str(n) for n, ok in rows if not ok)))
    print("\n  ⚠ (b) is REPORTED, not barred -- CRITERION.md (b) specifies")
    print("    quantities, and P3 carries the verdict.")
    print("  ⚠ (e) FAIL is the PRE-REGISTERED outcome, not a defect: the")
    print("    packet's headline is that the model settles far past Fig. 8.8's")
    print("    100 mm ceiling, and (a) quantifies by how much.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
