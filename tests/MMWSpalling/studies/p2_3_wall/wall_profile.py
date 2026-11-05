#!/usr/bin/env python3
"""P2-3 (c1) -- the gas temperature up the bore, from the wall march's profile.

Reads `<plot_file>_wall_profile.csv` (time, k, z, n_faces, T_gas, P_k), written
by WallEnthalpyMarch at `jet_profile_interval`.

⚠⚠ EVERY COLUMN HERE IS PER-STEP STATE, SO IT IS AVERAGED OVER A WINDOW AND
NEVER SUMMED (standing rule 7(i)). `jet_P_wall` is the enthalpy the stream lost
to the wall ON THAT STEP, not a running total; summing it would multiply by the
step count. The csv rows are per-interval SNAPSHOTS, so the window this script
averages is stated in its output rather than left implicit.

⚠ The cooling length is fitted on ln(T_gas - T_ent) vs z by least squares, the
same form the integrator fits online into `wall_L_fit`; both are printed so the
offline and online fits can be compared rather than one trusted.

Usage: wall_profile.py <leg-dir> [t_lo t_hi]
"""
import math
import os
import sys

T_ENT = 293.15      # surface_patch.jet_T_ent on the Meier deck


def rows(leg_dir):
    hits = []
    for dp, _, fn in os.walk(leg_dir):
        hits += [os.path.join(dp, f) for f in fn if f.endswith("wall_profile.csv")]
    if not hits:
        return None
    lines = open(sorted(hits)[0]).read().splitlines()
    if len(lines) < 2:
        return []
    head = lines[0].split(",")
    out = []
    for ln in lines[1:]:
        f = ln.split(",")
        if len(f) != len(head):
            continue
        try:
            out.append({h: float(v) for h, v in zip(head, f)})
        except ValueError:
            continue
    return out


def fit_L(zs, Ts):
    """cooling length from ln(T - T_ent) = a - z/L, least squares."""
    xs, ys = [], []
    for z, T in zip(zs, Ts):
        d = T - T_ENT
        if d > 0.0:
            xs.append(z)
            ys.append(math.log(d))
    n = len(xs)
    if n < 2:
        return None
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    if abs(den) <= 0.0:
        return None
    slope = (n * sxy - sx * sy) / den
    return (-1.0 / slope) if slope < 0.0 else None


def report(leg_dir, t_lo=None, t_hi=None, label=""):
    rs = rows(leg_dir)
    if rs is None:
        print("  not computed -- no *_wall_profile.csv under %s" % leg_dir)
        return None
    if not rs:
        print("  not computed -- the profile csv is empty (march off, or no "
              "above-plane faces)")
        return None
    times = sorted({r["time"] for r in rs})
    if t_lo is None:
        t_lo = times[0]
    if t_hi is None:
        t_hi = times[-1]
    use = [t for t in times if t_lo <= t <= t_hi]
    print("  %s: %d snapshots, window t = %.1f-%.1f s (%d in window)"
          % (label or os.path.basename(leg_dir), len(times), t_lo, t_hi, len(use)))
    if not use:
        return None
    # the LAST snapshot in the window is the shape; the window mean is the level.
    last = sorted([r for r in rs if r["time"] == use[-1]], key=lambda x: x["k"])
    zs = [r["z"] for r in last]
    Ts = [r["T_gas"] for r in last]
    Pk = [r["P_k"] for r in last]
    nf = [r["n_faces"] for r in last]
    print("     bins %d, z = %.3f-%.3f m, faces %d-%d"
          % (len(last), min(zs), max(zs), int(min(nf)), int(max(nf))))
    print("     T_gas floor %.1f K -> mouth %.1f K  (drop %.1f K)"
          % (Ts[0], Ts[-1], Ts[0] - Ts[-1]))
    L = fit_L(zs, Ts)
    print("     cooling length L_fit (offline, ln(T-T_ent) vs z): %s"
          % ("%.4f m" % L if L else "not fitted (<2 cooling bins)"))
    # the window mean of the per-bin power, as a per-step state average
    tot = {}
    for r in rs:
        if t_lo <= r["time"] <= t_hi:
            tot[r["time"]] = tot.get(r["time"], 0.0) + r["P_k"]
    pw = sum(tot.values()) / len(tot) if tot else 0.0
    print("     wall power, MEAN over the window (per-step state, NOT summed): "
          "%.1f W = %.3f kW (quarter)" % (pw, pw / 1e3))
    print("     per-bin: z[m] n_faces T_gas[K] P_k[W]")
    step = max(1, len(last) // 12)
    for i in range(0, len(last), step):
        print("       %6.3f  %6d  %8.1f  %8.2f"
              % (zs[i], int(nf[i]), Ts[i], Pk[i]))
    return {"L_fit": L, "T_mouth": Ts[-1], "T_floor": Ts[0], "P_wall": pw,
            "bins": len(last)}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    d = sys.argv[1]
    lo = float(sys.argv[2]) if len(sys.argv) > 3 else None
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else None
    print("P2-3 (c1) -- the gas up the bore")
    print("=" * 72)
    report(d, lo, hi)
    return 0


if __name__ == "__main__":
    sys.exit(main())
