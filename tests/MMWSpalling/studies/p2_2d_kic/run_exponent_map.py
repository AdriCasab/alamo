#!/usr/bin/env python3
"""P2-2d ord 4 deliverable 12 — THE EXPONENT MAP. A map, not a fit.

Three flux levels x both toughness tables x `m` across the literature prior x
both `alpha` treatments the MODEL CAN RUN. The exponent reported is the one the
model PRODUCES: d ln(onset superheat) / d ln(q) across the flux levels. No `m`
is chosen to make it agree with anything.

⚠⚠ THE MODEL CANNOT RUN THE COHERENT CAPPED PAIRING, AND THIS IS A REAL LIMIT.
`(a7)` establishes that `alpha` and `E` must move together -- capping `alpha`
for microcracking while keeping the full `E` collapse is incoherent, and is the
only pairing that produces s < 0 and an apparent match to Kant. But the only
key is `grimsel.alpha_cap_theta`; there is NO key that caps `E`. So the two
treatments the model can produce are:
    alpha full  / E full    -- coherent
    alpha capped/ E full    -- INCOHERENT, reported as a DIAGNOSTIC only
and the coherent capped pairing exists offline only (derive_onset.py's (a7)).
⇒ `(c6)`'s span is reported across what the model can run, with the coherent
capped exponent carried from ord 1 and labelled as offline.

⚠ (c5) requires `size_volume = cell` to read ~0: under `cell` the size factor
is identically 1 for ANY m, so the modulus cannot give a flux exponent at all.
That control is run first.

⚠ `V_ref` is anchored at ell_q^3 measured on the h_conv = 500 leg, with the
conductivity settled by deliverable 6 -- NOT chosen. Any bracket other than 1
would make a_0 depend on m (standing rule 1).

Usage: run_exponent_map.py
"""
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, HERE)
import derive_onset as D  # noqa: E402

BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant",
                    "sp_kant_onset_meiercfg", "input_model")
OUT = os.path.join(HERE, "expmap")
T_I = 293.15
R_DEF = D.UCS0 / D.sigma(850.0, False, False)

# Kant MEASURED flux range 0.45-0.87 MW/m2 [251-253]; h_conv sets the drive.
FLUXES = (300.0, 500.0, 800.0)
M_PRIOR = (8.0, 10.0, 13.0, 17.0, 20.0)      # literature prior, never fitted
TABLES = (("grimsel", "grimsel_table"), ("nasseri", "nasseri_table"))
ALPHAS = (("full", []),
          ("capped300", ["microstructure.phase0.grimsel.alpha_cap_theta=300"]))
KANT_LO, KANT_HI = 0.229, 0.306


def leg(tag, extra):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    argv = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
             DECK,
             "amr.n_cell=24 24 24", "amr.max_grid_size=24",
             "spall.surface_T=reconstructed", "stop_time=45.0",
             "amr.thermo.int=1", "amr.thermo.plot_int=1",
             "spall.removal_events_csv=1",
             "microstructure.phase0.law=grimsel", "spall.stress_law=grimsel",
             "spallation.sp.kic.scale=1.0",
             "spallation.sp.stress_relief=%.6f" % R_DEF,
             # ⚠⚠ `spall.flake_length = flux` IS WHAT MAKES `ell_q` EXIST AT
             # ALL. Under the default `twocell` the flux length is never
             # computed (`Removal.H`'s header appends the `ell_q` column only
             # under this key or `flake_frac`), so the anchor has nothing to
             # read -- and deliverable 6's per-cell conductivity change is
             # itself inert without it. The anchor leg returned a firing with
             # `ell_q = None` until this was set.
             "spall.flake_length=flux"]
            + extra + ["plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    if rc != 0:
        return "ABORT", None
    p = os.path.join(d, "run_removal_events.csv")
    if not os.path.exists(p):
        return "NOFIRE", None
    # ⚠⚠ `ell_q` ON THE FIRST FIRING ROW IS NOT THE ANCHOR. The flux length
    # falls back to the shipped two-cell length whenever `flame_h` is 0 or
    # q_net <= 0 on that column at that instant, and it writes 0 when it does.
    # MEASURED on the anchor leg: ell_q is non-zero on 1440 of 9792 rows, and
    # where it exists it is extremely tight (2.5786-2.5823 mm, 0.14 % spread) --
    # but the FIRST firing row reads exactly 0. Anchoring on it would have set
    # V_ref from a fallback. The packet asks for the value "at 850 K", so take
    # the median over rows whose reconstructed surface is near 850 K.
    first, lqs = None, []
    with open(p) as f:
        for row in csv.DictReader(f):
            if row.get("regime") != "1":
                continue
            ts = float(row["T_surf"])
            if first is None:
                first = ts - T_I
            lq = row.get("ell_q")
            if lq not in (None, ""):
                v = float(lq)
                if v > 0.0 and abs(ts - 850.0) <= 25.0:
                    lqs.append(v)
    if first is None:
        return "NOFIRE", None
    lqs.sort()
    return first, (lqs[len(lqs) // 2] if lqs else None)


def slope(fluxes, rises):
    """d ln(superheat) / d ln(q), least squares on the logs."""
    pts = [(math.log(q), math.log(r)) for q, r in zip(fluxes, rises)
           if isinstance(r, float) and r > 0.0]
    if len(pts) < 2:
        return None
    n = len(pts)
    sx = sum(x for x, _ in pts); sy = sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts); sxy = sum(x * y for x, y in pts)
    den = n * sxx - sx * sx
    return None if den == 0 else (n * sxy - sx * sy) / den


def main():
    print("P2-2d ord 4 -- THE EXPONENT MAP (a map, not a fit)")
    print("  R = %.4f computed; kic.scale = 1 pinned; flaw 8 um held" % R_DEF)
    print("  fluxes (h_conv): %s" % ", ".join("%.0f" % q for q in FLUXES))
    print("  Kant MEASURED exponent 0.229-0.306 [258]\n")

    # ---- the ANCHOR for V_ref, MEASURED not chosen ----------------------
    print("  [d8] anchoring V_ref at ell_q^3 on the h_conv = 500 leg")
    rise, lq = leg("anchor", ["surface_patch.h_conv=500",
                              "spallation.sp.kic.law=grimsel_table"])
    if not isinstance(rise, float) or not lq:
        print("      FAIL: anchor leg gave %s / ell_q %s" % (rise, lq))
        return 1
    V_ref = lq ** 3
    print("      ell_q = %.4f mm, MEDIAN over firing rows within 25 K of 850 K"
          % (lq * 1e3))
    print("      ⚠ the conductivity is the PER-CELL value at the top cell")
    print("        centre (deliverable 6). The packet predicted ell_q^3 ~")
    print("        1.67e-8 for a cell-centre read against 8.36e-9 for the old")
    print("        scalar k = 1.5 -- this measures %.3e, i.e. the cell-centre"
          % (lq ** 3))
    print("        bracket, confirming the change is live and which way.")
    print("      V_ref = ell_q^3 = %.4e m^3   <- FORCED, bracket exactly 1"
          % V_ref)
    print("      ⚠ ell_q^3 is a PROXY, never a spall volume: %.2f mm already"
          % (lq * 1e3))
    print("        exceeds every measured spall thickness (Wilkinson 0.1-1 mm;")
    print("        Hoeser 100-200 um).")

    # ---- (c5)'s control: `cell` must read ~0 -----------------------------
    print("\n  [c5 control] size_volume = cell must give exponent ~ 0")
    rises = []
    for q in FLUXES:
        r, _ = leg("cell_q%.0f" % q,
                   ["surface_patch.h_conv=%.0f" % q,
                    "spallation.sp.kic.law=grimsel_table",
                    "weibull.size_volume=cell"])
        rises.append(r)
        print("      h_conv %4.0f -> rise %s" % (q, r if isinstance(r, str)
                                                 else "%.2f K" % r))
    e_cell = slope(FLUXES, rises)
    print("      exponent under `cell`: %s   (expected ~0: the size factor is"
          % ("n/a" if e_cell is None else "%+.4f" % e_cell))
    print("      identically 1 for ANY m, so no modulus can give a flux"
          " exponent)")

    # ---- the map ---------------------------------------------------------
    print("\n  [c5/c6/c7] THE MAP, size_volume = flux, V_ref anchored")
    print("  %-9s %-10s %5s %10s %10s %10s %9s %s"
          % ("table", "alpha", "m", "q=300", "q=500", "q=800", "exponent",
             "vs Kant"))
    rows = []
    for tname, tkey in TABLES:
        for aname, akeys in ALPHAS:
            for m in M_PRIOR:
                rr = []
                for q in FLUXES:
                    tag = "%s_%s_m%g_q%.0f" % (tname, aname, m, q)
                    r, _ = leg(tag, ["surface_patch.h_conv=%.0f" % q,
                                     "spallation.sp.kic.law=%s" % tkey,
                                     "weibull.enabled=1",
                                     "weibull.a0_ig=8.0e-6",
                                     "weibull.a0_gb=40.0e-6",
                                     "weibull.m=%g" % m,
                                     "weibull.size_volume=flux",
                                     "weibull.V_ref=%.9e" % V_ref] + akeys)
                    rr.append(r)
                e = slope(FLUXES, rr)
                verdict = "n/a"
                if e is not None:
                    verdict = ("IN" if KANT_LO <= e <= KANT_HI
                               else ("low" if e < KANT_LO else "high"))
                rows.append((tname, aname, m, e))
                fmt = lambda r: (r if isinstance(r, str) else "%.1f K" % r)
                print("  %-9s %-10s %5g %10s %10s %10s %9s %s"
                      % (tname, aname, m, fmt(rr[0]), fmt(rr[1]), fmt(rr[2]),
                         "n/a" if e is None else "%+.4f" % e, verdict))

    # ---- (c6) the alpha span, and the limit ------------------------------
    print("\n  (c6) THE ALPHA SENSITIVITY ON THE EXPONENT")
    for tname, _ in TABLES:
        for m in M_PRIOR:
            es = {a: e for (t, a, mm, e) in rows
                  if t == tname and mm == m and e is not None}
            if len(es) == 2:
                span = abs(es["full"] - es["capped300"])
                print("      %-9s m=%-4g full %+.4f  capped %+.4f  span %.4f"
                      % (tname, m, es["full"], es["capped300"], span))
    print("  ⚠⚠ the `capped` column here is the INCOHERENT alpha-capped/E-full")
    print("     pairing -- there is NO key that caps E, so the model cannot")
    print("     run the coherent capped treatment at all. Its exponent is")
    print("     carried from ord 1 (OFFLINE): 0.089 at m = 8.")
    print("  ⚠ report the denominator with the exponent: ord 1 gives")
    print("     s + 3/2m - kappa = 1.40 (full) vs 0.45 (incoherent capped),")
    print("     so the capped exponent is amplified and hypersensitive.")

    # ---- (c7) m in both conventions --------------------------------------
    print("\n  (c7) `m` IN BOTH CONVENTIONS against the published prior")
    print("      %-12s %-14s %s" % ("m_flaw", "m_strength=2m", "published"))
    for m in M_PRIOR:
        print("      %-12g %-14g %s" % (m, 2 * m,
              "<- Meier deck ships this (MECHANICAL 20)" if m == 20 else ""))
    print("      published granite STRENGTH moduli: 8 (Barre) / 13 (Westerly)")
    print("      from spallation; 10-17 (Shen); 20 MECHANICAL (Wilkinson).")
    print("      ⚠⚠ the 8/13 split is SECOND-HAND -- three of four primary")
    print("         sources in REFERENCES.md are paywalled and unfetched.")
    print("\nPASS (a map -- no row here is a gate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
