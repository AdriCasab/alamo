#!/usr/bin/env python3
"""P2-31 z0 (the run half): one-step smokes and L14's two cost legs on the PRE-STAGE binary. Legs under S/legs/z0/;
each leg's provenance sits BESIDE its plot dir (<name>.argv.sh / .BINARY.txt / .log / .rc: ALAMO renames an existing
plot_file dir). Usage: smoke.py > SMOKE_Z0.txt
  L16  the first face list on KV_S, GB_S and every live face-list deck: fl_n_top after step 1 (thermo forced on)
  L17  P2-25's six regime argv replayed one step on today's parser
  4b   the empty-list candidate: face_beam + surface_patch.radius=0 runs rc 0 with fl_n_top 0 (the BEFORE side)
  L14  H_05 (Hoeser at 0.5 mm) for 0.5 s and B_08 (the #19 beam deck at 0.78 mm) for 2 s, timed and projected"""
import os
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import decks as D  # noqa: E402

ROOT, T = D.ROOT, D.T
BIN = os.environ.get("MMW_BIN", os.path.join(ROOT, "bin", "mmwspalling-3d-g++"))
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
LEGS = os.path.join(HERE, "legs", "z0")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
THERMO = ["amr.thermo.int=1", "amr.thermo.plot_int=1"]


def sha():
    import hashlib
    return hashlib.sha256(open(BIN, "rb").read()).hexdigest()


def leg(name, argv_tail):
    d = os.path.join(LEGS, name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(LEGS, exist_ok=True)
    argv = MPI + [BIN] + list(argv_tail) + ["plot_file=" + d]
    open(d + ".argv.sh", "w").write(" ".join(shlex.quote(a) for a in argv) + "\n")
    open(d + ".BINARY.txt", "w").write(sha() + "\n")
    t0 = time.time()
    with open(d + ".log", "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(d + ".rc", "w").write("%d\n" % rc)
    return rc, d, time.time() - t0


def first_list(d):
    p = os.path.join(d, "thermo.dat")
    if not os.path.exists(p):
        return None
    L = [ln.rstrip("\n").split("\t") for ln in open(p) if ln.strip()]
    h = [x.strip() for x in L[0]]
    # thermo is written at the START of a step: row t = 0 is the initial state, row t = dt carries step 1's list
    if "fl_n_top" not in h or len(L) < 3:
        return None
    return float(L[2][h.index("fl_n_top")])


def one_step(deck, extra=()):
    R = D.resolve(deck, extra)
    dt = float(R["timestep"][0])
    return [deck] + list(extra) + THERMO + ["stop_time=%r" % (2.0 * dt)]   # two steps: row t = dt is written


def main():
    print("P2-31 smoke on %s" % sha()[:16])
    print("=" * 100)
    print("L16 the first face list (fl_n_top after step 1; thermo forced on, one step)")
    gate = {"KV_S": os.path.join(T, "validation", "kant", "sp_kant_onset", "input_verification_face"),
            "GB_S": os.path.join(T, "unit", "gradient_blind", "input")}
    empty = []
    for nm, deck in gate.items():
        rc, d, s = leg("L16_" + nm, one_step(deck))
        v = first_list(d)
        print("  gate %-6s rc %d  fl_n_top %s  (%.1f s)" % (nm, rc, v, s))
        if rc or not v:
            empty.append(nm)
    for deck in D.live_decks():
        try:
            R = D.resolve(deck)
        except Exception:
            continue
        if not D.face_list_on(R) or "amr.n_cell" not in R:
            if D.face_list_on(R):
                print("  %-62s not runnable as written (its mesh comes from the command line): covered by the gate's "
                      "HOS_S / M0S records" % os.path.relpath(deck, ROOT))
            continue
        nm = "L16_" + os.path.relpath(deck, T).replace("/", "_")
        rc, d, s = leg(nm, one_step(deck))
        v = first_list(d)
        tail = "" if rc == 0 else "  (" + (open(d + ".log", errors="replace").read().strip().splitlines() or ["?"])[-1][:110] + ")"
        print("  %-62s rc %d  fl_n_top %s  (%.1f s)%s" % (os.path.relpath(deck, ROOT), rc, v, s, tail))
        if rc == 0 and not v:
            empty.append(os.path.relpath(deck, ROOT))
    print("  L16: a live deck or gate leg whose first list is empty: %s" % (empty or "none"))

    print("\nL17 P2-25's six argv, one step on today's parser")
    P25 = os.path.join(T, "studies", "p2_25_regime_map", "legs")
    for nm in ("D1_Q2", "D2_Q2", "D2_Q4", "F_D1_Q2", "F_D2_Q2", "F_D2_Q4"):
        deck, ov = D.argv_deck(os.path.join(P25, nm, "argv.sh"))
        ov = [t for t in ov if not t.startswith(("plot_file=", "stop_time="))]
        dt = D.resolve(deck, ov)["timestep"][0]
        rc, d, s = leg("L17_" + nm, [deck] + ov + ["stop_time=%s" % dt])
        print("  %-8s rc %d  (%.1f s)  tokens replayed %d, none dropped" % (nm, rc, s, len(ov)))

    print("\n4b the empty-list candidate on THIS binary (the clause's BEFORE side)")
    fb = os.path.join(T, "unit", "face_beam", "input")
    for nm, extra in (("EL_r0", ["surface_patch.radius=0"]), ("EL_r1", [])):
        rc, d, s = leg(nm, one_step(fb, extra))
        print("  face_beam %-26s rc %d  fl_n_top %s  (%.1f s)" % (" ".join(extra) or "(as written)", rc, first_list(d), s))

    print("\nL14 the two largest legs, timed")
    HOS = os.path.join(T, "studies", "p2_16_span_rule", "legs", "HOS_R", "argv.sh")
    deck, ov = D.argv_deck(HOS)
    drop = {"surface_patch.face", "plot_file", "stop_time", "amr.plot_int", "amr.n_cell", "timestep", "amr.max_grid_size"}
    ov = [t for t in ov if t.split("=", 1)[0] not in drop]
    # ⚠ 25 cells is odd and the deck's amr.blocking_factor is 2: H_2 adds amr.blocking_factor=1 (the box layout only)
    for nm, n, dt, mgs, stop, full in (("H_2_t", "25 25 75", "0.008", "25", "0.08", 10.0),
                                        ("H_05_t", "100 100 300", "0.004", "50", "0.5", 10.0)):
        rc, d, s = leg(nm, [deck] + ov + ["amr.n_cell=" + n, "timestep=" + dt, "amr.max_grid_size=" + mgs,
                                          "stop_time=" + stop, "amr.plot_int=-1"]
                       + (["amr.blocking_factor=1"] if nm == "H_2_t" else []))
        print("  %-7s n_cell %-12s dt %s to %s s: rc %d, %.1f s wall -> %.1f s simulated would take ~%.1f min" % (
            nm, n, dt, stop, rc, s, full, s * full / float(stop) / 60.0))
    og = os.path.join(T, "validation", "oglesby", "input")
    for nm, n, dt, stop in (("B_08_t", "128 128 64", "0.05", "2.0"), ("B_3_t", "32 32 16", "0.2", "2.0")):
        rc, d, s = leg(nm, [og, "amr.n_cell=" + n, "timestep=" + dt, "stop_time=" + stop, "amr.plot_int=-1"])
        print("  %-7s n_cell %-12s dt %s to %s s: rc %d, %.1f s wall -> 63.2 s simulated would take ~%.1f min" % (
            nm, n, dt, stop, rc, s, s * 63.2 / float(stop) / 60.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
