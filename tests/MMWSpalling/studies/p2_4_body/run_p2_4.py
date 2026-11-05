#!/usr/bin/env python3
"""P2-4 -- the cold body in the annulus.

Legs (SERIAL; five fast cores, four ranks each -- nothing concurrent):
  S0   W1's configuration with the WALL OFF and the march off, 500 s. The
       control P2-3 never ran: the floor scale ALONE. It settles P5 (which rate
       branch holds) before any wall or body number is interpreted, which is why
       it runs FIRST.
  B1L  W1's configuration + the cold body, 800 s. ONE run: its first 500 s are
       the mechanism leg (b1-b5, c1), scored AT THE 500 s MARK from the live
       logs and KILLED on any failure; its 500-800 s are the finding (c2-c9).

⚠⚠ `--gate-at 500` SCORES WHILE THE RUN IS STILL GOING and terminates it on a
failure. The alternative -- let 800 s finish and score afterwards -- spends
51 minutes to discover at minute 32 that the energy closure broke. The kill is
by PID on the Popen handle; `pkill -f` is forbidden here and would in any case
match this script's own command line.

⚠ The base configuration is built from P2-3's OWN harness
(`run_p2_3.base_keys() + COMMON + LEGS["W1"]`), de-duplicated LAST-WINS to match
ParmParse -- NOT parsed from `legs/W1/argv.sh`, which carries duplicate keys
(`side_face_factor_above` appears as both 0.0 and 1.0) and a stray sha on the
`side_face_h_law` line. A first-match parse of that file silently selects THE
WALL OFF. See `gate_a1.base_keys` for the full note.

⚠ Depth guard (b5): `stop_time` is 800 s because that is the longest stop that
keeps >= 60 mm of rock under BOTH rate branches. It is not raised without
redoing the projection.

Usage: run_p2_4.py S0 | B1L [--gate-at 500] | reduce [legs...]
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                    "sp_meier_pilot", "input_drilling")
OUT = os.path.join(HERE, "legs")
W1_DIR = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_3_wall",
                      "legs", "W1")

USABLE_MM = 500.0          # 600 mm of rock less the 40 mm guard and 60 mm margin

LEGS = {
    "S0":  (["surface_patch.side_face_factor_above=0.0",
             "surface_patch.jet_wall_march=0"], 500.0),
    "B1L": (["surface_patch.jet_body=1",
             "surface_patch.jet_body_T=320.0"], 800.0),
    # ⚠⚠ THE b1 DISCRIMINATOR, added after b1 failed by 2 rows of 30,784.
    # b1 compares B1L (P2-4 binary) against W1 (P2-3 binary), so a mismatch has
    # TWO possible causes and they demand opposite responses:
    #   (i)  the body acting before an annulus exists -- a real leak, must fix;
    #   (ii) the RECOMPILE, since this toolchain is not byte-reproducible
    #        (standing rule 9) and P2-4 changed code motion in Parse and in the
    #        bin loop.
    # NB is the SAME binary as B1L with `jet_body=0`, run just past the crossing.
    #   NB == B1L through the crossing  => the body is provably inert there and
    #                                      the W1 mismatch is the binary;
    #   NB != B1L                       => the body IS acting early: a real bug.
    #   NB vs W1 also differing by ~2 rows is the positive confirmation of (ii).
    # ⚠ This is a DIAGNOSTIC of the failed gate, not a continuation of the
    # campaign -- no further 500/800 s leg is launched.
    "NB":  (["surface_patch.jet_body=0"], 160.0),
}


def base_keys():
    """W1's configuration, from P2-3's harness, last-wins."""
    sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies",
                                    "p2_3_wall"))
    import run_p2_3 as R3
    raw = list(R3.base_keys()) + list(R3.COMMON) + list(R3.LEGS["W1"])
    seen = {}
    for k in raw:
        if "=" in k:
            seen[k.split("=", 1)[0]] = k
    return list(seen.values())


def argv_for(tag, extra, stop, d):
    keys = {}
    for k in base_keys() + list(extra):
        keys[k.split("=", 1)[0]] = k
    return (["caffeinate", "-dis", "mpirun", "--oversubscribe", "--bind-to",
             "none", "-np", "4", BIN, DECK]
            + list(keys.values())
            + ["stop_time=%g" % stop,
               "amr.thermo.int=1", "amr.thermo.plot_int=6",
               "surface_patch.jet_profile_interval=2.0",
               "energy_ledger.enabled=1", "spall.removal_events_csv=1",
               "plot_file=%s" % os.path.join(d, "run")])


def thermo(d):
    p = os.path.join(d, "run", "thermo.dat")
    if not os.path.exists(p):
        return None, []
    try:
        txt = open(p).read()
    except OSError:
        return None, []
    lines = [l.split() for l in txt.splitlines() if l.split()]
    if len(lines) < 2:
        return None, []
    return lines[0], lines[1:]


def col(head, rows, name):
    if not head or name not in head:
        return []
    i = head.index(name)
    out = []
    for r in rows:
        if len(r) > i:
            try:
                out.append(float(r[i]))
            except ValueError:
                pass
    return out


def deepest_void_mm(d):
    """from the event log: rock top 0.500, cell 2 mm. ⚠ k_top is the top BEFORE
    the removal, so the post-removal top is k_top - n_voided (standing rule 4)."""
    p = os.path.join(d, "run_removal_events.csv")
    if not os.path.exists(p):
        return None
    deep = 0.0
    try:
        with open(p) as f:
            for r in csv.DictReader(f):
                try:
                    k = float(r["k_top"]) - float(r["n_voided"])
                except (TypeError, ValueError, KeyError):
                    continue
                z = -0.1 + (k + 0.5) * 0.002
                deep = max(deep, (0.500 - z) * 1e3)
    except OSError:
        return None
    return deep


def ev_prefix_md5(d, t_max):
    """md5 of the event-log rows with time < t_max, for b1's identity."""
    import hashlib
    p = os.path.join(d, "run_removal_events.csv")
    if not os.path.exists(p):
        return None, 0
    h = hashlib.md5()
    n = 0
    with open(p) as f:
        head = f.readline()
        h.update(head.encode())
        for ln in f:
            try:
                t = float(ln.split(",", 1)[0])
            except ValueError:
                continue
            if t >= t_max:
                break
            h.update(ln.encode())
            n += 1
    return h.hexdigest(), n


def score_mark(d, tag, out):
    """b1-b5 and c1 at the 500 s mark. Returns (ok, lines)."""
    L, bad = [], []

    def need(ok, msg, det=""):
        L.append("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if det:
            L.append("        %s" % det)
        if not ok:
            bad.append(msg)

    head, rows = thermo(d)
    if not head:
        return False, ["  FAIL  no thermo.dat at the mark"]
    t = col(head, rows, "time")
    sel = [i for i, x in enumerate(t) if x <= 500.0]
    if not sel:
        return False, ["  FAIL  no rows at or before 500 s"]
    last = sel[-1]

    # ---- b2: the body was REACHED ----
    wb = col(head, rows, "wall_bins")
    bh = col(head, rows, "body_h_used")
    nl = col(head, rows, "body_noLaw_n")
    withbins = [i for i in sel if i < len(wb) and wb[i] > 0.0]
    reached = [i for i in withbins if i < len(bh) and bh[i] > 0.0]
    frac = (len(reached) / len(withbins)) if withbins else 0.0
    need(bool(withbins) and frac >= 0.95,
         "b2 body reached: body_h_used > 0 on >= 95 % of steps with bins",
         "%d of %d stepswith bins (%.2f %%); body_noLaw_n at the mark %.0f "
         "(CUMULATIVE)" % (len(reached), len(withbins), 100 * frac,
                           nl[last] if last < len(nl) else 0.0))

    # ---- b3: the ledger ----
    les = col(head, rows, "ledger_err_solid")
    worst = max((abs(les[i]) for i in sel[1:] if i < len(les)), default=0.0)
    need(worst <= 5e-3, "b3 ledger_err_solid <= 5e-3 at the mark",
         "max %.4e" % worst)

    # ---- b4: pads descend ----
    nz = col(head, rows, "nozzle_z")
    st = col(head, rows, "foot_stall_time")
    pts = [(t[i], nz[i]) for i in sel if i < len(nz) and 150.0 <= t[i] <= 350.0
           and nz[i] > 0]
    v = 0.0
    if len(pts) > 1 and pts[-1][0] > pts[0][0]:
        v = (pts[0][1] - pts[-1][1]) / (pts[-1][0] - pts[0][0]) * 3600.0
    stall = max((st[i] for i in sel if i < len(st)), default=0.0)
    need(v > 0.8 and stall < 10.0,
         "b4 pads descend: v > 0.8 m/h over 150-350 s and foot_stall_time < 10 s",
         "v = %.3f m/h, max foot_stall_time %.3f s" % (v, stall))

    # ---- b5: THE DEPTH GUARD ----
    deep = deepest_void_mm(d) or 0.0
    pts2 = [(t[i], nz[i]) for i in sel if i < len(nz) and 350.0 <= t[i] <= 500.0
            and nz[i] > 0]
    v2 = 0.0
    if len(pts2) > 1 and pts2[-1][0] > pts2[0][0]:
        v2 = (pts2[0][1] - pts2[-1][1]) / (pts2[-1][0] - pts2[0][0]) * 3600.0
    proj = deep + v2 / 3600.0 * 300.0 * 1e3
    need(proj <= USABLE_MM,
         "b5 DEPTH GUARD: deepest void at 500 s + v(350-500) x 300 s <= %.0f mm"
         % USABLE_MM,
         "deepest %.1f mm + %.3f m/h x 300 s = %.1f mm projected at 800 s"
         % (deep, v2, proj))

    # ---- c1: the mouth gas vs W1 at 500 s ----
    hw, rw = thermo(W1_DIR)
    tm = col(head, rows, "wall_T_mouth")
    mine = tm[last] if last < len(tm) else 0.0
    w1 = 0.0
    if hw:
        tw = col(hw, rw, "time")
        mw = col(hw, rw, "wall_T_mouth")
        if tw and mw:
            j = min(range(len(mw)), key=lambda k: abs(tw[k] - 500.0))
            w1 = mw[j]
    drop = w1 - mine
    L.append("  c1 wall_T_mouth at the 500 s mark: %.1f K vs W1's %.1f K "
             "=> %.1f K lower  (bar 100-130)" % (mine, w1, drop))
    L.append("     %s c1" % ("ok  " if 100.0 <= drop <= 130.0 else "MISS"))

    # ---- b1: the pre-crossing identity ----
    # ⚠⚠ TWO CLOCKS, AND MIXING THEM COST THIS GATE A FALSE FAILURE.
    # `thermo` stamps a step's END, the event log stamps its START (standing
    # rule 7(iii), which this packet's own guardrails cite). The first
    # `wall_bins > 0` thermo row is at 149.472 while that SAME step's event rows
    # are stamped 149.424 -- so a window of `event_time < 149.472` INCLUDES the
    # rows of the first bin-bearing step, which is exactly the step where the
    # body first legitimately acts. Scored that way, b1 failed by 2 rows of
    # 30,784 and the run was killed; the control (NB, same binary, body off) is
    # byte-identical to W1, so there was never a leak.
    # ⇒ bound the window by the PREVIOUS thermo row's time, which is strictly
    #   before the first bin-bearing step on EITHER clock.
    i_bin = next((i for i in sel if i < len(wb) and wb[i] > 0.0), None)
    first_bin = t[i_bin - 1] if (i_bin is not None and i_bin > 0) else None
    if first_bin is None:
        need(False, "b1 pre-crossing identity: no above-plane step found")
    else:
        a, na = ev_prefix_md5(d, first_bin)
        b, nb = ev_prefix_md5(W1_DIR, first_bin)
        need(a is not None and a == b,
             "b1 event log byte-identical to W1 strictly BEFORE the first "
             "above-plane step (event rows t < %.3f s; the step itself ends at "
             "%.3f s in thermo -- two clocks, see the note above)"
             % (first_bin, t[i_bin]),
             "%d rows here vs %d in W1; md5 %s vs %s"
             % (na, nb, (a or "-")[:12], (b or "-")[:12]))
    with open(out, "w") as f:
        f.write("P2-4 -- %s scored AT THE 500 s MARK (live)\n" % tag)
        f.write("=" * 72 + "\n")
        f.write("\n".join(L) + "\n")
        f.write("\n" + ("ALL GATES PASS -- the run continues to 800 s\n"
                        if not bad else
                        "GATE FAILURE: %s\n  => the run is KILLED at the mark\n"
                        % "; ".join(bad)))
    return (not bad), L


def run(tag, gate_at=None):
    extra, stop = LEGS[tag]
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    done = os.path.join(d, ".done")
    if os.path.exists(done):
        os.remove(done)
    argv = argv_for(tag, extra, stop, d)
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join("'%s'" % a if " " in a else a for a in argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    t0 = time.time()
    lg = open(os.path.join(d, "run.log"), "w")
    proc = subprocess.Popen(argv, cwd=ROOT, env=env, stdout=lg,
                            stderr=subprocess.STDOUT)
    scored, killed = False, False
    while proc.poll() is None:
        time.sleep(20)
        if gate_at and not scored:
            head, rows = thermo(d)
            tt = col(head, rows, "time") if head else []
            if tt and tt[-1] >= gate_at:
                ok, lines = score_mark(
                    d, tag, os.path.join(HERE, "MARK_%s.txt" % tag))
                print("\n  --- %s scored at the %.0f s mark ---" % (tag, gate_at))
                print("\n".join(lines))
                scored = True
                if not ok:
                    print("  => GATE FAILURE at the mark: terminating the run "
                          "(by PID, not pkill)")
                    proc.terminate()
                    killed = True
    rc = proc.wait()
    lg.close()
    wall = time.time() - t0
    with open(done, "w") as f:
        f.write("rc=%d wall=%.1f killed=%d scored=%d\n"
                % (rc, wall, int(killed), int(scored)))
    return rc, wall, killed


def main():
    if not sys.argv[1:]:
        print(__doc__)
        return 2
    if sys.argv[1] == "reduce":
        import reduce_p2_4
        return reduce_p2_4.main(sys.argv[2:] or None)
    gate_at = None
    if "--gate-at" in sys.argv:
        gate_at = float(sys.argv[sys.argv.index("--gate-at") + 1])
    tags = [a for a in sys.argv[1:] if a in LEGS]
    if not tags:
        print("  no known leg named; have: %s" % ", ".join(LEGS))
        return 2
    print("P2-4 -- legs %s" % ", ".join(tags))
    print("  binary %s" % subprocess.run(["shasum", "-a", "256", BIN],
                                         capture_output=True,
                                         text=True).stdout[:24])
    print("  base config: %d keys from P2-3's harness (last-wins)"
          % len(base_keys()))
    bad = []
    for tg in tags:
        rc, wall, killed = run(tg, gate_at if tg == "B1L" else None)
        print("  %-4s rc=%d wall %.1f s (%.2f wall-s/s)%s"
              % (tg, rc, wall, wall / LEGS[tg][1],
                 "  KILLED AT THE MARK" if killed else ""))
        if rc != 0 and not killed:
            bad.append(tg)
    if bad:
        print("  ABORTED: %s" % ", ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
