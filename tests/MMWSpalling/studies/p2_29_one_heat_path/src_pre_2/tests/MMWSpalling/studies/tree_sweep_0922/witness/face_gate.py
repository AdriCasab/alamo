#!/usr/bin/env python3
"""THE FACE-LIST GATE (P2-18, A1). The 2,639-file tree gate (hash.sh) never runs the face
list: no gate deck sets surface_patch.face_list = 1, so a change to the production removal
code passes it unseen. This runner is the gate for that code. Every leg is SHORT (seconds
to a minute) and every leg runs the face list.

  usage: face_gate.py run <label>        run every leg into witness/face_<label>/ and write
                                         face_<label>.md5 (+ face_<label>_headers.txt)
         face_gate.py cmp <a> <b>        face_<a>.md5 vs face_<b>.md5: paths only on one
                                         side, differing paths; differing thermo.dat column
                                         by column (the wall-clock timers excluded)
         face_gate.py rankcheck <label>  re-run CAV and T00_1_k1 at -np 1 into
                                         face_<label>_np1/ and compare to face_<label>/

Legs (all -np 4, the test's own runner where one exists):
  fr/T00_1_k1 fr/TL35_1 fr/TL45_1 fr/SPH_1 fr/BORE_1   unit/face_removal fixtures, 8 s (P2-9/P2-10's
                                                       W without its key-0 member; firings from 2.4 s)
  fl/CAV fl/STAIR                                      unit/face_list builder fixtures (1 step)
  fb/F_S                                               unit/face_buckle: span ON (0.015 m), 8 s
  fh/FLAT_P fh/BORE_P                                  unit/face_heat: face_heat = path, 3 steps
  M0S                                                  P2-9's Hoeser leg M0 minus P2-10's 13 retired
                                                       tokens (P2-10's M0T), stop_time 3 s
  HOS_S                                                P2-16's production Hoeser leg HOS_R (enthalpy
                                                       closure, buckle span), stop_time 2 s
  CAL_S                                                validation/kant/grimsel_calibration/input
                                                       (the onset-calibration deck), stop_time 10 s
  BEAM_S                                               unit/face_beam/input (P2-23): the face list WITHOUT
                                                       a flame on the #19 beam deck, stop_time 1 s
  FEET_S FEET_0                                        unit/face_feet/input (P2-27, census G): the burner on its
                                                       feet with the pad clearance ON (regime 4 on the face
                                                       list) and OFF, 10 s
  KV_S HU_S GB_S                                       P2-28: validation/kant/sp_kant_onset/input_verification_face
                                                       (EVALUATE-ONLY, 10 s), validation/hu/hu_spall_onset/
                                                       input_granite2_lefm (the PRESCRIBED patch, buckling on, 10 s),
                                                       unit/gradient_blind/input (P2-28 stage 1: input_face; removal
                                                       on, its 0.04 s)
The three argv-built legs drop `surface_patch.face=zhi` (a dead key with no reader; P2-18 A2). From P2-26 to P2-27 they
also dropped `spall.flake_coherence_length=0.0` / `spall.surface_normal=0` (the table then compared strings); since
P2-27 the retired-at-default table compares a one-token number as a double, so those tokens are replayed as written.

Hashed: Header, Cell_D_*, Cell_H, Checkpoint, thermo.dat, *.csv. Not hashed: run.log, argv.sh,
rc, BINARY.txt, .done, metadata, diff.patch, *.visit (timestamps / argv / git state).
"""
import csv
import hashlib
import importlib.machinery
import importlib.util
import math
import os
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none"]
M0 = os.path.join(T, "studies", "p2_9_face_removal", "legs", "M0", "argv.sh")
HOS_R = os.path.join(T, "studies", "p2_16_span_rule", "legs", "HOS_R", "argv.sh")
CAL = os.path.join(T, "validation", "kant", "grimsel_calibration", "input")
# P2-10's 13 retired tokens (studies/p2_10_retire/run_p2_10.py M0T_DELETE) + the dead key
DROP = {"surface_patch.side_face_flux", "surface_patch.side_face_step_height",
        "surface_patch.side_face_area_law", "surface_patch.side_face_h", "surface_patch.side_face_T",
        "surface_patch.side_face_staircase", "spall.lateral_drive", "spall.lateral_removal",
        "spall.lateral_subcell", "spall.lateral_thermal", "surface_patch.side_face_h_datum",
        "surface_patch.side_face_h_r_min", "surface_patch.side_face_h_law",
        "surface_patch.face",
        # P2-28: the bit's OFF switch (bit.* aborts by presence since the column path retired); inert at 0
        "bit.enabled"}
# P2-27: the two C2 tokens (spall.flake_coherence_length = 0.0, spall.surface_normal = 0) are no longer dropped: the
# retired-at-default table compares a one-token number as a double, so M0S replays its record as written.
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}
SKIP = {"run.log", "argv.sh", "rc", "BINARY.txt", ".done", "metadata", "diff.patch", "DROPPED.txt"}
FR_LEGS = ["T00_1_k1", "TL35_1", "TL45_1", "SPH_1", "BORE_1"]
FR_STOP = "8.0"


def load(path, name):
    """the test files have no .py suffix: name the loader explicitly"""
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stamp(d, rc, t0):
    open(os.path.join(d, "BINARY.txt"), "w").write("bin/mmwspalling-3d-g++ sha256 %s\n" % sha(BIN))
    if rc == 0:
        open(os.path.join(d, ".done"), "w").write("rc 0 wall %.1f\n" % (time.time() - t0))
    print("  %-16s rc=%d wall %6.1f s" % (os.path.relpath(d, os.path.dirname(os.path.dirname(d))), rc,
                                           time.time() - t0), flush=True)


def run_argv(d, argv):
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(a) for a in argv) + "\n")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    stamp(d, rc, t0)
    return rc


def from_record(argv_sh, outdir, override):
    """a recorded argv.sh with the retired/dead tokens dropped and `override` keys set"""
    toks = shlex.split(open(argv_sh).read())
    out, seen = [], set()
    for t in toks:
        k = t.split("=", 1)[0] if "=" in t else None
        if k in DROP:
            continue
        if k in override:
            t = "%s=%s" % (k, override[k])
            seen.add(k)
        if k == "plot_file":
            t = "plot_file=" + os.path.join(outdir, "run")
        out.append(t)
    out += ["%s=%s" % (k, v) for k, v in override.items() if k not in seen]
    return out


# ------------------------------------------------------------------------- run ----
def run(label):
    W = os.path.join(HERE, "face_" + label)
    if os.path.isdir(W):
        shutil.rmtree(W)
    os.makedirs(W)
    print("face gate %s on binary %s" % (label, sha(BIN)[:16]), flush=True)

    fr = load(os.path.join(T, "unit", "face_removal", "test"), "fr_face_gate")
    fr.OUT = os.path.join(W, "fr")
    os.makedirs(fr.OUT)
    k1 = {"surface_patch.face_list": "1", "stop_time": FR_STOP}
    for nm in FR_LEGS:
        t0 = time.time()
        if nm.startswith("T00"):
            rc = fr.run_leg(nm, dict(fr.mesh_keys("1"), **k1))
        elif nm.startswith("TL"):
            th, ph = fr.TILTS[nm.split("_")[0]]
            rc = fr.run_leg(nm, dict(fr.mesh_keys("1"), **k1), seeds=fr.tilt_seeds(th, ph, "1"))
        else:
            c = nm.split("_")[0]
            rc = fr.run_leg(nm, dict(fr.MI.curved_keys(c, 0.001), **k1), seeds=fr.MI.curved_seeds(c, 0.001))
        stamp(os.path.join(fr.OUT, nm), rc, t0)

    fl = load(os.path.join(T, "unit", "face_list", "test"), "fl_face_gate")
    fl.OUT = os.path.join(W, "fl")
    os.makedirs(fl.OUT)
    for nm in ("CAV", "STAIR"):
        t0 = time.time()
        rc = fl.run(nm)
        stamp(os.path.join(fl.OUT, nm), rc, t0)

    fb = load(os.path.join(T, "unit", "face_buckle", "test"), "fb_face_gate")
    fb.OUT = os.path.join(W, "fb")
    os.makedirs(fb.OUT)
    t0 = time.time()
    rc = fb.run("F_S", list(fb.LEGS["F_S"]) + ["stop_time=" + FR_STOP])
    stamp(os.path.join(fb.OUT, "F_S"), rc, t0)

    fh = load(os.path.join(T, "unit", "face_heat", "test"), "fh_face_gate")
    fh.OUT = os.path.join(W, "fh")
    os.makedirs(fh.OUT)
    fh.BASE = [t for t in fh.BASE if t.split("=", 1)[0] not in DROP]
    for nm in ("FLAT_P", "BORE_P"):
        keys, seeds, stop = fh.LEGS[nm]
        t0 = time.time()
        rc = fh.run(nm, keys, seeds, stop)
        stamp(os.path.join(fh.OUT, nm), rc, t0)

    d = os.path.join(W, "M0S")
    run_argv(d, from_record(M0, d, {"stop_time": "3.0", "amr.plot_int": "125"}))
    d = os.path.join(W, "HOS_S")
    run_argv(d, from_record(HOS_R, d, {"stop_time": "2.0", "amr.plot_int": "125"}))
    d = os.path.join(W, "CAL_S")
    run_argv(d, MPI + ["-np", "4", BIN, CAL, "stop_time=10.0", "amr.plot_int=250",
                       "plot_file=" + os.path.join(d, "run")])
    # P2-23: a FLAMELESS list (the beam port), so the beam deck's face list is gated from here on
    d = os.path.join(W, "BEAM_S")
    run_argv(d, MPI + ["-np", "4", BIN, os.path.join(T, "unit", "face_beam", "input"), "stop_time=1.0",
                       "plot_file=" + os.path.join(d, "run")])
    # P2-27 (census G): the pad clearance on the face list, on and off
    for nm, extra in (("FEET_S", []), ("FEET_0", ["surface_patch.foot_body_clearance=0"])):
        d = os.path.join(W, nm)
        run_argv(d, MPI + ["-np", "4", BIN, os.path.join(T, "unit", "face_feet", "input"), "stop_time=10.0",
                           "plot_file=" + os.path.join(d, "run")] + extra)
    # P2-28: what site A served, on the face list -- the evaluate-only mode (Kant's onset deck), the prescribed patch
    # with buckling on (Hu's deck) and the uniform slab with removal on (the locks' stage-2 fixture)
    for nm, deck, stop in (("KV_S", os.path.join(T, "validation", "kant", "sp_kant_onset", "input_verification_face"),
                            "10.0"),
                           ("HU_S", os.path.join(T, "validation", "hu", "hu_spall_onset", "input_granite2_lefm"),
                            "10.0"),
                           ("GB_S", os.path.join(T, "unit", "gradient_blind", "input"), None)):
        d = os.path.join(W, nm)
        run_argv(d, MPI + ["-np", "4", BIN, deck] + (["stop_time=" + stop] if stop else []) +
                 ["plot_file=" + os.path.join(d, "run")])

    n = write_hashes(W, os.path.join(HERE, "face_%s.md5" % label), os.path.join(HERE, "face_%s_headers.txt" % label))
    print("face gate %s: %d files hashed -> face_%s.md5" % (label, n, label), flush=True)


def write_hashes(W, md5_out, hdr_out):
    rows, hdrs = [], []
    for dp, _, fs in os.walk(W):
        for fn in sorted(fs):
            if fn in SKIP or fn.endswith(".visit"):
                continue
            p = os.path.join(dp, fn)
            rel = os.path.relpath(p, W)
            if ".old" in rel:
                continue
            rows.append("%s  %s" % (md5(p), rel))
            if fn == "thermo.dat":
                hdrs.append(rel + "\t" + open(p).readline().rstrip("\n"))
    rows.sort(key=lambda r: r.split("  ", 1)[1])
    open(md5_out, "w").write("\n".join(rows) + "\n")
    open(hdr_out, "w").write("\n".join(sorted(hdrs)) + "\n")
    return len(rows)


# ------------------------------------------------------------------------- cmp ----
def read_md5(label):
    out = {}
    for ln in open(os.path.join(HERE, "face_%s.md5" % label)):
        h, rel = ln.rstrip("\n").split("  ", 1)
        out[rel] = h
    return out


def thermo_cols(path):
    rows = [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]
    hdr = [h.strip() for h in rows[0]]
    cols = {h: [] for h in hdr}
    for r in rows[1:]:
        for h, v in zip(hdr, r):
            try:
                cols[h].append(float(v))
            except ValueError:
                cols[h].append(float("nan"))
    return cols


def col_moved(xa, xb):
    """P2-25 (DR-1): one thermo column against another, element-wise with NaN EQUAL to NaN (a logged NaN -- an
    undefined depth, a beam-off record -- is a value, not a difference). Returns (moved, max|d| over the rows where
    both are finite, the count of rows where exactly one side is NaN)."""
    if len(xa) != len(xb):
        return True, float("inf"), 0
    same = all(x == y or (x != x and y != y) for x, y in zip(xa, xb))
    fin = [abs(x - y) for x, y in zip(xa, xb) if math.isfinite(x) and math.isfinite(y)]
    nan_mis = sum(1 for x, y in zip(xa, xb) if (x != x) != (y != y))
    return (not same), (max(fin) if fin else 0.0), nan_mis


def cmp(a, b):
    A, B = read_md5(a), read_md5(b)
    only_a = sorted(set(A) - set(B))
    only_b = sorted(set(B) - set(A))
    common = sorted(set(A) & set(B))
    diff = [r for r in common if A[r] != B[r]]
    L = ["face_%s vs face_%s: %d / %d files; %d common; only in %s: %d; only in %s: %d; DIFFERING: %d"
         % (a, b, len(A), len(B), len(common), a, len(only_a), b, len(only_b), len(diff))]
    for r in only_a[:20]:
        L.append("  only %s: %s" % (a, r))
    for r in only_b[:20]:
        L.append("  only %s: %s" % (b, r))
    real = 0
    for r in diff:
        if r.endswith("thermo.dat"):
            ca = thermo_cols(os.path.join(HERE, "face_" + a, r))
            cb = thermo_cols(os.path.join(HERE, "face_" + b, r))
            moved = []
            for h in ca:
                if h in TIMERS or h not in cb:
                    continue
                mv, m, nan_mis = col_moved(ca[h], cb[h])
                if mv:
                    moved.append("%s(max|d| %.3g%s)" % (h, m, ", %d NaN-mismatched rows" % nan_mis if nan_mis else ""))
            gone = sorted(set(ca) - set(cb))
            new = sorted(set(cb) - set(ca))
            tag = "timers only" if not (moved or gone or new) else "MOVED " + " ".join(moved + ["-" + g for g in gone] + ["+" + n for n in new])
            if moved or gone or new:
                real += 1
            L.append("  thermo %-40s %s" % (r, tag))
        else:
            real += 1
            L.append("  DIFF   %s" % r)
    L.append("VERDICT: %s (%d substantive difference(s); thermo wall-clock timers %s are never compared)"
             % ("IDENTICAL" if real == 0 and not only_a and not only_b else "DIFFERENT", real, sorted(TIMERS)))
    txt = "\n".join(L)
    print(txt)
    return real == 0 and not only_a and not only_b


# ------------------------------------------------------------------- rankcheck ----
def rankcheck(label):
    W = os.path.join(HERE, "face_" + label)
    W1 = W + "_np1"
    if os.path.isdir(W1):
        shutil.rmtree(W1)
    os.makedirs(W1)
    fl = load(os.path.join(T, "unit", "face_list", "test"), "fl_face_gate_np1")
    fl.OUT = os.path.join(W1, "fl")
    os.makedirs(fl.OUT)
    d = os.path.join(fl.OUT, "CAV")
    os.makedirs(d)
    boxes, keys = fl.LEGS["CAV"]
    sv = " ".join(str(v) for b in boxes for v in b)
    run_argv(d, MPI + ["-np", "1", BIN, fl.INPUT, "plot_file=%s/run" % d, "spall.seed_void=" + sv] + keys)
    fr = load(os.path.join(T, "unit", "face_removal", "test"), "fr_face_gate_np1")
    fr.OUT = os.path.join(W1, "fr")
    os.makedirs(fr.OUT)
    t0 = time.time()
    rc = fr.run_leg("T00_1_k1", dict(fr.mesh_keys("1"), **{"surface_patch.face_list": "1", "stop_time": FR_STOP}), np_=1)
    stamp(os.path.join(fr.OUT, "T00_1_k1"), rc, t0)
    write_hashes(W1, os.path.join(HERE, "face_%s_np1.md5" % label), os.path.join(HERE, "face_%s_np1_headers.txt" % label))
    # Plotfile bytes (Cell_D / Cell_H / Header) depend on the box-to-rank layout, so they can
    # never be byte-identical across rank counts; the PHYSICS is compared: every thermo.dat
    # column numerically (max |d|, the wall-clock timers excluded) and every CSV byte for byte.
    A, B = read_md5(label), read_md5(label + "_np1")
    L = ["rank check %s: -np 1 vs -np 4, legs fl/CAV and fr/T00_1_k1 (plotfile bytes are layout-dependent"
         " and are not compared; thermo columns numerically, CSVs byte for byte)" % label]
    bad = 0
    worst = 0.0
    for r in sorted(B):
        if r not in A:
            L.append("  missing at np4: %s" % r)
            bad += 1
        elif r.endswith("thermo.dat"):
            ca, cb = thermo_cols(os.path.join(W, r)), thermo_cols(os.path.join(W1, r))
            if len(next(iter(ca.values()))) != len(next(iter(cb.values()))):
                L.append("  thermo %-30s ROW COUNT differs" % r)
                bad += 1
                continue
            moved = []
            for h in ca:
                if h in TIMERS or h not in cb or ca[h] == cb[h]:
                    continue
                m = max(abs(x - y) for x, y in zip(ca[h], cb[h]))
                scale = max(abs(x) for x in ca[h])
                # a column that IS a round-off residual (the ledger closure errors, ~1e-15) has no
                # meaningful relative error: score it on the absolute difference instead
                score = m / scale if scale > 1e-9 else m
                moved.append("%s(max|d| %.2e, %s %.1e)" % (h, m, "rel" if scale > 1e-9 else "abs", score))
                worst = max(worst, score)
            L.append("  thermo %-30s %s" % (r, "identical but the timers" if not moved else " ".join(moved)))
        elif r.endswith(".csv"):
            if A[r] != B[r]:
                L.append("  CSV DIFFERS %s" % r)
                bad += 1
            else:
                L.append("  csv    %-30s identical" % r)
    L.append("VERDICT: %s" % ("decomposition-stable: counts and CSVs identical, thermo columns agree to rel %.1e (round-off)" % worst
                               if bad == 0 and worst < 1e-9 else "%d difference(s), worst rel %.1e" % (bad, worst)))
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "face_%s_RANKCHECK.txt" % label), "w").write(txt + "\n")
    return bad == 0


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "run":
        run(sys.argv[2])
    elif len(sys.argv) >= 4 and sys.argv[1] == "cmp":
        sys.exit(0 if cmp(sys.argv[2], sys.argv[3]) else 1)
    elif len(sys.argv) >= 3 and sys.argv[1] == "rankcheck":
        sys.exit(0 if rankcheck(sys.argv[2]) else 1)
    else:
        print(__doc__)
        sys.exit(1)
