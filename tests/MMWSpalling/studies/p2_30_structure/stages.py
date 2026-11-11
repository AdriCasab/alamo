#!/usr/bin/env python3
"""P2-30: stages 2, 3, 4a, 4b of the driver (check.sh calls this). Runs a stage's Criteria rows in table order; a red GATE row
stops the stage with exit 1; a red verdict / report row is recorded and the stage continues. Output: S/CHECK_<stage>.log
(every command's output) and S/RESULTS.txt (one line per row)."""
import collections
import glob
import hashlib
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
GATE = os.path.join(W, "face_gate.py")
ENV = dict(os.environ, EIGEN=os.path.join(ROOT, "ext"), CPLUS_INCLUDE_PATH="/opt/homebrew/include",
           LIBRARY_PATH="/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current")
STAGE = sys.argv[1]
LOG = os.path.join(HERE, "CHECK_%s.log" % STAGE)
RES = os.path.join(HERE, "RESULTS.txt")


def say(s):
    print(s, flush=True)
    with open(RES, "a") as f:
        f.write(s + "\n")


def run(cmd, tag=None, env=None, out=None):
    """cmd's output appended to the stage log (and to `out`, when given); returns its rc"""
    with open(LOG, "a") as lg:
        lg.write("\n$ %s\n" % " ".join(cmd))
        lg.flush()
        r = subprocess.run(cmd, cwd=ROOT, env=env or os.environ, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        lg.write(r.stdout)
    if out:
        with open(out, "w") as f:
            f.write(r.stdout)
    return r.returncode


def row(kind, name, cmd, env=None, out=None):
    t0 = time.time()
    rc = cmd() if callable(cmd) else run(cmd, env=env, out=out)
    say("  %s: %s (rc %d, %.0f s) [%s]" % (name, "PASS" if rc == 0 else "RED", rc, time.time() - t0, kind))
    if rc != 0 and kind == "gate":
        say("== RED GATE %s: stopped at %s" % (name, time.strftime("%H:%M:%S")))
        sys.exit(1)
    return rc


def note_tail(pattern, n=20):
    """the lines of the stage log matching `pattern`, into RESULTS (the numbers behind a row)"""
    hits = [ln for ln in open(LOG, errors="replace").read().splitlines() if re.search(pattern, ln)]
    for h in hits[-n:]:
        say("    " + h.strip()[:200])


def sha():
    return hashlib.sha256(open(BIN, "rb").read()).hexdigest()


def build(label):
    def go():
        out = os.path.join(HERE, "build_%s.log" % label)
        rc = run(["make", "-j8"], env=ENV, out=out)
        with open(os.path.join(HERE, "BINARY_%s.txt" % label), "w") as f:
            f.write("%s  bin/mmwspalling-3d-g++\n" % sha())
        say("    binary %s: %s" % (label, sha()[:16]))
        return rc
    return go


def warnings_multiset(path):
    c = collections.Counter()
    for ln in open(path, errors="replace").read().splitlines():
        m = re.search(r"warning: (.*)$", ln)
        if m:
            c[m.group(1).strip()] += 1
    return c


def no_new_warnings(label, base):
    def go():
        new = warnings_multiset(os.path.join(HERE, "build_%s.log" % label))
        old = warnings_multiset(base)
        extra = new - old
        print_lines = ["warnings %s: %d (distinct %d); reference %s: %d; new messages: %d" % (
            label, sum(new.values()), len(new), os.path.basename(base), sum(old.values()), sum(extra.values()))]
        print_lines += ["  NEW: " + k[:180] for k in extra]
        with open(LOG, "a") as f:
            f.write("\n".join(print_lines) + "\n")
        say("    " + print_lines[0])
        return 1 if extra else 0
    return go


def lines(p):
    return open(p, errors="replace").read().count("\n")


def m1():
    d = os.path.join(ROOT, "src", "Integrator", "MMWSpalling")
    p = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
    say("  m1 (report): parent %d lines (prediction <= 600); slices: %s" % (lines(p), ", ".join(
        "%s %d" % (os.path.basename(x), lines(x)) for x in sorted(glob.glob(os.path.join(d, "*.H"))))))
    moved = sum(1 for _ in open(os.path.join(HERE, "COMMENT_MOVED.txt")))
    say("  m1: comment lines moved %d of the 586 candidate lines (prediction 200-450 of ~610)" % moved)
    u = os.path.join(T, "unit")
    old = sum(lines(os.path.join(HERE, "src_pre_1", "tests", n, "test")) for n in os.listdir(os.path.join(HERE, "src_pre_1", "tests")) if os.path.isdir(os.path.join(HERE, "src_pre_1", "tests", n)))
    new = sum(lines(os.path.join(u, n, "test")) for n in ("lateral_flux_lock", "retired_keys", "buckle_lock", "profile_depth", "face_heat", "face_removal"))
    say("  m1: the six source-reading tests %d -> %d lines (prediction -10 to -30 %%: %.1f %%)" % (old, new, 100.0 * (new - old) / old))


def stage2():
    row("gate", "a3_identity", [PY, os.path.join(HERE, "pp_identity.py"), "src_pre_2"])
    note_tail(r"tokens compared|VERDICT", 2)
    row("gate", "a4_sort", [PY, os.path.join(HERE, "sort_check.py"), "sort"])
    note_tail(r"candidates;|multisets|non-comment lines|flagged moved|VERDICT", 5)
    row("gate", "a5_capture", [PY, os.path.join(HERE, "locks.py"), "capture", "stage2"])
    row("gate", "a5_trim", [PY, os.path.join(HERE, "sort_check.py"), "trim", "stage1", "stage2"])
    note_tail(r"clause line\(s\) declared|UNDECLARED|declared but|lateral_flux_lock  .*removed", 4)
    row("gate", "a6_build", build("p2_30a"))
    row("gate", "a6_warnings", no_new_warnings("p2_30a", os.path.join(T, "studies", "p2_29_one_heat_path", "build_2.log")))
    row("gate", "a6_gate_run", [PY, GATE, "run", "p2_30a"])
    row("gate", "a6_gate_cmp", [PY, GATE, "cmp", "p2_28", "p2_30a"])
    note_tail(r"^face_p2_28 vs|VERDICT", 2)
    m1()


# ---- stage 3: the tests (no source change) --------------------------------------------------------------------------------
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
LEGS = os.path.join(HERE, "legs")


def leg(name, deck, extra=()):
    """a run into S/legs/<name>/ (never over a test's own output/: list N's baselines live there). The provenance sits BESIDE
    that directory (<name>.argv.sh / .BINARY.txt / .log / .rc): ALAMO owns plot_file and renames an existing one (SR 8)"""
    d = os.path.join(LEGS, name)
    if os.path.isdir(d):
        import shutil
        shutil.rmtree(d)
    os.makedirs(LEGS, exist_ok=True)
    argv = MPI + [BIN, deck] + list(extra) + ["plot_file=" + d]
    open(d + ".argv.sh", "w").write(" ".join(argv) + "\n")
    open(d + ".BINARY.txt", "w").write(sha() + "\n")
    with open(d + ".log", "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(d + ".rc", "w").write("%d\n" % rc)
    return rc, d


def c1():
    rc, d = leg("zhang_w0.01", os.path.join(T, "validation", "zhang", "zhang_oglesby", "input"))
    if rc:
        say("    c1: the Zhang run failed (rc %d)" % rc)
        return 1
    r = subprocess.run([PY, os.path.join(T, "validation", "zhang", "zhang_oglesby", "test"), d], cwd=ROOT, capture_output=True, text=True)
    open(os.path.join(HERE, "c1_zhang_test.log"), "w").write(r.stdout + r.stderr)
    with open(LOG, "a") as f:
        f.write(r.stdout + r.stderr)
    lines_ = r.stdout.strip().splitlines()
    last = lines_[-1] if lines_ else ""
    want = {50: 823.0, 100: 2184.8, 150: 2526.4}
    got = {}
    for ln in lines_:
        m = re.match(r"\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s", ln)
        if m and int(m.group(1)) in want and "t [s]" not in ln:
            got[int(m.group(1))] = float(m.group(2))
    onset = next((float(m.group(1)) for m in (re.search(r"c6 onset ([\d.]+) s", ln) for ln in lines_) if m), None)
    b12 = all(("b1" in ln or "b2" in ln) and "PASS" in ln for ln in lines_ if re.search(r"PASS|FAIL", ln) and ("b1 " in ln or "b2 " in ln))
    same = all(abs(got.get(k, -1) - v) < 0.05 for k, v in want.items()) and onset is not None and abs(onset - 72.0) < 0.05
    say("    c1: raw T at 50 / 100 / 150 s = %s (want 823.0 / 2184.8 / 2526.4); onset %s s (want 72.0); test rc %d" % (
        " / ".join("%.1f" % got.get(k, float("nan")) for k in want), onset, r.returncode))
    say("    c1: b1 and b2 PASS: %s; last line: %s" % (b12, last[:170]))
    return 0 if r.returncode == 0 and b12 and same and last.startswith("PASS: beam energy only") else 1


def c2():
    rc, d = leg("stefan", os.path.join(T, "unit", "stefan", "input"))
    if rc:
        say("    c2: the Stefan run failed (rc %d)" % rc)
        return 1
    t = os.path.join(T, "unit", "stefan", "test")
    r1 = subprocess.run([PY, t, d], cwd=ROOT, capture_output=True, text=True)
    r2 = subprocess.run([PY, t, d, "--plant"], cwd=ROOT, capture_output=True, text=True)
    r3 = subprocess.run([PY, t, d, "--plant-hot"], cwd=ROOT, capture_output=True, text=True)
    with open(LOG, "a") as f:
        f.write(r1.stdout + r1.stderr + "\n--- planted (median) ---\n" + r2.stdout + r2.stderr + "\n--- planted (hottest) ---\n" + r3.stdout + r3.stderr)
    for ln in r1.stdout.splitlines():
        if ln.startswith(("liquid region", "H(T)", "Stefan front", "PASS", "FAIL")):
            say("    c2 run: " + ln[:190])
    for tag, r in (("median +5 K", r2), ("hottest +6 K", r3)):
        for ln in r.stdout.splitlines():
            if ln.startswith(("PLANTED", "FAIL")):
                say("    c2 planted (%s): %s" % (tag, ln[:170]))
    return 0 if r1.returncode == 0 and r2.returncode == 1 and "H(T) not monotone" in r2.stdout \
        and r3.returncode == 1 and "the liquid region's maximum T" in r3.stdout else 1


def c3():
    rc, d = leg("equilibrium", os.path.join(T, "unit", "equilibrium", "input"))
    if rc:
        say("    c3: the equilibrium run failed (rc %d)" % rc)
        return 1
    r = subprocess.run([PY, os.path.join(T, "unit", "equilibrium", "test"), d], cwd=ROOT, capture_output=True, text=True)
    with open(LOG, "a") as f:
        f.write(r.stdout + r.stderr)
    for ln in r.stdout.splitlines():
        if ln.startswith(("alpha dz", "steady", "Surface equilibrium", "L8", "PASS", "FAIL", "Temp range")):
            say("    c3: " + ln[:190])
    return r.returncode


def stage3():
    row("verdict", "c1_zhang", c1)
    row("verdict", "c2_stefan", c2)
    row("verdict", "c3_equilibrium", c3)


# ---- stage 4a: one column list per log file --------------------------------------------------------------------------------
COVERS = os.path.join(HERE, "covers.py")
FR_H = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "FaceRemoval.H")


def gate_dones(label):
    """every leg of the face gate ran to rc 0 on this binary (face_gate stamps .done only on rc 0)"""
    def go():
        W_ = os.path.join(W, "face_" + label)
        legs_ = [dp for dp, _, fs in os.walk(W_) if "BINARY.txt" in fs]
        done = [d for d in legs_ if os.path.exists(os.path.join(d, ".done"))]
        say("    b2: %d gate legs ran, %d of them to rc 0 (the field-count check ran on every row those legs wrote)" % (len(legs_), len(done)))
        return 0 if len(legs_) == 19 and len(done) == len(legs_) else 1
    return go


def moved_backtraces(t0, d):
    import shutil
    for b in os.listdir(ROOT):
        if b.startswith("Backtrace.") and os.path.getmtime(os.path.join(ROOT, b)) >= t0 - 1.0:
            shutil.move(os.path.join(ROOT, b), os.path.join(d, b))


def b2_plant():
    """a planted extra field in the vapour row must ABORT (naming the file and the two counts); the source is restored"""
    orig = open(FR_H).read()
    old = '<< h_v << "," << h_v << ",0,0,0," << c.sheet << ",0,0,0";'
    if orig.count(old) != 1:
        say("    b2 plant: the planted target was not found once")
        return 2
    t0 = time.time()
    try:
        open(FR_H, "w").write(orig.replace(old, '<< h_v << "," << h_v << ",0,0,0," << c.sheet << ",0,0,0,0";'))
        rc_b = build("p2_30b_plant")()
        if rc_b:
            say("    b2 plant: the planted build failed")
            return 3
        rc, d = leg("plant4a", os.path.join(T, "unit", "face_vapor", "input"))
        log = open(d + ".log", errors="replace").read()
        moved_backtraces(t0, LEGS)
        hit = [ln for ln in log.splitlines() if "_removal_events.csv" in ln and "fields" in ln]
        say("    b2 plant: an extra field in the vapour row: rc %d; abort text: %s" % (rc, hit[0][:170] if hit else "(none)"))
        ok = rc != 0 and bool(hit)
    finally:
        open(FR_H, "w").write(orig)
        rc_r = build("p2_30b_restored")()
    say("    b2 plant: source restored; rebuilt binary %s (p2_30b: %s)" % (sha()[:16], open(os.path.join(HERE, "BINARY_p2_30b.txt")).read()[:16]))
    return 0 if ok and rc_r == 0 else 1


def stage4a():
    row("gate", "b1_covers_before", [PY, COVERS, "pre4a"])
    note_tail(r"^  [VB] run_", 3)
    row("gate", "b1_build", build("p2_30b"))
    row("report", "b1_warnings", no_new_warnings("p2_30b", os.path.join(HERE, "build_p2_30a.log")))
    row("gate", "b1_gate_run", [PY, GATE, "run", "p2_30b"])
    row("gate", "b1_gate_cmp", [PY, GATE, "cmp", "p2_28", "p2_30b"])
    note_tail(r"^face_p2_28 vs|VERDICT", 2)
    row("gate", "b1_covers_after", [PY, COVERS, "post4a"])
    row("gate", "b1_covers_bytes", [PY, COVERS, "cmp", "pre4a", "post4a"])
    note_tail(r"^  [VB] run_.*bytes", 3)
    row("verdict", "b2_every_row", gate_dones("p2_30b"))
    row("verdict", "b2_planted", b2_plant)


# ---- stage 4b: the patch's three constants, one definition ------------------------------------------------------------------
import importlib.machinery  # noqa: E402
import importlib.util  # noqa: E402
import shutil  # noqa: E402

PRE_4B = os.path.join(HERE, "src_pre_4b", "Integrator")
SRC_INT = os.path.join(ROOT, "src", "Integrator")
P29 = os.path.join(T, "studies", "p2_29_one_heat_path", "legs2")
_ldr = importlib.machinery.SourceFileLoader("fg_mod", GATE)
_spc = importlib.util.spec_from_loader("fg_mod", _ldr)
FG = importlib.util.module_from_spec(_spc)
_ldr.exec_module(FG)


def compared_files(d):
    out = {}
    for dp, _, fs in os.walk(d):
        if ".old" in dp:
            continue
        for fn in fs:
            if fn in ("Header", "Cell_H", "thermo.dat") or fn.startswith("Cell_D_") or fn.endswith(".csv"):
                out[os.path.relpath(os.path.join(dp, fn), d)] = os.path.join(dp, fn)
    return out


def same_file(rel, a, b):
    ba, bb = open(a, "rb").read(), open(b, "rb").read()
    if ba == bb:
        return True, ""
    if rel.endswith("thermo.dat"):
        ca, cb = FG.thermo_cols(a), FG.thermo_cols(b)
        moved = [h for h in ca if h not in FG.TIMERS and (h not in cb or FG.col_moved(ca[h], cb[h])[0])]
        extra = [h for h in cb if h not in ca]
        if not moved and not extra:
            return True, "timers only"
        return False, "thermo moved %s %s" % (moved[:6], extra[:3])
    return False, "bytes differ"


def compare_dirs(tag, base, new):
    A, B = compared_files(base), compared_files(new)
    only_a, only_b = sorted(set(A) - set(B)), sorted(set(B) - set(A))
    diffs = []
    for rel in sorted(set(A) & set(B)):
        ok, why = same_file(rel, A[rel], B[rel])
        if not ok:
            diffs.append((rel, why))
    ok = not only_a and not only_b and not diffs and len(A) > 0
    say("    list N %-24s %5d files %s%s%s" % (tag, len(A), "byte-identical" if ok else "DIFFERENT",
                                               " (only baseline %d: %s)" % (len(only_a), only_a[:3]) if only_a else "",
                                               " (only new %d: %s)" % (len(only_b), only_b[:3]) if only_b else ""))
    for rel, why in diffs[:6]:
        say("        %s: %s" % (rel, why))
    return ok


def list_n():
    """after 4b only: each output against the outputs of the unchanged code; baselines under unit/<n>/output are copied first
    (the tests rewrite their output/; SR 9)"""
    all_ok = True
    U = os.path.join(T, "unit")
    V = os.path.join(T, "validation")
    # decks into S/legs/N_<name>
    decks = [("zhang_quarter", os.path.join(V, "zhang", "zhang_oglesby", "input"), os.path.join(LEGS, "zhang_w0.01")),
             ("hu_conduction_granite2", os.path.join(V, "hu", "hu_conduction", "input_granite2"), os.path.join(P29, "hu_conduction", "granite2")),
             ("hu_conduction_sandstone2", os.path.join(V, "hu", "hu_conduction", "input_sandstone2"), os.path.join(P29, "hu_conduction", "sandstone2")),
             ("convective_patch", os.path.join(U, "convective_patch", "input"), os.path.join(LEGS, "N_base_convective"))]
    for name, deck, base in decks:
        rc, d = leg("N_" + name, deck)
        ok = rc == 0 and compare_dirs(name, base, d)
        all_ok = all_ok and ok
    # beam, stefan: the deck into its own output/ (baseline copied first); their outputs were written by ccd15b91
    for name in ("beam", "stefan"):
        out = os.path.join(U, name, "output")
        base = os.path.join(HERE, "base_N", name)
        if not os.path.isdir(base):
            shutil.copytree(out, base)
        rc, d = leg("N_" + name, os.path.join(U, name, "input"))
        ok = rc == 0 and compare_dirs(name, base, d)
        all_ok = all_ok and ok
    # the tests that launch their own legs: run them (each rewrites output/), compare with the copy
    for name, args in (("latent_heat", []), ("beam_tau", []), ("face_dirichlet", [])):
        out = os.path.join(U, name, "output")
        base = os.path.join(HERE, "base_N", name)
        if not os.path.isdir(base):
            shutil.copytree(out, base)
        rc = run([PY, os.path.join(U, name, "test")] + args)
        ok = rc == 0 and compare_dirs(name, base, out)
        all_ok = all_ok and ok
    return 0 if all_ok else 1


def b4_l7():
    r1 = subprocess.run([PY, os.path.join(T, "unit", "lateral_flux_lock", "test"), "--check"], cwd=ROOT, capture_output=True, text=True)
    r2 = subprocess.run([PY, os.path.join(T, "unit", "lateral_flux_lock", "test"), "--check", "--plant-patch-def"], cwd=ROOT, capture_output=True, text=True)
    with open(LOG, "a") as f:
        f.write(r1.stdout + r2.stdout)
    l7 = [ln for ln in r1.stdout.splitlines() if "L7 `PatchConstants`" in ln]
    l7p = [ln for ln in r2.stdout.splitlines() if "L7 `PatchConstants`" in ln]
    say("    b4: L7 %s" % (l7[0].strip()[:150] if l7 else "(no clause line)"))
    say("    b4: planted second definition: %s; test rc %d (exactly the L7 clause RED: %s)" % (
        l7p[0].strip()[:90] if l7p else "(no clause line)", r2.returncode, r2.stdout.count("  FAIL  ") == 1))
    return 0 if r1.returncode == 0 and l7 and "ok" in l7[0] and r2.returncode == 1 and r2.stdout.count("  FAIL  ") == 1 and "FAIL" in l7p[0] else 1


def b4_locks():
    rc = run([PY, os.path.join(HERE, "locks.py"), "capture", "stage4b"])
    bad = 1 if rc else 0
    for name, _ in [("retired_keys", 0), ("lateral_flux_lock", 0), ("buckle_lock", 0), ("profile_depth", 0), ("face_heat", 0), ("face_removal", 0)]:
        a = open(os.path.join(HERE, "logs", "stage2", name + ".log")).read().splitlines()
        b = open(os.path.join(HERE, "logs", "stage4b", name + ".log")).read().splitlines()
        ca, cb = collections.Counter(a), collections.Counter(b)
        gone, new = ca - cb, cb - ca
        odd = [ln for ln in new if ln.strip() and "L7" not in ln and "PATCH CONSTANTS" not in ln]
        say("    b4 locks %-20s %4d -> %4d lines; removed %d, added %d, unexplained %d" % (name, len(a), len(b), sum(gone.values()), sum(new.values()), len(odd)))
        for ln in odd[:4] + list(gone)[:4]:
            say("        %s" % ln[:150])
        bad += 1 if gone or odd else 0
    return 1 if bad else 0


def restore_4b():
    for rel in ["MMWSpalling.H"] + [os.path.join("MMWSpalling", n) for n in sorted(os.listdir(os.path.join(PRE_4B, "MMWSpalling")))]:
        shutil.copy2(os.path.join(PRE_4B, rel), os.path.join(SRC_INT, rel))
    return build("p2_30_restored")()


def stage4b():
    rc0, _ = leg("N_base_convective", os.path.join(T, "unit", "convective_patch", "input"))
    say("    b3: convective_patch's baseline (no ccd15b91 output on disk): run on the 4a binary %s, rc %d" % (sha()[:16], rc0))
    rc_b = row("verdict", "b3_build", build("p2_30"))
    row("report", "b3_warnings", no_new_warnings("p2_30", os.path.join(HERE, "build_p2_30a.log")))
    bad = rc_b != 0
    if not bad:
        bad = row("verdict", "b3_gate_run", [PY, GATE, "run", "p2_30"]) != 0
    if not bad:
        bad = row("verdict", "b3_gate_cmp", [PY, GATE, "cmp", "p2_28", "p2_30"]) != 0
        note_tail(r"^face_p2_28 vs|VERDICT", 2)
    if not bad:
        bad = row("verdict", "b3_list_N", list_n) != 0
    if bad:
        say("== b3 RED: the shared values changed an expression's shape; restoring src_pre_4b, rebuilding (4a stands)")
        rc_r = restore_4b()
        say("== restored, rebuilt rc %d binary %s" % (rc_r, sha()[:16]))
        return
    row("verdict", "b4_L7", b4_l7)
    row("verdict", "b4_locks", b4_locks)


def main():
    {"2": stage2, "3": stage3, "4a": stage4a, "4b": stage4b,
     "4a_plant": lambda: row("verdict", "b2_planted", b2_plant),
     "4b_locks": lambda: row("verdict", "b4_locks", b4_locks)}[STAGE]()


if __name__ == "__main__":
    main()
