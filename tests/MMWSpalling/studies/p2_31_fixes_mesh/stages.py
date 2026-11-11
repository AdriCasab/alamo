#!/usr/bin/env python3
"""P2-31: the stages of the driver (check.sh calls this). Runs a stage's Criteria rows in table order; a red GATE row
stops the stage with exit 1; a red verdict / report row is recorded and the stage continues. Output: S/CHECK_<stage>.log
(every command's output) and S/RESULTS.txt (one line per row, and the numbers behind it)."""
import collections
import glob
import hashlib
import importlib.machinery
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
U = os.path.join(T, "unit")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
S30 = os.path.join(T, "studies", "p2_30_structure")
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
GATE = os.path.join(W, "face_gate.py")
ENV = dict(os.environ, EIGEN=os.path.join(ROOT, "ext"), CPLUS_INCLUDE_PATH="/opt/homebrew/include",
           LIBRARY_PATH="/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current")
STAGE = sys.argv[1]
LOG = os.path.join(HERE, "CHECK_%s.log" % STAGE)
RES = os.path.join(HERE, "RESULTS.txt")
LEGS = os.path.join(HERE, "legs")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
sys.path.insert(0, HERE)
import components as C  # noqa: E402


def say(s):
    print(s, flush=True)
    with open(RES, "a") as f:
        f.write(s + "\n")


def log(s):
    with open(LOG, "a") as f:
        f.write(s + "\n")


def run(cmd, env=None, out=None):
    with open(LOG, "a") as lg:
        lg.write("\n$ %s\n" % " ".join(cmd))
        lg.flush()
        r = subprocess.run(cmd, cwd=ROOT, env=env or os.environ, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, errors="replace")
        lg.write(r.stdout)
    if out:
        open(out, "w").write(r.stdout)
    return r.returncode, r.stdout


def row(kind, name, cmd, env=None):
    t0 = time.time()
    rc = cmd() if callable(cmd) else run(cmd, env=env)[0]
    say("  %s: %s (rc %d, %.0f s) [%s]" % (name, "PASS" if rc == 0 else "RED", rc, time.time() - t0, kind))
    if rc != 0 and kind == "gate":
        say("== RED GATE %s: stopped at %s" % (name, time.strftime("%H:%M:%S")))
        sys.exit(1)
    return rc


def sha(p=BIN):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def build(label):
    """make with mmwspalling.cc taken as changed (-W: no file is touched), so every warning is emitted again"""
    def go():
        out = os.path.join(HERE, "build_%s.log" % label)
        rc, _ = run(["make", "-j8", "-W", "src/mmwspalling.cc"], env=ENV, out=out)
        open(os.path.join(HERE, "BINARY_%s.txt" % label), "w").write("%s  bin/mmwspalling-3d-g++\n" % sha())
        say("    binary %s: %s" % (label, sha()[:16]))
        return rc
    return go


def warnings_multiset(path):
    c = collections.Counter()
    for ln in open(path, errors="replace").read().splitlines():
        m = re.search(r"warning: (.*)$", ln)
        if m and "warnings generated" not in ln:
            c[m.group(1).strip()] += 1
    return c


def no_new_warnings(label, base):
    def go():
        new = warnings_multiset(os.path.join(HERE, "build_%s.log" % label))
        old = warnings_multiset(base)
        extra = new - old
        say("    warnings %s: %d; reference %s: %d; new: %d %s" % (label, sum(new.values()), os.path.relpath(base, T),
                                                                  sum(old.values()), sum(extra.values()),
                                                                  [k[:90] for k in extra]))
        return 1 if extra else 0
    return go


def load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    sp = importlib.util.spec_from_loader(name, ld)
    m = importlib.util.module_from_spec(sp)
    ld.exec_module(m)
    return m


def tagged(txt):
    """{tag: ok} from a test's `  ok    <tag> ...` / `  FAIL  <tag> ...` lines"""
    out = {}
    for ln in txt.splitlines():
        m = re.match(r"^\s+(ok|FAIL)\s+(\S+)", ln)
        if m:
            out.setdefault(m.group(2), m.group(1) == "ok")
    return out


# ---- the face gate and list N ------------------------------------------------------------------------------------------
def gate_cmp(a, b, allow=()):
    """face_<a> vs face_<b> component by component: substantive differences only on the legs in `allow`"""
    def go():
        D = C.diff(a, b)
        for lg in sorted(D):
            log(C.describe(lg, D[lg]))
        sub = sorted(lg for lg, e in D.items() if C.substantive(e))
        bad = [lg for lg in sub if lg not in allow]
        n = len(C.md5s(b))
        say("    face_%s vs face_%s: %d files; legs with a substantive difference %s (allowed %s); outside: %s" % (
            a, b, n, sub or "none", list(allow) or "none", bad or "none"))
        for lg in sub:
            say("    " + C.describe(lg, D[lg]).strip()[:400])
        return 1 if bad else 0
    return go


def leg(name, argv_tail, binary=None):
    d = os.path.join(LEGS, name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(LEGS, exist_ok=True)
    b = binary or BIN
    argv = MPI + [b] + list(argv_tail) + ["plot_file=" + d]
    open(d + ".argv.sh", "w").write(" ".join(argv) + "\n")
    open(d + ".BINARY.txt", "w").write(sha(b) + "\n")
    with open(d + ".log", "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, stdout=lg, stderr=subprocess.STDOUT)
    open(d + ".rc", "w").write("%d\n" % rc)
    return rc, d


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
    if open(a, "rb").read() == open(b, "rb").read():
        return True, ""
    if rel.endswith("thermo.dat"):
        ca, cb = C.thermo_cols(a), C.thermo_cols(b)
        moved = [h for h in ca if h not in C.TIMERS and (h not in cb or ca[h] != cb[h])]
        if not moved and set(ca) == set(cb):
            return True, "timers only"
        return False, "thermo moved %s" % moved[:6]
    return False, "bytes differ"


def compare_dirs(tag, base, new):
    A, B = compared_files(base), compared_files(new)
    oa, ob = sorted(set(A) - set(B)), sorted(set(B) - set(A))
    diffs = [(r, w) for r in sorted(set(A) & set(B)) for ok, w in [same_file(r, A[r], B[r])] if not ok]
    ok = not oa and not ob and not diffs and len(A) > 0
    say("    list N %-26s %5d files %s%s%s" % (tag, len(A), "byte-identical" if ok else "DIFFERENT",
                                              " (only base %d)" % len(oa) if oa else "", " (only new %d)" % len(ob) if ob else ""))
    for r, w in diffs[:6]:
        say("        %s: %s" % (r, w))
    return ok


def list_n():
    """P2-30's nine non-gate decks against P2-30's outputs (S30/base_N and S30/legs): byte-identical"""
    ok_all = True
    V = os.path.join(T, "validation")
    P29 = os.path.join(T, "studies", "p2_29_one_heat_path", "legs2")
    for name, deck, base in (("zhang_quarter", os.path.join(V, "zhang", "zhang_oglesby", "input"),
                              os.path.join(S30, "legs", "zhang_w0.01")),
                             ("hu_conduction_granite2", os.path.join(V, "hu", "hu_conduction", "input_granite2"),
                              os.path.join(P29, "hu_conduction", "granite2")),
                             ("hu_conduction_sandstone2", os.path.join(V, "hu", "hu_conduction", "input_sandstone2"),
                              os.path.join(P29, "hu_conduction", "sandstone2")),
                             ("convective_patch", os.path.join(U, "convective_patch", "input"),
                              os.path.join(S30, "legs", "N_base_convective")),
                             ("beam", os.path.join(U, "beam", "input"), os.path.join(S30, "base_N", "beam")),
                             ("stefan", os.path.join(U, "stefan", "input"), os.path.join(S30, "base_N", "stefan"))):
        rc, d = leg("N_%s_%s" % (STAGE, name), [deck])
        ok = rc == 0 and compare_dirs(name, base, d)
        if ok:                    # DISK: byte-identical to the baseline, which stays -- keep the manifest, drop plotfiles
            run([PY, os.path.join(HERE, "prune_legs.py"), d])
        ok_all = ok and ok_all
    for name in ("latent_heat", "beam_tau", "face_dirichlet"):
        rc, _ = run([PY, os.path.join(U, name, "test")])
        ok_all = (rc == 0 and compare_dirs(name, os.path.join(S30, "base_N", name), os.path.join(U, name, "output"))) and ok_all
    return 0 if ok_all else 1


# ---- stage 1 ----------------------------------------------------------------------------------------------------------
PRE1 = os.path.join(HERE, "bins", "pre_1")


def fixture_clauses(binary, out):
    fr = load(os.path.join(U, "face_removal", "test"), "fr_p231_%d" % abs(hash(out)))
    fr.BIN, fr.OUT = binary, out
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    with open(LOG, "a") as lg:
        lg.write("\n# face_removal --p2_31 legs on %s into %s\n" % (sha(binary)[:16], out))
    fr.run_p2_31()
    return fr.check_tref(), fr.check_seedvoid()


def a1_clauses():
    bad = 0
    say("    a1 goal 8's clauses: each PASS on p2_31a (%s), RED on the pre-stage binary bins/pre_1 (%s)" % (
        sha()[:16], sha(PRE1)[:16]))
    (tn, sn) = fixture_clauses(BIN, os.path.join(LEGS, "a1_new"))
    (to, so) = fixture_clauses(PRE1, os.path.join(LEGS, "a1_old"))
    for nm, (okn, ln), (oko, lo) in (("TREF", tn, to), ("SV", sn, so)):
        say("      %-5s new %s: %s" % (nm, "PASS" if okn else "RED", ln))
        say("      %-5s old %s: %s" % (nm, "PASS" if oko else "RED", lo))
        bad += 0 if (okn and not oko) else 1
    il = os.path.join(U, "interlocks", "test")
    rn, tx_n = run([PY, il], env=dict(os.environ, MMW_BIN=BIN))
    ro, tx_o = run([PY, il], env=dict(os.environ, MMW_BIN=PRE1))
    gn, go = tagged(tx_n), tagged(tx_o)
    new_legs, twins = ("N4", "Z1", "S1", "E1"), ("Z0", "E0")
    others = [k for k in gn if k not in new_legs + twins]
    ok_il = rn == 0 and all(gn.get(k) for k in new_legs + twins) and \
        not any(go.get(k, True) for k in new_legs) and all(go.get(k) for k in twins) and all(go.get(k) for k in others)
    say("      interlocks new rc %d (%s); old: %s; old's other %d legs all ok: %s" % (
        rn, " ".join("%s %s" % (k, "ok" if gn.get(k) else "FAIL") for k in new_legs + twins),
        " ".join("%s %s" % (k, "ok" if go.get(k) else "RED") for k in new_legs + twins),
        len(others), all(go.get(k) for k in others)))
    bad += 0 if ok_il else 1
    rk = os.path.join(U, "retired_keys", "test")
    rn2, tk_n = run([PY, rk, "--run"], env=dict(os.environ, MMW_BIN=BIN))
    ro2, tk_o = run([PY, rk, "--run"], env=dict(os.environ, MMW_BIN=PRE1))
    kn, ko = tagged(tk_n), tagged(tk_o)
    red_old = sorted(k for k, v in ko.items() if not v)
    say("      retired_keys --run: new rc %d (%d legs, K31 %s); old: RED on %s (want exactly K31)" % (
        rn2, len(kn), "ok" if kn.get("K31") else "FAIL", red_old))
    bad += 0 if (rn2 == 0 and red_old == ["K31"]) else 1
    return 1 if bad else 0


def a3_gb():
    d_new = os.path.join(W, "face_p2_31a", "GB_S")
    d_old = os.path.join(W, "face_p2_30c", "GB_S")
    rows = []
    for d in (d_old, d_new):
        p = os.path.join(d, "run_removal_events.csv")
        n = max(0, sum(1 for _ in open(p)) - 1) if os.path.exists(p) else -1
        rc = open(os.path.join(d, "rc")).read().strip()
        nan = len(re.findall(r"\bnan\b", open(os.path.join(d, "run.log"), errors="replace").read(), re.I))
        rows.append((rc, n, nan))
    D = C.diff("p2_30c", "p2_31a", legs={"GB_S"})
    e = D.get("GB_S")
    say("    a3 GB_S (T_ref 293.15, IC 820 K): p2_30c rc %s, %d event rows, %d 'nan' in run.log; p2_31a rc %s, %d event "
        "rows, %d 'nan'; differing: %s" % (rows[0] + rows[1] + (C.describe("GB_S", e).strip() if e else "nothing",)))
    return 0 if rows[0][1] == 0 and rows[1][1] == 0 and rows[1][0] == "0" else 1


def stage1():
    row("gate", "a1_build", build("p2_31a"))
    row("gate", "a1_warnings", no_new_warnings("p2_31a", os.path.join(S30, "build_comment_review.log")))
    row("gate", "a1_locks", lambda: run([PY, os.path.join(HERE, "locks.py"), "capture", "stage1"])[0]
        or run([PY, os.path.join(HERE, "locks.py"), "cmp", "pre_1", "stage1"])[0])
    row("gate", "a1_clauses", a1_clauses)
    row("gate", "a2_gate_run", [PY, GATE, "run", "p2_31a"])
    row("gate", "a2_gate_cmp", gate_cmp("p2_30c", "p2_31a", allow=("GB_S",)))
    row("gate", "a2_list_N", list_n)
    row("report", "a3_GB_S", a3_gb)


def locks(a, b, allow=None):
    def go():
        rc = run([PY, os.path.join(HERE, "locks.py"), "capture", b])[0]
        rc2, out = run([PY, os.path.join(HERE, "locks.py"), "cmp", a, b] + ([allow] if allow else []))
        for ln in out.splitlines():
            if ln.startswith("  ") and ("->" in ln or "every lock" in ln or "changed" in ln):
                say("    " + ln.strip()[:160])
        return rc or rc2
    return go


def drop_gate(label, binary_note):
    """DISK (16 GiB free at z0; a gate run writes 3.3 GB): an intermediate gate's output directory is deleted once the
    next comparison has used it. Its md5 list and header list stay (the record); the outputs are regenerable exactly:
    face_gate.py run <label> on the kept binary named here."""
    def go():
        d = os.path.join(W, "face_" + label)
        if not os.path.isdir(d):
            say("    disk: face_%s already absent" % label)
            return 0
        n = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(d) for f in fs)
        shutil.rmtree(d)
        say("    disk: deleted witness/face_%s/ (%.2f GB); kept face_%s.md5 and face_%s_headers.txt; regenerable with "
            "face_gate.py run %s on %s" % (label, n / 1e9, label, label, label, binary_note))
        return 0
    return go


# ---- stage 2a: dead code out ------------------------------------------------------------------------------------------
PRE2 = os.path.join(HERE, "bins", "pre_2")
PRE3 = os.path.join(HERE, "bins", "pre_3")


def b1_rk():
    rk = os.path.join(U, "retired_keys", "test")
    rn, tx = run([PY, rk])
    ro, txo = run([PY, rk, "--run"], env=dict(os.environ, MMW_BIN=PRE2))
    ko = tagged(txo)
    red = sorted(k for k, v in ko.items() if not v)
    r13 = [ln.strip() for ln in tx.splitlines() if " R13 " in ln]
    mut = [ln.strip() for ln in tx.splitlines() if "microstructure_enabled = false;` back" in ln]
    kk = {k: v for k, v in tagged(tx).items() if k in ("K61", "K62", "K63")}
    say("    b1 retired_keys (both halves) on p2_31b0: rc %d; %s; mutation: %s; K61-K63 %s" % (
        rn, r13[0][:150] if r13 else "(no R13 line)", mut[0][-40:] if mut else "(none)", kk))
    say("    b1 retired_keys --run on the pre-stage binary (bins/pre_2): RED on %s (want exactly K62)" % red)
    return 0 if rn == 0 and r13 and r13[0].startswith("ok") and mut and "goes RED" in mut[0] and red == ["K62"] else 1


def stage2a():
    row("gate", "b1_build", build("p2_31b0"))
    row("gate", "b1_warnings", no_new_warnings("p2_31b0", os.path.join(HERE, "build_p2_31a.log")))
    row("gate", "b1_locks", locks("stage1", "stage2a", r"R13|microstructure_enabled"))
    row("gate", "b1_retired_keys", b1_rk)
    say("    b1: the touched decks parse -- 2a touches no deck (Z0 L9: no live deck sets a vapour key)")
    row("gate", "b2_gate_run", [PY, GATE, "run", "p2_31b0"])
    row("gate", "b2_gate_cmp", gate_cmp("p2_31a", "p2_31b0"))
    row("report", "b2_vs_p2_30c", gate_cmp("p2_30c", "p2_31b0", allow=("GB_S",)))
    row("gate", "b2_list_N", list_n)
    row("report", "disk", drop_gate("p2_31a", "bins/pre_2 (p2_31a's binary)"))


# ---- stage 2b: regime_field out ---------------------------------------------------------------------------------------
def b3_audit():
    D = C.diff("p2_31b0", "p2_31b")
    bad, lost = [], []
    for lg, e in sorted(D.items()):
        extra = (e["plot"] or e["plot_new"] or (e["plot_gone"] - {"regime_field"}) or e["thermo"] or e["thermo_gone"]
                 or e["thermo_new"] or e["csv"] or e["other"] or e["only_a"] or e["only_b"])
        if extra:
            bad.append(lg)
            say("    " + C.describe(lg, e).strip()[:300])
        if e["plot_gone"]:
            lost.append(lg)
    say("    b3 p2_31b0 -> p2_31b: legs whose plotfiles lost exactly regime_field: %d %s; every other component, every "
        "thermo column and every CSV identical: %s" % (len(lost), lost, "yes" if not bad else "NO on %s" % bad))
    return 1 if bad or not lost else 0


def stage2b():
    row("gate", "b3_build", build("p2_31b"))
    row("gate", "b3_warnings", no_new_warnings("p2_31b", os.path.join(HERE, "build_p2_31b0.log")))
    row("gate", "b3_locks", locks("stage2a", "stage2b"))
    row("gate", "b3_gate_run", [PY, GATE, "run", "p2_31b"])
    row("gate", "b3_audit", b3_audit)
    row("report", "disk", drop_gate("p2_31b0", "bins/pre_2b (p2_31b0's binary)"))


# ---- stage 3: the overlap rule ----------------------------------------------------------------------------------------
L10_PLOT = {"phi", "surface", "spall_thickness", "RoP_spall", "RoP"}
L10_THERMO = {"ledger_surface_missing_cols"}
OVL_LEGS = ("fr/TL45_1", "fr/SPH_1", "fr/BORE_1")


def c1_ovl():
    bad = 0
    res = {}
    for tag, b in (("new", BIN), ("old", PRE3)):
        fr = load(os.path.join(U, "face_removal", "test"), "fr_ovl_" + tag)
        out = os.path.join(LEGS, "c1_" + tag)
        if os.path.isdir(out):
            shutil.rmtree(out)
        os.makedirs(out)
        fr.BIN, fr.OUT = b, out
        fr.run_ovl()
        res[tag] = fr.check_ovl()
        say("      OVL %s (%s) %s: %s" % (tag, sha(b)[:16], "PASS" if res[tag][0] else "RED", res[tag][1]))
    return 0 if res["new"][0] and not res["old"][0] else 1


def c2_cmp():
    D = C.diff("p2_31b", "p2_31")
    bad = []
    for lg, e in sorted(D.items()):
        if not C.substantive(e):
            continue
        say("    " + C.describe(lg, e).strip()[:400])
        if lg not in OVL_LEGS:
            bad.append(lg)
            continue
        if (set(e["plot"]) - L10_PLOT) or e["plot_gone"] or e["plot_new"] or (set(e["thermo"]) - L10_THERMO) \
                or e["thermo_gone"] or e["thermo_new"] or e["csv"] or e["other"] or e["only_a"] or e["only_b"]:
            bad.append(lg)
    say("    c2 p2_31b -> p2_31: outside the three overlap legs, substantive on %s; inside them, outside L10's set "
        "(plot %s, thermo %s): %s" % ([lg for lg in bad if lg not in OVL_LEGS] or "none", sorted(L10_PLOT),
                                     sorted(L10_THERMO), [lg for lg in bad if lg in OVL_LEGS] or "none"))
    return 1 if bad else 0


def c3_cells():
    for lg in OVL_LEGS:
        da, db = os.path.join(W, "face_p2_31b", lg, "run"), os.path.join(W, "face_p2_31", lg, "run")
        n_cells, n_pf = 0, 0
        for pf in sorted(glob.glob(os.path.join(da, "[0-9]" * 5 + "cell"))):
            pb = os.path.join(db, os.path.basename(pf))
            _, no = C.AC.read_header(pf)
            _, nn = C.AC.read_header(pb)
            fo, fn = C.AC.read_cell_h(pf), C.AC.read_cell_h(pb)
            co, cn, moved = {}, {}, 0
            for (a_f, a_o), (b_f, b_o) in zip(fo, fn):
                co.setdefault(a_f, C.AC.read_fabs(pf, a_f))
                cn.setdefault(b_f, C.AC.read_fabs(pb, b_f))
                A, B = co[a_f][a_o][2], cn[b_f][b_o][2]
                mask = None
                for c in L10_PLOT & set(no):
                    m = A[no.index(c)].view("u8") != B[nn.index(c)].view("u8")
                    mask = m if mask is None else (mask | m)
                moved += int(mask.sum()) if mask is not None else 0
            n_cells += moved
            n_pf += 1 if moved else 0
        th = C.thermo_cols(os.path.join(db, "thermo.dat"))
        ov = th.get("fr_overlap", [])
        say("    c3 %-11s overlap census (fr_overlap): %g cells in %d steps; plotfile cells moved in L10's components: %d "
            "over %d plotfiles" % (lg, sum(ov), sum(1 for x in ov if x > 0), n_cells, n_pf))
    return 0


def stage3():
    row("gate", "c1_build", build("p2_31"))
    row("gate", "c1_warnings", no_new_warnings("p2_31", os.path.join(HERE, "build_p2_31b.log")))
    row("gate", "c1_locks", locks("stage2b", "stage3"))
    row("gate", "c1_OVL", c1_ovl)
    row("gate", "c2_gate_run", [PY, GATE, "run", "p2_31"])
    row("gate", "c2_gate_cmp", c2_cmp)
    row("gate", "c2_list_N", list_n)
    row("report", "c3_cells", c3_cells)
    row("report", "disk", drop_gate("p2_31b", "bins/pre_3 (p2_31b's binary)"))


# ---- stage 4: the mesh test -------------------------------------------------------------------------------------------
MESH = os.path.join(HERE, "mesh.py")


def m_reduce():
    rc, out = run([PY, MESH, "reduce"], out=os.path.join(HERE, "MESH_RESULTS.txt"))
    for ln in out.splitlines():
        say(ln)
    return rc


def m_conditional():
    txt = open(os.path.join(HERE, "MESH_RESULTS.txt")).read()
    todo = re.findall(r"=> run (R\d+_08) ", txt)
    say("    m8 conditional fine legs to run (a 1.5625 mm replay that did not reproduce P2-25): %s" % (todo or "none"))
    if not todo:
        return 0
    return run([PY, MESH, "run"] + todo)[0]


def stage4():
    row("gate", "m_keys", [PY, MESH, "keys"])
    row("verdict", "m0_run", [PY, MESH, "run"])
    row("report", "m_reduce", m_reduce)
    if row("report", "m8_conditional", m_conditional) == 0 and \
            re.search(r"=> run R\d+_08 ", open(os.path.join(HERE, "MESH_RESULTS.txt")).read()):
        row("report", "m_reduce_after_conditional", m_reduce)


def main():
    {"1": stage1, "2a": stage2a, "2b": stage2b, "3": stage3, "4": stage4}[STAGE]()


if __name__ == "__main__":
    main()
