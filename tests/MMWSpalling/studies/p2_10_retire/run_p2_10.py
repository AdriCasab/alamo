#!/usr/bin/env python3
"""P2-10 -- the references and gates of the deletion step (ACTIVE_STEP.md `## Criteria`).

MODES
  refs <label>     M0T + the witness set W into W_<label>/ on the CURRENT binary;
                   writes W_<label>.md5 (every output file, path relative to
                   W_<label>/) and W_<label>_headers.txt (each leg's thermo header).
                     W = unit/face_removal T00_1_k0 T00_1_k1 TL35_1 TL45_1 SPH_1 BORE_1
                         (P2-9's fixture and keys) + unit/face_list CAV STAIR OFF0 CAVU
                         (its trimmed input) + M0T
  units <label>    every KEPT unit test (`test`, no flags) -> UNIT_<label>/<name>.log,
                   UNIT_<label>/rc.txt. The test's own output/ is archived (cp -c,
                   an APFS clone) to the scratchpad first (SR 9).
  a1               M0T (W_bB) vs studies/p2_9_face_removal/legs/M0 on disk.
  cmp <label>      W_<label> vs W_bB: per file, byte-identical; thermo.dat column by
                   column (columns only in W_bB listed with max |value|).
  units_cmp <label> UNIT_<label> vs UNIT_bB: rc and the PASS/FAIL lines.
  retired          a7: the six retired-key legs (M0T + one token, stop_time 1e-9) and
                   their rc + the log lines naming the key -> RETIRED.txt.

Every leg writes argv.sh and the binary's sha (SR 8); runs are serial, -np 4, and the
harness holds the host awake itself (SR 11).
"""
import glob
import hashlib
import importlib.machinery
import importlib.util
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
M0 = os.path.join(T, "studies", "p2_9_face_removal", "legs", "M0")
SCRATCH = "/private/tmp/claude-501/-Users-tzetze20-amr-tools-alamo/ac0a0f5a-4af6-48bf-b61d-9009cf0c8681/scratchpad"

M0T_DELETE = [  # the packet's 13 tokens, by key
    "surface_patch.side_face_flux", "surface_patch.side_face_step_height",
    "surface_patch.side_face_area_law", "surface_patch.side_face_h", "surface_patch.side_face_T",
    "surface_patch.side_face_staircase", "spall.lateral_drive", "spall.lateral_removal",
    "spall.lateral_subcell", "spall.lateral_thermal", "surface_patch.side_face_h_datum",
    "surface_patch.side_face_h_r_min", "surface_patch.side_face_h_law"]
RETIRED_LEGS = {
    "R_side_face_flux": "surface_patch.side_face_flux=1",
    "R_side_face_h": "surface_patch.side_face_h=169.3",
    "R_lateral_drive": "spall.lateral_drive=1",
    "R_per_cell_removal": "spall.per_cell_removal=1",
    "R_jet_wall_march": "surface_patch.jet_wall_march=1",
    "R_jet_body": "surface_patch.jet_body=1",
}
DELETED_UNITS = {"side_face_flux", "side_face_area", "side_face_step_height", "side_staircase",
                 "lateral_drive", "wall_march"}
FR_LEGS = ["T00_1_k0", "T00_1_k1", "TL35_1", "TL45_1", "SPH_1", "BORE_1"]
FL_LEGS = ["CAV", "STAIR", "OFF0", "CAVU"]
# wall-clock timers (P2-9 rev 2 item 7): non-deterministic by construction, never compared
TIMERS = ("fr_t_removal", "fr_t_step")
# the face list's event log appends these under spall.lateral_drive (FaceRemoval.H)
LAT_CSV = ["riser_depth", "ell_lat", "T_face_lat", "K_I_lat", "K_Ic_lat", "q_lat", "sigma_lat",
           "h_f_lat", "h_f_unc_lat", "lat_accum", "lat_nfire"]


def load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


def sha(path, n=24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def m0t_argv(outdir, extra=()):
    toks = shlex.split(open(os.path.join(M0, "argv.sh")).read())
    out, gone = [], []
    for t in toks:
        k = t.split("=", 1)[0] if "=" in t else None
        if k in M0T_DELETE:
            gone.append(t)
            continue
        out.append("plot_file=" + os.path.join(outdir, "run") if k == "plot_file" else t)
    assert sorted(t.split("=", 1)[0] for t in gone) == sorted(M0T_DELETE), gone
    return out + list(extra), gone


def run_cmd(d, cmd):
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    open(os.path.join(d, "BINARY.txt"), "w").write("bin/mmwspalling-3d-g++ sha256 %s\n" % sha(BIN, 64))
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as f:
        rc = subprocess.call(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    if rc == 0:
        open(os.path.join(d, ".done"), "w").write("rc 0 wall %.1f\n" % (time.time() - t0))
    print("  %-22s rc=%d wall %7.1f s" % (os.path.relpath(d, HERE), rc, time.time() - t0), flush=True)
    return rc


def hash_tree(top, out_md5):
    rows = []
    for dp, _, fs in os.walk(top):
        for fn in fs:
            p = os.path.join(dp, fn)
            rel = os.path.relpath(p, top)
            if fn in ("run.log", "argv.sh", "BINARY.txt", ".done"):
                continue
            rows.append("%s  %s" % (md5(p), rel))
    rows.sort(key=lambda r: r.split("  ", 1)[1])
    open(out_md5, "w").write("\n".join(rows) + "\n")
    return len(rows)


# ------------------------------------------------------------------ refs ------
def refs(label):
    W = os.path.join(HERE, "W_" + label)
    os.makedirs(W, exist_ok=True)
    # M0T
    d = os.path.join(W, "M0T")
    cmd, gone = m0t_argv(d)
    rc = run_cmd(d, cmd)
    open(os.path.join(d, "DROPPED.txt"), "w").write(
        "M0T: studies/p2_9_face_removal/legs/M0/argv.sh (sha256 %s) token for token, plot_file "
        "re-pointed, these %d tokens DELETED, nothing added:\n%s\n"
        % (sha(os.path.join(M0, "argv.sh")), len(gone), "\n".join("   " + g for g in gone)))
    # unit/face_removal legs, P2-9's fixture and keys
    fr = load(os.path.join(T, "unit", "face_removal", "test"), "fr_p2_10")
    fr.OUT = os.path.join(W, "fr")
    os.makedirs(fr.OUT, exist_ok=True)
    k1 = {"surface_patch.face_list": "1"}
    for nm in FR_LEGS:
        t0 = time.time()
        if nm.startswith("T00"):
            rc = fr.run_leg(nm, dict(fr.mesh_keys("1"), **{"surface_patch.face_list": nm[-1]}))
        elif nm.startswith("TL"):
            th, ph = fr.TILTS[nm.split("_")[0]]
            rc = fr.run_leg(nm, dict(fr.mesh_keys("1"), **k1), seeds=fr.tilt_seeds(th, ph, "1"))
        else:
            c = nm.split("_")[0]
            rc = fr.run_leg(nm, dict(fr.MI.curved_keys(c, 0.001), **k1, **{"stop_time": "20.0"}),
                            seeds=fr.MI.curved_seeds(c, 0.001))
        dd = os.path.join(fr.OUT, nm)
        open(os.path.join(dd, "BINARY.txt"), "w").write("bin/mmwspalling-3d-g++ sha256 %s\n" % sha(BIN, 64))
        if rc == 0:
            open(os.path.join(dd, ".done"), "w").write("rc 0 wall %.1f\n" % (time.time() - t0))
    # unit/face_list legs, its (trimmed) input
    fl = load(os.path.join(T, "unit", "face_list", "test"), "fl_p2_10")
    fl.OUT = os.path.join(W, "fl")
    os.makedirs(fl.OUT, exist_ok=True)
    for nm in FL_LEGS:
        rc = fl.run(nm)
        print("  fl/%-18s rc=%d" % (nm, rc), flush=True)
        open(os.path.join(fl.OUT, nm, "BINARY.txt"), "w").write("bin/mmwspalling-3d-g++ sha256 %s\n" % sha(BIN, 64))
    n = hash_tree(W, os.path.join(HERE, "W_%s.md5" % label))
    with open(os.path.join(HERE, "W_%s_headers.txt" % label), "w") as f:
        for p in sorted(glob.glob(os.path.join(W, "**", "thermo.dat"), recursive=True)):
            f.write("%s\t%s" % (os.path.relpath(p, W), open(p).readline()))
    print("refs %s: hashed %d files -> W_%s.md5" % (label, n, label))


# ------------------------------------------------------------------ units -----
def kept_units():
    return sorted(d for d in os.listdir(os.path.join(T, "unit"))
                  if os.path.isfile(os.path.join(T, "unit", d, "test")) and d not in DELETED_UNITS)


# check-only tests: `test` re-reads tests/.../<name>/output, so the input runs first
CHECK_ONLY = ["alpha_beta_transition", "amr_microstructure_regrid", "beam", "convective_patch",
              "equilibrium", "grain_topology", "heterogeneous_enthalpy", "heterogeneous_kappa",
              "stefan", "thermal_stress", "voronoi"]


def units(label, only=None):
    """every kept unit test on the current binary. A check-only test (CHECK_ONLY)
    first runs its `input` (-np 4, cwd = repo root, the input's own plot_file).
    Each test's output* dirs are then hashed to UNIT_<label>/<name>.md5 (AMReX
    `metadata` excluded), so a8 compares files as well as verdict lines."""
    U = os.path.join(HERE, "UNIT_" + label)
    os.makedirs(U, exist_ok=True)
    arch = os.path.join(SCRATCH, "unit_output_pre_%s" % label)
    rcp = os.path.join(U, "rc.txt")
    rcs = {}
    if os.path.exists(rcp):
        for l in open(rcp).read().splitlines()[1:]:
            rcs[l.split()[0]] = l
    for name in kept_units():
        if only and name not in only:
            continue
        src = os.path.join(T, "unit", name)
        if name != "face_removal":      # its output/ is P2-9's record and no-flag runs nothing
            for o in glob.glob(os.path.join(src, "output*")):
                dst = os.path.join(arch, name, os.path.basename(o))
                if not os.path.exists(dst):
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    subprocess.call(["cp", "-cR", o, dst])
        t0 = time.time()
        pre = ""
        if name in CHECK_ONLY:
            with open(os.path.join(U, name + ".run.log"), "w") as f:
                rr = subprocess.call(["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4",
                                      BIN, os.path.join(src, "input")],
                                     cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
            pre = " input rc=%d" % rr
        with open(os.path.join(U, name + ".log"), "w") as f:
            rc = subprocess.call([PY, os.path.join(src, "test")], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        rows = []
        if name != "face_removal":
            for o in sorted(glob.glob(os.path.join(src, "output*"))):
                for dp, _, fs in os.walk(o):
                    for fn in fs:
                        if fn == "metadata":
                            continue
                        pth = os.path.join(dp, fn)
                        rows.append("%s  %s" % (md5(pth), os.path.relpath(pth, src)))
        open(os.path.join(U, name + ".md5"), "w").write("\n".join(sorted(rows, key=lambda r: r[34:])) + "\n")
        rcs[name] = "%-28s rc=%d  %7.1f s%s" % (name, rc, time.time() - t0, pre)
        print("  " + rcs[name], flush=True)
        open(rcp, "w").write("binary sha256 %s\n%s\n" % (sha(BIN, 64), "\n".join(rcs[k] for k in sorted(rcs))))


VERDICT = re.compile(r"\b(PASS|FAIL|OK|ERROR|ABORT)\b")


def verdict_lines(path):
    out = []
    for ln in open(path, errors="replace"):
        if VERDICT.search(ln):
            # strip wall times and paths so only the verdicts and the measured values compare
            ln = re.sub(r"\b\d+(\.\d+)?\s*s\b(?=[^a-zA-Z]|$)", "<t>s", ln.rstrip())
            out.append(ln)
    return out


def units_cmp(label):
    A, B = os.path.join(HERE, "UNIT_bB"), os.path.join(HERE, "UNIT_" + label)
    rcA = {l.split()[0]: l.split()[1] for l in open(os.path.join(A, "rc.txt")).read().splitlines()[1:]}
    rcB = {l.split()[0]: l.split()[1] for l in open(os.path.join(B, "rc.txt")).read().splitlines()[1:]}
    rows = []
    for n in sorted(set(rcA) | set(rcB)):
        a = verdict_lines(os.path.join(A, n + ".log")) if n in rcA else None
        b = verdict_lines(os.path.join(B, n + ".log")) if n in rcB else None
        same = a == b
        def md5s(D):
            p = os.path.join(D, n + ".md5")
            return dict((l[34:], l[:32]) for l in open(p).read().splitlines() if l) if os.path.exists(p) else {}
        ma, mb = md5s(A), md5s(B)
        nd = sum(1 for k in set(ma) | set(mb) if ma.get(k) != mb.get(k))
        rows.append("%-28s rc %s -> %s   verdict lines %s   output files %d, differing %d" % (
            n, rcA.get(n, "-"), rcB.get(n, "-"),
            "identical (%d)" % len(a) if same else "DIFFER", len(mb), nd))
        if nd:
            for k in sorted(set(ma) | set(mb)):
                if ma.get(k) != mb.get(k):
                    rows.append("      file differs: " + k)
                    if sum(1 for r in rows if r.startswith("      file differs")) > 400:
                        break
        if not same and a is not None and b is not None:
            import difflib
            for d in list(difflib.unified_diff(a, b, lineterm="", n=0))[2:40]:
                rows.append("      " + d)
    # a8: same rc and identical verdict lines on every kept test; the two the packet
    # edits (face_list's input, lateral_flux_lock's checks) are reported, not gated
    DECL = ("lateral_flux_lock",)
    bad = []
    for n in sorted(set(rcA) | set(rcB)):
        if n in DECL:
            continue
        a = verdict_lines(os.path.join(A, n + ".log")) if n in rcA else None
        b = verdict_lines(os.path.join(B, n + ".log")) if n in rcB else None
        if rcA.get(n) != rcB.get(n) or a != b:
            bad.append(n)
    rows.append("a8 (rc + verdict lines, %s reported only): %s" % (
        ", ".join(DECL), "PASS" if not bad else "FAIL: %s" % ", ".join(bad)))
    if "lateral_flux_lock" in rcB:
        rows.append("lateral_flux_lock on %s: rc %s, last line: %s" % (label, rcB["lateral_flux_lock"], (
            open(os.path.join(B, "lateral_flux_lock.log")).read().strip().splitlines() or ["-"])[-1][:100]))
    txt = "\n".join(rows)
    print(txt)
    open(os.path.join(HERE, "UNITS_%s_vs_bB.txt" % label), "w").write(txt + "\n")


# ------------------------------------------------------------------ compare ---
def read_thermo(p):
    L = open(p).read().splitlines()
    hdr = [h.strip() for h in L[0].split("\t")]
    rows = [r.split("\t") for r in L[1:]]
    return hdr, rows


def thermo_cmp(pa, pb):
    """common columns token-identical; columns only in a, with max |value|."""
    ha, ra = read_thermo(pa)
    hb, rb = read_thermo(pb)
    out = {"rows": (len(ra), len(rb)), "only_a": {}, "only_b": [c for c in hb if c not in ha],
           "diff_cols": []}
    if len(ra) != len(rb):
        return out
    ia = {c: i for i, c in enumerate(ha)}
    ib = {c: i for i, c in enumerate(hb)}
    for c in ha:
        if c not in ib:
            vals = [abs(float(r[ia[c]])) for r in ra if len(r) > ia[c] and r[ia[c]].strip()]
            out["only_a"][c] = max(vals) if vals else 0.0
            continue
        if c in TIMERS:
            continue
        if any(x[ia[c]] != y[ib[c]] for x, y in zip(ra, rb)):
            out["diff_cols"].append(c)
    # order of the common columns
    out["order_same"] = [c for c in ha if c in ib] == [c for c in hb if c in ia]
    return out


def tree_files(top):
    r = {}
    for dp, _, fs in os.walk(top):
        for fn in fs:
            # `metadata` is AMReX's run record (argv, date, wall time) and `diff.patch`
            # the build's own git diff (GIT_DIFF_PATCH_OUTPUT): both differ between any
            # two builds/runs and carry no result (declared; hash.sh excludes them too)
            if fn in ("run.log", "argv.sh", "BINARY.txt", ".done", "DROPPED.txt", "metadata",
                      "diff.patch"):
                continue
            p = os.path.join(dp, fn)
            r[os.path.relpath(p, top)] = p
    return r


def compare_dirs(a, b, lines, tag):
    fa, fb = tree_files(a), tree_files(b)
    only_a = sorted(set(fa) - set(fb))
    only_b = sorted(set(fb) - set(fa))
    diff, same, th = [], 0, []
    for rel in sorted(set(fa) & set(fb)):
        if os.path.basename(rel) == "thermo.dat":
            th.append(rel)
            continue
        if md5(fa[rel]) == md5(fb[rel]):
            same += 1
        else:
            diff.append(rel)
    lines.append("%s: %d files byte-identical, %d differ, %d only in ref, %d only in new"
                 % (tag, same, len(diff), len(only_a), len(only_b)))
    for rel in diff[:20]:
        lines.append("     DIFFERS  " + rel)
    for rel in only_a[:10]:
        lines.append("     ONLY REF " + rel)
    for rel in only_b[:10]:
        lines.append("     ONLY NEW " + rel)
    ok = not diff and not only_a and not only_b
    for rel in th:
        c = thermo_cmp(fa[rel], fb[rel])
        byte = md5(fa[rel]) == md5(fb[rel])
        lines.append("   %s: rows %s, byte-identical %s, common columns differing %s, order same %s, "
                     "only in new %s" % (rel, c["rows"], byte, c["diff_cols"] or "none",
                                         c.get("order_same"), c["only_b"] or "none"))
        if c["only_a"]:
            nz = {k: v for k, v in c["only_a"].items() if v != 0.0}
            lines.append("       %d columns only in ref; non-zero among them: %s"
                         % (len(c["only_a"]), nz or "none"))
        ok = ok and c["rows"][0] == c["rows"][1] and not c["diff_cols"] and not c["only_b"] \
            and c.get("order_same", False)
    return ok, th


def a1():
    lines = ["a1 -- binary B: M0T (W_bB/M0T) vs M0 on disk (studies/p2_9_face_removal/legs/M0)", "=" * 100]
    ok, th = compare_dirs(M0, os.path.join(HERE, "W_bB", "M0T"), lines, "M0 vs M0T")
    c = thermo_cmp(os.path.join(M0, "run", "thermo.dat"), os.path.join(HERE, "W_bB", "M0T", "run", "thermo.dat"))
    lines.append("   columns only in M0 (%d): %s" % (len(c["only_a"]), " ".join(c["only_a"])))
    la = open(os.path.join(M0, "run_removal_events.csv")).read().splitlines()
    lb = open(os.path.join(HERE, "W_bB", "M0T", "run_removal_events.csv")).read().splitlines()
    ha = la[0].split(",")
    keep = [i for i, h in enumerate(ha) if h not in LAT_CSV]
    hdr_ok = [ha[i] for i in keep] == lb[0].split(",")
    nbad = sum(1 for x, y in zip(la[1:], lb[1:])
               if ",".join(x.split(",")[i] for i in keep) != y) + abs(len(la) - len(lb))
    nz = sum(1 for x in la[1:] for i, h in enumerate(ha) if h in LAT_CSV and float(x.split(",")[i]) != 0.0)
    lines.append("   run_removal_events.csv: M0 carries %d lateral columns (FaceRemoval.H `lat_csv = "
                 "sp_lateral_drive`); with them stripped the header matches %s and %d of %d rows "
                 "differ; non-zero lateral cells %d" % (len(ha) - len(keep), hdr_ok, nbad, len(lb) - 1, nz))
    strict = ok
    sub = hdr_ok and nbad == 0 and nz == 0 and all(
        d.strip().endswith("run_removal_events.csv") for d in lines if "DIFFERS" in d)
    lines.append("a1  as written (every file but thermo byte-identical): %s; on substance (event log "
                 "identical once the 11 all-zero lateral columns M0's key appends are stripped; every "
                 "plotfile identical; thermo common columns identical but the two wall-clock timers; "
                 "AMReX `metadata` excluded): %s" % ("PASS" if strict else "FAIL", "PASS" if sub else "FAIL"))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "GATE_a1.txt"), "w").write(txt + "\n")


# the thermo columns each binary deletes (DELETION_MAP.txt), cumulative
DECLARED = {"c1": {"fl_lat_top", "fl_lat_nfb", "lat_len_fb_h0", "lat_len_fb_q0"}}
DECLARED["c2"] = DECLARED["c1"] | {"surf_side_P_top", "surf_side_P_buried", "surf_side_P_rin",
                                   "surf_side_P_rout"}


def cmp(label):
    lines = ["W_%s vs W_bB (a4/a6): every file byte-identical but thermo.dat; thermo: common columns "
             "identical (timers excluded), no new column, every dropped column DECLARED and identically "
             "0 in the reference; AMReX `metadata` excluded" % label, "=" * 100]
    ok = True
    A, B = os.path.join(HERE, "W_bB"), os.path.join(HERE, "W_" + label)
    dropped = {}
    for sub in ["M0T"] + ["fr/" + x for x in FR_LEGS] + ["fl/" + x for x in FL_LEGS]:
        o, th = compare_dirs(os.path.join(A, sub), os.path.join(B, sub), lines, sub)
        ok = ok and o
        for rel in th:
            c = thermo_cmp(os.path.join(A, sub, rel), os.path.join(B, sub, rel))
            for col, mx in c["only_a"].items():
                dropped.setdefault(col, []).append((sub, mx))
    bad = [(c, v) for c, v in dropped.items()
           if c not in DECLARED[label] or any(mx != 0.0 for _, mx in v)]
    lines.append("dropped columns: %s" % ", ".join("%s (%d legs, max |v| %g)" % (
        c, len(v), max(mx for _, mx in v)) for c, v in sorted(dropped.items())) or "none")
    lines.append("undeclared or non-zero dropped columns: %s" % (bad or "none"))
    ok = ok and not bad
    lines.append("W_%s vs W_bB: %s" % (label, "PASS" if ok else "FAIL"))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "W_%s_vs_bB.txt" % label), "w").write(txt + "\n")


def retired():
    R = os.path.join(HERE, "retired")
    out = ["a7 -- each retired key given on M0T (stop_time 1e-9): rc and the log lines naming it", "=" * 100]
    for nm, tok in RETIRED_LEGS.items():
        d = os.path.join(R, nm)
        cmd, _ = m0t_argv(d, [tok, "stop_time=1e-9"])
        rc = run_cmd(d, cmd)
        key = tok.split("=", 1)[0]
        hits = [l.strip() for l in open(os.path.join(d, "run.log"), errors="replace")
                if key in l and "P2-10" in l]
        out.append("  %-20s %-40s rc=%d  %s" % (nm, tok, rc, ("log: " + hits[0][:160]) if hits else
                                                    "NO log line names the key and P2-10"))
        out.append("      %s" % ("PASS" if rc != 0 and hits else "FAIL"))
    txt = "\n".join(out)
    print(txt)
    open(os.path.join(HERE, "RETIRED.txt"), "w").write(txt + "\n")


def a9():
    """lines removed per file (src_bB -> src/), and the src/ grep for deleted names."""
    L = ["a9 -- lines removed per file (binary B -> the shipped source) and deleted names left in src/",
         "=" * 100]
    tot = 0
    for f, p in (("MMWSpalling.H", "src/Integrator/MMWSpalling.H"),
                 ("Removal.H", "src/Integrator/MMWSpalling/Removal.H"),
                 ("FaceRemoval.H", "src/Integrator/MMWSpalling/FaceRemoval.H"),
                 ("HeatConduction.H", "src/Integrator/HeatConduction.H")):
        a = len(open(os.path.join(HERE, "src_bB", f)).read().split("\n"))
        b = len(open(os.path.join(ROOT, p)).read().split("\n"))
        tot += a - b
        L.append("  %-18s %6d -> %6d   removed %5d" % (f, a, b, a - b))
    L.append("  %-18s %22s   removed %5d" % ("total", "", tot))
    sys.path.insert(0, HERE)
    import retire_tests as RT
    pat = re.compile(RT.DELETED_NAMES)
    left = []
    for dp, _d, fs in os.walk(os.path.join(ROOT, "src")):
        for fn in fs:
            if not fn.endswith((".H", ".h", ".cpp", ".cc")):
                continue
            pth = os.path.join(dp, fn)
            inl = False
            for ln, line in enumerate(open(pth, errors="replace"), 1):
                if "P2-10: THE RETIRED KEYS" in line:
                    inl = True
                if inl and "remove it from the input" in line:
                    inl = False
                    continue
                if not inl and pat.search(line):
                    left.append("%s:%d %s" % (os.path.relpath(pth, ROOT), ln, line.strip()[:80]))
    L.append("deleted names left in src/ outside the retired-key list: %d" % len(left))
    L += ["  " + x for x in left]
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "A9.txt"), "w").write(txt + "\n")


def main():
    a = sys.argv[1:]
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    if a[0] == "ord1":                  # the packet's ord 1 on binary B, one launch
        refs("bB")
        a1()
        units("bB")
        print("ORD1 DONE", flush=True)
    elif a[0] == "refs":
        refs(a[1])
    elif a[0] == "units":
        units(a[1], a[2].split(",") if len(a) > 2 else None)
    elif a[0] == "units_cmp":
        units_cmp(a[1])
    elif a[0] == "a1":
        a1()
    elif a[0] == "cmp":
        cmp(a[1])
    elif a[0] == "retired":
        retired()
    elif a[0] == "a9":
        a9()
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
