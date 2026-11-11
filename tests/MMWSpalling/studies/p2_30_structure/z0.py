#!/usr/bin/env python3
"""P2-30 z0 (offline) -- the packet's ledger L1, L4-L6, L9, the git state, the baselines on disk. L2 is pp_identity.py --selftest."""
import collections
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
M = open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")).read().split("\n")
FAIL = []


def need(ok, msg):
    print(("  ok    " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def live_tests():
    out = []
    for pat in ("unit/*/test", "unit/*/*.py", "validation/**/test*", "validation/**/*.py", "extensions/**/test*",
                "studies/tree_sweep_0922/witness/*.py"):
        out += [p for p in glob.glob(os.path.join(T, pat), recursive=True) if os.path.isfile(p)
                and "/output" not in p and "/reference" not in p and "__pycache__" not in p]
    return sorted(set(out))


def main():
    print("L1 an in-class #include of a slice is legal and already used")
    need('#include "Integrator/MMWSpalling/Shared.H"' in M[4739] and '#include "Integrator/MMWSpalling/FaceRemoval.H"' in M[4740],
         "M:4740-4741 are the two existing in-class includes")

    print("L4 no kept test compares an abort's file:line text")
    hits, doc = [], 0
    for p in live_tests():
        in_doc = False
        for ln, s in enumerate(open(p, errors="replace"), 1):
            n3 = s.count('"""')
            was_doc = in_doc
            if n3 % 2:
                in_doc = not in_doc
            if re.search(r"\.H:[0-9]", s) and not s.lstrip().startswith("#"):
                if was_doc or in_doc or n3:
                    doc += 1                     # prose in a docstring compares nothing
                else:
                    hits.append("%s:%d %s" % (os.path.relpath(p, ROOT), ln, s.strip()[:110]))
    for h in hits:
        print("       ", h)
    need(not hits, "`.H:<digits>` appears in no executable test line (%d hit(s); %d more inside docstrings, which compare nothing)" % (len(hits), doc))

    print("L5 source-reading tests (`\"src\"`, src/ or Integrator inside open( / os.path.join( on a code line)")
    rd = collections.defaultdict(list)
    pat = re.compile(r'(open|os\.path\.join|os\.walk)\([^\n]*("src"|src/|"Integrator")')
    for p in live_tests():
        for ln, s in enumerate(open(p, errors="replace"), 1):
            if s.lstrip().startswith("#"):
                continue
            if pat.search(s):
                rd[os.path.relpath(p, ROOT)].append(ln)
    for k, v in sorted(rd.items()):
        print("        %-60s lines %s" % (k, v))
    six = {"tests/MMWSpalling/unit/" + n + "/test" for n in
           ("lateral_flux_lock", "retired_keys", "buckle_lock", "profile_depth", "face_heat", "face_removal")}
    others = sorted(set(rd) - six)
    need(six <= set(rd), "the packet's six are all readers")
    print("        other readers: %s" % (others or "none"))
    need(not others, "no seventh reader of the source tree among the live tests")
    # a reader may also build the path across lines, or point its own variable at it: what it would read is the question
    wide = collections.defaultdict(list)
    for p in live_tests():
        if os.path.relpath(p, ROOT) in six:
            continue
        for ln, s in enumerate(open(p, errors="replace"), 1):
            if not s.lstrip().startswith("#") and re.search(r"MMWSpalling\.H|FaceRemoval\.H|Shared\.H|\"Integrator\"", s):
                wide[os.path.relpath(p, ROOT)].append((ln, s.strip()[:100]))
    for k, v in sorted(wide.items()):
        print("        names a header on a code line: %s %s" % (k, v))

    print("L6 the log writers in the gate reference, by kind")
    cnt = collections.Counter()
    for ln in open(os.path.join(W, "face_p2_28.md5")):
        r = ln.rstrip("\n").split("  ", 1)[1]
        if r.endswith(".csv"):
            cnt[re.sub(r".*run", "", r)] += 1
    print("        reference file count per kind: %s" % dict(cnt))
    need(cnt == {"_removal_events.csv": 19, "_jet_profile.csv": 6, "_side_power.csv": 6, "_buckle_log.csv": 2},
         "19 / 2 / 6 / 6 files as the packet says")
    regimes, sizes = collections.Counter(), collections.Counter()
    gate = os.path.join(W, "face_p2_29")
    for dp, _, fs in os.walk(gate):
        if ".old" in dp:
            continue
        for f in fs:
            if not f.endswith(".csv"):
                continue
            kind = re.sub(r".*run", "", f)
            lines = open(os.path.join(dp, f)).read().splitlines()
            if len(lines) > 1:
                sizes[kind] += 1
            if kind == "_removal_events.csv":
                for row in lines[1:]:
                    regimes[row.split(",")[3]] += 1
    print("        files with at least one data row: %s" % dict(sizes))
    print("        event rows per regime across the gate: %s" % dict(regimes))
    need(all(sizes[k] >= 1 for k in cnt if k != "_buckle_log.csv"), "events, jet profile, side power carry data rows in at least one leg")
    print("        ⚠ FINDING (not in the packet): `_buckle_log.csv` is HEADER-ONLY in both gate legs (BEAM_S, HU_S end before any")
    print("          line fires), so the gate never hashes the buckle ROW writer (F:1011). 4a's bytes claim needs a cover for it:")
    print("          ADDED, in the packet's own vapour pattern: unit/face_buckle's F_S + spall.buckle_log_csv=1 (covers.py B), before / after 4a.")
    smoke = os.path.join(HERE, "covers", "z0smoke", "B", "run_buckle_log.csv")
    nrow = len(open(smoke).read().splitlines()) - 1 if os.path.exists(smoke) else 0
    need(nrow > 0, "the added cover writes buckle-log data rows (%d, smoke on this binary)" % nrow)
    need(regimes["1"] > 0 and regimes["4"] > 0, "firing rows (regime 1) and feet rows (regime 4) are in the gate")
    need(regimes["2"] == 0, "regime 2 absent (the packet expected 0): unit/face_vapor's event log is the vapour writer's only cover")

    print("L9 sp_on and void_aware differ between consumers by design")
    # the packet cites M:1672 for the kernel's sp_on; that line is sp_T's tail. The claim is checked at the definition (M:1664).
    for ln, want in ((1664, "sp_on    = surface_patch.enabled && owner"), (3436, "sp_on    = surface_patch.enabled;"),
                     (1602, "list_state()"), (3423, "list_state()"), (4212, "list_state()"),
                     (3922, "!k_top_col.empty()"), (4129, "!k_top_col.empty()")):
        need(want in M[ln - 1], "M:%d: %s" % (ln, M[ln - 1].strip()[:90]))

    print("git state")
    print("        HEAD %s; %d paths differ from HEAD (the working tree carries the uncommitted packets; pre-images are taken either way)" % (
        subprocess.run(["git", "log", "--oneline", "-1"], cwd=ROOT, capture_output=True, text=True).stdout.strip()[:40],
        len(subprocess.run(["git", "status", "--short"], cwd=ROOT, capture_output=True, text=True).stdout.splitlines())))

    print("the binary and the baselines on disk (list N, face_vapor)")
    import hashlib
    b = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
    sha = hashlib.sha256(open(b, "rb").read()).hexdigest()
    need(sha.startswith("ccd15b91"), "bin/mmwspalling-3d-g++ is the shipped P2-29 binary ccd15b91 (%s)" % sha[:16])
    built = os.path.getmtime(b)
    p29 = os.path.join(T, "studies", "p2_29_one_heat_path")
    for name, path in (("hu_conduction", p29 + "/legs2/hu_conduction"), ("ZQ (Zhang at 0.02, not list N's baseline)", p29 + "/legs2/ZQ"),
                       ("convective_patch", T + "/unit/convective_patch/output"), ("beam", T + "/unit/beam/output"),
                       ("beam_tau", T + "/unit/beam_tau/output"), ("stefan", T + "/unit/stefan/output"),
                       ("latent_heat", T + "/unit/latent_heat/output"), ("face_dirichlet", T + "/unit/face_dirichlet/output"),
                       ("face_vapor", T + "/unit/face_vapor/output")):
        there = os.path.isdir(path)
        age = ""
        if there:
            newest = max((os.path.getmtime(os.path.join(dp, f)) for dp, _, fs in os.walk(path) for f in fs), default=0)
            age = "newest file %+.0f s from the binary's build (%s)" % (newest - built, "after" if newest > built else "BEFORE")
        print("        %-44s %s %s" % (name, "on disk" if there else "ABSENT", age))
    print("preprocessor lines in the class body outside the two existing includes: %s" % (
        [i + 1 for i, s in enumerate(M[64:6821], 64) if s.lstrip().startswith("#") and "Integrator/MMWSpalling/" not in s]))
    print("FAIL (%d)" % len(FAIL) if FAIL else "z0 rows PASS")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
