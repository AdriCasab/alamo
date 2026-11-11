#!/usr/bin/env python3
"""P2-29 -- THE KEPT-TEST RUNNER: every kept test, run from its own deck(s), then scored. Never a check-only score of an
old output (P2-28's b4 scored the `beam` and `convective_patch` outputs of 11-06 / 11-03: its runner never launched
their decks).

For each MANIFEST entry, in order:
  1. archive (MOVE) every output directory the entry writes to <W>/kept_<label>_archive/<path> (SR 9);
  2. launch its deck(s) on bin/ (-np 4) -- `deck` entries; a `self` entry's test launches its own runs;
  3. run the test command; its outcome is PASS (rc 0) or must match the entry's declared RECORD (rc and a signature
     string in the output);
  4. assert every plotfile under the entry's outputs is NEWER than this runner's start (a stale output fails);
then the DISCOVERY CROSS-CHECK: every `unit/*/test` and `validation/**/test*` is either in the manifest or in EXCLUDED
with a reason, or the runner fails. One summary file, <W>/kept_<label>.txt, one line per test.

⚠ CADENCE (the user, 2026-11-10): the FULL run is the user's periodic check, every few days -- not a packet step.
A packet proves the runner on `--subset` (latent_heat, beam, stefan, face_dirichlet) with `--plant-stale`: after the
run, a copy of one fresh plotfile Header is planted under beam's output with an mtime BEFORE the runner's start, and the
staleness assertion must FAIL on it (then the plant is removed).

Usage: run_kept.py run <label> [--subset] [--plant-stale] [name ...]   (the whole manifest, the subset, or names)
       run_kept.py discover                  (the cross-check alone; no run)
Exit 1 on any test outside its declared outcome, a stale output, an unlisted test, or a planted stale output that
passes.
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import time

W = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(W, "..", "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
U, V = "tests/MMWSpalling/unit/", "tests/MMWSpalling/validation/"


def self_test(name, path, args=(), outputs=(), record=None):
    return dict(name=name, kind="self", test=path, args=list(args), decks=[], outputs=list(outputs), record=record)


def deck_test(name, path, decks, args=(), outputs=(), record=None):
    """decks: [(deck, [argv ...], plot_file)]"""
    return dict(name=name, kind="deck", test=path, args=list(args), decks=list(decks), outputs=list(outputs),
                record=record)


def unit_self(n, outputs=True, args=()):
    return self_test(n, U + n + "/test", args, [U + n + "/output"] if outputs else [])


def unit_deck(n, deck="input", out="output", record=None):
    pf = U + n + "/" + out
    return deck_test(n, U + n + "/test", [(U + n + "/" + deck, [], pf)], [pf], [pf], record)


def hu_pair(n, out_g, out_s, deck_g, deck_s):
    base = V + "hu/" + n + "/"
    g, s = base + out_g, base + out_s
    return deck_test(n, base + "test", [(base + deck_g, [], g), (base + deck_s, [], s)], [g, s],
                     [base + os.path.dirname(out_g)])


# P2-29 stage 2: the phase decks took the deck names (the old decks are records under retired_*/).
MANIFEST = [
    # --- source locks and argument checks (they read src/ and run their own short decks)
    unit_self("retired_keys", outputs=False),
    unit_self("lateral_flux_lock", outputs=False),
    unit_self("buckle_lock", outputs=False),
    unit_self("interlocks", outputs=False),
    # --- the face list and its mechanisms (each test launches its own legs)
    unit_self("face_list"), unit_self("face_heat"), unit_self("face_beam"), unit_self("face_buckle"),
    unit_self("face_dirichlet"), unit_self("face_feet"), unit_self("face_melt"), unit_self("face_vapor"),
    dict(unit_self("face_removal"), record=(1, "FAIL (4)")),   # P2-19: the four tilted-line rows (TL20/35/45, TR)
    unit_self("gradient_blind"), unit_self("robin_face", outputs=False), unit_self("profile_depth", outputs=False),
    unit_self("latent_heat"), unit_self("beam_schedule"),
    # P2-29 stage 3: the DR-5 gaps' tests (the march recursion, the beam across a void; face_feet carries the stall
    # and cap legs)
    unit_self("jet_march"), unit_self("beam_void"),
    # --- the fixtures, on the production kernel (P2-29)
    unit_self("beam_tau"),
    unit_deck("beam"),
    unit_deck("equilibrium"),
    unit_deck("stefan"),
    deck_test("convective_patch", U + "convective_patch/test",
              [(U + "convective_patch/input", [], U + "convective_patch/output/granite")],
              [U + "convective_patch/output/granite"], [U + "convective_patch/output"]),
    # --- the microstructure unit decks (score tests)
    unit_deck("alpha_beta_transition"), unit_deck("grain_topology"), unit_deck("heterogeneous_enthalpy"),
    unit_deck("thermal_stress"), unit_deck("voronoi"),
    # --- validation
    deck_test("zhang_oglesby", V + "zhang/zhang_oglesby/test",
              [(V + "zhang/zhang_oglesby/input", [], V + "zhang/zhang_oglesby/output")],
              [V + "zhang/zhang_oglesby/output"], [V + "zhang/zhang_oglesby/output"]),
    hu_pair("hu_conduction", "output/granite2", "output/sandstone2", "input_granite2", "input_sandstone2"),
    hu_pair("hu_thermoelastic", "output/granite2", "output/sandstone2", "input_granite2", "input_sandstone2"),
    deck_test("hu_breakage_index", V + "hu/hu_breakage_index/test",
              [(V + "hu/hu_breakage_index/input_%s" % c, [], V + "hu/hu_breakage_index/output/" + c)
               for c in ("granite2_37", "sandstone2_37", "sandstone2_89")],
              [V + "hu/hu_breakage_index/output"], [V + "hu/hu_breakage_index/output"]),
    self_test("hu_spall_onset_lefm", V + "hu/hu_spall_onset/test_lefm", [], [V + "hu/hu_spall_onset/output_lefm"]),
    self_test("kant_onset_face", V + "kant/sp_kant_onset/test_face", [], [V + "kant/sp_kant_onset/output_face"]),
    self_test("kant_pressure_face", V + "kant/sp_kant_pressure_sweep/test_face", [],
              [V + "kant/sp_kant_pressure_sweep/output_face"]),
    # Oglesby #19 on a SHORT B0 (WHICH_TESTS.md: the face list without a flame); the full 1,800 s is a study leg
    deck_test("oglesby", V + "oglesby/test", [(V + "oglesby/input", ["stop_time=60.0"], V + "oglesby/output_kept/run")],
              [V + "oglesby/output_kept/run"], [V + "oglesby/output_kept"]),
]
SUBSET = ["latent_heat", "beam", "stefan", "face_dirichlet"]
EXCLUDED = {
    "tests/MMWSpalling/retired_units/amr_microstructure_regrid/test": "RETIRED in P2-29 (the user, 2026-11-10: AMR is parked; it fails on a level-0 "
                                          "is_grain_boundary mismatch since before P2-27); stage 2 moves it to "
                                          "retired_units/, re-instated with the AMR-repair packet",
    "tests/MMWSpalling/retired_validation/separate_path/hu_thermoelastic_amr/test": "RETIRED in P2-29 (the user's trim, 2026-11-10: AMR is parked; hu_thermoelastic "
                                        "covers the physics); stage 2 moves it to retired_validation/separate_path/",
    V + "meier/sp_meier_pilot/test": "a RECORD since P2-28: the Meier column decks ran the deleted column path",
    V + "meier/sp_meier_pilot/test_feet": "a RECORD since P2-28 (the column path's feet)",
    V + "meier/sp_meier_pilot/test_feet_d2c": "a RECORD since P2-28 (the column path's D2c feet)",
}


def discover():
    found = sorted(set(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(T, "unit", "*", "test"))) |
                   set(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(T, "validation", "**", "test*"),
                                                                   recursive=True)
                       if os.path.isfile(p) and not p.endswith(".pyc") and "/output" not in p and "/legs/" not in p
                       and "__pycache__" not in p))
    listed = {e["test"] for e in MANIFEST}
    missing = [p for p in found if p not in listed and p not in EXCLUDED]
    ghost = [p for p in listed if not os.path.isfile(os.path.join(ROOT, p))]
    print("discovery: %d tests on disk; %d in the manifest, %d excluded with a reason" % (
        len(found), len(listed), len(EXCLUDED)))
    for p in missing:
        print("  UNLISTED  %s (neither in the manifest nor in EXCLUDED)" % p)
    for p in ghost:
        print("  MISSING   %s (in the manifest, not on disk)" % p)
    return missing, ghost


def plotfiles_under(paths):
    out = []
    for p in paths:
        out += glob.glob(os.path.join(ROOT, p, "**", "Header"), recursive=True)
    return out


def stale_of(e, start):
    return [p for p in plotfiles_under(e["outputs"]) if os.path.getmtime(p) < start]


def run_entry(e, label, start, log):
    arch = os.path.join(W, "kept_%s_archive" % label)
    for o in e["outputs"]:
        src = os.path.join(ROOT, o)
        if os.path.exists(src):
            dst = os.path.join(arch, o)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if os.path.exists(dst):
                shutil.rmtree(dst)
            shutil.move(src, dst)
    txt = ""
    for deck, extra, pf in e["decks"]:
        os.makedirs(os.path.dirname(os.path.join(ROOT, pf)), exist_ok=True)
        argv = MPI + [BIN, deck] + extra + ["plot_file=" + pf]
        r = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
        txt += "$ " + " ".join(argv) + "\n" + (r.stdout or "")[-2000:] + (r.stderr or "")[-2000:] + "\n"
        if r.returncode != 0:
            return "DECK rc %d (%s)" % (r.returncode, deck), txt
    args = [os.path.join(ROOT, a) if a.startswith("tests/") else a for a in e["args"]]
    r = subprocess.run([PY, os.path.join(ROOT, e["test"])] + args, cwd=ROOT, capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    txt += "$ test %s %s\n%s" % (e["test"], " ".join(args), out)
    stale = stale_of(e, start)
    fresh = [p for p in plotfiles_under(e["outputs"]) if os.path.getmtime(p) >= start]
    if stale:
        return "STALE: %d plotfile(s) older than the runner (%s)" % (len(stale), os.path.relpath(stale[0], ROOT)), txt
    if e["kind"] == "deck" and not fresh:
        return "NO PLOTFILE written under %s" % e["outputs"], txt
    last = [ln for ln in out.strip().splitlines() if ln.strip()][-1:] or [""]
    if e["record"]:
        rc_want, sig = e["record"]
        ok = r.returncode == rc_want and sig in out
        return ("RECORD %s (rc %d, '%s')" % ("as declared" if ok else "NOT AS DECLARED", r.returncode, sig),
                txt) if ok else ("UNEXPECTED: rc %d, signature '%s' %s" % (r.returncode, sig,
                                                                          "present" if sig in out else "absent"), txt)
    return ("PASS (%d fresh plotfile files)" % len(fresh) if r.returncode == 0 else
            "FAIL rc %d: %s" % (r.returncode, last[0][:120])), txt


def main(argv):
    if not argv or argv[0] not in ("run", "discover"):
        print(__doc__)
        return 1
    missing, ghost = discover()
    if argv[0] == "discover":
        return 1 if (missing or ghost) else 0
    label = argv[1]
    names = [a for a in argv[2:] if not a.startswith("--")]
    if "--subset" in argv:
        names = SUBSET
    entries = [e for e in MANIFEST if not names or e["name"] in names]
    start = time.time()
    rows, bad = [], 0
    logd = os.path.join(W, "kept_%s_logs" % label)
    os.makedirs(logd, exist_ok=True)
    for e in entries:
        t0 = time.time()
        res, txt = run_entry(e, label, start, None)
        open(os.path.join(logd, e["name"] + ".log"), "w").write(txt)
        ok = res.startswith("PASS") or res.startswith("RECORD as declared")
        bad += not ok
        rows.append("%-26s %-6s %6.0f s  %s" % (e["name"], "ok" if ok else "BAD", time.time() - t0, res))
        print(rows[-1], flush=True)
    if "--plant-stale" in argv:
        # the mtime assertion's planted control: one stale Header under the first deck entry's output
        e = next(x for x in entries if x["kind"] == "deck" and plotfiles_under(x["outputs"]))
        src = plotfiles_under(e["outputs"])[0]
        pdir = os.path.join(ROOT, e["outputs"][0], "99999cell")
        os.makedirs(pdir, exist_ok=True)
        dst = os.path.join(pdir, "Header")
        shutil.copy(src, dst)
        os.utime(dst, (start - 3600.0, start - 3600.0))
        caught = any(os.path.samefile(p, dst) for p in stale_of(e, start))
        shutil.rmtree(pdir)
        bad += not caught
        rows.append("%-26s %-6s %6s    %s" % ("(planted stale output)", "ok" if caught else "BAD", "-",
                                            "the mtime assertion %s a Header aged 1 h under %s" % (
                                                "FAILS as it must on" if caught else "MISSED", e["outputs"][0])))
        print(rows[-1], flush=True)
    sha = subprocess.run(["shasum", "-a", "256", BIN], capture_output=True, text=True).stdout.split()[0]
    head = ["run_kept.py %s: %d tests on bin %s; started %s; %d outside their declared outcome; discovery: %d unlisted, "
            "%d missing" % (label, len(entries), sha[:16], time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start)),
                            bad, len(missing), len(ghost))]
    open(os.path.join(W, "kept_%s.txt" % label), "w").write("\n".join(head + rows) + "\n")
    print(head[0])
    return 1 if (bad or missing or ghost) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
