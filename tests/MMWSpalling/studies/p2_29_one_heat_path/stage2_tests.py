#!/usr/bin/env python3
"""P2-29 stage 2, Goal 10 -- THE PHASE DECKS BECOME THE DECKS; the old ones and the AMR tests become records.
Moves (never overwrites; every move logged to S/STAGE2_TESTS.txt):
  * the old fixture decks + their outputs/logs -> retired_units/<name>_separate_path/ (beam, beam_tau, equilibrium,
    stefan, convective_patch); the old Zhang, Hu conduction / thermoelastic / breakage decks + outputs/logs ->
    retired_validation/separate_path/<test>/;
  * hu_thermoelastic_amr (the whole test) -> retired_validation/separate_path/hu_thermoelastic_amr/ and
    unit/amr_microstructure_regrid (the whole test) -> retired_units/amr_microstructure_regrid/ (the user: AMR parked);
  * the phase decks take the deck names (`input_phase` -> `input`, `input_<x>_phase` -> `input_<x>`) and write to
    output/ (plot_file `/output_phase` -> `/output`); the stage-1 `unit/beam_tau/output_phase` -> S/archive/.
Edits: unit/beam_tau/test (U3 retired, --deck default `input`), unit/interlocks/test (D4/D5 retired with AMR: the
user, 2026-11-10), witness/run_kept.py (the manifest's deck names), the two READMEs.
Usage: stage2_tests.py"""
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
LOG = []


def log(msg):
    LOG.append(msg)
    print(msg)


def mv(src, dst):
    s, d = os.path.join(T, src), os.path.join(T, dst)
    if not os.path.exists(s):
        log("  (absent) %s" % src)
        return
    assert not os.path.exists(d), "refusing to overwrite " + dst
    os.makedirs(os.path.dirname(d), exist_ok=True)
    shutil.move(s, d)
    log("  moved %s -> %s" % (src, dst))


def edit(rel, pairs):
    p = os.path.join(T, rel)
    s = open(p).read()
    for a, b in pairs:
        assert s.count(a) == 1, "%s: anchor %r" % (rel, a[:70])
        s = s.replace(a, b)
    open(p, "w").write(s)
    log("  edited %s (%d change(s))" % (rel, len(pairs)))


log("== the old decks and their outputs -> records")
for n in ("beam", "beam_tau", "equilibrium", "stefan", "convective_patch"):
    base = "unit/%s/" % n
    for x in ("input", "output", "sweep_run.log", "sweep_sim.log"):
        mv(base + x, "retired_units/%s_separate_path/%s" % (n, x))
mv("unit/beam_tau/output_phase", "studies/p2_29_one_heat_path/archive/unit_beam_tau_output_phase_stage1")
mv("validation/zhang/zhang_oglesby/input", "retired_validation/separate_path/zhang_oglesby/input")
mv("validation/zhang/zhang_oglesby/output", "retired_validation/separate_path/zhang_oglesby/output")
for t, files in (("hu_conduction", ["input_granite2", "input_sandstone2", "output", "granite.log", "sandstone.log",
                                    "sweep_run.log"]),
                 ("hu_thermoelastic", ["input_granite2", "input_sandstone2", "output", "granite.log",
                                       "granite_topfix.log", "sandstone.log", "sandstone_topfix.log", "sweep_run.log"]),
                 ("hu_breakage_index", ["input_granite2_37", "input_sandstone2_37", "input_sandstone2_89",
                                        "granite37_roller.log", "granite37_run.log", "granite37_topfix.log",
                                        "sandstone37_roller.log", "sandstone37_run.log", "sandstone37_topfix.log",
                                        "sandstone89_roller.log", "sandstone89_run.log", "sandstone89_topfix.log",
                                        "sweep_run.log"])):
    for x in files:
        mv("validation/hu/%s/%s" % (t, x), "retired_validation/separate_path/%s/%s" % (t, x))
log("== the AMR tests -> records (the user: AMR parked)")
mv("validation/hu/hu_thermoelastic_amr", "retired_validation/separate_path/hu_thermoelastic_amr")
mv("unit/amr_microstructure_regrid", "retired_units/amr_microstructure_regrid")

log("== the phase decks become the decks")
REN = [("unit/%s/input_phase" % n, "unit/%s/input" % n) for n in
       ("beam", "beam_tau", "equilibrium", "stefan", "convective_patch")] + \
      [("validation/zhang/zhang_oglesby/input_phase", "validation/zhang/zhang_oglesby/input")] + \
      [("validation/hu/hu_conduction/input_%s_phase" % r, "validation/hu/hu_conduction/input_%s" % r)
       for r in ("granite2", "sandstone2")] + \
      [("validation/hu/hu_thermoelastic/input_%s_phase" % r, "validation/hu/hu_thermoelastic/input_%s" % r)
       for r in ("granite2", "sandstone2")] + \
      [("validation/hu/hu_breakage_index/input_%s_phase" % r, "validation/hu/hu_breakage_index/input_%s" % r)
       for r in ("granite2_37", "sandstone2_37", "sandstone2_89")]
for a, b in REN:
    mv(a, b)
    p = os.path.join(T, b)
    s = open(p).read()
    n_pf = len(re.findall(r"/output_phase\b", s))
    s = re.sub(r"/output_phase\b", "/output", s)
    s = s.replace("# retires in P2-29 stage 2. Changes:", "# retired in P2-29 stage 2 (the old deck is a record under "
                  "retired_*/; this deck took its name). Changes:")
    open(p, "w").write(s)
    log("    %s: plot_file output_phase -> output (%d)" % (b, n_pf))

log("== tests")
edit("unit/beam_tau/test", [
    ('''  U3  unit/beam's OWN deck (no table) against studies/p2_20_beam_deposition/unit_beam_output_pre/: every plotfile
      (Header, Cell_H, Cell_D_*) byte-identical -- the closed form untouched
''', '''  (U3, a byte lock on unit/beam's const-alpha output, retired in P2-29 with the const-alpha path:
      retired_units/beam_tau_separate_path/README line)
'''),
    ('''  P2-29: --deck input_phase runs U1 and U2 on the production-kernel deck (one constant phase, kappa 0, rho = c_p = 1:
  H [J/m^3] = the old H [K]) into output_phase/; U3 is a byte lock on the const-alpha path's output and retires with
  it (it runs only on the old deck). --bin PATH runs every leg on that binary (P2-29 m2: the old deck on bin_pre).''',
     '''  P2-29: the deck is on the production kernel (one constant phase, kappa 0, rho = c_p = 1: H [J/m^3] = the old
  H [K]). --deck NAME runs another deck of this directory into output_<name>/; --bin PATH runs every leg on that
  binary.'''),
    ('''BEAM_DECK = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "beam", "input")
PRE = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_20_beam_deposition", "unit_beam_output_pre")
''', ''),
])
p = os.path.join(T, "unit/beam_tau/test")
s = open(p).read()
a = s.index("    # ---- U3 (the const-alpha byte lock")
b = s.index('    print("=" * 100)\n    print("FAIL (%d)" % len(FAIL) if FAIL else "PASS")\n    return 1 if FAIL else 0\n',
            a + 1)
s = s[:a] + s[b:]
s = s.replace("import filecmp\n", "")
open(p, "w").write(s)
log("  edited unit/beam_tau/test: U3 removed")
edit("unit/interlocks/test", [
    ('''    rg = os.path.join(U, "amr_microstructure_regrid", "input")
    for tag, extra, nlev, warn_lev in (("D4", [], 2, []), ("D5", ["timestep=2.0e-4"], 2, [1])):
        rc, txt = run(BIN, rg, extra, np_=4)
        msgs = sorted(set(int(m) for m in re.findall(r"explicit time-step limit, level (\\d+):", txt)))
        warns = sorted(set(int(m.group(1)) for m in re.finditer(r"level (\\d+): timestep dt = \\S+ s exceeds the "
                                                                 r"explicit limit dt_lim", txt)))
        need(rc == 0 and msgs == list(range(nlev)) and warns == warn_lev,
             "%s amr_microstructure_regrid %s runs (rc %d); per-level messages %s; warnings on levels %s (want %s)"
             % (tag, " ".join(extra) or "(as written)", rc, msgs, warns, warn_lev))
''', '''    # P2-29 (the user, 2026-11-10): D4 / D5 -- the per-level limit on the 2-level amr_microstructure_regrid deck --
    # are RETIRED with AMR (parked; that deck moved to retired_units/). CheckTimeStepLimit's level >= 1 branch is
    # untested until the AMR-repair packet re-instates them.
'''),
])
p = os.path.join(T, "unit/interlocks/test")
s = open(p).read()
m = re.search(r"\n  D4 [^\n]*\n(  D5 [^\n]*\n)?", s)
if m:
    s = s[:m.start()] + "\n  (D4 / D5, the per-level limit on an AMR deck, retired with AMR in P2-29: the user, 2026-11-10)\n" + \
        s[m.end():]
    open(p, "w").write(s)
    log("  edited unit/interlocks/test: D4/D5 docstring lines replaced")
edit("studies/tree_sweep_0922/witness/run_kept.py", [
    ('''# ⚠ The P2-29 stage-1 deck names (`input_phase`, `*_phase`): stage 2 renames the phase decks to the deck names and
# edits these rows with them.
''', '''# P2-29 stage 2: the phase decks took the deck names (the old decks are records under retired_*/).
'''),
    ('''    self_test("beam_tau", U + "beam_tau/test", ["--deck", "input_phase"], [U + "beam_tau/output_phase"]),
    unit_deck("beam", "input_phase", "output_phase"),
    unit_deck("equilibrium", "input_phase", "output_phase"),
    unit_deck("stefan", "input_phase", "output_phase"),
    deck_test("convective_patch", U + "convective_patch/test",
              [(U + "convective_patch/input_phase", [], U + "convective_patch/output_phase/granite")],
              [U + "convective_patch/output_phase/granite"], [U + "convective_patch/output_phase"]),''',
     '''    unit_self("beam_tau"),
    unit_deck("beam"),
    unit_deck("equilibrium"),
    unit_deck("stefan"),
    deck_test("convective_patch", U + "convective_patch/test",
              [(U + "convective_patch/input", [], U + "convective_patch/output/granite")],
              [U + "convective_patch/output/granite"], [U + "convective_patch/output"]),'''),
    ('''              [(V + "zhang/zhang_oglesby/input_phase", [], V + "zhang/zhang_oglesby/output_phase")],
              [V + "zhang/zhang_oglesby/output_phase"], [V + "zhang/zhang_oglesby/output_phase"]),
    hu_pair("hu_conduction", "output_phase/granite2", "output_phase/sandstone2",
            "input_granite2_phase", "input_sandstone2_phase"),
    hu_pair("hu_thermoelastic", "output_phase/granite2", "output_phase/sandstone2",
            "input_granite2_phase", "input_sandstone2_phase"),''',
     '''              [(V + "zhang/zhang_oglesby/input", [], V + "zhang/zhang_oglesby/output")],
              [V + "zhang/zhang_oglesby/output"], [V + "zhang/zhang_oglesby/output"]),
    hu_pair("hu_conduction", "output/granite2", "output/sandstone2", "input_granite2", "input_sandstone2"),
    hu_pair("hu_thermoelastic", "output/granite2", "output/sandstone2", "input_granite2", "input_sandstone2"),'''),
    ('''              [(V + "hu/hu_breakage_index/input_%s_phase" % c, [], V + "hu/hu_breakage_index/output_phase/" + c)
               for c in ("granite2_37", "sandstone2_37", "sandstone2_89")],
              [V + "hu/hu_breakage_index/output_phase"], [V + "hu/hu_breakage_index/output_phase"]),''',
     '''              [(V + "hu/hu_breakage_index/input_%s" % c, [], V + "hu/hu_breakage_index/output/" + c)
               for c in ("granite2_37", "sandstone2_37", "sandstone2_89")],
              [V + "hu/hu_breakage_index/output"], [V + "hu/hu_breakage_index/output"]),'''),
    ('''    V + "hu/hu_thermoelastic_amr/test": "RETIRED in P2-29''', '''    "tests/MMWSpalling/retired_validation/separate_path/hu_thermoelastic_amr/test": "RETIRED in P2-29'''),
    ('''    U + "amr_microstructure_regrid/test": "RETIRED in P2-29''', '''    "tests/MMWSpalling/retired_units/amr_microstructure_regrid/test": "RETIRED in P2-29'''),
])

log("== READMEs")
R = os.path.join(T, "retired_units", "README.md")
s = open(R).read()
assert s.count("- `jet_d2c` -- P2-28: the D2c steady closures on the column path (feed cap, stall, foot_rule modes)") == 1
assert s.count("- `robin_feet` -- P2-28: the burner on its feet on the column path [face_feet]") == 1
s = s.replace("- `jet_d2c` -- P2-28: the D2c steady closures on the column path (feed cap, stall, foot_rule modes) [face_feet; "
              "the feet descent's feed cap / stall / foot_rule modes are a named coverage gap]",
              "- `jet_d2c` -- P2-28: the D2c steady closures on the column path [face_feet] (P2-29: the feed cap / stall / "
              "foot_rule specs were `robin_feet`'s (FM/FA/FC/SM/SP), not this test's -- see that line)")
s = s.replace("- `robin_feet` -- P2-28: the burner on its feet on the column path [face_feet]",
              "- `robin_feet` -- P2-28: the burner on its feet on the column path [face_feet]; its FM/FA/FC/SM/SP legs are "
              "the spec of the feet descent's feed cap / stall / foot_rule modes (P2-29 stage 3: face_feet's stall and cap "
              "legs; foot_rule max / mean deleted)")
s = s.rstrip("\n") + '''

P2-29 (census D): THE SEPARATE HEAT KERNELS are deleted (one heat path, AdvanceMicrostructure). Each test below kept
its name and took a production-kernel deck (one constant phase, physical units); the OLD deck and its outputs are
the record here (`<name>_separate_path/`). Not re-run.
- `beam_separate_path`, `equilibrium_separate_path`, `stefan_separate_path`, `convective_patch_separate_path`,
  `beam_tau_separate_path` -- the const-alpha / material decks (H in kelvin, rho c_p = 1) and their last outputs
  [the same test on `unit/<name>/input`]; `beam_tau`'s U3 (a byte lock on unit/beam's const-alpha output) retired
  with them
- `amr_microstructure_regrid` -- RETIRED (the user, 2026-11-10): AMR is parked (SR 12); red since before P2-27 (a
  level-0 is_grain_boundary mismatch after regrid); AMR + mechanics also fails on the production path (MLMG at the
  first regrid). Re-instate with the AMR-repair packet. `unit/interlocks` D4/D5 (the per-level limit, on this deck)
  retired with it.
'''
open(R, "w").write(s)
log("  retired_units/README.md: the P2-29 lines + the robin_feet / jet_d2c attribution corrected")
RV = os.path.join(T, "retired_validation", "separate_path", "README.md")
open(RV, "w").write('''# The separate heat kernels' validation decks (retired in P2-29)

P2-29 (census D) deleted `AdvanceConstAlpha`, `AdvanceMaterial` and `AdvanceMaterialImplicit`: one heat path,
`AdvanceMicrostructure`. Each test below kept its name and took a production-kernel deck (one constant phase, or
`microstructure.phase0.law = zhang`); the OLD decks, their logs and their last outputs are the record here. A deck here
would now ABORT (`material.*` / `phase.*` / `implicit.*` retired by presence, `microstructure.nphases` required).

- `zhang_oglesby/` -- the global Zhang material + the implicit solver (2 s) [validation/zhang/zhang_oglesby/input:
  `phase0.law = zhang` on a QUARTER block, explicit 0.05 s; P2-29 a3/a4: the production kernel = the old explicit
  kernel = the full block to 2.3e-13 K]
- `hu_conduction/`, `hu_thermoelastic/` -- the material kernel with the HALF-conductance prescribed term and the
  homogeneous mechanics; `output/` is the 09-25 record (P2-29 m4/m5's baseline) [the same tests on the phase decks]
- `hu_breakage_index/` -- the AMR decks (max_level 1, no output since 09-17) [the same test on single-level 40^3 phase
  decks: AMR + mechanics fails on the production path, the user 2026-11-10]
- `hu_thermoelastic_amr/` -- RETIRED whole (decks + test; the user, 2026-11-10): AMR is parked (SR 12);
  `hu_thermoelastic` covers the thermoelastic physics. Re-instate with the AMR-repair packet.
''')
log("  retired_validation/separate_path/README.md written")
open(os.path.join(HERE, "STAGE2_TESTS.txt"), "w").write("\n".join(LOG) + "\n")
