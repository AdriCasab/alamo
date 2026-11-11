#!/usr/bin/env python3
"""P2-28 stage 2 -- the test side of deleting the column path (items 8, 9): moves, README lines, record headers, deck
edits. Idempotent: a moved directory is skipped, a header already present is not added twice, a deleted line stays
deleted. Writes S/STAGE2_TESTS.txt (every move, header and deleted line).

  MOVE -> retired_units/        the unit tests of the column path (each needs it or its keys)
  MOVE -> retired_validation/column_path/   site A's onset decks and their tests, the meiercfg onset, p0_steady_column,
                                rossi -- an onset deck would RUN silently as a thermal deck once Sp_field is gone (SR 1)
  RECORD (in place, a header)   the Meier column decks and tests, the Hu LEFM extension decks, the s1 study input: they
                                abort loudly now
  DECKS                         live decks lose the lines retired at their default (spall.surface_T = reconstructed,
                                spall.removal_energy = flake, surface.follow_mask = 1, weibull.detach_mode = per_face,
                                surface_patch.pinned_idle_cycles = any: inert) and the declared LIVE lines
                                (unit/face_feet/input spall.h_col_events_csv; input_feet_face bit.*)
"""
import glob
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
T = os.path.abspath(os.path.join(HERE, "..", ".."))
LOG = []

UNITS = {  # name -> (why, successor)
    "robin_pinned": ("the pinned Robin closure (robin_form = pinned) and its pin rule", "face_heat, face_list"),
    "robin_jet": ("the jet map on the column path's Robin forms", "face_heat (v)/(vi), the face gate"),
    "jet_enthalpy": ("the jet enthalpy march's column flux (per-bin recursion vs the strip closed form)",
                     "face_heat (v)/(vi) + the face gate; the per-bin recursion is a named coverage gap"),
    "robin_feet": ("the burner on its feet on the column path", "face_feet"),
    "jet_d2c": ("the D2c steady closures on the column path (feed cap, stall, foot_rule modes)",
                "face_feet; the feet descent's feed cap / stall / foot_rule modes are a named coverage gap"),
    "pad_clear_ledger": ("the column path's pad-clearance ledger", "face_feet (e)"),
    "pin_clock": ("the pinned closure's idle clock", "none (robin_form = pinned is retired)"),
    "removal_increment": ("the column removal's whole-cell increment and its prescribed-T leg", "face_removal"),
    "scalar_flaw": ("the deterministic flaw on the column criterion", "face_removal, face_buckle"),
    "beam_void_closure": ("the beam path through the column path's voided columns",
                          "beam, beam_tau, face_beam; a beam path across an interior void is a named coverage gap"),
    "seed_void": ("seeded interior voids for the column scans", "face_list (CAV)"),
    "surface_T_recon": ("the column reconstruction (spall.surface_T = cell_average | diagnostic | reconstructed)",
                        "face_removal, face_heat (the face list's R1)"),
    "sp_v_n_regime": ("the column path's spall / vapour regime field", "face_vapor"),
    "sp_lefm": ("site A's LEFM checks (prescribed stress and T hooks)",
                "validation/kant/sp_kant_onset/test_face, sp_kant_pressure_sweep/test_face"),
}
VALID = {  # path under validation/ -> destination name under retired_validation/column_path/
    "kant/sp_kant_onset_meiercfg": "sp_kant_onset_meiercfg",
    "p0_steady_column": "p0_steady_column",
    "rossi": "rossi",
}
VALID_FILES = {  # the site-A files of the two kept Kant dirs (their face successors stay)
    "kant/sp_kant_onset": ["input_verification", "input_verification_hfl", "input_model", "input_weibull",
                           "input_weibull_cluster", "test", "output", "output.p1g0_archive", "sweep_run.log"],
    "kant/sp_kant_pressure_sweep": ["input", "test", "output", "output.p1g0_archive", "sweep_run.log"],
}
RECORDS = ["validation/meier/sp_meier_pilot/" + n for n in
           ("input_2d_dev", "input_2d_dev_melt", "input_energy", "input_drilling", "input_drilling_25mm",
            "input_feet", "input_feet_d2c", "test", "test_feet", "test_feet_d2c")] + \
          ["studies/s1_surface_resolution/input"] + \
          sorted(os.path.relpath(p, T) for p in glob.glob(os.path.join(T, "extensions", "hu_spall_onset_mmwbeam_lefm",
                                                                       "input*")))
HEADER = ("# P2-28 RECORD: this file ran the COLUMN removal path, retired in P2-28 (one removal path, the face list); "
          "it now aborts at parse. Kept as the record of its last runs, not run.")
DEFAULT = {"spall.surface_T": "reconstructed", "spall.removal_energy": "flake", "surface.follow_mask": "1",
           "weibull.detach_mode": "per_face"}
LINE = re.compile(r"^\s*([A-Za-z0-9_.]+)\s*=\s*(.*?)\s*(#.*)?$")


def mv(src, dst):
    if not os.path.exists(src):
        LOG.append("  (absent, skipped) %s" % os.path.relpath(src, T))
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    assert not os.path.exists(dst), dst
    shutil.move(src, dst)
    LOG.append("  moved %s -> %s" % (os.path.relpath(src, T), os.path.relpath(dst, T)))


def main():
    LOG.append("== moves")
    RU = os.path.join(T, "retired_units")
    for n in UNITS:
        mv(os.path.join(T, "unit", n), os.path.join(RU, n))
    os.makedirs(os.path.join(RU, "robin_face_ledger"), exist_ok=True)
    mv(os.path.join(T, "unit", "robin_face", "input_ledger"), os.path.join(RU, "robin_face_ledger", "input_ledger"))
    gb = os.path.join(T, "unit", "gradient_blind")
    if os.path.exists(os.path.join(gb, "input_face")):
        mv(os.path.join(gb, "input"), os.path.join(RU, "gradient_blind_column", "input"))
        shutil.move(os.path.join(gb, "input_face"), os.path.join(gb, "input"))
        LOG.append("  renamed unit/gradient_blind/input_face -> input")
    RV = os.path.join(T, "retired_validation", "column_path")
    for src, dst in VALID.items():
        mv(os.path.join(T, "validation", src), os.path.join(RV, dst))
    for d, files in VALID_FILES.items():
        for f in files:
            mv(os.path.join(T, "validation", d, f), os.path.join(RV, os.path.basename(d), f))

    LOG.append("== READMEs")
    rp = os.path.join(RU, "README.md")
    txt = open(rp).read()
    if "P2-28 (census H)" not in txt:
        txt = txt.rstrip("\n") + "\n\nP2-28 (census H): THE COLUMN REMOVAL PATH is deleted (one removal path, the face list). " \
              "Each test below needed it or one of its keys (which now abort); not re-run. Successor in brackets.\n"
        for n, (why, succ) in UNITS.items():
            txt += "- `%s` -- P2-28: %s [%s]\n" % (n, why, succ)
        txt += "- `robin_face_ledger` -- P2-28: robin_face's (b) legs, the ledger on the column path " \
               "(`input_ledger`) [face_feet (e), face_vapor, face_heat (v)]; robin_face's (a) form legs stay\n"
        txt += "- `gradient_blind_column` -- P2-28: the uniform slab on the column path (removal_energy = cell) " \
               "[unit/gradient_blind, now the face-list deck]\n"
        open(rp, "w").write(txt)
        LOG.append("  retired_units/README.md: %d lines added" % (len(UNITS) + 3))
    rv = os.path.join(RV, "README.md")
    if not os.path.exists(rv):
        open(rv, "w").write(
            "# The column path's validation decks (retired in P2-28)\n\n"
            "Site A -- the column path's no-removal onset criterion (`Sp_field`) -- is deleted in P2-28 (census H, decision 8 "
            "of 2026-11-05). These decks and tests ran it; an onset deck would now RUN silently as a thermal-only deck (its "
            "criterion field is gone, standing rule 1), so they are moved here as records, not run. Their last runs are "
            "P2-28 stage 1's (`studies/p2_28_column_path/K1_kant_siteA.log`, `K1_pressure_siteA.log`), on the stage-1 binary.\n\n"
            "- `sp_kant_onset/` -- the five site-A onset decks (`input_verification`, `_hfl`, `input_model`, `input_weibull`, "
            "`_cluster`) and `test` (targets 1-6 and the 25-case sweep) [validation/kant/sp_kant_onset/test_face: targets 3, "
            "5 on the face list's evaluate-only mode]\n"
            "- `sp_kant_pressure_sweep/` -- the site-A base deck and `test` (R1-R4) [validation/kant/sp_kant_pressure_sweep/"
            "test_face]\n"
            "- `sp_kant_onset_meiercfg/` -- the onset at Meier's configuration on the column path [validation/kant/"
            "grimsel_calibration (c6)]\n"
            "- `p0_steady_column/` -- the steady-front reconstruction check on the column path [unit/face_removal --floor]\n"
            "- `rossi/` -- the Rossi damage-profile deck (decision 3 of 2026-11-05: retire)\n"
            "The decks write to their old `tests/MMWSpalling/validation/...` paths: records, not runnable as moved.\n")
        LOG.append("  retired_validation/column_path/README.md written")

    LOG.append("== record headers")
    for rel in RECORDS:
        p = os.path.join(T, rel)
        if not os.path.exists(p):
            LOG.append("  (absent) %s" % rel)
            continue
        txt = open(p).read()
        if HEADER in txt:
            continue
        if rel.endswith(("test", "test_feet", "test_feet_d2c")) and txt.startswith("#!"):
            first, rest = txt.split("\n", 1)
            txt = first + "\n" + HEADER + "\n" + rest
        else:
            txt = HEADER + "\n" + txt
        open(p, "w").write(txt)
        LOG.append("  header: %s" % rel)

    LOG.append("== deck edits (live decks: unit/, validation/, extensions/ minus the records; shared includes)")
    decks = []
    for base in ("unit", "validation"):
        for p in glob.glob(os.path.join(T, base, "**", "*"), recursive=True):
            if not os.path.isfile(p) or "/output" in p or ".old" in p or "/legs/" in p:
                continue
            b = os.path.basename(p)
            if not (b.startswith("input") or b.endswith(".inputs")):
                continue
            if os.path.relpath(p, T) in RECORDS:
                continue
            decks.append(p)
    n_del = 0
    for p in sorted(decks):
        rel = os.path.relpath(p, T)
        L = open(p).read().split("\n")
        keep = []
        for i, ln in enumerate(L):
            m = LINE.match(ln)
            drop = None
            if m:
                k, v = m.group(1), m.group(2).strip()
                if k in DEFAULT and v == DEFAULT[k]:
                    drop = "DEFAULT"
                elif k == "surface_patch.pinned_idle_cycles":
                    drop = "INERT"
                elif k == "spall.h_col_events_csv" and rel == "unit/face_feet/input":
                    drop = "LIVE (declared: FEET_S's h_col CSV leaves)"
                elif k.startswith("bit.") and rel == "validation/meier/sp_meier_pilot/input_feet_face":
                    drop = "LIVE (declared: the bit was off on the face list)"
            if drop:
                LOG.append("  %-8s %s:%d  %s" % (drop.split()[0], rel, i + 1, ln.strip()))
                n_del += 1
                continue
            keep.append(ln)
        if len(keep) != len(L):
            open(p, "w").write("\n".join(keep))
    LOG.append("  (%d lines deleted)" % n_del)
    out = "\n".join(LOG) + "\n"
    open(os.path.join(HERE, "STAGE2_TESTS.txt"), "a").write(out)
    print(out)


main()
