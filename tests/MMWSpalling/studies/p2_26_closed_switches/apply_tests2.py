#!/usr/bin/env python3
"""P2-26 stage 2 -- the test edits that go with C2 + C3 (exact pairs, each required to match ONCE), then the move of
unit/sp_melt_skeleton to retired_units/ with its README line.

  retired_keys      R8 (the four C2/C3 keys, once each, in the table) + its mutation; K10 (the four at their old
                    defaults run), K11 surface_normal=1, K12 flake_coherence_length=0.004, K13 melt_aware=1 (abort)
  lateral_flux_lock the depth scan's stress regex loses `f_skel \\* \\(` (ScanSigmaT's scalar return is now
                    `return sp_relief * phase_E_loc ...`); the P2-26 name clause gains the C2 and C3 identifiers
  face_heat         the BASE tokens spall.flake_coherence_length=0.0 / spall.surface_normal=0 are dropped (the table
                    compares strings: 0.0 != 0, so the argv must not write them)
  removal_increment the SN leg and checks (e)-(e4), (f) retire with spall.surface_normal
  sp_melt_skeleton  moves to tests/MMWSpalling/retired_units/ (it A/B-tests the deleted skeleton)
Usage: apply_tests2.py --dry | --apply
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
U = os.path.join(ROOT, "tests", "MMWSpalling", "unit")
RK = os.path.join(U, "retired_keys", "test")
LF = os.path.join(U, "lateral_flux_lock", "test")
FH = os.path.join(U, "face_heat", "test")
RI = os.path.join(U, "removal_increment", "test")
RU = os.path.join(ROOT, "tests", "MMWSpalling", "retired_units")
EDITS = []


def E(path, what, old, new):
    EDITS.append((path, what, old, new))


# ---------------------------------------------------------------- retired_keys ----
E(RK, "docstring R8", '''     exactly ONCE (comments included), as an entry of the retired-at-default table with that value
Each clause has a MUTATION''', '''     exactly ONCE (comments included), as an entry of the retired-at-default table with that value
  R8 (P2-26, census C2 + C3) the same for spall.surface_normal (0), spall.flake_coherence_length (0),
     spall.melt_aware (0) and melt.lambda_crit (*): the "0" is a STRING, so a deck writing 0.0 aborts
Each clause has a MUTATION''')
E(RK, "docstring K10-K13", '''  K9 + surface_patch.pad_shield_open=0.5 -> rc != 0, the strict surface_patch scope names the key (unparsed)
''', '''  K9 + surface_patch.pad_shield_open=0.5 -> rc != 0, the strict surface_patch scope names the key (unparsed)
  P2-26 (census C2 + C3):
  K10 + the four C2/C3 keys at their old defaults (0, 0, 0, 0.4) -> rc == 0
  K11 + spall.surface_normal=1           -> rc != 0, output names the key and P2-26
  K12 + spall.flake_coherence_length=0.004 -> rc != 0, output names the key and P2-26
  K13 + spall.melt_aware=1               -> rc != 0, output names the key and P2-26
''')
E(RK, "P2_26 R8", '''                "weibull.size_volume": "cell", "weibull.V_ref": "*"}}
''', '''                "weibull.size_volume": "cell", "weibull.V_ref": "*"},
         "R8": {"spall.surface_normal": "0", "spall.flake_coherence_length": "0",
                "spall.melt_aware": "0", "melt.lambda_crit": "*"}}
''')
E(RK, "MUTATION R8", '''                              '        pp.query_default("spall.flake_c_f"', 1))),
]''', '''                              '        pp.query_default("spall.flake_c_f"', 1))),
    ("R8", "re-add the parse of spall.melt_aware (a second occurrence of the key string)",
     _mmw(lambda t: t.replace('        pp.query_default("spall.flake_c_f"',
                              '        { int mo = 0; pp.query_default("spall.melt_aware", mo, 0); }\\n'
                              '        pp.query_default("spall.flake_c_f"', 1))),
]''')
E(RK, "K10-K13", '''            ("K9", ["surface_patch.pad_shield_open=0.5"], True, ["surface_patch.pad_shield_open",
                                                                 "specified but not used"])):''',
  '''            ("K9", ["surface_patch.pad_shield_open=0.5"], True, ["surface_patch.pad_shield_open",
                                                                 "specified but not used"]),
            ("K10", ["spall.surface_normal=0", "spall.flake_coherence_length=0", "spall.melt_aware=0",
                     "melt.lambda_crit=0.4"], False, []),
            ("K11", ["spall.surface_normal=1"], True, ["spall.surface_normal", "P2-26"]),
            ("K12", ["spall.flake_coherence_length=0.004"], True, ["spall.flake_coherence_length", "P2-26"]),
            ("K13", ["spall.melt_aware=1"], True, ["spall.melt_aware", "P2-26"])):''')

# ----------------------------------------------------------- lateral_flux_lock ----
E(LF, "the depth scan's stress regex", '''    STRESS = [("depth scan", rem, r"f_skel \\* \\(sp_relief \\* phase_E_loc"),''',
  '''    # P2-26: ScanSigmaT lost its f_skel factor (C3), so the scalar return is `sp_relief * phase_E_loc ...`
    STRESS = [("depth scan", rem, r"return sp_relief \\* phase_E_loc"),''')
E(LF, "the P2-26 name clause gains C2 + C3", '''               r"\\brc_ellq\\b|\\bWeibullSizeFactor\\b|\\bsp_size_flux\\b|\\bweibull_V_ref\\b"),
    }''', '''               r"\\brc_ellq\\b|\\bWeibullSizeFactor\\b|\\bsp_size_flux\\b|\\bweibull_V_ref\\b"),
        "C2": (r"\\bspall_surface_normal\\b|\\bsn_recession_superseded\\b|\\bSP_SN_COS_MIN\\b|\\bColumnSurfaceCos\\b|"
               r"\\bsn_cos_p\\b|\\bsn_cos_col\\b|\\bsn_cth\\b|\\bsn_icth\\b|\\bicth\\b|\\bsn_active\\b|"
               r"\\bspall_flake_coh_len\\b|\\bcoh_cells\\b"),
        # `melt_aware` also ends the key string "spall.melt_aware", which survives once in the table
        "C3": (r"\\bMeltSkeleton\\b|(?<!\\.)\\bmelt_aware\\b|\\bmelt_lambda_crit\\b|\\bf_skel\\b|\\bLL_col_prof\\b|"
               r"\\bLL_col_p\\b|\\bLLcol_i\\b|\\bLL_at_didx\\b"),
    }''')

# ------------------------------------------------------------------- face_heat ----
E(FH, "BASE tokens", '''    "surface_patch.pinned_idle_cycles=1.0e9", "surface.follow_mask=1", "spall.flake_coherence_length=0.0",
    "spall.surface_normal=0", "spall.surface_T=reconstructed", "spall.removal_energy=flake", "spall.flake_c_f=0.4",''',
  '''    # P2-26: G1's spall.flake_coherence_length=0.0 and spall.surface_normal=0 are dropped (both keys retired at
    # the string "0"; the argv must not write 0.0)
    "surface_patch.pinned_idle_cycles=1.0e9", "surface.follow_mask=1",
    "spall.surface_T=reconstructed", "spall.removal_energy=flake", "spall.flake_c_f=0.4",''')

# ----------------------------------------------------------- removal_increment ----
E(RI, "header variants", '''#   PRESC   ./input_presc                      prescribed_T, scan RETAINED
#   SN      ./input + spall.surface_normal=1   the superseded-recession witness
''', '''#   PRESC   ./input_presc                      prescribed_T, scan RETAINED
#   (P2-26: the SN leg -- spall.surface_normal=1, the superseded-recession witness -- and checks (e)-(e4), (f)
#   retire with the key, census C2.)
''')
E(RI, "header (e), (f)", '''#   (e) THE STEP 20 LOSS IS VISIBLE ON DISK, NOT JUST IN A LOG: SN carries the
#       thermo column sn_recession_superseded and W does not, and the parse-time
#       message names the superseded scaling. ⚠ A MESSAGE, NOT AN ABORT: two of
#       the four inputs that set surface_normal (unit/beam_void_closure and the
#       meier dev2d leg) are inside the 2639-file identity sweep, so aborting
#       would stop the sweep running them at all.
#   (f) THE WITNESS IS MEANINGFUL: sn_recession_superseded is a mean 1/cos(theta)
#       over firing columns, so it is >= 1 wherever it fires. Reading exactly 1
#       means nothing was lost on that step; < 1 would mean the estimator is
#       wrong.
''', "")
E(RI, "VARIANTS", '''    "PRESC": (INPUT_P, []),
    "SN":    (INPUT,   ["spall.surface_normal=1"]),
}

# The parse-time message must name what is superseded, not just that something is.
SN_MSG = "RECESSION scaling is superseded by whole-cell"
''', '''    "PRESC": (INPUT_P, []),
}
''')
E(RI, "checks (e)-(f)", '''    # ---- (e) the Step 20 loss is on disk ---------------------------------
    col = "sn_recession_superseded"
    check("(e) SN registers the superseded-recession witness", col in T["SN"],
          "present" if col in T["SN"] else "MISSING")
    check("(e2) a run without surface_normal does NOT register it",
          col not in T["W"], "absent on W" if col not in T["W"] else "PRESENT on W")
    log = open(os.path.join(OUT, "SN.log")).read()
    check("(e3) the parse-time message names the superseded scaling",
          SN_MSG in log, "found" if SN_MSG in log else f"MISSING: {SN_MSG!r}")
    check("(e4) surface_normal still RUNS (message, not abort)",
          len(E["SN"]) > 0, f"{len(E['SN'])} regime-1 rows")

    # ---- (f) the witness is meaningful -----------------------------------
    if col in T["SN"]:
        y = T["SN"][col]
        y = y[y > 0.0]
        check("(f) sn_recession_superseded is a 1/cos(theta) mean, so >= 1",
              y.size > 0 and np.all(y >= 1.0 - 1e-12),
              f"{y.size} firing rows, range [{y.min():.5f}, {y.max():.5f}]"
              if y.size else "no firing rows")

''', "")


def main():
    dry = "--apply" not in sys.argv
    text = {}
    bad = 0
    for path, what, old, new in EDITS:
        s = text.setdefault(path, open(path).read())
        n = s.count(old)
        print("  %-4s %-28s %-44s [%d match]" % ("ok" if n == 1 else "FAIL", os.path.relpath(path, U), what, n))
        if n == 1:
            text[path] = s.replace(old, new, 1)
        else:
            bad += 1
    src = os.path.join(U, "sp_melt_skeleton")
    dst = os.path.join(RU, "sp_melt_skeleton")
    mv_ok = os.path.isdir(src) and not os.path.exists(dst)
    print("  %-4s move unit/sp_melt_skeleton -> retired_units/sp_melt_skeleton" % ("ok" if mv_ok else "FAIL"))
    bad += 0 if mv_ok else 1
    if bad:
        print("FAIL: %d edit(s) did not match; nothing written" % bad)
        return 1
    if dry:
        print("DRY RUN: %d edits match exactly once; the move is possible" % len(EDITS))
        return 0
    for p, t in text.items():
        open(p, "w").write(t)
    shutil.move(src, dst)
    rd = os.path.join(RU, "README.md")
    lines = open(rd).read() if os.path.exists(rd) else (
        "# Retired unit tests\n\nTests of deleted features, kept as records (not run, not in any gate). One line each.\n\n")
    for nm, why in (("constraint_ref_lock", "P2-16: locked the constraint depth d (spallation.sp.constraint_depth), "
                     "deleted in P2-16"),
                    ("face_cref", "P2-16: the face list's constraint reference (P2-14a), deleted with d in P2-16"),
                    ("sp_melt_skeleton", "P2-26: A/B test of spall.melt_aware's skeleton f(Lambda_L) (P1), deleted "
                     "in P2-26 (census C3, user Q6); last on-disk run sweep_run.log (2026-09-17) PASS, the census "
                     "(2026-11-05) recorded rc 1 -- not re-run")):
        if "`%s`" % nm not in lines:
            lines += "- `%s` -- %s\n" % (nm, why)
    open(rd, "w").write(lines)
    print("APPLIED %d edits; moved sp_melt_skeleton; README.md updated" % len(EDITS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
