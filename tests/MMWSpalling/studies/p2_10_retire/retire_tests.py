#!/usr/bin/env python3
"""P2-10 item 6 -- the test edits, as a script (provenance, like retire.py).

  retire_tests.py <repo_root> <out_dir>
writes <out_dir>/{lateral_flux_lock,constraint_ref_lock}/test from the binary-B copies in
tests_bB/, with every edit an
asserted replacement, and <out_dir>/TEST_EDITS.txt: the declared list of the lock's
checks DELETED (their code is gone) and RESTATED (a count moved because the deleted
lateral site carried one of its members). Nothing is written into the repo.
"""
import os
import re
import sys

LOG = []
HERE = os.path.dirname(os.path.abspath(__file__))


def rep(t, what, old, new, count=1):
    n = t.count(old)
    assert n == count, "%s: %d occurrences of %r" % (what, n, old[:90])
    LOG.append("  %-8s %s" % ("RESTATED" if new.strip() else "DELETED", what))
    return t.replace(old, new)


def cut(t, what, start, end, keep_end=True):
    """delete from `start` (inclusive) to `end` (exclusive if keep_end)."""
    i = t.index(start)
    assert t.count(start) == 1, (what, t.count(start))
    j = t.index(end, i)
    if not keep_end:
        j += len(end)
    LOG.append("  %-8s %s" % ("DELETED", what))
    return t[:i] + t[j:]


# the names P2-10 deletes: every key, every function, every member a reader would grep for
DELETED_NAMES = (
    r"side_face_(?:flux|h|T|h_r_min|h_law|h_datum|factor_above|area_law|step_height|staircase)\b"
    r"|lateral_(?:drive|removal|subcell|thermal)\b|per_cell_removal|jet_wall_march|jet_wall_mdot"
    r"|jet_body|step_diag_radius|PerCellRemoval|FaceListFlux|ReconstructSurfaceT|fl_lat_top"
    r"|fl_lat_nfb|SideFaceH|SideFaceTgas|ColumnWallCorner|PatchAzimuthFraction|UpdateWallLawState"
    r"|WallJetH|WallJetBranch|WallZFill|WallMu\b|UpdateColumnPhase|AssertStepSplit"
    r"|WallStepDiagnostic|StepConductance|StepCellWeight|StaircaseFactor|AreaOfDelta|ColPhi"
    r"|PairDelta|side_T_eff|groupB_given|surf_side_P_(?:top|buried|rin|rout)|lat_accum_col"
    r"|lat_shift|ledger_E_flake_wall|sp_lateral_|pcr_on|wall_T_gas|jet_P_wall|wall_T_mouth")


def lock(t):
    t = rep(t, "docstring: the P2-10 scope note",
            '"""P1g-0b (a1)/(a2) -- ONE CRITERION, ONE KERNEL: the lateral flux argument lock.\n',
            '"""P1g-0b (a1)/(a2) -- ONE CRITERION, ONE KERNEL: the lateral flux argument lock.\n'
            "\n"
            "⚠⚠ P2-10: THE LATERAL SITES ARE DELETED. The face list's one breaking rule\n"
            "(FaceRemoval.H) replaced the lateral criterion, the side-face heat family and the\n"
            "wall march, and P2-10 removed them from src/. The blocks that locked those sites\n"
            "argument by argument are deleted with them (the declared list is\n"
            "studies/p2_10_retire/TEST_EDITS.txt); what stays locks the criterion, the\n"
            "reference, the property route and the face list -- and one new clause asserts\n"
            "that no deleted name survives in src/ outside the retired-key list.\n")
    t = cut(t, "SIG / LOCK: the four-site argument table (the lateral sites are deleted)",
            "# The signature, from MMWSpalling.H's declaration of SurfaceCellFlux. Positions",
            "def strip_comments(text):")
    t = cut(t, "calls(): the SurfaceCellFlux call parser (its only reader, the four-site lock, is deleted)",
            "def calls(path, text):", "def main():")
    t = cut(t, "P1g-0b (a2) the four lateral sites; (a1) the cap; the march's reason; (a3) the "
               "criterion's runtime aborts; P2-1 (b3) the SideFaceH resolution lock; the "
               "ColumnWallCorner datum",
            "    try:\n        cm, cr = calls(MMW, mmw), calls(REM, rem)",
            "    # ---- P2-2a (b4) + (b8): ONE RESTRAINT REFERENCE, NO SECOND FORM ------")
    # ---- P2-2a: the stress sites, 4 -> 3
    t = rep(t, "P2-2a: FOUR stress expressions -> THREE (the lateral criterion's is deleted)",
            '    STRESS = [("lateral criterion", rem, r"relief_loc \\* phase_E_loc \\* phase_beta_loc"),\n'
            '              ("depth scan", rem, r"f_skel \\* \\(sp_relief \\* phase_E_loc"),\n',
            '    # P2-10: the lateral criterion\'s expression is deleted with its site.\n'
            '    STRESS = [("depth scan", rem, r"f_skel \\* \\(sp_relief \\* phase_E_loc"),\n')
    t = rep(t, "P2-2a: found == 4 -> 3",
            '    need(found == 4,\n'
            '         "all FOUR surviving stress expressions carry a relief factor (%d) -- a "\n',
            '    need(found == 3,\n'
            '         "all THREE surviving stress expressions carry a relief factor (%d) -- a "\n')
    t = rep(t, "P2-2a: four E*beta products -> three",
            '    need(prod == 4,\n'
            '         "there are exactly FOUR E*beta products, so none is left un-relieved "\n',
            '    need(prod == 3,\n'
            '         "there are exactly THREE E*beta products, so none is left un-relieved "\n')
    t = rep(t, "P2-2a: `relief_loc` (Removal.H's alias lived in the lateral criterion) -> it is gone",
            '    for alias, txt, nm in (("te_relief", mmw, "MMWSpalling.H"),\n'
            '                           ("relief_loc", rem, "Removal.H")):\n',
            '    need(not re.search(r"\\brelief_loc\\b", rem),\n'
            '         "Removal.H has no `relief_loc` alias left (P2-10: it lived in the deleted "\n'
            '         "lateral criterion; the depth scan reads sp_relief itself)")\n'
            '    for alias, txt, nm in (("te_relief", mmw, "MMWSpalling.H"),):\n')
    # ---- P2-2d: the toughness sites 3 -> 2, the flaw reads 3 -> 2
    t = rep(t, "P2-2d clause 1: THREE toughness sites -> TWO",
            '    need(n_sites == 3,\n'
            '         "all THREE toughness sites EVALUATE the route (%d) -- MMWSpalling.H\'s "\n'
            '         "Sp_field site, Removal.H\'s lateral face, and the depth scan\'s ONE "\n'
            '         "ScanKIcPa (P2-9 joined its two textually identical blocks)" % n_sites)\n',
            '    # ⚠ P2-10 RESTATED: 3 -> 2. The lateral face was the third site; it is\n'
            '    # deleted with the lateral criterion.\n'
            '    need(n_sites == 2,\n'
            '         "both toughness sites EVALUATE the route (%d) -- MMWSpalling.H\'s "\n'
            '         "Sp_field site and the depth scan\'s ONE ScanKIcPa" % n_sites)\n')
    t = rep(t, "P2-2d clause 2: THREE flaw reads -> TWO (the lateral face's is deleted)",
            '        ("lateral face",         rem, "a_f *= WeibullSizeFactor(ell * ell * ell)"),\n', "")
    # ---- P2-2b: the LawSigma sites, the verbatim lateral expression
    t = rep(t, "P2-2b: `Removal.H lateral` routes through LawSigma (the site is deleted)",
            '                            ("Removal.H lateral",    rem, "LawSigma(lat_g"),\n', "")
    t = rep(t, "P2-2b: `C_th * jump` kept verbatim at the lateral site (the site is deleted)",
            '    need("C_th * jump" in rem,\n'
            '         "the lateral site keeps `C_th * jump` verbatim for `law = constant`")\n', "")
    # ---- the runtime block: the live control
    t = rep(t, "P2-8 (a2): the live control `side_face_flux=0` -> `face_list=0` (side_face_flux is "
               "retired by P2-10 and now aborts)",
            '        ab, txt = aborts_for(["surface_patch.side_face_flux=0"])\n'
            '        print("       %-42s %s" % ("surface_patch.side_face_flux=0 (live)",\n',
            '        # ⚠ P2-10: side_face_flux is itself retired now, so the live key is face_list.\n'
            '        ab, txt = aborts_for(["surface_patch.face_list=0"])\n'
            '        print("       %-42s %s" % ("surface_patch.face_list=0 (live)",\n')
    t = cut(t, "P2-3 (a3) the gas-temperature lock (SideFaceTgas, five readers, its negative "
               "control); P2-4 (a2) the cold body (jet_body_h_scale x wall_h_duct, its negative control)",
            "    # ---- P2-3 (a3): THE GAS-TEMPERATURE LOCK, FIVE READERS ---------------",
            "    # ---- P2-8 (a3): THE FACE LIST -- one builder, one heat rule, one area rule")
    # ---- P2-8 (a3): FaceListFlux, the march, the criterion, group B
    t = rep(t, "P2-8 (a3): FaceListFlux declared once, 0 calls -> deleted (0 declarations)",
            '            ("FaceWeight", r"static\\s+Set::Scalar\\s+FaceWeight\\s*\\(", 1,\n'
            '             "only inside BuildFaceList -- the ONE area rule"),\n'
            '            # ⚠ P2-9 RESTATED: F4\'s Newton BC branch is NOT CALLED under the key\n'
            '            # -- every listed face takes the cell\'s T_s through FrFaceFlux\n'
            '            # (checked below). P2-10 deletes FaceListFlux.\n'
            '            ("FaceListFlux", r"static\\s+void\\s+FaceListFlux\\s*\\(", 0,\n'
            '             "none -- P2-9 replaced it at all four sites by FrFaceFlux")):\n',
            '            ("FaceWeight", r"static\\s+Set::Scalar\\s+FaceWeight\\s*\\(", 1,\n'
            '             "only inside BuildFaceList -- the ONE area rule")):\n')
    t = rep(t, "P2-8 (a3): (FaceListFlux) -- now asserted DELETED",
            '    # the heat rule is the ONLY place h_expr is evaluated per face\n'
            if False else
            '    # the heat resolver is the ONLY place h_expr is evaluated per face\n',
            '    # P2-10: P2-8\'s face Newton (F4) had no caller after P2-9 and is deleted.\n'
            '    need(ndecl(r"void\\s+FaceListFlux\\s*\\(") == 0,\n'
            '         "`FaceListFlux` is DELETED (P2-10; FrFaceFlux heats every listed face)")\n'
            '    # the heat resolver is the ONLY place h_expr is evaluated per face\n')
    t = rep(t, "P2-8 (a3): the march builds no slot gather UNDER THE KEY -> it has none at all",
            '    need(mar is not None and "fl_col_first[c]" in mar.group(0)\n'
            '         and re.search(r"side_on\\s*=\\s*surface_patch\\.side_face_flux\\s*&&\\s*!fl_on_m",\n'
            '                       mar.group(0)) is not None,\n'
            '         "the MARCH bills the list\'s faces per column and builds no slot gather "\n'
            '         "under the key (slot builder and side term both replaced)")\n',
            '    need(mar is not None and "fl_col_first[c]" in mar.group(0)\n'
            '         and "side_nslot" not in mar.group(0) and "side_on" not in mar.group(0),\n'
            '         "the MARCH bills the list\'s faces per column and has no slot gather and "\n'
            '         "no side term at all (P2-10 deleted both)")\n')
    t = rep(t, "P2-8 (a3): the CRITERION takes exposure and h from the list; its q_net_l (the "
               "lateral criterion is deleted)",
            '    need("fla(i, j, k, 0)" in code_r and "h_dir[e.dir] = e.h" in code_r\n'
            '         and re.search(r"if\\s*\\(k\\s*==\\s*kT\\)\\s*continue;", code_r) is not None,\n'
            '         "the CRITERION takes exposure and h from the list and leaves the column-top "\n'
            '         "cell to the vertical criterion (the transitional rule)")\n'
            '    need(re.search(r"h_q\\s*\\*\\s*\\(T_face_gas\\s*-\\s*T_face\\)", code_r) is not None,\n'
            '         "the criterion\'s q_net_l uses the FACE\'s h and gas under the key")\n', "")
    t = rep(t, "P2-8 (a3): the kernel's and ledger's side blocks are OFF under the key -> they "
               "do not exist",
            '    # group B is not read under the key at the old side sites\n'
            '    need(re.search(r"r2_col <= sp_r2 && !fl_on;", code_m) is not None\n'
            '         and re.search(r"side_on\\s*=\\s*sp_on\\s*&&\\s*surface_patch\\.side_face_flux\\s*&&\\s*!fl_on_l",\n'
            '                       code_m) is not None,\n'
            '         "the kernel\'s and the ledger\'s old side-face blocks are OFF under the key")\n',
            '    # P2-10: the old side sites are deleted, not merely off under the key\n'
            '    need(re.search(r\'(?<!")surface_patch\\.side_face_flux\', code_m) is None\n'
            '         and re.search(r"\\bside_on\\b", code_m) is None,\n'
            '         "the kernel\'s and the ledger\'s old side-face blocks are DELETED (P2-10)")\n')
    # ---- P2-9: the dispatch order (the lateral criterion is gone)
    t = rep(t, "P2-9: the dispatch precedes the lateral criterion and its transitional skip -> "
               "precedes the column depth scan, and no lateral criterion exists",
            '    i_lat = ur.find("const bool lat_on = sp_lateral_drive")\n'
            '    i_skip = re.search(r"if\\s*\\(k\\s*==\\s*kT\\)\\s*continue;", ur)\n'
            '    i_scan = ur.find("h_col = DepthScanH(")\n'
            '    need(0 <= i_disp < i_ret < i_lat and i_skip is not None and i_ret < i_skip.start()\n'
            '         and i_ret < i_scan and "RecedeAlongLines(" in ur[i_disp:i_ret],\n'
            '         "the face-list DISPATCH (RecedeAlongLines; return) precedes the lateral criterion, "\n'
            '         "its transitional skip `if (k == kT) continue;` and the column depth scan")\n',
            '    i_scan = ur.find("h_col = DepthScanH(")\n'
            '    need(0 <= i_disp < i_ret < i_scan and "RecedeAlongLines(" in ur[i_disp:i_ret]\n'
            '         and "lat_on" not in ur,\n'
            '         "the face-list DISPATCH (RecedeAlongLines; return) precedes the column depth "\n'
            '         "scan, and no lateral criterion is left after it (P2-10)")\n')
    t = rep(t, "P2-9 rev 2: the column tail `if (remap_any || lat_shift_any)` -> `if (remap_any)`",
            '    tail = re.search(r"if\\s*\\(remap_any\\s*\\|\\|\\s*lat_shift_any\\)\\s*\\{(.*?)\\}", ur2, re.S)\n',
            '    tail = re.search(r"if\\s*\\(remap_any\\)\\s*\\{(.*?)\\}", ur2, re.S)\n')
    # ---- P2-10: the new clause
    t = rep(t, "P2-10 (new): no deleted name survives in src/ outside the retired-key list",
            '    print()\n    if fails:\n',
            '    # ---- P2-10: NOTHING DELETED SURVIVES ------------------------------------\n'
            '    # The retired-key list (MMWSpalling::Parse) is the ONE place the deleted keys\n'
            '    # may be named -- it is what makes them abort. Everywhere else in src/, code\n'
            '    # and comments alike, a deleted key, function or member is a leftover.\n'
            '    print("\\nP2-10 NOTHING DELETED SURVIVES IN src/ (outside the retired-key list)")\n'
            '    DELETED = re.compile(\n'
            '        ' + repr(DELETED_NAMES) + ')\n'
            '    left = []\n'
            '    for dp, _d, fs in os.walk(os.path.join(ROOT, "src")):\n'
            '        for fn in fs:\n'
            '            if not fn.endswith((".H", ".h", ".cpp", ".cc")):\n'
            '                continue\n'
            '            p = os.path.join(dp, fn)\n'
            '            in_list = False\n'
            '            for ln, line in enumerate(open(p, errors="replace"), 1):\n'
            '                if "P2-10: THE RETIRED KEYS" in line:\n'
            '                    in_list = True\n'
            '                if in_list and "remove it from the input" in line:\n'
            '                    in_list = False\n'
            '                    continue\n'
            '                if in_list:\n'
            '                    continue\n'
            '                if DELETED.search(line):\n'
            '                    left.append("%s:%d %s" % (os.path.relpath(p, ROOT), ln, line.strip()[:70]))\n'
            '    for o in left[:20]:\n'
            '        print("       ** %s" % o)\n'
            '    need(not left, "no deleted key, function or member is named anywhere in src/ outside "\n'
            '                   "the retired-key list (%d line(s))" % len(left))\n'
            '    nlist = sum(open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")).read()\n'
            '                .count(k) for k in ("\\"spall.lateral_drive\\"", "\\"surface_patch.side_face_flux\\"",\n'
            '                                    "\\"surface_patch.jet_wall_march\\""))\n'
            '    need(nlist == 3, "the retired-key list itself exists and names them (%d of 3 probes)" % nlist)\n'
            '\n'
            '    print()\n    if fails:\n')
    t = rep(t, "the PASS message",
            '    print("PASS: the four lateral flux sites agree argument by argument, all four\\n"\n'
            '          "      resolve `h` PER FACE through the one resolver, the stress-free\\n"\n'
            '          "      reference is the STORED INITIAL FIELD -- write-once, unplotted,\\n"\n'
            '          "      read at BOTH firing sites -- the\\n"\n'
            '          "      second driving form and the support rule are DELETED rather than\\n"\n'
            '          "      defaulted off, their keys go UNCONSUMED while live keys still\\n"\n'
            '          "      bind, no strain table has entered the criterion, and all five\\n"\n'
            '          "      gas-temperature readers -- the criterion among them -- resolve\\n"\n'
            '          "      through the one height-aware resolver.")\n',
            '    print("PASS: the stress-free reference is the STORED INITIAL FIELD -- write-once,\\n"\n'
            '          "      unplotted, read at BOTH firing sites -- the second driving form and\\n"\n'
            '          "      the support rule are DELETED rather than defaulted off, no strain\\n"\n'
            '          "      table has entered the criterion, the face list is one builder, one\\n"\n'
            '          "      heat rule and one breaking rule, and nothing P2-10 deleted survives\\n"\n'
            '          "      in src/ outside the retired-key list.")\n')
    return t


def cref_lock(t):
    """unit/constraint_ref_lock (P2-2e): found at a8, not in the packet's item 6 list.
    It is a source lock like lateral_flux_lock, and three of its clauses name the lateral
    criterion's cref branch and the Sp_field resolver's lateral path, both deleted."""
    LOG.append("")
    LOG.append("unit/constraint_ref_lock (found at a8 -- a source lock the packet's item 6 did not list)")
    t = rep(t, "docstring: the lateral criterion's two reads (the criterion is deleted)",
            "  Removal.H  lateral criterion   law  LawSigma(lat_g, T0 + jump, ...)\n"
            "                                 verbatim  C_th * (...)\n",
            "  (Removal.H  lateral criterion  law + verbatim: DELETED in P2-10 with the\n"
            "   criterion; its clause C4 and mutation control went with it)\n")
    t = cut(t, "C4 the lateral criterion's cref branch (law + verbatim): the criterion is deleted",
            "    # C4 the lateral criterion\n", "    # C5 the Sp_field site")
    t = rep(t, "C5 restated: the resolver's LATERAL path (the exp(-(depth+d)/ell_s) profile model "
               "for nx_in/ny_in faces) is deleted; C5 now asserts it is GONE (vertical-only resolver)",
            "its LATERAL path the profile model d further in; d = 0 is T0_at(...)\",",
            "no LATERAL path is left (deleted in P2-10); d = 0 is T0_at(...)\",")
    t = rep(t, "C5 restated (condition)",
            '         and "std::exp(-(depth_n * icth + cref_d_sp) / ell_s)" in tb)',
            '         and "ell_s" not in tb and "nx_in" not in tb)')
    t = rep(t, "C10 restated: its anchor pcr_lateral_cells (the per-cell census, deleted) -> the "
               "registration that now precedes the cref block, sn_recession_superseded",
            "    i_p = mmw.find('\"pcr_lateral_cells\",  false);')",
            "    i_p = mmw.find('\"sn_recession_superseded\", false);')")
    t = cut(t, "mutation m_lat (deletes the lateral criterion's cref branch): its target is gone",
            "def m_lat(mmw, rem):", "def m_sp(mmw, rem):")
    t = rep(t, "MUTATIONS: the m_lat row",
            '    ("delete the lateral criterion\'s cref branch", m_lat, "C4"),\n', "")
    return t


def main():
    root, out = sys.argv[1], sys.argv[2]
    src = os.path.join(HERE, "tests_bB", "lateral_flux_lock", "test")  # the binary-B copy
    t = lock(open(src).read())
    os.makedirs(os.path.join(out, "lateral_flux_lock"), exist_ok=True)
    p = os.path.join(out, "lateral_flux_lock", "test")
    open(p, "w").write(t)
    os.chmod(p, 0o755)
    src = os.path.join(HERE, "tests_bB", "constraint_ref_lock", "test")  # the binary-B copy
    t = cref_lock(open(src).read())
    os.makedirs(os.path.join(out, "constraint_ref_lock"), exist_ok=True)
    p = os.path.join(out, "constraint_ref_lock", "test")
    open(p, "w").write(t)
    os.chmod(p, 0o755)
    hdr = ("P2-10 item 6 -- unit/lateral_flux_lock and unit/constraint_ref_lock: the checks of "
           "deleted code, and the counts that moved\n" + "=" * 100 + "\n")
    open(os.path.join(out, "TEST_EDITS.txt"), "w").write(hdr + "\n".join(LOG) + "\n")
    print(hdr + "\n".join(LOG))


if __name__ == "__main__":
    main()
