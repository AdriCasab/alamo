#!/usr/bin/env python3
"""P2-10 -- the deletion, as a script (provenance: every cut is an anchored range).

  retire.py c1 <in_dir> <out_dir>   binary C1: the lateral path, per-cell removal, dead code
  retire.py c2 <in_dir> <out_dir>   binary C2: side-face heat (group B + side_face_flux) and
                                    the wall march + cold body
<in_dir>/<out_dir> hold MMWSpalling.H, Removal.H, FaceRemoval.H (flat). Every cut asserts
its anchors (unique, or the first after a unique anchor) and the line counts it removes, so
a source that moved under the script stops it instead of cutting the wrong lines.
Writes <out_dir>/CUTS.txt: one row per cut, lines removed, and the first/last line cut.
"""
import os
import sys

FILES = {"M": "MMWSpalling.H", "R": "Removal.H", "F": "FaceRemoval.H"}


class Doc:
    def __init__(self, path, tag, log):
        self.L = open(path).read().split("\n")
        self.tag, self.log = tag, log

    def idx(self, anchor, after=0):
        hits = [i for i, l in enumerate(self.L) if anchor in l]
        if after:
            hits = [i for i in hits if i >= after]
            assert hits, "%s: anchor not found after %d: %r" % (self.tag, after, anchor)
            return hits[0]
        assert len(hits) == 1, "%s: anchor %r found %d times %s" % (self.tag, anchor, len(hits), hits[:5])
        return hits[0]

    def cut(self, what, start, end=None, end_off=0, after=None, keep_end=False,
            start_off=0, insert=None):
        """delete lines [start + start_off .. end + end_off] (end found at or after
        start); with keep_end, stop BEFORE the end anchor line (end_off then counts
        back). `insert` (a string) replaces the cut lines."""
        a = self.idx(after) if after else 0
        s = (self.idx(start, a) if after else self.idx(start)) + start_off
        if end is None:
            e = s + end_off
        else:
            e = self.idx(end, s) + end_off
            if keep_end:
                e -= 1
        assert e >= s, (what, s, e)
        gone = self.L[s:e + 1]
        ins = insert.split("\n") if insert is not None else []
        self.log.append("%s  %-60s -%4d  [%s ... %s]%s" % (self.tag, what, len(gone),
                                                       gone[0].strip()[:50], gone[-1].strip()[:50],
                                                       ("  +%d inserted" % len(ins)) if ins else ""))
        self.L[s:e + 1] = ins

    def rep(self, what, old, new, count=1):
        t = "\n".join(self.L)
        n = t.count(old)
        assert n == count, "%s: %s: %d occurrences of %r" % (self.tag, what, n, old[:80])
        n0 = len(self.L)
        self.L = t.replace(old, new).split("\n")
        self.log.append("%s  %-60s %+5d  (replace)" % (self.tag, what, len(self.L) - n0))

    def text(self):
        return "\n".join(self.L)


def load(d, log):
    return {k: Doc(os.path.join(d, f), k, log) for k, f in FILES.items()}


def save(D, out, log):
    os.makedirs(out, exist_ok=True)
    for k, f in FILES.items():
        open(os.path.join(out, f), "w").write(D[k].text())
    open(os.path.join(out, "CUTS.txt"), "w").write("\n".join(log) + "\n")


# ======================================================================= C1 ==
RETIRED_C1 = ["spall.lateral_drive", "spall.lateral_removal", "spall.lateral_subcell",
              "spall.lateral_thermal", "spall.per_cell_removal"]


def retired_block(keys):
    names = "".join('                                "%s",\n' % k for k in keys)
    return (
        "        // P2-10: THE RETIRED KEYS. The flame wall machinery the face list made\n"
        "        // redundant -- the side-face heat family, the lateral criterion (regime 6,\n"
        "        // banking, the thermal shift), per-cell removal and the wall march with\n"
        "        // its cold body -- is deleted. A deck that still names one of its keys\n"
        "        // aborts here, naming the key, instead of running with it ignored\n"
        "        // (`spall.*` sits outside any strict ParmParse scope, so this list is\n"
        "        // the only guard those keys have).\n"
        "        for (const char* k : {\n" + names.rstrip(",\n") + "})\n"
        "            if (pp.contains(k))\n"
        "                Util::Abort(INFO, k, \" was retired in P2-10 (the flame wall machinery: \"\n"
        "                            \"side-face heat, the lateral criterion, per-cell removal and \"\n"
        "                            \"the wall march were deleted); remove it from the input\");\n")


def c1(D):
    R, M, F = D["R"], D["M"], D["F"]
    # ---------------------------------------------------------------- Removal.H
    R.cut("lateral criterion: declarations",
          "// ================= P1f: THE LATERAL CRITERION (deliverable 2) =========",
          "long lat_len_fb_h0 = 0, lat_len_fb_q0 = 0;")
    R.rep("cref: the lateral site's counter",
          "        // ⚠⚠ P2-2e: the constraint depth at BOTH Removal.H sites (the lateral\n"
          "        // criterion and the depth scan). `sp_cref_d > 0` alone is the switch:\n"
          "        // Parse aborts unless d > 0 comes with the T-derived criterion path, and\n"
          "        // `sp_weibull_use_T_face` is not formed until much later in this\n"
          "        // function -- AFTER the lateral criterion, which must read it too.\n"
          "        // Witnesses are accumulated per BOX (the scan runs a column only on the\n"
          "        // rank owning its top, :1497; the lateral loop is per box) and reduced\n"
          "        // ONCE below the scan.\n"
          "        const bool        cref_on    = sp_cref_d > 0.0;\n"
          "        const Set::Scalar cref_d_loc = sp_cref_d;\n"
          "        long cref_vert = 0, cref_lat = 0, cref_clamp = 0, cref_neg = 0;\n",
          "        // ⚠⚠ P2-2e: the constraint depth at the depth scan. `sp_cref_d > 0`\n"
          "        // alone is the switch: Parse aborts unless d > 0 comes with the\n"
          "        // T-derived criterion path. Witnesses are accumulated per BOX (the scan\n"
          "        // runs a column only on the rank owning its top) and reduced ONCE below\n"
          "        // the scan. (P2-10 deleted the lateral criterion, the second site here,\n"
          "        // so cref_n_lat reads 0.)\n"
          "        const bool        cref_on    = sp_cref_d > 0.0;\n"
          "        const Set::Scalar cref_d_loc = sp_cref_d;\n"
          "        long cref_vert = 0, cref_clamp = 0, cref_neg = 0;\n")
    R.cut("lateral criterion: the pass and its reductions",
          "// ⚠⚠ deliverable 12: taken over EVERY candidate, before the firing",
          "// Step 19: sigma_xx-reconstruction selector for k_i_at_zcrack.",
          end_off=-1, keep_end=True)
    R.rep("toughness locals only the lateral criterion read",
          "        const KIcLaw      law_loc          = kic_law;\n"
          "        // ⚠ P2-2d d7: the law object for `grimsel_table`. Hoisted once -- the\n"
          "        // three sites below must not each decide where the toughness comes\n"
          "        // from, which is how the Pa-vs-MPa unit bug and the h_side cap bug\n"
          "        // both arose in this file.\n"
          "        const Numeric::Material::Grimsel *kic_g_loc =\n"
          "            phases.empty() ? nullptr : &phases[0].grimsel;\n"
          "        const Set::Scalar kic_scale_loc    = kic_scale;\n"
          "        const Set::Scalar kic_K_Ic0_loc    = kic_K_Ic0;\n",
          "")
    R.rep("col_rho: its consumer is the pin now",
          "        // ⚠ Gathered UNCONDITIONALLY (not under `rev_on`): site 3's consumer is\n"
          "        // the lateral banking energy estimate, which runs in a per-column loop\n"
          "        // with no MFIter in scope and under a DIFFERENT gate. Coupling it to\n"
          "        // `rev_on` would tie two unrelated conditions together.\n",
          "        // Gathered unconditionally (not under `rev_on`); the pin reads it below.\n")
    R.cut("lat_csv (event-log gate)",
          "// P1f deliverable 4: gates the three APPENDED event-log columns. ⚠ Keyed",
          "const bool lat_csv = sp_lateral_drive;")
    R.rep("the flaw reads: the lateral one is gone",
          "                            // THE THREE READS, and their equality BY\n"
          "                            // CONFIGURATION:\n"
          "                            //   MMWSpalling.H:7403  criterion / Sp_field  scaled\n"
          "                            //   Removal.H:823       lateral face          scaled\n"
          "                            //   HERE                the thickness scan    scaled now\n",
          "                            // THE READS, and their equality BY\n"
          "                            // CONFIGURATION (P2-10 deleted the third, the\n"
          "                            // lateral face's):\n"
          "                            //   MMWSpalling.H:7403  criterion / Sp_field  scaled\n"
          "                            //   HERE                the thickness scan    scaled now\n")
    R.rep("cref: the reduction carries 0 for the lateral share",
          "        // P2-2e: the depth scan's and the lateral criterion's share of the\n"
          "        // constraint-reference witnesses, added to what the Sp_field site\n"
          "        // contributed this step. ⚠ `cref_on` is global, so every rank reaches\n"
          "        // this collective or none does.\n"
          "        if (cref_on)\n"
          "        {\n"
          "            long cnt[4] = {cref_vert, cref_lat, cref_clamp, cref_neg};\n",
          "        // P2-2e: the depth scan's share of the constraint-reference witnesses,\n"
          "        // added to what the Sp_field site contributed this step. ⚠ `cref_on`\n"
          "        // is global, so every rank reaches this collective or none does.\n"
          "        if (cref_on)\n"
          "        {\n"
          "            long cnt[4] = {cref_vert, 0, cref_clamp, cref_neg};\n")
    R.cut("regime 6 in the winner loop",
          "// ⚠⚠ P1f: REGIME 6 -- a LATERAL firing. h_c is deliberately ZERO:",
          "h_c        = 0.0;   // ⚠ NOT dz: phi must not move for regime 6", end_off=1)
    R.cut("regime 6 ledger tallies",
          "// ⚠⚠ P1f: REGIME 6 OWES ITS OWN LEDGER TERM. Regime 5's omission took a",
          "long lat_voided = 0;")
    R.cut("banking (P1g-1)",
          "// P1g-1: banking tallies. lat_void_col_local is PER COLUMN because the",
          "lat_bank_step[idx] = hf;   // P1g-2: what the shift must displace", end_off=2)
    R.cut("regime 6 in PASS 3",
          "// ⚠⚠ P1f: REGIME 6 -- void EXACTLY the cell the lateral",
          "// P0b deliverable 4: write the remapped profile. SECOND loop,", keep_end=True)
    R.cut("regime 6 void-count assertion",
          "// ⚠ P1f: regime 6 voids exactly one cell, so D2z-1's",
          "\"instead bank and void none.\");", end_off=1)
    R.cut("lateral ledger + the thermal shift (P1g-2) + banking diagnostics",
          "// ⚠⚠ P1f: THE LATERAL LEDGER TERM IS UNCONDITIONAL, NOT INSIDE `remap_on`.",
          "// A1: every-removal event log, written by the I/O rank after the void",
          end_off=-1, keep_end=True)
    R.cut("event log header: the lateral columns",
          "// P1f: riser_depth is how far BELOW the column's top the",
          "\"lat_nfire\" : \"\")")
    R.rep("event log support-rule row: the lateral zeros",
          "                                << (recon_on ? \",0,0,0,0,0,0,0,0\" : \"\")\n"
          "                                << (lat_csv ? \",0,0,0,0,0,0,0,0\" : \"\") << \"\\n\";\n",
          "                                << (recon_on ? \",0,0,0,0,0,0,0,0\" : \"\") << \"\\n\";\n")
    R.cut("event log rows: the lateral columns",
          "                        if (lat_csv)",
          "<< \",\" << (lat_nf_log.empty()  ? 0.0 : lat_nf_log[idx]);")
    R.cut("the shift's InvertHtoT note",
          "// ⚠⚠ P1g-2 deliverable 7 rides on EXACTLY this hook, and that is the",
          "lat_shift_any is built from reduced", end_off=1)
    R.rep("the remap hook without the shift", "        if (remap_any || lat_shift_any)\n",
          "        if (remap_any)\n")
    # ---------------------------------------------------------- FaceRemoval.H
    F.cut("fr_lateral_given", "std::vector<std::string> fr_lateral_given;")
    F.cut("lat_csv", "const bool lat_csv = sp_lateral_drive;")
    F.cut("event log header: lateral columns",
          "<< (lat_csv ? \",riser_depth,ell_lat,T_face_lat,\"", "\"lat_nfire\" : \"\")")
    F.cut("event log rows: lateral zeros", "if (lat_csv) ofs << \",0,0,0,0,0,0,0,0,0,0,0\";")
    # ---------------------------------------------------------- MMWSpalling.H
    M.rep("the retired-key list (C1), first in Parse",
          "    static void Parse(MMWSpalling& value, IO::ParmParse& pp)\n    {\n",
          "    static void Parse(MMWSpalling& value, IO::ParmParse& pp)\n    {\n"
          + retired_block(RETIRED_C1))
    M.rep("LEN counters' placement note",
          "        //   - inside `if (sp_lateral_drive)` (within the side-face-flux guard):\n",
          "        //   - inside the lateral-drive guard (deleted in P2-10):\n")
    M.cut("lat_len_fb registrations", "value.RegisterIntegratedVariable(&value.lat_len_fb_h0_n,",
          end_off=3)
    M.cut("lateral registrations (P1f/P1g-1/P1g-2)",
          "// P1f: the lateral mode's own counters, registered only under the",
          "\"lat_void_colmax\", false);", end_off=1)
    M.cut("per-cell removal preconditions + registrations",
          "// D2q-2b item 5: per-cell removal preconditions, cross-checked here",
          "\"pcr_lateral_cells\",  false);", end_off=3)
    M.cut("fl_lat_top registration", "value.RegisterIntegratedVariable(&value.fl_lat_top, \"fl_lat_top\", false);")
    M.cut("fl_lat_nfb registration", "value.RegisterIntegratedVariable(&value.fl_lat_nfb, \"fl_lat_nfb\", false);")
    M.cut("face-list abort for per_cell_removal", "if (value.spall_per_cell_removal) fr_abort(")
    M.cut("fr_lateral_given in the unread-keys line", "for (const auto& k : value.fr_lateral_given)")
    M.cut("PerCellRemoval call", "// D2q-2b item 6: lateral per-cell removal, AFTER the column pass so the",
          "PerCellRemoval(lev, time, dt);")
    M.rep("F3-F5 (stale since P2-9 rev 2)",
          "    // F3 (FaceWeight)  w = 1/(|nx|+|ny|+|nz|), n the unit outward normal from\n"
          "    //     the central-difference gradient of the solid indicator\n"
          "    //     (Numeric::Gradient) over the cell's 3x3x3 neighbourhood. Flat floor:\n"
          "    //     w = 1 exactly. 45 deg staircase: 1/sqrt2 per face, summing to the\n"
          "    //     plane. Zero gradient => w = 1, counted (fl_w_fallback). One w per\n"
          "    //     CELL, applied to each of its faces.\n"
          "    // F4  the column-top +z face keeps the surface form (SurfaceRobinT) x w;\n"
          "    //     every other face takes the face Newton x w (FaceListFlux). ⚠ The ONE\n"
          "    //     transitional rule of P2-8, retired by P2-9.\n"
          "    // F5  the kernel, the ledger, the march (billing each face to its column's\n"
          "    //     bin) and the lateral criterion read F1-F4 and nothing else.\n",
          "    // F3 (FaceWeight)  w = 1/(|nx|+|ny|+|nz|), n the unit outward normal of the\n"
          "    //     SMOOTHED solid fraction (FrSmoothedGradient, P2-9 rev 2). Flat floor:\n"
          "    //     w = 1 exactly. 45 deg staircase: 1/sqrt2 per face, summing to the\n"
          "    //     plane. Zero gradient => w = 1, counted (fl_w_fallback). One w per\n"
          "    //     CELL, applied to each of its faces.\n"
          "    // F4  every listed face -- the column top included -- is heated at its\n"
          "    //     cell's one T_s, x w (FrFaceFlux, P2-9).\n"
          "    // F5  the kernel, the ledger, the march (billing each face to its column's\n"
          "    //     bin) and the line rule read F1-F4 and nothing else.\n")
    M.rep("fl_lat members",
          "    Set::Scalar fl_P_march = 0.0, fl_w_fallback = 0.0, fl_lat_top = 0.0;\n"
          "    Set::Scalar fl_h_top_mismatch = 0.0, fl_lat_nfb = 0.0;\n",
          "    Set::Scalar fl_P_march = 0.0, fl_w_fallback = 0.0;\n"
          "    Set::Scalar fl_h_top_mismatch = 0.0;\n")
    M.cut("FaceListFlux (0 callers)",
          "// F4 for every face that is NOT the column-top +z face: the face Newton.",
          "q_in, q_out, q_rob, q_los, 0.0, nullptr, 0.0);", end_off=2)
    M.cut("lateral key parse (1 of 2)", "// P1f: the lateral mode. Default 0 = INERT.",
          "pp.query_default(\"spall.lateral_thermal\", value.sp_lateral_thermal, false);")
    M.cut("lateral key parse (2 of 2) + the P1f self-test + the staircase abort",
          "// ⚠⚠ DELIBERATELY NOT an abort-if-unaccompanied, unlike lateral_removal",
          "\"measuring the grid.\");")
    M.cut("per_cell_removal parse",
          "// D2q-2b item 5: per-cell removal. Default 0 = off and byte-identical",
          "value.spall_per_cell_removal = (pcr == 1);", end_off=1)
    # UpdateSpAfterMechanics: the per-cell branch
    M.rep("Sp site: no lateral cref counter",
          "        long sp_cref_vert = 0, sp_cref_lat = 0, sp_cref_clamp = 0, sp_cref_neg = 0;\n",
          "        long sp_cref_vert = 0, sp_cref_clamp = 0, sp_cref_neg = 0;\n")
    M.cut("Sp site: pcr_on + its ghost refill",
          "// D2q-2b item 6: per-cell removal scores every exposed face, so the",
          "if (pcr_on) Util::RealFillBoundary(*removed_mf[lev], geom[lev]);", end_off=1)
    M.cut("Sp site: rem_pc placeholder",
          "// Bound to Sp_field_mf as a harmless placeholder when the key is",
          ": (*Sp_field_mf[lev]).const_array(mfi);")
    M.rep("Sp site: the top face only",
          "                // D2q-2b item 6: with per-cell removal the criterion is\n"
          "                // evaluated on EVERY exposed face, not only the column top.\n"
          "                // `nx_in, ny_in, nz_in` is the inward normal of the face being\n"
          "                // scored -- (0,0,-1) is the top face, i.e. exactly the\n"
          "                // pre-D2q-2b behaviour, which is why the key-off path below is\n"
          "                // the original code unchanged.\n"
          "                int nx_in = 0, ny_in = 0, nz_in = -1;\n"
          "                if (!pcr_on)\n"
          "                {\n"
          "                    if (k != ksurf) return;\n"
          "                }\n",
          "                if (k != ksurf) return;\n"
          "                // (P2-10: the per-cell lateral branch is deleted.)\n"
          "                if (false)\n"
          "                {\n"
          "                    if (k != ksurf) return;\n"
          "                }\n")
    M.cut("Sp site: the per-cell lateral branch",
          "                // (P2-10: the per-cell lateral branch is deleted.)",
          "else           { ny_in = -1; }", end_off=2)
    M.rep("Sp site: recon on the top face (now the only face)",
          "                // P0: reconstruct the SURFACE temperature from the top two\n"
          "                // cell averages. ⚠ Applied to the TOP FACE only: under\n"
          "                // spall.per_cell_removal a lateral face has no column profile,\n"
          "                // and per-cell removal stays off in this programme -- so a\n"
          "                // lateral face keeps the cell average and that is stated, not\n"
          "                // silently assumed. Needs a solid second cell; when the\n"
          "                // contiguous solid run below the top is a single cell there is\n"
          "                // no profile to read and the cell average stands (counted).\n"
          "                SurfaceRecon rec;\n"
          "                bool rec_ok = false;\n"
          "                if (recon_on && nz_in == -1 && nx_in == 0 && ny_in == 0)\n",
          "                // P0: reconstruct the SURFACE temperature from the top two\n"
          "                // cell averages. Needs a solid second cell; when the\n"
          "                // contiguous solid run below the top is a single cell there is\n"
          "                // no profile to read and the cell average stands (counted).\n"
          "                SurfaceRecon rec;\n"
          "                bool rec_ok = false;\n"
          "                if (recon_on)\n")
    M.cut("Sp site: Tref's LATERAL-face note",
          "//  * LATERAL face (reached only under spall.per_cell_removal,",
          "//    one model (this one is isothermal across the flake).")
    M.rep("Sp site: Tref, the top face only",
          "                    if (!cref_on_sp) return T0_at(ii, jj, kk);\n"
          "                    if (nx_in == 0 && ny_in == 0)\n"
          "                    {\n"
          "                        const std::size_t c =\n"
          "                            static_cast<std::size_t>(col_index_sp(i, j));\n"
          "                        return CrefTAtDepth(c, (depth_n + cref_d_sp) * icth\n"
          "                                               + Set::Scalar(cref_ktop[c] - ksurf) * dz,\n"
          "                                            dz, sp_cref_vert, sp_cref_clamp);\n"
          "                    }\n"
          "                    ++sp_cref_lat;\n"
          "                    const Set::Scalar T0c = T0_at(i, j, k);\n"
          "                    const Set::Scalar dTc = T_a(i, j, k) - T0c;\n"
          "                    const Set::Scalar dTb = T_a(i + nx_in, j + ny_in, k) - T0c;\n"
          "                    if (!(dTc > 0.0) || !(dTb > 0.0) || dTc <= dTb)\n"
          "                        return T_a(i, j, k);\n"
          "                    const Set::Scalar ell_s = DX[0] / std::log(dTc / dTb);\n"
          "                    if (!(ell_s > 0.0) || !std::isfinite(ell_s))\n"
          "                        return T_a(i, j, k);\n"
          "                    return T0c + dTc * std::exp(-(depth_n * icth + cref_d_sp) / ell_s);\n",
          "                    if (!cref_on_sp) return T0_at(ii, jj, kk);\n"
          "                    const std::size_t c =\n"
          "                        static_cast<std::size_t>(col_index_sp(i, j));\n"
          "                    return CrefTAtDepth(c, (depth_n + cref_d_sp) * icth\n"
          "                                           + Set::Scalar(cref_ktop[c] - ksurf) * dz,\n"
          "                                        dz, sp_cref_vert, sp_cref_clamp);\n")
    M.rep("Sp site: the sample steps down the column only",
          "                        // containing the requested depth along the INWARD FACE\n"
          "                        // NORMAL, clamped into the box.\n"
          "                        //\n"
          "                        // D2q-2b: for the top face (nz_in = -1, nx_in = ny_in\n"
          "                        // = 0) this is bit-for-bit the pre-D2q-2b expression.\n"
          "                        // For a lateral face it steps sideways into the rock\n"
          "                        // instead of downwards -- the same rule, evaluated\n"
          "                        // where it currently is not. At a0 = 20 um against\n"
          "                        // dx = 2 mm the step count is 0 for every sample, so\n"
          "                        // the criterion is a function of THIS CELL alone.\n"
          "                        const int nstep = static_cast<int>(std::floor(dv / dz));\n"
          "                        int i_d = i + nx_in * nstep;\n"
          "                        int j_d = j + ny_in * nstep;\n"
          "                        int k_d = k + nz_in * nstep;\n"
          "                        if (k_d > ksurf)    k_d = ksurf;\n"
          "                        if (k_d < vbx_lo_k) k_d = vbx_lo_k;\n"
          "                        if (i_d < dlo_pc[0]) i_d = dlo_pc[0];\n"
          "                        if (i_d > dhi_pc[0]) i_d = dhi_pc[0];\n"
          "                        if (j_d < dlo_pc[1]) j_d = dlo_pc[1];\n"
          "                        if (j_d > dhi_pc[1]) j_d = dhi_pc[1];\n",
          "                        // containing the requested depth below the top face,\n"
          "                        // clamped into the box.\n"
          "                        const int nstep = static_cast<int>(std::floor(dv / dz));\n"
          "                        const int i_d = i, j_d = j;\n"
          "                        int k_d = k - nstep;\n"
          "                        if (k_d > ksurf)    k_d = ksurf;\n"
          "                        if (k_d < vbx_lo_k) k_d = vbx_lo_k;\n")
    M.rep("Sp site: the reference note",
          "                        // cell at d = 0 (bit-identical); under d > 0 the\n"
          "                        // gathered column (top face) or the profile model\n"
          "                        // (lateral face). ⚠ The SAMPLE T_a(i_d, j_d, k_d) is\n",
          "                        // cell at d = 0 (bit-identical); under d > 0 the\n"
          "                        // gathered column. ⚠ The SAMPLE T_a(i_d, j_d, k_d) is\n")
    M.rep("Sp site: cref reduction carries 0 for the lateral share",
          "        // Removal.H adds the depth scan's and the lateral criterion's share.\n"
          "        if (cref_on_sp)\n"
          "        {\n"
          "            long cnt[4] = {sp_cref_vert, sp_cref_lat, sp_cref_clamp, sp_cref_neg};\n",
          "        // Removal.H adds the depth scan's share. (P2-10: no lateral site is\n"
          "        // left, so cref_n_lat stays 0.)\n"
          "        if (cref_on_sp)\n"
          "        {\n"
          "            long cnt[4] = {sp_cref_vert, 0, sp_cref_clamp, sp_cref_neg};\n")
    M.cut("PerCellRemoval", "// D2q-2b items 6 + 8: LATERAL per-cell removal.",
          "pcr_overhang_cells = n_over;", end_off=2)
    M.cut("per-cell removal members", "// D2q-2b item 5: per-cell removal (default off). The two preconditions",
          "Set::Scalar pcr_lateral_cells  = 0.0;", end_off=1)
    M.cut("lateral key members", "// ⚠⚠ P1f: the LATERAL mode. Both default OFF.",
          "bool sp_lateral_thermal = false;  // spall.lateral_thermal (P1g-2 the sink)")
    M.cut("lateral diagnostic members", "// the LATERAL site's own pair -- reported separately, because the wall's",
          "Set::Scalar lat_Tface_max  = 0.0;  // T_face of the candidate attaining it [K]")
    M.rep("ReconFromX note", "    // Everything from a given x. Split out from ReconstructSurfaceT so the\n",
          "    // Everything from a given x. Split out of the reconstruction so the\n")
    M.cut("ReconstructSurfaceT (0 callers)",
          "static SurfaceRecon ReconstructSurfaceT(Set::Scalar T_top, Set::Scalar T_1,",
          "RawSurfaceX(T_top, T_1, T0), T_gas_cap, c_f, h_f_max);", end_off=2)


def main():
    stage, src, out = sys.argv[1:4]
    log = []
    D = load(src, log)
    {"c1": c1, "c2": c2}[stage](D)
    save(D, out, log)
    print("\n".join(log))


# ======================================================================= C2 ==
RETIRED_C2 = RETIRED_C1 + [
    "surface_patch.side_face_flux", "surface_patch.side_face_h", "surface_patch.side_face_T",
    "surface_patch.side_face_h_r_min", "surface_patch.side_face_h_law",
    "surface_patch.side_face_h_datum", "surface_patch.side_face_factor_above",
    "surface_patch.side_face_area_law", "surface_patch.side_face_step_height",
    "surface_patch.side_face_staircase", "surface_patch.step_diag_radius",
    "surface_patch.jet_wall_march", "surface_patch.jet_wall_mdot", "surface_patch.jet_body",
    "surface_patch.jet_body_T", "surface_patch.jet_body_h_scale"]

PAD_PARSE = (
    "            // D2aa deliverable 1: the pad ring's OPEN FRACTION. P2-10: read here,\n"
    "            // outside the deleted side-face block that held it, so a face-list deck\n"
    "            // keeps the key (FaceH and BuildFlameColumns read it).\n"
    "            // ⚠ contains() BEFORE query_default, which ADDS the key.\n"
    "            {\n"
    "                const bool pad_given = pp.contains(\"pad_shield_open\");\n"
    "                pp.query_default(\"pad_shield_open\", v.pad_shield_open, -1.0);\n"
    "                if (pad_given\n"
    "                    && !(v.pad_shield_open >= 0.0 && v.pad_shield_open <= 1.0))\n"
    "                    Util::Abort(INFO, \"surface_patch.pad_shield_open is the OPEN FRACTION \"\n"
    "                                      \"of the pad ring and must lie in [0, 1]; got \",\n"
    "                                v.pad_shield_open);\n"
    "            }")


def c2(D):
    R, M, F = D["R"], D["M"], D["F"]
    M.rep("the retired-key list (C2: + the side-face and wall-march keys)",
          retired_block(RETIRED_C1), retired_block(RETIRED_C2))
    M.cut("<cstdlib> (std::strtod, side_face_T = live only)", "#include <cstdlib>   // D2ab: std::strtod")
    # ---- registrations
    M.cut("surf_side_P_{top,buried,rin,rout} registrations",
          "value.RegisterIntegratedVariable(&value.surf_side_P_top,", "\"surf_side_P_rout\", false);")
    M.cut("wall-march registrations", "// P2-3: gated on the march's own key and appended after every",
          "\"wall_L_fit\", false);", end_off=1)
    M.cut("side_face_flux registrations (group B) + jet_body + jet_wall_mdot",
          "// D2k: side-face diagnostics, registered LAST and only when the key",
          "\"jet_P_retained\", false);", end_off=1)
    M.cut("side_face_step_height material abort",
          "// D2z: alpha = kappa/(rho*cp) must come FROM THE RUN. Without a",
          "\"this closure exists to remove\");", end_off=1)
    M.cut("face-list abort for jet_wall_march + the group-B message",
          "            if (sp.jet_wall_march)",
          "\"set but NOT READ under the face list:\", s);", end_off=1)
    for old, new in [
            ('fr_abort("regime 3 (bit.enabled)", "P2-10");',
             'fr_abort("regime 3 (bit.enabled)", "a later step");'),
            ('fr_abort("regime 4 (surface_patch.foot_body_clearance)", "P2-10");',
             'fr_abort("regime 4 (surface_patch.foot_body_clearance)", "P2-11");'),
            ('"it returns)", "P2-10");', '"it returns)", "a later step");'),
            ('fr_abort("spall.surface_normal = 1", "P2-10");',
             'fr_abort("spall.surface_normal = 1", "a later step");'),
            ('fr_abort("amr.max_level > 0", "P2-10");', 'fr_abort("amr.max_level > 0", "a later step");'),
            ('"reconstruction)", "P2-10");', '"reconstruction)", "a later step");'),
            ('(R3 IS the remap)", "P2-10");', '(R3 IS the remap)", "a later step");'),
            ('"weibull.size_volume = flux (per-column flux lengths)", "P2-10");',
             '"weibull.size_volume = flux (per-column flux lengths)", "a later step");'),
            ('"spall.flake_coherence_length > 0 (a per-column cap)", "P2-10");',
             '"spall.flake_coherence_length > 0 (a per-column cap)", "a later step");'),
            ('"(a top-plane labeller)", "P2-10");', '"(a top-plane labeller)", "a later step");')]:
        M.rep("face-list abort owner: " + new.split('"')[1][:40], old, new)
    M.rep("the unread-keys line (no group B, no lateral keys left)",
          "            {\n"
          "                // unread under the key; no abort, so a control differs by\n"
          "                // this key only\n"
          "                std::string s;\n"
          "                for (const auto& k : sp.groupB_given) s += \" surface_patch.\" + k;\n"
          "                if (!s.empty())\n"
          "                    Util::Message(INFO, \"surface_patch.face_list = 1 (P2-9): these keys are \"\n"
          "                                  \"set but NOT READ under the line rule:\", s);\n"
          "            }\n", "")
    M.cut("UpdateWallLawState call", "// P2-1 deliverable 2: the wall law's per-step state and the ONE height",
          "UpdateWallLawState(lev);")
    # ---- the conduction kernel
    M.cut("kernel: the step-height column phase",
          "// D2z: advance the per-column recession phase off this same mask,",
          "const int*         c_nv = col_nvoid.empty()  ? nullptr : col_nvoid.data();")
    M.cut("kernel: side-face locals",
          "// D2k: side-face flux. One SurfaceCellFlux call serves both lateral",
          "* surface_patch.side_face_h_r_min;")
    M.rep("kernel: inv_dx only (inv_dy and the cubic-cell abort were side-only)",
          "        const Set::Scalar inv_dx  = 1.0 / DX[0];\n"
          "        const Set::Scalar inv_dy  = 1.0 / DX[1];\n"
          "        if (side_on && !(std::abs(DX[0] - DX[1]) <= 1.0e-12 * DX[0] &&\n"
          "                         std::abs(DX[0] - DX[2]) <= 1.0e-12 * DX[0]))\n"
          "            Util::Abort(INFO, \"surface_patch.side_face_flux = 1 requires cubic cells \"\n"
          "                        \"(dx = dy = dz); got \", DX[0], \", \", DX[1], \", \", DX[2]);\n",
          "        const Set::Scalar inv_dx  = 1.0 / DX[0];\n")
    M.cut("kernel: side_T_ex", "// D2ab deliverable 1: the RESOLVED wall temperature, not",
          "const Set::Scalar side_T_ex = side_T_eff;")
    M.cut("kernel: flame_s pointer (side factor only)",
          "// D2r-1: per-column stand-off, so the side factor can differ above and",
          "const Set::Scalar* fs  = sp_on ? flame_s.data()  : nullptr;")
    M.cut("kernel: sp_surfform (the lateral rule's form test)",
          "// P0b rev5 deliverable 7: the lateral rule fires only under `surface`,",
          "SurfacePatch::RobinForm::Surface;")
    M.cut("kernel: the step conductance in k_face",
          "// D2z Step 3: correct the LATERAL face conductance across a",
          "kf *= StepConductance(d, kT - kS, a_diff);", end_off=2)
    M.cut("kernel: the side-face flux (D2k ... P2-3)",
          "// D2k: flux on EXPOSED VERTICAL FACES of a staircase wall.",
          "* (Set::Scalar(nfx) * inv_dx + Set::Scalar(nfy) * inv_dy);", end_off=2)
    # ---- BuildFlameColumns
    M.cut("BuildFlameColumns: side_T_eff",
          "// D2ab deliverable 1: resolve the wall gas temperature ONCE per step,",
          "side_T_eff = surface_patch.side_face_T;")
    M.cut("BuildFlameColumns: WallStepDiagnostic call",
          "// D2r-1 item 3: gated on the same key as its thermo columns, so an",
          "WallStepDiagnostic(lev, k_top_col, time);")
    M.cut("UpdateColumnPhase + AssertStepSplit + WallStepDiagnostic",
          "// D2z: per-column recession phase, advanced once per step from the",
          "// D2b: burner-on-feet descent (surface_patch.nozzle_descent = feet), lev 0,",
          keep_end=True)
    # ---- JetEnthalpyMarch
    M.cut("march: side locals", "// D2ab deliverable 2: the annulus pair's radius gate (which faces take",
          "(surface_patch.side_face_area_law == SurfacePatch::SideAreaLaw::Wetted);")
    M.cut("march: surf_form_m", "// P0b rev5 deliverable 7, march copy of the lateral rule (rev7 d6 switch).",
          "SurfacePatch::RobinForm::Surface;")
    M.cut("march: the side slot gather (note)",
          "// D2k: exposed vertical faces, so the side power comes OUT OF THE JET",
          "// below), so the slot gather is not built.")
    M.cut("march: the side slot gather",
          "const bool side_on = surface_patch.side_face_flux && !fl_on_m;",
          "ReduceRealSum(sN.data(), static_cast<int>(side_total));", end_off=2)
    M.cut("march: q_col's side note + side_f_m",
          "// D2k: with side_face_flux the column draws on the jet through its",
          "const Set::Scalar side_f_m = Set::Scalar(1.0);")
    M.cut("march: q_col's side term",
          "// ⚠⚠ P2-3: `== 1` and not a truth test. The gather now also holds",
          "q_tot += side_f_m * a_rob * sN[idx] * w;", end_off=2)
    M.cut("the wall march + the cold body", "// P2-3 deliverable 2 -- THE WALL MARCH.",
          "jet_P_wall = P_wall;", start_off=-1, end_off=1)
    M.rep("march: the exhaust with no wall section",
          "        // ⚠⚠ `jet_T_exhaust` KEEPS ITS MEANING: the gas state leaving the FLOOR.\n"
          "        // With the wall march on, `T_gas` has since been marched up the bore, so\n"
          "        // the floor-exit value is recovered from the stream's own bookkeeping\n"
          "        // rather than by reordering the assignments -- P_wall is by construction\n"
          "        // mcp*(T_floor_exit - T_mouth).\n"
          "        // P2-4: the floor-exit value is recovered from BOTH wall-section sinks.\n"
          "        // ⚠ mcp_w here must be the same mass the bin loop marched, or the\n"
          "        // recovered floor exit is wrong under `fresh`.\n"
          "        const Set::Scalar mcp_w =\n"
          "            (surface_patch.jet_wall_mdot == SurfacePatch::WallMdot::Fresh)\n"
          "                ? surface_patch.jet_mdot * surface_patch.jet_cp : mcp;\n"
          "        jet_T_exhaust = (mcp_w > 0.0)\n"
          "            ? (T_gas + (P_wall + jet_P_body) / mcp_w) : T_gas;\n"
          "        // the stream that actually leaves the domain is at the MOUTH, so the\n"
          "        // exhaust term is taken there and `face + wall + exhaust == cap` holds.\n"
          "        // ⚠⚠ UNDER `fresh` ONLY THE FRESH MASS LEAVES BY THE MOUTH. The\n"
          "        // entrained share never climbed the bore, so it is booked as RETAINED\n"
          "        // and the identity becomes\n"
          "        //   face + wall + body + exhaust_mouth + retained == cap,\n"
          "        // which holds because\n"
          "        //   mcp_jet (T_fe - T_ent) = mcp_fresh (T_fe - T_ent) + retained\n"
          "        // and the fresh part then pays wall + body + exhaust_mouth.\n"
          "        jet_P_exhaust = mcp_w * (T_gas - T_ent);\n"
          "        jet_P_retained = (mcp_w == mcp)\n"
          "            ? Set::Scalar(0.0)\n"
          "            : (mcp - mcp_w) * (jet_T_exhaust - T_ent);\n",
          "        // the gas leaving the floor, which is the gas leaving the domain (P2-10:\n"
          "        // the wall march that carried it on up the bore is deleted).\n"
          "        jet_T_exhaust = T_gas;\n"
          "        jet_P_exhaust = mcp * (T_gas - T_ent);\n")
    M.rep("march: the budget identity without the wall sinks",
          "              P_face <= jet_P_cap + tol &&\n"
          "              // P2-3: the wall is a third sink on the SAME stream.\n"
          "              // P2-4: the body is a fourth, and under `fresh` the retained\n"
          "              // entrained share is a fifth term that never leaves.\n"
          "              P_wall >= -tol && jet_P_body >= -tol && jet_P_retained >= -tol &&\n"
          "              std::abs(P_face + P_wall + jet_P_body + jet_P_exhaust\n"
          "                       + jet_P_retained - jet_P_cap) <= tol &&\n",
          "              P_face <= jet_P_cap + tol &&\n"
          "              std::abs(P_face + jet_P_exhaust - jet_P_cap) <= tol &&\n")
    M.rep("march: the invariant's message",
          "                        \", jet_P_neg = \", P_neg, \", jet_P_wall = \", P_wall,\n"
          "                        \", jet_P_body = \", jet_P_body,\n"
          "                        \", jet_P_retained = \", jet_P_retained,\n",
          "                        \", jet_P_neg = \", P_neg,\n")
    # ---- EnergyLedger
    M.cut("ledger: sp_surfform_l", "// P0b rev5 deliverable 7, ledger copy of the lateral rule. ⚠ It must",
          "SurfacePatch::RobinForm::Surface;")
    M.rep("ledger: side_P_topcell / side_P_buried",
          "        Set::Scalar P_robin = 0.0, P_loss = 0.0, P_robin_side = 0.0;\n"
          "        Set::Scalar side_P_topcell = 0.0, side_P_buried = 0.0;\n",
          "        Set::Scalar P_robin = 0.0, P_loss = 0.0, P_robin_side = 0.0;\n")
    M.cut("ledger: side_P_rin / side_P_rout",
          "// ⚠⚠ REVISION 7 DELIVERABLE 2 -- RE-BIN AT THE LIVE GATE.",
          "Set::Scalar side_P_rin = 0.0, side_P_rout = 0.0;")
    M.rep("ledger: the stray `// D2k`",
          "        std::vector<Set::Scalar> sidecol;\n        // D2k\n",
          "        std::vector<Set::Scalar> sidecol;\n")
    M.cut("ledger: side locals",
          "const bool        side_on = sp_on && surface_patch.side_face_flux && !fl_on_l;",
          "// report. (P2-8 group A: the side_w_mean sums are deleted with it.)")
    M.cut("ledger: the side-face term (D2k ... P2-3)",
          "// D2k: the identical side-face term, or ledger_err breaks.",
          "side_ab_n += Set::Scalar(nfx + nfy);", start_off=-1, end_off=3)
    M.cut("ledger: the side reductions",
          "amrex::ParallelDescriptor::ReduceRealSum(side_rob);",
          "side_cols_above   = side_ab_n;", start_off=-2, end_off=2)
    M.rep("ledger: surf_side_P reductions",
          "            amrex::ParallelDescriptor::ReduceRealSum(P_robin_side);\n"
          "            amrex::ParallelDescriptor::ReduceRealSum(side_P_topcell);\n"
          "            amrex::ParallelDescriptor::ReduceRealSum(side_P_buried);\n"
          "            amrex::ParallelDescriptor::ReduceRealSum(side_P_rin);\n"
          "            amrex::ParallelDescriptor::ReduceRealSum(side_P_rout);\n",
          "            amrex::ParallelDescriptor::ReduceRealSum(P_robin_side);\n")
    M.rep("ledger: surf_side_P assignments",
          "            surf_side_P_top    = side_P_topcell;\n"
          "            surf_side_P_buried = side_P_buried;\n"
          "            surf_side_P_rin    = side_P_rin;\n"
          "            surf_side_P_rout   = side_P_rout;\n", "")
    # ---- the wall law (P2-1/P2-3)
    M.cut("the composite wall law: constants, resolvers, datum, UpdateWallLawState",
          "// P2-1: THE COMPOSITE WALL LAW", "if (!wall_law_ready) wall_corner.clear();",
          start_off=-1, end_off=2)
    # ---- members
    M.cut("side-face members + the step helpers (ColPhi ... StepConductance)",
          "// D2k: side-face flux on exposed vertical faces (surface_patch.side_face_flux).",
          "Set::Scalar side_faces  = 0.0;  // thermo: number of exposed lateral faces heated", end_off=1)
    M.cut("wall-law, wall-march and cold-body members",
          "// ---- P2-1: THE DERIVED WALL COEFFICIENT ----",
          "Set::Scalar wall_law_on     = 0.0;  // 1 while the law is actually applied")
    M.cut("side_P_col (no reader)",
          "// Per-column side power at the current flame vectors, for the enthalpy",
          "std::vector<Set::Scalar> side_P_col;")
    M.cut("surf_side_P members", "// P0b rev5 deliverable 9 bin 2: the side Robin power split by whether the",
          "Set::Scalar surf_side_P_rout   = 0.0;   // r >= side_face_h_r_min")
    # ---- SurfacePatch
    M.cut("SurfacePatch: side_face_flux", "// D2k: flux on exposed vertical faces of a staircase wall (default-off;",
          "bool side_face_flux = false;")
    M.cut("SurfacePatch: groupB_given", "// P2-8: the GROUP-B keys the deck set, captured before any",
          "std::vector<std::string> groupB_given;")
    M.cut("SurfacePatch: group B, wall march, cold body (up to pad_shield_open)",
          "// P2-8 group A: `side_face_factor` (1.0) and `side_face_h_max` (80,",
          "bool side_face_h_derived = false;")
    M.cut("SurfacePatch: step height, staircase, area law",
          "// P2-8 group A: `side_face_area_exact` (and RiserWeight) DELETED.",
          "SideAreaLaw side_face_area_law = SideAreaLaw::Step;")
    M.cut("CheckPadRing: the side_face_T / h_law / wall march / body guards",
          "// D2ab deliverable 1. Same both-exits home, for the same reason:",
          "\"surface_patch.jet_closure = enthalpy.\");")
    M.cut("Parse: the group-B capture",
          "// P2-8: group-B keys present on the command line or in the deck,",
          "if (pp.contains(kB)) v.groupB_given.push_back(kB);")
    M.cut("Parse: the side_face_flux block -> the pad_shield_open parse alone",
          "// D2k: side-face flux on exposed vertical faces of a staircase",
          "else if (pp.contains(\"side_face_h_r_min\"))", end_off=3, insert=PAD_PARSE)
    # ---- comments that name deleted code (a9)
    R.rep("comment: ReconstructSurfaceT",
          "        // over the WRONG thickness, and ReconstructSurfaceT -- which assumes\n",
          "        // over the WRONG thickness, and the reconstruction -- which assumes\n")
    M.rep("comment: UpdateWallLawState",
          "    // Called FIRST in Advance, beside P2-1's UpdateWallLawState, so the\n",
          "    // Called FIRST in Advance (where P2-1's wall-law state also ran), so the\n")
    M.rep("comment: unit/lateral_drive",
          "        // exactly what `unit/lateral_drive/input:120` does, while Meier's\n",
          "        // exactly what P1f's lateral-drive fixture did, while Meier's\n")
    # ---- FaceRemoval.H: the wall share of the flake ledger was printed only under jet_body
    F.rep("E_wall", "Set::Scalar E_left_tot = 0.0, dH_book = 0.0, E_wall = 0.0;",
          "Set::Scalar E_left_tot = 0.0, dH_book = 0.0;")
    F.cut("E_wall accumulation", "if (cls != 0) E_wall += E_left;")
    F.cut("ledger_E_flake_wall (registered only under jet_body)", "ledger_E_flake_wall += E_wall * Vcell;")
    del R


if __name__ == "__main__":
    main()
