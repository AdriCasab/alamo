#!/usr/bin/env python3
"""P2-28 stage 2 -- THE COLUMN PATH DELETED (census H). A transform of the stage-2 PRE-IMAGES (src_pre_2/) into src/:
re-runnable, every anchor asserted unique, so the edit is one reviewable artefact rather than a hand session.

  MMWSpalling.H   site A + the cluster labeller, the column reconstruction (BuildSurfaceTColumns & co.), the pinned
                  closure, the column patch diagnostics, the column gathers, the AMR column repair and surface band, the
                  bit, the march's column flux, the P0 keys' flags, DetachMode, the column CSVs, the prescribed
                  verification hooks -- deleted; robin_form {cell, face, surface} (pinned aborts); one-removal-path
                  interlocks; ledger_err -> ledger_err_raw; RemovalStep in Advance
  FaceRemoval.H   RemovalStep (the step's three calls); SurfaceCellFlux's narrowed signature
  Shared.H        comments
  Removal.H       deleted
Then R11: every column-only name 0 times in src/ (comments scrubbed), the retired key strings once each.
Usage: stage2_src.py [--check]   (--check: transform in memory, report, write nothing)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PRE = os.path.join(HERE, "src_pre_2", "src", "Integrator")
SRC = os.path.join(ROOT, "src", "Integrator")
LOG = []


def sub(s, old, new, count=1):
    n = s.count(old)
    assert n == count, "expected %d, found %d: %r" % (count, n, old[:120])
    return s.replace(old, new)


def between(s, start, end, new):
    """replace [start .. end] (both markers included; each unique)"""
    assert s.count(start) == 1, "start x%d: %r" % (s.count(start), start[:100])
    i = s.index(start)
    j = s.index(end, i)
    return s[:i] + new + s[j + len(end):]


def del_func(s, sig):
    """delete the class-scope function whose signature line contains `sig` (unique), with the comment block
    directly above it and one trailing blank line"""
    L = s.split("\n")
    idx = [n for n, l in enumerate(L) if sig in l]
    assert len(idx) == 1, "sig x%d: %r" % (len(idx), sig)
    i = idx[0]
    j = i
    while j > 0 and L[j - 1].strip().startswith("//"):
        j -= 1
    k = i
    while L[k] != "    }":
        k += 1
    end = k + 1
    if end < len(L) and L[end].strip() == "":
        end += 1
    LOG.append("deleted %-34s %5d lines (from line %d)" % (sig.split("(")[0].split()[-1], end - j, j + 1))
    return "\n".join(L[:j] + L[end:])


R11_NAMES = ["UpdateRemovalAfterCohesive", "BuildSurfaceTColumns", "UpdateSpAfterMechanics", "UpdateSpClusters",
             "PatchDiagnostics", "GatherColumnSurface", "ColumnSolidRunBelowTop", "PinRule", "BuildPinEffective",
             "SurfaceRobinT", "SharedSurfaceRecon", "RepairRemovalStateAfterRegrid", "CheckSurfaceCovered",
             "rc_ell_prev", "DetachMode", "connected_cluster", "E_sup_local", "Removal.H"]
SCRUB = {"UpdateRemovalAfterCohesive": "the column removal step (deleted, P2-28)",
         "BuildSurfaceTColumns": "the column reconstruction (deleted, P2-28)",
         "SharedSurfaceRecon": "the column reconstruction (deleted, P2-28)",
         "SurfaceRobinT": "the column reconstruction (deleted, P2-28)",
         "UpdateSpAfterMechanics": "site A (deleted, P2-28)",
         "UpdateSpClusters": "the cluster labeller (deleted, P2-28)",
         "PatchDiagnostics": "the column patch diagnostics (deleted, P2-28)",
         "GatherColumnSurface": "a column gather (deleted, P2-28)",
         "ColumnSolidRunBelowTop": "a column gather (deleted, P2-28)",
         "BuildPinEffective": "the pin rule (deleted, P2-28)",
         "PinRule": "the pin rule (deleted, P2-28)",
         "RepairRemovalStateAfterRegrid": "the AMR column repair (deleted, P2-28)",
         "CheckSurfaceCovered": "the AMR column repair (deleted, P2-28)",
         "rc_ell_prev": "a column diagnostic (deleted, P2-28)",
         "DetachMode": "the cluster detach mode (deleted, P2-28)",
         "connected_cluster": "the cluster detach mode (deleted, P2-28)",
         "Removal.H": "the column path (deleted, P2-28)"}


def scrub(s):
    """the R11 names in COMMENTS -> a neutral phrase (code references are deleted above, so any left in code
    is a compile error, not a scrub)"""
    out = []
    for ln in s.split("\n"):
        if "//" in ln:
            a, b = ln.split("//", 1)
            for k, v in SCRUB.items():
                b = re.sub(r"(?<![A-Za-z_])" + re.escape(k), v, b)
            ln = a + "//" + b
        out.append(ln)
    return "\n".join(out)


# =========================================================================================== MMWSpalling.H ====
def mmw(s):
    # ---- the deleted functions ----------------------------------------------------------------------------------
    for sig in ["    void BuildSurfaceTColumns(int lev)", "    bool SharedSurfaceRecon(std::size_t c, SurfaceRecon& r) const",
                "    Set::Scalar SurfaceRobinT(std::size_t c) const", "    void BuildPinEffective(int lev, Set::Scalar time)",
                "    Set::Scalar PinRule(std::size_t c, Set::Scalar h, Set::Scalar T_gas,",
                "    void PatchDiagnostics(int lev, const amrex::Vector<int>& k_top_col)",
                "    void RepairRemovalStateAfterRegrid(int lev)", "    void CheckSurfaceCovered(int lev, Set::Scalar time)",
                "    amrex::Vector<int> ColumnSolidRunBelowTop(int lev, const amrex::Vector<int>& k_top)",
                "    static void ParseBitSettings(MMWSpalling& value, IO::ParmParse& pp)",
                "    void GatherColumnSurface(int lev, amrex::Vector<int>& k_top_col,",
                "    void UpdateSpAfterMechanics(int lev, Set::Scalar time, Set::Scalar dt)",
                "    void UpdateSpClusters(int lev, Set::Scalar time, Set::Scalar dt)"]:
        s = del_func(s, sig)

    s = sub(s, """// P1g-2c deliverable 1: the lateral thermal shift's write-back collects
// per-cell deltas in a small ordered map before touching the field, so a cell
// written by two rows is composed instead of overwritten. Removal.H is included
// INSIDE the class body, so this cannot live there.""", """// The face list's ONE WRITE collects per-cell deltas in a small ordered map before
// touching the field (FaceRemoval.H, included INSIDE the class body, so this
// cannot live there).""")

    # ---- Parse: the P2-28 abort lists and retired-at-default rows ------------------------------------------------
    s = sub(s, """                    Util::Abort(INFO, k, " was retired in P2-27 (the Drucker-Prager damage "
                                "law and the grain-boundary cohesive zone are deleted: LEFM "
                                "is the model); remove it from the input");
""", """                    Util::Abort(INFO, k, " was retired in P2-27 (the Drucker-Prager damage "
                                "law and the grain-boundary cohesive zone are deleted: LEFM "
                                "is the model); remove it from the input");
        // P2-28 (census H): THE COLUMN REMOVAL PATH is retired -- one removal path, the face list.
        // Its own keys abort by PRESENCE, whatever the value: the bit (regime 3), the column CSVs,
        // the cluster labeller's area, the AMR surface band and site A's verification hooks.
        for (const char* k : {"spall.clusters_csv", "spall.h_col_events_csv", "weibull.A_crit",
                              "surface.refine_depth", "surface.refine_height"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-28 (the column removal path is deleted: "
                            "the face list is the one removal path); remove it from the input");
        for (const char* family : {"bit", "spallation.sp.prescribed_stress",
                                   "spallation.sp.prescribed_T_surface"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                Util::Abort(INFO, k, " was retired in P2-28 (the column removal path is deleted: "
                            "the face list is the one removal path); remove it from the input");
""")
    s = sub(s, """            {"spall.sample",                    "top_cell",
             "P2-27: the damage-law removal's D-sampling window is deleted with it"},
        };""", """            {"spall.sample",                    "top_cell",
             "P2-27: the damage-law removal's D-sampling window is deleted with it"},
            // P2-28 (census H): the column path's selectors, each retired at the ONE value the face list runs
            {"spall.surface_T",      "reconstructed",
             "P2-28: the face list's R1 IS the reconstruction (site A's cell_average / diagnostic forms are deleted)"},
            {"spall.removal_energy", "flake",
             "P2-28: the face list's R3 IS the flake remap (the column path's whole-cell removal is deleted)"},
            {"weibull.detach_mode",  "per_face",
             "P2-28: the other mode was the cluster labeller, deleted with site A"},
            {"surface.follow_mask",  "1",
             "P2-28: the face list's `top` class IS the mask column top, and removal requires the list"},
        };""")
    s = sub(s, """        // Step 19: optional mechanical bit feed (drawworks). Default-off.
        // Parsed after spall/vapor so it can require one of them.
        ParseBitSettings(value, pp);

""", "")

    # ---- registrations ------------------------------------------------------------------------------------------
    s = between(s, "            // Step 12: surface Sp diagnostic for the sp_weibull branch only.\n",
                """                    value.RegisterIntegratedVariable(&value.recon_dT_mean,
                                                     "recon_dT_mean", false);
                }
            }
""", """            // Step 15b: the per-cell Weibull flaw length (the face list's R2 reads it at each line's
            // surface cell). P2-28: site A's fields (Sp_field, Sp_cluster_id, h_spall_field,
            // crack_event_count, T_surf_field, ell_eff_field) and the recon_* columns are deleted with it.
            if (value.spallation_model == SpallationModel::SpWeibull && value.weibull_enabled)
                value.RegisterNewFab(value.flaw_a_mf, nullptr, 1, 0, "flaw_a", true);
""")
    s = between(s, "            int fm = 0, le = 0;\n",
                """            pp.query_default("surface.refine_height", value.surface_refine_height, 0.0);
""", """            int le = 0;
            pp.query_default("energy_ledger.enabled", le, 0);
            // P2-28: surface.follow_mask is retired at 1 (the table above): the mask column top is the
            // face list's `top` class, and removal requires the list, so the switch IS the list's.
            value.surface_follow_mask = value.surface_patch.face_list;
            value.ledger_enabled      = (le != 0);
""")
    s = sub(s, """            // C1: pinned-surface closure (the pin comes from the removal
            // pipeline's per-column firing T at the mask top).
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Pinned)
            {
                value.RegisterIntegratedVariable(&value.pinned_cols_pinned, "pinned_cols_pinned", false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face,   "pinned_cols_face",   false);
                value.RegisterIntegratedVariable(&value.pinned_q_mean,      "pinned_q_mean",      false);
            }
""", "")
    s = sub(s, """                // Deliverable 1.4: pinned_idle_cycles governs the pin rule,
                // which this form replaces. Say so rather than ignore it.
                if (value.surface_patch.pinned_idle_cycles != 2.0)
                    Util::Message(INFO, "surface_patch.pinned_idle_cycles = ",
                                  value.surface_patch.pinned_idle_cycles,
                                  " is INERT under robin_form = surface: the pin rule it "
                                  "governs is replaced by the per-column reconstruction.");
                value.RegisterIntegratedVariable(&value.surf_cols_robin, "surf_cols_robin", false);
                value.RegisterIntegratedVariable(&value.surf_q_mean,     "surf_q_mean",     false);
                value.RegisterIntegratedVariable(&value.surf_side_share, "surf_side_share", false);""",
            """                // P2-28: surf_cols_robin / surf_q_mean (the column reconstruction's witnesses) are deleted
                value.RegisterIntegratedVariable(&value.surf_side_share, "surf_side_share", false);""")
    s = sub(s, """                value.RegisterIntegratedVariable(&value.ledger_err,     "ledger_err",     false);""",
            """                // P2-28 (P2-27 DR-1): RENAMED -- the all-cell residual omits the flake credit by
                // design, so it is NOT a closure on a spalling deck; ledger_err_solid closes
                value.RegisterIntegratedVariable(&value.ledger_err_raw, "ledger_err_raw", false);""")
    s = sub(s, """                // header alone whatever the numbers do.
                if (value.sp_surface_T_recon)
                {
                    value.RegisterIntegratedVariable(&value.ledger_E_flake,""",
            """                // header alone whatever the numbers do. P2-28: the reconstruction's one home is the face list.
                if (value.surface_patch.face_list)
                {
                    value.RegisterIntegratedVariable(&value.ledger_E_flake,""")
    s = sub(s, """            // D2a: pinned branch-reason counters (sum = pinned_cols_face) and the
            // spatial-flame diagnostics. Registered after every pre-D2a variable
            // so existing thermo.dat columns keep their positions.
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Pinned)
            {
                value.RegisterIntegratedVariable(&value.pinned_cols_face_unset,   "pinned_cols_face_unset",   false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face_idle,    "pinned_cols_face_idle",    false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face_qneg,    "pinned_cols_face_qneg",    false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face_minrule, "pinned_cols_face_minrule", false);
            }
""", """            // D2a: the spatial-flame diagnostic, after every pre-D2a variable (P2-28: the pinned
            // branch counters and patch_P_robin -- the column patch diagnostics -- are deleted).
""")
    s = sub(s, """                value.RegisterIntegratedVariable(&value.patch_min_standoff, "patch_min_standoff", false);
                value.RegisterIntegratedVariable(&value.patch_P_robin,      "patch_P_robin",      false);""",
            """                value.RegisterIntegratedVariable(&value.patch_min_standoff, "patch_min_standoff", false);""")

    # ---- CheckInterlocks ----------------------------------------------------------------------------------------
    s = between(s, "        // flake recession recedes phi by h_f = min(ell_eff, dz). ell_eff has\n",
                """                              "(the bit reuses the phi/removed surface machinery)");
""", "")
    s = sub(s, """                              "spallation.model = none nothing would fire");
""", """                              "spallation.model = none nothing would fire");
        // P2-28 (census H): ONE REMOVAL PATH. The column path is deleted, so removal runs on the face
        // list or not at all (the list's evaluate-only mode, spall and vapour off, is the onset decks').
        if ((value.spall_enabled || value.vapor_enabled) && !value.surface_patch.face_list)
            Util::Abort(INFO, "spall.enabled / vapor.enabled need surface_patch.face_list = 1: the column "
                              "removal path was retired in P2-28 (one removal path, the face list)");
        // P2-28: the jet march bills the face list's faces; its column flux is deleted
        if (value.surface_patch.enabled && value.surface_patch.jet_enthalpy &&
            !value.surface_patch.face_list)
            Util::Abort(INFO, "surface_patch.jet_closure = enthalpy needs surface_patch.face_list = 1: the "
                              "march's column flux was retired in P2-28");
""")
    s = between(s, "        // S1: surface-resolution study switches (microstructure path, single level).\n",
                """                Util::Abort(INFO, "surface.follow_mask = 1 is not supported with "
                                  "surface_patch.mode = prescribed_T");
        }
""", "")
    s = sub(s, """        if ((value.surface_follow_mask || value.ledger_enabled) &&
            !value.microstructure_enabled)
            Util::Abort(INFO, "surface.follow_mask / energy_ledger.enabled require "
                              "the microstructure path (microstructure.nphases)");""",
            """        if (value.ledger_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "energy_ledger.enabled requires "
                              "the microstructure path (microstructure.nphases)");""")
    s = between(s, "        // AMR: the surface band must live on the finest level.\n",
                """                                      "once per step on the finest level)");
        }
""", "")
    s = between(s, "        // C1: pinned-surface closure prerequisites (the pin comes from the\n",
                """                                  " < ", value.eps_rad, "): the face-temperature Newton assumes a "
                                  "monotone residual with a unique root");
        }
        // P0b deliverable 1""", """        // P0b deliverable 1""")
    s = sub(s, """            if (!value.surface_follow_mask)
                Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                  "surface.follow_mask = 1");""",
            """            // P2-28: the form's surface temperature is the face list's (R1): the column
            // reconstruction it used without the list is deleted
            if (!value.surface_patch.face_list)
                Util::Abort(INFO, "surface_patch.robin_form = surface requires surface_patch.face_list = 1: "
                                  "the column reconstruction it read without the list was retired in P2-28");""")
    s = sub(s, """            if (!value.sp_surface_T_apply)
                Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                  "spall.surface_T = reconstructed: the form evaluates "
                                  "the flux AND the losses at the reconstructed surface "
                                  "temperature, and without the reconstruction there is "
                                  "no surface temperature to evaluate them at.");
""", "")
    s = between(s, "        // P0 SCOPE, part 2: spallation.model is known here (part 1 is above).\n",
                """            Util::Abort(INFO, "spall.removal_energy = flake requires "
                              "spall.enabled = 1");
""", "")
    s = sub(s, """        if (value.surface_patch.enabled && value.surface_patch.use_expr)
        {
            if (!value.surface_follow_mask)
                Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr require "
                                  "surface.follow_mask = 1");
        }""", """        // P2-28: the flame map is evaluated per FACE (F2): the face list
        if (value.surface_patch.enabled && value.surface_patch.use_expr && !value.surface_patch.face_list)
            Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr need surface_patch.face_list = 1 (F2 "
                              "evaluates the map at every face centre; the column path was retired in P2-28)");""")
    s = sub(s, """            if (!value.sp_surface_T_apply)
                Util::Abort(INFO, "spall.buckle_span > 0 (P2-15) reads the reconstructed "
                            "heated-layer length: it needs spall.surface_T = reconstructed");
""", "")
    s = sub(s, """            if (!value.surface_follow_mask)
                Util::Abort(INFO, "surface_patch.face_list = 1 requires "
                            "surface.follow_mask = 1: the `top` class is the mask "
                            "column top");
""", "")
    s = sub(s, """            if (value.bit_enabled)   fr_abort("regime 3 (bit.enabled)", "a later step");
""", "")
    s = sub(s, """            if (!value.sp_surface_T_apply)
                fr_abort("spall.surface_T other than reconstructed (R1 IS the "
                         "reconstruction)", "a later step");
            if (!value.sp_removal_energy_flake)
                fr_abort("spall.removal_energy other than flake (R3 IS the remap)", "a later step");
            if (value.weibull_enabled &&
                value.weibull_detach_mode == DetachMode::ConnectedCluster &&
                value.weibull_A_crit > 0.0)
                fr_abort("weibull.detach_mode = connected_cluster with A_crit > 0 "
                         "(a top-plane labeller)", "a later step");
""", "")

    # ---- Regrid, Advance ----------------------------------------------------------------------------------------
    s = sub(s, """        // AMR: a remade level keeps exact phi/removed where it already
        // existed and gets conservatively interpolated values where it is
        // new; the lateral slope limiter smears the per-column height at
        // steps. Rebuild the column-linear phi from each column's top cell.
        if (removal_enabled() && lev > 0) RepairRemovalStateAfterRegrid(lev);
""", "")
    s = sub(s, """        // P0b deliverable 1.1: ONE surface temperature per column for this
        // whole step -- criterion, BC and removal. Must be first: every
        // consumer below reads it, and temp_mf is swapped into temp_old_mf
        // partway through.
        BuildSurfaceTColumns(lev);
""", "")
    s = between(s, "        // Step 12: Sp = K_I/K_Ic(T) surface diagnostic for the sp_weibull\n",
                "        UpdateRemovalAfterCohesive(lev, time, dt);\n",
                """        // Step 10/11, P2-28: THE removal step -- the face list's (FaceRemoval.H), the only one: the
        // criterion on every listed line and the recession, BEFORE the H/Temp swaps below, so what it
        // writes is the "old" state the thermal advance consumes.
        RemovalStep(lev, time, dt);
""")

    # ---- the kernel ---------------------------------------------------------------------------------------------
    s = sub(s, """        // every step. D is refilled every step (removal and damage write it).""",
            """        // every step. (D, refilled every step until P2-27, is deleted with the damage lineage.)""")
    s = sub(s, """        // C1: pinned runs the face Newton (face_form = true) plus the pin.
        const bool        sp_pinned = sp_on && surface_patch.robin_form ==
                                      SurfacePatch::RobinForm::Pinned;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face || sp_pinned;
        if (sp_pinned) BuildPinEffective(lev, time);
        const Set::Scalar* tpe = sp_pinned ? pin_T_eff.data() : nullptr;
""", """        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face;
""")
    s = sub(s, """                                    q_in, q_out, q_rob, q_los,
                                    tpe ? tpe[col_c] : Set::Scalar(0.0),
                                    nullptr,
                                    SurfaceRobinT(
                                        static_cast<std::size_t>(col_c)));""",
            """                                    q_in, q_out, q_rob, q_los);""")
    s = sub(s, """        if (ledger_enabled) EnergyLedger(lev, time, dt, k_top_col);
        // C1/D2a: pinned-closure branch counts and spatial-flame diagnostics.
        if (owner && (sp_pinned || (sp_on && surface_patch.use_expr)))
            PatchDiagnostics(lev, k_top_col);
""", """        if (ledger_enabled) EnergyLedger(lev, time, dt, k_top_col);
""")

    # ---- the march: the face list's faces only ------------------------------------------------------------------
    s = between(s, """        const bool pinned = surface_patch.robin_form == SurfacePatch::RobinForm::Pinned;
        const bool face   = surface_patch.robin_form == SurfacePatch::RobinForm::Face || pinned;
""", """        const bool fl_on_m = surface_patch.face_list && fl_ready;
""", """        // P2-28: the march bills the face list's faces (its column flux, and the top-cell gather it
        // read, are deleted with the column path; Parse requires the list under jet_closure = enthalpy)
        if (!(surface_patch.face_list && fl_ready))
            Util::Abort(INFO, "the jet march bills the face list's faces: the list is not built (P2-28)");
""")
    s = sub(s, """        const bool have_pin = pinned && pin_T.size() == ncols;
""", "")
    s = between(s, "        // q_robin of column c at gas temperature Tg (the kernel's inputs).\n",
                """                return q_tot;
            }
            return q_tot;
        };
""", """        // the robin power of column c's listed faces at gas temperature Tg, billed to its bin (P2-9 R1:
        // every listed face of the column -- the top included -- at its cell's T_s, FrFaceFlux, x w; the
        // face area dx*dz against the bin's dA = dx*dy, cubic)
        auto q_col = [&](std::size_t c, Set::Scalar Tg) -> Set::Scalar
        {
            Set::Scalar q_tot = 0.0;
            for (int q = fl_col_first[c]; q < fl_col_first[c] + fl_col_count[c]; ++q)
            {
                const FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                FrFaceFlux(e, Tg, e.K, e.T, DX[0], a_in, a_out, a_rob, a_los);
                q_tot += e.w * a_rob * (DX[2] / DX[0]);
            }
            return q_tot;
        };
""")
    s = sub(s, "        const bool path = surface_patch.face_heat_path && fl_on_m;",
            "        const bool path = surface_patch.face_heat_path;")
    s = sub(s, "                if (fl_on_m) fl_Pm += q_rob * dA;", "                fl_Pm += q_rob * dA;")
    s = sub(s, "        if (fl_on_m) fl_P_march = fl_Pm;", "        fl_P_march = fl_Pm;")

    # ---- SurfaceCellFlux: the pinned branch and its two parameters ----------------------------------------------
    s = between(s, "    // C1 pinned form (face_form = true and T_pin > 0, robin_form = pinned):\n",
                "    // so the cell and face arithmetic are untouched.\n", "")
    s = sub(s, """                                Set::Scalar& q_robin, Set::Scalar& q_loss,
                                Set::Scalar T_pin = 0.0, int* branch = nullptr,
                                Set::Scalar T_surf_in = 0.0,""", """                                Set::Scalar& q_robin, Set::Scalar& q_loss,
                                Set::Scalar T_surf_in = 0.0,""")
    s = sub(s, """            q_out   = 0.0;
            if (branch) *branch = 2;
            if (T_face_out) *T_face_out = Ts;   // P1f: the reconstructed surface T""",
            """            q_out   = 0.0;
            if (T_face_out) *T_face_out = Ts;   // P1f: the reconstructed surface T""")
    s = between(s, "            if (!(T_pin > 0.0 && Tf > T_pin))\n",
                """            if (T_face_out) *T_face_out = Ts;   // P1f: Ts == T_pin on this branch
            return;
        }
""", """            q_in    = G * (Tf - Tc);
            q_out   = 0.0;
            q_robin = h * (T_gas - Tf);
            q_loss  = q_robin - q_in;
            // P1f: the Newton's own root -- the temperature this q_robin used.
            // ⚠ Inferred on a LINEAR half-cell (G = 2k/dz) while the criterion
            // reads an exponential; the two disagree once ell_lat < dz/2, which
            // criterion (k2) scores rather than assumes away.
            if (T_face_out) *T_face_out = Tf;
            return;
        }
""")
    s = sub(s, """        // ⚠ Every existing caller passes nothing, so T_surf_in defaults to 0
        // and `cell`, `face` and `pinned` stay byte-identical. They are this
        // packet's CONTROLS and every frozen study from d2v to d2ad used
        // `pinned`, so that is not a convenience -- it is scored.""",
            """        // ⚠ A caller that passes nothing gets T_surf_in = 0, so `cell` and
        // `face` stay byte-identical (the face list passes its cell's T_s).""")

    # ---- the ledger ---------------------------------------------------------------------------------------------
    s = sub(s, """    //   ledger_err = (E_H - E_H0 - E_src)/max(|E_src|, 1e-30)""",
            """    //   ledger_err_raw = (E_H - E_H0 - E_src)/max(|E_src|, 1e-30)""")
    s = sub(s, """    // step. Likewise the C1 pin state (pin_T/pin_t/pin_rho/pin_Cp) is not
    // checkpointed: after a restart every column uses the face form until it
    // fires again. E_H is a compensated (Neumaier) sum so the""",
            """    // step. E_H is a compensated (Neumaier) sum so the""")
    s = sub(s, """        const bool        sp_pinned = sp_on && surface_patch.robin_form ==
                                      SurfacePatch::RobinForm::Pinned && lev == 0;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face || sp_pinned;""",
            """        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face;""")
    s = sub(s, """                                q_in, q_out, q_rob, q_los,
                                sp_pinned ? pin_T_eff[col] : Set::Scalar(0.0),
                                nullptr,
                                SurfaceRobinT(static_cast<std::size_t>(col)));""",
            """                                q_in, q_out, q_rob, q_los);""")
    s = sub(s, """            Set::Scalar nsurf = 0.0, qsum = 0.0;
            for (std::size_t c = 0; c < surf_T_col.size(); ++c)
                if (surf_T_col[c] > 0.0)
                { nsurf += 1.0; qsum += flame_h[c] * (flame_Tg[c] - surf_T_col[c]); }
            surf_cols_robin = nsurf;
            surf_q_mean     = (nsurf > 0.0) ? qsum / nsurf : 0.0;
""", "")
    s = sub(s, "        ledger_err = (dE - ledger_E_src) / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));",
            "        ledger_err_raw = (dE - ledger_E_src) / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));")
    s = sub(s, "        // gap between ledger_err and ledger_err_solid is the flake credit the",
            "        // gap between ledger_err_raw and ledger_err_solid is the flake credit the")
    s = sub(s, "    Set::Scalar ledger_err                   = 0.0;", "    Set::Scalar ledger_err_raw               = 0.0;")
    s = sub(s, """    // The per-class robin power, from the replicated list with the kernel's
    // own inputs (T_old, kc, flame_h/flame_Tg, SurfaceRobinT) -- the kernel's""",
            """    // The per-class robin power, from the replicated list with the kernel's
    // own inputs (T_old, kc, each face's h and gas temperature) -- the kernel's""")

    # ---- TagCellsForRefinement: the surface band ----------------------------------------------------------------
    s = between(s, "        // AMR surface band: every column's top solid cell (on this level's\n",
                """                if (k >= kt - n_below && k <= kt + n_above) tags(i, j, k) = amrex::TagBox::SET;
            });
        }
""", "        // P2-28: the AMR surface band (surface.refine_depth / refine_height) is retired with the column path\n")

    # ---- ParseSpallSettings: the column CSVs and the two P0 keys ------------------------------------------------
    s = between(s, "        // CSV diagnostic toggles (default off). Parsed before the early-return\n",
                """        // A1: every-removal event log (validated at the end of Parse, once
        // spallation.model is known).
        csv_on = 0;
""", """        // A1: every-removal event log (validated at the end of Parse, once
        // spallation.model is known). P2-28: the clusters and h_col CSVs are retired with the column path.
        int csv_on = 0;
""")
    s = between(s, "        // P0. BOTH keys are parsed BEFORE the spall.enabled early-return, for\n",
                """            Util::Abort(INFO, "spall.surface_T must be 'cell_average', "
                              "'diagnostic' or 'reconstructed', got: ", surfT);
""", """        // P2-28: spall.surface_T is retired at `reconstructed` (the table in Parse): the face list's
        // R1 IS the reconstruction; site A's cell-average and diagnostic forms are deleted.
""")
    s = between(s, "        std::string remE;\n",
                """        // reconstruction) and of the scope (both keys need a temperature-derived
        // stress source, and sp_weibull) are in CheckInterlocks.
""", """        // P2-28: spall.removal_energy is retired at `flake` (the table in Parse): the face list's R3
        // IS the flake remap.
""")

    # ---- members ------------------------------------------------------------------------------------------------
    s = sub(s, """    // AMR surface band (m): cells kept on the finest level below / above each
    // column's top solid cell. depth > 0 is required when max_level > 0.
    Set::Scalar surface_refine_depth  = 0.0;
    Set::Scalar surface_refine_height = 0.0;
""", "")
    s = between(s, "    // C1: pinned-surface closure state (lev 0, ncols-sized host vectors,\n",
                "    Set::Scalar pinned_cols_face_minrule = 0.0;\n", "")
    s = sub(s, """    // T_flame everywhere. Read by the kernel, BuildPinEffective, EnergyLedger
    // and PatchDiagnostics.""", """    // T_flame everywhere. Read by the kernel, EnergyLedger and the face list.""")
    s = sub(s, """    Set::Scalar patch_P_robin      = 0.0; // thermo: sum q_robin*dx*dy over patch top cells [W]
""", "")
    s = between(s, "    // Cell-centred per-column h_spall diagnostic, written at (i, j, k_top)\n",
                "    Set::Field<Set::Scalar> crack_event_count_mf;\n", "")
    s = between(s, "    // Step 19: mechanical bit feed (drawworks). When bit.enabled, a rigid flat\n",
                "    Set::Scalar bit_start_time = 0.0;   // [s] time the feed begins\n", "")
    s = between(s, "    // Optional prescribed sigma_xx(x, y, z, t) override for verification.\n",
                "    Set::Scalar recon_dT_mean     = 0.0;  // mean (T_surf - T_top) [K]\n", "")
    s = sub(s, """    // P0b witnesses, registered only under robin_form = surface.
    // surf_side_share is deliverable 1.2's requirement: side faces keep the
    // cell average, so the Meier rerun is only interpretable if the share of
    // P_robin they carry is on the record.
    Set::Scalar surf_cols_robin = 0.0;
    Set::Scalar surf_q_mean     = 0.0;
    Set::Scalar surf_side_share = 0.0;""", """    // The share of P_robin the side faces carry (the ledger; registered only
    // under robin_form = surface). P2-28: surf_cols_robin / surf_q_mean, the
    // column reconstruction's witnesses, are deleted.
    Set::Scalar surf_side_share = 0.0;""")
    s = between(s, "    // P0b deliverable 1.1: ONE surface temperature per column per step, built\n",
                "    amrex::Vector<Set::Scalar> rc_ell_prev;\n", "")
    s = between(s, "    // Step 15c: connected-component spall labeller + detach mode gate.\n",
                "    bool sp_h_col_events_csv_header_written = false;\n", "")
    s = sub(s, """    // edge-triggered h_col CSV above, it is a complete removal record.""",
            """    // column path's retired edge-triggered CSVs, it is a complete removal record.""")
    s = between(s, "    // Compute + log. False => every path below is skipped and the run is\n",
                "    bool sp_removal_energy_flake = false;\n",
                "    // P2-28: the three P0 flags (sp_surface_T_recon / _apply, sp_removal_energy_flake) are deleted:\n"
                "    // the face list always reconstructs (R1) and always removes a flake (R3).\n")

    # ---- ParseSpallationModel: site A's hooks, the detach mode ----------------------------------------------------
    s = between(s, "        int presc_on = 0;\n",
                """            value.sp_prescribed_T_surface_f =
                value.sp_prescribed_T_surface_parser.compile<1>();
        }

""", "")
    s = between(s, "            // Step 15c: detach-mode gate + cluster-area threshold. Default\n",
                """                        "detach_mode = connected_cluster, got ",
                        value.weibull_A_crit);
            }
""", "            // P2-28: weibull.detach_mode is retired at per_face, weibull.A_crit aborts (Parse)\n")

    # ---- InitializeFlaws: the A_crit warning --------------------------------------------------------------------
    s = between(s, "        // Resolvability check for the connected_cluster lateral gate. A_crit is\n",
                """            weibull_A_crit_warned = true;
        }

""", "")

    # ---- SurfacePatch -------------------------------------------------------------------------------------------
    s = sub(s, "        enum class RobinForm { Cell, Face, Pinned, Surface };",
            "        // P2-28: `pinned` (the column pin rule) is retired -- it aborts by value.\n"
            "        enum class RobinForm { Cell, Face, Surface };")
    s = sub(s, "        Set::Scalar pinned_idle_cycles = 2.0; // C1: robin_form = pinned only\n", "")
    s = sub(s, """            else if (robin_str == "pinned") v.robin_form = RobinForm::Pinned;
            else if (robin_str == "surface") v.robin_form = RobinForm::Surface;
            else Util::Abort(INFO, "surface_patch.robin_form must be 'cell', 'face', "
                                   "'pinned' or 'surface', got: ", robin_str);""",
            """            else if (robin_str == "surface") v.robin_form = RobinForm::Surface;
            else if (robin_str == "pinned")
                Util::Abort(INFO, "surface_patch.robin_form = pinned was retired in P2-28 (the pin rule "
                                  "and the column path are deleted): use surface (the face list) or "
                                  "cell / face (a deck without removal)");
            else Util::Abort(INFO, "surface_patch.robin_form must be 'cell', 'face' "
                                   "or 'surface', got: ", robin_str);""")
    s = between(s, "            // C1: pinned-surface closure (see SurfaceCellFlux). A column that\n",
                """                    Util::Abort(INFO, "surface_patch.pinned_idle_cycles must be finite "
                                      "and >= 0, got: ", v.pinned_idle_cycles);
            }
""", """            if (v.robin_form == RobinForm::Surface)
            {
                // P2-28: a PRESCRIBED patch on the face list keeps the list's form (its non-
                // prescribed faces are billed by FrFaceFlux at h = 0: their losses only)
                if (v.mode != Mode::ConvectiveFlame && !v.face_list)
                    Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                      "surface_patch.mode = convective_flame");
                // P2-28: pinned_idle_cycles governed the pin rule, retired with robin_form = pinned.
                // Read here, where the strict scope would otherwise abort, and reported inert.
                if (pp.contains("pinned_idle_cycles"))
                {
                    Set::Scalar pic = 0.0;
                    pp.query_required("pinned_idle_cycles", pic);
                    Util::Message(INFO, "surface_patch.pinned_idle_cycles = ", pic, " is INERT: the pin "
                                  "rule it governed was retired in P2-28 (robin_form = pinned)");
                }
            }
""")

    # ---- the include --------------------------------------------------------------------------------------------
    s = sub(s, """    // The removal code lives in three class-body fragments, textually included
    // here at the original position of UpdateRemovalAfterCohesive (R1 refactor):
    // Shared.H (P2-18: the helpers and the step prologue BOTH paths call),
    // Removal.H (the column path, face_list = 0) and FaceRemoval.H (the line
    // path, face_list = 1, production).
    #include "Integrator/MMWSpalling/Shared.H"
    #include "Integrator/MMWSpalling/Removal.H"
    #include "Integrator/MMWSpalling/FaceRemoval.H\"""", """    // The removal code lives in two class-body fragments, textually included
    // here: Shared.H (the helpers and the step prologue) and FaceRemoval.H (the
    // face list, the ONE removal path since P2-28 deleted the column path).
    #include "Integrator/MMWSpalling/Shared.H"
    #include "Integrator/MMWSpalling/FaceRemoval.H\"""")
    return scrub(s)


# =========================================================================================== FaceRemoval.H ====
def face(s):
    s = sub(s, """            SurfaceCellFlux(true, true, 1.0, el.h, Tg, kc, dx, Tc,
                            eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                            a, b, c, d, 0.0, nullptr, 0.0, &Tf);""",
            """            SurfaceCellFlux(true, true, 1.0, el.h, Tg, kc, dx, Tc,
                            eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                            a, b, c, d, 0.0, &Tf);""")
    s = sub(s, """                        q_in, q_out, q_rob, q_los, 0.0, nullptr, Ts);""",
            """                        q_in, q_out, q_rob, q_los, Ts);""")
    s = sub(s, """                                    a1, a2, a3, a4, 0.0, nullptr, 0.0, &Tf);""",
            """                                    a1, a2, a3, a4, 0.0, &Tf);""")
    s = sub(s, """    // ---------------------------------------------------------------------
    // The dispatch target: under the key this IS the removal step.
    // ---------------------------------------------------------------------
    void RecedeAlongLines(""", """    // ---------------------------------------------------------------------
    // ⚠ P2-28 (census H): THE removal step -- the face list is the ONE removal
    // path (the column path is deleted). The prologue both paths shared
    // (Shared.H), the column tops, then the lines. Called by Advance once per
    // step, BEFORE the thermal swap.
    // ---------------------------------------------------------------------
    void RemovalStep(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (!RemovalStepBegins(lev)) return;
        ResetRemovalStepFields(lev);
        RecedeAlongLines(lev, time, dt, ColumnTopSolid(lev));
    }

    // ---------------------------------------------------------------------
    // The lines: R0-R3, the vapour rule, the ONE WRITE and the event log.
    // ---------------------------------------------------------------------
    void RecedeAlongLines(""")
    return scrub(s)


def shared(s):
    return scrub(s)


def main():
    check = "--check" in sys.argv
    out = {}
    for name, fn in (("MMWSpalling.H", mmw), ("MMWSpalling/FaceRemoval.H", face), ("MMWSpalling/Shared.H", shared)):
        out[name] = fn(open(os.path.join(PRE, name)).read())
    for l in LOG:
        print("  " + l)
    # R11 on the transformed text (+ any other src/ file)
    others = []
    for dp, _, fs in os.walk(os.path.join(ROOT, "src")):
        for f in fs:
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, SRC)
            if rel in out or rel == os.path.join("MMWSpalling", "Removal.H"):
                continue
            if f.endswith((".H", ".cpp", ".cc")):
                others.append(open(p, errors="replace").read())
    allsrc = "\n".join(list(out.values()) + others)
    cnt = lambda n: len(re.findall(r"(?<![A-Za-z_])" + re.escape(n), allsrc))
    bad = [(n, cnt(n)) for n in R11_NAMES if cnt(n)]
    print("R11 (column-only names in src/ after the transform): %s" % ("none" if not bad else bad))
    lines = {k: v.count("\n") for k, v in out.items()}
    pre = {k: open(os.path.join(PRE, k)).read().count("\n") for k in out}
    pre_r = open(os.path.join(PRE, "MMWSpalling", "Removal.H")).read().count("\n")
    for k in out:
        print("  %-26s %6d -> %6d  (%+d)" % (k, pre[k], lines[k], lines[k] - pre[k]))
    print("  %-26s %6d -> deleted" % ("MMWSpalling/Removal.H", pre_r))
    print("  total change %+d lines" % (sum(lines.values()) - sum(pre.values()) - pre_r))
    if bad:
        sys.exit(1)
    if check:
        return
    for k, v in out.items():
        open(os.path.join(SRC, k), "w").write(v)
    rp = os.path.join(SRC, "MMWSpalling", "Removal.H")
    if os.path.exists(rp):
        os.remove(rp)
    print("written; Removal.H deleted (its pre-image is src_pre_2/src/Integrator/MMWSpalling/Removal.H)")


main()
