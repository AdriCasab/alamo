#!/usr/bin/env python3
"""P2-26 stage 2 (census C2 + C3) -- the source edits, as exact (old, new) pairs, each required to match ONCE.

  C2  spall.surface_normal (parse, message, member, the conditional thermo column sn_recession_superseded,
      ColumnSurfaceCos, SP_SN_COS_MIN, site A's /cos(theta) sampling, the column path's sn_cos_col / sn_cth / sn_icth
      and the column's h_col/cos recession, DepthScanH's sn_cth / sn_icth parameters and both callers' arguments,
      the face list's literal 1.0, 1.0) and spall.flake_coherence_length (parse, member, the column path's cap,
      the face-list refusal); both join the retired-at-default table at "0"
  C3  spall.melt_aware + melt.lambda_crit (parse, members, MeltSkeleton, the face-list refusal, the Lambda_L gather,
      LL_at_didx, f_skel, ScanSigmaT's f_skel parameter and every caller's argument); melt_aware joins the table at
      "0", lambda_crit at "*"
Bit-safety of what is deleted (ledger L6, L7): with the key off sn_cth = sn_icth = 1.0 and f_skel = 1.0, so every
deleted multiplication was x * 1.0 == x; deleting `1.0 * (a*b + c)` leaves the inner expression (and its FMA
contraction) unchanged.
Usage: apply_stage2.py --dry | --apply
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
M = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
R = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "Removal.H")
F = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "FaceRemoval.H")
S = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "Shared.H")

EDITS = []


def E(path, what, old, new=""):
    EDITS.append((path, what, old, new))


# ============================================================== MMWSpalling.H ====
E(M, "table: the C2/C3 keys retired at a default", '''            {"weibull.V_ref",        "*",
             "P2-26: read only under the deleted size_volume = flux"},
        };
''', '''            {"weibull.V_ref",        "*",
             "P2-26: read only under the deleted size_volume = flux"},
            // P2-26 (census C2, C3): closed, and deleted with one declared re-pin
            {"spall.surface_normal",          "0",
             "P2-26: the Step 20 surface-normal frame; its recession half was superseded at D2z-1 (user Q7: retire)"},
            {"spall.flake_coherence_length",  "0",
             "P2-26: the column path's coherence cap, a needle regulariser of the vertical-only law"},
            {"spall.melt_aware",              "0",
             "P2-26: the P1 melt skeleton f(Lambda_L) (user Q6: delete); the face list's melt rule is P2-24's"},
            {"melt.lambda_crit",              "*",
             "P2-26: read only under the deleted spall.melt_aware"},
        };
''')
E(M, "C2: the conditional thermo column", '''            // D2z-1 Goal 4 (revised): NO firings_per_void column. Under
            // whole-cell it is 1.000 by construction, so it carries no
            // information, and it appended after the eight jet_* columns and
            // broke jet_enthalpy criterion (a). The invariant lives in the
            // n_voided == 1 assertion in Removal.H instead. Registered here is
            // only the conditional Step 20 witness, so the tree audit is
            // absolute: no column appears or disappears anywhere else.
            if (value.spall_surface_normal)
                value.RegisterIntegratedVariable(&value.sn_recession_superseded,
                                                 "sn_recession_superseded", false);
        }
''', '''            // D2z-1 Goal 4 (revised): NO firings_per_void column. Under
            // whole-cell it is 1.000 by construction, so it carries no
            // information, and it appended after the eight jet_* columns and
            // broke jet_enthalpy criterion (a). The invariant lives in the
            // n_voided == 1 assertion in Removal.H instead.
        }
''')
E(M, "C3 + C2: the face-list refusals", '''            if (value.melt_aware) fr_abort("spall.melt_aware", "the MMW step");
            if (value.spall_surface_normal) fr_abort("spall.surface_normal = 1", "a later step");
''')
E(M, "C2: the coherence refusal", '''            if (value.spall_flake_coh_len > 0.0)
                fr_abort("spall.flake_coherence_length > 0 (a per-column cap)", "a later step");
''')
_sn_start = '''        // Step 20: surface-normal (slope-aware) spall CRITERION.
'''
_sn_end = '''                                "thermo column sn_recession_superseded (1.0 = "
                                "nothing lost on that step).");

'''
E(M, "C2: the surface_normal parse and message", ("SPAN", _sn_start, _sn_end))
_melt_start = '''        // P1: melt-aware spallation (default off, byte-identical when off).
'''
_melt_end = '''            Util::Abort(INFO, "melt.lambda_crit must be in (0, 1], got ",
                        value.melt_lambda_crit);

'''
E(M, "C3: the melt_aware / lambda_crit parse", ("SPAN", _melt_start, _melt_end))
E(M, "C2: the coherence parse", '''        pp.query_default("spall.flake_coherence_length", value.spall_flake_coh_len, 0.0);
        if (!(std::isfinite(value.spall_flake_coh_len) && value.spall_flake_coh_len >= 0.0))
            Util::Abort(INFO, "spall.flake_coherence_length must be finite and >= 0, got ",
                        value.spall_flake_coh_len);

''')
E(M, "C2: ColumnSurfaceCos", ("SPAN", '''    // Step 20: per-column cos(theta) of the local cavity-wall inclination,
''', '''            cos_out[static_cast<std::size_t>(cj) * nx + ci] = c;
        }
    }

'''))
E(M, "C2: site A's per-column cos", '''        // Step 20: surface-normal frame -- per-column cos(theta) so the K_I
        // scan below samples the vertical T/stress profile at depth/cos(theta)
        // (the crack runs along the local surface normal). nullptr when off:
        // the conversion multiplies by exactly 1.0 and stays bit-identical.
        amrex::Vector<Set::Scalar> sn_cos_col;
        if (spall_surface_normal)
            ColumnSurfaceCos(k_top_col.data(), nx_dom, ny_dom,
                             DX[0], DX[1], dz, domain.smallEnd(2), sn_cos_col);
        const Set::Scalar* sn_cos_p =
            sn_cos_col.empty() ? nullptr : sn_cos_col.data();

''')
E(M, "C2: site A's icth", '''                // Step 20: under surface_normal the crack frame is the local
                // surface normal; `depth` below is a NORMAL distance and the
                // vertical profile is sampled at depth/cos(theta). icth = 1
                // exactly when the flag is off (bit-identical).
                const Set::Scalar icth = sn_cos_p
                    ? 1.0 / sn_cos_p[col_index_sp(i, j)] : 1.0;

''')
E(M, "C2: site A's vertical sample depth", '''                    const Set::Scalar dv = depth * icth;  // vertical sample depth
''', '''                    const Set::Scalar dv = depth;   // the sample depth below the surface
''')
E(M, "C2: site A's K_Ic tip depth (comment)", '''                // The tip of the surface-breaking flaw sits at depth a_loc
                // (normal), i.e. vertical depth a_loc*icth, where the rock is
''', '''                // The tip of the surface-breaking flaw sits at depth a_loc,
                // where the rock is
''')
E(M, "C2: site A's K_Ic tip depth", '''                    rec_ok ? T_recon_at(a_loc * icth) : Tsurf;
''', '''                    rec_ok ? T_recon_at(a_loc) : Tsurf;
''')
E(M, "C2: the coherence member", '''    // Flake coherence length [m] (rock-mechanics cavity-wall coherence). Default
    // 0 = off. When > 0, a firing column may recede at most this far below its
    // shallowest live 4-neighbour, so the spall cavity floor has a bounded wall
    // slope: a smooth bowl can descend gradually (<= L_coh per lateral cell) but
    // a free-standing high-aspect-ratio needle cannot form (the runaway column
    // is capped until its neighbours catch up). Decomposition-independent
    // (operates on the gathered per-column surface heights).
    Set::Scalar spall_flake_coh_len = 0.0;

''')
E(M, "C2: the sn_recession_superseded member", '''    // Registered ONLY when spall.surface_normal = 1, so it appears in exactly
    // the four inputs it describes and changes no other thermo.dat header.
    // Mean 1/cos(theta) over each step's spall-firing columns = the factor by
    // which recession WOULD have been scaled before whole-cell discarded the
    // scan's result. 1.0 means nothing was lost.
    Set::Scalar sn_recession_superseded = 0.0;
''')
E(M, "C2 + C3: members, MeltSkeleton, SP_SN_COS_MIN", ("SPAN", '''    // Step 20: surface-normal (slope-aware) spall criterion + removal. When
''', '''    static constexpr Set::Scalar SP_SN_COS_MIN = 0.15;
'''))

# ================================================================== Removal.H ====
E(R, "C3: the Lambda_L gather", ("SPAN", '''        // P1 (melt-aware): gather each column's liquid-fraction profile the
''', '''        const Set::Scalar* LL_col_p =
            LL_col_prof.empty() ? nullptr : LL_col_prof.data();

'''))
E(R, "C2: the column path's per-column cos", '''        // Step 20: per-column cos(theta) of the cavity-wall inclination for
        // the surface-normal spall frame (computed from the gathered tops,
        // so identical on every rank). nullptr when the flag is off -- every
        // conversion below then multiplies by exactly 1.0 (bit-identical).
        amrex::Vector<Set::Scalar> sn_cos_col;
        if (spall_surface_normal && sp_weibull_spall)
            ColumnSurfaceCos(k_top_global.data(), nx_dom, ny_dom,
                             DX[0], DX[1], dz, domain.smallEnd(2), sn_cos_col);
        const Set::Scalar* sn_cos_p =
            sn_cos_col.empty() ? nullptr : sn_cos_col.data();

''')
E(R, "C2: sn_cth / sn_icth", '''                                // Step 20: surface-normal frame factors.
                                // Depth variables in this scan (z_crack, d,
                                // h_col) are NORMAL distances when the flag
                                // is on; the cell at k_top - s sits at
                                // normal depth s*dz*cos(theta). Exactly 1.0
                                // when off (bit-identical arithmetic).
                                const Set::Scalar sn_cth =
                                    sn_cos_p ? sn_cos_p[idx] : 1.0;
                                const Set::Scalar sn_icth = 1.0 / sn_cth;
''')
E(R, "C3: LL_at_didx", ("SPAN", '''                                // P1: liquid fraction Lambda_L by depth-below-
''', '''                                    return 0.0;
                                };

'''))
E(R, "C2 + C3: z_v and f_skel in the scan's stress", '''                                        // z_abs is a NORMAL distance under
                                        // Step 20; z_v is the vertical depth
                                        // at which the column profile holds
                                        // the same material point.
                                        const Set::Scalar z_abs =
                                            z_crack + d;
                                        const Set::Scalar z_v =
                                            z_abs * sn_icth;
                                        // P1: load-bearing skeleton factor at
                                        // this depth's cell. f == 1.0 (off)
                                        // leaves every return byte-identical;
                                        // in melt f -> 0 kills the elastic
                                        // driving stress (equivalently
                                        // E_eff = f*E), so a liquid cell carries
                                        // no crack-driving stress and the scan
                                        // breaks there instead of spalling it.
                                        Set::Scalar f_skel = 1.0;
                                        if (melt_aware)
                                        {
                                            int d_ll = static_cast<int>(
                                                std::floor(z_v / dz));
                                            f_skel = MeltSkeleton(
                                                LL_at_didx(d_ll));
                                        }
''', '''                                        // z_v: the depth below the surface at
                                        // which the column profile is sampled.
                                        const Set::Scalar z_v =
                                            z_crack + d;
''')
E(R, "C3: the stress_mf return", '''                                            return f_skel * (sx(z_v)
                                                   + sp_conf_offset);
''', '''                                            return sx(z_v) + sp_conf_offset;
''')
E(R, "C3: the ScanSigmaT call", '''                                        return ScanSigmaT(scan_g, T_z, T0_z, f_skel,
                                                          sp_conf_offset, phase_E_loc,
''', '''                                        return ScanSigmaT(scan_g, T_z, T0_z,
                                                          sp_conf_offset, phase_E_loc,
''')
E(R, "C2: the DepthScanH call", '''                                    max_steps, k_top, scan_kbot, dz, sn_cth, sn_icth,
''', '''                                    max_steps, k_top, scan_kbot, dz,
''')
E(R, "C2: the h_col/cos recession", '''                                    // Step 20: h_col is the flake thickness
                                    // along the surface normal; the column
                                    // (height-function) representation must
                                    // recede VERTICALLY by h_col/cos(theta)
                                    // to remove that flake. The graph form
                                    // of dphi/dt = -v_n*|grad phi|. The CSV
                                    // and h_spall_field diagnostics keep the
                                    // physical (normal) thickness.
                                    const Set::Scalar h_vert =
                                        h_col * sn_icth;
''', '''                                    const Set::Scalar h_vert = h_col;
''')
E(R, "C2: the coherence cap", ("SPAN", '''        // Flake coherence-length cap (rock-mechanics cavity-wall coherence):
''', '''                    if (cap_h <= 0.0) spall_fires_col[idx] = 0;
                }
            }
        }

'''))
E(R, "C2: the D2z-1 note on Step 20's recession", '''            //
            // ⚠ STEP 20 LOSES ITS RECESSION SCALING HERE, NOT ITS CRITERION.
            // h_c = dz discards the scan's result, and the slope factor lives
            // inside the scan (sn_icth, :691 and :992), while recession is a
            // plain phi shift with no cos(theta) at its own site (:1444). So a
            // wall column drops ONE CELL PER FIRING whatever its slope. The
            // slope-aware FIRING DECISION is untouched -- it samples K_I at
            // depth/cos(theta) in the criterion update (MMWSpalling.H:4668),
            // which never goes through the scan. The capability degrades; it
            // does not die. sn_recession_superseded reports it on every run
            // with spall.surface_normal = 1.
''')
E(R, "C2: the sn_recession_superseded value", '''                // disappears, with no exception carved out for our own.
                //
                // sn_recession_superseded: the mean 1/cos(theta) over this
                // step's spall-firing columns -- the factor by which recession
                // WOULD have been scaled at :992 before whole-cell discarded
                // the scan's result. 1.0 means nothing was lost on this step.
                // Registered only when spall.surface_normal = 1, so it appears
                // in exactly the runs it describes. sn_cos_p is non-null under
                // precisely that key (see :256-261).
                if (sn_cos_p)
                {
                    Set::Scalar icth_sum = 0.0, icth_n = 0.0;
                    for (int idx = 0; idx < ncols; ++idx)
                    {
                        if (regime_col[idx] != 1) continue;
                        if (sn_cos_p[idx] > 0.0)
                        {
                            icth_sum += 1.0 / sn_cos_p[idx];
                            icth_n   += 1.0;
                        }
                    }
                    sn_recession_superseded =
                        (icth_n > 0.0) ? icth_sum / icth_n : 0.0;
                }
''', '''                // disappears, with no exception carved out for our own.
''')

# =================================================================== Shared.H ====
E(S, "C3: ScanSigmaT's f_skel parameter", '''    // The T-derived driving stress at one sample (1-D confinement,
    // compression-positive): the rock law through LawSigma when one is
    // selected, else the scalar thermoelastic form. f_skel = 1 off melt.
    Set::Scalar ScanSigmaT(const Numeric::Material::Grimsel *scan_g,
                           Set::Scalar T_z, Set::Scalar T0_z, Set::Scalar f_skel,
                           Set::Scalar sp_conf_offset, Set::Scalar phase_E_loc,
                           Set::Scalar phase_beta_loc, Set::Scalar phase_nu_loc) const
    {
        if (scan_g)
            return f_skel
                 * (sp_relief
                    * LawSigma(scan_g, T_z, T0_z)
                    + sp_conf_offset);
        // ⚠ pre-P2-2b expression, VERBATIM.
        return f_skel * (sp_relief * phase_E_loc
               * phase_beta_loc
               * (T_z - T0_z)
               / (1.0 - phase_nu_loc)
               + sp_conf_offset);
    }
''', '''    // The T-derived driving stress at one sample (1-D confinement,
    // compression-positive): the rock law through LawSigma when one is
    // selected, else the scalar thermoelastic form.
    Set::Scalar ScanSigmaT(const Numeric::Material::Grimsel *scan_g,
                           Set::Scalar T_z, Set::Scalar T0_z,
                           Set::Scalar sp_conf_offset, Set::Scalar phase_E_loc,
                           Set::Scalar phase_beta_loc, Set::Scalar phase_nu_loc) const
    {
        if (scan_g)
            return sp_relief
                   * LawSigma(scan_g, T_z, T0_z)
                   + sp_conf_offset;
        // ⚠ pre-P2-2b expression, VERBATIM.
        return sp_relief * phase_E_loc
               * phase_beta_loc
               * (T_z - T0_z)
               / (1.0 - phase_nu_loc)
               + sp_conf_offset;
    }
''')
E(S, "C2: DepthScanH's signature", '''    static Set::Scalar DepthScanH(int max_steps, int k_top, int scan_kbot,
                                  Set::Scalar dz, Set::Scalar sn_cth,
                                  Set::Scalar sn_icth, Set::Scalar a_f, bool rec_ok,
''', '''    static Set::Scalar DepthScanH(int max_steps, int k_top, int scan_kbot,
                                  Set::Scalar dz, Set::Scalar a_f, bool rec_ok,
''')
E(S, "C2: DepthScanH z_crack", '''            const Set::Scalar z_crack =
                Set::Scalar(s) * dz * sn_cth;
''', '''            const Set::Scalar z_crack =
                Set::Scalar(s) * dz;
''')
E(S, "C2: DepthScanH K_Ic at the tip", '''                ? kic_pa_at_z((z_crack + a_f) * sn_icth)
''', '''                ? kic_pa_at_z(z_crack + a_f)
''')
E(S, "C2: DepthScanH bracket", '''            Set::Scalar z_lo = Set::Scalar(s_last_pass)
                               * dz * sn_cth; // pass
            Set::Scalar z_hi = Set::Scalar(s_last_pass + 1)
                               * dz * sn_cth; // fail
''', '''            Set::Scalar z_lo = Set::Scalar(s_last_pass)
                               * dz; // pass
            Set::Scalar z_hi = Set::Scalar(s_last_pass + 1)
                               * dz; // fail
''')
E(S, "C2: DepthScanH bisection cell", '''                        std::floor(z_mid * sn_icth / dz));
''', '''                        std::floor(z_mid / dz));
''')
E(S, "C2: DepthScanH bisection K_Ic", '''                    ? kic_pa_at_z((z_mid + a_f) * sn_icth)
''', '''                    ? kic_pa_at_z(z_mid + a_f)
''')

# ============================================================== FaceRemoval.H ====
E(F, "C3: the face list's ScanSigmaT (scan)", '''                        return ScanSigmaT(scan_g, T_z, T0_z, 1.0, conf_off, pE, pB, pNu);
''', '''                        return ScanSigmaT(scan_g, T_z, T0_z, conf_off, pE, pB, pNu);
''')
E(F, "C2: the face list's DepthScanH literals", '''                    o.h_col = DepthScanH(nprof, 0, -(L.Lg + 1), dx, 1.0, 1.0, a_f, true,
''', '''                    o.h_col = DepthScanH(nprof, 0, -(L.Lg + 1), dx, a_f, true,
''')
E(F, "C3: the face list's ScanSigmaT (surface stress)", '''                    const Set::Scalar sig_s = ScanSigmaT(scan_g, T_at_zv(0.0), T0_at_d(0), 1.0,
                                                         conf_off, pE, pB, pNu);
''', '''                    const Set::Scalar sig_s = ScanSigmaT(scan_g, T_at_zv(0.0), T0_at_d(0),
                                                         conf_off, pE, pB, pNu);
''')


def main():
    dry = "--apply" not in sys.argv
    text = {p: open(p).read() for p in (M, R, F, S)}
    bad = 0
    for path, what, old, new in EDITS:
        s = text[path]
        if isinstance(old, tuple):            # ("SPAN", start, end): delete start..end inclusive
            _, a, b = old
            na = s.count(a)
            i = s.find(a)
            j = s.find(b, i) if i >= 0 else -1
            nb = s.count(b)
            ok = na == 1 and nb == 1 and j > i
            print("  %-4s %-14s %-58s span %d line(s)" % ("ok" if ok else "FAIL", os.path.basename(path), what,
                                                        s[i:j + len(b)].count("\n") if ok else -1))
            if ok:
                text[path] = s[:i] + new + s[j + len(b):]
        else:
            n = s.count(old)
            ok = n == 1
            print("  %-4s %-14s %-58s [%d match]" % ("ok" if ok else "FAIL", os.path.basename(path), what, n))
            if ok:
                text[path] = s.replace(old, new, 1)
        bad += 0 if ok else 1
    if bad:
        print("FAIL: %d edit(s) did not match exactly once; nothing written" % bad)
        return 1
    if not dry:
        for p, t in text.items():
            open(p, "w").write(t)
        print("APPLIED %d edits" % len(EDITS))
    else:
        print("DRY RUN: %d edits match exactly once" % len(EDITS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
