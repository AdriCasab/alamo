    // S1: whole-domain energy ledger at level 0, called after the H update.
    //   E_H   = sum_all H*V (void cells included; their H is frozen)
    //   E_src += dt*(P_beam + P_robin - P_loss)
    //   ledger_err_raw = (E_H - E_H0 - E_src)/max(|E_src|, 1e-30)
    // P_beam reuses the E1 BeamEnergyBalance sum (same BeamSource calls);
    // P_robin/P_loss re-evaluate SurfaceCellFlux with the H update's inputs
    // (prescribed_T patch flux, if any, is counted in P_robin).
    // EXACT only with Neumann-0 outer temperature BCs on every face: a
    // Dirichlet/Robin domain BC puts a flux through the boundary faces that
    // the ledger does not count. The ledger state (ledger_started, E_H0,
    // E_src) is NOT checkpointed: a restarted run re-baselines at its first
    // step. E_H is a compensated (Neumaier) sum so the
    // difference E_H - E_H0 keeps the small per-step increments.
    // surface_missing_cols counts columns that have a mask-top solid cell but
    // no surface_mf cell (the phi-window dropout at phi_top ~ dz ties).
    void EnergyLedger(int lev, Set::Scalar time, Set::Scalar dt,
                      const amrex::Vector<int>& k_top_col)
    {
        if (lev != 0) return;
        // P2-28: the kernel's own switch (list_state: the evaluate-only list included)
        const bool void_aware  = list_state();
        const bool follow_mask = surface_follow_mask;
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const int khi_dom = domain.bigEnd(2);
        const int ncx = domain.length(0), ncy = domain.length(1);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz = DX[2];
        const Set::Scalar A  = DX[0] * DX[1];
        const Set::Scalar V  = A * dz;
        const PatchConst patch_c = PatchConstants(DX);
        const Set::Scalar inv_dz2 = patch_c.inv_dz2;

        const bool        sp_on    = surface_patch.enabled;
        const bool        sp_presc = patch_c.presc;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face;
        const Set::Scalar sp_r2    = patch_c.r2;
        const Set::Scalar sp_T     = (sp_on && sp_presc) ? surface_patch.T_at(time) : 0.0;

        auto neumaier = [](Set::Scalar& s, Set::Scalar& c, Set::Scalar x) {
            const Set::Scalar t = s + x;
            if (std::abs(s) >= std::abs(x)) c += (s - t) + x;
            else                            c += (x - t) + s;
            s = t;
        };

        Set::Scalar EH_s = 0.0, EH_c = 0.0, EH0_s = 0.0, EH0_c = 0.0;

        Set::Scalar EHsol_s = 0.0, EHsol_c = 0.0;

        Set::Scalar EH0sol_s = 0.0, EH0sol_c = 0.0;
        Set::Scalar P_robin = 0.0, P_loss = 0.0, P_robin_side = 0.0;
        // Per-column side power, for the wall/riser split that no plotfile can
        // give (no plotfile field is a side flux). Sized only when the key
        // is on.
        const bool sidecsv = surface_patch.side_power_interval > 0.0;
        std::vector<Set::Scalar> sidecol;
        // P2-8: under the face list the ledger mirrors the LIST, face for face.
        // P2-23: the kernel's switch (the ledger is single-level, so `owner` is
        // true). P2-28: a prescribed patch is on the list too (its top faces' T_presc)
        const bool        fl_on_l = surface_patch.face_list && void_aware && fl_ready;
        std::vector<int> has_surf(static_cast<std::size_t>(ncx) * ncy, 0);
        // rev7 d2: sized here, where ncx/ncy are known, not at the declaration.
        if (sidecsv) sidecol.assign(static_cast<std::size_t>(ncx) * ncy, 0.0);
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& H     = (*H_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& H_old = (*H_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp    = (*kappa_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& surf  = (*surface_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& rem   =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
            amrex::Array4<const int> fli_l;
            if (fl_on_l) fli_l = fl_index->const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                neumaier(EH_s, EH_c, H(i, j, k) * V);
                if (!ledger_started) neumaier(EH0_s, EH0_c, H_old(i, j, k) * V);
                const int col = (j - jlo_dom) * ncx + (i - ilo_dom);
                if (surf(i, j, k) > 0.5) has_surf[col] = 1;
                const bool self_void = void_aware && rem(i, j, k) > 0.5;
                // P0b (c''): the solid-only sum. See ledger_E_H_solid.
                if (!self_void)
                {
                    neumaier(EHsol_s, EHsol_c, H(i, j, k) * V);
                    if (!ledger_started)
                        neumaier(EH0sol_s, EH0sol_c, H_old(i, j, k) * V);
                }
                if (self_void) continue;

                const Set::Scalar Tc = T_old(i, j, k);
                const Set::Scalar kc = kp(i, j, k);
                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const int k_top_c = (follow_mask && void_aware) ? k_top_col[col] : khi_dom;
                const Set::Scalar surf_c = follow_mask
                                         ? Set::Scalar(k == k_top_c ? 1.0 : 0.0)
                                         : surf(i, j, k);
                bool in_patch = false;
                if (sp_on && k == k_top_c)
                {
                    const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                    if (dxp * dxp + dyp * dyp <= sp_r2)
                    {
                        // P2-28: the kernel's half-cell Dirichlet term, mirrored; under the
                        // face list it is the list's (T_presc, below), not added twice
                        if (sp_presc) { if (!fl_on_l) P_robin += 2.0 * kc * (sp_T - Tc) * inv_dz2 * V; }
                        else          in_patch = true;
                    }
                }
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                SurfaceCellFlux(in_patch, sp_face, surf_c,
                                sp_on ? flame_h[col]  : Set::Scalar(0.0),
                                sp_on ? flame_Tg[col] : Set::Scalar(0.0), kc, dz, Tc,
                                eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                                q_in, q_out, q_rob, q_los);
                // P2-8: the cell's listed faces (first entry, count), as the kernel reads them.
                int fl_f = 0, fl_n = 0;
                if (fl_on_l) { fl_f = fli_l(i, j, k, 0); fl_n = fli_l(i, j, k, 1); }
                // P2-9: a listed cell's faces all come through FrFaceFlux below.
                if (fl_n > 0) { q_rob = 0.0; q_los = 0.0; }
                P_robin += q_rob * A;
                P_loss  += q_los * A;
                // P2-9 R1: every listed face, mirrored from the kernel.
                for (int q = fl_f; q < fl_f + fl_n; ++q)
                {
                    const FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                    // P2-28: a PRESCRIBED face -- the kernel's term x V (a top face: no side share),
                    // in the column mirror's product order (w first)
                    if (e.T_presc > 0.0)
                    {
                        P_robin += e.w * 2.0 * kc * (e.T_presc - Tc) * inv_dz2 * V;
                        continue;
                    }
                    Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                    FrFaceFlux(e, e.Tg, kc, Tc, DX[0], a_in, a_out, a_rob, a_los);
                    const Set::Scalar Af = e.w * DX[0] * dz;
                    P_robin      += a_rob * Af;
                    P_loss       += a_los * Af;
                    if (e.cls == FL_TOP) continue;
                    P_robin_side += a_rob * Af;
                    if (sidecsv && col >= 0 && col < static_cast<int>(sidecol.size()))
                        sidecol[col] += a_rob * Af;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(EH_s);
        amrex::ParallelDescriptor::ReduceRealSum(EH_c);
        amrex::ParallelDescriptor::ReduceRealSum(EHsol_s);
        amrex::ParallelDescriptor::ReduceRealSum(EHsol_c);
        amrex::ParallelDescriptor::ReduceRealSum(P_robin);
        amrex::ParallelDescriptor::ReduceRealSum(P_loss);
        amrex::ParallelDescriptor::ReduceIntMax(has_surf.data(), static_cast<int>(has_surf.size()));
        if (!ledger_started)
        {
            amrex::ParallelDescriptor::ReduceRealSum(EH0_s);
            amrex::ParallelDescriptor::ReduceRealSum(EH0_c);
            ledger_E_H0_s  = EH0_s;
            ledger_E_H0_c  = EH0_c;
            amrex::ParallelDescriptor::ReduceRealSum(EH0sol_s);
            amrex::ParallelDescriptor::ReduceRealSum(EH0sol_c);
            ledger_E_H0_sol_s = EH0sol_s;
            ledger_E_H0_sol_c = EH0sol_c;
            ledger_started = true;
        }

        int n_missing = 0;
        for (std::size_t c = 0; c < has_surf.size(); ++c)
        {
            const bool has_solid = void_aware ? (k_top_col[c] != std::numeric_limits<int>::min())
                                              : true;
            if (has_solid && !has_surf[c]) ++n_missing;
        }

        const Set::Scalar P_beam = (beam.power_at(time) > 0.0) ? beam_P_dep : 0.0;
        ledger_P_beam  = P_beam;
        ledger_P_robin = P_robin;
        if (fl_on_l)
        {
            const Set::Scalar Ps = fl_P_top + fl_P_up + fl_P_side + fl_P_down;
            fl_id_ledger = (std::abs(Ps) > 0.0) ? std::abs(P_robin - Ps) / std::abs(Ps) : 0.0;
        }
        // P0b witness for robin_form = surface: the non-top listed faces' share
        // of the robin power (registered only with the patch on: P2-23)
        if (sp_on && surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
        {
            amrex::ParallelDescriptor::ReduceRealSum(P_robin_side);
            surf_side_share = (std::abs(P_robin) > 0.0)
                              ? P_robin_side / P_robin : 0.0;
            // rev7 d2: the per-column side-power log. The wall/riser split
            // cannot be a re-analysis -- no plotfile field is a side flux, and
            // the thermo holds only aggregates. One row per column that carries
            // side power, at side_power_interval seconds of simulated time.
            if (sidecsv)
            {
                amrex::ParallelDescriptor::ReduceRealSum(sidecol.data(),
                                                         static_cast<int>(sidecol.size()));
                const long widx = static_cast<long>(
                    std::floor(time / surface_patch.side_power_interval));
                if (widx > side_power_last)
                {
                    side_power_last = widx;
                    if (amrex::ParallelDescriptor::IOProcessor())
                    {
                        const std::string path = plot_file + "_side_power.csv";
                        const bool first = !side_power_header_written;
                        std::ofstream ofs(path,
                            (first ? std::ios::trunc : std::ios::app) | std::ios::out);
                        if (ofs.is_open())
                        {
                            if (first) ofs << CsvHeader(SidePowerColumns());
                            ofs.setf(std::ios::scientific);
                            ofs.precision(10);
                            const std::size_t n_side = SidePowerColumns().size();
                            for (int cj = 0; cj < ncy; ++cj)
                            for (int ci = 0; ci < ncx; ++ci)
                            {
                                const std::size_t c =
                                    static_cast<std::size_t>(cj) * ncx + ci;
                                if (!(std::abs(sidecol[c]) > 0.0)) continue;
                                const Set::Scalar cx = PLO[0]
                                    + (Set::Scalar(ci + ilo_dom) + 0.5) * DX[0]
                                    - surface_patch.x0;
                                const Set::Scalar cy = PLO[1]
                                    + (Set::Scalar(cj + jlo_dom) + 0.5) * DX[1]
                                    - surface_patch.y0;
                                std::ostringstream row;
                                row.copyfmt(ofs);
                                row << time << "," << (ci + ilo_dom) << ","
                                    << (cj + jlo_dom) << ","
                                    << std::sqrt(cx * cx + cy * cy) << ","
                                    << sidecol[c] << "\n";
                                CsvCheckRow(row.str(), n_side, "_side_power.csv");
                                ofs << row.str();
                            }
                        }
                        side_power_header_written = true;
                    }
                }
            }
        }
        ledger_P_loss  = P_loss;
        ledger_E_src  += dt * (P_beam + P_robin - P_loss);
        ledger_E_H     = EH_s + EH_c;
        const Set::Scalar dE = (EH_s - ledger_E_H0_s) + (EH_c - ledger_E_H0_c);
        ledger_err_raw = (dE - ledger_E_src) / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));
        // P0b (c''): with the flake sink, the SOLID enthalpy must close against
        // the source minus what left with the rock.
        ledger_E_H_solid = EHsol_s + EHsol_c;
        const Set::Scalar dE_sol = (EHsol_s - ledger_E_H0_sol_s)
                                 + (EHsol_c - ledger_E_H0_sol_c);
        // P2-2f: foot_mech_E is the pad-clearance sink -- enthalpy that left the
        // SOLID set mechanically. Without it this row cannot close on any leg
        // with feet + foot_body_clearance, which is one cause and not two: the
        // gap between ledger_err_raw and ledger_err_solid is the flake credit the
        // TOTAL row omits by design, and the accepted CTRL5 shows a larger gap
        // while closing here to 7.78e-04.
        // P2-24: and the vapour's enthalpy (the face list's regime 2: H V of each vaporised cell)
        ledger_err_solid = (dE_sol - (ledger_E_src - ledger_E_flake - foot_mech_E - ledger_E_vapor))
                           / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));
        ledger_surface_missing_cols = Set::Scalar(n_missing);
    }
