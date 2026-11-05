#ifndef INTEGRATOR_MMWSPALLING_REMOVAL_H
#define INTEGRATOR_MMWSPALLING_REMOVAL_H

// Fragment of class Integrator::MMWSpalling (R1 refactor). This file is
// textually #included INSIDE the class body of src/Integrator/MMWSpalling.H
// at the original location of UpdateRemovalAfterCohesive. It is a class-body
// fragment, NOT a standalone header: it must not be compiled or #included
// anywhere else, and intentionally has no #include / namespace / class of its
// own (everything it needs is already included by MMWSpalling.H before the
// class). The definition stays an inline member of MMWSpalling.

    // pass populates each rank's local maximum k_top via ReduceIntMax;
    // a second pass lets the rank that owns (i, j, k_top_global)
    // compute the spall and vapor candidates and write into the
    // per-column buffers, which are summed across ranks. The final
    // pass applies the column-wide phi shift and per-cell detach
    // bookkeeping locally on every rank using the reduced data. This
    // makes the update correct for arbitrary AMReX z-decompositions,
    // not just column-coherent partitions.
    void UpdateRemovalAfterCohesive(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (!removal_enabled()) return;
        Util::Assert(INFO, TEST(microstructure_enabled));
        // AMR: removal is decided and applied on the surface-owning level
        // only; coarser levels receive phi/removed by average-down.
        if (!SurfaceOwner(lev)) return;

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz   = DX[2];

        const amrex::Box& domain = geom[lev].Domain();
        const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
        const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
        const int nx_dom = dom_ihi - dom_ilo + 1;
        const int ny_dom = dom_jhi - dom_jlo + 1;
        const int ncols  = nx_dom * ny_dom;

        // Spall calibration ---------------------------------------------
        const bool spall_on = spall_enabled;
        const Set::Scalar D_th = spall_damage_threshold;
        const Set::Scalar xi_s = spall_xi_s;
        Set::Scalar h_spall = 0.0;
        if (spall_on)
        {
            const Set::Scalar C_h  = spall_C_h;
            const Set::Scalar alpha_eff   = spall_alpha_override ? spall_alpha   : therm_alpha;
            const Set::Scalar t_spall_eff = spall_t_spall_override ? spall_t_spall : (time + dt);
            if (!(std::isfinite(alpha_eff) && alpha_eff > 0.0))
                Util::Abort(INFO, "spall: thermal diffusivity for h_spall is non-positive (",
                            alpha_eff, "); set spall.alpha or hc.heat.alpha");
            if (!(std::isfinite(t_spall_eff) && t_spall_eff > 0.0))
                Util::Abort(INFO, "spall: t_spall for h_spall is non-positive (", t_spall_eff,
                            "); set spall.t_spall");
            h_spall = C_h * std::sqrt(alpha_eff * t_spall_eff);
            if (!(std::isfinite(h_spall) && h_spall >= 0.0))
                Util::Abort(INFO, "spall: computed h_spall is non-finite (", h_spall, ")");
        }
        const bool spall_prescribed = spall_prescribed_ratio_enabled;
        auto ratio_f = spall_ratio_f;
        // stress_mf is only registered/populated when the elastic solve runs.
        // With el.type = disable (e.g. spall.stress_source =
        // local_thermoelastic), it must not be touched; the depth scan still
        // runs via the T-derived sigma_xx path. mech_disabled gates every
        // stress_mf access; have_stress still gates whether the scan runs at
        // all (a stress source — T-derived or stress_mf — exists).
        const bool mech_disabled =
            (m_type == Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable);
        const bool have_stress = spall_on && !spall_prescribed;
        if (have_stress && !mech_disabled)
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

        // Vapor calibration ---------------------------------------------
        const bool vapor_on = vapor_enabled;
        const Set::Scalar A_beam = vapor_A_beam;
        const Set::Scalar T0_v   = vapor_T0;
        const Set::Scalar P_now  = vapor_on ? beam.power_at(time) : 0.0;
        const Set::Scalar oneR   = 1.0 - beam.R;
        // Step 17: local Gaussian-beam parameters for per-cell v_vap
        // under sp_weibull. Identical to the values used by the heat
        // source path; `beam.power_at(time)` is already captured in
        // `P_now`. At the surface (z = z0) the beam waist is w0 with
        // no spreading; per-area surface power density is then
        // (2·P / (π·w0²)) · exp(-2·r² / w0²). Captured outside the
        // MFIter so the inner block can read them.
        const Set::Scalar bw0_v  = beam.omega0;
        const Set::Scalar bx0_v  = beam.x0;
        const Set::Scalar by0_v  = beam.y0;
        const bool vapor_prescribed = vapor_prescribed_rop_enabled;
        auto vapor_rop_eval = vapor_rop_f;
        if (vapor_on && !(std::isfinite(A_beam) && A_beam > 0.0))
            Util::Abort(INFO, "vapor: A_beam is non-positive (", A_beam,
                        "); set vapor.A_beam");

        // Per-step diagnostic fields are reset here. `phi` and `removed`
        // carry state across steps and are not reset.
        regime_mf[lev]->setVal(0.0);
        RoP_mf[lev]->setVal(0.0);
        if (spall_on)
        {
            spall_event_mf[lev]->setVal(0.0);
            spall_thickness_mf[lev]->setVal(0.0);
            RoP_spall_mf[lev]->setVal(0.0);
        }
        if (vapor_on)
        {
            vapor_event_mf[lev]->setVal(0.0);
            vapor_thickness_mf[lev]->setVal(0.0);
            RoP_vap_mf[lev]->setVal(0.0);
        }
        // Step 17: per-step reset for the per-cell regime field under
        // sp_weibull. The field is registered only when
        // spallation_model == SpWeibull && (spall_on || vapor_on).
        if (spallation_model == SpallationModel::SpWeibull &&
            (spall_on || vapor_on))
        {
            regime_field_mf[lev]->setVal(0.0);
        }
        // Step 16: sp_weibull spall path flag. The h_spall_field_mf
        // diagnostic accumulates as a running max per column (written
        // at the INITIAL top-z slice z = top_k_dom, which doesn't move
        // as the surface erodes) — see PASS 2 write. The field is NOT
        // reset per step, so the last plotfile records the deepest h_col
        // every column has ever produced.
        const bool sp_weibull_spall =
            spall_on && spallation_model == SpallationModel::SpWeibull;

        const Numeric::Material::Table *ptabs =
            phase_tables.empty() ? nullptr : phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        const Set::Scalar dt_local = dt;

        auto col_index = [=](int i, int j) {
            return (j - dom_jlo) * nx_dom + (i - dom_ilo);
        };

        // D2c: burner-body clearance (surface_patch.foot_body_clearance).
        // FeetDescent (inside the previous step's thermal advance, which runs
        // after this pass) recorded, per annulus column standing proud of its
        // pad, the pad height's cell index k_c (foot_clear_k). Here, before
        // PASS 1, every solid cell of that column with k > k_c is flipped to
        // removed -- regime 4 (mechanical), D = 0, thermal state kept (void H
        // is frozen and still counted by EnergyLedger) -- and the column's phi
        // is shifted by n*dz, n = k_top - k_c, so the new top is k_c with its
        // sub-cell phi unchanged. Cells are flipped by index, not by the sign
        // of phi (round-off ties phi_top = 0 / dz would otherwise miss or add
        // one); a flipped cell left at phi >= 0 by a tie gets phi = -tiny,
        // and k_c is held at phi >= 0. The pin state (pin_T/pin_t) is left
        // untouched. The thermal passes below then see the cleared surface.
        // foot_mech_vol / foot_mech_cols log the flips; removal_events.csv
        // gets one regime-4 row per cleared column (k_top = the top before
        // the cut, h_applied = n*dz, n_voided = cells flipped, other fields 0).
        //
        // ⚠⚠ P2-2f: PAD CLEARANCE OWES ITS OWN LEDGER TERM, AND IT IS **NOT**
        // ledger_E_flake. This is the THIRD instance of the class that regime 5
        // (:2557) and regime 6 (:2571) already record: a path that moves a cell
        // out of the SOLID set without crediting the sink runs ledger_err_solid
        // negative, here to -3.6e-02 by 500 s on the Meier feet leg, tracking
        // foot_mech_vol exactly. Neither Kant nor CTRL5 has feet, which is why
        // it never appeared before.
        //
        // WHERE THE ENERGY GOES, and why it is a separate counter: this removal
        // is MECHANICAL -- the burner body shears pad-height rock out of the
        // way, it is not heat-driven spallation. The specific-energy price is
        // E / V over the HEAT-DRIVEN paths only (floor regime 1 + wall regime
        // 6), so folding this into ledger_E_flake would put pad energy in the
        // numerator while the denominator excludes pad volume -- a ratio scaled
        // on one side only. foot_mech_E is therefore its own column, credited
        // to ledger_err_solid alongside ledger_E_flake but reported apart from
        // it, so a price can exclude it by construction rather than by
        // subtraction after the fact.
        //
        // ⚠⚠ BOOKING ONLY -- THE FIELD IS NOT TOUCHED, and that is forced, not
        // chosen. Void H is frozen and still counted by EnergyLedger's TOTAL
        // sum, so the enthalpy never leaves the domain: jet_d2c case CLD fires
        // this path on 620 steps and asserts |ledger_err| <= 1e-10, measuring
        // 4.3e-15. Zeroing H here would destroy enthalpy the total ledger is
        // counting and break that passing assertion. The three sinks above
        // (:2790, :2877, :2931) all book H*Vcell without writing H, for the
        // same reason; this one mirrors them exactly.
        const bool mech_on = surface_patch.enabled && surface_patch.feet &&
                             surface_patch.foot_body_clearance;
        Set::Scalar E_mech_local = 0.0;
        amrex::Vector<int> mech_n(mech_on ? ncols : 0, 0);
        amrex::Vector<int> mech_ktop(mech_on ? ncols : 0, std::numeric_limits<int>::min());
        amrex::Vector<int> mech_flip(mech_on ? ncols : 0, 0);
        if (mech_on)
        {
            if (!SurfaceOwner(lev))
                Util::Abort(INFO, "surface_patch.foot_body_clearance runs on the surface-owning level only");
            foot_mech_cols = 0.0;
            if (foot_clear_k.size() == static_cast<std::size_t>(ncols))
            {
                const amrex::Vector<int> kt = ColumnTopSolid(lev);
                for (int idx = 0; idx < ncols; ++idx)
                {
                    const int kc = foot_clear_k[idx];
                    if (kc == std::numeric_limits<int>::min() || kt[idx] <= kc) continue;
                    mech_n[idx]    = kt[idx] - kc;
                    mech_ktop[idx] = kt[idx];
                }
                for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    amrex::Array4<Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& reg = (*regime_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& D   = (*D_mf[lev]).array(mfi);
                    amrex::Array4<const Set::Scalar> const& Hc =
                        (*H_mf[lev]).const_array(mfi);
                    amrex::Array4<Set::Scalar> regf;
                    if (sp_weibull_spall) regf = (*regime_field_mf[lev]).array(mfi);
                    const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                    for (int j = lo.y; j <= hi.y; ++j)
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        const int idx = col_index(i, j);
                        const int n = mech_n[idx];
                        if (n == 0) continue;
                        const int kc = foot_clear_k[idx];
                        const Set::Scalar h = Set::Scalar(n) * dz;
                        for (int k = lo.z; k <= hi.z; ++k)
                        {
                            phi(i, j, k) -= h;
                            if (k > kc && rem(i, j, k) <= 0.5)
                            {
                                rem(i, j, k) = 1.0;
                                reg(i, j, k) = 4.0;
                                if (sp_weibull_spall) regf(i, j, k) = 4.0;
                                D(i, j, k)   = 0.0;
                                if (phi(i, j, k) >= 0.0)
                                    phi(i, j, k) = -std::numeric_limits<Set::Scalar>::min();
                                // per-rank: each rank owns distinct cells here,
                                // unlike mech_flip's per-column tally. Reduced
                                // once below, with it.
                                E_mech_local += Hc(i, j, k) * DX[0] * DX[1] * dz;
                                ++mech_flip[idx];
                            }
                            else if (k == kc && phi(i, j, k) < 0.0)
                                phi(i, j, k) = 0.0;
                        }
                    }
                }
                amrex::ParallelDescriptor::ReduceIntSum(mech_flip.data(), ncols);
                amrex::ParallelDescriptor::ReduceRealSum(E_mech_local);
                foot_mech_E += E_mech_local;
                long nflip = 0;
                for (int idx = 0; idx < ncols; ++idx)
                {
                    nflip += mech_flip[idx];
                    if (mech_flip[idx] > 0) foot_mech_cols += 1.0;
                }
                foot_mech_vol += Set::Scalar(nflip) * DX[0] * DX[1] * dz;
            }
            foot_clear_k.clear();   // consumed; FeetDescent re-arms it
        }

        // PASS 1: each rank fills its local maximum k_top per column,
        // sentinel INT_MIN for columns it does not own at all. ReduceIntMax
        // collapses to the global topmost material cell per (i, j). After
        // the reduce, every rank knows the same k_top_global.
        amrex::Vector<int> k_top_global(ncols, std::numeric_limits<int>::min());
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                int k_top_local = std::numeric_limits<int>::min();
                for (int k = khi; k >= klo; --k)
                {
                    // D2q-1: removed_mf is authoritative (see ColumnTopSolid).
                    if (rem(i, j, k) <= 0.5) { k_top_local = k; break; }
                }
                if (k_top_local > std::numeric_limits<int>::min())
                {
                    int idx = col_index(i, j);
                    if (k_top_local > k_top_global[idx]) k_top_global[idx] = k_top_local;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(k_top_global.data(), ncols);

        // Step 20: per-column cos(theta) of the cavity-wall inclination for
        // the surface-normal spall frame (computed from the gathered tops,
        // so identical on every rank). nullptr when the flag is off -- every
        // conversion below then multiplies by exactly 1.0 (bit-identical).
        amrex::Vector<Set::Scalar> sn_cos_col;
        if (spall_surface_normal && sp_weibull_spall)
            ColumnSurfaceCos(k_top_global.data(), nx_dom, ny_dom,
                             DX[0], DX[1], dz, domain.smallEnd(2), sn_cos_col);
        const Set::Scalar* sn_cos_p =
            sn_cos_col.empty() ? nullptr : sn_cos_col.data();

        // PASS 1.5: per-column max D over `spall.sample`'s scan window
        // (top_cell | column_max | shallow_max). MPI z-decomposition is
        // collapsed by ReduceRealMax; empty columns get 0.0 (no-fire).
        amrex::Vector<Set::Scalar> D_max_col(ncols, 0.0);
        if (spall_on)
        {
            const SpallSampleMode sample_mode = spall_sample_mode;
            const Set::Scalar shallow_d       = spall_shallow_depth;
            const Set::Scalar dz_local        = DX[2];
            const int k_band_shallow =
                static_cast<int>(std::floor(shallow_d / dz_local));

            for (amrex::MFIter mfi(*D_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& Dfab =
                    (*D_mf[lev]).array(mfi);
                const int klo_box = bx.smallEnd(2), khi_box = bx.bigEnd(2);
                const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
                const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);

                for (int j = jlo; j <= jhi; ++j)
                for (int i = ilo; i <= ihi; ++i)
                {
                    int idx = col_index(i, j);
                    int k_top = k_top_global[idx];
                    if (k_top == std::numeric_limits<int>::min()) continue;

                    int k_hi_scan, k_lo_scan;
                    switch (sample_mode)
                    {
                        case SpallSampleMode::TopCell:
                            // Only the rank that owns k_top contributes;
                            // every other rank skips so the reduction
                            // picks up exactly D(i, j, k_top).
                            if (k_top < klo_box || k_top > khi_box) continue;
                            k_hi_scan = k_top;
                            k_lo_scan = k_top;
                            break;
                        case SpallSampleMode::ColumnMax:
                            k_hi_scan = std::min(k_top, khi_box);
                            if (k_hi_scan < klo_box) continue;
                            k_lo_scan = klo_box;
                            break;
                        case SpallSampleMode::ShallowMax:
                            k_hi_scan = std::min(k_top, khi_box);
                            if (k_hi_scan < klo_box) continue;
                            k_lo_scan = std::max(klo_box,
                                                 k_top - k_band_shallow);
                            if (k_lo_scan > k_hi_scan) continue;
                            break;
                        default:
                            continue;  // unreachable
                    }

                    Set::Scalar local_max = 0.0;
                    for (int k = k_lo_scan; k <= k_hi_scan; ++k)
                    {
                        Set::Scalar v = Dfab(i, j, k);
                        if (!std::isfinite(v)) continue;
                        if (v > local_max) local_max = v;
                    }
                    if (local_max > D_max_col[idx]) D_max_col[idx] = local_max;
                }
            }
            amrex::ParallelDescriptor::ReduceRealMax(D_max_col.data(), ncols);
        }

        // PASS 2: only the rank that owns (i, j, k_top_global) writes
        // into the per-column candidate buffers; others contribute 0.
        // ReduceRealSum / ReduceIntSum collapse the per-rank
        // contributions onto every rank.
        amrex::Vector<int> spall_fires_col(ncols, 0);
        amrex::Vector<int> vapor_fires_col(ncols, 0);
        amrex::Vector<Set::Scalar> h_s_cand_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> h_v_cand_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> rop_s_cand_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> rop_v_cand_col(ncols, 0.0);

        // ------------------------------------------------------------------
        // P1b: the support rule's GEOMETRY, computed once, before the scan.
        //
        // n_exposed = # of the 4 lateral neighbours whose top solid cell sits
        // at least spall.support_step_min BELOW this column's. Read from
        // k_top_global -- the gathered, ReduceIntMax'd surface -- and NEVER
        // from a recession sum: k_top is the authoritative void geometry
        // (D2q-1), and a recession sum is both the wrong quantity and
        // definition-sensitive (the census measured 69/73/54/23 % for the live
        // fraction depending on which was used, against a stable 78 % for the
        // frozen one). Because the input is already reduced, every rank
        // computes the identical array and the count is decomposition-
        // independent by construction -- the P0b (e) failure cannot recur here.
        //
        // A domain-edge neighbour is NOT exposure. On the quarter domain the
        // x = 0 and y = 0 faces are symmetry planes, so the absent neighbour is
        // this column's own mirror image and its step is exactly zero; the far
        // faces are undrilled rock, which is shallower, not deeper. Skipping is
        // correct in both cases, and it is the same guard the flake coherence
        // cap already uses below.
        // Step 16: capture sp_weibull-specific K_Ic + scan parameters
        // outside the MFIter so the lambda in the depth scan can close over
        // them. These values are unused when spallation_model == DamageLaw.
        const KIcLaw      law_loc          = kic_law;
        // ⚠ P2-2d d7: the law object for `grimsel_table`. Hoisted once -- the
        // three sites below must not each decide where the toughness comes
        // from, which is how the Pa-vs-MPa unit bug and the h_side cap bug
        // both arose in this file.
        const Numeric::Material::Grimsel *kic_g_loc =
            phases.empty() ? nullptr : &phases[0].grimsel;
        const Set::Scalar kic_scale_loc    = kic_scale;
        const Set::Scalar kic_K_Ic0_loc    = kic_K_Ic0;
        const Set::Scalar kic_T_ref_loc    = kic_T_ref;
        const Set::Scalar kic_T_melt_loc   = kic_T_melt;
        const Set::Scalar h_spall_max_loc  = spall_h_spall_max;
        const int         h_spall_n_seg_loc =
            (spall_h_spall_n_segments > 0) ? spall_h_spall_n_segments
                                           : sp_n_segments;

        // Step 16d: capture phase 0 material params + the prescribed-T
        // face callable for the T-derived sigma_xx reconstruction
        // (option (a) in rossi-validation-diagnostic-design.md §8).
        // The bisection sample uses
        //   sigma_xx(z) = -E * beta * (T(z) - T_ref) / (1 - nu)
        // with T(z) linearly interpolated between sp_T_face_loc(time) at
        // z = 0 and cell-centred temp_mf values at z = (n + 0.5) * dz.
        // Single-phase Voronoi only — multi-phase requires per-cell phase
        // lookup which is out of scope per the active step's Guardrails.
        // Fail fast if surface_patch isn't in prescribed_T mode (the face
        // callable would be invalid).
        Set::Scalar phase_E_loc     = 0.0;
        Set::Scalar phase_beta_loc  = 0.0;
        Set::Scalar phase_nu_loc    = 0.0;
        Set::Scalar phase_T_ref_loc = 0.0;
        Set::Scalar phase_kappa_loc = 0.0;   // P1 witness: the q_fit datum
        auto sp_T_face_loc          = surface_patch.T_f;  // ParserExecutor<1>
        if (sp_weibull_spall)
        {
            if (phases.empty())
                Util::Abort(INFO, "Step 16d (sp_weibull + spall.enabled): "
                                  "phases vector is empty; need at least "
                                  "phase 0 for T-derived sigma_xx");
            phase_E_loc    = phases[0].E;
            phase_beta_loc = phases[0].beta;
            if (!(phases[0].E > 0.0 && phases[0].mu > 0.0))
                Util::Abort(INFO, "Step 16d: phase 0 must have E > 0 and "
                                  "mu > 0 for the T-derived sigma_xx path; "
                                  "got E=", phases[0].E, ", mu=",
                                  phases[0].mu);
            phase_nu_loc    = phases[0].E / (2.0 * phases[0].mu) - 1.0;
            phase_T_ref_loc = phases[0].T_ref;
            phase_kappa_loc = phases[0].kappa;
            if (!(phase_nu_loc > -1.0 && phase_nu_loc < 0.5))
                Util::Abort(INFO, "Step 16d: phase 0 Poisson's ratio out "
                                  "of physical range; got nu=", phase_nu_loc);
        }

        // ================= P1f: THE LATERAL CRITERION (deliverable 2) =========
        // Per (cell, FACE), not per column.
        //
        // WHY A FACE AND NOT A COLUMN -- a theorem, not a measurement: after P1e's
        // staircase correction side heat is isotropic per unit of wall ARC, but a
        // 45 deg wall carries the same arc with sqrt(2) FEWER columns, so a diagonal
        // wall column still carries 1.4142 the face-equivalents of an axis column.
        // ⇒ no per-column 1-D criterion can be azimuthally isotropic on a round
        // wall however the area is corrected, which is why P1e's (f) got WORSE.
        //
        // ⚠⚠ WHY HERE AND NOT IN THE KERNEL: MMWSpalling.H:2070-2084 enumerates the
        // same (cell, face) pairs, but it sits inside amrex::ParallelFor with a [=]
        // capture -- a device-capable lambda, where SpCriterion::K_I's callable and
        // these per-column vectors cannot go. Removal.H's MFIter is a plain HOST
        // loop; that is why the vertical scan closes over locals by reference.
        amrex::Vector<Set::Scalar> lat_score, lat_ell, lat_Tf, lat_depth;
        // ⚠⚠ P1f deliverable 10: THE QUANTITIES THE BRANCH ACTUALLY DECIDES ON.
        // They were never output, so the Pa-vs-MPa unit error at the comparison
        // was invisible to every bar in the packet: 224,359 firings on rock a
        // factor of 9 BELOW threshold, while (h4) went green and (r) read a
        // BETTER isotropy. A branch's deciding quantity must be loggable in
        // stated units. lat_KI/lat_KIc are Pa*sqrt(m); lat_sig is Pa.
        // ⚠ P1g-0 deliverable 2: lat_hf is the flake thickness the lateral
        // stress was referenced to, in METRES. Its median is 0 today, so the
        // quantity that enters sigma under `differential` -- P1g's form -- was
        // not an output at all. It is also what settles the cap share: P1f
        // quoted "99.8 % vs 61.9 %" from two estimators on two legs, and
        // h_f AS LOGGED on regime-6 rows is at the cap on only 4.9 %.
        amrex::Vector<Set::Scalar> lat_KI, lat_KIc, lat_q, lat_sig, lat_hf;
        // P1g-1 deliverable 7: lat_hfu is the UNCAPPED c_f*ell -- the
        // thickness the criterion earned before spall.flake_h_max clipped
        // it. It is logged so a later step can re-derive the bound FROM
        // EVIDENCE rather than from the cap that currently sets it.
        // lat_acc_log / lat_nf_log are the accumulator and that cell's own
        // cumulative firing count, both read back AFTER banking.
        amrex::Vector<Set::Scalar> lat_hfu, lat_acc_log, lat_nf_log;
        // ⚠⚠ P1g-2 deliverable 3: THE FACE NORMAL, CARRIED INTERNALLY AND NEVER
        // LOGGED. `lat_csv = sp_lateral_drive` (:1012) -- the DRIVE key alone --
        // so adding a CSV column would move the header of every lateral leg
        // INCLUDING GEOM's, and (a2) (GEOM byte-identical to P1g-1's archive)
        // would fail by construction rather than by physics. The shift needs the
        // normal in-register, so it stays here.
        // ⚠ `lat_sx/sy` are encoded as sign+2 (so 1, 2 or 3, never negative)
        // because the per-column winner is selected with ReduceRealMax and a
        // negative payload would lose to a non-winner's zero.
        amrex::Vector<Set::Scalar> lat_sx_enc, lat_sy_enc, lat_wx, lat_wy;
        amrex::Vector<int>         lat_kfire;
        long lat_no_faceT = 0, lat_cand = 0, lat_skip_hside = 0;
        // P1g-3: the LATERAL site's own two LEN fallback counters.
        long lat_len_fb_h0 = 0, lat_len_fb_q0 = 0;
        // ⚠⚠ P2-2e: the constraint depth at BOTH Removal.H sites (the lateral
        // criterion and the depth scan). `sp_cref_d > 0` alone is the switch:
        // Parse aborts unless d > 0 comes with the T-derived criterion path, and
        // `sp_weibull_use_T_face` is not formed until much later in this
        // function -- AFTER the lateral criterion, which must read it too.
        // Witnesses are accumulated per BOX (the scan runs a column only on the
        // rank owning its top, :1497; the lateral loop is per box) and reduced
        // ONCE below the scan.
        const bool        cref_on    = sp_cref_d > 0.0;
        const Set::Scalar cref_d_loc = sp_cref_d;
        long cref_vert = 0, cref_lat = 0, cref_clamp = 0, cref_neg = 0;
        // ⚠⚠ deliverable 12: taken over EVERY candidate, before the firing
        // branch -- the sub-threshold score (w)'s upper bound has to read.
        Set::Scalar lat_qmax_l = 0.0, lat_Tfq_l = 0.0;
        const bool lat_on = sp_lateral_drive && sp_weibull_spall;
        // P1g-1: the banking law's reachability, as the ledger states it. Three
        // deep on purpose -- subcell alone, or subcell without removal, is inert
        // rather than an error, because the NEG leg scores exactly that.
        const bool bank_on = lat_on && sp_lateral_removal && sp_lateral_subcell;
        // ⚠⚠ P1g-0b deliverable 3 -- THE LATENT ROWS, CLOSED AT RUNTIME.
        // The criterion hard-codes two SurfaceCellFlux arguments that the kernel
        // computes. Each is equal only under a condition nobody asserted, so each
        // is asserted here. These forbid configurations, they are not new physics,
        // and both are silent on every run in the programme.
        if (lat_on)
        {
            // `face_form`: the criterion passes a hard-coded `true`; the kernel
            // passes `sp_surfform ? true : sp_face` (MMWSpalling.H:2075).
            // ⚠⚠ `sp_surfform` is TWO conditions (MMWSpalling.H:1825-1828) --
            // robin_form == Surface AND side_face_closure == follow -- while the
            // criterion's own scope test at :544 checks only the FIRST. So with
            // `surface_patch.side_face_closure = cell` the criterion would run
            // with face_form = true against a kernel using face_form = false:
            // a live mismatch of the same family, and NOT covered by the scope
            // test. Both conditions are required here.
            const bool surfform_eq =
                surface_patch.robin_form == SurfacePatch::RobinForm::Surface
             && surface_patch.side_T_form == SurfacePatch::SideTForm::Follow;
            if (!surfform_eq)
                Util::Abort(INFO, "spall.lateral_drive requires "
                            "surface_patch.robin_form = surface AND "
                            "surface_patch.side_face_closure = follow: the lateral "
                            "criterion hard-codes face_form = true, which equals "
                            "the kernel's `sp_surfform ? true : sp_face` only "
                            "under BOTH. See MMWSpalling.H:1825-1828 and :2075.");
            // `T_pin`: the criterion passes a hard-coded 0.0; the kernel passes
            // `tpe ? tpe[col_c] : 0.0`, and `tpe` is non-null only when
            // `sp_pinned` (MMWSpalling.H:1830).
            // ⚠ This cannot fire independently of the abort above -- sp_pinned
            // means robin_form == Pinned, which already fails surfform_eq -- so it
            // is defence in depth against that coupling being broken, not a
            // second reachable configuration. Kept because the cost is one
            // comparison and the coupling is not obvious from either site.
            if (surface_patch.robin_form == SurfacePatch::RobinForm::Pinned)
                Util::Abort(INFO, "spall.lateral_drive is incompatible with "
                            "surface_patch.robin_form = pinned: the lateral "
                            "criterion hard-codes T_pin = 0.0 while the kernel "
                            "passes the thermostat array (MMWSpalling.H:1830).");
        }
        if (lat_on)
        {
            lat_score.assign(ncols, 0.0);
            lat_ell  .assign(ncols, 0.0);
            lat_Tf   .assign(ncols, 0.0);
            lat_depth.assign(ncols, 0.0);
            lat_KI   .assign(ncols, 0.0);
            lat_KIc  .assign(ncols, 0.0);
            lat_q    .assign(ncols, 0.0);
            lat_sig  .assign(ncols, 0.0);
            lat_hf   .assign(ncols, 0.0);
            lat_hfu    .assign(ncols, 0.0);
            lat_sx_enc .assign(ncols, 0.0);
            lat_sy_enc .assign(ncols, 0.0);
            lat_wx     .assign(ncols, 0.0);
            lat_wy     .assign(ncols, 0.0);
            lat_acc_log.assign(ncols, 0.0);
            lat_nf_log .assign(ncols, 0.0);
            lat_kfire.assign(ncols, std::numeric_limits<int>::min());
            // ⚠ Guards asserted, not assumed (deliverable 2).
            const bool lat_face_form =
                (surface_patch.robin_form == SurfacePatch::RobinForm::Surface);
            const Set::Scalar rmin2 = surface_patch.side_face_h_r_min
                                    * surface_patch.side_face_h_r_min;
            // ⚠⚠ P1g-0b deliverable 1 -- THE COEFFICIENT CAP, now resolved
            // through the KERNEL'S OWN CONDITION. This line read
            //     side_hmax = surface_patch.side_face_h_max;
            // i.e. it capped UNCONDITIONALLY at the 80.0 W/m^2K default, while
            // the kernel (MMWSpalling.H:1768-1770) and the ledger (:3937-3939)
            // cap only under `side_face_flux && spall_per_cell_removal` -- a
            // conjunction that is FALSE on every live run, because
            // `spall.per_cell_removal` defaults to 0 (MMWSpalling.H:5964) and no
            // live argv sets it. So the criterion judged with h = 80 while the
            // kernel heated the same face with 169.3: a face tens of kelvin
            // COOLER than the one that heated the rock, and a lift well under
            // half what the flux produced.
            //
            // ⚠⚠ THIS IS THE THIRD INSTANCE OF ONE BUG FAMILY IN THIS ONE CALL:
            // Pa vs MPa (P1f d9), the Newton LENGTH d_n vs DX[0] (P1g-0), and
            // now the COEFFICIENT. Each time the criterion read something the
            // flux never used. The argument-list lock in
            // tests/MMWSpalling/unit/lateral_flux_lock/test exists so the fourth
            // is caught by a test rather than by a campaign.
            //
            // ⚠ -1.0 is how the kernel and the ledger spell "no cap": the
            // application site below tests `side_hmax > 0.0`, so -1.0 disables
            // it. Nothing there changes.
            const Set::Scalar side_hmax  =
                (surface_patch.side_face_flux && spall_per_cell_removal)
                ? surface_patch.side_face_h_max : Set::Scalar(-1.0);
            const Set::Scalar side_h_ex  = surface_patch.side_face_h;
            // P2-3 site 2 of 5 (the LATERAL CRITERION). ⚠⚠ THIS IS THE SITE
            // THAT GETS MISSED (P2-0 lesson): it is the earliest of the five, it
            // lives in a different file from the other four, and a criterion
            // reading a different gas temperature from the kernel decides
            // whether a face fires against a flux the kernel never applied.
            // Resolved PER FACE at :678, not hoisted, because the march's whole
            // point is that the value depends on the face's height.
            const Set::Scalar side_T_ex  = side_T_eff;   // fallback only
            // P2-2 change 4, site 4 of 7 (the LATERAL criterion).
            const Set::Scalar relief_loc = sp_relief;
            const Set::Scalar C_th = relief_loc * phase_E_loc * phase_beta_loc
                                   / (1.0 - phase_nu_loc);
            // P2-2b deliverable 4: the rock law, or nullptr under `constant`.
            const Numeric::Material::Grimsel *lat_g = RockLaw();
            const Set::Scalar T0 = phase_T_ref_loc;
            const Set::Scalar sp_x0 = surface_patch.x0, sp_y0 = surface_patch.y0;
            for (amrex::MFIter mfi(*temp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& T  = (*temp_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& rm = (*removed_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> fa;
                if (weibull_enabled) fa = (*flaw_a_mf[lev]).const_array(mfi);
                // ⚠⚠ P1g-0b deliverable 3: THE PER-CELL CONDUCTIVITY, read the
                // way the kernel reads it (MMWSpalling.H:1840, :1854). This
                // replaced the SCALAR `phases[0].kappa`, which agreed with the
                // kernel only while the rock was single-phase AND
                // `conductivity.kappa_damage_alpha` was 0 -- two configurations
                // nobody asserted. Reading the field retires BOTH at once,
                // instead of forbidding each in turn.
                // ⚠ The kernel also zeroes kc on a void cell; here the cell is
                // solid by construction (`rm(i,j,k) > 0.5` is skipped below), so
                // that branch cannot be reached and is deliberately absent.
                amrex::Array4<const Set::Scalar> const& kpa =
                    (*kappa_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& Dfa =
                    (*D_mf[lev]).const_array(mfi);
                const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
                const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
                const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
                for (int j = jlo; j <= jhi; ++j)
                for (int i = ilo; i <= ihi; ++i)
                {
                    const int idx = col_index(i, j);
                    const int kT  = k_top_global[idx];
                    if (kT == std::numeric_limits<int>::min()) continue;
                    const Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    const Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    const Set::Scalar r2 = (xc - sp_x0) * (xc - sp_x0)
                                         + (yc - sp_y0) * (yc - sp_y0);
                    // ⚠⚠ SCOPE: the 7a HOLE-WALL branch only. A lateral firing on a
                    // one-cell floor riser would substitute whole-cell removal (dz)
                    // for the vertical flake path (<= flake_h_max), cutting the
                    // floor rate to 0.75 per event for no physical reason. (h4) is
                    // the witness. ⚠ Deliverable 1 measured that this cut also
                    // excludes 88.9 % of P1e's viable faces -- which is why Leg A
                    // is inert BY SCOPE and only Leg B is launched.
                    if (!(lat_face_form && r2 >= rmin2)) continue;
                    // the staircase riser heights -- the same inputs P1e uses
                    const int nb[4][2] = {{i-1,j},{i+1,j},{i,j-1},{i,j+1}};
                    int rn[4] = {0,0,0,0};
                    for (int d = 0; d < 4; ++d)
                    {
                        const int ni = nb[d][0], nj = nb[d][1];
                        if (ni < dom_ilo || ni > dom_ihi ||
                            nj < dom_jlo || nj > dom_jhi) continue;
                        const int kS = k_top_global[col_index(ni, nj)];
                        if (kS == std::numeric_limits<int>::min()) continue;
                        const int n = kT - kS;
                        if (n >= 1) rn[d] = n;
                    }
                    const Set::Scalar hx = Set::Scalar(rn[0] + rn[1]);
                    const Set::Scalar hy = Set::Scalar(rn[2] + rn[3]);
                    const Set::Scalar hyp = std::hypot(hx, hy);
                    if (!(hyp > 0.0)) continue;
                    // side-face h and T, resolved exactly as the kernel does
                    Set::Scalar h_side = (idx < static_cast<int>(flame_h.size()))
                                       ? flame_h[idx] : Set::Scalar(0.0);
                    Set::Scalar T_side = (idx < static_cast<int>(flame_Tg.size()))
                                       ? flame_Tg[idx] : Set::Scalar(0.0);
                    if (side_h_ex >= 0.0) { h_side = side_h_ex; T_side = side_T_ex; }
                    if (side_hmax > 0.0 && h_side > side_hmax) h_side = side_hmax;
                    if (!(h_side > 0.0)) { ++lat_skip_hside; continue; }
                    for (int k = klo; k <= khi && k <= kT; ++k)
                    {
                        if (rm(i, j, k) > 0.5) continue;          // this cell is void
                        // P2-1 deliverable 4, ROLE 4 of 4 (the CRITERION).
                        // ⚠⚠ RESOLVED HERE, INSIDE THE k LOOP, because the law is
                        // a function of HEIGHT and the column-level block above
                        // has no k. Leaving the criterion on the column scalar
                        // while the kernel went per-height would have been the
                        // FOURTH instance of the family this call pre-registers
                        // at :589-594 -- the criterion reading something the flux
                        // never used. unit/lateral_flux_lock's P2-1 block is what
                        // makes that impossible rather than merely unlikely.
                        // ⚠ The cap is re-applied to the resolved value, through
                        // the same `side_hmax` conjunction, exactly as the kernel
                        // and the ledger do.
                        Set::Scalar h_face = SideFaceH(h_side, r2, idx, k);
                        if (side_hmax > 0.0 && h_face > side_hmax) h_face = side_hmax;
                        if (!(h_face > 0.0)) continue;
                        // ⚠⚠ P2-3 site 2 of 5 (the LATERAL CRITERION) -- AND IT
                        // HAD TO MOVE. `h_side`/`T_side` above are resolved ONCE
                        // PER COLUMN, outside this k loop, which is correct for
                        // a single scalar but wrong the moment the gas
                        // temperature depends on height. `h_face` was already
                        // per-k (SideFaceH takes k); the gas temperature now is
                        // too, resolved on the same line of reasoning and in the
                        // same loop, so the criterion and the kernel cannot
                        // disagree about which gas a face is standing in.
                        const Set::Scalar T_face_gas = (side_h_ex >= 0.0)
                            ? SideFaceTgas(idx, k) : T_side;
                        for (int d = 0; d < 4; ++d)
                        {
                            if (rn[d] == 0) continue;
                            const int ni = nb[d][0], nj = nb[d][1];
                            // the face is open AT THIS k only if the neighbour cell
                            // at the same k is void
                            if (rm(ni, nj, k) <= 0.5) continue;
                            // the solid cell BEHIND the face, opposite the void
                            const int bi = i - (ni - i), bj = j - (nj - j);
                            if (bi < dom_ilo || bi > dom_ihi ||
                                bj < dom_jlo || bj > dom_jhi) continue;
                            if (rm(bi, bj, k) > 0.5) continue;
                            // ⚠⚠ THE NORMAL DISTANCE, NOT dx. A plain dT/dx
                            // under-reads the wall-normal gradient by 1/|n_x| =
                            // 1.4142 at 45 deg -- a 41 % artefact that would fail
                            // (c) with nothing wrong with the physics. ⚠ NOT the
                            // staircase factor f = 1/(|n_x|+|n_y|): that equals
                            // |n_x| at 0 and 45 deg ONLY -- the two angles (c)
                            // tests -- and is 15.5 % off at 15 and 30 deg.
                            const Set::Scalar comp = (d < 2) ? std::abs(hx / hyp)
                                                             : std::abs(hy / hyp);
                            if (!(comp > 1.0e-12)) continue;
                            const Set::Scalar d_n = DX[0] * comp;
                            // ---- the FACE temperature, from the flux itself ----
                            // ⚠⚠ deliverable 0's out-parameter. NOT an
                            // extrapolation of the vertical half-cell formula: on a
                            // 45 deg apex the surface plane is ambiguous by 0-0.707
                            // cells, worth x2.03 at ell = 1 cell. Flux and criterion
                            // must read ONE temperature, which is P0b's rule
                            // sideways -- and until P1g-0 they did not.
                            //
                            // ⚠⚠ P1g-0: TWO LENGTHS, TWO QUESTIONS. They are not
                            // interchangeable and this site had them crossed.
                            //
                            //   the NEWTON's length, HERE -> DX[0].
                            //     G = 2k/dz is the conduction half-distance from
                            //     the cell CENTRE to the cell FACE. A cube face is
                            //     DX[0]/2 away whatever the wall's inclination, so
                            //     the length is DX[0] and carries NO |n_x|. It must
                            //     equal what the KERNEL passes, because the face
                            //     this Newton describes is the face the kernel's
                            //     flux heated: MMWSpalling.H:2082, :3192 and :4066
                            //     all pass DX[0]. Passing d_n here made G 41 %
                            //     larger on a 45 deg corner, pulling the face
                            //     towards the cell, so the criterion read a face
                            //     ~9 K COOLER THAN THE FLUX THAT HEATED IT -- and
                            //     a corner then needed its cell 9 K hotter than an
                            //     axis face to reach the same threshold, which the
                            //     vertical mode never allowed. Measured on P1f's
                            //     Leg B: regime 6 fired 232 times on the axes and
                            //     6 on the diagonals below 150 mm, while the
                            //     VERTICAL mode was isotropic to 1-2 % in the same
                            //     bands -- so it was never a flux anisotropy.
                            //   the SLOPE's length, at :611 -> d_n.
                            //     ell = d_n/log(dT_c/dT_b) IS wall-normal: the two
                            //     cells are d_n apart along the normal, and using
                            //     DX[0] there would under-read the gradient by
                            //     1/|n_x| = 1.4142 -- the 41 % artefact the comment
                            //     above warns about. d_n STAYS there.
                            //
                            // ⚠ Whether the staircase should carry inclination in
                            // the Newton as WELL as in the area factor is a P1g
                            // question. What is settled here is only the
                            // CONSISTENCY claim: criterion and kernel read one
                            // temperature.
                            Set::Scalar s_in = 0, s_out = 0, s_rob = 0, s_los = 0;
                            Set::Scalar T_face = -1.0;      // sentinel
                            // P1g-0b deliverable 3: the conductivity the KERNEL
                            // would use on this same cell, expression for
                            // expression (MMWSpalling.H:1854).
                            const Set::Scalar kc_cell =
                                kappa_damage_alpha != 0.0
                                ? kpa(i, j, k)
                                  * std::exp(-kappa_damage_alpha * Dfa(i, j, k))
                                : kpa(i, j, k);
                            if (!(kc_cell > 0.0)) continue;
                            SurfaceCellFlux(true, true, 1.0, h_face, T_face_gas,
                                            kc_cell, DX[0], T(i, j, k),
                                            eps_rad, eps_rad_hi, T_eps_break,
                                            sigma_SB, h_conv, T_amb,
                                            s_in, s_out, s_rob, s_los,
                                            0.0, nullptr, 0.0, &T_face);
                            if (!(T_face > 0.0)) { ++lat_no_faceT; continue; }
                            // ---- the SLOPE, from the two cells ----
                            const Set::Scalar dT_c = T(i, j, k)    - T0;
                            const Set::Scalar dT_b = T(bi, bj, k)  - T0;
                            if (!(dT_c > 0.0) || !(dT_b > 0.0) || dT_c <= dT_b)
                                continue;   // no inward decay: soaked or inverted
                            const Set::Scalar ell = d_n / std::log(dT_c / dT_b);
                            if (!(ell > 0.0) || !std::isfinite(ell)) continue;
                            const Set::Scalar dT_s = T_face - T0;
                            if (!(dT_s > 0.0)) continue;
                            // ⚠⚠ ONE FORM, ONE DRIVE FRACTION, ONE FLAW. The same
                            // two-cell exponential and the same c_f as the vertical
                            // mode. A LINEAR lateral profile would give drive
                            // fraction c_f = 0.400 against 1-e^{-c_f} = 0.32968 --
                            // a threshold 17.6 % lower on the same flaw, i.e. a
                            // second calibration by the back door.
                            // ⚠⚠ P1g-3 deliverables 5 + 6 AT THE LATERAL SITE.
                            // Both deliverables say "at BOTH sites" and the first
                            // implementation wired only the vertical one, which
                            // left h_f_lat pinned at the 1.500 mm cap on a Meier
                            // leg that was supposed to be testing a derived wall
                            // thickness -- a false negative on the width, for an
                            // implementation reason rather than a physical one.
                            //
                            // ⚠⚠ THE WALL NEEDS ITS OWN FLUX AND CANNOT BORROW
                            // THE FLOOR'S. `q_ratio` (Removal.H's floor route) is
                            // zero or non-finite on 96.0 % of regime-6 events,
                            // because `flame_h` is 0 outside the jet patch and
                            // above the nozzle plane: for the wall the floor flux
                            // is not merely different, it is ABSENT. So q_net is
                            // built from the side-face pair the criterion already
                            // resolved (h_side, T_side) and the SAME net form the
                            // floor uses -- net of radiation and convection, or
                            // the two sites would not be the same quantity:
                            //   q_net = h_side (T_side - T_face) - loss(T_face)
                            //   ell_q = kc_cell (T_face - T0) / q_net
                            // which is exactly the (a3) offline route, now in the
                            // kernel. ⚠ kc_cell, not a scalar kappa: P1g-0b
                            // retired that for the per-cell field.
                            Set::Scalar ell_use = ell;
                            if (sp_flake_len_flux)
                            {
                                const Set::Scalar ep_l =
                                    (T_face <= T_eps_break) ? eps_rad : eps_rad_hi;
                                const Set::Scalar loss_l =
                                    ep_l * sigma_SB
                                    * (T_face*T_face*T_face*T_face
                                       - T_amb*T_amb*T_amb*T_amb)
                                    + h_conv * (T_face - T_amb);
                                const Set::Scalar q_net_l =
                                    h_side * (T_side - T_face) - loss_l;
                                // ⚠ the fallback is to the two-cell length, never
                                // to a clamp, and each reason is counted
                                // separately so an inert key cannot hide.
                                if (!(h_side > 0.0))          ++lat_len_fb_h0;
                                else if (!(q_net_l > 0.0)
                                         || !(T_face > T0))   ++lat_len_fb_q0;
                                else ell_use = kc_cell * (T_face - T0) / q_net_l;
                            }
                            // ⚠ w/m replaces c_f, with w the cone-weighted MEDIAN
                            // (see W_CONE_MEDIAN's provenance note).
                            const Set::Scalar frac_l = sp_flake_frac_weibull
                                ? (W_CONE_MEDIAN / weibull_m) : sp_flake_c_f;
                            Set::Scalar h_f = frac_l * ell_use;
                            if (sp_flake_h_max > 0.0 && h_f > sp_flake_h_max)
                                h_f = sp_flake_h_max;
                            if (!(h_f > 0.0)) continue;
                            // P2-2 change 3, READ 2 of 2 (the LATERAL flaw).
                            // ⚠⚠ ONE factor across both reads -- a size term on
                            // the vertical flaw and not the lateral one is the
                            // same wall-vs-floor split `Removal.H:589-594`
                            // pre-registers, and it would show up as an
                            // anisotropy nobody could attribute.
                            // ⚠ `ell` here is the LATERAL flux length already in
                            // scope, the same quantity surf_ellq_col holds for a
                            // column -- so both reads scale on ell_q^3.
                            Set::Scalar a_f = (weibull_enabled && fa)
                                            ? fa(i, j, k) : sp_a0;
                            if (sp_size_flux && ell > 0.0)
                                a_f *= WeibullSizeFactor(ell * ell * ell);
                            if (!(std::isfinite(a_f) && a_f > 0.0)) continue;
                            // sigma across the flake at normal depth xi, from the
                            // ONE exponential -- the vertical form, sideways.
                            // ⚠ P2-2a: under the retired `differential`
                            // form the drive was referenced to the
                            // flake's own back face and carried a
                            // (1 - e^{-h_f/ell}) factor. `absolute` has
                            // no reference subtraction, so the factor
                            // must NOT appear -- and with one form left
                            // it cannot be selected back on.
                            // ⚠ h_f never enters the stress here, so the
                            // criterion is ell-blind apart from
                            // exp(-a_f/ell) -- 0.16 % at 8 um / 5 mm.
                            auto sigma_lat = [&] (Set::Scalar xi) -> Set::Scalar {
                                Set::Scalar jump = dT_s * std::exp(-xi / ell);
                                if (!(jump > 0.0)) return Set::Scalar(0.0);
                                // ⚠⚠ P2-2e: the reference continues THIS site's
                                // own exponential d further behind the face,
                                // T_ref(xi) = T0 + dT_s exp(-(xi + d)/ell_use),
                                // with ell_use the FLUX length (or its counted
                                // two-cell fallbacks, lat_len_fb_h0/_q0). At the
                                // face the bracket is dT_s (1 - e^{-d/ell_q}) ~
                                // q d / k: the same minimum-flux law the
                                // vertical field read gives.
                                // ⚠ RESOLVED AGAINST THE CODE, NOT THE CITATION
                                // (standing rule 14): the jump above uses `ell`,
                                // the two-cell SLOPE length (:803), not the flux
                                // length -- `ell_use` feeds only h_f -- although
                                // the comment above a_f calls `ell` the flux
                                // length. The jump is left as it is (the stress
                                // formula does not change); only the reference
                                // uses ell_use, as the packet specifies. On every
                                // leg without spall.flake_length = flux the two
                                // are the SAME number.
                                // ⚠ ell_use > 0 here: h_f = frac_l * ell_use and
                                // a non-positive h_f `continue`d above.
                                if (cref_on)
                                {
                                    const Set::Scalar jref =
                                        dT_s * std::exp(-(xi + cref_d_loc) / ell_use);
                                    ++cref_lat;
                                    if (jump < jref) ++cref_neg;
                                    if (lat_g)
                                        return relief_loc
                                             * LawSigma(lat_g, T0 + jump, T0 + jref);
                                    return C_th * (jump - jref);
                                }
                                // P2-2b deliverable 4: strain AND stiffness at
                                // the LOCAL temperature T0 + jump, not a fixed
                                // coefficient times dT.
                                if (lat_g)
                                    return relief_loc
                                         * LawSigma(lat_g, T0 + jump, T0);
                                // ⚠⚠ the pre-P2-2b expression, VERBATIM. Not
                                // `relief_loc * (E*beta*jump/(1-nu))`: that
                                // re-associates the multiply and is NOT bit
                                // identical, which would break (b1) on the
                                // shipped default.
                                return C_th * jump;
                            };
                            const Set::Scalar K_I = Numeric::SpCriterion::K_I(
                                sigma_lat, a_f, h_spall_n_seg_loc);
                            const Set::Scalar T_tip =
                                T0 + dT_s * std::exp(-a_f / ell);
                            // ⚠⚠ P1f deliverable 9: UNITS. KIcNasseri returns the
                            // Nasseri table DIRECTLY, in MPa*sqrt(m) (1.43 … 0.22),
                            // and kic.K_Ic0 is MPa*sqrt(m) too; K_I is built from
                            // C_th = E*beta/(1-nu) with E in Pa, so it is
                            // Pa*sqrt(m). Comparing them unconverted inflated the
                            // score 1e6x -- 224,359 firings a factor of 9 below
                            // threshold. The vertical path converts at :1333 and
                            // :1357; this was the ONLY site in src/ that did not
                            // (swept: KIcNasseri / KIcLinear / MPa_sqrtm_to_Pa_sqrtm).
                            // ⚠ BOTH branches are MPa -- the table AND the fallback.
                            Set::Scalar K_Ic_MPa = kic_K_Ic0_loc;
                            // P2-2d d7, toughness site 2 of 4 (the LATERAL face).
                            if (law_loc == KIcLaw::GrimselTable && kic_g_loc)
                                K_Ic_MPa = kic_scale_loc
                                         * kic_g_loc->K_Ic_of(T_tip);
                            else if (law_loc == KIcLaw::NasseriTable)
                                K_Ic_MPa = Numeric::SpCriterion::KIcNasseri(
                                               T_tip, kic_scale_loc);
                            if (!(K_Ic_MPa > 0.0)) continue;
                            const Set::Scalar K_Ic =
                                Numeric::SpCriterion::MPa_sqrtm_to_Pa_sqrtm(
                                    K_Ic_MPa);
                            ++lat_cand;
                            const Set::Scalar score = K_I / K_Ic;
                            // deliverable 12: BEFORE the branch, so a leg that
                            // fires nothing still reports how close it came.
                            if (score > lat_qmax_l)
                            { lat_qmax_l = score; lat_Tfq_l = T_face; }
                            // ⚠ An APEX cell offers two candidates with the same
                            // flaw and nearly the same drive: evaluate both, keep
                            // the larger, void the cell exactly ONCE.
                            if (score >= 1.0 && score > lat_score[idx])
                            {
                                lat_score[idx] = score;
                                lat_kfire[idx] = k;
                                lat_ell[idx]   = ell;
                                lat_Tf[idx]    = T_face;
                                lat_depth[idx] = Set::Scalar(kT - k) * dz;
                                // deliverable 10: in stated units -- Pa*sqrt(m),
                                // Pa*sqrt(m), dimensionless, Pa. sigma at the
                                // FACE (xi = 0) is the drive's magnitude.
                                lat_KI[idx]    = K_I;
                                lat_KIc[idx]   = K_Ic;
                                lat_q[idx]     = score;
                                lat_sig[idx]   = sigma_lat(0.0);
                                // P1g-0 deliverable 2: metres, capped value
                                lat_hf[idx]    = h_f;
                                // ⚠ BEFORE the cap, by construction: h_f
                                // above has already been clipped, so the
                                // uncapped value must be recomputed here
                                // rather than read back from h_f.
                                lat_hfu[idx]   = sp_flake_c_f * ell;
                                // ⚠⚠ P1g-2: THE SIGNED INWARD DIRECTION.
                                // nb[] is {-x, +x, -y, +y} and rn[d] > 0 means
                                // the neighbour in direction d sits LOWER, i.e.
                                // the face is exposed TOWARD d. The rock is
                                // therefore OPPOSITE d, so inward = -d.
                                // ⚠ rn[0] and rn[1] both set is a one-cell
                                // ridge exposed on both sides: there is no
                                // inward direction, so the sign is left 0 and
                                // the row-length rule below refuses it.
                                int sx = 0, sy = 0;
                                if (rn[0] > 0 && rn[1] == 0) sx = +1;
                                else if (rn[1] > 0 && rn[0] == 0) sx = -1;
                                if (rn[2] > 0 && rn[3] == 0) sy = +1;
                                else if (rn[3] > 0 && rn[2] == 0) sy = -1;
                                lat_sx_enc[idx] = Set::Scalar(sx + 2);
                                lat_sy_enc[idx] = Set::Scalar(sy + 2);
                                // ⚠⚠ cos SQUARED, and the reason is arithmetic,
                                // not taste. The cos weights hx/hyp and hy/hyp
                                // sum to sqrt(2) at a 45 deg corner, so shifting
                                // BOTH rows by cos*h would displace 1.414*h of
                                // row material against a bank of ONE h -- the
                                // accumulator banks one thickness per firing.
                                // cos^2 weights sum to EXACTLY 1, so displaced
                                // volume is h*A. (b5) scores this ratio at
                                // 1.000, and neither (b1) nor (b2) would see the
                                // error: geometry would still conserve and
                                // energy would still balance while the thermal
                                // half removed sqrt(2) too much.
                                // ⚠ This is P1e's sqrt(2) theorem reappearing as
                                // a geometry-vs-thermal mismatch inside ONE step.
                                const Set::Scalar h2 = hyp * hyp;
                                lat_wx[idx] = (h2 > 0.0) ? (hx * hx) / h2 : 0.0;
                                lat_wy[idx] = (h2 > 0.0) ? (hy * hy) / h2 : 0.0;
                            }
                        }
                    }
                }
            }
            // ⚠ Each cell is owned by exactly one rank, but a column spans ranks in
            // z, so the winner needs a reduction. ReduceRealMax on the score, then
            // ReduceIntMax on the firing k among ranks holding that score -- ties
            // (which require bit-identical scores) resolve to the larger k, so the
            // choice is deterministic rather than decomposition-dependent.
            if (ncols > 0)
            {
                amrex::Vector<Set::Scalar> sc_max = lat_score;
                amrex::ParallelDescriptor::ReduceRealMax(sc_max.data(), ncols);
                for (int c = 0; c < ncols; ++c)
                    if (!(lat_score[c] >= sc_max[c] && sc_max[c] >= 1.0))
                        lat_kfire[c] = std::numeric_limits<int>::min();
                amrex::ParallelDescriptor::ReduceIntMax(lat_kfire.data(), ncols);
                for (int c = 0; c < ncols; ++c)
                {
                    if (lat_kfire[c] == std::numeric_limits<int>::min()
                        || !(lat_score[c] >= sc_max[c]))
                    {
                        lat_ell[c] = 0.0; lat_Tf[c] = 0.0; lat_depth[c] = 0.0;
                        lat_KI[c]  = 0.0; lat_KIc[c] = 0.0;
                        lat_q[c]   = 0.0; lat_sig[c] = 0.0;
                        lat_hf[c]  = 0.0;
                        lat_hfu[c] = 0.0;
                        lat_sx_enc[c] = 0.0; lat_sy_enc[c] = 0.0;
                        lat_wx[c] = 0.0;     lat_wy[c] = 0.0;
                    }
                }
                amrex::ParallelDescriptor::ReduceRealMax(lat_ell.data(),   ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_Tf.data(),    ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_depth.data(), ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_KI.data(),    ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_KIc.data(),   ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_q.data(),     ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_sig.data(),   ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_hf.data(),    ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_hfu.data(),   ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_sx_enc.data(), ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_sy_enc.data(), ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_wx.data(),     ncols);
                amrex::ParallelDescriptor::ReduceRealMax(lat_wy.data(),     ncols);
                lat_score = sc_max;
                amrex::ParallelDescriptor::ReduceLongSum(lat_cand);
                amrex::ParallelDescriptor::ReduceLongSum(lat_no_faceT);
                amrex::ParallelDescriptor::ReduceLongSum(lat_skip_hside);
                // ⚠ reduced here because these accumulate inside the per-BOX
                // criterion loop -- the rule P1g-1 and P1g-2 both paid for:
                // reduce what a per-box loop accumulated, never what a
                // per-column loop did.
                amrex::ParallelDescriptor::ReduceLongSum(lat_len_fb_h0);
                amrex::ParallelDescriptor::ReduceLongSum(lat_len_fb_q0);
                lat_len_fb_h0_n = static_cast<Set::Scalar>(lat_len_fb_h0);
                lat_len_fb_q0_n = static_cast<Set::Scalar>(lat_len_fb_q0);
                int nfire = 0; Set::Scalar esum = 0.0;
                for (int c = 0; c < ncols; ++c)
                    if (lat_kfire[c] != std::numeric_limits<int>::min())
                    { ++nfire; esum += lat_ell[c]; }
                // deliverable 12: one ReduceRealMax each, but the PAIR is kept
                // together -- reduce the score, then take the T only from the
                // ranks that attain it (ties resolve to the larger T, so the
                // answer is deterministic rather than decomposition-dependent).
                // This mirrors the lat_score / lat_kfire pattern above.
                {
                    Set::Scalar q_g = lat_qmax_l;
                    amrex::ParallelDescriptor::ReduceRealMax(q_g);
                    Set::Scalar t_g = (lat_qmax_l >= q_g) ? lat_Tfq_l : 0.0;
                    amrex::ParallelDescriptor::ReduceRealMax(t_g);
                    lat_qmax      = q_g;
                    lat_Tface_max = t_g;
                }
                lat_cand_n     = static_cast<Set::Scalar>(lat_cand);
                lat_fires      = static_cast<Set::Scalar>(nfire);
                lat_no_faceT_n = static_cast<Set::Scalar>(lat_no_faceT);
                lat_ell_med    = (nfire > 0) ? esum / Set::Scalar(nfire) : 0.0;
            }
        }

        // Step 19: sigma_xx-reconstruction selector for k_i_at_zcrack.
        // The T-derived form is only valid when surface_patch is ENABLED
        // AND in prescribed_T mode (only then is surface_patch.T_f
        // compiled to a real parser callable; an uncompiled T_f returns
        // std::numeric_limits<double>::max() which silently produces
        // garbage K_I values). Under prescribed_T + enabled we keep the
        // Step 16d byte-identical path. Under any other configuration
        // (surface_patch disabled, convective_flame, or no surface_patch
        // at all — e.g. the MMW beam case) sigma_xx is sampled directly
        // from the mechanics solve's stress_mf via NodeToCellAverage,
        // mirroring UpdateSpAfterMechanics's Sp_field diagnostic path.
        // Mechanics must be active for the stress_mf branch (have_stress
        // = el.time_evolving = 1); fail-fast otherwise.
        // Surface temperature for the T-face surface region: the prescribed-T
        // parser when surface_patch is in prescribed_T mode, else (e.g. the
        // local_thermoelastic flag under an MMW beam / computed-T source) the
        // top-cell temp_mf value. sp_weibull_use_T_face fires for EITHER
        // prescribed_T OR spall.stress_source = local_thermoelastic, so the
        // same T-derived sigma_xx path serves both without an MLMG solve.
        const bool sp_weibull_presc_T_face =
            surface_patch.enabled &&
            surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        const bool sp_weibull_use_T_face =
            sp_weibull_spall &&
            (sp_weibull_presc_T_face || sp_local_thermoelastic);
        if (sp_weibull_spall && !sp_weibull_use_T_face &&
            (!have_stress || mech_disabled))
            Util::Abort(INFO, "Step 19 (sp_weibull + spall.enabled under "
                              "non-prescribed_T heat source): mechanics "
                              "solve must be active to source sigma_xx "
                              "from stress_mf — set el.type = static and "
                              "el.time_evolving = 1, switch to "
                              "surface_patch.mode = prescribed_T, or set "
                              "spall.stress_source = local_thermoelastic");

        // Confining-pressure offset +p·nu/(1-nu) for the k_i_at_zcrack
        // lambda (Kant Eq. 15). Default-zero short-circuit preserves
        // byte-identity when confining.p is unset. Applies in BOTH the
        // Step 16d prescribed_T branch AND the Step 19 stress_mf branch
        // (compression-positive sigma_xx convention is the same).
        Set::Scalar sp_conf_offset = 0.0;
        if (sp_weibull_spall && confining_p != 0.0 &&
            phase_nu_loc > -1.0 && phase_nu_loc < 1.0)
        {
            sp_conf_offset = confining_p * phase_nu_loc / (1.0 - phase_nu_loc);
        }

        // Per-rank buffer of first-firing CSV events. Each column emits
        // once per run when cec(i, j, k_top) flips 0 -> 1; rows are
        // appended to <plot_file>_h_col_events.csv after PASS 2.
        struct HColEvent {
            Set::Scalar t, a_f, h_col, sp;
            int col_i, col_j, is_gb;
        };
        std::vector<HColEvent> local_h_col_events;

        // A1: every-removal event log (spall.removal_events_csv). Per-column
        // buffers, decomposition-independent: owner-only values reduced by
        // ReduceRealSum (non-owners contribute 0), n_voided by ReduceIntSum
        // after the void flip; the I/O rank writes one row per column with
        // regime_col != 0, ordered by (col_j, col_i). Level 0 only.
        // C1: the same per-column gather also feeds the pinned-surface closure
        // (surface_patch.robin_form = pinned), even with the CSV off; the CSV
        // itself stays gated on spall.removal_events_csv.
        const bool pin_on = sp_weibull_spall && SurfaceOwner(lev) && surface_patch.enabled &&
                            surface_patch.robin_form == SurfacePatch::RobinForm::Pinned;
        const bool rev_on = sp_weibull_spall && SurfaceOwner(lev) &&
                            (sp_emit_removal_events_csv || pin_on);
        // D2z-1: whole-cell removal is now the only mode, so these are no longer
        // keyed on anything -- `inc_on` is just "this level owns the surface"
        // and `inc_active` additionally excludes prescribed_T, where the scan is
        // retained (see the override below). The void counter is INDEPENDENT of
        // rev_nvoid, which only exists when the event CSV or the pin is on --
        // the assertion must hold on every run, not just logged ones.
        const bool inc_on = SurfaceOwner(lev);
        // ⚠ The assertion must be scoped EXACTLY like the override below
        // (sp_weibull_spall && !prescribed_T), or a damage-law run that
        // legitimately removes several cells in one event trips it. That is not
        // hypothetical: unit/spall_event removes two.
        const bool inc_active =
            inc_on && sp_weibull_spall && !sp_weibull_presc_T_face;
        amrex::Vector<int> inc_nvoid(inc_on ? ncols : 0, 0);
        amrex::Vector<Set::Scalar> rev_T_top, rev_a_f, rev_Sp, rev_h_scan;
        // ⚠⚠ P2-2b THE RULE: no site in the rock path reads `rho` or `Cp` as a
        // per-phase constant, and NO STORED PRE-MULTIPLIED PRODUCT. The retired
        // `rev_rhoCp`/`pin_rhoCp` held `rho*Cp`, which cannot carry temperature
        // dependence without being recomputed anyway.
        // ⚠ Gathered UNCONDITIONALLY (not under `rev_on`): site 3's consumer is
        // the lateral banking energy estimate, which runs in a per-column loop
        // with no MFIter in scope and under a DIFFERENT gate. Coupling it to
        // `rev_on` would tie two unrelated conditions together.
        amrex::Vector<Set::Scalar> col_rho, col_Cp;
        if (microstructure_enabled && ncols > 0)
        {
            col_rho.assign(ncols, 0.0);
            col_Cp.assign(ncols, 0.0);
        }
        amrex::Vector<int>         rev_nvoid;
        // ------------------------------------------------------------------
        // P0. The depth scan's surface endpoint (`T_face` below) is the SECOND
        // criterion site; UpdateSpAfterMechanics is the first. Both must read
        // the same surface temperature, or the Kant anchor and the Meier
        // firing run different formulas -- the disease this packet exists to
        // cure. prescribed_T is excluded: Step 16d already gives that path a
        // real face temperature.
        const bool recon_on = sp_surface_T_recon && sp_weibull_use_T_face &&
                              !sp_weibull_presc_T_face;
        // P1f deliverable 4: gates the three APPENDED event-log columns. ⚠ Keyed
        // off by default, so the frozen tree is untouched and no re-pin is owed.
        const bool lat_csv = sp_lateral_drive;
        // Per-column, owner-written, ReduceRealSum-collapsed (the rev_* /
        // PASS 1b pattern). rc_h_f is READ BY EVERY RANK in the winner loop,
        // so it must be reduced before that loop -- not merely logged.
        amrex::Vector<Set::Scalar> rc_Tsurf, rc_ell, rc_h_f, rc_Tbar, rc_flag;
        amrex::Vector<Set::Scalar> rc_ellq;   // P1g-3: the flux-based length used
        // P1 deliverable 5: rc_qrat = q_fit / q_real -- the flux the
        // two-cell fit IMPLIES, k*(T_surf - T0)/ell_eff, over the flux
        // the Robin BC DELIVERS. A receding front reads ~1; a thin hot
        // skin on rock pre-warmed from the side reads far below.
        // ⚠ LOGGED, NOT ACTED ON.
        // ⚠ P2-2a: rc_dfrac retired with `differential` -- it measured the
        // fraction of the absolute drive that form applied, and there is
        // now one form. Its own finding is kept in the packet: on a
        // self-similar front it was a CONSTANT 1 - e^{-c_f} = 0.330, i.e.
        // no flux information at all.
        amrex::Vector<Set::Scalar> rc_qrat;
        if (recon_on)
        {
            rc_Tsurf.assign(ncols, 0.0);
            rc_ell.assign(ncols, 0.0);
            rc_ellq.assign(ncols, 0.0);
            rc_h_f.assign(ncols, 0.0);
            rc_Tbar.assign(ncols, 0.0);
            rc_flag.assign(ncols, 0.0);   // 1 flat | 2 ell clamp | 4 T clamp
            rc_qrat.assign(ncols, 0.0);
            // Step-to-step variation of ell_eff: the two-cell ratio is exact
            // only in steady state and a firing event disturbs it, so the
            // packet asks for the disturbance to be MEASURED. ncols is the
            // domain column count and is regrid-invariant, so no repair.
            if (static_cast<int>(rc_ell_prev.size()) != ncols)
                rc_ell_prev.assign(ncols, 0.0);
        }
        if (rev_on)
        {
            rev_T_top.assign(ncols, 0.0);
            rev_a_f.assign(ncols, 0.0);
            rev_Sp.assign(ncols, 0.0);
            rev_h_scan.assign(ncols, 0.0);
            rev_nvoid.assign(ncols, 0);
        }

        // PASS 1b (depth-scan z-gather): collect each column's temperature
        // profile for a window of up to max_steps_cap cells below its surface,
        // collapsed across MPI z-decomposition via ReduceRealSum (each (i,j,k)
        // is owned by exactly one validbox; non-owning ranks contribute 0). The
        // K_I-vs-depth scan below then walks the FULL physical column instead of
        // clamping at the local z-box floor (k_top - klo) — that clamp otherwise
        // makes the detected flake thickness, and hence crater depth, depend on
        // amr.max_grid_size / rank count. Only the T-derived (local_thermoelastic
        // / prescribed_T) scan path consumes this; the stress_mf path keeps its
        // box-local read (sig_node is not gathered).
        const int dom_klo_z = domain.smallEnd(2);
        // D2q-1: the K_I depth scan strides downward from k_top in 1-cell
        // strides and reads each cell as material. This is how far that holds:
        // the contiguous solid run below each column's top. While void is
        // contiguous it reaches the domain floor, so the scan is unchanged and
        // the clamp below is inert -- but once stage 2 can leave an interior
        // void, the scan must not read through it as if it were rock.
        const amrex::Vector<int> solid_run_col =
            ColumnSolidRunBelowTop(lev, k_top_global);
        int Dmax_col = static_cast<int>(std::floor(h_spall_max_loc / dz));
        if (Dmax_col < 1) Dmax_col = 1;
        if (Dmax_col > domain.bigEnd(2) - dom_klo_z)
            Dmax_col = domain.bigEnd(2) - dom_klo_z;
        const int colw = Dmax_col + 1;
        amrex::Vector<Set::Scalar> T_col_prof;
        if (sp_weibull_use_T_face)
        {
            T_col_prof.assign(static_cast<std::size_t>(ncols) * colw, 0.0);
            for (amrex::MFIter mfi(*temp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).array(mfi);
                const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
                const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
                const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
                for (int j = jlo; j <= jhi; ++j)
                for (int i = ilo; i <= ihi; ++i)
                {
                    const int idx = col_index(i, j);
                    const int ktop = k_top_global[idx];
                    if (ktop == std::numeric_limits<int>::min()) continue;
                    for (int d = 0; d <= Dmax_col; ++d)
                    {
                        const int k = ktop - d;
                        if (k < klo || k > khi) continue;
                        T_col_prof[static_cast<std::size_t>(idx) * colw + d] = T(i, j, k);
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(
                T_col_prof.data(), static_cast<int>(T_col_prof.size()));
        }
        const Set::Scalar* T_col_p = T_col_prof.empty() ? nullptr : T_col_prof.data();

        // ⚠⚠ P2-2a change A: THE STRESS-FREE REFERENCE, gathered by the SAME
        // depth index as the temperature above.
        //
        // ⚠ Why a gather and not a box-local read: the depth scan walks the
        // FULL physical column (PASS 1b) and a column spans ranks under
        // z-decomposition, so the temperature is gathered. The reference has to
        // be addressed identically or the two would be read at different cells
        // at a box seam, and the bracket would carry a decomposition artefact
        // -- the D2q-1 defect shape.
        // ⚠ It is a READ of a never-updated field, so unlike the temperature
        // this profile is the same every step; it is re-gathered rather than
        // cached because caching needs invalidation on regrid and the
        // reduction is the same size as the one already being paid.
        amrex::Vector<Set::Scalar> T0_col_prof;
        if (sp_weibull_use_T_face &&
            T_init_mf.size() > static_cast<std::size_t>(lev) && T_init_mf[lev])
        {
            T0_col_prof.assign(static_cast<std::size_t>(ncols) * colw, 0.0);
            for (amrex::MFIter mfi(*T_init_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& T0 =
                    (*T_init_mf[lev]).const_array(mfi);
                const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                {
                    const int idx = col_index(i, j);
                    const int ktop = k_top_global[idx];
                    if (ktop == std::numeric_limits<int>::min()) continue;
                    for (int d = 0; d <= Dmax_col; ++d)
                    {
                        const int k = ktop - d;
                        if (k < klo || k > khi) continue;
                        T0_col_prof[static_cast<std::size_t>(idx) * colw + d] =
                            T0(i, j, k);
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(
                T0_col_prof.data(), static_cast<int>(T0_col_prof.size()));
        }
        const Set::Scalar* T0_col_p =
            T0_col_prof.empty() ? nullptr : T0_col_prof.data();

        // ------------------------------------------------------------------
        // P0b deliverable 4: THE CONSERVATIVE COLUMN REMAP.
        //
        // ⚠⚠ WHY NOT REVISION 1'S SUBTRACTION. `phi` is bound in the conduction
        // kernel but NEVER READ inside its ParallelFor, so a receded cell is
        // still a full dz of rock thermally. Subtracting the flake's energy
        // from the firing cell therefore leaves the remainder's energy spread
        // over the WRONG thickness, and ReconstructSurfaceT -- which assumes
        // two full cells dz apart -- misreads it: planner-verified at
        // 1.32x / 3.81x / FLAT for 2 / 1 / 0.5 mm. Sub-cell recession needs a
        // sub-cell thermal frame, and that frame is what this remap supplies.
        //
        // THE SCHEME. Removing a flake of thickness h_f shifts the profile up:
        // the new profile at depth z is the old one at z + h_f. Writing
        // h_f = m*dz + r and integrating the per-cell EXPONENTIAL sub-cell
        // profile that the reconstruction already supplies,
        //
        //     new[d] = (1 - b)*old[d + m] + b*old[d + m + 1],
        //     b      = (1 - e^{-r/ell}) / (1 - e^{-dz/ell}).
        //
        // Three things fall out, and they are the reason this is the right
        // scheme rather than merely a workable one:
        //  * ⚠ IT IS EXACTLY CONSERVATIVE. The sum telescopes, so the energy
        //    removed is b*(H_top - H_deep) = rho*c*h_f*(T_bar_f - T0) -- the
        //    flake's own enthalpy, EXACTLY, with no separate subtraction.
        //    Verified algebraically against T_bar = T0 + (T_surf-T0)(1-e^-cf)/cf.
        //  * The two-cell ratio is INVARIANT under the shift for an exponential
        //    (the profile keeps its shape, scaled by e^{-c_f}), so ell_eff is
        //    continuous across firings -- P3.
        //  * P0's "the surface is at the top face of the top cell" becomes true
        //    at all times, so the P0 algebra needs no change.
        //
        // ⚠ It is applied to H, not T: H is the conserved state, the blend is
        // linear, and the telescoping is then exact in the variable that must
        // be conserved.
        const bool remap_on = recon_on && sp_surface_T_apply &&
                              sp_removal_energy_flake;
        amrex::Vector<Set::Scalar> H_col_prof, phi_top_col;
        if (remap_on)
        {
            H_col_prof.assign(static_cast<std::size_t>(ncols) * colw, 0.0);
            phi_top_col.assign(ncols, 0.0);
            for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& Hh = (*H_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& pp_ = (*phi_mf[lev]).array(mfi);
                const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                {
                    const int idx = col_index(i, j);
                    const int ktop = k_top_global[idx];
                    if (ktop == std::numeric_limits<int>::min()) continue;
                    if (ktop >= klo && ktop <= khi) phi_top_col[idx] = pp_(i, j, ktop);
                    for (int d = 0; d <= Dmax_col; ++d)
                    {
                        const int k = ktop - d;
                        if (k < klo || k > khi) continue;
                        H_col_prof[static_cast<std::size_t>(idx) * colw + d] = Hh(i, j, k);
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(
                H_col_prof.data(), static_cast<int>(H_col_prof.size()));
            amrex::ParallelDescriptor::ReduceRealSum(phi_top_col.data(), ncols);
        }

        // P1 (melt-aware): gather each column's liquid-fraction profile the
        // same z-decomposition-safe way as the temperature profile above, so
        // the sp_weibull driving stress can be scaled by the load-bearing
        // skeleton f(Lambda_L) at the crack-tip depth on the box-independent
        // full column. Only when melt_aware && sp_weibull_spall; off => empty
        // buffer and every skeleton factor below is exactly 1.0 (byte-
        // identical). Lambda_L_mf has 0 ghosts; the gather reads validbox
        // cells only (each (i,j,k) owned by exactly one box), ReduceRealSum
        // collapses the z-decomposition.
        amrex::Vector<Set::Scalar> LL_col_prof;
        if (melt_aware && sp_weibull_spall)
        {
            LL_col_prof.assign(static_cast<std::size_t>(ncols) * colw, 0.0);
            for (amrex::MFIter mfi(*Lambda_L_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& LLa =
                    (*Lambda_L_mf[lev]).array(mfi);
                const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
                const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
                const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
                for (int j = jlo; j <= jhi; ++j)
                for (int i = ilo; i <= ihi; ++i)
                {
                    const int idx = col_index(i, j);
                    const int ktop = k_top_global[idx];
                    if (ktop == std::numeric_limits<int>::min()) continue;
                    for (int d = 0; d <= Dmax_col; ++d)
                    {
                        const int k = ktop - d;
                        if (k < klo || k > khi) continue;
                        LL_col_prof[static_cast<std::size_t>(idx) * colw + d] =
                            LLa(i, j, k);
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(
                LL_col_prof.data(), static_cast<int>(LL_col_prof.size()));
        }
        const Set::Scalar* LL_col_p =
            LL_col_prof.empty() ? nullptr : LL_col_prof.data();

        // P2-2b deliverable 4: the rock law for the DEPTH SCAN -- the firing
        // site the Meier run actually fires on. ⚠ Hoisted OUTSIDE the MFIter:
        // the selection is a run-level fact, and deciding it per box invites
        // the two-firing-site invariant to be broken box by box.
        const Numeric::Material::Grimsel *scan_g = RockLaw();

        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);
            // P2-2b: the per-cell rho and c_p. ⚠ READS, never written here.
            amrex::Array4<const Set::Scalar> rho_pc, Cp_pc;
            if (microstructure_enabled && rho_mf[lev] && Cp_mf[lev])
            {
                rho_pc = (*rho_mf[lev]).const_array(mfi);
                Cp_pc  = (*Cp_mf[lev]).const_array(mfi);
            }
            amrex::Array4<const Set::Matrix> sig_node;
            if (have_stress && !mech_disabled) sig_node = (*stress_mf[lev]).array(mfi);
            // Step 16: per-cell h_spall diagnostic write target + the
            // Sp_cluster_id detach gate / flaw_a / temperature reads for
            // the depth scan. All four fields exist iff sp_weibull_spall
            // is true (Register block in Parse()).
            // Step 16b: crack_event_count_mf write target (edge-triggered
            // 3-D first-crossing count field, design note §2).
            amrex::Array4<Set::Scalar>      hsf;
            amrex::Array4<Set::Scalar>      cec;
            amrex::Array4<const Set::Scalar> clu, fa_a, T_a;
            // Step 16c: bind Sp_field + is_grain_boundary on the top
            // slice to populate the per-column h_col CSV event row.
            amrex::Array4<const Set::Scalar> Sp_top, is_gb_top;
            if (sp_weibull_spall)
            {
                hsf       = (*h_spall_field_mf[lev]).array(mfi);
                cec       = (*crack_event_count_mf[lev]).array(mfi);
                // A1: Sp_cluster_id / flaw_a exist only with weibull.enabled
                // (Parse register block); the scalar-flaw path reads neither.
                if (weibull_enabled)
                {
                    clu   = (*Sp_cluster_id_mf[lev]).array(mfi);
                    fa_a  = (*flaw_a_mf[lev]).array(mfi);
                }
                T_a       = (*temp_mf[lev]).array(mfi);
                Sp_top    = (*Sp_field_mf[lev]).array(mfi);
                is_gb_top = (*is_grain_boundary_mf[lev]).array(mfi);
            }

            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);

            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                int idx = col_index(i, j);
                int k_top = k_top_global[idx];
                if (k_top == std::numeric_limits<int>::min()) continue;
                if (k_top < klo || k_top > khi) continue; // not owned here

                // A1: owner-only per-column values for the removal-event log
                // (non-owners leave 0; ReduceRealSum after PASS 2).
                if (rev_on)
                {
                    rev_T_top[idx] = T_a(i, j, k_top);
                    rev_a_f[idx]   = weibull_enabled ? fa_a(i, j, k_top) : sp_a0;
                    rev_Sp[idx]    = Sp_top(i, j, k_top);
                }
                // P2-2b site 1: the pin reference, from the PER-CELL fields at
                // the column's own top cell, carried as two numbers.
                // ⚠ Outside the `rev_on` block on purpose -- site 3 needs these
                // under a different gate. Under `law = constant` the fields hold
                // the phase scalars, so this is arithmetically what it replaced.
                if (!col_rho.empty() && rho_pc && Cp_pc)
                {
                    col_rho[idx] = rho_pc(i, j, k_top);
                    col_Cp[idx]  = Cp_pc (i, j, k_top);
                }

                int p_top = static_cast<int>(std::round(ph(i, j, k_top)));
                if (p_top < 0 || p_top >= n_phases) p_top = 0;

                Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar zc = PLO[2] + (Set::Scalar(k_top) + 0.5) * DX[2];

                if (spall_on)
                {
                    if (sp_weibull_spall)
                    {
                        // Step 16: Sp-vs-depth zero-crossing h_spall.
                        // Detach gate is Sp_cluster_id_mf > 0 (consumes
                        // the §15c labeller — propagates per_face vs
                        // connected_cluster + A_crit cleanly to material
                        // removal). h_spall is the deepest crack-tip depth
                        // z_crack at which K_I(z_crack; a_f) >= K_Ic(T(z))
                        // for the hypothetical edge crack of fixed length
                        // a_f starting at z_crack below the surface.
                        // A1: without weibull.enabled there is no cluster
                        // labeller; its per_face gate is Sp_field >= 1 at the
                        // column's surface cell, and the flaw is the scalar
                        // spallation.sp.a0 (as in UpdateSpAfterMechanics).
                        const bool detach_true = weibull_enabled
                            ? (clu(i, j, k_top) > 0.0)
                            : (Sp_top(i, j, k_top) >= 1.0);
                        if (detach_true)
                        {
                            Set::Scalar a_f = weibull_enabled ? fa_a(i, j, k_top) : sp_a0;
                            // ⚠⚠ P2-2d rev 6 DEFECT 2 -- THE THIRD FLAW READ,
                            // AND THE ONE THAT HAD NO SIZE TERM.
                            //
                            // A removal needs BOTH `detach_true` (the cluster
                            // gate, whose flaw IS size-scaled at
                            // MMWSpalling.H:7403) AND `h_col > 0` from this
                            // thickness scan (:2001). This read fed the scan
                            // and took the flaw RAW, so the effective flaw was
                            // the SMALLER of the two: the size term could
                            // DELAY an onset at high flux but never ADVANCE
                            // one at low flux. Measured first-row firing
                            // ratios were 1.04-1.08 at q = 300 against
                            // ~1.00 at q = 500/800, so the low-flux end was
                            // pinned to the unscaled level and every
                            // three-point slope read about half the lever.
                            //
                            // THE THREE READS, and their equality BY
                            // CONFIGURATION:
                            //   MMWSpalling.H:7403  criterion / Sp_field  scaled
                            //   Removal.H:823       lateral face          scaled
                            //   HERE                the thickness scan    scaled now
                            // They are EQUAL under `size_volume = cell` (the
                            // factor is identically 1) and under a
                            // deterministic flaw with `cell`; they are NOT
                            // equal under `flux`, where each read scales by
                            // its own column's length. That by-configuration
                            // statement is what the lock asserts -- its old
                            // "exactly 2 sites" clause had MECHANISED the
                            // omission.
                            if (sp_size_flux)
                            {
                                const Set::Scalar lq_sc =
                                    (idx < static_cast<int>(surf_ellq_col.size()))
                                    ? surf_ellq_col[idx] : Set::Scalar(0.0);
                                if (lq_sc > 0.0)
                                    a_f *= WeibullSizeFactor(lq_sc * lq_sc * lq_sc);
                            }
                            // Step 16c: capture cec(i, j, k_top) BEFORE the
                            // scan modifies it, so the post-scan check can
                            // detect a 0 -> 1 transition (i.e. this column's
                            // first firing this run).
                            const Set::Scalar cec_at_top_pre =
                                cec(i, j, k_top);
                            if (std::isfinite(a_f) && a_f > 0.0 && have_stress)
                            {
                                // Walk z_crack downward from the surface in
                                // 1-cell strides; stop at the first
                                // K_I < K_Ic crossing, h_spall_max, or the
                                // bottom of the local box. Also drives the
                                // edge-triggered crack_event_count_mf (each
                                // (i, j, k_tip) flips 0 -> 1 on first pass,
                                // never re-tested — freeze, option (ii))
                                // and the sub-cell bisection that refines
                                // h_col below dz.
                                Set::Scalar h_col = 0.0;
                                // Step 20: surface-normal frame factors.
                                // Depth variables in this scan (z_crack, d,
                                // h_col) are NORMAL distances when the flag
                                // is on; the cell at k_top - s sits at
                                // normal depth s*dz*cos(theta). Exactly 1.0
                                // when off (bit-identical arithmetic).
                                const Set::Scalar sn_cth =
                                    sn_cos_p ? sn_cos_p[idx] : 1.0;
                                const Set::Scalar sn_icth = 1.0 / sn_cth;
                                const int max_steps_cap = static_cast<int>(
                                    std::floor(h_spall_max_loc / dz));
                                const int klo_loc = klo;
                                const int khi_loc = khi;
                                // T-derived path walks the gathered full column
                                // down to the domain floor (box-independent);
                                // stress_mf path stays box-local. See PASS 1b.
                                const int scan_kbot =
                                    sp_weibull_use_T_face ? dom_klo_z : klo;
                                // D2q-1: also clamp to the contiguous solid run
                                // below the top (inert while void is contiguous).
                                const int max_steps =
                                    std::min(std::min(max_steps_cap, k_top - scan_kbot),
                                             solid_run_col[idx] - 1);
                                // Column temperature by depth-below-surface index
                                // (cells): gathered profile for the T-path (box-
                                // independent), box-local T_a otherwise.
                                const Set::Scalar* Tcol_i =
                                    (T_col_p && sp_weibull_use_T_face)
                                      ? (T_col_p + static_cast<std::size_t>(idx) * colw)
                                      : nullptr;
                                // P2-2a change A: this column's stress-free
                                // reference profile, same depth indexing.
                                const Set::Scalar* T0col_i =
                                    (T0_col_p && sp_weibull_use_T_face)
                                      ? (T0_col_p + static_cast<std::size_t>(idx) * colw)
                                      : nullptr;
                                auto T_at_didx = [&] (int d) -> Set::Scalar {
                                    if (Tcol_i)
                                    {
                                        if (d < 0)         d = 0;
                                        if (d > Dmax_col)  d = Dmax_col;
                                        return Tcol_i[d];
                                    }
                                    int kk = k_top - d;
                                    if (kk < klo_loc) kk = klo_loc;
                                    if (kk > khi_loc) kk = khi_loc;
                                    return T_a(i, j, kk);
                                };
                                // ⚠⚠ P2-2a change A: THE STRESS-FREE REFERENCE
                                // by depth index -- the SAME clamping and the
                                // SAME gathered addressing as T_at_didx above,
                                // so the bracket is a difference of two reads
                                // at ONE cell.
                                // ⚠ A CELL READ, not a lerp. The temperature is
                                // interpolated inside the top half-cell (the P0
                                // surface lerp); the reference is not, because
                                // it is the state the cell STARTED in and has no
                                // sub-cell profile. On every shipped deck the IC
                                // is `constant`, so the two agree identically;
                                // the distinction only bites a graded IC, and
                                // then the cell value is the defensible one.
                                // ⚠ Falls back to the DECLARED constant when the
                                // field is absent or non-finite -- the pre-P2-2a
                                // behaviour exactly, so a cell with no stored
                                // value cannot silently acquire a zero bracket.
                                auto T0_at_didx = [&] (int d) -> Set::Scalar {
                                    Set::Scalar v;
                                    if (T0col_i)
                                    {
                                        if (d < 0)         d = 0;
                                        if (d > Dmax_col)  d = Dmax_col;
                                        v = T0col_i[d];
                                    }
                                    else return phase_T_ref_loc;
                                    return std::isfinite(v) ? v : phase_T_ref_loc;
                                };

                                // P0: THE SURFACE TEMPERATURE. Recovered from
                                // the top two cell AVERAGES, which for an
                                // exponential front have the exact ratio
                                // e^{dz/l_eff}. Needs a solid second cell --
                                // the same contiguous-solid-run test that
                                // already bounds the scan.
                                SurfaceRecon rec;
                                bool rec_ok = false;
                                // P0b: read the ONE reconstruction built at
                                // the top of the step (BuildSurfaceTColumns),
                                // not a second gather here. P0 had two copies
                                // of the same estimate; they agreed, but
                                // neither could carry the smoothing or the
                                // flat-hold, which are per-column state ACROSS
                                // steps.
                                if (recon_on &&
                                    SharedSurfaceRecon(
                                        static_cast<std::size_t>(idx), rec))
                                {
                                    // `diagnostic` logs it without
                                    // feeding it to the criterion.
                                    rec_ok = sp_surface_T_apply;
                                    rc_Tsurf[idx] = rec.T_surf;
                                    rc_ellq[idx]  = std::isfinite(rec.ell_q)
                                                  ? rec.ell_q : 0.0;
                                    rc_ell[idx]   = std::isfinite(rec.ell_eff)
                                                    ? rec.ell_eff : 0.0;
                                    rc_h_f[idx]   = rec.h_f;
                                    rc_Tbar[idx]  = rec.T_bar;
                                    // ⚠ rev7: bits 8 and 16 were MISSING here.
                                    // The event log's rc_flag encoded only
                                    // flat/ell/T, so "the h_f cap binds on 0 %
                                    // of firings" -- read off this column in
                                    // revision 2 -- was reading a bit that was
                                    // never written. Both thickness bounds are
                                    // now reported where the firings are.
                                    rc_flag[idx]  = (rec.flat        ? 1.0 : 0.0)
                                                  + (rec.ell_clamped ? 2.0 : 0.0)
                                                  + (rec.T_clamped   ? 4.0 : 0.0)
                                                  + (rec.h_f_capped  ? 8.0 : 0.0)
                                                  + (rec.h_f_bounded ? 16.0 : 0.0);
                                    {
                                        const Set::Scalar dT =
                                            rec.T_surf - phase_T_ref_loc;
                                        if (dT > 0.0 && rec.ell_eff > 0.0 &&
                                            std::isfinite(rec.ell_eff) &&
                                            idx < static_cast<int>(flame_h.size()))
                                        {
                                            const Set::Scalar q_fit =
                                                phase_kappa_loc * dT / rec.ell_eff;
                                            const Set::Scalar Ts = rec.T_surf;
                                            const Set::Scalar eps =
                                                (Ts <= T_eps_break) ? eps_rad
                                                                    : eps_rad_hi;
                                            const Set::Scalar q_real =
                                                flame_h[idx]
                                                * (flame_Tg[idx] - Ts)
                                                - eps * sigma_SB
                                                  * (Ts*Ts*Ts*Ts
                                                     - T_amb*T_amb*T_amb*T_amb)
                                                - h_conv * (Ts - T_amb);
                                            if (std::abs(q_real) > 0.0)
                                                rc_qrat[idx] = q_fit / q_real;
                                        }
                                    }
                                }

                                // P1: liquid fraction Lambda_L by depth-below-
                                // surface index from the gathered column (box-
                                // independent). nullptr unless melt_aware, in
                                // which case LL_at_didx is never called and the
                                // skeleton factor stays 1.0.
                                const Set::Scalar* LLcol_i =
                                    (LL_col_p && melt_aware)
                                      ? (LL_col_p +
                                         static_cast<std::size_t>(idx) * colw)
                                      : nullptr;
                                auto LL_at_didx = [&] (int d) -> Set::Scalar {
                                    if (LLcol_i)
                                    {
                                        if (d < 0)        d = 0;
                                        if (d > Dmax_col) d = Dmax_col;
                                        return LLcol_i[d];
                                    }
                                    return 0.0;
                                };

                                // P0 deliverable 3: T(z_v) extracted from sigma_at_depth so
                                // that the driving stress AND K_Ic read the SAME profile --
                                // "by construction", not as two copies that can drift. The
                                // arithmetic is the pre-P0 arithmetic verbatim; only the
                                // z = 0 endpoint changes, and only under the key.
                                auto T_at_zv = [&] (Set::Scalar z_v) -> Set::Scalar {
                                    const Set::Scalar half_dz =
                                        0.5 * dz;
                                    Set::Scalar T_z;
                                    if (z_v <= half_dz)
                                    {
                                        // Surface region: lerp between the
                                        // surface T and the top-cell-centre
                                        // T. Surface T is the prescribed-T
                                        // parser under prescribed_T; else
                                        // the P0 reconstruction when
                                        // spall.surface_T = reconstructed;
                                        // else (pre-P0) the top cell
                                        // itself -- which is what made the
                                        // whole top half-cell isothermal
                                        // at the CELL AVERAGE.
                                        const Set::Scalar T_top =
                                            T_at_didx(0);
                                        const Set::Scalar T_face =
                                            sp_weibull_presc_T_face
                                            ? sp_T_face_loc(time)
                                            : (rec_ok ? rec.T_surf
                                                      : T_top);
                                        Set::Scalar alpha =
                                            z_v / half_dz;
                                        if (alpha < 0.0) alpha = 0.0;
                                        if (alpha > 1.0) alpha = 1.0;
                                        T_z = (1.0 - alpha) * T_face
                                            +        alpha  * T_top;
                                    }
                                    else
                                    {
                                        // Interior: linear T interp
                                        // between adjacent cell
                                        // centres at z = (k + 0.5)*dz.
                                        const Set::Scalar u =
                                            z_v / dz - 0.5;
                                        const Set::Scalar uf =
                                            std::floor(u);
                                        int n = static_cast<int>(uf);
                                        if (n < 0) n = 0;
                                        const Set::Scalar alpha =
                                            u - uf;  // in [0, 1)
                                        // Depth-below-surface indices n,
                                        // n+1; T_at_didx clamps and reads the
                                        // gathered column (box-independent).
                                        T_z = (1.0 - alpha)
                                              * T_at_didx(n)
                                            +        alpha
                                              * T_at_didx(n + 1);
                                    }
                                    return T_z;
                                };

                                // Local helper: K_Ic(T) under the
                                // configured law, in Pa·m^0.5.
                                auto kic_pa_at_cell =
                                    [&] (int k_cell) -> Set::Scalar {
                                    Set::Scalar T_at = T_at_didx(k_top - k_cell);
                                    Set::Scalar K_Ic_MPa;
                                    // P2-2d d7, toughness sites 3 and 4 of 4
                                    // (the DEPTH SCAN -- the site Kant fires
                                    // on). ⚠ These two blocks are textually
                                    // identical and BOTH must carry the route;
                                    // converting one would make the scan
                                    // disagree with itself by branch.
                                    if (law_loc == KIcLaw::GrimselTable && kic_g_loc)
                                        K_Ic_MPa = kic_scale_loc
                                                 * kic_g_loc->K_Ic_of(T_at);
                                    else if (law_loc == KIcLaw::NasseriTable)
                                        K_Ic_MPa =
                                            Numeric::SpCriterion::KIcNasseri(
                                                T_at, kic_scale_loc);
                                    else
                                        K_Ic_MPa =
                                            Numeric::SpCriterion::KIcLinear(
                                                T_at, kic_K_Ic0_loc,
                                                kic_T_ref_loc, kic_T_melt_loc);
                                    return Numeric::SpCriterion
                                        ::MPa_sqrtm_to_Pa_sqrtm(K_Ic_MPa);
                                };

                                // P0 deliverable 3: K_Ic AT THE CRACK TIP.
                                // The pre-P0 form evaluates K_Ic at the tip
                                // CELL AVERAGE while sigma_xx is integrated
                                // over a continuous T(z), so at a marginal
                                // firing the two sides of Sp = K_I/K_Ic are
                                // read off different profiles. Under the key
                                // both come from T_at_zv.
                                auto kic_pa_at_z =
                                    [&] (Set::Scalar z_v) -> Set::Scalar {
                                    const Set::Scalar T_at = T_at_zv(z_v);
                                    Set::Scalar K_Ic_MPa;
                                    // P2-2d d7, toughness sites 3 and 4 of 4
                                    // (the DEPTH SCAN -- the site Kant fires
                                    // on). ⚠ These two blocks are textually
                                    // identical and BOTH must carry the route;
                                    // converting one would make the scan
                                    // disagree with itself by branch.
                                    if (law_loc == KIcLaw::GrimselTable && kic_g_loc)
                                        K_Ic_MPa = kic_scale_loc
                                                 * kic_g_loc->K_Ic_of(T_at);
                                    else if (law_loc == KIcLaw::NasseriTable)
                                        K_Ic_MPa =
                                            Numeric::SpCriterion::KIcNasseri(
                                                T_at, kic_scale_loc);
                                    else
                                        K_Ic_MPa =
                                            Numeric::SpCriterion::KIcLinear(
                                                T_at, kic_K_Ic0_loc,
                                                kic_T_ref_loc, kic_T_melt_loc);
                                    return Numeric::SpCriterion
                                        ::MPa_sqrtm_to_Pa_sqrtm(K_Ic_MPa);
                                };

                                // K_I(z_crack; a_f). Two sigma_xx
                                // reconstructions, selected by
                                // sp_weibull_use_T_face:
                                //   (Step 16d, prescribed_T)  T-derived:
                                //     T(z) lerp between surface_patch.T_f
                                //     at z=0 and cell-centred T at
                                //     z=(n+0.5)*dz; sigma_xx = -E*beta*
                                //     dT/(1-nu) under 1-D confinement
                                //     (single-phase, phases[0] only).
                                //   (Step 19, non-prescribed_T)  stress_mf:
                                //     sigma_xx sampled from the mechanics
                                //     solve at the cell containing
                                //     z_crack + d, same pattern as
                                //     UpdateSpAfterMechanics's Sp_field
                                //     diagnostic. Captures the full FEM
                                //     stress state (lateral free-BC
                                //     effects, heterogeneous Voronoi,
                                //     etc.) — no 1-D-confinement
                                //     assumption.
                                // sp_conf_offset (Step 18) applies to
                                // both branches (compression-positive
                                // convention is the same).
                                // P1b: the support rule's FORCED DRIVE
                                // FRACTION, <= 0 = off. The ordinary criterion
                                // never sets it, so every branch below is
                                // byte-identical on every existing leg.
                                auto k_i_at_zcrack =
                                    [&] (Set::Scalar z_crack) -> Set::Scalar {
                                    auto sigma_at_depth =
                                        [&] (Set::Scalar d) -> Set::Scalar {
                                        // z_abs is a NORMAL distance under
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
                                        if (!sp_weibull_use_T_face)
                                        {
                                            // Step 19: stress_mf path.
                                            auto sx = [&] (Set::Scalar zz)
                                                      -> Set::Scalar {
                                                int k_d = k_top -
                                                    static_cast<int>(
                                                        std::floor(zz / dz));
                                                if (k_d > k_top)   k_d = k_top;
                                                if (k_d < klo_loc) k_d = klo_loc;
                                                if (k_d > khi_loc) k_d = khi_loc;
                                                Set::Matrix sg =
                                                    Numeric::Interpolate
                                                    ::NodeToCellAverage(
                                                        sig_node, i, j, k_d, 0);
                                                return -sg(0, 0);
                                            };
                                            // P1: the stress JUMP across the
                                            // flake. sp_conf_offset cancels --
                                            // correct, a uniformly confined
                                            // layer has no delamination drive.
                                            return f_skel * (sx(z_v)
                                                   + sp_conf_offset);
                                        }
                                        const Set::Scalar T_z = T_at_zv(z_v);
                                        // P2-2a change A: the reference at the
                                        // cell holding this sample depth.
                                        // ⚠⚠ P2-2e: under constraint_depth > 0
                                        // the LIVE temperature d further
                                        // behind the face, from the column
                                        // gathered before EITHER site ran --
                                        // the same values the Sp_field site
                                        // reads. Normal depth z_abs + d, so the
                                        // vertical read is (z_abs + d)*sn_icth
                                        // (identical to z_v + d with Step 20
                                        // off). METRES, interpolated: never a
                                        // cell index. d = 0 is the pre-P2-2e
                                        // call, bit for bit.
                                        const Set::Scalar T0_z = cref_on
                                            ? CrefTAtDepth(
                                                  static_cast<std::size_t>(idx),
                                                  (z_abs + cref_d_loc) * sn_icth
                                                  + Set::Scalar(cref_ktop[idx] - k_top) * dz,
                                                  dz, cref_vert, cref_clamp)
                                            : T0_at_didx(
                                                  static_cast<int>(
                                                      std::floor(z_v / dz)));
                                        if (cref_on && T_z < T0_z) ++cref_neg;
                                        // ⚠⚠ P1: the profile comes from the
                                        // EXPONENTIAL reconstruction, never
                                        // from T_at_zv -- that lambda
                                        // interpolates between CELL CENTRES,
                                        // so the same physical depth would read
                                        // differently at 2 / 1 / 0.5 mm and the
                                        // mesh criterion would catch it first.
                                        // 1-D confinement: physical sigma_xx
                                        // = -E*beta*dT/(1-nu). K_I integrator
                                        // wants compression-positive, so the
                                        // double-negation returns +E*beta*dT/
                                        // (1-nu). + sp_conf_offset adds the
                                        // confining-pressure contribution.
                                        // P2-2b deliverable 4: the law, at the
                                        // LOCAL temperature of the criterion's
                                        // own reconstructed profile.
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
                                    };
                                    return Numeric::SpCriterion::K_I(
                                        sigma_at_depth, a_f,
                                        h_spall_n_seg_loc);
                                };

                                int  s_last_pass = -1;
                                bool found_failure = false;
                                for (int s = 0; s <= max_steps; ++s)
                                {
                                    // Cell k_top - s sits at normal depth
                                    // s*dz*cos(theta) below the inclined
                                    // surface (s*dz vertical).
                                    // E1: z_crack = 0 is the TOP FACE of
                                    // k_top (cell-quantised; T lerps to the
                                    // k_top centre at dz/2), while h_col is
                                    // applied to the level set phi, which
                                    // sits within dz of that face. The beam
                                    // path now also starts at that face. The
                                    // <= dz offset between phi and the crack
                                    // origin is DELIBERATE: do not harmonise.
                                    const Set::Scalar z_crack =
                                        Set::Scalar(s) * dz * sn_cth;
                                    const int k_tip = k_top - s;
                                    if (k_tip < scan_kbot) break;

                                    Set::Scalar K_Ic_Pa = rec_ok
                                        ? kic_pa_at_z((z_crack + a_f) * sn_icth)
                                        : kic_pa_at_cell(k_tip);
                                    Set::Scalar K_I_val =
                                        k_i_at_zcrack(z_crack);

                                    if (K_I_val >= K_Ic_Pa)
                                    {
                                        s_last_pass = s;
                                        h_col       = z_crack;
                                        // Step 16b: edge-triggered event
                                        // count at (i, j, k_tip). Only writable
                                        // for cells in the surface-owning box;
                                        // under the depth-scan z-gather k_tip may
                                        // descend into a box this rank does not
                                        // own, so guard against an OOB write
                                        // (deep crossings simply go uncounted in
                                        // this secondary diagnostic).
                                        if (k_tip >= klo && k_tip <= khi &&
                                            cec(i, j, k_tip) < 0.5)
                                            cec(i, j, k_tip) = 1.0;
                                    }
                                    else
                                    {
                                        found_failure = true;
                                        break;
                                    }
                                }

                                // Step 16b: sub-cell binary search if the
                                // integer scan brackets the crossing
                                // between s_last_pass and s_last_pass + 1
                                // AND we have an in-box deeper cell to
                                // evaluate. Refines h_col to dz / 256
                                // (~6 um at dz = 1.5625 mm) — the GB/IG
                                // quantization fix from design note §4.
                                if (found_failure && s_last_pass >= 0 &&
                                    s_last_pass < max_steps &&
                                    (k_top - (s_last_pass + 1)) >= scan_kbot)
                                {
                                    Set::Scalar z_lo = Set::Scalar(s_last_pass)
                                                       * dz * sn_cth; // pass
                                    Set::Scalar z_hi = Set::Scalar(s_last_pass + 1)
                                                       * dz * sn_cth; // fail
                                    for (int b = 0; b < 8; ++b)
                                    {
                                        Set::Scalar z_mid =
                                            0.5 * (z_lo + z_hi);
                                        int k_mid = k_top -
                                            static_cast<int>(
                                                std::floor(z_mid * sn_icth / dz));
                                        // kic_pa_at_cell reads via T_at_didx,
                                        // which is box-independent on the T-path;
                                        // clamp to [scan_kbot, k_top] (== box on
                                        // the stress_mf path).
                                        if (k_mid < scan_kbot) k_mid = scan_kbot;
                                        if (k_mid > k_top)      k_mid = k_top;
                                        Set::Scalar K_Ic_mid = rec_ok
                                            ? kic_pa_at_z((z_mid + a_f) * sn_icth)
                                            : kic_pa_at_cell(k_mid);
                                        Set::Scalar K_I_mid =
                                            k_i_at_zcrack(z_mid);
                                        if (K_I_mid >= K_Ic_mid)
                                            z_lo = z_mid;
                                        else
                                            z_hi = z_mid;
                                    }
                                    h_col = z_lo;
                                }

                                if (h_col > 0.0 && detach_true)
                                {
                                    // Step 20: h_col is the flake thickness
                                    // along the surface normal; the column
                                    // (height-function) representation must
                                    // recede VERTICALLY by h_col/cos(theta)
                                    // to remove that flake. The graph form
                                    // of dphi/dt = -v_n*|grad phi|. The CSV
                                    // and h_spall_field diagnostics keep the
                                    // physical (normal) thickness.
                                    const Set::Scalar h_vert =
                                        h_col * sn_icth;
                                    spall_fires_col[idx] = 1;
                                    h_s_cand_col[idx]    = h_vert;
                                    if (rev_on) rev_h_scan[idx] = h_col;
                                    rop_s_cand_col[idx]  =
                                        (dt_local > 0.0) ? h_vert / dt_local : 0.0;
                                    // Running max at the INITIAL top-z
                                    // slice (immovable; column anchor for
                                    // the Rossi depth-profile post-process).
                                    // Each column records the deepest h_col
                                    // ever observed across the simulation.
                                    const int top_k_dom = domain.bigEnd(2);
                                    if (top_k_dom >= klo && top_k_dom <= khi)
                                    {
                                        Set::Scalar prev =
                                            hsf(i, j, top_k_dom);
                                        if (h_col > prev) hsf(i, j, top_k_dom) = h_col;
                                    }
                                    // Step 16c: first-firing CSV event.
                                    // Detect the cec(i, j, k_top) 0 -> 1
                                    // transition from this step's integer
                                    // scan. The freeze guard inside the
                                    // scan already set cec to 1 if any
                                    // crossing passed; we use the pre-loop
                                    // snapshot to decide whether THIS step
                                    // is the column's first-firing step.
                                    if (cec_at_top_pre < 0.5)
                                    {
                                        HColEvent ev;
                                        ev.t      = time;
                                        ev.col_i  = i;
                                        ev.col_j  = j;
                                        ev.a_f    = a_f;
                                        ev.h_col  = h_col;
                                        ev.is_gb  = (is_gb_top(i, j, k_top)
                                                     >= 0.5) ? 1 : 0;
                                        ev.sp     = Sp_top(i, j, k_top);
                                        local_h_col_events.push_back(ev);
                                    }
                                }

                                // ----------------------------------------
                                // P1b: UNSUPPORTED HOT ROCK COMES OFF.
                                //
                                // Only for a column that did NOT already fire
                                // this step -- the ordinary criterion keeps
                                // priority, so the rule can never double-remove
                                // and n_voided stays 1 (D2z-1's assertion,
                                // extended to this regime below).
                                //
                                // One K_I evaluation, at the surface, through
                                // the SAME k_i_at_zcrack the criterion just
                                // used: same flaw a_f, same toughness law, same
                                // integral, only the drive fraction forced to
                                // its self-similar value. There is no second
                                // temperature constant anywhere in this block
                                // -- the threshold it implies is the
                                // criterion's own, ~768 K at the shipped Kant
                                // calibration, and the parser refuses the
                                // absolute form precisely so it stays there.
                            }
                        }
                    }
                    else
                    {
                        // Step 14b: D sampled per the configured spall.sample
                        // mode (PASS 1.5 above). top_cell mode preserves the
                        // pre-Step-14b behaviour bit-for-bit because the scan
                        // window collapses to k = k_top.
                        Set::Scalar Dv = D_max_col[idx];
                        if (!std::isfinite(Dv)) Dv = 0.0;
                        if (Dv >= D_th)
                        {
                            Set::Scalar ratio = 0.0;
                            if (spall_prescribed)
                            {
                                ratio = ratio_f(xc, yc, zc, time);
                                if (!std::isfinite(ratio)) ratio = 0.0;
                            }
                            else
                            {
                                Set::Matrix sigma = Numeric::Interpolate::NodeToCellAverage(
                                    sig_node, i, j, k_top, 0);
                                ratio = PrincipalStressRatio(sigma);
                            }
                            if (ratio >= xi_s)
                            {
                                spall_fires_col[idx] = 1;
                                h_s_cand_col[idx]    = h_spall;
                                rop_s_cand_col[idx]  = (dt_local > 0.0) ? h_spall / dt_local : 0.0;
                            }
                        }
                    }
                }

                if (vapor_on && ptabs)
                {
                    Set::Scalar H_top = H(i, j, k_top);
                    Set::Scalar H_vap = ptabs[p_top].H_of_T(phases[p_top].T_vap_lo);
                    // P2-2b site 6: the per-cell values at the top cell. ⚠ The
                    // worst place for a constant c_p -- L_tot integrates one
                    // across the whole span up to vaporisation, over which
                    // Grimsel's c_p rises ~1.5x. ⚠ Inert on both production
                    // decks (`vapor.enabled = 0`); its witness is unit/latent_heat.
                    Set::Scalar Cp_p  = (Cp_pc  ? Cp_pc (i, j, k_top) : 0.0);
                    Set::Scalar rho_p = (rho_pc ? rho_pc(i, j, k_top) : 0.0);
                    Set::Scalar L_tot = Cp_p * (phases[p_top].T_vap_lo - T0_v)
                                      + phases[p_top].L_m + phases[p_top].L_v;
                    Set::Scalar rop_v = 0.0;
                    if (vapor_prescribed)
                    {
                        // Prescribed-RoP path stays as-is for synthetic
                        // unit tests (e.g. regime_low_high_power).
                        rop_v = vapor_rop_eval(xc, yc, zc, time);
                        if (!std::isfinite(rop_v) || rop_v < 0.0) rop_v = 0.0;
                    }
                    else if (sp_weibull_spall)
                    {
                        // Step 17: LOCAL v_vap from per-cell Gaussian
                        // beam power flux at the surface (W/m²). At the
                        // top face z = z0, the beam waist is w0 with no
                        // spreading:
                        //   Q_local(x, y, t) = (2·P(t) / (π·w0²))
                        //                       · exp(-2·r² / w0²)
                        // r² = (x - x0)² + (y - y0²). Compared to the
                        // bulk averaged form below (used under
                        // damage_law), the local v_vap at the beam
                        // centre equals the bulk value, and decays
                        // exponentially off-centre — matching the
                        // revised-approach §5 "interpreted locally
                        // rather than averaged" prescription.
                        if (std::isfinite(L_tot) && L_tot > 0.0 &&
                            rho_p > 0.0 && bw0_v > 0.0)
                        {
                            const Set::Scalar pi  = 3.141592653589793238;
                            const Set::Scalar r2  = (xc - bx0_v) * (xc - bx0_v)
                                                  + (yc - by0_v) * (yc - by0_v);
                            const Set::Scalar Q_local =
                                (2.0 * P_now / (pi * bw0_v * bw0_v))
                                * std::exp(-2.0 * r2 / (bw0_v * bw0_v));
                            rop_v = oneR * Q_local / (rho_p * L_tot);
                        }
                        if (!std::isfinite(rop_v) || rop_v < 0.0) rop_v = 0.0;
                    }
                    else
                    {
                        // damage_law: BULK-averaged RoP_vap path
                        // (byte-identical to Step 11). At the beam
                        // centre this matches the local Q above; off
                        // centre it incorrectly applies the centre
                        // value uniformly. damage_law tests
                        // (hu_end_to_end etc.) lock this behaviour in.
                        if (std::isfinite(L_tot) && L_tot > 0.0 && rho_p > 0.0)
                            rop_v = oneR * P_now / (rho_p * L_tot * A_beam);
                        if (!std::isfinite(rop_v) || rop_v < 0.0) rop_v = 0.0;
                    }
                    if (H_top >= H_vap && rop_v > 0.0)
                    {
                        vapor_fires_col[idx] = 1;
                        h_v_cand_col[idx]    = rop_v * dt_local;
                        rop_v_cand_col[idx]  = rop_v;
                    }
                }
            }
        }
        // P2-2e: the depth scan's and the lateral criterion's share of the
        // constraint-reference witnesses, added to what the Sp_field site
        // contributed this step. ⚠ `cref_on` is global, so every rank reaches
        // this collective or none does.
        if (cref_on)
        {
            long cnt[4] = {cref_vert, cref_lat, cref_clamp, cref_neg};
            amrex::ParallelDescriptor::ReduceLongSum(cnt, 4);
            cref_n_vert  += static_cast<Set::Scalar>(cnt[0]);
            cref_n_lat   += static_cast<Set::Scalar>(cnt[1]);
            cref_n_clamp += static_cast<Set::Scalar>(cnt[2]);
            cref_n_neg   += static_cast<Set::Scalar>(cnt[3]);
        }
        if (ncols > 0)
        {
            amrex::ParallelDescriptor::ReduceIntSum(spall_fires_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceIntSum(vapor_fires_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(h_s_cand_col.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(h_v_cand_col.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rop_s_cand_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rop_v_cand_col.data(), ncols);
        }
        if (rev_on && ncols > 0)
        {
            amrex::ParallelDescriptor::ReduceRealSum(rev_T_top.data(),  ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rev_a_f.data(),    ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rev_Sp.data(),     ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rev_h_scan.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(col_rho.data(),    ncols);
            amrex::ParallelDescriptor::ReduceRealSum(col_Cp.data(),     ncols);
        }

        // P0: collapse the owner-only reconstruction onto every rank BEFORE
        // the winner loop -- `flake` recession reads rc_h_f there, so this is
        // load-bearing, not a log.
        amrex::Vector<Set::Scalar> rc_dell;
        if (recon_on && ncols > 0)
        {
            amrex::ParallelDescriptor::ReduceRealSum(rc_Tsurf.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rc_ell.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rc_h_f.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rc_Tbar.data(),  ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rc_flag.data(),  ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rc_qrat.data(),  ncols);
            // ⚠⚠ P2-2d rev 6 DEFECT 1, ONE LINE. `rc_ellq` was filled
            // owner-only at :1610 and then written from the IO RANK'S LOCAL
            // COPY at :3651 while being absent from this block, so the CSV's
            // `ell_q` read true only where rank 0 happened to own the firing
            // column -- positive iff `col_j <= 11 && k_top <= 11`, on 100 % of
            // rows across four legs. A per-column diagnostic that is non-zero
            // on a GEOMETRIC SUBSET is a decomposition signature before it is
            // physics, and it is what made the flux lever look inert.
            // ⚠ The lever itself was never inert: `surf_ellq_col` is built in
            // BuildSurfaceTColumns AFTER the Ttop/T1/Ktop reduce, so it is
            // identical on every rank. Only the REPORTING was broken.
            amrex::ParallelDescriptor::ReduceRealSum(rc_ellq.data(), ncols);
            rc_dell.assign(ncols, 0.0);
            for (int idx = 0; idx < ncols; ++idx)
            {
                if (rc_ell[idx] <= 0.0) continue;
                if (rc_ell_prev[idx] > 0.0)
                    rc_dell[idx] = (rc_ell[idx] - rc_ell_prev[idx])
                                   / rc_ell_prev[idx];
                rc_ell_prev[idx] = rc_ell[idx];
            }
        }

        // Flake coherence-length cap (rock-mechanics cavity-wall coherence):
        // a firing column may recede at most coh_cells below its shallowest
        // live 4-neighbour. Bounds the cavity-wall slope so a smooth bowl can
        // descend gradually but a free-standing high-aspect-ratio needle cannot
        // form. Uses k_top_global (the gathered per-column live surface) so it
        // is decomposition-independent. k_top is the topmost live cell index
        // (higher k = shallower); a column at k_top removing thickness h lands
        // at k_top - h/dz, which must stay >= (shallowest-neighbour k_top) -
        // coh_cells.
        // Step 20: under the surface-normal kernel the coherence cap is
        // BYPASSED -- it is superseded, not just relaxed. The cap was a
        // needle regularizer for the vertical-only law; it also hard-limits
        // the cavity-wall slope to coh_len/dx, which is exactly what froze
        // the cone (every wall column cap-saturated against the cold rim).
        // Slope-aware normal motion is self-regularizing: steeper walls see
        // a thinner normal skin (fire more easily) AND recede further per
        // flake (h/cos(theta)), so slopes erode preferentially and needles
        // heal without a hard cap. Per-step removal stays bounded by the
        // scan window: h_vert = h_col/cos <= max_steps*dz = h_spall_max.
        const bool sn_active = (sn_cos_p != nullptr);
        const int coh_cells = (spall_flake_coh_len > 0.0 && !sn_active)
            ? std::max(1, static_cast<int>(std::lround(spall_flake_coh_len / dz)))
            : 0;
        if (coh_cells > 0)
        {
            for (int j = dom_jlo; j <= dom_jhi; ++j)
            for (int i = dom_ilo; i <= dom_ihi; ++i)
            {
                const int idx = col_index(i, j);
                if (spall_fires_col[idx] == 0 || h_s_cand_col[idx] <= 0.0) continue;
                const int kt = k_top_global[idx];
                if (kt == std::numeric_limits<int>::min()) continue;
                int nbr_max = std::numeric_limits<int>::min();
                const int nb[4][2] = {{i-1,j},{i+1,j},{i,j-1},{i,j+1}};
                for (int n = 0; n < 4; ++n)
                {
                    const int ni = nb[n][0], nj = nb[n][1];
                    if (ni < dom_ilo || ni > dom_ihi || nj < dom_jlo || nj > dom_jhi)
                        continue;
                    const int nkt = k_top_global[col_index(ni, nj)];
                    if (nkt > nbr_max) nbr_max = nkt;
                }
                if (nbr_max == std::numeric_limits<int>::min()) continue;
                const int floor_k = nbr_max - coh_cells;     // deepest allowed surface
                Set::Scalar cap_h = static_cast<Set::Scalar>(kt - floor_k) * dz;
                if (cap_h < 0.0) cap_h = 0.0;
                if (h_s_cand_col[idx] > cap_h)
                {
                    h_s_cand_col[idx]   = cap_h;
                    rop_s_cand_col[idx] = (dt_local > 0.0 && cap_h > 0.0)
                                        ? cap_h / dt_local : 0.0;
                    if (cap_h <= 0.0) spall_fires_col[idx] = 0;
                }
            }
        }

        // Step 16c: rank-ordered append of per-column first-firing events
        // to <plot_file>_h_col_events.csv. Each rank writes its own
        // subset; barriers serialize the file writes so the header (from
        // the I/O proc on the first emission) and subsequent rows don't
        // interleave. Append mode preserves cross-step accumulation; the
        // I/O proc latch ensures the header is emitted exactly once per
        // run (re-runs of the same plot_file truncate-and-rewrite).
        if (sp_weibull_spall && SurfaceOwner(lev) && sp_emit_h_col_events_csv)
        {
            const std::string csv_path =
                plot_file + "_h_col_events.csv";
            const int io_rank =
                amrex::ParallelDescriptor::IOProcessorNumber();
            const int my_rank = amrex::ParallelDescriptor::MyProc();
            const int n_ranks = amrex::ParallelDescriptor::NProcs();
            for (int r = 0; r < n_ranks; ++r)
            {
                if (my_rank == r)
                {
                    const bool need_header =
                        (r == io_rank) &&
                        !sp_h_col_events_csv_header_written;
                    std::ios::openmode mode = need_header
                                                  ? std::ios::trunc
                                                  : std::ios::app;
                    std::ofstream ofs(csv_path,
                                      mode | std::ios::out);
                    if (ofs.is_open())
                    {
                        if (need_header)
                        {
                            ofs << "time,col_i,col_j,a_f,h_col,"
                                   "is_gb_at_top,Sp_field_at_top\n";
                            sp_h_col_events_csv_header_written = true;
                        }
                        ofs.setf(std::ios::scientific);
                        ofs.precision(8);
                        for (auto const& ev : local_h_col_events)
                        {
                            ofs << ev.t << "," << ev.col_i << ","
                                << ev.col_j << "," << ev.a_f << ","
                                << ev.h_col << "," << ev.is_gb << ","
                                << ev.sp << "\n";
                        }
                    }
                }
                amrex::ParallelDescriptor::Barrier();
            }
        }

        // Step 19: mechanical bit feed. Within bit_radius of the axis, force
        // every column's surface down to the descending flat bit floor
        // z_bit(t) = bit_z0 - feed*(t - start), quantized to the cell whose
        // index is target_k. The recession needed is (k_top - target_k) cells.
        // NOT coherence-capped (the rigid bit cuts a flat floor); folded into
        // the winner loop below as max(spall, vapor, bit), so spallation may
        // still remove material AHEAD of (below) the bit. Uses k_top_global so
        // it is decomposition-independent.
        amrex::Vector<Set::Scalar> h_bit_col(ncols, 0.0);
        if (bit_enabled && time >= bit_start_time)
        {
            const Set::Scalar z_bit  = bit_z0 - bit_feed_rate * (time - bit_start_time);
            const int target_k = static_cast<int>(std::floor((z_bit - PLO[2]) / dz));
            const Set::Scalar r2_max = bit_radius * bit_radius;
            for (int j = dom_jlo; j <= dom_jhi; ++j)
            for (int i = dom_ilo; i <= dom_ihi; ++i)
            {
                const int idx = col_index(i, j);
                const int kt  = k_top_global[idx];
                if (kt == std::numeric_limits<int>::min()) continue;
                if (kt <= target_k) continue;            // already at/below bit
                // GATE: the bit advances only into spall-ready rock (the
                // spallation criterion fired this column this step). If the
                // beam has not heated the floor to Sp >= 1, the bit waits, so
                // ROP falls below the feed -- making ROP an emergent
                // prediction rather than an assumption. vapor counts too.
                if (spall_fires_col[idx] == 0 && vapor_fires_col[idx] == 0) continue;
                const Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const Set::Scalar r2 = (xc - bit_x0) * (xc - bit_x0)
                                     + (yc - bit_y0) * (yc - bit_y0);
                if (r2 > r2_max) continue;
                // Rate cap: the bit floor descends at the feed, so the
                // recession is bounded by (k_top - target_k) cells even when
                // the rock could spall faster -- the drawworks controls ROP.
                h_bit_col[idx] = static_cast<Set::Scalar>(kt - target_k) * dz;
            }
        }

        // Per-column winner from reduced data (deterministic on every rank).
        // Step 17: under sp_weibull, spall has priority over vap when both
        // gates fire at the same column. Under damage_law, the
        // winner-takes-all on h (Step 11 behaviour) is preserved
        // byte-identical via the else-branch.
        amrex::Vector<int> regime_col(ncols, 0);
        amrex::Vector<Set::Scalar> h_chosen_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> rop_chosen_col(ncols, 0.0);
        for (int idx = 0; idx < ncols; ++idx)
        {
            const bool sf = spall_fires_col[idx] != 0;
            const bool vf = vapor_fires_col[idx] != 0;
            const Set::Scalar hs = h_s_cand_col[idx];
            const Set::Scalar hv = h_v_cand_col[idx];
            int reg_choice = 0;
            Set::Scalar h_c = 0.0;
            if (sp_weibull_spall)
            {
                // Spall priority at the same (x, y) per revised-approach
                // §5: v_n = v_spall if Sp fires, else v_vap if H >= H_vap.
                if      (sf) { reg_choice = 1; h_c = hs; }
                else if (vf) { reg_choice = 2; h_c = hv; }
            }
            else
            {
                if (sf && vf)
                {
                    if (hv > hs) { reg_choice = 2; h_c = hv; }
                    else         { reg_choice = 1; h_c = hs; }
                }
                else if (sf) { reg_choice = 1; h_c = hs; }
                else if (vf) { reg_choice = 2; h_c = hv; }
            }
            // Step 19: when the bit is enabled it CONTROLS the surface (the
            // drawworks sets the rate; the flame rides with the bit, so the
            // surface cannot run ahead of it). The bit is the sole remover;
            // spallation/vapor remain as DIAGNOSTICS (Sp_field, computed
            // elsewhere) confirming the floor is spall-ready at this feed.
            if (bit_enabled)
            {
                reg_choice = (h_bit_col[idx] > 0.0) ? 3 : 0;
                h_c        = h_bit_col[idx];
            }
            // ⚠⚠ P1f: REGIME 6 -- a LATERAL firing. h_c is deliberately ZERO:
            // phi recession voids a column from its TOP, and a lateral firing must
            // void a MID-COLUMN cell (its riser is 75-150 cells tall with the hot
            // part at the bottom, so removing the top cell would take rock 150-300
            // mm away from the heat). The cell is voided explicitly through the
            // per-cell removed_mf path below -- D2q-1 made removed_mf
            // authoritative, and ColumnSolidRunBelowTop already exists precisely so
            // the depth scan does not read through an interior void as if it were
            // rock. ⚠ The interior cavity this leaves is resolved by COLLAPSE,
            // which is P1g; see the takeaways.
            // ⚠ PRECEDENCE: vertical > support (5) > lateral (6). reg_choice == 0
            // here means neither of the first two fired on this column.
            if (sp_lateral_removal && !lat_kfire.empty()
                && lat_kfire[idx] != std::numeric_limits<int>::min()
                && reg_choice == 0 && !bit_enabled)
            {
                reg_choice = 6;
                h_c        = 0.0;   // ⚠ NOT dz: phi must not move for regime 6
            }
            // D2z-0: retire the half-cell increment. Under computed-T the
            // depth scan's surface region is a lerp between T_face and T_top
            // with T_face = T_top (the documented fallback above), so T --
            // hence sigma and K_I -- is FLAT across the whole top half-cell.
            // At a marginal firing (Sp - 1 median 0.003-0.006) the crossing
            // therefore lands just past dz/2 and the column fires TWICE per
            // cell (measured 1.977 firings/void). full_cell voids the top cell
            // in one firing instead.
            //
            // SPALL REGIME ONLY. Vapor sets its own h_v from the enthalpy
            // budget and the bit is drawworks-set; neither carries the
            // artefact, and overriding them would change what fires rather
            // than when it is voided.
            //
            // prescribed_T is EXCLUDED: Step 16d gives that path a real face
            // temperature, so the flat half-cell does not arise and the scan
            // stays as it is.
            //
            // D2z-1: UNCONDITIONAL. There is no half-cell mode any more, so an
            // event in removal_events.csv IS a removal -- no n_voided > 0.5
            // filter, no twice-per-cell caveat. (That caveat had already
            // produced one wrong estimator: a step-height phase built on firing
            // intervals read 0.716 where the void-only answer was 0.034.)
            //
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
            // ⚠ sp_weibull_spall IS LOAD-BEARING, not decoration. The half-cell
            // artefact is a property of the sp_weibull DEPTH SCAN (the K_I
            // bisection landing just past dz/2). The damage-law /
            // spall.prescribed_ratio path computes h_s from a calibrated ratio
            // instead, where a multi-cell thickness is the intended answer --
            // unit/spall_event deliberately removes TWO cells per event at
            // h_spall = 2*dz. Without this term the override flattens that to
            // one cell and the test fails on five separate assertions. D2z-0
            // never saw it because the key was off on every input but its own.
            //
            // P0b (spall.removal_energy = flake): recede by the PHYSICAL flake
            // thickness h_f = min(ell_eff, dz) instead of forcing a whole cell,
            // and let the EXISTING mechanism void the cell when phi crosses
            // zero. ⚠ This deliberately reopens D2z-0/-1, and the record
            // matters: those packets retired the sub-cell increment for a
            // PLUMBING reason (firings_per_void 1.957 -> exactly 1.000,
            // n_voided singleton over 79,308 firings), with a rate effect of
            // +0.35/+0.44 %. So firings_per_void > 1 and h/dz < 1 now RETURN BY
            // DESIGN -- firing twice to shed a 1 mm flake from a 2 mm cell is
            // correct -- and the assertion below is restated accordingly.
            // ⚠ NOT an enthalpy transfer. You cannot conserve the rock's heat
            // while deleting the rock's geometry: the un-flaked remainder has
            // nowhere to be, and dumping it into the neighbour raises that cell
            // by ~dTbar, fires it next step and cascades one cell per timestep.
            // Here the heat leaves exactly when the geometry does, at the void.
            // Falls back to dz when the column has no reconstruction (no solid
            // second cell), so a firing always removes something.
            if (sp_weibull_spall && reg_choice == 1 && !sp_weibull_presc_T_face)
                h_c = (sp_removal_energy_flake && sp_surface_T_apply &&
                       !rc_h_f.empty() && rc_h_f[idx] > 0.0)
                      ? rc_h_f[idx] : dz;
            regime_col[idx]     = reg_choice;
            h_chosen_col[idx]   = h_c;
            rop_chosen_col[idx] = (dt_local > 0.0) ? h_c / dt_local : 0.0;
        }

        // P0b deliverable 4: build each firing column's remapped H profile and
        // its void count, from the REDUCED per-column buffers, so every rank
        // computes the same thing. H_mf has 0 ghost cells; PASS 3 writes each
        // cell from this array and never from a neighbour read.
        amrex::Vector<Set::Scalar> H_new_col;
        amrex::Vector<int>         nvoid_pred(remap_on ? ncols : 0, 0);
        Set::Scalar E_flake_step = 0.0;
        bool remap_any = false;
        // rev7 d5: per-rank, because each rank owns different cells -- unlike
        // E_flake_step, which every rank computes identically from reduced
        // per-column buffers. Reduced once after PASS 3.
        Set::Scalar E_bnd_local = 0.0;
        // ⚠⚠ P1b: THE SUPPORT REMOVAL'S OWN SINK. ledger_E_flake is what left
        // the SOLID set with the rock, and it was accumulated for regime 1
        // only, so every regime-5 void took a whole cell's enthalpy out of the
        // books unrecorded. Measured before the fix: ledger_err_solid ran to
        // -0.258 against P1's 7.6e-4, i.e. -8.95 J per voided cell against the
        // 8.69 J one cell holds at 500 K of superheat -- the defect was exactly
        // one cell, every time. It is a DIAGNOSTIC sum and feeds back into no
        // physics, which is why the fixed binary must reproduce the event log
        // byte for byte; that reproduction is the witness that this is
        // bookkeeping and not a change of behaviour.
        Set::Scalar E_sup_local = 0.0;
        // ⚠⚠ P1f: REGIME 6 OWES ITS OWN LEDGER TERM. Regime 5's omission took a
        // cell's enthalpy out of the books unrecorded and ran ledger_err_solid to
        // -0.258, pinned at -8.95 J per voided cell. (e) is the witness.
        Set::Scalar E_lat_local = 0.0;
        long lat_voided = 0;
        const Set::Scalar Vcell = dz * DX[0] * DX[1];
        // P1g-1: banking tallies. lat_void_col_local is PER COLUMN because the
        // cascade witness is "voids per step per column must never exceed 1" --
        // a domain-wide count cannot see a single column running away.
        long        lat_bank_local    = 0;
        Set::Scalar lat_bank_hf_local = 0.0;
        amrex::Vector<int> lat_void_col_local;
        if (ncols > 0) lat_void_col_local.assign(ncols, 0);
        // ⚠⚠ P1g-1: THE BANKING DECISION IS MADE HERE, NOT IN THE MFIter, AND
        // THAT IS WHAT MAKES (b4) TRUE BY CONSTRUCTION RATHER THAN BY LUCK.
        // `regime_col` (:2449) and `lat_hf` (reduced at :863) are already
        // rank-consistent, so every rank can run this identical per-column loop
        // and reach the same accumulator state with NO reduction at all. Doing
        // the arithmetic inside the MFIter instead would make the result depend
        // on which rank owned the firing cell, i.e. on the decomposition.
        // ⚠ Consequently lat_bank_local / lat_bank_hf_local must NOT be reduced:
        // every rank computes the same value, so a sum would multiply by the
        // rank count. Only lat_voided and E_lat_local, which accumulate on the
        // OWNING rank inside the MFIter, still need their reductions.
        amrex::Vector<char> lat_will_void;
        amrex::Vector<Set::Scalar> lat_bank_step;
        if (bank_on && ncols > 0)
        {
            if (static_cast<int>(lat_accum_col.size()) != ncols)
            {
                lat_accum_col.assign(ncols, 0.0);
                lat_nfire_col.assign(ncols, 0.0);
            }
            lat_will_void.assign(ncols, 0);
            lat_bank_step.assign(ncols, 0.0);
            const Set::Scalar eps = 1.0e-12 * DX[0];
            for (int idx = 0; idx < ncols; ++idx)
            {
                if (regime_col[idx] != 6) continue;
                const Set::Scalar hf =
                    lat_hf.empty() ? Set::Scalar(0.0) : lat_hf[idx];
                if (!(hf > 0.0)) continue;
                Set::Scalar acc = lat_accum_col[idx] + hf;
                lat_nfire_col[idx] += 1.0;
                ++lat_bank_local;
                lat_bank_hf_local += hf;
                if (acc >= DX[0] - eps)
                {
                    acc -= DX[0];
                    if (acc < 0.0) acc = 0.0;       // only ever by <= eps
                    lat_will_void[idx] = 1;
                }
                lat_accum_col[idx] = acc;
                lat_acc_log[idx] = acc;
                lat_nf_log[idx]  = lat_nfire_col[idx];
                lat_bank_step[idx] = hf;   // P1g-2: what the shift must displace
            }
        }
        if (remap_on)
        {
            H_new_col.assign(static_cast<std::size_t>(ncols) * colw, 0.0);
            for (int idx = 0; idx < ncols; ++idx)
            {
                if (regime_col[idx] != 1) continue;
                const Set::Scalar h_f = h_chosen_col[idx];
                const Set::Scalar ell = (idx < static_cast<int>(surf_ell_col.size()))
                                        ? surf_ell_col[idx] : 0.0;
                if (!(h_f > 0.0) || !(ell > 0.0)) continue;
                remap_any = true;
                const Set::Scalar x = dz / ell;
                const int m = static_cast<int>(std::floor(h_f / dz));
                const Set::Scalar r = h_f - Set::Scalar(m) * dz;
                // b in [0, 1): the fraction of the NEXT cell that the shift
                // pulls in, weighted by the exponential sub-cell profile.
                const Set::Scalar den = 1.0 - std::exp(-x);
                const Set::Scalar b = (den > 0.0)
                    ? (1.0 - std::exp(-r / ell)) / den : (r / dz);
                Set::Scalar* Hn = H_new_col.data()
                                + static_cast<std::size_t>(idx) * colw;
                const Set::Scalar* Ho = H_col_prof.data()
                                + static_cast<std::size_t>(idx) * colw;
                // P0b rev5 deliverable 6 -- THE (e) FIX. Revision 2 ALSO
                // remapped T with these weights and wrote it straight into
                // temp_mf. Two things went wrong with that, and together they
                // are the decomposition dependence:
                //   1. T is not a linear function of H. blend(T) and
                //      Invert(blend(H)) differ, so the pair (H, T) left the
                //      cell in a state that satisfies no equation of state --
                //      the only place in the code where T != Invert(H).
                //   2. temp_mf CARRIES GHOST CELLS and H_mf does not. Writing
                //      T into valid cells left every ghost copy of a remapped
                //      cell stale until the next FillPatch, so the conduction
                //      kernel read the PRE-remap temperature across a lateral
                //      box seam -- and only across a lateral seam, which is
                //      exactly the layout (mgs = 6) that differed while every
                //      layout keeping x and y whole agreed.
                // So: remap H, the conserved state, and NOTHING else. T is
                // recovered by the ordinary InvertHtoT below and its ghosts
                // are refilled, which is what every other H write in this
                // integrator already does.
                // ⚠ Ambient enters at depth: indices past the gathered window
                // clamp to its deepest cell, which sits 9-19 decay lengths down
                // (spall.h_spall_max = 0.015 m) and is at ambient.
                auto Hat = [&] (int d) -> Set::Scalar {
                    if (d < 0)         d = 0;
                    if (d > Dmax_col)  d = Dmax_col;
                    return Ho[d];
                };
                for (int d = 0; d <= Dmax_col; ++d)
                    Hn[d] = (1.0 - b) * Hat(d + m) + b * Hat(d + m + 1);
                // The energy the flake carries away: the telescoped difference,
                // which IS rho*c*h_f*(T_bar_f - T0) for an exponential profile.
                Set::Scalar dE = 0.0;
                for (int d = 0; d <= Dmax_col; ++d) dE += Ho[d] - Hn[d];
                E_flake_step += dE * dz * DX[0] * DX[1];
                // ⚠ h_f > dz is just several crossings -- no clamp. A cell at
                // depth index d holds phi = phi_top + d*dz and crosses when
                // phi - h_f < 0, so n = ceil((h_f - phi_top)/dz).
                const Set::Scalar over = h_f - phi_top_col[idx];
                nvoid_pred[idx] = (over > 0.0)
                    ? static_cast<int>(std::ceil(over / dz)) : 0;
            }
        }

        // C1: record each spall firing's top-cell T as the column's pin
        // (identical on every rank: all inputs are reduced per-column buffers).
        if (pin_on)
        {
            if (pin_T.size() != static_cast<std::size_t>(ncols))
            {
                pin_T.assign(ncols, 0.0);
                pin_t.assign(ncols, 0.0);
                pin_rho.assign(ncols, 0.0);
                pin_Cp.assign(ncols, 0.0);
            }
            for (int idx = 0; idx < ncols; ++idx)
            {
                if (regime_col[idx] != 1) continue;
                pin_T[idx]     = rev_T_top[idx];
                pin_t[idx]     = time;
                pin_rho[idx]   = col_rho.empty() ? 0.0 : col_rho[idx];
                pin_Cp[idx]    = col_Cp.empty()  ? 0.0 : col_Cp[idx];
            }
        }

        // PASS 3: every rank applies the per-column phi shift to its
        // local cells and writes per-cell detach diagnostics. Both
        // tasks read from the reduced per-column buffers, so the
        // result is identical regardless of how the BoxArray splits z.
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar>       const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& rem = (*removed_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& reg = (*regime_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& ropU= (*RoP_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& T   = (*temp_mf[lev]).array(mfi);
            (void) T;   // rev5: the remap writes H only; T comes back from InvertHtoT
            amrex::Array4<Set::Scalar>       const& LS  = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LL  = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LV  = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);

            amrex::Array4<Set::Scalar> ev_s, thk_s, rop_s;
            if (spall_on)
            {
                ev_s  = (*spall_event_mf[lev]).array(mfi);
                thk_s = (*spall_thickness_mf[lev]).array(mfi);
                rop_s = (*RoP_spall_mf[lev]).array(mfi);
            }
            amrex::Array4<Set::Scalar> ev_v, thk_v, rop_v;
            if (vapor_on)
            {
                ev_v  = (*vapor_event_mf[lev]).array(mfi);
                thk_v = (*vapor_thickness_mf[lev]).array(mfi);
                rop_v = (*RoP_vap_mf[lev]).array(mfi);
            }
            // Step 17: per-cell regime field under sp_weibull. Registered
            // only when sp_weibull_spall (gated above in Parse()), so
            // binding is conditional too.
            amrex::Array4<Set::Scalar> regf;
            if (sp_weibull_spall)
                regf = (*regime_field_mf[lev]).array(mfi);

            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);

            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                const int idx = col_index(i, j);
                const int reg_choice = regime_col[idx];
                if (reg_choice == 0) continue;

                const Set::Scalar h_c     = h_chosen_col[idx];
                const Set::Scalar rop_c   = rop_chosen_col[idx];
                const bool sf             = spall_fires_col[idx] != 0;
                const bool vf             = vapor_fires_col[idx] != 0;
                const Set::Scalar h_s_c   = h_s_cand_col[idx];
                const Set::Scalar h_v_c   = h_v_cand_col[idx];
                const Set::Scalar rop_s_c = rop_s_cand_col[idx];
                const Set::Scalar rop_v_c = rop_v_cand_col[idx];

                for (int k = klo; k <= khi; ++k)
                {
                    Set::Scalar phi_old = phi(i, j, k);
                    Set::Scalar phi_new = phi_old - h_c;
                    phi(i, j, k) = phi_new;
                    if (rem(i, j, k) > 0.5) continue;
                    if (phi_old >= 0.0 && phi_new < 0.0)
                    {
                        if (rev_on) ++rev_nvoid[idx];
                        if (inc_on) ++inc_nvoid[idx];
                        rem(i, j, k)  = 1.0;
                        reg(i, j, k)  = static_cast<Set::Scalar>(reg_choice);
                        // P1b: read H BEFORE anything else can touch it. An
                        // unsupported cell is hot through, so the cell's own
                        // enthalpy is the whole sink -- there is no flake
                        // remainder to remap and no telescoped sum to take.
                        if (reg_choice == 5) E_sup_local += H(i, j, k) * Vcell;
                        ropU(i, j, k) = rop_c;
                        // Step 17: also write the cell-resolved
                        // regime_field under sp_weibull. Same value as
                        // the scalar `reg(i, j, k)` here — the
                        // distinction is the registered field name
                        // (`regime` vs `regime_field`) so plotfile
                        // consumers can pick the canonical per-cell
                        // diagnostic under sp_weibull.
                        if (sp_weibull_spall)
                            regf(i, j, k) =
                                static_cast<Set::Scalar>(reg_choice);
                        if (spall_on)
                        {
                            ev_s (i, j, k) = sf ? 1.0 : 0.0;
                            thk_s(i, j, k) = sf ? h_s_c : 0.0;
                            rop_s(i, j, k) = rop_s_c;
                        }
                        if (vapor_on)
                        {
                            ev_v (i, j, k) = vf ? 1.0 : 0.0;
                            thk_v(i, j, k) = vf ? h_v_c : 0.0;
                            rop_v(i, j, k) = rop_v_c;
                        }
                        D(i, j, k) = 0.0;
                        // Do NOT reset the voided cell's thermal state to
                        // ambient. The cell now holds no material and is
                        // insulating (the void-aware conduction stencil carries
                        // zero flux across any void face and zeroes the void's
                        // own conductivity), so its stored T/H/latent state is
                        // physically inert. Resetting it to T_amb only created a
                        // spurious cold (280 K) value that showed up as a
                        // checkerboard in temperature plotfiles next to the
                        // un-spalled hot cells at the rough drill front. Keeping
                        // the prior thermal state leaves the temperature field
                        // continuous; the cavity is delineated for visualization
                        // by the `removed` mask (threshold removed < 0.5 to
                        // carve the hole). NaN-masking is avoided because
                        // UpdateSpAfterMechanics reads the absolute domain-top
                        // cell, which becomes a void after drilling — a NaN
                        // there would poison the Sp = K_I/K_Ic gate.
                    }
                }

                // ⚠⚠ P1f: REGIME 6 -- void EXACTLY the cell the lateral
                // criterion fired, and nothing else. Not through phi.
                // ⚠ n_voided stays 1, so D2z-1's invariant is intact and needs no
                // relaxation -- that is the whole reason collapse is P1g.
                if (reg_choice == 6 && !lat_kfire.empty())
                {
                    const int kf = lat_kfire[idx];
                    if (kf >= klo && kf <= khi && rem(i, j, kf) <= 0.5)
                    {
                        // P1g-1 deliverables 2-4: BANK THE THICKNESS, then void
                        // only once a whole cell's worth has accumulated.
                        //
                        // ⚠⚠ GEOMETRY BANKS; ENTHALPY LEAVES ONLY ON VOID. The
                        // accumulator takes NO enthalpy out of a cell that does
                        // not void: `H` is untouched on the banking path, which is
                        // what (b3) asserts to 1e-14 against the OFF leg. The
                        // alternative -- subtracting a fraction of the cell's
                        // enthalpy sub-cell -- is KNOWN WRONG: it reads `ell` back
                        // 1.32x too long at dz/ell = 2 and 3.81x at 1, corrupting
                        // the very two cells the lateral criterion reads.
                        //
                        // ⚠ The honest cost: a banked cell holds up to
                        // dx/h_f_lat = 1.33 flakes' worth of enthalpy it has
                        // nominally shed. That is a BOUNDED approximation, and it
                        // is not a loss of accuracy that previously existed --
                        // `phi` is not read by the conduction kernel either, so
                        // this step does not claim a sub-cell thermal frame.
                        if (bank_on)
                        {
                            // ⚠⚠ THE REMAINDER CARRIES PAST THE VOID -- it lives
                            // with the recession TRACK (the column), not the cell.
                            // A voided cell is skipped by the criterion forever
                            // (it tests rm > 0.5 first), so a per-cell remainder
                            // could never be spent and the realised law would read
                            // 0.5 cells per firing instead of 0.75. The decision
                            // itself was taken in the rank-invariant pre-loop
                            // above; this branch only EXECUTES it.
                            if (!lat_will_void.empty() && lat_will_void[idx])
                            {
                                rem(i, j, kf)  = 1.0;
                                reg(i, j, kf)  = 6.0;
                                ropU(i, j, kf) = 0.0;
                                if (sp_weibull_spall) regf(i, j, kf) = 6.0;
                                E_lat_local += H(i, j, kf) * Vcell;
                                if (rev_on) ++rev_nvoid[idx];
                                if (inc_on) ++inc_nvoid[idx];
                                ++lat_voided;
                                ++lat_void_col_local[idx];
                            }
                        }
                        else
                        {
                            // The shipped whole-cell law, unchanged.
                            rem(i, j, kf)  = 1.0;
                            reg(i, j, kf)  = 6.0;
                            ropU(i, j, kf) = 0.0;
                            if (sp_weibull_spall) regf(i, j, kf) = 6.0;
                            // the sink is the cell's own enthalpy, as regime 5 does
                            E_lat_local += H(i, j, kf) * Vcell;
                            if (rev_on) ++rev_nvoid[idx];
                            if (inc_on) ++inc_nvoid[idx];
                            ++lat_voided;
                            ++lat_void_col_local[idx];
                        }
                    }
                }
                // P0b deliverable 4: write the remapped profile. SECOND loop,
                // so `rem` is settled -- a cell that flipped above is void and
                // must keep its frozen state.
                // ⚠ THE INDEX SHIFT IS THE WHOLE BOOKKEEPING. After n cells
                // void, the new top cell is the one that was at depth index n,
                // and it must hold the profile's depth 0. So the cell at OLD
                // depth index jd takes new[jd - n].
                if (remap_on && reg_choice == 1 && !H_new_col.empty())
                {
                    const int n_v  = nvoid_pred[idx];
                    const int ktop = k_top_global[idx];
                    const std::size_t base =
                        static_cast<std::size_t>(idx) * colw;
                    for (int k = klo; k <= khi; ++k)
                    {
                        if (rem(i, j, k) > 0.5) continue;
                        const int jd = ktop - k;
                        const int d  = jd - n_v;
                        if (jd < 0 || d < 0 || d > Dmax_col) continue;
                        // ⚠⚠ REVISION 7 DELIVERABLE 5 -- THE WINDOW-BOUNDARY
                        // TERM. E_flake telescopes Ho - Hn over the GATHERED
                        // window only (depths 0..Dmax). The solid set also
                        // loses the n cells at depths Dmax+1..Dmax+n, which
                        // the shift moves up into the window and which no
                        // telescoped sum inside it can see. That term is zero
                        // only if the rock at h_spall_max is at ambient -- true
                        // on the 1-D column (exp(-15/1.3) = 1e-5, and the
                        // ledger closed to 4e-4 there) and false on Meier,
                        // where the hole wall heats the rock laterally and it
                        // sits ~82 K up, which is the whole of the -0.16
                        // residual. Read the OLD H, before the overwrite.
                        if (jd > Dmax_col) E_bnd_local += H(i, j, k) * Vcell;
                        H(i, j, k) = H_new_col[base + d];
                    }
                }
            }
        }

        // D2z-0 Goal 2: h = dz lands EXACTLY on the phi tie. The top cell holds
        // phi_top in [0, dz), so phi_top - dz is in [-dz, 0) and flips, while
        // the cell below holds phi in [dz, 2dz) and lands in [0, dz), which is
        // NOT < 0 and must not flip. That is a one-ulp question at phi_top = 0,
        // and Removal.H:145-148 records that ties here "miss or add one".
        // So this is an ASSERTION, not a diagnostic: a whole-cell spall firing
        // that voids anything other than exactly one cell is a hard stop.
        // D2z-0 ran 79,308 such firings across two meshes with no violation;
        // D2z-1 keeps the assertion permanently because the tie is structural,
        // not a property of that campaign.
        if (inc_on && ncols > 0)
        {
            // Collective: ReduceIntSum deposits the per-column totals on the
            // IO rank. Thermo variables registered `extensive = false` are NOT
            // reduced by the Integrator (Integrator.cpp:1141-1144) and only the
            // IO rank writes thermo.dat, so the IO rank is where both the
            // assertion and the counters below must be evaluated.
            amrex::ParallelDescriptor::ReduceIntSum(inc_nvoid.data(), ncols);
            if (amrex::ParallelDescriptor::IOProcessor())
            {
                if (inc_active)
                {
                    // P1b: a support firing recedes by exactly dz, so D2z-1's
                    // original invariant applies to it UNCONDITIONALLY -- one
                    // firing, one cell -- with none of the `flake` sub-cell
                    // caveat below. This is witness 3: if the rule ever fired a
                    // column that had already fired, this is where it stops.
                    // ⚠ P1f: regime 6 voids exactly one cell, so D2z-1's
                    // invariant applies to it unchanged -- no relaxation, which is
                    // why collapse (n_voided > 1 by design) is a separate step.
                    if (lat_on)
                        for (int idx = 0; idx < ncols; ++idx)
                        {
                            if (regime_col[idx] != 6) continue;
                            // ⚠⚠ P1g-1 deliverable 5: RELAXED BY EXACTLY ONE CASE.
                            // Under banking a lateral firing legitimately voids
                            // ZERO cells -- it banks instead -- so the permitted
                            // set becomes {0, 1} and ONLY when the key is on.
                            // ⚠ The upper bound is deliberately untouched: 2 or
                            // more is still an abort on every leg, because
                            // multi-cell lateral removal is collapse (P1h) and
                            // must not arrive silently through this path.
                            // (a2) scores both sides -- banking runs with voids =
                            // 0, and a forced 2-cell lateral void still aborts.
                            const int lo_ok = bank_on ? 0 : 1;
                            if (inc_nvoid[idx] < lo_ok || inc_nvoid[idx] > 1)
                                Util::Abort(INFO, "spall.lateral_removal: a lateral "
                                            "firing at column ", idx, " voided ",
                                            inc_nvoid[idx], " cells, not ",
                                            (bank_on ? "0 or 1" : "exactly 1"),
                                            ". Regime 6 voids at most ONE mid-column "
                                            "cell through removed_mf and never moves "
                                            "phi; under spall.lateral_subcell it may "
                                            "instead bank and void none.");
                        }
                    for (int idx = 0; idx < ncols; ++idx)
                    {
                        if (regime_col[idx] != 1) continue;
                        // P0b: under `flake` a firing recedes phi by
                        // h_f <= dz, so a firing voids ZERO OR ONE cell and
                        // several firings may precede one void. The invariant
                        // that still holds -- and the one D2z-0's assertion was
                        // really protecting -- is "exactly one void per
                        // zero-crossing", i.e. never more than one. The strict
                        // form is kept on the default path.
                        // P0b: under the remap the invariant is "voids ==
                        // crossings": h_f is not clamped to dz, so one firing
                        // may cross several cell boundaries, and the predicted
                        // count is ceil((h_f - phi_top)/dz). Without the remap
                        // (P0's `flake`) it is 0 or 1; on the default path the
                        // strict D2z-0 form is kept.
                        const int nv_want = remap_on ? nvoid_pred[idx] : -1;
                        const int nv_min = sp_removal_energy_flake ? 0 : 1;
                        const bool bad = remap_on
                            ? (inc_nvoid[idx] != nv_want)
                            : (inc_nvoid[idx] > 1 || inc_nvoid[idx] < nv_min);
                        if (bad)
                        {
                            // Revision 7 deliverable 4. THE MESSAGE USED TO
                            // NAME ONE CAUSE AND THERE ARE TWO. Revision 5's
                            // Meier run aborted here because the column had
                            // drilled to within 6 cells of the domain floor
                            // and a 15 mm flake wanted 8 -- nothing to do with
                            // a tie at a cell boundary, which is what the
                            // message said, and the wrong lead cost time.
                            // A column can only void as many cells as it has
                            // above the floor, so compare the shortfall with
                            // what was actually available before blaming the
                            // sign test.
                            const int avail = (idx < static_cast<int>(k_top_global.size()))
                                ? (k_top_global[idx] - dom_klo_z + 1) : -1;
                            const bool out_of_rock =
                                remap_on && avail >= 0 && nv_want > avail;
                            if (out_of_rock)
                                Util::Abort(INFO,
                                    "spall removal RAN OUT OF ROCK at column (",
                                    idx % nx_dom + dom_ilo, ", ",
                                    idx / nx_dom + dom_jlo, ") at t = ", time,
                                    ": the flake asked for ", nv_want,
                                    " cells and the column has only ", avail,
                                    " above the domain floor (it voided ",
                                    inc_nvoid[idx], "). THE BLOCK IS TOO "
                                    "SHALLOW FOR THIS RUN -- deepen the domain "
                                    "or stop earlier. This is NOT the phi tie.");
                            Util::Abort(INFO,
                                "spall removal voided ",
                                inc_nvoid[idx], " cells at column (",
                                idx % nx_dom + dom_ilo, ", ",
                                idx / nx_dom + dom_jlo, ") at t = ", time,
                                " -- expected ",
                                (remap_on ? nv_want
                                          : (sp_removal_energy_flake ? 1 : 1)),
                                ", with ", avail, " cells above the domain "
                                "floor (so the column did NOT run out of rock). "
                                "This is the phi tie at "
                                "Removal.H:145-148: h = dz landed on a cell "
                                "boundary and the sign test missed or doubled "
                                "a flip.");
                        }
                    }
                }

                // D2z-1 Goal 4 (revised): firings_per_void is NOT a thermo
                // column. Under whole-cell its value is 1.000 by construction,
                // so it carries no information, and registering it appended
                // after the eight jet_* columns and broke jet_enthalpy's
                // criterion (a), which asserts those eight are the LAST eight.
                // The invariant is already enforced where it belongs -- the
                // n_voided == 1 assertion above, on every firing. Dropping the
                // column lets the tree audit be absolute: no column appears or
                // disappears, with no exception carved out for our own.
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
            }
        }

        // P0b (c''): the sink, for the domain ledger. Summed across ranks
        // because each column's contribution was computed on every rank from
        // reduced buffers -- so take it from the IO rank's copy only.
        if (remap_on)
        {
            amrex::ParallelDescriptor::ReduceRealSum(E_bnd_local);
            amrex::ParallelDescriptor::ReduceRealSum(E_sup_local);
            ledger_E_flake += E_flake_step + E_bnd_local + E_sup_local;
        }
        // ⚠⚠ P1f: THE LATERAL LEDGER TERM IS UNCONDITIONAL, NOT INSIDE `remap_on`.
        // Regime 6 voids a cell whatever spall.removal_energy is, so its enthalpy
        // leaves the solid set on a `cell`-energy leg too. Folding it into the
        // remap-only block above cost ledger_err_solid -1.0 on the unit test and
        // would have been INVISIBLE on both Meier legs, which run
        // removal_energy = flake and therefore take that branch -- i.e. it would
        // have PASSED (e) and been silently wrong on every cell-energy run.
        // The energy invariant is what found it.
        // ==================================================================
        // P1g-2 deliverables 4-10: THE LATERAL FLAKE SINK, THERMAL HALF
        // ==================================================================
        // ⚠⚠ WHY THIS STEP EXISTS. P1g-1's banking is geometrically exact but
        // PHYSICALLY INERT: on matched rungs both legs voided 56 cells at void
        // times agreeing to the digit, because a bank changes nothing the
        // conduction kernel can see. The vertical remap delays the next firing
        // by TAKING THE HOT SKIN AWAY so the surface must re-heat. This supplies
        // that half for the lateral face.
        //
        // ⚠⚠ THIS IS NEW INFRASTRUCTURE, NOT A MIRROR OF THE VERTICAL GATHER.
        // The vertical gather (:1147) is one vector per column over DEPTH, and
        // depth is the direction every column recedes in. Here the direction is
        // the FACE NORMAL and it differs per column, so this is a ROW AT FIXED
        // HEIGHT per in-plane direction, built only for the faces that banked
        // this step. It is reduced so every rank holds the same row, then
        // written back per box, which is how the vertical shift does it.
        //
        // ⚠⚠ WRITE ORDER, STATED (deliverable 8). The vertical remap writes a
        // COLUMN's enthalpy inside the removal MFIter above. This block runs
        // AFTER that loop has closed, so on any cell both writers touch in one
        // step, THE LATERAL SHIFT WINS. That precedence is deliberate: the
        // lateral face is the surface that just receded, so its profile is the
        // more recent physical statement about that cell. Overlaps are COUNTED
        // (`lat_overlap`) rather than assumed absent -- on this fixture the
        // count must be 0 because there are no vertical firings at all
        // (regimes {6: 80}), which is what makes the count a real witness
        // instead of a tautology.
        bool lat_shift_any = false;
        if (bank_on && sp_lateral_thermal && ncols > 0)
        {
            // ⚠ ROWMIN is deliverable 4's refusal rule. The shared fixture box
            // leaves only 3-4 cells behind the face, so without this the shift
            // would clamp at the row end and effectively pull material IN FROM
            // OUTSIDE THE DOMAIN. The vertical remap pulls COLD rock from below,
            // which is physical; a mirrored hot cell is not. (c4) is the witness
            // and the WIDE leg is the configuration where it passes.
            const int ROWW   = 32;
            const int ROWMIN = 6;
            const std::size_t RW = static_cast<std::size_t>(ROWW);
            amrex::Vector<Set::Scalar> rH(static_cast<std::size_t>(ncols) * RW * 2, 0.0);
            amrex::Vector<Set::Scalar> rOK(static_cast<std::size_t>(ncols) * RW * 2, 0.0);

            auto decode = [](Set::Scalar e) {
                return (e > 0.0) ? static_cast<int>(std::lround(e)) - 2 : 0;
            };

            // ---- gather the two rows (x then y) for every banking column ----
            for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& Hc =
                    (*H_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& rmc =
                    (*removed_mf[lev]).const_array(mfi);
                for (int idx = 0; idx < ncols; ++idx)
                {
                    if (!(lat_bank_step[idx] > 0.0)) continue;
                    const int kf = lat_kfire[idx];
                    if (kf < bx.smallEnd(2) || kf > bx.bigEnd(2)) continue;
                    const int ci = idx % nx_dom + dom_ilo;
                    const int cj = idx / nx_dom + dom_jlo;
                    const int sx = decode(lat_sx_enc[idx]);
                    const int sy = decode(lat_sy_enc[idx]);
                    for (int ax = 0; ax < 2; ++ax)
                    {
                        const int sgx = (ax == 0) ? sx : 0;
                        const int sgy = (ax == 0) ? 0  : sy;
                        if (sgx == 0 && sgy == 0) continue;
                        const std::size_t base =
                            (static_cast<std::size_t>(idx) * 2 + ax) * RW;
                        for (int m = 0; m < ROWW; ++m)
                        {
                            const int ii = ci + m * sgx, jj = cj + m * sgy;
                            if (ii < dom_ilo || ii > dom_ihi ||
                                jj < dom_jlo || jj > dom_jhi) break;
                            if (ii < bx.smallEnd(0) || ii > bx.bigEnd(0) ||
                                jj < bx.smallEnd(1) || jj > bx.bigEnd(1)) continue;
                            if (rmc(ii, jj, kf) > 0.5) continue;   // void: not row
                            rH [base + m] = Hc(ii, jj, kf);
                            rOK[base + m] = 1.0;
                        }
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(rH.data(),
                static_cast<int>(rH.size()));
            amrex::ParallelDescriptor::ReduceRealSum(rOK.data(),
                static_cast<int>(rOK.size()));

            // ---- build the shifted rows, identically on every rank ----------
            amrex::Vector<Set::Scalar> rNew(rH.size(), 0.0);
            amrex::Vector<char>        rDo(static_cast<std::size_t>(ncols) * 2, 0);
            Set::Scalar E_shift_local = 0.0, E_est_local = 0.0;
            Set::Scalar vol_disp = 0.0, vol_bank = 0.0;
            long n_shift = 0, n_refuse = 0;
            // P1g-2c deliverable 6: the three outcomes n_refuse conflated.
            // ⚠⚠ "refusal" read as a failure and sent a reader looking for a
            // broken shift; (a2) proved shifts + refusals == banks EXACTLY and
            // refusals/voids == 1.01, i.e. almost every "refusal" is simply a
            // bank that VOIDED its own face cell, so the row gathers length 0.
            // n_refuse itself is kept, unchanged, so lat_shift_refuse and the
            // identity shifts+refusals==banks both still hold.
            long n_ref_ell = 0, n_ref_void = 0, n_ref_short = 0;
            // P1g-2c deliverable 2: two-axis coverage. Per-COLUMN, so it must
            // NOT be rank-reduced (see the reduction note below).
            long n_twoaxis = 0;
            for (int idx = 0; idx < ncols; ++idx)
            {
                const Set::Scalar hb = lat_bank_step[idx];
                if (!(hb > 0.0)) continue;
                const Set::Scalar ell = lat_ell[idx];
                if (!(ell > 0.0)) { ++n_refuse; ++n_ref_ell; continue; }
                const Set::Scalar wgt[2] = {lat_wx[idx], lat_wy[idx]};
                // the contiguous solid run available in each direction
                int runlen[2] = {0, 0};
                for (int ax = 0; ax < 2; ++ax)
                {
                    const std::size_t base =
                        (static_cast<std::size_t>(idx) * 2 + ax) * RW;
                    int L = 0;
                    while (L < ROWW && rOK[base + L] > 0.5) ++L;
                    runlen[ax] = L;
                }
                // ⚠ refuse PER AXIS: a corner may have a long row one way and a
                // stub the other, and shifting only the long one would displace
                // less than it banked -- (b5) would then read below 1.
                bool refuse = false;
                for (int ax = 0; ax < 2; ++ax)
                    if (wgt[ax] > 0.0 && runlen[ax] < ROWMIN) refuse = true;
                if (refuse)
                {
                    ++n_refuse;
                    // P1g-2c deliverable 6: which kind of short row was it?
                    // runlen == 0 on an active axis means the face cell itself
                    // is already void -- this bank voided it -- and that is the
                    // ~25,000 case, not a failure of anything.
                    bool zero = false;
                    for (int ax = 0; ax < 2; ++ax)
                        if (wgt[ax] > 0.0 && runlen[ax] == 0) zero = true;
                    if (zero) ++n_ref_void; else ++n_ref_short;
                    continue;
                }
                // deliverable 2: counted BEFORE the per-axis loop, so it counts
                // faces (columns), not row-writes.
                if (wgt[0] > 0.0 && wgt[1] > 0.0) ++n_twoaxis;
                bool did = false;
                for (int ax = 0; ax < 2; ++ax)
                {
                    if (!(wgt[ax] > 0.0)) continue;
                    const int L = runlen[ax];
                    const std::size_t base =
                        (static_cast<std::size_t>(idx) * 2 + ax) * RW;
                    // ⚠⚠ cos^2 SPLIT. Each row is shifted by its OWN share
                    // wgt*hb, and the two shares sum to hb exactly, so the row
                    // material displaced is hb*A however the corner is oriented.
                    const Set::Scalar hsh = wgt[ax] * hb;
                    const Set::Scalar dl  = DX[(ax == 0) ? 0 : 1];
                    // the vertical remap's own exponential-weighted rule
                    // (:2538-2557), applied along the normal instead of depth
                    const Set::Scalar x   = dl / ell;
                    const int         mm  = static_cast<int>(std::floor(hsh / dl));
                    const Set::Scalar r   = hsh - Set::Scalar(mm) * dl;
                    const Set::Scalar den = 1.0 - std::exp(-x);
                    const Set::Scalar b   = (den > 0.0)
                        ? (1.0 - std::exp(-r / ell)) / den : (r / dl);
                    for (int m = 0; m < L; ++m)
                    {
                        // ⚠ clamp at the row end: indices past the gathered run
                        // hold the deepest cell, which the ROWMIN rule has
                        // guaranteed is at least 6 cells from the face.
                        const int i0 = std::min(m + mm,     L - 1);
                        const int i1 = std::min(m + mm + 1, L - 1);
                        rNew[base + m] = (1.0 - b) * rH[base + i0]
                                       +        b  * rH[base + i1];
                        E_shift_local += (rH[base + m] - rNew[base + m]) * Vcell;
                    }
                    // ⚠ THE FACE AREA OF A LATERAL FACE IS NOT dx*dy. For a
                    // row running along axis `ax` the face normal is `ax`, so
                    // the area is the product of the OTHER two spacings. On this
                    // fixture dx=dy=dz=2 mm so every choice is numerically
                    // 4e-6 m^2 and the error would be invisible -- which is
                    // exactly why it is worth writing correctly now.
                    const Set::Scalar A_face = (ax == 0)
                        ? DX[1] * DX[2] : DX[0] * DX[2];
                    vol_disp += hsh * A_face;
                    // ⚠⚠ deliverable 9, THE INDEPENDENT ESTIMATE, in the
                    // DISCRETE form the shift actually realises. With a sub-cell
                    // shift (mm = 0) the telescoped sum collapses exactly to
                    //      b * (Ho[0] - Ho[L-1]) * Vcell
                    // so the matching estimate is the GEOMETRIC series over
                    // cell centres, not a continuum integral:
                    //      sum_m exp(-(m+1/2) dl/ell) = e^{-dl/2ell}/(1-e^{-dl/ell})
                    // A continuum integral over [0,h] instead differs by
                    // dl/(ell(1-e^{-dl/ell})) = 1.31 here, and anchoring it at
                    // the face rather than the first cell centre differs by a
                    // further e^{-dl/2ell} = 0.75. Both were measured before this
                    // form was adopted; see the takeaway.
                    // ⚠ Built from the LOGGED face temperature and decay length
                    // and nothing else, so it stays independent of the shift's
                    // own sum, which is what makes (b3) not a tautology.
                    const Set::Scalar dTf = lat_Tf[idx] - phase_T_ref_loc;
                    // P2-2b site 3: rho and c_p read SEPARATELY from the
                    // per-column gather and multiplied at the point of use.
                    // ⚠ This loop has no MFIter, which is why the gather exists.
                    const Set::Scalar rho_c =
                        (col_rho.empty() || col_Cp.empty())
                        ? Set::Scalar(0.0) : col_rho[idx] * col_Cp[idx];
                    const Set::Scalar dd = 1.0 - std::exp(-dl / ell);
                    if (dTf > 0.0 && rho_c > 0.0 && dd > 0.0)
                        E_est_local += rho_c * dTf
                                     * (1.0 - std::exp(-hsh / ell))
                                     * std::exp(-0.5 * dl / ell) / dd
                                     * dl * A_face;
                    rDo[static_cast<std::size_t>(idx) * 2 + ax] = 1;
                    did = true;
                }
                if (did)
                {
                    ++n_shift;
                    // the banked volume uses the SAME face area convention as
                    // the displaced volume, so (b5)'s ratio is a pure geometry
                    // check and cannot be rescued by an area mismatch.
                    vol_bank += hb * ((lat_wx[idx] >= lat_wy[idx])
                                      ? DX[1] * DX[2] : DX[0] * DX[2]);
                    // ⚠ The estimate itself is accumulated PER AXIS above, in
                    // the DISCRETE geometric form. The continuum expression
                    //     rho_c (Tf - T0) ell (1 - e^{-h/ell}) A
                    // was measured and REJECTED (it reads 0.843 face-anchored and
                    // 1.152 cell-anchored); do not reinstate it here.
                    // ⚠ The LEDGER takes the telescoped sum, never the estimate:
                    // the ledger must equal the enthalpy the field actually lost
                    // or (b2) cannot close. P1g-2's takeaway 3 records that
                    // deliverable 9's wording and (b2)'s bar cannot both be met.
                    // ⚠⚠ AND THE ESTIMATOR ASSUMES DECAY TO AMBIENT, which is why
                    // (b3) misses 10 % at 14.2 % run-integrated: the P1g-2 review
                    // showed the per-step ratio DRIFTS 0.95 -> 0.78 as the far
                    // field soaks, and the first steps sit inside 4.7 %. A
                    // mis-chosen anchor would be a CONSTANT bias, not a drift.
                    // Any successor needs a substrate-base third parameter or a
                    // bar on the early window.
                }
            }

            // ---- write the rows back, per box -------------------------------
            // ⚠⚠ P1g-2c deliverable 1: THE TRUTH DIAGNOSTIC, AND WHY IT MUST
            // READ THE FIELD. lat_shift_ratio reads 1.0000 at p05/p50/p95 and
            // lat_shift_E_est/lat_shift_E reads 0.947, and BOTH ARE BLIND to the
            // overwrite because both are built from the INTENDED shift
            // (vol_disp, E_shift_local) rather than from H after the write. So
            // the realised energy is summed from Hw itself, before and after.
            // Anything derived from the row buffers would read 1.00 by
            // construction -- the tautology that let this through two reviews.
            //
            // The per-cell record also carries the touch count, which is what
            // deliverable 5's fix needs: a cell written by two rows must compose
            // its two deltas instead of keeping only the last.
            struct ShiftCell
            {
                int ii, jj, kf, cnt;
                Set::Scalar H_old, sum_delta, last_new;
            };
            // key packs (ii,jj,kf); the domain is <= a few hundred cells per
            // axis so a 16-bit field with a 32768 offset cannot collide.
            auto pack = [](int a, int b, int c) -> long long
            {
                const long long O = 32768, M = 65536;
                return ((static_cast<long long>(a) + O) * M
                        + (static_cast<long long>(b) + O)) * M
                       + (static_cast<long long>(c) + O);
            };
            Set::Scalar E_real_local = 0.0;
            long n_multi = 0;
            long n_overlap = 0;
            if (n_shift > 0)
            {
                lat_shift_any = true;
                for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    // ⚠ One map per BOX: a cell belongs to exactly one box
                    // (MFIter(..., false) + validbox), so two boxes can never
                    // hold the same cell and the counts cannot double-apply.
                    std::map<long long, ShiftCell> touched;
                    const amrex::Box& bx = mfi.validbox();
                    amrex::Array4<Set::Scalar> const& Hw = (*H_mf[lev]).array(mfi);
                    amrex::Array4<const Set::Scalar> const& rmc =
                        (*removed_mf[lev]).const_array(mfi);
                    for (int idx = 0; idx < ncols; ++idx)
                    {
                        if (!(lat_bank_step[idx] > 0.0)) continue;
                        const int kf = lat_kfire[idx];
                        if (kf < bx.smallEnd(2) || kf > bx.bigEnd(2)) continue;
                        const int ci = idx % nx_dom + dom_ilo;
                        const int cj = idx / nx_dom + dom_jlo;
                        const int sx = decode(lat_sx_enc[idx]);
                        const int sy = decode(lat_sy_enc[idx]);
                        for (int ax = 0; ax < 2; ++ax)
                        {
                            if (!rDo[static_cast<std::size_t>(idx) * 2 + ax])
                                continue;
                            const int sgx = (ax == 0) ? sx : 0;
                            const int sgy = (ax == 0) ? 0  : sy;
                            const std::size_t base =
                                (static_cast<std::size_t>(idx) * 2 + ax) * RW;
                            for (int m = 0; m < ROWW; ++m)
                            {
                                if (!(rOK[base + m] > 0.5)) break;
                                const int ii = ci + m * sgx, jj = cj + m * sgy;
                                if (ii < bx.smallEnd(0) || ii > bx.bigEnd(0) ||
                                    jj < bx.smallEnd(1) || jj > bx.bigEnd(1))
                                    continue;
                                if (rmc(ii, jj, kf) > 0.5) continue;
                                // deliverable 8: did the vertical remap write
                                // this same cell this step?
                                const int oc = col_index(ii, jj);
                                if (remap_on && oc >= 0 && oc < ncols
                                    && regime_col[oc] == 1) ++n_overlap;
                                // P1g-2c: record the touch BEFORE writing, so
                                // H_old is the field's own value and the
                                // diagnostic cannot be a tautology.
                                const long long key = pack(ii, jj, kf);
                                auto it = touched.find(key);
                                if (it == touched.end())
                                    touched.emplace(key, ShiftCell{
                                        ii, jj, kf, 1, Hw(ii, jj, kf),
                                        rNew[base + m] - Hw(ii, jj, kf),
                                        rNew[base + m]});
                                else
                                {
                                    ++it->second.cnt;
                                    it->second.sum_delta +=
                                        rNew[base + m] - it->second.H_old;
                                    it->second.last_new = rNew[base + m];
                                }
                                // ⚠⚠ P1g-2c deliverable 5: THE WRITE HAS MOVED
                                // OUT OF THIS LOOP. It used to be
                                //     Hw(ii, jj, kf) = rNew[base + m];
                                // right here -- an ASSIGNMENT. A face exposed on
                                // two sides runs two rows under the cos^2 split;
                                // the rows share the face cell at m = 0, and
                                // CROSS-FACE rows collide at m > 0 too (a face
                                // at (i,j) shifting +x and one at (i+2,j-2)
                                // shifting +y both write (i+2,j)). Both shares
                                // were booked into E_shift_local but the second
                                // assignment DISCARDED the first, so the field
                                // lost only one of them: measured 0.669 realised
                                // / booked on the 45 deg fixture and 0.848 at
                                // Meier scale (a 15.2 % leak, 7,616 J).
                            }
                        }
                    }
                    // ---- compose the deltas, THEN write ---------------------
                    for (const auto& kv : touched)
                    {
                        const ShiftCell& c = kv.second;
                        // ⚠⚠ ASSIGN WHEN TOUCHED ONCE. This is not an
                        // optimisation -- it is what keeps every frozen
                        // flat-wall anchor BIT-IDENTICAL. A blanket
                        //     Hw += (rNew - rH)
                        // is algebraically identity-preserving on a one-row wall
                        // but `old + (new - old) != new` in floating point, so
                        // it would drift all 70 anchor files at round-off and
                        // the identity gate would have to be loosened -- which
                        // is how a real change hides. (c1) is the witness.
                        Hw(c.ii, c.jj, c.kf) =
                            (c.cnt == 1) ? c.last_new
                                         : (c.H_old + c.sum_delta);
                        // the realised energy, read back FROM THE FIELD
                        E_real_local +=
                            (c.H_old - Hw(c.ii, c.jj, c.kf)) * Vcell;
                        if (c.cnt > 1) ++n_multi;
                    }
                }
            }

            // ---- deliverable 10: the orphaned-bank policy, with its counter --
            // ⚠⚠ POLICY: CARRY, DO NOT RESET. Across steps the vertical
            // criterion can remove a lip cell this column has partially banked,
            // leaving a remainder that refers to no cell. Resetting it would
            // DESTROY earned depth and (b1)'s conservation -- the one exact
            // claim this line owns -- would fail by exactly that amount. So the
            // remainder is carried and the occurrence is counted. It is ZERO on
            // this fixture (no vertical firings at all) and REAL on Leg A.
            long n_orphan = 0;
            for (int idx = 0; idx < ncols; ++idx)
                if (lat_accum_col[idx] > 0.0 && regime_col[idx] == 1) ++n_orphan;

            // ⚠⚠ E_shift_local IS DELIBERATELY NOT REDUCED, and getting this
            // wrong cost a (b2) failure at 0.0898 against a 5e-3 bar -- exactly
            // 4x the truth on 4 ranks. The row-build loop above runs over
            // COLUMNS using the already-reduced rows, so every rank computes the
            // same total; a sum reduction multiplies it by the rank count and
            // the ledger then books four times the enthalpy the field lost.
            // n_overlap is different -- it is accumulated inside the per-BOX
            // write-back, so it genuinely needs the reduction.
            // ⚠ Same trap as P1g-1's lat_bank / lat_bank_hf. The rule: reduce
            // what a per-box loop accumulated, never what a per-column loop did.
            amrex::ParallelDescriptor::ReduceLongSum(n_overlap);
            // ⚠⚠ P1g-2c: E_real_local and n_multi ARE accumulated per BOX, so
            // unlike E_shift_local they DO need the reduction -- the same rule,
            // applied the other way. Getting this backwards is what produced
            // exactly 4x on 4 ranks twice before (P1g-1's lat_bank, P1g-2's
            // E_shift_local). The rule: reduce what a per-BOX loop accumulated,
            // never what a per-COLUMN loop did.
            amrex::ParallelDescriptor::ReduceRealSum(E_real_local);
            amrex::ParallelDescriptor::ReduceLongSum(n_multi);
            lat_shift_n      = static_cast<Set::Scalar>(n_shift);
            lat_shift_refuse = static_cast<Set::Scalar>(n_refuse);
            lat_shift_E      = E_shift_local;
            lat_shift_E_est  = E_est_local;
            lat_shift_ratio  = (vol_bank > 0.0) ? (vol_disp / vol_bank) : 0.0;
            lat_overlap_n    = static_cast<Set::Scalar>(n_overlap);
            lat_orphan_n     = static_cast<Set::Scalar>(n_orphan);
            // P1g-2c deliverables 1, 2 and 6.
            lat_shift_real   = E_real_local;
            lat_shift_realfrac = (std::abs(E_shift_local) > 0.0)
                               ? (E_real_local / E_shift_local) : 0.0;
            lat_touch_multi  = static_cast<Set::Scalar>(n_multi);
            lat_twoaxis_n    = static_cast<Set::Scalar>(n_twoaxis);
            lat_ref_ell      = static_cast<Set::Scalar>(n_ref_ell);
            lat_ref_void     = static_cast<Set::Scalar>(n_ref_void);
            lat_ref_short    = static_cast<Set::Scalar>(n_ref_short);
            // the shift took this enthalpy out of the solid set
            ledger_E_flake  += E_shift_local;
        }

        if (lat_on)
        {
            amrex::ParallelDescriptor::ReduceRealSum(E_lat_local);
            amrex::ParallelDescriptor::ReduceLongSum(lat_voided);
            lat_voided_n    = static_cast<Set::Scalar>(lat_voided);
            ledger_E_flake += E_lat_local;
            // P2-4: and the same quantity booked SEPARATELY, so the price can be
            // split by regime. ⚠ Not a second sink -- it is a share of the one
            // above, so it must NOT also be credited in ledger_err_solid or the
            // wall's energy would be counted twice.
            ledger_E_flake_wall += E_lat_local;
        }
        // P1g-1: the banking diagnostics. ⚠ Guarded on the KEY ALONE, not on
        // bank_on, and that is deliberate: the NEG leg sets the key with the
        // drive off, and (c3)'s bar is that the banking counters READ ZERO.
        // Guarding these on bank_on would remove the columns entirely and (c3)
        // would have nothing to read -- it would pass by absence instead of by
        // measurement, which is the vacuity trap.
        if (sp_lateral_subcell)
        {
            // ⚠ deliberately NOT reduced -- see the pre-loop's note: the
            // banking loop runs identically on every rank, so a sum reduction
            // would report rank_count x the true figure.
            lat_bank_n  = static_cast<Set::Scalar>(lat_bank_local);
            lat_bank_hf = lat_bank_hf_local;
            // the per-column max, reduced as a MAX so one runaway column
            // survives the reduction instead of being averaged away
            int colmax = 0;
            for (int c = 0; c < static_cast<int>(lat_void_col_local.size()); ++c)
                if (lat_void_col_local[c] > colmax) colmax = lat_void_col_local[c];
            amrex::ParallelDescriptor::ReduceIntMax(colmax);
            lat_void_colmax = static_cast<Set::Scalar>(colmax);
            // Sum of the accumulator over every column. ⚠ P2 predicts this is
            // NON-ZERO at run end: with columns stopping part-way through a flake,
            // a zero would mean the accumulator is being flushed somewhere.
            // ⚠ NOT reduced across ranks -- lat_accum_col is per-column host state
            // already made rank-consistent below, so a sum reduction would
            // multiply it by the rank count.
            Set::Scalar acc_sum = 0.0;
            for (int c = 0; c < static_cast<int>(lat_accum_col.size()); ++c)
                acc_sum += lat_accum_col[c];
            lat_accum_sum = acc_sum;
            // ⚠ no reduction: both log arrays are filled by the rank-invariant
            // pre-loop, so every rank already holds the same values.
        }

        // A1: every-removal event log, written by the I/O rank after the void
        // flip. Rows ordered by column index (col_j, then col_i); the file
        // is identical for any rank count / box layout. Reals are printed
        // with 17 significant digits (exact round trip).
        //   T_top, a_f, Sp_top : at k_top (k_top_global, `removed` mask +
        //                        phi >= 0) before removal
        //   h_scan    : depth-scan h_col (0 if the scan did not fire)
        //   h_applied : h_chosen_col, the phi shift actually applied
        //   n_voided  : cells flipped to removed in the column this step
        if (rev_on && sp_emit_removal_events_csv)
        {
            if (ncols > 0)
                amrex::ParallelDescriptor::ReduceIntSum(rev_nvoid.data(), ncols);
            if (amrex::ParallelDescriptor::IOProcessor())
            {
                const std::string csv_path = plot_file + "_removal_events.csv";
                const bool need_header = !sp_removal_events_csv_header_written;
                std::ofstream ofs(csv_path, (need_header ? std::ios::trunc : std::ios::app)
                                            | std::ios::out);
                if (ofs.is_open())
                {
                    // P0: the reconstruction columns are APPENDED and
                    // exist only under spall.surface_T = reconstructed. A new
                    // field breaks the 2639-file identity on the header alone,
                    // whatever the numbers do, so the default header is the
                    // pre-P0 header character for character.
                    //   T_surf     reconstructed surface temperature [K]
                    //   ell_eff    recovered decay length [m] (0 = flat)
                    //   h_f        flake thickness min(ell_eff, dz) [m]
                    //   T_bar_flake flake-average temperature [K]
                    //   rc_flag    1 flat | 2 ell_eff clamped | 4 T_surf capped
                    //              | 8 h_f at spall.h_spall_max (the SEARCH
                    //                depth) | 16 h_f at spall.flake_h_max (the
                    //                physical flake bound, rev7 deliverable 1)
                    //   d_ell_rel  (ell_eff - previous)/previous for this column
                    if (need_header)
                        ofs << "time,col_i,col_j,regime,k_top,T_top,a_f,Sp_top,"
                               "h_scan,h_applied,n_voided"
                            << (recon_on ? ",T_surf,ell_eff,h_f,T_bar_flake,"
                                           "rc_flag,d_ell_rel,ell_raw,delta,"
                                           "q_ratio" : "")
                            // ⚠⚠ P1g-3 deliverable 7: APPENDED, and ONLY under
                            // the two new keys, so a keys-off run's header is
                            // character for character what it was and (b1)'s
                            // 0-of-2639 cannot break on the header. This is the
                            // third step in a row where an unguarded new column
                            // would have failed the identity gate for a
                            // reporting reason rather than a physics one.
                            //   ell_q  the FLUX-based heated depth actually used
                            //          (0 when the shipped two-cell length was)
                            << ((sp_flake_len_flux || sp_flake_frac_weibull)
                                ? ",ell_q" : "")
                            // P1b: APPENDED, and only under spall.support_rule,
                            // so the default header stays character for
                            // character what it was and the 2639-file identity
                            // cannot break on the header alone.
                            //   n_exposed  lateral neighbours >= support_step_min
                            //              below this column (0-4)
                            // P1f: riser_depth is how far BELOW the column's top the
                            // lateral firing occurred -- the quantity (h1) needs to
                            // ask whether firing propagates UPWARD from the corner,
                            // which is where deliverable 1 says ignition is.
                            // ⚠⚠ deliverable 10: K_I_lat and K_Ic_lat are
                            // Pa*sqrt(m), q_lat = K_I/K_Ic is dimensionless and
                            // fires at >= 1, sigma_lat is Pa AT THE FACE
                            // (xi = 0). These four are what the branch decides
                            // on; none of them was an output before, which is
                            // why the Pa-vs-MPa comparison error was invisible.
                            // ⚠ Still inside `lat_csv`, so the default header
                            // stays character for character what it was and the
                            // 2639-file identity cannot break on the header.
                            << (lat_csv ? ",riser_depth,ell_lat,T_face_lat,"
                                          "K_I_lat,K_Ic_lat,q_lat,sigma_lat,"
                                          "h_f_lat,h_f_unc_lat,lat_accum,"
                                          "lat_nfire" : "")
                            << "\n";
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    for (int idx = 0; idx < ncols; ++idx)
                    {
                        if (mech_on && mech_n[idx] > 0)
                            ofs << time << "," << (idx % nx_dom + dom_ilo) << ","
                                << (idx / nx_dom + dom_jlo) << ",4,"
                                << mech_ktop[idx] << ",0,0,0,0,"
                                << Set::Scalar(mech_n[idx]) * dz << ","
                                << mech_flip[idx]
                                << (recon_on ? ",0,0,0,0,0,0,0,0" : "")
                                << (lat_csv ? ",0,0,0,0,0,0,0,0" : "") << "\n";
                        if (regime_col[idx] == 0) continue;
                        ofs << time << "," << (idx % nx_dom + dom_ilo) << ","
                            << (idx / nx_dom + dom_jlo) << "," << regime_col[idx] << ","
                            << k_top_global[idx] << "," << rev_T_top[idx] << ","
                            << rev_a_f[idx] << "," << rev_Sp[idx] << ","
                            << rev_h_scan[idx] << "," << h_chosen_col[idx] << ","
                            << rev_nvoid[idx];
                        if (recon_on)
                            ofs << "," << rc_Tsurf[idx] << "," << rc_ell[idx]
                                << "," << rc_h_f[idx]   << "," << rc_Tbar[idx]
                                << "," << rc_flag[idx]  << "," << rc_dell[idx]
                                // P0b deliverable 6: the RAW (unsmoothed)
                                // ell_eff beside the smoothed one -- report
                                // both, so smoothing cannot hide a transient --
                                // and `delta`, the geometric offset of the
                                // thermal frame from the cell's top face after
                                // this firing.
                                << "," << (recon_on && idx < static_cast<int>(surf_ell_raw.size())
                                           ? surf_ell_raw[idx] : 0.0)
                                << "," << (remap_on && !phi_top_col.empty()
                                           ? phi_top_col[idx] - h_chosen_col[idx]
                                             + Set::Scalar(nvoid_pred[idx]) * dz
                                           : 0.0)
                                // P1 deliverable 5.
                                << "," << (idx < static_cast<int>(rc_qrat.size())
                                           ? rc_qrat[idx] : 0.0);
                        // ⚠⚠ P1g-3: WRITTEN HERE, AFTER the recon block, because
                        // the header appends ",ell_q" AFTER the recon columns.
                        // An earlier version wrote it BEFORE them, which would
                        // have silently misaligned every column after it -- the
                        // header and the row must be appended in the same order,
                        // and nothing in a column-count check would catch it.
                        if (sp_flake_len_flux || sp_flake_frac_weibull)
                            ofs << "," << (rc_ellq.empty() ? 0.0 : rc_ellq[idx]);
                        if (lat_csv)
                            ofs << "," << (lat_depth.empty() ? 0.0 : lat_depth[idx])
                                << "," << (lat_ell.empty()   ? 0.0 : lat_ell[idx])
                                << "," << (lat_Tf.empty()    ? 0.0 : lat_Tf[idx])
                                << "," << (lat_KI.empty()    ? 0.0 : lat_KI[idx])
                                << "," << (lat_KIc.empty()   ? 0.0 : lat_KIc[idx])
                                << "," << (lat_q.empty()     ? 0.0 : lat_q[idx])
                                << "," << (lat_sig.empty()   ? 0.0 : lat_sig[idx])
                                << "," << (lat_hf.empty()    ? 0.0 : lat_hf[idx])
                                << "," << (lat_hfu.empty()     ? 0.0 : lat_hfu[idx])
                                << "," << (lat_acc_log.empty() ? 0.0 : lat_acc_log[idx])
                                << "," << (lat_nf_log.empty()  ? 0.0 : lat_nf_log[idx]);
                        ofs << "\n";
                    }
                }
            }
            sp_removal_events_csv_header_written = true;
        }

        // AMR: the band must still hold the surface and the scan window.
        if (max_level > 0) CheckSurfaceCovered(lev, time);

        // P0b rev5 deliverable 6 -- the second half of the (e) fix. The remap
        // wrote H and nothing else, so T (and the Lambda fractions) are one
        // step stale on exactly the cells it touched. Recover them the way the
        // rest of the integrator does -- the SAME InvertHtoT the thermal
        // advance calls -- and then refill the ghost layer, because temp_mf
        // carries ghosts, the conduction kernel reads them laterally, and a
        // stale ghost is visible only at a box seam. Both calls are collective
        // and remap_any is built from reduced per-column buffers, so every
        // rank takes the same branch.
        // ⚠ InvertHtoT over the whole level is bit-idempotent on every cell
        // the remap did not touch: T there already IS Invert(H) from the
        // previous step's call, and Invert is a deterministic cell-local
        // function of H alone. Void cells included -- their H is frozen, so
        // re-inverting returns the frozen T.
        // ⚠⚠ P1g-2 deliverable 7 rides on EXACTLY this hook, and that is the
        // point: the lateral shift writes H and NOTHING else, so T and the
        // Lambda fractions come back through the same InvertHtoT the vertical
        // remap uses, and the ghost layer is refilled by the same call. The P0b
        // lesson in full -- writing T directly leaves (H, T) satisfying no
        // equation of state AND leaves every ghost copy stale until the next
        // fill, visible only at a box seam. The lateral shift writes rows that
        // run straight into box boundaries, so it is MORE exposed to that bug
        // than the vertical remap was, not less.
        // ⚠ Both calls are collective, and lat_shift_any is built from reduced
        // per-column buffers, so every rank takes the same branch.
        if (remap_any || lat_shift_any)
        {
            InvertHtoT(lev);
            Util::RealFillBoundary(*temp_mf[lev], geom[lev]);
        }

        // Rebuild diagnostics that depend on the new mask/D field.
        UpdateSurfaceMaskFromPhi(lev);
        RefreshKappaEff(lev);
        (void)dz;
    }

#endif
