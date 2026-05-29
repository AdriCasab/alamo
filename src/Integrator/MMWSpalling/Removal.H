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
        const bool have_stress = spall_on && !spall_prescribed;
        if (have_stress) Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

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
        const Set::Scalar T_amb_local = T_amb;
        const Set::Scalar dt_local = dt;

        auto col_index = [=](int i, int j) {
            return (j - dom_jlo) * nx_dom + (i - dom_ilo);
        };

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
                    if (rem(i, j, k) <= 0.5 && phi(i, j, k) >= 0.0) { k_top_local = k; break; }
                }
                if (k_top_local > std::numeric_limits<int>::min())
                {
                    int idx = col_index(i, j);
                    if (k_top_local > k_top_global[idx]) k_top_global[idx] = k_top_local;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(k_top_global.data(), ncols);

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

        // Step 16: capture sp_weibull-specific K_Ic + scan parameters
        // outside the MFIter so the lambda in the depth scan can close over
        // them. These values are unused when spallation_model == DamageLaw.
        const KIcLaw      law_loc          = kic_law;
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
            if (!(phase_nu_loc > -1.0 && phase_nu_loc < 0.5))
                Util::Abort(INFO, "Step 16d: phase 0 Poisson's ratio out "
                                  "of physical range; got nu=", phase_nu_loc);
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
        const bool sp_weibull_use_T_face =
            sp_weibull_spall &&
            surface_patch.enabled &&
            surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        if (sp_weibull_spall && !sp_weibull_use_T_face && !have_stress)
            Util::Abort(INFO, "Step 19 (sp_weibull + spall.enabled under "
                              "non-prescribed_T heat source): mechanics "
                              "solve must be active to source sigma_xx "
                              "from stress_mf — set el.type = static and "
                              "el.time_evolving = 1, or switch to "
                              "surface_patch.mode = prescribed_T for the "
                              "Step 16d T-derived path");

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

        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> sig_node;
            if (have_stress) sig_node = (*stress_mf[lev]).array(mfi);
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
                clu       = (*Sp_cluster_id_mf[lev]).array(mfi);
                fa_a      = (*flaw_a_mf[lev]).array(mfi);
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
                        Set::Scalar cluster_id = clu(i, j, k_top);
                        if (cluster_id > 0.0)
                        {
                            Set::Scalar a_f = fa_a(i, j, k_top);
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
                                const int max_steps_cap = static_cast<int>(
                                    std::floor(h_spall_max_loc / dz));
                                const int max_steps_box = k_top - klo;
                                const int max_steps =
                                    std::min(max_steps_cap, max_steps_box);
                                const int klo_loc = klo;
                                const int khi_loc = khi;

                                // Local helper: K_Ic(T) under the
                                // configured law, in Pa·m^0.5.
                                auto kic_pa_at_cell =
                                    [&] (int k_cell) -> Set::Scalar {
                                    Set::Scalar T_at = T_a(i, j, k_cell);
                                    Set::Scalar K_Ic_MPa;
                                    if (law_loc == KIcLaw::NasseriTable)
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
                                auto k_i_at_zcrack =
                                    [&] (Set::Scalar z_crack) -> Set::Scalar {
                                    auto sigma_at_depth =
                                        [&] (Set::Scalar d) -> Set::Scalar {
                                        const Set::Scalar z_abs =
                                            z_crack + d;
                                        if (!sp_weibull_use_T_face)
                                        {
                                            // Step 19: stress_mf path.
                                            int k_d = k_top -
                                                static_cast<int>(
                                                    std::floor(z_abs / dz));
                                            if (k_d > k_top)   k_d = k_top;
                                            if (k_d < klo_loc) k_d = klo_loc;
                                            if (k_d > khi_loc) k_d = khi_loc;
                                            Set::Matrix sigma =
                                                Numeric::Interpolate
                                                ::NodeToCellAverage(
                                                    sig_node, i, j, k_d, 0);
                                            return -sigma(0, 0)
                                                   + sp_conf_offset;
                                        }
                                        const Set::Scalar half_dz =
                                            0.5 * dz;
                                        Set::Scalar T_z;
                                        if (z_abs <= half_dz)
                                        {
                                            // Surface region: lerp
                                            // between T_face (prescribed)
                                            // and top-cell-centre T.
                                            const Set::Scalar T_face =
                                                sp_T_face_loc(time);
                                            const Set::Scalar T_top =
                                                T_a(i, j, k_top);
                                            Set::Scalar alpha =
                                                z_abs / half_dz;
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
                                                z_abs / dz - 0.5;
                                            const Set::Scalar uf =
                                                std::floor(u);
                                            int n = static_cast<int>(uf);
                                            int k_above = k_top - n;
                                            int k_below = k_above - 1;
                                            if (k_above > k_top)   k_above = k_top;
                                            if (k_above < klo_loc) k_above = klo_loc;
                                            if (k_above > khi_loc) k_above = khi_loc;
                                            if (k_below > k_top)   k_below = k_top;
                                            if (k_below < klo_loc) k_below = klo_loc;
                                            if (k_below > khi_loc) k_below = khi_loc;
                                            const Set::Scalar alpha =
                                                u - uf;  // in [0, 1)
                                            T_z = (1.0 - alpha)
                                                  * T_a(i, j, k_above)
                                                +        alpha
                                                  * T_a(i, j, k_below);
                                        }
                                        // 1-D confinement: physical sigma_xx
                                        // = -E*beta*dT/(1-nu). K_I integrator
                                        // wants compression-positive, so the
                                        // double-negation returns +E*beta*dT/
                                        // (1-nu). + sp_conf_offset adds the
                                        // confining-pressure contribution.
                                        return phase_E_loc * phase_beta_loc
                                               * (T_z - phase_T_ref_loc)
                                               / (1.0 - phase_nu_loc)
                                               + sp_conf_offset;
                                    };
                                    return Numeric::SpCriterion::K_I(
                                        sigma_at_depth, a_f,
                                        h_spall_n_seg_loc);
                                };

                                int  s_last_pass = -1;
                                bool found_failure = false;
                                for (int s = 0; s <= max_steps; ++s)
                                {
                                    const Set::Scalar z_crack =
                                        Set::Scalar(s) * dz;
                                    const int k_tip = k_top - s;
                                    if (k_tip < klo) break;

                                    Set::Scalar K_Ic_Pa =
                                        kic_pa_at_cell(k_tip);
                                    Set::Scalar K_I_val =
                                        k_i_at_zcrack(z_crack);

                                    if (K_I_val >= K_Ic_Pa)
                                    {
                                        s_last_pass = s;
                                        h_col       = z_crack;
                                        // Step 16b: edge-triggered event
                                        // count at (i, j, k_tip).
                                        if (cec(i, j, k_tip) < 0.5)
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
                                    (k_top - (s_last_pass + 1)) >= klo)
                                {
                                    Set::Scalar z_lo = Set::Scalar(s_last_pass)
                                                       * dz;       // pass
                                    Set::Scalar z_hi = Set::Scalar(s_last_pass + 1)
                                                       * dz;       // fail
                                    for (int b = 0; b < 8; ++b)
                                    {
                                        Set::Scalar z_mid =
                                            0.5 * (z_lo + z_hi);
                                        int k_mid = k_top -
                                            static_cast<int>(
                                                std::floor(z_mid / dz));
                                        if (k_mid < klo) k_mid = klo;
                                        if (k_mid > khi) k_mid = khi;
                                        Set::Scalar K_Ic_mid =
                                            kic_pa_at_cell(k_mid);
                                        Set::Scalar K_I_mid =
                                            k_i_at_zcrack(z_mid);
                                        if (K_I_mid >= K_Ic_mid)
                                            z_lo = z_mid;
                                        else
                                            z_hi = z_mid;
                                    }
                                    h_col = z_lo;
                                }

                                if (h_col > 0.0)
                                {
                                    spall_fires_col[idx] = 1;
                                    h_s_cand_col[idx]    = h_col;
                                    rop_s_cand_col[idx]  =
                                        (dt_local > 0.0) ? h_col / dt_local : 0.0;
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
                    Set::Scalar Cp_p  = phases[p_top].Cp;
                    Set::Scalar rho_p = phases[p_top].rho;
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
        if (ncols > 0)
        {
            amrex::ParallelDescriptor::ReduceIntSum(spall_fires_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceIntSum(vapor_fires_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(h_s_cand_col.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(h_v_cand_col.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rop_s_cand_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rop_v_cand_col.data(), ncols);
        }

        // Step 16c: rank-ordered append of per-column first-firing events
        // to <plot_file>_h_col_events.csv. Each rank writes its own
        // subset; barriers serialize the file writes so the header (from
        // the I/O proc on the first emission) and subsequent rows don't
        // interleave. Append mode preserves cross-step accumulation; the
        // I/O proc latch ensures the header is emitted exactly once per
        // run (re-runs of the same plot_file truncate-and-rewrite).
        if (sp_weibull_spall && lev == 0)
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
            regime_col[idx]     = reg_choice;
            h_chosen_col[idx]   = h_c;
            rop_chosen_col[idx] = (dt_local > 0.0) ? h_c / dt_local : 0.0;
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
                        rem(i, j, k)  = 1.0;
                        reg(i, j, k)  = static_cast<Set::Scalar>(reg_choice);
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
                        T(i, j, k) = T_amb_local;
                        if (ptabs)
                        {
                            int p = static_cast<int>(std::round(ph(i, j, k)));
                            if (p < 0 || p >= n_phases) p = 0;
                            Set::Scalar Hv = ptabs[p].H_of_T(T_amb_local);
                            Set::Scalar Tdummy, ls, ll, lv;
                            ptabs[p].Invert(Hv, Tdummy, ls, ll, lv);
                            H(i, j, k)  = Hv;
                            LS(i, j, k) = ls;
                            LL(i, j, k) = ll;
                            LV(i, j, k) = lv;
                        }
                        else
                        {
                            H(i, j, k)  = 0.0;
                            LS(i, j, k) = 1.0;
                            LL(i, j, k) = 0.0;
                            LV(i, j, k) = 0.0;
                        }
                    }
                }
            }
        }

        // Rebuild diagnostics that depend on the new mask/D field.
        UpdateSurfaceMaskFromPhi(lev);
        RefreshKappaEff(lev);
        (void)dz;
    }

#endif
