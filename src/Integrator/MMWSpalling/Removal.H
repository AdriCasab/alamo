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

        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);
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
                                const int max_steps =
                                    std::min(max_steps_cap, k_top - scan_kbot);
                                // Column temperature by depth-below-surface index
                                // (cells): gathered profile for the T-path (box-
                                // independent), box-local T_a otherwise.
                                const Set::Scalar* Tcol_i =
                                    (T_col_p && sp_weibull_use_T_face)
                                      ? (T_col_p + static_cast<std::size_t>(idx) * colw)
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

                                // Local helper: K_Ic(T) under the
                                // configured law, in Pa·m^0.5.
                                auto kic_pa_at_cell =
                                    [&] (int k_cell) -> Set::Scalar {
                                    Set::Scalar T_at = T_at_didx(k_top - k_cell);
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
                                        // z_abs is a NORMAL distance under
                                        // Step 20; z_v is the vertical depth
                                        // at which the column profile holds
                                        // the same material point.
                                        const Set::Scalar z_abs =
                                            z_crack + d;
                                        const Set::Scalar z_v =
                                            z_abs * sn_icth;
                                        if (!sp_weibull_use_T_face)
                                        {
                                            // Step 19: stress_mf path.
                                            int k_d = k_top -
                                                static_cast<int>(
                                                    std::floor(z_v / dz));
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
                                        if (z_v <= half_dz)
                                        {
                                            // Surface region: lerp between the
                                            // surface T and the top-cell-centre
                                            // T. Surface T is the prescribed-T
                                            // parser when in prescribed_T mode,
                                            // else (local_thermoelastic under a
                                            // computed-T source) the top cell.
                                            const Set::Scalar T_top =
                                                T_at_didx(0);
                                            const Set::Scalar T_face =
                                                sp_weibull_presc_T_face
                                                ? sp_T_face_loc(time)
                                                : T_top;
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
                                    // Cell k_top - s sits at normal depth
                                    // s*dz*cos(theta) below the inclined
                                    // surface (s*dz vertical).
                                    const Set::Scalar z_crack =
                                        Set::Scalar(s) * dz * sn_cth;
                                    const int k_tip = k_top - s;
                                    if (k_tip < scan_kbot) break;

                                    Set::Scalar K_Ic_Pa =
                                        kic_pa_at_cell(k_tip);
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
        if (sp_weibull_spall && lev == 0 && sp_emit_h_col_events_csv)
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
            }
        }

        // Rebuild diagnostics that depend on the new mask/D field.
        UpdateSurfaceMaskFromPhi(lev);
        RefreshKappaEff(lev);
        (void)dz;
    }

#endif
