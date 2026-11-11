    void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::UpdateModel(a_step, a_time);

        // Step 8b: prescribed-temperature surface_patch acts as a flux source
        // on the cell-centred heat conduction (AdvanceMicrostructure), but the
        // elastic eigenstrain below reads cell-centred T through
        // CellToNodeAverage. With Neumann zhi on temp_mf, the ghost row equals
        // the topmost real cell, so the surface node sees a half-cell-lagged
        // temperature instead of the prescribed Dirichlet sp_T(t). Override
        // the top-node temperature inside the patch with sp_T(a_time) so that
        // the eigenstrain at z=z_top reflects the BC Hu prescribes.
        // Step 15a: the top-node eigenstrain override below is meaningful
        // only under mode = prescribed_T. Under convective_flame the surface
        // temperature is whatever the heat solve produces, so the top node
        // reads solved temp_mf through CellToNodeAverage with no override.
        const bool        sp_on  = surface_patch.enabled &&
            surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        const Set::Scalar sp_x0c = surface_patch.x0;
        const Set::Scalar sp_y0c = surface_patch.y0;
        const Set::Scalar sp_r2c = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_Tnow = sp_on ? surface_patch.T_at(a_time) : 0.0;

        if (microstructure_enabled)
        {
            // Step 7/alpha-beta: heterogeneous mechanics. Build each nodal
            // model from cell-centred phase stiffness and cell-centred free
            // strain, both averaged onto the node. The free strain is
            // beta*(T-T_ref) plus an optional phase-local alpha-beta
            // transition strain, so nonlinear mineral transitions are
            // evaluated before nodal averaging.
            // Stiffness uses the (E, mu) → lambda relation from
            // Affine::Isotropic::Parse:
            //     lambda = mu * (E - 2*mu) / (3*mu - E).
            // F0 = eps_free_node * I.
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                Util::RealFillBoundary(*beta_mf[lev],        geom[lev]);
                Util::RealFillBoundary(*E_phase_mf[lev],     geom[lev]);
                Util::RealFillBoundary(*mu_phase_mf[lev],    geom[lev]);
                Util::RealFillBoundary(*T_ref_phase_mf[lev], geom[lev]);
                Util::RealFillBoundary(*alpha_beta_linear_strain_mf[lev], geom[lev]);
                Util::RealFillBoundary(*alpha_beta_T0_mf[lev],            geom[lev]);
                Util::RealFillBoundary(*alpha_beta_width_mf[lev],         geom[lev]);
                UpdateMicrostructureFreeStrain(lev);

                const Set::Scalar* DXl  = geom[lev].CellSize();
                const Set::Scalar* PLOl = geom[lev].ProbLo();
                const amrex::Box domain = geom[lev].Domain();
                const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
                const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
                const int dom_klo = domain.smallEnd(2), dom_khi = domain.bigEnd(2);
                const int sp_topk_cell = domain.bigEnd(2);
                const int sp_topk_node = geom[lev].Domain().bigEnd(2) + 1;
                for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.grownnodaltilebox();
                    amrex::Array4<Model::Solid::Affine::Isotropic> const& model =
                        model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& epsc =
                        thermal_free_strain_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& Ec   = E_phase_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& muc  = mu_phase_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& betc = beta_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& Trc  = T_ref_phase_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& abec = alpha_beta_linear_strain_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& abTc = alpha_beta_T0_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& abWc = alpha_beta_width_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Set::Scalar epsn = CellToNodeAverageInDomain(
                            epsc, i, j, k, 0,
                            dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                            0.1);
                        if (sp_on && k == sp_topk_node) {
                            Set::Scalar xn = PLOl[0] + Set::Scalar(i) * DXl[0];
                            Set::Scalar yn = PLOl[1] + Set::Scalar(j) * DXl[1];
                            Set::Scalar dxp = xn - sp_x0c, dyp = yn - sp_y0c;
                            if (dxp*dxp + dyp*dyp <= sp_r2c) {
                                Set::Scalar sum = 0.0;
                                int count = 0;
                                for (int di = 0; di <= 1; ++di) {
                                    int ci = i - di;
                                    if (ci < dom_ilo || ci > dom_ihi) continue;
                                    for (int dj = 0; dj <= 1; ++dj) {
                                        int cj = j - dj;
                                        if (cj < dom_jlo || cj > dom_jhi) continue;
                                        sum += ThermalFreeStrain(
                                            sp_Tnow,
                                            betc(ci, cj, sp_topk_cell, 0),
                                            Trc (ci, cj, sp_topk_cell, 0),
                                            abec(ci, cj, sp_topk_cell, 0),
                                            abTc(ci, cj, sp_topk_cell, 0),
                                            abWc(ci, cj, sp_topk_cell, 0));
                                        ++count;
                                    }
                                }
                                if (count > 0) epsn = sum / Set::Scalar(count);
                            }
                        }
                        Set::Scalar En = CellToNodeAverageInDomain(
                            Ec, i, j, k, 0,
                            dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                            1.0e13);
                        Set::Scalar mun = CellToNodeAverageInDomain(
                            muc, i, j, k, 0,
                            dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                            1.0e13);
                        Set::Scalar denom = 3.0 * mun - En;
                        // Guard against the degenerate (lambda denominator
                        // small / non-finite moduli) input combinations.
                        Set::Scalar lam = (std::abs(denom) > 1.0e-30)
                                        ? mun * (En - 2.0 * mun) / denom
                                        : 0.0;
                        if (!std::isfinite(lam) || !std::isfinite(mun)) {
                            lam = 0.0; mun = 0.0;
                        }
                        Set::Matrix F0 = epsn * Set::Matrix::Identity();
                        model(i, j, k).Define(mun, lam, F0);
                    });
                }
                Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            }
            return;
        }
    }

    // P0b: declared here, ahead of its first use as a parameter type in
    // the column reconstruction (deleted, P2-28). An in-class member function BODY is a
    // complete-class context, but its DECLARATION is not.
    struct SurfaceRecon
    {
        Set::Scalar T_surf  = 0.0;   // reconstructed surface temperature [K]
        Set::Scalar ell_eff = 0.0;   // recovered decay length [m]
        // the face list's old tangent depth k_c (T_s - T0) / q_net, logged
        // only (FaceRemoval.H); 0 everywhere else
        Set::Scalar ell_q   = 0.0;
        Set::Scalar h_f     = 0.0;   // flake thickness min(ell_eff, dz) [m]
        Set::Scalar T_bar   = 0.0;   // flake-average temperature [K]
        bool ell_clamped    = false; // the ell_eff >= dz/4 floor bound
        bool h_f_capped     = false; // h_f hit the spall.h_spall_max SEARCH depth
        bool h_f_bounded    = false; // h_f hit spall.flake_h_max (deliverable 1)
        bool T_clamped      = false; // the gas-temperature ceiling bound
        bool flat           = false; // no usable gradient; T_surf == T_top
    };

    // P2-21: the column optical depth (ColumnOpticalDepth), built once per step
    // per path when the loss table varies and the beam is on. Indexed
    // [col*nz + (k - klo)] as ColumnSolidAboveCount; tau_top is the optical depth
    // from the column's top solid face to the cell's top face, tau_bot to its
    // bottom face (tau_top + a_k dz, the SAME a_k the gather used). k_last is the
    // column's lowest solid cell (the transmitted power's exit), INT_MIN if none.
    // Declared here, ahead of its first use as a parameter type (as SurfaceRecon).
    struct BeamTau
    {
        amrex::Vector<Set::Scalar> top, bot;
        amrex::Vector<int>         k_last;
        // P2-24: the column's reflection R_col = R + (R_melt - R) f_L(T_top), one per column (ncols); used at
        // P_inc, P_trans and every cell's (1 - R) prefactor in that column
        amrex::Vector<Set::Scalar> Rcol;
        bool built = false;
    };
    // P2-21: a cell's (tau_top, tau_bot) for the kernels, captured by value in a
    // ParallelFor. Unbuilt => -1, -1 => BeamSource's closed form.
    struct TauView
    {
        const Set::Scalar* top = nullptr;
        const Set::Scalar* bot = nullptr;
        int ilo = 0, jlo = 0, klo = 0, nx = 0;
        std::size_t nz = 0;
        AMREX_GPU_HOST_DEVICE std::size_t at(int i, int j, int k) const
        {
            return (static_cast<std::size_t>(j - jlo) * nx + (i - ilo)) * nz + (k - klo);
        }
        AMREX_GPU_HOST_DEVICE Set::Scalar t(int i, int j, int k) const { return top ? top[at(i, j, k)] : -1.0; }
        AMREX_GPU_HOST_DEVICE Set::Scalar b(int i, int j, int k) const { return bot ? bot[at(i, j, k)] : -1.0; }
        // P2-24: the column's R when built, else the deck's scalar R (byte for byte)
        const Set::Scalar* rc = nullptr;
        AMREX_GPU_HOST_DEVICE Set::Scalar R(int i, int j, Set::Scalar R_scalar) const
        {
            return rc ? rc[static_cast<std::size_t>(j - jlo) * nx + (i - ilo)] : R_scalar;
        }
    };
    TauView ViewOf(const BeamTau& bt, int lev) const
    {
        TauView v;
        if (!bt.built) return v;
        const amrex::Box& d = geom[lev].Domain();
        v.top = bt.top.data();
        v.bot = bt.bot.data();
        v.rc  = bt.Rcol.empty() ? nullptr : bt.Rcol.data();
        v.ilo = d.smallEnd(0); v.jlo = d.smallEnd(1); v.klo = d.smallEnd(2);
        v.nx  = d.length(0);   v.nz  = static_cast<std::size_t>(d.length(2));
        return v;
    }
