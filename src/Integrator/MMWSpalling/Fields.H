    // The H<->T tables are per phase, built in InitializeMicrostructure(); see
    // phase_tables. A cell reads the table of its phase ID.

    void ConvertTtoH(int lev)
    {
        const Numeric::Material::Table *ptabs = phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& H    = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Hold = (*H_old_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph = (*phase_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tv = T(i, j, k);
                Set::Scalar Hv = 0.0;
                Set::Scalar ls = 1.0, ll = 0.0, lv = 0.0;
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                Hv = ptabs[p].H_of_T(Tv);
                Set::Scalar dummyT;
                ptabs[p].Invert(Hv, dummyT, ls, ll, lv);
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
                H(i,j,k)    = Hv;
                Hold(i,j,k) = Hv;
            });
        }
    }

    void InvertHtoT(int lev)
    {
        // P2-29: one path -- the per-cell phase table
        const Numeric::Material::Table *ptabs = phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& H = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& T  = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph = (*phase_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Hv = H(i, j, k);
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                Set::Scalar Tv, ls, ll, lv;
                ptabs[p].Invert(Hv, Tv, ls, ll, lv);
                T(i,j,k) = Tv;
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
            });
        }
    }

    void InitSurfaceMask(int lev)
    {
        if (list_state())
        {
            UpdateSurfaceMaskFromPhi(lev);
            return;
        }
        const int top_k = geom[lev].Domain().bigEnd(2);
        for (amrex::MFIter mfi(*surface_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& mask = (*surface_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                mask(i, j, k) = (k == top_k) ? 1.0 : 0.0;
            });
        }
    }

    // The phi-window surface mask: 1 where 0 <= phi < dz, else 0 -- a pointwise
    // test on the level set, no column scan.
    //
    // D2q-1: this mask is deliberately left on phi and is NOT rebuilt from the
    // `removed` mask: "solid here, void directly above" is a different mask
    // (see studies/d2q1_percellvoid/RESULTS.md section "Site 7").
    // It is built only under list_state(), i.e. under the face list, where the
    // kernel and EnergyLedger take the surface cell from the mask column top
    // (surface_follow_mask) and do not read it for the flux: its readers there
    // are the plotfile (`surface`) and the ledger's surface_missing_cols count.
    void UpdateSurfaceMaskFromPhi(int lev)
    {
        const Set::Scalar dz = geom[lev].CellSize()[2];
        for (amrex::MFIter mfi(*surface_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& mask = (*surface_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar p = phi(i, j, k);
                mask(i, j, k) = (p >= 0.0 && p < dz) ? 1.0 : 0.0;
            });
        }
    }

    // D2q-2b item 0: flip the seeded index-space boxes to void, once, right
    // after InitPhi. Sets `removed` ONLY -- phi is deliberately left alone, so
    // the two genuinely disagree and every consumer is forced through the
    // authoritative mask that D2q-1 installed. That disagreement is the whole
    // point of the probe: with phi still saying "solid", any path that quietly
    // fell back to phi shows up immediately.
    //
    // No-op unless `spall.seed_void` was given.
    void ApplySeedVoid(int lev)
    {
        if (spall_seed_void.size() < 6 || spall_seed_void.size() % 6 != 0)
            return;
        long n = 0;
        for (std::size_t sb = 0; sb < spall_seed_void.size(); sb += 6)
        {
            const amrex::Box seed(
                amrex::IntVect{AMREX_D_DECL(spall_seed_void[sb + 0],
                                            spall_seed_void[sb + 1],
                                            spall_seed_void[sb + 2])},
                amrex::IntVect{AMREX_D_DECL(spall_seed_void[sb + 3],
                                            spall_seed_void[sb + 4],
                                            spall_seed_void[sb + 5])});
            for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & seed;
                if (!bx.ok()) continue;
                amrex::Array4<Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
                const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                    if (rem(i, j, k) <= 0.5) { rem(i, j, k) = 1.0; ++n; }
            }
        }
        amrex::ParallelDescriptor::ReduceLongSum(n);
        Util::Message(INFO, "spall.seed_void: level ", lev, ", voided ", n, " cells");
    }

    // Initialize phi to a flat top surface at z = prob_hi[2]. The cell-
    // centred value is z_top - z_center; with cell index k counting from
    // the bottom, that is (Nz - 0.5 - k) * dz.
    void InitPhi(int lev)
    {
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PHI = geom[lev].ProbHi();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar z_top = PHI[2];
        for (amrex::MFIter mfi(*phi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar zc = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                phi(i, j, k) = z_top - zc;
                rem(i, j, k) = 0.0;
            });
        }
        regime_mf[lev]->setVal(0.0);
        RoP_mf[lev]->setVal(0.0);
        if (spall_enabled)
        {
            spall_event_mf[lev]->setVal(0.0);
            spall_thickness_mf[lev]->setVal(0.0);
            RoP_spall_mf[lev]->setVal(0.0);
        }
        if (vapor_enabled)
        {
            vapor_event_mf[lev]->setVal(0.0);
            vapor_thickness_mf[lev]->setVal(0.0);
            RoP_vap_mf[lev]->setVal(0.0);
        }
    }

    void RefreshMicrostructureDerivedFields(bool announce)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            InitSurfaceMask(lev);
            InitializeMicrostructure(lev, announce);
        }
    }

    void RepairMechanicsAfterRegrid(Set::Scalar time)
    {
        if (m_type == Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            return;

        if (max_level > 0)
            RefreshMicrostructureDerivedFields(/*announce=*/false);

        UpdateModel(/*a_step=*/0, time);

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (m_zero_out_displacement)
                disp_mf[lev]->setVal(Set::Vector::Zero());
            stress_mf[lev]->setVal(Set::Matrix::Zero());
            strain_mf[lev]->setVal(Set::Matrix::Zero());

            Util::RealFillBoundary(*disp_mf[lev], geom[lev]);
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);
            Util::RealFillBoundary(*strain_mf[lev], geom[lev]);
        }
    }


    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar AlphaBetaFraction(Set::Scalar T, Set::Scalar T0,
                                         Set::Scalar width)
    {
        if (width <= 0.0) return (T >= T0) ? 1.0 : 0.0;
        Set::Scalar lo = T0 - 0.5 * width;
        Set::Scalar s = (T - lo) / width;
        if (s <= 0.0) return 0.0;
        if (s >= 1.0) return 1.0;
        return s * s * (3.0 - 2.0 * s);
    }

    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ThermalFreeStrain(Set::Scalar T,
                                         Set::Scalar beta,
                                         Set::Scalar T_ref,
                                         Set::Scalar alpha_beta_eps,
                                         Set::Scalar alpha_beta_T0,
                                         Set::Scalar alpha_beta_width)
    {
        return beta * (T - T_ref)
             + alpha_beta_eps * AlphaBetaFraction(T, alpha_beta_T0,
                                                  alpha_beta_width);
    }

    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar CellToNodeAverageInDomain(const amrex::Array4<const Set::Scalar>& f,
                                          int i, int j, int k, int comp,
                                          int ilo, int ihi,
                                          int jlo, int jhi,
                                          int klo, int khi,
                                          Set::Scalar max_abs)
    {
        Set::Scalar sum = 0.0;
        int count = 0;
        for (int dk = 0; dk <= 1; ++dk)
        {
            int kk = k - dk;
            if (kk < klo || kk > khi) continue;
            for (int dj = 0; dj <= 1; ++dj)
            {
                int jj = j - dj;
                if (jj < jlo || jj > jhi) continue;
                for (int di = 0; di <= 1; ++di)
                {
                    int ii = i - di;
                    if (ii < ilo || ii > ihi) continue;
                    Set::Scalar v = f(ii, jj, kk, comp);
                    if (!std::isfinite(v) || std::abs(v) > max_abs) continue;
                    sum += v;
                    ++count;
                }
            }
        }
        return (count > 0) ? sum / Set::Scalar(count) : 0.0;
    }

    void UpdateMicrostructureFreeStrain(int lev)
    {
        for (MFIter mfi(*thermal_free_strain_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& temp = temp_old_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& beta = beta_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& Tref = T_ref_phase_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& abe  = alpha_beta_linear_strain_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& abT0 = alpha_beta_T0_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& abW  = alpha_beta_width_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& eps = thermal_free_strain_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                eps(i, j, k) = ThermalFreeStrain(
                    temp(i, j, k, 0),
                    beta(i, j, k, 0),
                    Tref(i, j, k, 0),
                    abe (i, j, k, 0),
                    abT0(i, j, k, 0),
                    abW (i, j, k, 0));
            });
        }
        Util::RealFillBoundary(*thermal_free_strain_mf[lev], geom[lev]);
    }

    // Step 5: build the Voronoi grain field, mineral assignment, heterogeneous
    // property fields, and grain-boundary flag for level `lev`. Grain centres
    // and the grain→phase mapping are generated once (on the first call) using
    // a deterministic mt19937 seeded from `microstructure.seed`, so the result
    // is reproducible and consistent across AMR levels.
    void InitializeMicrostructure(int lev, bool announce = true)
    {
        // Static-field ghost layers must be refilled on the next advance.
        if (lev >= static_cast<int>(micro_ghost_filled.size()))
            micro_ghost_filled.resize(lev + 1, 0);
        micro_ghost_filled[lev] = 0;
        // Step 6b: build per-phase Constant material and H<->T table once
        // (lazy on first level call). Each phase's Constant captures the
        // user-specified rho, Cp, kappa plus the optional phase-change
        // scalars; the Table integrates rho*cp_eff to give H(T) and the
        // inverse, with phase fractions inferred from H. (A law phase's
        // Table is built from its law instead: see SITE 4 below.)
        if (phase_tables.size() != phases.size())
        {
            phase_materials.clear();
            phase_tables.clear();
            phase_materials.resize(phases.size());
            phase_tables.resize(phases.size());
            for (size_t p = 0; p < phases.size(); ++p)
            {
                phase_materials[p].m_rho   = phases[p].rho;
                phase_materials[p].m_cp    = phases[p].Cp;
                phase_materials[p].m_kappa = phases[p].kappa;
                phase_materials[p].T_ref     = phases[p].T_ref;
                phase_materials[p].T_m       = phases[p].T_m;
                phase_materials[p].fus_width = phases[p].fus_width;
                phase_materials[p].L_m       = phases[p].L_m;
                phase_materials[p].T_vap_lo  = phases[p].T_vap_lo;
                phase_materials[p].T_vap_hi  = phases[p].T_vap_hi;
                phase_materials[p].L_v       = phases[p].L_v;
                // ⚠⚠ P2-2b SITE 4, resolution (ii): BUILD FROM THE LAW, never
                // from the scalar. `Table::Build` takes `const Material&` and is
                // polymorphic (`Table.H:30`), so a Grimsel object goes straight
                // in and the table then integrates rho(T)*c_p(T) instead of a
                // constant -- which is what routes rho and c_p into conduction
                // AND into the H->T inversion, from ONE curve.
                //
                // ⚠⚠ THE LATENT HEAT MUST BE COPIED ACROSS FIRST. Grimsel parses
                // none of Material's reference / latent scalars, so `T_ref`/
                // `T_m`/`fus_width`/`L_m`/`T_vap_lo`/`T_vap_hi`/`L_v` would
                // sit at Material's base defaults and `effective_rhocp`'s
                // boxcars would contribute NOTHING -- the latent heat vanishes
                // silently. `unit/latent_heat` is the witness (the shared rock
                // file, validation/shared/rock_grimsel.inputs, declares both
                // bands under law = grimsel).
                if (phases[p].grimsel_law)
                {
                    Numeric::Material::Grimsel g = phases[p].grimsel;
                    g.T_ref     = phases[p].T_ref;
                    g.T_m       = phases[p].T_m;
                    g.fus_width = phases[p].fus_width;
                    g.L_m       = phases[p].L_m;
                    g.T_vap_lo  = phases[p].T_vap_lo;
                    g.T_vap_hi  = phases[p].T_vap_hi;
                    g.L_v       = phases[p].L_v;
                    phase_tables[p].Build(g);
                }
                // P2-29: a zhang phase -- the same rule, one rule for both laws (Zhang's own latent
                // scalars are refused at the parse; the phase's are the ones that count)
                else if (phases[p].zhang_law)
                {
                    Numeric::Material::Zhang z = phases[p].zhang;
                    z.T_ref     = phases[p].T_ref;
                    z.T_m       = phases[p].T_m;
                    z.fus_width = phases[p].fus_width;
                    z.L_m       = phases[p].L_m;
                    z.T_vap_lo  = phases[p].T_vap_lo;
                    z.T_vap_hi  = phases[p].T_vap_hi;
                    z.L_v       = phases[p].L_v;
                    phase_tables[p].Build(z);
                }
                else
                    phase_tables[p].Build(phase_materials[p]);
                // P2-24: THE TABLE THE RUN INVERTS, printed once per build when a latent band is declared: its
                // band temperatures and enthalpies [J/m^3, H(T_ref) = 0] (b1 reads the fusion jump from here)
                if ((phases[p].T_m < 1.0e29 || phases[p].T_vap_lo < 1.0e29)
                    && amrex::ParallelDescriptor::IOProcessor())
                {
                    const Numeric::Material::Table& tb = phase_tables[p];
                    const Set::Scalar fl = phases[p].T_m - 0.5 * phases[p].fus_width;
                    const Set::Scalar fh = phases[p].T_m + 0.5 * phases[p].fus_width;
                    std::cout << std::setprecision(17)
                              << "P2-24 phase " << p << " enthalpy table [J/m^3]: H(fus_lo " << fl << " K) = "
                              << tb.H_fus_lo() << ", H(fus_hi " << fh << " K) = " << tb.H_fus_hi()
                              << ", H(vap_lo " << phases[p].T_vap_lo << " K) = " << tb.H_vap_lo()
                              << ", H(vap_hi " << phases[p].T_vap_hi << " K) = " << tb.H_vap_hi()
                              << std::setprecision(6) << std::endl;
                }
            }
        }

        // Lazy one-shot generation of grain centres + per-grain mineral
        // assignment. We do the tessellation locally rather than going through
        // IC::Voronoi because IC::Voronoi::Parse reads its `seed` member only
        // after Define() runs, so its grain placement isn't actually under
        // user-seed control. Doing it inline keeps the microstructure fully
        // reproducible from `microstructure.seed`. Skipped entirely in
        // expression mode — Step 6 prescribes the phase IDs from a
        // user-supplied amrex::Parser.
        if (micro_mode == MicroMode::Voronoi && grain_centres.empty())
        {
            std::mt19937 rng(static_cast<unsigned int>(micro_seed));
            std::uniform_real_distribution<Set::Scalar> uni(0.0, 1.0);
            const Set::Scalar* lo = geom[0].ProbLo();
            const Set::Scalar* hi = geom[0].ProbHi();
            grain_centres.assign(3 * num_grains, 0.0);
            for (int g = 0; g < num_grains; ++g)
            {
                grain_centres[3*g + 0] = lo[0] + (hi[0] - lo[0]) * uni(rng);
                grain_centres[3*g + 1] = lo[1] + (hi[1] - lo[1]) * uni(rng);
                grain_centres[3*g + 2] = lo[2] + (hi[2] - lo[2]) * uni(rng);
            }
            std::vector<Set::Scalar> cdf(phases.size(), 0.0);
            Set::Scalar acc = 0.0;
            for (size_t p = 0; p < phases.size(); ++p)
            {
                acc += phases[p].fraction;
                cdf[p] = acc;
            }
            grain_phase_ids.assign(num_grains, 0);
            for (int g = 0; g < num_grains; ++g)
            {
                Set::Scalar u = uni(rng);
                int pid = (int)phases.size() - 1;
                for (size_t p = 0; p < phases.size(); ++p)
                    if (u <= cdf[p]) { pid = (int)p; break; }
                grain_phase_ids[g] = pid;
            }
        }

        // Assign each cell to a phase ID. Ghost cells get a -1 sentinel so
        // out-of-domain ghosts are identifiable in the is_gb pass. The kernel
        // branches on `micro_mode`: Voronoi looks up the nearest grain centre;
        // Expression evaluates the user parser at the cell centre and rounds
        // to the nearest integer phase index.
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const int n_phases = static_cast<int>(phases.size());
        phase_mf[lev]->setVal(-1.0, 0, 1, phase_mf[lev]->nGrow());
        grain_id_mf[lev]->setVal(-1.0, 0, 1, grain_id_mf[lev]->nGrow());
        if (micro_mode == MicroMode::Voronoi)
        {
            const int n_grains = num_grains;
            const Set::Scalar* gx = grain_centres.data(); // packed [num_grains*3]
            const int* gpid = grain_phase_ids.data();
            for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& gid = (*grain_id_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                    Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                    Set::Scalar best_d2 = std::numeric_limits<Set::Scalar>::infinity();
                    int best_g = 0;
                    for (int g = 0; g < n_grains; ++g)
                    {
                        Set::Scalar dx_ = x - gx[3*g + 0];
                        Set::Scalar dy_ = y - gx[3*g + 1];
                        Set::Scalar dz_ = z - gx[3*g + 2];
                        Set::Scalar d2 = dx_*dx_ + dy_*dy_ + dz_*dz_;
                        if (d2 < best_d2) { best_d2 = d2; best_g = g; }
                    }
                    // Step 6c: grain_id is the unique grain index; phase is
                    // the per-grain mineral assignment, so two adjacent
                    // same-mineral grains have grain_id != grain_id but
                    // identical phase.
                    gid(i, j, k) = static_cast<Set::Scalar>(best_g);
                    ph (i, j, k) = static_cast<Set::Scalar>(gpid[best_g]);
                });
            }
        }
        else // Expression
        {
            auto pf = phase_f;
            const bool has_gexpr = has_grain_expr;
            auto gf = grain_f;
            for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& gid = (*grain_id_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                    Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                    int p = static_cast<int>(std::round(pf(x, y, z)));
                    if (p < 0)         p = 0;
                    if (p >= n_phases) p = n_phases - 1;
                    ph(i, j, k) = static_cast<Set::Scalar>(p);
                    // Step 6c: grain_id either follows a user expression or
                    // defaults to the phase id. Expression values are
                    // rounded to a non-negative integer label and never used
                    // as a vector index.
                    int g = p;
                    if (has_gexpr)
                    {
                        int gv = static_cast<int>(std::round(gf(x, y, z)));
                        if (gv < 0) gv = 0;
                        g = gv;
                    }
                    gid(i, j, k) = static_cast<Set::Scalar>(g);
                });
            }
        }

        // Property lookup tables snapshotted onto raw pointers so the
        // ParallelFor capture works with std::vector storage outside the
        // lambda. (Lambdas here intentionally run host-side; properties are
        // a per-cell scalar lookup with no fast-path benefit from GPU.)
        std::vector<Set::Scalar> kap(n_phases), bet(n_phases), Ev(n_phases),
                                 muv(n_phases), Trv(n_phases),
                                 ab_eps(n_phases), ab_T0(n_phases),
                                 ab_width(n_phases), rh(n_phases),
                                 Cpv(n_phases);
        for (int p = 0; p < n_phases; ++p)
        {
            kap[p] = phases[p].kappa;
            bet[p] = phases[p].beta;
            Ev[p]  = phases[p].E;
            muv[p] = phases[p].mu;
            Trv[p] = phases[p].T_ref;
            ab_eps[p]   = phases[p].alpha_beta_linear_strain;
            ab_T0[p]    = phases[p].alpha_beta_T0;
            ab_width[p] = phases[p].alpha_beta_width;
            // ⚠⚠ P2-2b SITE 5 -- THE ONE LEGITIMATE SCALAR -> FIELD WRITER for
            // rho and Cp: this is where the per-phase values become per-cell
            // ones. unit/lateral_flux_lock permits a `phases[..]` rho / Cp
            // scalar read only here and at the constant-law table object above
            // (phase_materials); CheckTimeStepLimit also reads a constant-law
            // phase's scalars, for the explicit limit.
            // ⚠ rho_mf / Cp_mf are written and plotted (`rho_phase`,
            // `Cp_phase`) and read back by no code: conduction takes rho c_p
            // through the phase's H<->T table.
            //
            // ⚠⚠ DO NOT DELETE THIS IN FAVOUR OF `UpdatePropertyFields` ALONE.
            // That function returns early unless SOME phase sets
            // `law = grimsel` or `zhang`, and skips constant-law phases.
            // Under the shipped default the fields would then be NEVER
            // WRITTEN. Together the two cover every phase under every law:
            // this at init, that at init and per step.
            //
            // ⚠ LAW-AWARE: a law phase's rho / Cp scalars are refused at the
            // parse and stay 0, so its init values here are the law's at
            // T_ref. UpdatePropertyFields(lev, false) below then rewrites them
            // at the cell temperature wherever that is positive and finite, so
            // the fields hold the law's values from Initialize on.
            if (phases[p].grimsel_law)
            {
                rh[p]  = phases[p].grimsel.rho(phases[p].T_ref);
                Cpv[p] = phases[p].grimsel.cp (phases[p].T_ref);
            }
            else if (phases[p].zhang_law)   // P2-29
            {
                rh[p]  = phases[p].zhang.rho(phases[p].T_ref);
                Cpv[p] = phases[p].zhang.cp (phases[p].T_ref);
            }
            else
            {
                rh[p]  = phases[p].rho;
                Cpv[p] = phases[p].Cp;
            }
        }
        const Set::Scalar* d_kap = kap.data();
        const Set::Scalar* d_bet = bet.data();
        const Set::Scalar* d_Ev  = Ev.data();
        const Set::Scalar* d_mu  = muv.data();
        const Set::Scalar* d_Tr  = Trv.data();
        const Set::Scalar* d_ab_eps = ab_eps.data();
        const Set::Scalar* d_ab_T0  = ab_T0.data();
        const Set::Scalar* d_ab_w   = ab_width.data();
        const Set::Scalar* d_rh  = rh.data();
        const Set::Scalar* d_Cpv = Cpv.data();

        for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& kapa = (*kappa_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& beta = (*beta_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& EE   = (*E_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& muu  = (*mu_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Trf  = (*T_ref_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& abe  = (*alpha_beta_linear_strain_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& abT0 = (*alpha_beta_T0_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& abW  = (*alpha_beta_width_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rho_a= (*rho_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Cp_a = (*Cp_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                kapa (i, j, k) = d_kap[p];
                beta (i, j, k) = d_bet[p];
                EE   (i, j, k) = d_Ev [p];
                muu  (i, j, k) = d_mu [p];
                Trf  (i, j, k) = d_Tr [p];
                abe  (i, j, k) = d_ab_eps[p];
                abT0 (i, j, k) = d_ab_T0 [p];
                abW  (i, j, k) = d_ab_w  [p];
                rho_a(i, j, k) = d_rh [p];
                Cp_a (i, j, k) = d_Cpv[p];
            });
        }
        // Step 7: fill ghost cells of mechanical-property fields so cell-to-node
        // averaging at box boundaries reads coherent neighbour values.
        Util::RealFillBoundary(*beta_mf[lev],        geom[lev]);
        Util::RealFillBoundary(*E_phase_mf[lev],     geom[lev]);
        Util::RealFillBoundary(*mu_phase_mf[lev],    geom[lev]);
        Util::RealFillBoundary(*T_ref_phase_mf[lev], geom[lev]);
        Util::RealFillBoundary(*alpha_beta_linear_strain_mf[lev], geom[lev]);
        Util::RealFillBoundary(*alpha_beta_T0_mf[lev],            geom[lev]);
        Util::RealFillBoundary(*alpha_beta_width_mf[lev],         geom[lev]);
        // P2-27 stage 4: a LAW phase's scalars are refused (Parse), so its property
        // fields come from the one per-cell writer -- the law at the start
        // temperature, what the first Advance writes -- not from the scalars
        // above (no-op unless some phase sets law = grimsel or zhang; no counter moves).
        UpdatePropertyFields(lev, /*step=*/false);
        UpdateMicrostructureFreeStrain(lev);

        // Fill phase_mf and grain_id_mf ghost cells from neighbouring boxes;
        // out-of-domain ghosts retain the -1 sentinel set above.
        Util::RealFillBoundary(*phase_mf[lev],    geom[lev]);
        Util::RealFillBoundary(*grain_id_mf[lev], geom[lev]);

        const amrex::Box& domain = geom[lev].Domain();
        const int dlo[3] = {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
        const int dhi[3] = {domain.bigEnd(0),   domain.bigEnd(1),   domain.bigEnd(2)};
        // Step 6c: build phase- and grain-boundary flags from the same face
        // sweep. is_gb is a phase-boundary alias (the same value as
        // is_phase_boundary; unit/voronoi and unit/grain_topology read it);
        // the flaw sampler (InitializeFlaws) consumes is_grain_boundary.
        for (amrex::MFIter mfi(*is_gb_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gid = (*grain_id_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& gb  = (*is_gb_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& pbf = (*is_phase_boundary_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& gbf = (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                int p_self = static_cast<int>(ph(i, j, k));
                int g_self = static_cast<int>(gid(i, j, k));
                int p_flag = 0, g_flag = 0;
                const int di[6] = {-1, 1, 0, 0, 0, 0};
                const int dj[6] = { 0, 0,-1, 1, 0, 0};
                const int dk[6] = { 0, 0, 0, 0,-1, 1};
                for (int n = 0; n < 6; ++n)
                {
                    int ii = i + di[n], jj = j + dj[n], kk = k + dk[n];
                    // Skip neighbours that fall outside the physical domain.
                    if (ii < dlo[0] || ii > dhi[0] ||
                        jj < dlo[1] || jj > dhi[1] ||
                        kk < dlo[2] || kk > dhi[2]) continue;
                    int p_n = static_cast<int>(ph (ii, jj, kk));
                    int g_n = static_cast<int>(gid(ii, jj, kk));
                    if (p_n >= 0 && p_n != p_self) p_flag = 1;
                    if (g_n >= 0 && g_n != g_self) g_flag = 1;
                    if (p_flag && g_flag) break;
                }
                pbf(i, j, k) = static_cast<Set::Scalar>(p_flag);
                gbf(i, j, k) = static_cast<Set::Scalar>(g_flag);
                gb (i, j, k) = static_cast<Set::Scalar>(p_flag);
            });
        }


        // Step 15b: sample the per-cell Weibull flaw-length field. Runs
        // after is_grain_boundary_mf is populated above; default-off, only
        // active under spallation.model = sp_weibull + weibull.enabled = 1.
        InitializeFlaws(lev);
    }

    // Step 15b: deterministic per-cell Weibull flaw-length sampling.
    // a_f = a_scale * (-ln U)^(1/m), a 2-parameter Weibull with scale
    // a_scale and shape m; a_scale = a0_gb on is_grain_boundary cells,
    // a0_ig (times the size factor) otherwise. U is a splitmix64 hash of
    // weibull.seed and the GLOBAL cell index (i, j, k) on an intragranular
    // cell, of weibull.seed and the sorted grain pair of each adjacent facet
    // on a grain-boundary cell (see below), so flaw_a_mf is
    // bit-identical regardless of MPI rank count or BoxArray partition.
    // AMR regrid-repair of flaw_a_mf is a §15c deliverable; this resamples
    // deterministically on any call, which is single-level-safe.
    void InitializeFlaws(int lev)
    {
        if (!weibull_enabled) return;

        const Set::Scalar a0_gb = weibull_a0_gb;
        const Set::Scalar a0_ig = weibull_a0_ig;
        const Set::Scalar inv_m = 1.0 / weibull_m;
        const std::uint64_t seed =
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(weibull_seed));

        // Mesh-objective flaw assignment, split by physical population:
        //
        //   GB cells -> per-FACET. One Weibull draw per adjacent grain pair
        //     (gA, gB), keyed by a hash of the SORTED pair + seed. Grains are
        //     seed-fixed, so a facet draws the same flaw at any resolution and
        //     the facet population (hence the max surface flaw that sets onset,
        //     ΔT_onset ∝ a^(-1/2)) is mesh-invariant. A cell touching several
        //     grains takes the MAX (weakest-link) facet flaw. No vol_factor:
        //     the facet, not the cell, is the physical unit.
        //   IG cells -> per-CELL Weibull (hash of i,j,k + seed) scaled by the
        //     weakest-link size factor (V_cell/V_0)^(1/m); intragranular cracks
        //     are a volumetric population. V_0 (default (1 mm)^3) makes this
        //     exactly 1 at the calibration mesh, so 1 mm IG draws are identical.
        //
        // This removes the per-cell GB tag's mesh-dependence (tagged GB
        // fraction ~ dx/grain -> 0 as dx->0, which drove onset toward the wrong
        // all-IG limit). grain_id carries 1 ghost (ngrow=1) and is ghost-filled
        // below, so the 6-face neighbour grain lookup is valid across ranks;
        // the facet hash needs no MPI communication.
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar V_cell = DX[0] * DX[1] * DX[2];
        const Set::Scalar vol_factor = std::pow(V_cell / weibull_V0, inv_m);

        // grain_id ghosts must be current for the 6-face neighbour lookup.
        Util::RealFillBoundary(*grain_id_mf[lev], geom[lev]);
        const amrex::Box& domain = geom[lev].Domain();
        const auto dlo = domain.smallEnd();
        const auto dhi = domain.bigEnd();

        // Uniform U in (0, 1) from a 64-bit hash (53-bit mantissa).
        auto u_of = [] (std::uint64_t key) -> Set::Scalar {
            return (static_cast<Set::Scalar>(key >> 11) + 0.5)
                   * (1.0 / 9007199254740992.0);
        };
        auto u32 = [] (int v) -> std::uint64_t {
            return static_cast<std::uint64_t>(static_cast<std::uint32_t>(v));
        };

        for (amrex::MFIter mfi(*flaw_a_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& fa = (*flaw_a_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gb =
                (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gid =
                (*grain_id_mf[lev]).array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (gb(i, j, k) > 0.5)
                {
                    // GRAIN-BOUNDARY cell -> per-facet, weakest adjacent facet.
                    const int g_self = static_cast<int>(gid(i, j, k));
                    const int di[6] = {-1, 1, 0, 0, 0, 0};
                    const int dj[6] = { 0, 0,-1, 1, 0, 0};
                    const int dk[6] = { 0, 0, 0, 0,-1, 1};
                    Set::Scalar a_max = 0.0;
                    for (int n = 0; n < 6; ++n)
                    {
                        int ii = i + di[n], jj = j + dj[n], kk = k + dk[n];
                        if (ii < dlo[0] || ii > dhi[0] ||
                            jj < dlo[1] || jj > dhi[1] ||
                            kk < dlo[2] || kk > dhi[2]) continue;
                        int g_n = static_cast<int>(gid(ii, jj, kk));
                        if (g_n < 0 || g_n == g_self) continue;
                        // Sorted pair -> both sides of the facet draw alike.
                        int gA = g_self < g_n ? g_self : g_n;
                        int gB = g_self < g_n ? g_n : g_self;
                        std::uint64_t key = SplitMix64(seed);
                        key = SplitMix64(key ^ u32(gA));
                        key = SplitMix64(key ^ u32(gB));
                        Set::Scalar a_f =
                            a0_gb * std::pow(-std::log(u_of(key)), inv_m);
                        if (a_f > a_max) a_max = a_f;
                    }
                    if (a_max <= 0.0)
                    {
                        // A GB-tagged cell should always have a different-grain
                        // face neighbour; if not (shouldn't happen), fall back
                        // to a per-cell GB draw so flaw_a is never left at 0.
                        std::uint64_t key = SplitMix64(seed);
                        key = SplitMix64(key ^ u32(i));
                        key = SplitMix64(key ^ u32(j));
                        key = SplitMix64(key ^ u32(k));
                        a_max = a0_gb * std::pow(-std::log(u_of(key)), inv_m);
                    }
                    fa(i, j, k) = a_max;
                }
                else
                {
                    // INTRAGRANULAR cell -> per-cell volumetric (size effect).
                    std::uint64_t key = SplitMix64(seed);
                    key = SplitMix64(key ^ u32(i));
                    key = SplitMix64(key ^ u32(j));
                    key = SplitMix64(key ^ u32(k));
                    fa(i, j, k) =
                        a0_ig * vol_factor * std::pow(-std::log(u_of(key)), inv_m);
                }
            });
        }
    }

    // Conserved enthalpy and phase fractions (own to MMWSpalling; T comes from
    // HeatConduction's temp_mf / temp_old_mf).
    Set::Field<Set::Scalar> H_mf, H_old_mf;
    Set::Field<Set::Scalar> Lambda_S_mf, Lambda_L_mf, Lambda_V_mf;
    Set::Field<Set::Scalar> surface_mf; // 1.0 on surface cells, 0.0 otherwise

    Set::Scalar eps_rad      = 0.0;     // emissivity below T_eps_break  (-)
    Set::Scalar eps_rad_hi   = 0.0;     // emissivity above T_eps_break  (-)
    Set::Scalar T_eps_break  = 1.0e30;  // breakpoint for ε(T)           (K)
    Set::Scalar sigma_SB = 0.0; // Stefan-Boltzmann constant       (W/m²K⁴)
    Set::Scalar h_conv   = 0.0; // convective heat-transfer coeff. (W/m²K)
    Set::Scalar T_amb    = 0.0; // ambient temperature             (K)

    Numeric::MMWBeam beam; // Beer-Lambert + Gaussian volumetric source (step 3)

    // E1: beam-deposition closure counter (BeamEnergyBalance), written to
    // thermo.dat as non-extensive vars (reduced by BeamEnergyBalance itself).
    // Powers [W] are per step; energies [J] and the violation count are
    // cumulative. The surface-owning level only (level 0 on a single-level run).
    Set::Scalar beam_P_inc          = 0.0;
    Set::Scalar beam_P_dep          = 0.0;
    Set::Scalar beam_P_trans        = 0.0;
    Set::Scalar beam_E_inc          = 0.0;
    Set::Scalar beam_E_dep          = 0.0;
    Set::Scalar beam_E_trans        = 0.0;
    Set::Scalar beam_closure_err    = 0.0;
    Set::Scalar beam_invariant_viol = 0.0;
    long        beam_viol_steps     = 0;

    // P2-21: the deciders, thermo columns registered ONLY when
    // beam.alpha_varies() && beam.can_deposit() (every flat and every beam-OFF
    // deck keeps its header). Per step, non-extensive; set by ColumnOpticalDepth.
    Set::Scalar beam_tau_max     = 0.0;  // max over lit columns of tau_bot at the last solid cell
    Set::Scalar beam_alpha_max   = 0.0;  // max a_k over lit solid cells [1/m]
    Set::Scalar beam_col_err_max = 0.0;  // max over lit columns |sum_k (e^-tt - e^-tb) + e^-tN - 1|
    Set::Scalar bc_T_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();  // central column [K]
    Set::Scalar bc_z_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();  // [m below the top face]
    Set::Scalar bc_z_melt = std::numeric_limits<Set::Scalar>::quiet_NaN();  // shallowest T >= T_m [m]
    Set::Scalar bc_tau    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // central tau_bot at k_last
    Set::Scalar bc_ell_p  = std::numeric_limits<Set::Scalar>::quiet_NaN();  // P2-19 face e-fold [m]
    Set::Scalar bc_ell_pk = std::numeric_limits<Set::Scalar>::quiet_NaN();  // P2-19 peak e-fold [m]
    // P2-24: the melt and vapour record (same registration condition as the bc_* above)
    Set::Scalar beam_P_beam = 0.0;   // the schedule's power now [W] (the DOE's denominator; P_inc is after R)
    Set::Scalar bc_T_rad    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // synthetic radiometer eps T [K]
    Set::Scalar bc_T_view   = std::numeric_limits<Set::Scalar>::quiet_NaN();  // view mean of T_top [K]
    Set::Scalar bc_A_melt   = 0.0;   // top solid cells with f_L > 0 [m^2]
    Set::Scalar bc_f_L_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();  // f_L at the central peak [-]
    Set::Scalar bc_z_liq    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // deepest central f_L >= 1 [m]
    Set::Scalar bc_R_col    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // the central column's R_col [-]
    Set::Scalar bc_n_vapor     = 0.0;   // cumulative surface cells vaporised (the face list's vapour rule)
    Set::Scalar bc_n_vapor_int = 0.0;   // cumulative cell-steps an INTERIOR cell sat past H_vap_hi

    // S1: the mask-following surface -- the surface patch and the surface
    // losses sit at the `removed`-mask column top (ColumnTopSolid) instead of
    // the fixed domain top / the phi-window surface_mf. Not a key: Parse sets
    // it to surface_patch.face_list (on exactly under the face list).
    bool surface_follow_mask = false;

    // S1: whole-domain energy ledger (EnergyLedger), thermo.dat vars.
    // E_H [J] is sum(H*V) over all level-0 cells (void included), E_src [J]
    // the cumulative applied energy, P_* [W] per step. The compensated-sum
    // state (ledger_E_H0_s/c) is internal.
    bool        ledger_enabled               = false;
    // D2a: per-column flame (h, T_gas) and nozzle stand-off s for the current
    // AdvanceMicrostructure call (BuildFlameColumns). Constants path: h_conv /
    // T_flame everywhere. Read by the kernel, EnergyLedger and the face list.
    std::vector<Set::Scalar> flame_h, flame_Tg, flame_s;
    Set::Scalar patch_min_standoff = 0.0; // thermo: min s over in-patch live columns [m]
    // D2a2: jet_closure = enthalpy diagnostics (JetEnthalpyMarch).
    Set::Scalar jet_P_face    = 0.0; // thermo: enthalpy taken from the jet stream [W]
    Set::Scalar jet_P_cap     = 0.0; // thermo: m*cp*(T_stag - T_ent) [W]
    Set::Scalar jet_T_stag    = 0.0; // thermo: gas temperature entering the axis bin [K]
    Set::Scalar jet_s_c       = 0.0; // thermo: mean stand-off of the axis bin [m]
    Set::Scalar jet_P_decay   = 0.0; // thermo: nozzle excess discarded by the decay line [W]
    Set::Scalar jet_T_exhaust = 0.0; // thermo: T_gas after the outermost bin [K]
    Set::Scalar jet_P_exhaust = 0.0; // thermo: m*cp*(T_exhaust - T_ent) [W]
    Set::Scalar jet_r_reach   = 0.0; // thermo: max column r in a bin with P_bin > 0 [m]
    // D2b: exhaust recirculation (jet_T_ent_mode = exhaust). jet_T_rec = the
    // entrained-gas temperature used this step (= the previous step's
    // jet_T_mouth, which is its jet_T_exhaust when no free-surface bin exists;
    // jet_T_ent on the first step). Host state, NOT checkpointed.
    Set::Scalar jet_T_rec     = 0.0; // thermo [K]
    Set::Scalar jet_P_recirc  = 0.0; // thermo: (m - mdot)*cp*(T_rec - T_ent) [W]
    bool        jet_rec_init  = false;
    Set::Scalar jet_T_rec_next = 0.0;
    // D2b: burner-on-feet state (nozzle_descent = feet). nozzle_z is host
    // state, NOT checkpointed (a restart re-lands the burner on the feet).
    bool        feet_init           = false;
    Set::Scalar nozzle_z            = 0.0; // thermo: nozzle exit height [m]
    Set::Scalar foot_z              = 0.0; // thermo: foot (tripod) height [m]
    Set::Scalar foot_cols           = 0.0; // thermo: live annulus columns
    Set::Scalar foot_carry_cols     = 0.0; // thermo: columns at or above their pad height
    Set::Scalar foot_stall_time     = 0.0; // thermo: simulated time since z_n last decreased [s]
    // D2c: march options and burner-body clearance. Host state (jet_fs_z0,
    // foot_clear_k) is NOT checkpointed, like the D2b state above: a restart
    // re-takes the original surface from the restart step's faces.
    Set::Scalar jet_D_e      = 0.0; // thermo: decay diameter used this step [m] (momentum)
    Set::Scalar jet_fs_cols  = 0.0; // thermo: marching free-surface columns
    Set::Scalar jet_m_ratio  = 0.0; // thermo: m_exit/m_stag
    Set::Scalar jet_T_mouth  = 0.0; // thermo: inlet T_gas of the first free-surface bin [K]
    Set::Scalar jet_P_neg    = 0.0; // thermo: sum of negative q_robin*dA [W] (count)
    std::vector<Set::Scalar> jet_fs_z0;   // per column face height at the first march [m]
    Set::Scalar foot_mech_vol  = 0.0;     // thermo: cumulative mechanically removed volume [m^3]
    // P2-2f: cumulative enthalpy that left the SOLID set with pad-cleared rock
    // [J]. Credited to ledger_err_solid; kept OUT of ledger_E_flake so the
    // heat-driven price (floor + wall) excludes mechanical removal on both
    // sides of the ratio. See FrFootClearance (FaceRemoval.H) for why this is a
    // booking-only term and why the field must not be written.
    Set::Scalar foot_mech_E    = 0.0;     // thermo: cumulative pad-clearance enthalpy [J]
    Set::Scalar foot_mech_cols = 0.0;     // thermo: columns cleared this step
    std::vector<int> foot_clear_k;        // per column target top k (INT_MIN = none)
    // P2-2 change 1: how many times the property fields were refreshed.
    Set::Scalar prop_refresh_n  = 0.0;
    // (b5)'s two-temperature probe reads these: the MEAN refreshed
    // property over the level, beside the mean temperature it was
    // evaluated at. A refresh that never happens leaves them at 0.
    Set::Scalar prop_kappa_mean = 0.0, prop_beta_mean = 0.0;
    Set::Scalar prop_E_mean     = 0.0, prop_T_mean    = 0.0;
    // =====================================================================
    // P2-2 change 4: THE STRESS-RELIEF FACTOR R.
    // =====================================================================
    // `spallation.sp.stress_relief`, default 1.0 (identity, so every frozen
    // study stays byte-identical). It multiplies the thermoelastic driving
    // stress in ScanSigmaT (Shared.H), both branches -- the ONE T-derived
    // stress every criterion read goes through.
    //
    // ⚠⚠ WHY A RELIEF FACTOR AT ALL. With the measured properties the heated
    // layer carries 1.78x this rock's own compressive strength at Meier's onset
    // and 1.90x at Kant's, so an elastic criterion is outside its validity and
    // something relieves it. R stands in for that, bounded by the measurement:
    // R = UCS/sigma = 0.563 (Meier) / 0.527 (Kant), one material number to 7 %.
    //
    // ⚠⚠ RELIEVING ONE SITE RELIEVES THE WALL AND NOT THE FLOOR -- a family
    // that has cost this programme three packets. unit/lateral_flux_lock
    // asserts that every surviving stress expression takes exactly this
    // member, so no site can be missed and a new stress expression cannot
    // appear un-relieved.
    Set::Scalar sp_relief = 1.0;
    // P2-15: the buckling thickness condition (spall.buckle_span [m], 0 = off; face list only).
    // A fired line is withdrawn unless C_L x <= L_span (FaceRemoval.H, BuckleRatio). Witnesses
    // registered only under the key; per step, reset and added once in RecedeAlongLines.
    Set::Scalar sp_buckle_span = 0.0;
    // P2-23: spall.buckle_log_csv -- one row per line B evaluated, on
    // thermo-output steps, into <plot_file>_buckle_log.csv (FaceRemoval.H)
    bool sp_buckle_log_csv = false;
    bool sp_buckle_log_header_written = false;
    Set::Scalar bk_span_out = 0.0, bk_n_eval = 0.0, bk_n_block = 0.0, bk_n_inf = 0.0,
                bk_n_block_floor = 0.0, bk_n_block_wall = 0.0, bk_n_block_roof = 0.0,
                bk_max_pass = 0.0;
    long        jet_profile_last = -1; // last jet_profile_interval index written
    bool        jet_profile_header_written = false;

    bool        ledger_started               = false;
    Set::Scalar ledger_E_H0_s                = 0.0;
    Set::Scalar ledger_E_H0_c                = 0.0;
    Set::Scalar ledger_E_H                   = 0.0;
    Set::Scalar ledger_E_src                 = 0.0;
    Set::Scalar ledger_err_raw               = 0.0;
    Set::Scalar ledger_P_beam                = 0.0;
    Set::Scalar ledger_P_robin               = 0.0;
    Set::Scalar ledger_P_loss                = 0.0;
    Set::Scalar ledger_surface_missing_cols  = 0.0;
