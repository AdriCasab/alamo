    // Override HeatConduction::Advance with the enthalpy method: ONE heat path, AdvanceMicrostructure
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        // P2-9 rev 2 item 7: the step's wall time (a thermo column only under
        // the face list; measured always, read by nothing else)
        const double fr_step_t0 = amrex::second();
        // P2-2 change 1: the property fields, refreshed at the cell
        // temperature BEFORE the criterion and the conduction dispatch read
        // them. No-op unless a phase sets law = grimsel or (P2-29) zhang.
        UpdatePropertyFields(lev);
        Mechanics<Model::Solid::Affine::Isotropic>::Advance(lev, time, dt);
        // Step 10/11, P2-28: THE removal step -- the face list's (FaceRemoval.H), the only one: the
        // criterion on every listed line and the recession, BEFORE the H/Temp swaps below, so what it
        // writes is the "old" state the thermal advance consumes.
        RemovalStep(lev, time, dt);
        std::swap(*H_mf[lev], *H_old_mf[lev]);
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        AdvanceMicrostructure(lev, time, dt);
        // Step 6b: every path evolves H_mf as the conserved state, so
        // dispatch a single InvertHtoT here. The microstructure variant
        // selects the per-phase Material::Table by cell phase ID.
        InvertHtoT(lev);
        {
            double ts = amrex::second() - fr_step_t0;
            amrex::ParallelDescriptor::ReduceRealMax(ts);
            fr_t_step = static_cast<Set::Scalar>(ts);
        }
    }


    // Step 6b: heterogeneous-property explicit enthalpy update for the
    // microstructure path. H is the conserved variable:
    //
    //   ∂H/∂t = ∇·(κ_eff(x,D) ∇T) + Q − q_loss/dz
    //
    // where κ_eff(x,D) = κ_phase(x) · exp(−α_D · D). Face conductivities are
    // built from the harmonic mean of the cell-centred κ_eff values, with an
    // optional grain-boundary contact-conductance term `h_gb(D) = h_gb0·(1−D)`
    // added in series whenever the two cells across the face have different
    // phase IDs.
    //
    // T is recovered cell-by-cell in InvertHtoT using the per-phase
    // Numeric::Material::Table (built from the phase's Constant material).
    // For phases without phase change this reduces to T_new = T_old +
    // dt·(div+Q-loss)/(rho*Cp), matching the Step 6 behaviour exactly.
    void AdvanceMicrostructure(int lev, Set::Scalar time, Set::Scalar dt)
    {
        // Perf (2026-09-18): phase, grain_id and kappa_phase are written only
        // by InitializeMicrostructure (Initialize / Regrid), so their ghost
        // layers are filled once after each (re)initialisation instead of
        // every step. (D, refilled every step until P2-27, is deleted with the damage lineage.)
        if (lev >= static_cast<int>(micro_ghost_filled.size()))
            micro_ghost_filled.resize(lev + 1, 0);
        if (!micro_ghost_filled[lev])
        {
            Util::RealFillBoundary(*phase_mf[lev],     geom[lev]);
            Util::RealFillBoundary(*grain_id_mf[lev],  geom[lev]);
            Util::RealFillBoundary(*kappa_mf[lev],     geom[lev]);
            micro_ghost_filled[lev] = 1;
        }
        // Step 10: if spall is active, removed cells are void and must not
        // exchange heat with surviving material. Filling the ghost layer of
        // `removed` lets the face stencil below detect a void neighbour
        // through ParmParse-style ghost reads.
        // P2-28: and under the face list's evaluate-only mode (nothing is removed;
        // `removed` reads 0, the list and its column tops are built as with removal)
        const bool void_aware = list_state();
        if (void_aware) Util::RealFillBoundary(*removed_mf[lev], geom[lev]);

        // Domain extents for the face stencil's boundary-face fallback below.
        const amrex::Box& domain = geom[lev].Domain();
        const int dlo[3] = {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
        const int dhi[3] = {domain.bigEnd(0),   domain.bigEnd(1),   domain.bigEnd(2)};

        // E1: global top solid cell per (i, j) column for the beam path (see
        // the BeamSource call below). Empty unless void-aware with the beam on,
        // or (S1) with surface.follow_mask / the energy ledger on.
        // AMR: only the surface-owning level applies the beam, the surface
        // patch, the side faces and the losses (coarser levels: conduction
        // only, their covered cells are overwritten by average-down).
        const bool        owner = SurfaceOwner(lev);
        const Set::Scalar P_now = owner ? beam.power_at(time) : Set::Scalar(0.0);
        const bool follow_mask = surface_follow_mask;
        const amrex::Vector<int> k_top_col =
            (void_aware && (P_now > 0.0 || follow_mask || ledger_enabled))
            ? ColumnTopSolid(lev) : amrex::Vector<int>();
        const int* ktc = k_top_col.data();
        const int ncx = dhi[0] - dlo[0] + 1;
        // D2q-1: the beam path counts SOLID cells, not geometric distance.
        // Only needed when the beam is actually on (BeamSource returns 0 at
        // P = 0), which keeps the large beam-free runs free of the gather.
        const amrex::Vector<int> solid_above_col =
            (void_aware && P_now > 0.0 && !k_top_col.empty())
            ? ColumnSolidAboveCount(lev, k_top_col) : amrex::Vector<int>();
        const int* sab    = solid_above_col.data();
        const int  nz_dom = dhi[2] - dlo[2] + 1;
        const int  klo_d  = dlo[2];
        // P2-21: the column optical depth, through solid cells only from each
        // column's top solid cell when void-aware (the solid-count convention),
        // else from the domain top (unbuilt unless the loss table varies and the
        // beam is on: then the kernel keeps the closed form)
        const BeamTau btau = BuildBeamTau(lev, P_now, void_aware ? k_top_col : amrex::Vector<int>());
        const TauView tv = ViewOf(btau, lev);

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
        const Set::Scalar bn  = beam.n_idx;   // P2-22: in-rock spread
        const int         bprof = beam.profile;
        const Set::Scalar brad  = beam.radius;
        const Set::Scalar b_alpha = beam.alpha;   // P2-20: computed from beam.eps_r, beam.tan_delta, beam.lambda
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];
        const Set::Scalar inv_dx2 = 1.0 / (DX[0] * DX[0]);
        const Set::Scalar inv_dy2 = 1.0 / (DX[1] * DX[1]);
        const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);
        const Set::Scalar h_gb0   = h_gb_conductance;
        const Set::Scalar dx_ax[3] = {DX[0], DX[1], DX[2]};
        const Set::Scalar inv_dx  = 1.0 / DX[0];

        // Step 11b: the Hu onset validation uses a homogeneous one-phase
        // microstructure (until P2-27 so its DP damage fields existed). Apply the same
        // surface patch here as the deleted material kernel did (P2-29): prescribed_T injects the
        // Dirichlet face flux; convective_flame (Step 15a) injects the Robin
        // flame-jet flux q = h_fl*(T_flame - T_cell).
        const bool        sp_on    = surface_patch.enabled && owner;
        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
        const int         sp_topk  = geom[lev].Domain().bigEnd(2);
        const Set::Scalar sp_x0    = surface_patch.x0;
        const Set::Scalar sp_y0    = surface_patch.y0;
        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_T     = (sp_on && sp_presc)
                                   ? surface_patch.T_at(time) : 0.0;
        // D2a: one per-column (h, T_gas) evaluation per call, shared by the
        // kernel, the pin rule, the ledger and the diagnostics.
        // P2-8 F1: the face list, after the column tops and BEFORE the march
        // (BuildFlameColumns), which bills every listed face.
        // P2-23: the list runs WITHOUT a flame too. `owner` is kept (the finest
        // level owns the surface).
        // P2-28: and with a PRESCRIBED patch: its top faces carry the half-cell
        // Dirichlet term (FaceListHeat sets T_presc), billed in the face loop below.
        const bool fl_on = surface_patch.face_list && owner && void_aware;
        if (fl_on) BuildFaceList(lev, k_top_col);
        if (sp_on) BuildFlameColumns(lev, time, dt, k_top_col);
        // P2-23: with the patch OFF BuildFlameColumns -- which bills the list's
        // heat, gives every listed cell its T_s and every face its gas
        // temperature -- does not run, so its LIST half runs here, in its order,
        // without the flame: h = 0 on every face, T_s, T_gas = 0.
        else if (fl_on)
        {
            FaceListHeat(lev, 0.0, time, k_top_col);
            FrAssignTs(lev);
            FaceListTg(lev);
        }
        const Set::Scalar* fh  = sp_on ? flame_h.data()  : nullptr;
        const Set::Scalar* fTg = sp_on ? flame_Tg.data() : nullptr;
        // P2-8 F5: the list every kernel face reads.
        const FaceEntry* flp = fl_on ? fl_faces.data() : nullptr;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face;

        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> T_old = temp_old_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> H_old = H_old_mf   .Patch(lev, mfi);
            Set::Patch<Set::Scalar>       H     = H_mf       .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> gid   = grain_id_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> kp    = kappa_mf   .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> surf  = surface_mf.Patch(lev, mfi);
            // When void_aware is false, rem is bound to surface_mf as a
            // harmless placeholder that the kernel never reads.
            Set::Patch<const Set::Scalar> rem   =
                void_aware ? removed_mf.Patch(lev, mfi) : surface_mf.Patch(lev, mfi);
            amrex::Array4<const int> fli;
            if (fl_on) fli = fl_index->const_array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tc = T_old(i, j, k);
                Set::Scalar kc = kp(i, j, k);
                const bool self_void = void_aware && rem(i, j, k) > 0.5;
                if (self_void) kc = 0.0;
                int gc = static_cast<int>(gid(i, j, k));

                auto k_face = [&](int di, int dj, int dk, int axis) {
                    int ii = i + di, jj = j + dj, kk = k + dk;
                    // Domain-boundary face: kappa_mf has no BC
                    // handler so their ghost cells are uninitialised.
                    // Use the cell-centred κ_eff directly so the Dirichlet
                    // ghost in T_old still drives a flux of
                    //   q ≈ kc · (T_BC − T_cell) / dx
                    // (the standard FV boundary treatment when the BC value
                    // sits in the ghost cell at distance dx).
                    if (ii < dlo[0] || ii > dhi[0] ||
                        jj < dlo[1] || jj > dhi[1] ||
                        kk < dlo[2] || kk > dhi[2]) return kc;
                    if (void_aware && rem(ii, jj, kk) > 0.5) return Set::Scalar(0.0);
                    Set::Scalar k_n = kp(ii, jj, kk);
                    int g_n = static_cast<int>(gid(ii, jj, kk));   // Step 6c: the GB hook keys on grain_id
                    Set::Scalar inv = 0.0;
                    if (kc > 0.0 && k_n > 0.0)
                        inv = 0.5 / kc + 0.5 / k_n;
                    else
                        return Set::Scalar(0.0);
                    // Step 6c: physical grain-boundary contact resistance is
                    // keyed on true grain identity, not on mineral phase. Two
                    // adjacent same-mineral grains still see h_gb0.
                    if (h_gb0 > 0.0 && g_n != gc)
                    {
                        inv += 1.0 / (h_gb0 * dx_ax[axis]);
                    }
                    Set::Scalar kf = 1.0 / inv;
                    return kf;
                };

                Set::Scalar k_xp = k_face( 1, 0, 0, 0);
                Set::Scalar k_xm = k_face(-1, 0, 0, 0);
                Set::Scalar k_yp = k_face( 0, 1, 0, 1);
                Set::Scalar k_ym = k_face( 0,-1, 0, 1);
                Set::Scalar k_zp = k_face( 0, 0, 1, 2);
                Set::Scalar k_zm = k_face( 0, 0,-1, 2);

                Set::Scalar div_kgrad =
                    (k_xp * (T_old(i+1,j,k) - Tc) - k_xm * (Tc - T_old(i-1,j,k))) * inv_dx2 +
                    (k_yp * (T_old(i,j+1,k) - Tc) - k_ym * (Tc - T_old(i,j-1,k))) * inv_dy2 +
                    (k_zp * (T_old(i,j,k+1) - Tc) - k_zm * (Tc - T_old(i,j,k-1))) * inv_dz2;

                Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];

                // S1: with surface.follow_mask the patch and the losses sit on
                // the mask column top; otherwise the fixed domain top and the
                // phi-window surface_mf (unchanged).
                const int k_top_c = (follow_mask && ktc)
                                  ? ktc[(j - dlo[1]) * ncx + (i - dlo[0])] : sp_topk;
                const Set::Scalar surf_c = follow_mask
                                         ? Set::Scalar(k == k_top_c ? 1.0 : 0.0)
                                         : surf(i, j, k);
                bool in_patch = false;
                if (!self_void && sp_on && k == k_top_c) {
                    Set::Scalar dxp = x - sp_x0, dyp = y - sp_y0;
                    if (dxp * dxp + dyp * dyp <= sp_r2) {
                        // ⚠⚠ P2-28: THE HALF-CELL DIRICHLET TERM. The held temperature sits on the
                        // cell's top FACE, dz/2 from its centre: conductance 2 kc / dz, so the
                        // divergence term is 2 kc (T_patch - T_c) / dz^2. Until P2-28 this read
                        // kc (...) / dz^2 -- the held value a full cell away, half the boundary
                        // conductance (-12 to -20 % in the top cell's rise at 3.125 mm; first
                        // order, so refinement hid it). Under the face list the term is the
                        // LIST's (T_presc on the top face, the loop below), so it is not added
                        // here twice.
                        if (sp_presc)
                        {
                            if (!fl_on) div_kgrad += 2.0 * kc * (sp_T - Tc) * inv_dz2;
                        }
                        else
                            in_patch = true;
                    }
                }
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                const int col_c = (j - dlo[1]) * ncx + (i - dlo[0]);
                if (!self_void && owner)
                    SurfaceCellFlux(in_patch, sp_face, surf_c,
                                    fh  ? fh[col_c]  : Set::Scalar(0.0),
                                    fTg ? fTg[col_c] : Set::Scalar(0.0),
                                    kc, DX[2], Tc,
                                    eps_lo, eps_hi, T_eps_b, sig_SB, h_c, T_a,
                                    q_in, q_out, q_rob, q_los);
                // P2-8 F4: the column-top +z face keeps the surface form, x w
                // (gain AND losses: both scale with the face's area).
                int fl_f = 0, fl_n = 0;
                if (fl_on && !self_void) { fl_f = fli(i, j, k, 0); fl_n = fli(i, j, k, 1); }
                // P2-9: a LISTED cell's every face -- the top included -- is
                // heated below at the cell's one T_s, x w (FrFaceFlux).
                if (fl_n > 0) { q_in = 0.0; q_out = 0.0; }
                if (in_patch) div_kgrad += q_in * inv_dz;


                // P2-9 R1: EVERY listed face -- top, up, side, down -- at its
                // OWN h (F2), the CELL's T_s, x its weight (F3). Risers are
                // part of the surface, not surfaces of their own.
                for (int q = fl_f; q < fl_f + fl_n; ++q)
                {
                    // P2-28: a PRESCRIBED face -- the column term's own expression (2 kc
                    // (T_presc - T_c) / dz^2, the same kc, T_c and evaluation time), x w.
                    // ⚠ w FIRST, in the column term's product order: the compiler fuses the
                    // last multiply with the add (FMA), so w (P2-9's area weight; 1.0 exactly
                    // on a flat top) must not round the product before it -- a parenthesised
                    // w * (...) differed from the column term by one ulp in H at 8 of 76 cells.
                    if (flp[q].T_presc > 0.0)
                    {
                        div_kgrad += flp[q].w * 2.0 * kc * (flp[q].T_presc - Tc) * inv_dz2;
                        continue;
                    }
                    Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                    FrFaceFlux(flp[q], fTg ? flp[q].Tg : Set::Scalar(0.0), kc, Tc, DX[0],
                               a_in, a_out, a_rob, a_los);
                    div_kgrad += flp[q].w * (a_in - a_out) * inv_dx;
                }

                // E1: Beer-Lambert path is measured through material from the
                // TOP FACE OF THE FIRST SOLID CELL, so it matches the `removed`
                // mask. Removal keeps solid <=> phi >= 0 with phi linear down
                // the column, so the top solid cell has phi in [0, dz) and a
                // cell m cells below it has n_above = floor(phi/dz) == m; its
                // path interval [m*dz, (m+1)*dz] tiles the column with no gap
                // and no overlap for ANY phi_top. Passing raw phi instead left
                // the slice between the level set and the top solid face inside
                // a void cell (Q = 0): attenuated but never deposited. Before any
                // removal phi = (n + 1/2)*dz, so path == phi to round-off.
                // n_above is taken from the mask (k_top - k), which equals
                // floor(phi/dz) except at round-off ties phi_top ~ 0 / ~ dz,
                // where the floor picks the wrong cell (see SolidFacePath).
                Set::Scalar Q = self_void ? Set::Scalar(0.0)
                                          : BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                                       tv.R(i, j, bR), bprof, brad, b_alpha, x, y, z,
                                                       DX[2],
                                                       // ktc is null unless void-aware (Q = 0 when P_now == 0).
                                                       // D2q-1: `sab` is the solid-cell count above k;
                                                       // it is null when the beam is off, where the
                                                       // path is unused and k_top - k is kept.
                                                       (void_aware && ktc)
                                                         ? SolidFacePath(
                                                             sab ? sab[static_cast<std::size_t>(col_c) * nz_dom
                                                                       + (k - klo_d)]
                                                                 : ktc[col_c] - k, DX[2])
                                                         : Set::Scalar(-1.0),
                                                       tv.t(i, j, k), tv.b(i, j, k), bn);

                const Set::Scalar q_loss = q_out;

                // Conservative H update; T is recovered in InvertHtoT
                // using the per-phase Numeric::Material::Table. Detached
                // (void) cells exchange no heat and get no source, so their
                // H stays frozen at its value when voided (removal never
                // writes H).
                H(i, j, k) = H_old(i, j, k) + dt * (div_kgrad + Q - q_loss * inv_dz);
            });
        }

        // E1: standing beam-deposition closure counter (thermo.dat).
        if (P_now > 0.0) BeamEnergyBalance(lev, time, dt, k_top_col, btau);
        // S1: whole-domain energy ledger (thermo.dat), after the H update.
        // P2-8: the per-class robin power, the kernel's terms summed by class.
        if (fl_on) FaceListDiagnostics(lev);
        if (ledger_enabled) EnergyLedger(lev, time, dt, k_top_col);
    }

    // D2a: per-column flame for this AdvanceMicrostructure call.
    //   constants: flame_h = h_conv, flame_Tg = T_flame (bit-identical path).
    //   expressions: at the column centre (x, y), r = |(x, y) - (x0, y0)|,
    //     z_face = PLO_z + (k_top + 1)*dz (mask top face), z_n = FeetDescent
    //     (D2b, nozzle_descent = feet) or nozzle_z0 -
    //     nozzle_feed*time, s = z_n - z_face; h = h_expr(x,y,r,s,t),
    //     T_gas = T_flame_expr(...). Columns with no solid cell get 0.
    // Aborts if a live column with r <= nozzle_collision_radius (default: the
    // patch radius) has s <= 0 (nozzle in the rock), or if an in-patch column
    // gets a non-finite / non-positive h or T_gas. patch_min_standoff is the
    // min s over live columns with r <= nozzle_collision_radius (in the patch).
    // k_top_col is the reduced
    // ColumnTopSolid, so every rank builds identical vectors.
    void BuildFlameColumns(int lev, Set::Scalar time, Set::Scalar dt,
                           const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        if (!surface_patch.use_expr)
        {
            flame_h.assign(ncols, surface_patch.h_conv);
            flame_Tg.assign(ncols, surface_patch.T_flame);
            flame_s.assign(ncols, 0.0);
            // P2-28: a PRESCRIBED patch on the face list (the only constant-flame list: a
            // convective list needs the h_expr map, Parse) -- the list's half, in the
            // flameless branch's order: every face's h and T_presc, every cell's T_s,
            // every face's gas temperature
            if (surface_patch.face_list)
            {
                FaceListHeat(lev, 0.0, time, k_top_col);
                FrAssignTs(lev);
                FaceListTg(lev);
            }
            return;
        }
        flame_h.assign(ncols, 0.0);
        flame_Tg.assign(ncols, 0.0);
        flame_s.assign(ncols, 0.0);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;
        const Set::Scalar rc2   = surface_patch.nozzle_collision_radius *
                                  surface_patch.nozzle_collision_radius;
        const Set::Scalar z_n = surface_patch.feet
                              ? FeetDescent(lev, dt, k_top_col)
                              : surface_patch.nozzle_z0 - surface_patch.nozzle_feed * time;
        Set::Scalar s_min = std::numeric_limits<Set::Scalar>::infinity();
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t col = static_cast<std::size_t>(j) * nx + i;
            const int k_top = k_top_col[col];
            if (k_top == std::numeric_limits<int>::min()) continue;
            const Set::Scalar x = PLO[0] + (Set::Scalar(i + domain.smallEnd(0)) + 0.5) * DX[0];
            const Set::Scalar y = PLO[1] + (Set::Scalar(j + domain.smallEnd(1)) + 0.5) * DX[1];
            const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
            const Set::Scalar r2 = dxp * dxp + dyp * dyp;
            const Set::Scalar r = std::sqrt(r2);
            const Set::Scalar z_face = PLO[2] + Set::Scalar(k_top + 1) * DX[2];
            const Set::Scalar s = z_n - z_face;
            flame_s[col] = s;
            if (r2 > r2max) continue;
            const bool in_coll = r2 <= rc2;
            if (in_coll && !(s > 0.0))
                Util::Abort(INFO, "surface_patch nozzle collision: stand-off s = ", s,
                            " <= 0 at column (i, j) = (", i + domain.smallEnd(0), ", ",
                            j + domain.smallEnd(1), "), time = ", time,
                            " (z_nozzle = ", z_n, ", z_face = ", z_face, ")");
            if (surface_patch.jet_enthalpy)
            {
                // D2a2: no jet flux above the nozzle plane (s <= 0); T_gas is
                // set by JetEnthalpyMarch (T_ent here, a harmless finite value;
                // above the plane the march never visits the column, so it is
                // final there). P2-26 (census B3): the D2v above-plane pair and
                // the D2aa pad shield that reached this site are deleted.
                const Set::Scalar hj = (s > 0.0)
                                     ? surface_patch.h_f(x, y, r, s, time)
                                     : Set::Scalar(0.0);
                if (!(std::isfinite(hj) && hj >= 0.0))
                    Util::Abort(INFO, "surface_patch.h_expr must be finite and >= 0 in the "
                                "patch; got h = ", hj, " at column (", i + domain.smallEnd(0),
                                ", ", j + domain.smallEnd(1), "), r = ", r, ", s = ", s,
                                ", time = ", time);
                flame_h[col]  = hj;
                flame_Tg[col] = surface_patch.jet_T_ent;
                if (in_coll) s_min = std::min(s_min, s);
                continue;
            }
            const Set::Scalar h  = surface_patch.h_f(x, y, r, s, time);
            const Set::Scalar Tg = surface_patch.Tg_f(x, y, r, s, time);
            if (!(std::isfinite(h) && h > 0.0 && std::isfinite(Tg) && Tg > 0.0))
                Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr must be finite and > 0 "
                            "in the patch; got h = ", h, ", T_gas = ", Tg, " at column (",
                            i + domain.smallEnd(0), ", ", j + domain.smallEnd(1), "), r = ", r,
                            ", s = ", s, ", time = ", time);
            flame_h[col]  = h;
            flame_Tg[col] = Tg;
            if (in_coll) s_min = std::min(s_min, s);
        }
        patch_min_standoff = std::isfinite(s_min) ? s_min : 0.0;
        // P2-8 F2: every listed face's h, at its centre, BEFORE the march.
        if (surface_patch.face_list)
        {
            FaceListHeat(lev, z_n, time, k_top_col);
            // P2-9 R1: every listed cell's T_s, before the march bills them
            FrAssignTs(lev);
        }
        if (surface_patch.jet_enthalpy) JetEnthalpyMarch(lev, time, k_top_col);
        if (surface_patch.face_list) FaceListTg(lev);
    }
