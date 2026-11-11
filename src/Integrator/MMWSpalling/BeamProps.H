    // E1: optical path from the top face of the first solid cell to the
    // centre of a cell with n_above solid cells above it in its column (see
    // the BeamSource call in AdvanceMicrostructure). n_above is COUNTED from
    // the `removed` mask, which is exact by construction. It is NOT evaluated
    // from phi as floor(phi/dz): phi is not the void test (a top cell more
    // than half consumed is solid with phi < 0), and even on a column where
    // the two agree a round-off tie breaks the floor -- e.g. phi = -5e-18 (voided)
    // above phi = dz exactly (top solid) gives 1, and the whole top slice (98%
    // of the column power at alpha*dz = 4) is lost. Do NOT use
    // floor((phi - dz/2)/dz) either: for phi_top in [0, dz/2) the top and
    // second cells both get path dz/2 and the surface slice is deposited twice.
    AMREX_GPU_HOST_DEVICE
    static Set::Scalar SolidFacePath(int n_above, Set::Scalar dz)
    {
        return (Set::Scalar(n_above) + 0.5) * dz;
    }

    // AMR: the finest level owns the surface. Every surface closure (patch,
    // losses, beam, Sp, removal, the per-column state and the CSV logs) runs
    // on this level only; coarser levels advance conduction in their covered
    // cells and are overwritten by average-down. Single level: always true.
    bool SurfaceOwner(int lev) const { return max_level == 0 || lev == finest_level; }

    // =====================================================================
    // P2-2 change 1: UpdatePropertyFields — the per-step property refresh.
    // =====================================================================
    // ⚠⚠ THE PROPERTIES ARE PER-CELL FIELDS WITH TWO WRITERS.
    // `kappa_mf` and `beta_mf` are filled from the phase scalars at init
    // (`InitializeMicrostructure`) and written nowhere else but here, so
    // temperature dependence is NOT a parameter swap -- without this refresh
    // every reader sees the rock the fields were initialised with, whatever
    // the laws say.
    //
    // Called FIRST in Advance, so the criterion and the kernel/march/ledger
    // all see ONE set by construction -- the ordering P2-1 validated.
    //
    // ⚠ `beta_mf` carries the SECANT expansion coefficient, strain/(T - T_ref),
    // NOT alpha(T). Every reader of the field forms beta * (T - T_ref)
    // (ThermalFreeStrain: the eigenstrain of UpdateMicrostructureFreeStrain
    // and of UpdateModel), so writing the secant here makes the product the
    // accumulated strain identically. ⚠⚠ Writing alpha(T) here instead would
    // give every reader alpha(T)*dT. That is gate (b6).
    //
    // No-op, and no traversal, unless some phase sets `law = grimsel` or (P2-29) `zhang`.
    // P2-27 stage 4: also called from InitializeMicrostructure with
    // `step = false` -- a law phase's scalars are refused, so the t = 0 fields
    // come from THIS writer (the law at the start temperature) instead of the
    // scalars; the per-step counters (prop_refresh_n, prop_*_mean) are not
    // touched then, so no thermo column moves.
    void UpdatePropertyFields(int lev, bool step = true)
    {
        bool any = false;
        for (const auto& ph : phases) any = any || ph.grimsel_law || ph.zhang_law;
        if (!any) return;
        if (lev >= static_cast<int>(phase_mf.size()) || !phase_mf[lev]) return;

        const int n_phases = static_cast<int>(phases.size());
        std::vector<char> use(n_phases, 0);
        // P2-29: 2 = a zhang phase, refreshed in its THERMAL fields only (the mechanics are its scalars)
        for (int p = 0; p < n_phases; ++p) use[p] = phases[p].grimsel_law ? 1 : (phases[p].zhang_law ? 2 : 0);
        Set::Scalar acc_k = 0.0, acc_b = 0.0, acc_E = 0.0, acc_T = 0.0, acc_n = 0.0;

        for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T  = (*temp_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& Trf =
                (*T_ref_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& kapa = (*kappa_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& beta = (*beta_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& EE   = (*E_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& muu  = (*mu_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rho_a = (*rho_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Cp_a  = (*Cp_mf[lev]).array(mfi);
            // ⚠ the phase objects are not device-copyable, so this loop is on
            // the host. It is one pass over the level per step.
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                if (!use[p]) continue;
                if (use[p] == 2)
                {
                    // P2-29: Zhang -- kappa_phase is a CONDUCTIVITY [W/(m K)] (Zhang::conductivity)
                    const Numeric::Material::Zhang& z = phases[p].zhang;
                    const Set::Scalar Tz = T(i, j, k);
                    if (!(Tz > 0.0) || !std::isfinite(Tz)) continue;
                    kapa (i, j, k) = z.conductivity(Tz);
                    rho_a(i, j, k) = z.rho(Tz);
                    Cp_a (i, j, k) = z.cp(Tz);
                    continue;
                }
                const Numeric::Material::Grimsel& g = phases[p].grimsel;
                Set::Scalar Tc = T(i, j, k);
                if (!(Tc > 0.0) || !std::isfinite(Tc)) continue;
                const Set::Scalar T_r = Trf(i, j, k);
                kapa (i, j, k) = g.lambda(Tc);
                EE   (i, j, k) = g.E_of(Tc);
                muu  (i, j, k) = g.mu_of(Tc);
                rho_a(i, j, k) = g.rho(Tc);
                Cp_a (i, j, k) = g.cp(Tc);
                // THE SECANT COEFFICIENT -- see the warning above. ⚠ strain_of holds its
                // last node above 1000 degC (1273.15 K), so above it this beta is the held
                // strain over a growing dT (MODEL_CARD §11 item 11).
                const Set::Scalar dT = Tc - T_r;
                beta (i, j, k) = (std::abs(dT) > 1.0e-9)
                    ? (g.strain_of(Tc) - g.strain_of(T_r)) / dT
                    : g.alpha_of(Tc);
                acc_k += kapa(i, j, k); acc_b += beta(i, j, k);
                acc_E += EE(i, j, k);   acc_T += Tc; acc_n += 1.0;
            }
        }
        Util::RealFillBoundary(*kappa_mf[lev],       geom[lev]);
        Util::RealFillBoundary(*beta_mf[lev],        geom[lev]);
        Util::RealFillBoundary(*E_phase_mf[lev],     geom[lev]);
        Util::RealFillBoundary(*mu_phase_mf[lev],    geom[lev]);
        // the diagnostics the (b5) two-temperature probe reads. ⚠ per BOX, so
        // they MUST be reduced -- the P1g-2c lesson.
        if (!step) return;
        amrex::ParallelDescriptor::ReduceRealSum(acc_k);
        amrex::ParallelDescriptor::ReduceRealSum(acc_b);
        amrex::ParallelDescriptor::ReduceRealSum(acc_E);
        amrex::ParallelDescriptor::ReduceRealSum(acc_T);
        amrex::ParallelDescriptor::ReduceRealSum(acc_n);
        if (acc_n > 0.0)
        {
            prop_kappa_mean = acc_k / acc_n;
            prop_beta_mean  = acc_b / acc_n;
            prop_E_mean     = acc_E / acc_n;
            prop_T_mean     = acc_T / acc_n;
        }
        prop_refresh_n += 1.0;
    }

    // E1: global topmost solid (removed <= 0.5) cell index per (i, j)
    // column at level lev, INT_MIN for fully removed columns. Index is
    // (j - dlo_y)*nx + (i - dlo_x). MPI-reduced, so identical on every rank
    // and independent of the z box decomposition. D2q-1: the rule is `removed`
    // alone -- there is no phi >= 0 conjunct.
    amrex::Vector<int> ColumnTopSolid(int lev)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int nx = domain.length(0), ny = domain.length(1);
        amrex::Vector<int> k_top(static_cast<std::size_t>(nx) * ny,
                                 std::numeric_limits<int>::min());
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int& kt = k_top[(j - jlo) * nx + (i - ilo)];
                for (int k = hi.z; k >= lo.z && k > kt; --k)
                    if (rem(i, j, k) <= 0.5) { kt = k; break; }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(k_top.data(), static_cast<int>(k_top.size()));
        return k_top;
    }

    // D2q-1: number of SOLID cells strictly above cell k, up to and including
    // the column top, indexed [col*nz + (k - klo)] with col as in
    // ColumnTopSolid. The beam's optical path must be measured through
    // MATERIAL, not through geometric distance: with an interior void the beam
    // crosses it unattenuated, so `k_top - k` would overstate the path.
    //
    // While void is contiguous every cell in (k, k_top] is solid and this
    // returns exactly `k_top - k`, which is what keeps the beam runs
    // byte-identical.
    //
    // Cost: one ncols*nz int reduction. Built ONLY when the beam is on (the
    // path argument is unused otherwise -- BeamSource returns 0 at P = 0 and
    // BeamEnergyBalance skips columns with Pi == 0), so the Meier and arbiter
    // runs, which are the large ones, never pay for it.
    amrex::Vector<int> ColumnSolidAboveCount(int lev, const amrex::Vector<int>& k_top)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int klo = domain.smallEnd(2);
        const int nx = domain.length(0), ny = domain.length(1), nz = domain.length(2);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        // Each cell is owned by exactly one box at a level, so a sum over ranks
        // gathers the column without double counting.
        amrex::Vector<int> solid(ncols * static_cast<std::size_t>(nz), 0);
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
                if (rem(i, j, k) <= 0.5)
                    solid[static_cast<std::size_t>((j - jlo) * nx + (i - ilo)) * nz
                          + (k - klo)] = 1;
        }
        amrex::ParallelDescriptor::ReduceIntSum(solid.data(), static_cast<int>(solid.size()));

        amrex::Vector<int> cnt(ncols * static_cast<std::size_t>(nz), 0);
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = k_top[c];
            if (kt == std::numeric_limits<int>::min()) continue;
            int acc = 0;
            for (int k = kt; k >= klo; --k)
            {
                const std::size_t p = c * nz + (k - klo);
                cnt[p] = acc;                 // solid cells in (k, k_top]
                if (solid[p]) ++acc;
            }
        }
        return cnt;
    }

    // P2-21: THE COLUMN OPTICAL DEPTH for a temperature-dependent loss
    // (beam.tan_delta_T). Mirrors ColumnSolidAboveCount: every solid cell's
    // a_k*dz, a_k = beam.alpha_at(T_old), into an ncols*nz vector (-1 marks a
    // void cell; each cell is owned by exactly one box, so a sum over ranks
    // gathers the column), ONE ReduceRealSum, then per column from the top solid
    // cell down the running sums tau_top (0 at the top cell) and tau_bot =
    // tau_top + a_k*dz -- the SAME a_k*dz, never recomputed from another T.
    // The CENTRAL column's temperatures (the column containing (x0, y0); it can
    // span several z-boxes, so it is never assumed to sit on one rank) ride in
    // the same vector, so the deciders read the T_old the gather used.
    // k_top_col empty => every column's top is the domain top and every cell is
    // solid (a deck without the face list: nothing is ever removed);
    // the beam then enters at the domain top, so beam.z0 must be that top.
    // Also sets the thermo deciders beam_tau_max ... bc_ell_pk.
    // Cost: one (ncols*nz + nz + ncols)-double allreduce per step, ONLY on a deck whose
    // table varies and only while the beam is on (BuildBeamTau).
    BeamTau ColumnOpticalDepth(int lev, const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int klo = domain.smallEnd(2), khi = domain.bigEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t NZ = static_cast<std::size_t>(domain.length(2));
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz = DX[2];
        const bool have_tops = !k_top_col.empty();
        const int kmin = std::numeric_limits<int>::min();
        if (!have_tops && std::abs(beam.z0 - geom[lev].ProbHi()[2]) > 1.0e-9 * dz)
            Util::Abort(INFO, "beam.tan_delta_T: the column optical depth starts at the domain top (",
                        geom[lev].ProbHi()[2], " m); beam.z0 = ", beam.z0, " must equal it");
        auto top_of = [&](std::size_t c) { return have_tops ? k_top_col[c] : khi; };

        // a column is LIT when the beam reaches its top cell (Pi > 0 for any P > 0;
        // a Gaussian's waist only grows with depth, a uniform disk's radius is fixed)
        std::vector<char> lit(ncols, 0);
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = top_of(c);
            if (kt == kmin) continue;
            const int i = ilo + static_cast<int>(c % nx), j = jlo + static_cast<int>(c / nx);
            lit[c] = BeamIncidentIntensity(1.0, beam.omega0, beam.lambda, beam.x0, beam.y0, beam.z0,
                                           beam.profile, beam.radius,
                                           PLO[0] + (Set::Scalar(i) + 0.5) * DX[0],
                                           PLO[1] + (Set::Scalar(j) + 0.5) * DX[1],
                                           PLO[2] + (Set::Scalar(kt) + 0.5) * dz) > 0.0;
        }
        const int ic = std::min(nx - 1, std::max(0, static_cast<int>(std::floor((beam.x0 - PLO[0]) / DX[0]))));
        const int jc = std::min(ny - 1, std::max(0, static_cast<int>(std::floor((beam.y0 - PLO[1]) / DX[1]))));
        const std::size_t c0 = static_cast<std::size_t>(jc) * nx + ic;

        // [ncols*NZ: a_k dz | NZ: the central column's T | ncols: each column's TOP solid cell T (P2-24)]
        std::vector<Set::Scalar> g(ncols * NZ + NZ + ncols, 0.0);
        const std::size_t o_top = ncols * NZ + NZ;
        const Set::Scalar mf = beam.melt_factor;
        Set::Scalar amax = 0.0;
        for (amrex::MFIter mfi(*temp_old_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& rem =
                have_tops ? (*removed_mf[lev]).const_array(mfi) : (*temp_old_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(j - jlo) * nx + (i - ilo);
                const std::size_t p = c * NZ + (k - klo);
                if (c == c0) g[ncols * NZ + (k - klo)] = T(i, j, k);
                if (have_tops && rem(i, j, k) > 0.5) { g[p] = -1.0; continue; }
                if (k == top_of(c)) g[o_top + c] = T(i, j, k);
                // P2-24: the melt's absorption on the liquid fraction (1 + (m - 1) f_L); m = 1 skips it
                const Set::Scalar a = (mf != 1.0)
                    ? beam.alpha_at(T(i, j, k)) * (1.0 + (mf - 1.0) * MeltFraction(T(i, j, k)))
                    : beam.alpha_at(T(i, j, k));
                g[p] = a * dz;
                if (lit[c] && a > amax) amax = a;
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(g.data(), static_cast<int>(g.size()));
        amrex::ParallelDescriptor::ReduceRealMax(amax);

        BeamTau out;
        out.top.assign(ncols * NZ, -1.0);
        out.bot.assign(ncols * NZ, -1.0);
        out.k_last.assign(ncols, kmin);
        Set::Scalar tmax = 0.0, emax = 0.0;
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = top_of(c);
            if (kt == kmin) continue;
            Set::Scalar acc = 0.0, e_top = 1.0, s = 0.0;
            for (int k = kt; k >= klo; --k)
            {
                const std::size_t p = c * NZ + (k - klo);
                if (g[p] < 0.0) continue;              // void: no material, no attenuation
                out.top[p] = acc;
                out.bot[p] = acc + g[p];
                acc = out.bot[p];
                out.k_last[c] = k;
                // the kernel's own pair: e^{-tau_top} of this cell IS e^{-tau_bot}
                // of the solid cell above (the same double), so the sum telescopes
                const Set::Scalar e_bot = std::exp(-out.bot[p]);
                s += e_top - e_bot;
                e_top = e_bot;
            }
            if (!lit[c] || out.k_last[c] == kmin) continue;
            tmax = std::max(tmax, acc);
            emax = std::max(emax, std::abs(s + e_top - 1.0));
        }
        beam_tau_max     = tmax;
        beam_alpha_max   = amax;
        beam_col_err_max = emax;

        // P2-24: the column's reflection from its top solid cell's liquid fraction (only when the melt keys
        // differ from the no-op: unbuilt => every site reads the scalar R, byte for byte), and the surface
        // record of the melt: the synthetic radiometer over the DOE launcher's view, the melt spot's area
        const Set::Scalar R_s = beam.R, R_l = beam.R_melt;
        if (beam.melt_active())
        {
            out.Rcol.assign(ncols, R_s);
            for (std::size_t c = 0; c < ncols; ++c)
                if (top_of(c) != kmin) out.Rcol[c] = R_s + (R_l - R_s) * MeltFraction(g[o_top + c]);
        }
        {
            // the radiometer's view: the 32.5 mm launcher aperture (DOE DE-EE0005504 L992-1002), r <= 16.25 mm
            // from the beam axis. eps_137 = (1 - R_melt) f_L + (1 - R)(1 - f_L): the DOE's own melt absorptance
            // (0.7) and the solid's from its permittivity, both from the one dielectric description. A LOG.
            const Set::Scalar r_view = 0.01625;
            Set::Scalar n_v = 0.0, sTr = 0.0, sT = 0.0, A_m = 0.0;
            for (std::size_t c = 0; c < ncols; ++c)
            {
                if (top_of(c) == kmin) continue;
                const Set::Scalar Tt = g[o_top + c];
                const Set::Scalar fL = MeltFraction(Tt);
                if (fL > 0.0) A_m += DX[0] * DX[1];
                const int i = ilo + static_cast<int>(c % nx), j = jlo + static_cast<int>(c / nx);
                const Set::Scalar dxc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0] - beam.x0;
                const Set::Scalar dyc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1] - beam.y0;
                if (dxc * dxc + dyc * dyc > r_view * r_view) continue;
                n_v += 1.0;
                sT  += Tt;
                sTr += ((1.0 - R_l) * fL + (1.0 - R_s) * (1.0 - fL)) * Tt;
            }
            bc_T_rad  = (n_v > 0.0) ? sTr / n_v : std::numeric_limits<Set::Scalar>::quiet_NaN();
            bc_T_view = (n_v > 0.0) ? sT / n_v : std::numeric_limits<Set::Scalar>::quiet_NaN();
            bc_A_melt = A_m;
        }

        // the central column: its contiguous solid run from the top, T from the
        // gathered T_old; depths are cell centres below the top solid face
        const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
        bc_T_peak = bc_z_peak = bc_z_melt = bc_tau = bc_ell_p = bc_ell_pk = nan;
        bc_f_L_peak = bc_z_liq = bc_R_col = nan;
        const int kt0 = top_of(c0);
        if (kt0 != kmin)
        {
            // P2-24: the central column's reflection, from its top solid cell (the form every site reads)
            bc_R_col = R_s + (R_l - R_s) * MeltFraction(g[o_top + c0]);
            // P2-24: the ROCK PHASE's solidus (phase 0's T_m; 1e30 = no band when there is no phase)
            const Set::Scalar Tm = !phases.empty() ? phases[0].T_m : Set::Scalar(1.0e30);
            std::vector<Set::Scalar> Tcol;
            for (int k = kt0; k >= klo && g[c0 * NZ + (k - klo)] >= 0.0; --k)
                Tcol.push_back(g[ncols * NZ + (k - klo)]);
            const int nmax = static_cast<int>(Tcol.size()) - 1;
            int d_pk = 0;
            for (int d = 0; d <= nmax; ++d)
            {
                if (Tcol[d] > Tcol[d_pk]) d_pk = d;
                if (!std::isfinite(bc_z_melt) && Tcol[d] >= Tm) bc_z_melt = (Set::Scalar(d) + 0.5) * dz;
            }
            bc_T_peak = Tcol[d_pk];
            bc_z_peak = (Set::Scalar(d_pk) + 0.5) * dz;
            // P2-24: the liquid fraction at the peak, and the deepest fully liquid cell of the column
            bc_f_L_peak = MeltFraction(bc_T_peak);
            for (int d = 0; d <= nmax; ++d)
                if (MeltFraction(Tcol[d]) >= 1.0) bc_z_liq = (Set::Scalar(d) + 0.5) * dz;
            if (out.k_last[c0] != kmin) bc_tau = out.bot[c0 * NZ + (out.k_last[c0] - klo)];
            // P2-19's helpers on the rise above the ambient (the Zhang decks' IC);
            // the face sample is the top cell's rise (no surface reconstruction here)
            auto rise_at_d = [&](int d) { return Tcol[d] - T_amb; };
            Set::Scalar z_pk = nan;
            bc_ell_p = ProfileEFoldDepth(rise_at_d, Tcol[0] - T_amb, dz, nmax);
            ProfilePeak(rise_at_d, Tcol[0] - T_amb, dz, nmax, z_pk, bc_ell_pk);
        }
        out.built = true;
        return out;
    }

    // P2-21: the column optical depth when it applies (the table varies AND the
    // beam is on this step), else an unbuilt BeamTau (every kernel then keeps the
    // closed form) with the deciders cleared on a table deck.
    BeamTau BuildBeamTau(int lev, Set::Scalar P_now, const amrex::Vector<int>& k_top_col)
    {
        beam_P_beam = P_now;   // P2-24: the schedule's power now (the DOE's denominator)
        if (beam.alpha_varies() && P_now > 0.0) return ColumnOpticalDepth(lev, k_top_col);
        if (beam.alpha_varies())
        {
            const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
            beam_tau_max = beam_alpha_max = beam_col_err_max = 0.0;
            bc_T_peak = bc_z_peak = bc_z_melt = bc_tau = bc_ell_p = bc_ell_pk = nan;
            // P2-24: with the beam OFF the surface still cools; on a BEAM deck the radiometer and the melt spot
            // are read anyway (a top-cell gather, no optical depth). A flame deck (no beam power ever) skips it.
            if (beam.can_deposit()) SurfaceMeltRecord(lev, k_top_col);
        }
        return BeamTau();
    }

    // P2-24: THE ROCK'S LIQUID FRACTION f_L(T) -- ONE phase (SR 13): phase 0's fusion band on the microstructure
    // path (P2-29: the only path); 0 when no band is declared (T_m >= 1e29).
    // f_L = clamp((T - (T_m - w/2)) / w, 0, 1); a sharp band (w <= 0) is the step at T_m. The beam's melt factor
    // and reflection, the face list's molten / lid exclusions and the melt records all read THIS.
    Set::Scalar MeltFraction(Set::Scalar T) const
    {
        Set::Scalar Tm = 1.0e30, w = 0.0;
        if (!phases.empty()) { Tm = phases[0].T_m; w = phases[0].fus_width; }
        if (!(Tm < 1.0e29)) return 0.0;
        if (!(w > 0.0)) return (T >= Tm) ? 1.0 : 0.0;
        const Set::Scalar f = (T - (Tm - 0.5 * w)) / w;
        return (f <= 0.0) ? 0.0 : ((f >= 1.0) ? 1.0 : f);
    }

    // P2-24: the beam-OFF surface record (after the beam is switched off the melt still cools): bc_T_rad,
    // bc_T_view, bc_A_melt from each column's top solid cell, as ColumnOpticalDepth forms them
    void SurfaceMeltRecord(int lev, const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1), khi = domain.bigEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const bool have_tops = !k_top_col.empty();
        const int kmin = std::numeric_limits<int>::min();
        std::vector<Set::Scalar> Ttop(ncols, 0.0);
        for (amrex::MFIter mfi(*temp_old_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T = (*temp_old_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(j - jlo) * nx + (i - ilo);
                const int kt = have_tops ? k_top_col[c] : khi;
                if (kt != kmin && kt >= lo.z && kt <= hi.z) Ttop[c] = T(i, j, kt);
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(Ttop.data(), static_cast<int>(Ttop.size()));
        const Set::Scalar R_s = beam.R, R_l = beam.R_melt, r_view = 0.01625;   // as ColumnOpticalDepth
        Set::Scalar n_v = 0.0, sTr = 0.0, sT = 0.0, A_m = 0.0;
        for (std::size_t c = 0; c < ncols; ++c)
        {
            if ((have_tops ? k_top_col[c] : khi) == kmin) continue;
            const Set::Scalar fL = MeltFraction(Ttop[c]);
            if (fL > 0.0) A_m += DX[0] * DX[1];
            const int i = ilo + static_cast<int>(c % nx), j = jlo + static_cast<int>(c / nx);
            const Set::Scalar dxc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0] - beam.x0;
            const Set::Scalar dyc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1] - beam.y0;
            if (dxc * dxc + dyc * dyc > r_view * r_view) continue;
            n_v += 1.0;
            sT  += Ttop[c];
            sTr += ((1.0 - R_l) * fL + (1.0 - R_s) * (1.0 - fL)) * Ttop[c];
        }
        const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
        bc_T_rad  = (n_v > 0.0) ? sTr / n_v : nan;
        bc_T_view = (n_v > 0.0) ? sT / n_v : nan;
        bc_A_melt = A_m;
        bc_f_L_peak = bc_z_liq = bc_R_col = nan;
    }

    // E1: per-step closure of the microstructure-path beam deposition.
    //
    //   P_inc   = sum over lit columns of (1-R)*Pi*dx*dy at the top solid cell
    //   P_dep   = sum over solid cells of Q*dx*dy*dz  (same BeamSource call)
    //   P_trans = (1-R)*Pi*exp(-a*s_bot)*dx*dy at k == domain bottom
    //   closure : P_dep + P_trans == P_inc
    //
    // Exact (to round-off) only for CONSTANT alpha and the UNIFORM profile.
    // A temperature-dependent alpha deposits through the cumulative optical
    // depth tau = sum_k a_k*dz (P2-21: ColumnOpticalDepth), never a per-cell
    // a_k applied to the whole path above (which over-deposits). A Gaussian profile adds an
    // O((dz/z_R)^2) deviation because Pi evolves with the geometric depth, so
    // the P_dep <= P_inc double-count guard is applied to the uniform profile
    // only.
    //
    // Invariant (void-aware runs), per solid cell, where n = the count of SOLID
    // cells in (k, k_top] (D2q-1: a solid count, NOT k_top - k):
    //   (a) n == 0 <=> k == k_top.
    //       Exact at any topology because k_top is itself solid, so the count
    //       vanishes only at the top cell. Do NOT test "cell above is void or
    //       k is the domain top <=> n == 0" instead: that is a CONTIGUITY
    //       test, which an interior void breaks (the cell under the void
    //       satisfies the left side with solid still above it).
    // The P_inc tally below is gated on the SAME k == k_top test, for the same
    // reason: a contiguity test counts incident power twice per column once
    // an interior void exists (measured 14.29 % closure error).
    // Runs on the surface-owning level only (SurfaceOwner: level 0 on a
    // single-level run). Multi-level is unsupported: the sums cover that
    // level's own cells, not a composite.
    // P2-21: with a built BeamTau (a varying loss table) every cell's source is
    // the optical-depth pair the kernel used and P_trans leaves at each column's
    // LOWEST solid cell with e^{-tau_bot} there; the per-column pairs telescope
    // (beam_col_err_max), but the closure keeps a residual of the size
    // (d/z_R)^2 x the centre-to-edge contrast in alpha: Pi is the cell's own, and
    // a Gaussian's Pi changes with depth -- the lateral power a column sum cannot
    // follow (studies/p2_21_loss_table/step0). Unbuilt: the closed form.
    void BeamEnergyBalance(int lev, Set::Scalar time, Set::Scalar dt,
                           const amrex::Vector<int>& k_top_col,
                           const BeamTau& btau)
    {
        if (!SurfaceOwner(lev)) return;
        // P2-28: the kernel's own switch (list_state)
        const bool void_aware = list_state();
        const amrex::Box& domain = geom[lev].Domain();
        const int klo_dom = domain.smallEnd(2), khi_dom = domain.bigEnd(2);
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const int ncx = domain.length(0);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz = DX[2];
        const Set::Scalar A  = DX[0] * DX[1];
        const Set::Scalar V  = A * dz;
        const Set::Scalar P_now = beam.power_at(time);
        const Set::Scalar b_alpha = beam.alpha;   // P2-20: computed from beam.eps_r, beam.tan_delta, beam.lambda
        // D2q-1: solid-cell count above k, as in AdvanceMicrostructure.
        const amrex::Vector<int> solid_above =
            (void_aware && P_now > 0.0 && !k_top_col.empty())
            ? ColumnSolidAboveCount(lev, k_top_col) : amrex::Vector<int>();
        const int* sab    = solid_above.data();
        const int  nz_dom = domain.length(2);
        const TauView tv = ViewOf(btau, lev);

        Set::Scalar P_inc = 0.0, P_dep = 0.0, P_trans = 0.0, n_viol = 0.0;
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            // Placeholders as in AdvanceMicrostructure when not void-aware.
            // (No temperature is read here: the closed form's alpha is the
            // constant beam.alpha (P2-20), and a varying alpha arrives as the
            // BeamTau pairs the kernel used.)
            amrex::Array4<const Set::Scalar> const& rem =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (void_aware && rem(i, j, k) > 0.5) continue;
                const int col_b = (j - jlo_dom) * ncx + (i - ilo_dom);
                const int n_above = void_aware
                    ? (sab ? sab[static_cast<std::size_t>(col_b) * nz_dom + (k - klo_dom)]
                           : k_top_col[col_b] - k)
                    : 0;
                const Set::Scalar path = void_aware ? SolidFacePath(n_above, dz)
                                                    : Set::Scalar(-1.0);
                // D2q-2b: the column top, for BOTH invariant (a) and the P_inc
                // tally below. NOT "the cell directly above is void, or this
                // is the domain top": that is a contiguity test and counts
                // MORE THAN ONE cell per column the moment an interior void
                // exists -- the true top, and again the cell sitting just
                // under the void -- which double-counts the incident beam
                // power (measured: a 14.29 % closure error against a predicted
                // 140/120 over-count on the seeded probe). The two tests agree
                // while void is contiguous.
                const bool at_col_top = void_aware ? (k == k_top_col[col_b])
                                                   : (k == khi_dom);
                if (void_aware)
                {
                    // (a) (D2q-2b item 0.3): `n_above == 0` iff `k == k_top`.
                    //
                    // Exact at any topology, because n_above counts SOLID
                    // cells in (k, k_top] and k_top is itself solid: below the
                    // top there is always at least one, so the count vanishes
                    // exactly at the top cell. (A comparison against "the cell
                    // directly above is void" would be a CONTIGUITY test: a
                    // cell sitting just under an interior void satisfies it
                    // while still having solid above, and would be counted as
                    // a violation.)
                    if (at_col_top != (n_above == 0)) n_viol += 1.0;
                    // D2q-1: there is NO invariant tying phi to n (phi - n*dz
                    // in [0, dz)): that would assert phi and the `removed`
                    // mask agree, and phi is not a void oracle. n_viol feeds
                    // the `beam_invariant_viol` thermo.dat column.
                }

                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                // P2-22: `path` is -1 off the void-aware path (BeamSource then takes the
                // geometric z0 - z as the solid path); the incident intensity must use the
                // same solid path as the deposit or P_inc and P_dep disagree (smoke: 1.43).
                const Set::Scalar s_rock = (path >= 0.0) ? path : (beam.z0 - z);
                const Set::Scalar Pi = BeamIncidentIntensity(P_now, beam.omega0, beam.lambda,
                                                             beam.x0, beam.y0, beam.z0,
                                                             beam.profile, beam.radius, x, y, z,
                                                             beam.n_idx, s_rock);
                if (Pi == 0.0) continue;

                // P2-24: the column's R (the melt's reflection over a molten top), the deck's scalar R unbuilt
                const Set::Scalar Rc = tv.R(i, j, beam.R);
                const Set::Scalar Q = BeamSource(P_now, beam.omega0, beam.lambda,
                                                 beam.x0, beam.y0, beam.z0, Rc,
                                                 beam.profile, beam.radius, b_alpha,
                                                 x, y, z, dz, path, tv.t(i, j, k), tv.b(i, j, k),
                                                 beam.n_idx);
                P_dep += Q * V;
                if (at_col_top) P_inc += (1.0 - Rc) * Pi * A;
                if (btau.built)
                {
                    // P2-21: the gathered value at the column's exit cell
                    if (k == btau.k_last[static_cast<std::size_t>(col_b)])
                        P_trans += (1.0 - Rc) * Pi * std::exp(-tv.b(i, j, k)) * A;
                }
                else if (k == klo_dom)
                {
                    const Set::Scalar a = b_alpha;
                    const Set::Scalar s_bot = (path >= 0.0 ? path : (beam.z0 - z)) + 0.5 * dz;
                    P_trans += (1.0 - beam.R) * Pi * std::exp(-a * s_bot) * A;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(P_inc);
        amrex::ParallelDescriptor::ReduceRealSum(P_dep);
        amrex::ParallelDescriptor::ReduceRealSum(P_trans);
        amrex::ParallelDescriptor::ReduceRealSum(n_viol);

        // Double-count guard (uniform profile only; see header comment).
        if (beam.profile == Numeric::MMWBeam::Uniform &&
            P_dep > P_inc * (1.0 + 1.0e-9)) n_viol += 1.0;

        beam_P_inc   = P_inc;
        beam_P_dep   = P_dep;
        beam_P_trans = P_trans;
        beam_E_inc   += dt * P_inc;
        beam_E_dep   += dt * P_dep;
        beam_E_trans += dt * P_trans;
        beam_closure_err = (P_inc > 0.0) ? (P_dep + P_trans - P_inc) / P_inc : 0.0;

        if (n_viol > 0.0)
        {
#ifdef AMREX_DEBUG
            Util::Abort(INFO, "E1 beam closure violation at t = ", time,
                        ": invariant/double-count count = ", n_viol,
                        ", P_inc = ", P_inc, ", P_dep = ", P_dep, ", P_trans = ", P_trans);
#endif
            if (beam_invariant_viol == 0.0 || (beam_viol_steps % 1000) == 0)
                Util::Warning(INFO, "E1 beam closure violation at t = ", time,
                              ": invariant/double-count count = ", n_viol,
                              " (cumulative ", beam_invariant_viol + n_viol,
                              "), P_inc = ", P_inc, ", P_dep = ", P_dep,
                              ", P_trans = ", P_trans);
            beam_invariant_viol += n_viol;
            ++beam_viol_steps;
        }
    }

    // Cell-centred Beer-Lambert + Gaussian volumetric source (eqs. 1-3 of
    // main.tex), scaled by (1 − R) for surface reflectivity. Returns 0
    // above the surface (z > z_0).
    //
    // `dz_abs` is the optical path the beam travels through *material* to
    // reach this cell. When removal is active and the surface has receded,
    // the beam crosses the void cavity unattenuated, so the Beer-Lambert
    // attenuation must be measured from the current surface. E1: that surface
    // is the TOP FACE OF THE FIRST SOLID CELL (the `removed` mask), not the
    // level set -- pass dz_abs = SolidFacePath(n_above, dz) = (n_above+1/2)*dz,
    // with n_above = solid cells above in the column, counted from the mask
    // (ColumnSolidAboveCount; see SolidFacePath for why not floor(phi/dz)).
    // Raw phi leaks the slice between the level set and that face.
    // dz_abs < 0 (the default) falls back to the geometric depth z_0 - z,
    // which is correct whenever no material above this cell has been removed
    // (then z_surface == z_0 and phi == z_0 - z exactly). The Gaussian waist
    // evolves with the geometric distance z_0 - z through the void cavity,
    // less the in-rock reduction of the solid path (P2-22, see
    // BeamIncidentIntensity; none at n_rock = 1).
    AMREX_GPU_HOST_DEVICE
    // `dz_cell` (> 0) enables ENERGY-CONSERVING per-cell deposition: the
    // Beer-Lambert profile is integrated analytically over the cell's
    // material-path interval [path - dz/2, path + dz/2] and divided by dz,
    // instead of sampled pointwise at the cell centre. Pointwise sampling
    // under-deposits by a factor ~ (alpha*dz)*exp(-alpha*dz/2)/(1-e^{-alpha*dz})
    // once the absorption depth 1/alpha is smaller than the cell -- at
    // alpha*dz = 4 (Meier surrogate: alpha = 2000/m, dz = 2 mm) only 54% of
    // the beam power reached the rock; at alpha*dz = 40 essentially none.
    // The integrated form reduces to the pointwise one as alpha*dz -> 0.
    // dz_cell <= 0 keeps the legacy pointwise sample.
    static Set::Scalar BeamSource(Set::Scalar P, Set::Scalar w0, Set::Scalar lam,
                                  Set::Scalar x0, Set::Scalar y0, Set::Scalar z0,
                                  Set::Scalar R, int profile, Set::Scalar radius,
                                  Set::Scalar a,
                                  Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                  Set::Scalar dz_cell = 0.0,
                                  Set::Scalar dz_abs = -1.0,
                                  Set::Scalar tau_top = -1.0,
                                  Set::Scalar tau_bot = -1.0,
                                  Set::Scalar n_rock = 1.0)
    {
        if (P == 0.0) return 0.0; // perf: beam off (Pi would be 0 anyway)
        Set::Scalar dz = z0 - z;
        // P2-22: the solid path to this cell (dz when the beam origin is the
        // rock surface) sets the in-rock spread; n_rock 1 keeps the geometric dz.
        Set::Scalar path = (dz_abs >= 0.0) ? dz_abs : dz;
        Set::Scalar Pi = BeamIncidentIntensity(P, w0, lam, x0, y0, z0,
                                               profile, radius, x, y, z, n_rock, path);
        if (Pi == 0.0) return 0.0;
        // P2-21: a VARYING coefficient. The column's optical depth to this cell's
        // top and bottom faces (ColumnOpticalDepth) replaces a*s: the Beer-Lambert
        // integral of a(T) along the solid path, so the column's pairs telescope to
        // 1 - e^{-tau_N}. Pi stays the cell's own (the closed form's, unchanged).
        // tau_top < 0 (the default) is the closed form below, byte for byte.
        if (tau_top >= 0.0)
            return (1.0 - R) * Pi * (std::exp(-tau_top) - std::exp(-tau_bot)) / dz_cell;
        // P2-20: a = beam.alpha, the POWER coefficient computed from the
        // dielectric description (no temperature dependence: user, 2026-11-05)
        if (dz_cell > 0.0)
        {
            // Analytic cell integral of I*a*exp(-a*s) over the material
            // path spanned by this cell, / dz_cell. s_top clamps at 0 so
            // the surface cell absorbs from the true surface even when the
            // (possibly receded) surface lies inside the cell.
            Set::Scalar s_top = path - 0.5 * dz_cell;
            if (s_top < 0.0) s_top = 0.0;
            Set::Scalar s_bot = path + 0.5 * dz_cell;
            return (1.0 - R) * Pi
                   * (std::exp(-a * s_top) - std::exp(-a * s_bot)) / dz_cell;
        }
        return (1.0 - R) * Pi * a * std::exp(-a * path);
    }

    // E1: incident beam intensity Pi [W/m^2] at (x, y), before reflection and
    // Beer-Lambert attenuation. Factored out of BeamSource with the same
    // arithmetic in the same order (BeamSource output is bit-identical).
    // Returns 0 above the beam origin (z > z0), for P == 0, and outside the
    // uniform disk. The Gaussian waist evolves with the reduced distance
    // below (the geometric z0 - z at n_rock = 1).
    AMREX_GPU_HOST_DEVICE
    // P2-22: n_rock = sqrt(eps_r) and s_rock = the solid path from the surface to
    // this point. A Gaussian beam crossing a plane interface into a dielectric
    // keeps its waist and propagates with the Rayleigh range n*zR, i.e. the waist
    // law is evaluated at the REDUCED distance d_air + s_rock / n. Written as
    // dz - s_rock*(1 - 1/n) so that n = 1 reproduces the geometric dz BIT FOR BIT
    // (s_rock * 0.0 = 0.0; every eps_r 1 fixture is byte-identical).
    static Set::Scalar BeamIncidentIntensity(Set::Scalar P, Set::Scalar w0, Set::Scalar lam,
                                             Set::Scalar x0, Set::Scalar y0, Set::Scalar z0,
                                             int profile, Set::Scalar radius,
                                             Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                             Set::Scalar n_rock = 1.0, Set::Scalar s_rock = 0.0)
    {
        Set::Scalar dz = z0 - z;
        if (dz < 0.0 || P == 0.0) return 0.0;
        Set::Scalar dz_red = dz - s_rock * (1.0 - 1.0 / n_rock);
        Set::Scalar r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        Set::Scalar Pi;
        if (profile == Numeric::MMWBeam::Uniform)
        {
            // Flat-top disk: uniform intensity P/(π radius²) inside the disk,
            // zero outside. No Gaussian waist evolution (the beam propagates
            // geometrically through the void cavity as a collimated column).
            if (r2 > radius * radius) return 0.0;
            Pi = P / (M_PI * radius * radius);
        }
        else
        {
            Set::Scalar zR = M_PI * w0 * w0 / lam;
            Set::Scalar w  = w0 * std::sqrt(1.0 + (dz_red / zR) * (dz_red / zR));
            Pi = (2.0 * P / (M_PI * w0 * w0)) * (w0 / w) * (w0 / w)
               * std::exp(-2.0 * r2 / (w * w));
        }
        return Pi;
    }
