#ifndef INTEGRATOR_MMWSPALLING_FACEREMOVAL_H
#define INTEGRATOR_MMWSPALLING_FACEREMOVAL_H

// Fragment of class Integrator::MMWSpalling (P2-9). Textually #included INSIDE
// the class body of src/Integrator/MMWSpalling.H, right after Removal.H. A
// class-body fragment, NOT a standalone header: no #include / namespace /
// class of its own.

    // =====================================================================
    // P2-9 -- ONE BREAKING RULE FOR EVERY SURFACE CELL, READ AND RECEDED
    // ALONG ONE LINE INTO THE ROCK (surface_patch.face_list = 1; no new key).
    // =====================================================================
    // R0  A surface cell is a solid cell with >= 1 face in the list (F1). Its
    //     LINE runs into the rock opposite its exposed face e* with the largest
    //     n.e (n = F3's unit outward normal); mu = n.e*. Ties: the lowest
    //     direction index (the list is sorted by dir). A SHEET (F3 zero
    //     gradient: n undefined) evaluates a line behind EACH exposed face and
    //     takes the largest score, mu = 1. The normal is used for nothing else.
    // R1  ONE surface temperature per surface cell per step: ReconFromX on the
    //     cell and the next cell along the line, cell size dx NOT rescaled
    //     (along a straight line the tilt cancels: the profile is exponential
    //     with length ell/mu, which the two-cell ratio recovers exactly); only
    //     the two NORMAL lengths h_spall_max and flake_h_max become line
    //     lengths by / mu. No inward decay or no solid cell behind => the face
    //     heat balance on e* (the face Newton), always defined, no history.
    // R2  ONE criterion: the depth scan (Removal.H's DepthScanH and its
    //     pieces, the SAME functions the column path calls) along the line,
    //     over the cells gathered behind the surface cell up to
    //     floor(h_spall_max / (mu dx)), stopping at the first non-solid cell
    //     or the domain edge (`run`). Score = K_I/K_Ic at the surface.
    // R3  ONE removal function, RecedeAlongLine(line, h), blind to what
    //     produced h: a per-cell SOLID FRACTION f in (0, 1] drops by h/dx along
    //     the line, a cell reaching f <= eps voids and the remainder passes to
    //     the next cell; the heat is P0b's remap along the line (RemapProfile).
    //     Every line's change is a per-cell delta from the step's OLD field,
    //     summed, then one write (P1g-2c), then the existing InvertHtoT hook.

    // One surface cell (global, identical on every rank).
    struct FrCell
    {
        int i = 0, j = 0, k = 0, col = 0;
        int q0 = 0, nq = 0;        // its faces in fl_faces, [q0, q0 + nq)
        int qline = 0;             // the face e* (non-sheet)
        int sheet = 0;             // 1 = n undefined: a line behind every face
        int misaligned = 0;        // n defined but leans along no exposed face
        Set::Scalar mu = 1.0;
    };
    // One evaluated line (a non-sheet cell has one, a sheet one per face).
    struct FrLine
    {
        int cell = 0;              // index into the FrCell vector
        int dir = 0;               // e* (0 -x .. 5 +z); the line runs to -e*
        int qface = 0;             // fl_faces index of e*
        Set::Scalar mu = 1.0;
        int Dmax = 0;              // floor(h_spall_max / (mu dx))
        int Lg = 0;                // the scan / removal length (Dmax + the shift's reach)
        int Lr = 0;                // gathered length: Lg, or further for the constraint
                                   // reference (P2-14a; = Lg at d = 0)
        std::size_t off = 0;       // into the packed gather
    };
    struct FrCarry { Set::Scalar Ts = 0.0; int fb = 0; };

    // the line classes: floor (+z), wall (+-x, +-y), roof (-z)
    static int FrClass(int dir) { return dir == 5 ? 0 : (dir == 4 ? 2 : 1); }
    static Set::Scalar FrNdotE(const FaceEntry& e, int dir)
    {
        switch (dir)
        {
        case 0: return -e.nx; case 1: return e.nx;
        case 2: return -e.ny; case 3: return e.ny;
        case 4: return -e.nz; default: return e.nz;
        }
    }

    // state and witnesses
    Set::Field<Set::Scalar> fr_f_mf, fr_score_mf, fr_mu_mf;
    bool fr_f_init = false;
    std::map<long long, FrCarry> fr_carry;   // T_s handed from removal to the BC
    bool fr_carry_valid = false;
    Set::Scalar fr_n_c[3] = {0.0, 0.0, 0.0}, fr_v_c[3] = {0.0, 0.0, 0.0};
    Set::Scalar fr_fb_c[3] = {0.0, 0.0, 0.0}, fr_smax_c[3] = {0.0, 0.0, 0.0};
    Set::Scalar fr_clip = 0.0, fr_cap = 0.0, fr_sheet = 0.0, fr_multi = 0.0;
    Set::Scalar fr_overlap = 0.0, fr_id_remap = 0.0, fr_ncell = 0.0;
    Set::Scalar fr_through = 0.0, fr_E_drop = 0.0, fr_fresh = 0.0;
    // P2-9 rev 2 item 7
    Set::Scalar fr_misaligned = 0.0, fr_noskin_c[3] = {0.0, 0.0, 0.0};
    Set::Scalar fr_t_removal = 0.0, fr_t_step = 0.0;   // wall seconds

    long long FrKey(int lev, int i, int j, int k) const
    {
        const amrex::Box& d = geom[lev].Domain();
        const long long nx = d.length(0), ny = d.length(1);
        return (static_cast<long long>(k - d.smallEnd(2)) * ny
                + (j - d.smallEnd(1))) * nx + (i - d.smallEnd(0));
    }

    // R0 over the current list: the surface cells and the line of each.
    std::vector<FrCell> FrCells() const
    {
        std::vector<FrCell> cells;
        std::size_t q = 0;
        while (q < fl_faces.size())
        {
            const FaceEntry& e0 = fl_faces[q];
            FrCell c;
            c.i = e0.i; c.j = e0.j; c.k = e0.k; c.col = e0.col;
            c.q0 = static_cast<int>(q);
            std::size_t r = q;
            while (r < fl_faces.size() && fl_faces[r].col == e0.col
                   && fl_faces[r].k == e0.k) ++r;
            c.nq = static_cast<int>(r - q);
            c.sheet = e0.fb;
            c.qline = c.q0;
            c.mu = 1.0;
            if (!c.sheet)
            {
                Set::Scalar best = -2.0;
                for (std::size_t p = q; p < r; ++p)
                {
                    const Set::Scalar v = FrNdotE(fl_faces[p], fl_faces[p].dir);
                    if (v > best) { best = v; c.qline = static_cast<int>(p); }
                }
                // ⚠ a normal no exposed face leans along cannot carry a line
                // (mu -> 0 sends the window to infinity): it is read as a sheet
                // and counted with them.
                if (!(best > 1.0e-6)) { c.sheet = 1; c.misaligned = 1; }
                else c.mu = best;
            }
            cells.push_back(c);
            q = r;
        }
        return cells;
    }

    // F4 under P2-9: every listed face takes the surface form at its CELL's
    // T_s, x w (the caller weights). A fallback cell's T_s is the face Newton
    // on its e* at this consumer's inputs -- so the kernel, the ledger, the
    // diagnostics and the march evaluate the identical number. The ONE call
    // every heat consumer makes under the key.
    void FrFaceFlux(const FaceEntry& e, Set::Scalar Tg, Set::Scalar kc,
                    Set::Scalar Tc, Set::Scalar dx,
                    Set::Scalar& q_in, Set::Scalar& q_out,
                    Set::Scalar& q_rob, Set::Scalar& q_los) const
    {
        Set::Scalar Ts = e.Ts;
        if (e.tsfb || !(Ts > 0.0))
        {
            const FaceEntry& el = fl_faces[static_cast<std::size_t>(e.qline)];
            Set::Scalar a = 0.0, b = 0.0, c = 0.0, d = 0.0, Tf = Tc;
            SurfaceCellFlux(true, true, 1.0, el.h, Tg, kc, dx, Tc,
                            eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                            a, b, c, d, 0.0, nullptr, 0.0, &Tf);
            Ts = Tf;
        }
        SurfaceCellFlux(true, true, 1.0, e.h, Tg, kc, dx, Tc,
                        eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                        q_in, q_out, q_rob, q_los, 0.0, nullptr, Ts);
    }

    // R1 at the list rebuild (the kernel, after FaceListHeat, before the
    // march): every listed cell gets its T_s. A cell the removal handed one to
    // (itself, or the cell its line's face moved to) keeps it -- one T_s per
    // step; a cell first exposed this step gets its own from the current T.
    void FrAssignTs(int lev)
    {
        if (!fl_ready) return;
        const Set::Scalar dx = geom[lev].CellSize()[0];
        const amrex::Box& dom = geom[lev].Domain();
        const std::vector<FrCell> cells = FrCells();
        const std::size_t nc = cells.size();
        // fresh cells: gather T(c) and T(c - e*) with its solidity
        std::vector<int> fresh(nc, 0);
        long n_fresh = 0;
        for (std::size_t n = 0; n < nc; ++n)
        {
            const FrCell& c = cells[n];
            const long long key = FrKey(lev, c.i, c.j, c.k);
            const auto it = fr_carry_valid ? fr_carry.find(key) : fr_carry.end();
            if (c.sheet)
            {
                // a sheet: every face its own line, its own face Newton
                for (int p = c.q0; p < c.q0 + c.nq; ++p)
                {
                    FaceEntry& e = fl_faces[static_cast<std::size_t>(p)];
                    e.Ts = 0.0; e.tsfb = 1; e.qline = p;
                }
                continue;
            }
            for (int p = c.q0; p < c.q0 + c.nq; ++p)
                fl_faces[static_cast<std::size_t>(p)].qline = c.qline;
            if (it != fr_carry.end())
            {
                for (int p = c.q0; p < c.q0 + c.nq; ++p)
                {
                    FaceEntry& e = fl_faces[static_cast<std::size_t>(p)];
                    e.Ts = it->second.fb ? 0.0 : it->second.Ts;
                    e.tsfb = it->second.fb;
                }
            }
            else { fresh[n] = 1; ++n_fresh; }
        }
        fr_fresh = static_cast<Set::Scalar>(n_fresh);
        if (n_fresh > 0)
        {
            std::vector<Set::Scalar> g(3 * nc, 0.0);   // T(c), T(c-e*), solid(c-e*)
            for (amrex::MFIter mfi(*temp_old_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& T = (*temp_old_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& rm = (*removed_mf[lev]).const_array(mfi);
                for (std::size_t n = 0; n < nc; ++n)
                {
                    if (!fresh[n]) continue;
                    const FrCell& c = cells[n];
                    const int dir = fl_faces[static_cast<std::size_t>(c.qline)].dir;
                    int ib = c.i, jb = c.j, kb = c.k;
                    const int sg = (dir % 2 == 1) ? -1 : 1;   // a step AGAINST e*
                    if (dir / 2 == 0) ib += sg; else if (dir / 2 == 1) jb += sg; else kb += sg;
                    if (bx.contains(amrex::IntVect(AMREX_D_DECL(c.i, c.j, c.k))))
                        g[3 * n] = T(c.i, c.j, c.k);
                    if (dom.contains(amrex::IntVect(AMREX_D_DECL(ib, jb, kb))) &&
                        bx.contains(amrex::IntVect(AMREX_D_DECL(ib, jb, kb))))
                    {
                        g[3 * n + 1] = T(ib, jb, kb);
                        g[3 * n + 2] = (rm(ib, jb, kb) > 0.5) ? 0.0 : 1.0;
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(g.data(), static_cast<int>(g.size()));
            const Set::Scalar T0r = phases.empty() ? 0.0 : phases[0].T_ref;
            const Set::Scalar cap = SurfaceGasTempCap();
            for (std::size_t n = 0; n < nc; ++n)
            {
                if (!fresh[n]) continue;
                const FrCell& c = cells[n];
                const Set::Scalar x = (g[3 * n + 2] > 0.5)
                    ? RawSurfaceX(g[3 * n], g[3 * n + 1], T0r) : 0.0;
                Set::Scalar Ts = 0.0;
                int fb = 1;
                if (x > 0.0)
                {
                    Ts = ReconFromX(g[3 * n], T0r, dx, x, cap, sp_flake_c_f,
                                    spall_h_spall_max / c.mu, sp_flake_h_max / c.mu).T_surf;
                    fb = 0;
                }
                for (int p = c.q0; p < c.q0 + c.nq; ++p)
                {
                    FaceEntry& e = fl_faces[static_cast<std::size_t>(p)];
                    e.Ts = Ts; e.tsfb = fb;
                }
            }
        }
        fr_carry.clear();
        fr_carry_valid = false;
    }

    // ---------------------------------------------------------------------
    // R3 -- THE ONE REMOVAL FUNCTION. Blind to what produced h. Along one
    // line of `run` solid cells with fractions f[0..run-1] and enthalpy
    // Ho[0..Lg-1]: consume h/dx of solid fraction from the face inward (a
    // cell reaching f <= eps voids, the remainder passes on) and remap the
    // heat (RemapProfile, over the scan window W = min(Dmax, run - 1)). Out:
    // n_voided; df[d] (fraction taken from cell d); dH[d] (the cell's
    // enthalpy change if it stays solid, -Ho[d] for a voided one); E_left,
    // the enthalpy that left the solid set (per unit volume), by the
    // TELESCOPED route -- the per-cell dH must sum to -E_left (fr_id_remap).
    // ---------------------------------------------------------------------
    static int RecedeAlongLine(int run, int Dmax, const Set::Scalar* f,
                               const Set::Scalar* Ho, Set::Scalar h,
                               Set::Scalar ell, Set::Scalar dx,
                               std::vector<Set::Scalar>& df,
                               std::vector<Set::Scalar>& dH,
                               Set::Scalar& f_after, Set::Scalar& E_left)
    {
        const Set::Scalar eps = 1.0e-9;
        df.assign(static_cast<std::size_t>(run), 0.0);
        dH.assign(static_cast<std::size_t>(run), 0.0);
        Set::Scalar rem = h / dx;
        int n_v = 0;
        f_after = 0.0;
        for (int d = 0; d < run; ++d)
        {
            if (rem >= f[d] - eps) { df[d] = f[d]; rem -= f[d]; ++n_v; if (rem < 0.0) rem = 0.0; }
            else { df[d] = rem; f_after = f[d] - rem; rem = 0.0; break; }
        }
        // the heat: the remap over the scan window, then the index shift by
        // the voided count -- the cell at old depth jd takes new[jd - n_v]
        const int W = std::min(Dmax, run - 1);
        std::vector<Set::Scalar> Hn(static_cast<std::size_t>(W + 1), 0.0);
        const Set::Scalar dE = RemapProfile(Ho, W, h, ell, dx, Hn.data());
        const int Wr = std::min(W, run - 1 - n_v);   // remapped values that land
        Set::Scalar E_bnd = 0.0, E_lost = 0.0;
        for (int jd = W + 1; jd <= Wr + n_v; ++jd) E_bnd += Ho[jd];
        for (int d = Wr + 1; d <= W; ++d) E_lost += Hn[static_cast<std::size_t>(d)];
        E_left = dE + E_bnd + E_lost;
        for (int jd = 0; jd < run; ++jd)
        {
            if (jd < n_v) { dH[static_cast<std::size_t>(jd)] = -Ho[jd]; continue; }
            const int d = jd - n_v;
            if (d > Wr) break;
            dH[static_cast<std::size_t>(jd)] = Hn[static_cast<std::size_t>(d)] - Ho[jd];
        }
        return n_v;
    }

    // ---------------------------------------------------------------------
    // The dispatch target: under the key this IS the removal step.
    // ---------------------------------------------------------------------
    void RecedeAlongLines(int lev, Set::Scalar time, Set::Scalar dt,
                         const amrex::Vector<int>& k_top_global)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar dx = DX[0];
        const Set::Scalar Vcell = DX[0] * DX[1] * DX[2];
        const amrex::Box& dom = geom[lev].Domain();
        const int nx_dom = dom.length(0);
        for (int c3 = 0; c3 < 3; ++c3)
        { fr_n_c[c3] = fr_v_c[c3] = fr_fb_c[c3] = 0.0; fr_smax_c[c3] = 0.0; }
        fr_clip = fr_cap = fr_sheet = fr_multi = fr_overlap = 0.0;
        fr_id_remap = fr_ncell = fr_through = fr_E_drop = 0.0;
        fr_misaligned = 0.0;
        fr_noskin_c[0] = fr_noskin_c[1] = fr_noskin_c[2] = 0.0;
        const double fr_t0 = amrex::second();
        fr_carry.clear();
        fr_carry_valid = false;

        // the solid fraction: 1 on solid, 0 on void, at first use
        FrInitFraction(lev);
        fr_score_mf[lev]->setVal(0.0);
        fr_mu_mf[lev]->setVal(0.0);

        bool any_write = false;
        // the list is built by the kernel; on the first step there is none
        // yet (nothing has heated, nothing can fire)
        if (fl_ready && !fl_faces.empty())
        {
            const std::vector<FrCell> cells = FrCells();
            fr_ncell = static_cast<Set::Scalar>(cells.size());

            // ---- R0: the lines ------------------------------------------
            const Set::Scalar hcap = (sp_flake_h_max > 0.0)
                ? std::min(sp_flake_h_max, spall_h_spall_max) : spall_h_spall_max;
            const int span = std::max(dom.length(0), std::max(dom.length(1), dom.length(2)));
            std::vector<FrLine> lines;
            std::size_t off = 0;
            auto add_line = [&] (int ci, int q, Set::Scalar mu) {
                FrLine L;
                L.cell = ci; L.qface = q; L.dir = fl_faces[static_cast<std::size_t>(q)].dir;
                L.mu = mu;
                L.Dmax = std::max(1, std::min(span,
                         static_cast<int>(std::floor(spall_h_spall_max / (mu * dx)))));
                L.Lg = std::min(span + 1, L.Dmax + 2
                       + static_cast<int>(std::ceil(hcap / (mu * dx))));
                // P2-14a: the constraint reference reads the live T d/mu further
                // along the line, so the gather reaches h_spall_max + d (normal
                // lengths, / mu) plus the column path's two-cell margin. ONLY the
                // reference reads past Lg: `run` (scan window, through test, R3) is
                // still counted over Lg. d = 0: Lr = Lg, the gather is unchanged.
                L.Lr = L.Lg;
                if (sp_cref_d > 0.0)
                    L.Lr = std::min(span + 1, std::max(L.Lg,
                           static_cast<int>(std::floor((spall_h_spall_max + sp_cref_d) / (mu * dx))) + 2));
                L.off = off;
                off += static_cast<std::size_t>(L.Lr);
                lines.push_back(L);
            };
            for (std::size_t n = 0; n < cells.size(); ++n)
            {
                const FrCell& c = cells[n];
                if (c.sheet)
                {
                    // zero gradient (fr_sheet) and misaligned (fr_misaligned)
                    // are evaluated alike, counted apart
                    if (c.misaligned) fr_misaligned += 1.0; else fr_sheet += 1.0;
                    for (int p = c.q0; p < c.q0 + c.nq; ++p)
                        add_line(static_cast<int>(n), p, 1.0);
                }
                else add_line(static_cast<int>(n), c.qline, c.mu);
            }
            const std::size_t nl = lines.size();

            // ---- THE ONE GATHER: T, T0, H, solid, f along every line, plus
            // the surface cell's flaw and conductivity. Owner-written,
            // ReduceRealSum-collapsed (PASS 1b pattern): each cell has exactly
            // one owner, so the sum is exact and decomposition-free.
            constexpr int NF = 5;
            std::vector<Set::Scalar> G(off * NF, 0.0), Gc(2 * nl, 0.0);
            const bool have_T0 = T_init_mf.size() > static_cast<std::size_t>(lev)
                                 && T_init_mf[lev];
            for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& T  = (*temp_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& H  = (*H_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& rm = (*removed_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& ff = (*fr_f_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& kp = (*kappa_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& Df = (*D_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> T0a, fa;
                if (have_T0) T0a = (*T_init_mf[lev]).const_array(mfi);
                if (weibull_enabled) fa = (*flaw_a_mf[lev]).const_array(mfi);
                const int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
                const int hi[3] = {bx.bigEnd(0), bx.bigEnd(1), bx.bigEnd(2)};
                for (std::size_t l = 0; l < nl; ++l)
                {
                    const FrLine& L = lines[l];
                    const FrCell& c = cells[static_cast<std::size_t>(L.cell)];
                    const int a = L.dir / 2;
                    const int sg = (L.dir % 2 == 1) ? -1 : 1;   // a step AGAINST e*
                    int p0[3] = {c.i, c.j, c.k};
                    bool in = true;
                    for (int t = 0; t < 3; ++t)
                        if (t != a && (p0[t] < lo[t] || p0[t] > hi[t])) in = false;
                    if (!in) continue;
                    // positions p0[a] + sg*d inside [lo[a], hi[a]]
                    int dlo, dhi;
                    if (sg > 0) { dlo = lo[a] - p0[a]; dhi = hi[a] - p0[a]; }
                    else        { dlo = p0[a] - hi[a]; dhi = p0[a] - lo[a]; }
                    dlo = std::max(dlo, 0);
                    dhi = std::min(dhi, L.Lr - 1);
                    for (int d = dlo; d <= dhi; ++d)
                    {
                        int p[3] = {p0[0], p0[1], p0[2]};
                        p[a] += sg * d;
                        const std::size_t o = (L.off + static_cast<std::size_t>(d)) * NF;
                        G[o + 0] = T(p[0], p[1], p[2]);
                        G[o + 1] = have_T0 ? T0a(p[0], p[1], p[2])
                                           : std::numeric_limits<Set::Scalar>::quiet_NaN();
                        G[o + 2] = H(p[0], p[1], p[2]);
                        G[o + 3] = (rm(p[0], p[1], p[2]) > 0.5) ? 0.0 : 1.0;
                        G[o + 4] = ff(p[0], p[1], p[2]);
                        if (d == 0)
                        {
                            Gc[2 * l] = weibull_enabled ? fa(p[0], p[1], p[2]) : sp_a0;
                            Gc[2 * l + 1] = (kappa_damage_alpha != 0.0)
                                ? kp(p[0], p[1], p[2]) * std::exp(-kappa_damage_alpha * Df(p[0], p[1], p[2]))
                                : kp(p[0], p[1], p[2]);
                        }
                    }
                }
            }
            // ⚠ NaN marks "no reference field" -- it must not enter a sum
            if (!have_T0) for (std::size_t o = 0; o < off; ++o) G[o * NF + 1] = 0.0;
            if (!G.empty())
                amrex::ParallelDescriptor::ReduceRealSum(G.data(), static_cast<int>(G.size()));
            if (!Gc.empty())
                amrex::ParallelDescriptor::ReduceRealSum(Gc.data(), static_cast<int>(Gc.size()));

            // ---- R1 + R2 on every line (every rank, identical) ----------
            const Set::Scalar T0r = phases.empty() ? 0.0 : phases[0].T_ref;
            const Set::Scalar cap = SurfaceGasTempCap();
            const Numeric::Material::Grimsel *scan_g = RockLaw();
            const ElasticConsts ec = ElasticConstants();
            const Set::Scalar pE = ec.E, pB = ec.beta, pNu = ec.nu, conf_off = ec.conf;
            const int n_seg = (spall_h_spall_n_segments > 0) ? spall_h_spall_n_segments
                                                            : sp_n_segments;
            struct LineOut
            {
                int run = 0, fire = 0, through = 0, fb = 0, noskin = 0;
                Set::Scalar q_net = 0.0;   // the fallback's face balance (logged)
                Set::Scalar score = 0.0, h_col = 0.0;
                Set::Scalar T_cref = 0.0;  // P2-14a: the reference at z_v = 0 (logged)
                SurfaceRecon rec;
            };
            std::vector<LineOut> out(nl);
            // P2-14a: the constraint reference's witnesses on the LINES. ⚠ Not
            // reduced: every rank evaluates every line from the same reduced
            // gather (R1 + R2 below), so a ReduceSum would count each read once
            // per rank. [0] floor-class reads, [1] wall/roof-class reads,
            // [2] clamps, [3] T_z < T0_z.
            const bool cref_fl = sp_cref_d > 0.0;
            long cref_cnt[4] = {0, 0, 0, 0};
            for (std::size_t l = 0; l < nl; ++l)
            {
                const FrLine& L = lines[l];
                const Set::Scalar* Gl = G.data() + L.off * NF;
                LineOut& o = out[l];
                int run = 0;
                while (run < L.Lg && Gl[run * NF + 3] > 0.5) ++run;
                o.run = run;
                if (run < 1) continue;   // cannot happen: the cell is solid
                // P2-14a: the reference's own solid run, over the extended gather
                // (= run at d = 0, where it is not read)
                int run_ref = run;
                if (cref_fl)
                    while (run_ref < L.Lr && Gl[run_ref * NF + 3] > 0.5) ++run_ref;
                const Set::Scalar d_line = sp_cref_d / L.mu;   // the normal d along the line
                auto T_line = [Gl] (int m) -> Set::Scalar { return Gl[m * NF]; };
                long& n_ref = cref_cnt[FrClass(L.dir) == 0 ? 0 : 1];
                const Set::Scalar Tc = Gl[0];
                const Set::Scalar x = (run >= 2) ? RawSurfaceX(Tc, Gl[NF], T0r) : 0.0;
                if (x > 0.0)
                {
                    o.rec = ReconFromX(Tc, T0r, dx, x, cap, sp_flake_c_f,
                                       spall_h_spall_max / L.mu, sp_flake_h_max / L.mu);
                }
                else
                {
                    // P2-9 rev 2 item 5 -- THE FALLBACK, from the face heat balance
                    // on e*: T_s = the face Newton's root, q_net = h (T_g - T_s) -
                    // losses on e*, ell_q = k_c (T_s - T0) / q_net, the flake
                    // c_f ell_q under ReconFromX's two caps (/ mu). No heated skin
                    // (q_net <= 0 or T_s <= T0): no recession (fr_noskin). Nothing
                    // defaults to one cell -- ReconFromX's flat h_f = dz is NOT read.
                    o.rec = SurfaceRecon();
                    o.rec.flat = true;
                    o.rec.ell_eff = std::numeric_limits<Set::Scalar>::infinity();
                    const FaceEntry& ef = fl_faces[static_cast<std::size_t>(L.qface)];
                    Set::Scalar a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, Tf = Tc;
                    // P2-11b item 5: the face's own gas temperature (under
                    // face_heat = position it is flame_Tg[c.col] exactly)
                    const Set::Scalar Tg = ef.Tg;
                    SurfaceCellFlux(true, true, 1.0, ef.h, Tg, Gc[2 * l + 1], dx, Tc,
                                    eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                                    a1, a2, a3, a4, 0.0, nullptr, 0.0, &Tf);
                    o.rec.T_surf = Tf;
                    o.rec.T_bar = Tf;
                    o.fb = 1;
                    const Set::Scalar q_net = a3 - a4;
                    const Set::Scalar k_c = Gc[2 * l + 1];
                    if (q_net > 0.0 && Tf > T0r && k_c > 0.0)
                    {
                        const Set::Scalar ell_q = k_c * (Tf - T0r) / q_net;
                        Set::Scalar hf = sp_flake_c_f * ell_q;
                        const Set::Scalar h_max = spall_h_spall_max / L.mu;
                        const Set::Scalar h_bnd = sp_flake_h_max / L.mu;
                        if (h_max > 0.0 && hf > h_max) { hf = h_max; o.rec.h_f_capped = true; }
                        if (h_bnd > 0.0 && hf > h_bnd) { hf = h_bnd; o.rec.h_f_bounded = true; }
                        o.rec.h_f = hf;
                        o.rec.ell_eff = ell_q;
                        o.rec.ell_q = ell_q;
                        o.q_net = q_net;
                    }
                    else
                    {
                        o.rec.h_f = 0.0;
                        o.noskin = 1;
                        o.q_net = q_net;
                    }
                }
                // the scan window along the line, clamped to the solid run
                const int nprof = std::min(L.Dmax, run - 1);
                auto T_at_d = [&] (int d) -> Set::Scalar {
                    if (d < 0) d = 0;
                    if (d > nprof) d = nprof;
                    return Gl[d * NF];
                };
                auto T0_at_d = [&] (int d) -> Set::Scalar {
                    if (!have_T0) return T0r;
                    if (d < 0) d = 0;
                    if (d > nprof) d = nprof;
                    const Set::Scalar v = Gl[d * NF + 1];
                    return std::isfinite(v) ? v : T0r;
                };
                const Set::Scalar Ts = o.rec.T_surf;
                auto T_at_zv = [&] (Set::Scalar z_v) -> Set::Scalar {
                    return ScanTAtZv(T_at_d,
                        [&] (Set::Scalar) -> Set::Scalar { return Ts; }, dx, z_v);
                };
                const Set::Scalar a_f = Gc[2 * l];
                auto k_i_at_zcrack = [&] (Set::Scalar z_crack) -> Set::Scalar {
                    auto sigma_at_depth = [&] (Set::Scalar d) -> Set::Scalar {
                        const Set::Scalar z_v = z_crack + d;
                        const Set::Scalar T_z = T_at_zv(z_v);
                        // ⚠⚠ P2-14a: under constraint_depth > 0 the LIVE T a
                        // normal distance d further behind, i.e. d/mu along the
                        // line, through the column path's own arithmetic
                        // (CrefInterp), in METRES over the line's solid run.
                        // d = 0: the stored initial T at the sampled cell.
                        const Set::Scalar T0_z = cref_fl
                            ? CrefInterp(T_line, run_ref, z_v + d_line, dx, n_ref, cref_cnt[2])
                            : T0_at_d(static_cast<int>(std::floor(z_v / dx)));
                        if (cref_fl && T_z < T0_z) ++cref_cnt[3];
                        return ScanSigmaT(scan_g, T_z, T0_z, 1.0, conf_off, pE, pB, pNu);
                    };
                    return Numeric::SpCriterion::K_I(sigma_at_depth, a_f, n_seg);
                };
                auto kic_pa_at_z = [&] (Set::Scalar z_v) -> Set::Scalar {
                    return ScanKIcPa(T_at_zv(z_v));
                };
                auto kic_pa_at_cell = [&] (int k_tip) -> Set::Scalar {
                    return ScanKIcPa(T_at_d(-k_tip));
                };
                int s_last = -1;
                bool failed = false;
                Set::Scalar score = 0.0;
                if (std::isfinite(a_f) && a_f > 0.0)
                    o.h_col = DepthScanH(nprof, 0, -(L.Lg + 1), dx, 1.0, 1.0, a_f, true,
                                         k_i_at_zcrack, kic_pa_at_z, kic_pa_at_cell,
                                         [] (int) {}, s_last, failed, score);
                o.score = score;
                o.fire = (o.h_col > 0.0) ? 1 : 0;
                if (cref_fl)
                {
                    long nr = 0, nc = 0;   // the logged value is not a criterion read
                    o.T_cref = CrefInterp(T_line, run_ref, d_line, dx, nr, nc);
                }
                // the crack reaches the void behind the line (a sheet, or any
                // line whose run ends inside the window): no rock to stop it
                if (!o.fire && score >= 1.0 && !failed && s_last == nprof
                    && nprof == run - 1 && run - 1 < L.Dmax)
                { o.fire = 1; o.through = 1; }
            }

            // P2-14a: the lines' share of the witnesses (BuildConstraintRefColumns
            // zeroed them this step; the Sp_field site adds nothing under the list)
            if (cref_fl)
            {
                cref_n_vert  += static_cast<Set::Scalar>(cref_cnt[0]);
                cref_n_lat   += static_cast<Set::Scalar>(cref_cnt[1]);
                cref_n_clamp += static_cast<Set::Scalar>(cref_cnt[2]);
                cref_n_neg   += static_cast<Set::Scalar>(cref_cnt[3]);
            }

            // ---- one line per cell: a sheet keeps its largest score -------
            std::vector<int> win(cells.size(), -1);
            for (std::size_t l = 0; l < nl; ++l)
            {
                const int ci = lines[l].cell;
                if (win[static_cast<std::size_t>(ci)] < 0
                    || out[l].score > out[static_cast<std::size_t>(win[static_cast<std::size_t>(ci)])].score)
                    win[static_cast<std::size_t>(ci)] = static_cast<int>(l);
            }

            // ---- R3: recede every firing line, per-cell deltas -----------
            struct Upd { Set::Scalar dH = 0.0, df = 0.0, H0 = 0.0, f0 = 0.0;
                         Set::Scalar h_line = 0.0; int nl = 0; int phiz = 0; };
            std::map<long long, Upd> upd;
            std::map<long long, std::array<int, 3>> where;
            struct Ev { std::size_t l; int n_v; Set::Scalar h_line, f_after; };
            std::vector<Ev> evs;
            Set::Scalar E_left_tot = 0.0, dH_book = 0.0;
            std::vector<int> fired_cell(cells.size(), 0);
            for (std::size_t ci = 0; ci < cells.size(); ++ci)
            {
                const int l = win[ci];
                if (l < 0) continue;
                const FrLine& L = lines[static_cast<std::size_t>(l)];
                const LineOut& o = out[static_cast<std::size_t>(l)];
                const int cls = FrClass(L.dir);
                if (o.fb) fr_fb_c[cls] += 1.0;
                if (o.noskin) fr_noskin_c[cls] += 1.0;
                // no heated skin and no crack to the void: nothing to recede
                if (o.fire && o.noskin && !o.through)
                {
                    fr_smax_c[cls] = std::max(fr_smax_c[cls], o.score);
                    continue;
                }
                if (!o.fire)
                {
                    fr_smax_c[cls] = std::max(fr_smax_c[cls], o.score);
                    continue;
                }
                if (fired_cell[ci]++) fr_multi += 1.0;
                if (o.through) fr_through += 1.0;
                const Set::Scalar* Gl = G.data() + L.off * NF;
                std::vector<Set::Scalar> f(static_cast<std::size_t>(o.run)),
                                         Ho(static_cast<std::size_t>(o.run));
                Set::Scalar L_av = 0.0;
                for (int d = 0; d < o.run; ++d)
                {
                    f[static_cast<std::size_t>(d)] = Gl[d * NF + 4];
                    Ho[static_cast<std::size_t>(d)] = Gl[d * NF + 2];
                    L_av += f[static_cast<std::size_t>(d)] * dx;
                }
                // P2-9 rev 2 item 6: THROUGH MEANS THROUGH -- the crack reached
                // the void behind the line, so the whole solid along the run
                // comes off (<= h_spall_max / mu: the run ends inside the window).
                // Otherwise the flake, clipped to the solid along the line.
                Set::Scalar h = o.through ? L_av : o.rec.h_f;
                if (!o.through)
                {
                    if (o.rec.h_f_capped || o.rec.h_f_bounded) fr_cap += 1.0;
                    if (h > L_av) { h = L_av; fr_clip += 1.0; }
                }
                if (!(h > 0.0)) continue;
                const Set::Scalar ell = std::isfinite(o.rec.ell_eff) && o.rec.ell_eff > 0.0
                    ? o.rec.ell_eff : std::numeric_limits<Set::Scalar>::infinity();
                std::vector<Set::Scalar> df, dH;
                Set::Scalar f_after = 0.0, E_left = 0.0;
                const int n_v = RecedeAlongLine(o.run, L.Dmax, f.data(), Ho.data(), h, ell,
                                                dx, df, dH, f_after, E_left);
                fr_n_c[cls] += 1.0;
                fr_v_c[cls] += static_cast<Set::Scalar>(n_v);
                E_left_tot += E_left;
                Set::Scalar dsum = 0.0;
                const FrCell& c = cells[ci];
                const int a = L.dir / 2;
                const int sg = (L.dir % 2 == 1) ? -1 : 1;
                for (int d = 0; d < o.run; ++d)
                {
                    if (df[static_cast<std::size_t>(d)] == 0.0 && dH[static_cast<std::size_t>(d)] == 0.0)
                        continue;
                    int p[3] = {c.i, c.j, c.k};
                    p[a] += sg * d;
                    const long long key = FrKey(lev, p[0], p[1], p[2]);
                    Upd& u = upd[key];
                    u.dH += dH[static_cast<std::size_t>(d)];
                    u.df += df[static_cast<std::size_t>(d)];
                    u.H0 = Ho[static_cast<std::size_t>(d)];
                    u.f0 = f[static_cast<std::size_t>(d)];
                    u.h_line = h;
                    u.nl += 1;
                    where[key] = {p[0], p[1], p[2]};
                    dsum += dH[static_cast<std::size_t>(d)];
                }
                dH_book += dsum;
                // phi follows a +z line's column, as the column path moves it
                if (L.dir == 5)
                {
                    const long long key = FrKey(lev, c.i, c.j, c.k);
                    upd[key].phiz = 1;
                    upd[key].h_line = h;
                    where[key] = {c.i, c.j, c.k};
                }
                evs.push_back({static_cast<std::size_t>(l), n_v, h, f_after});
                // the value moves with the face: T_s to the cell the line's
                // face now sits on (or stays, if nothing voided)
                if (n_v < o.run)
                {
                    int p[3] = {c.i, c.j, c.k};
                    p[a] += sg * n_v;
                    const long long key = FrKey(lev, p[0], p[1], p[2]);
                    if (fr_carry.find(key) == fr_carry.end())
                        fr_carry[key] = {o.rec.T_surf, o.fb};
                }
            }
            // the cells that did not fire keep their own T_s
            for (std::size_t ci = 0; ci < cells.size(); ++ci)
            {
                const int l = win[ci];
                if (l < 0 || out[static_cast<std::size_t>(l)].fire) continue;
                const long long key = FrKey(lev, cells[ci].i, cells[ci].j, cells[ci].k);
                if (fr_carry.find(key) == fr_carry.end())
                    fr_carry[key] = {out[static_cast<std::size_t>(l)].rec.T_surf,
                                     cells[ci].sheet ? 1 : out[static_cast<std::size_t>(l)].fb};
            }
            fr_carry_valid = true;
            fr_id_remap = (E_left_tot > 0.0)
                ? std::abs(dH_book + E_left_tot) / E_left_tot : 0.0;

            // ---- overlap, and the energy a voided cell carries off another
            // line's books: a cell that ends void leaves the solid set with
            // its OLD enthalpy, whatever deltas the lines booked on it, so
            // the difference (sum of booked dH + H_old) is credited here.
            const Set::Scalar eps_f = 1.0e-9;
            Set::Scalar E_drop = 0.0;
            for (auto& kv : upd)
            {
                Upd& u = kv.second;
                if (u.nl > 1) fr_overlap += 1.0;
                if (u.nl > 0 && u.f0 - u.df <= eps_f) E_drop += u.dH + u.H0;
            }
            fr_E_drop = E_drop * Vcell;
            ledger_E_flake += (E_left_tot + E_drop) * Vcell;

            // ---- ONE WRITE: f, removed, H, phi and the per-step fields ----
            std::vector<long long> keys;
            keys.reserve(upd.size());
            for (const auto& kv : upd) keys.push_back(kv.first);
            for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<Set::Scalar> const& rm  = (*removed_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& H   = (*H_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& ff  = (*fr_f_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& reg = (*regime_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& ropU= (*RoP_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& D   = (*D_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& regf= (*regime_field_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& ev_s  = (*spall_event_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& thk_s = (*spall_thickness_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& rop_s = (*RoP_spall_mf[lev]).array(mfi);
                for (const long long key : keys)
                {
                    const Upd& u = upd[key];
                    const std::array<int, 3>& p = where[key];
                    if (u.phiz)
                    {
                        // the whole column's phi, the column path's PASS 3 shift
                        if (p[0] >= bx.smallEnd(0) && p[0] <= bx.bigEnd(0) &&
                            p[1] >= bx.smallEnd(1) && p[1] <= bx.bigEnd(1))
                            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
                                phi(p[0], p[1], k) -= u.h_line;
                    }
                    if (u.nl == 0) continue;
                    if (!bx.contains(amrex::IntVect(AMREX_D_DECL(p[0], p[1], p[2])))) continue;
                    const Set::Scalar fn = u.f0 - u.df;
                    if (fn <= eps_f)
                    {
                        // void: the thermal state stays frozen (never reset)
                        ff(p[0], p[1], p[2]) = 0.0;
                        rm(p[0], p[1], p[2]) = 1.0;
                        reg(p[0], p[1], p[2]) = 1.0;
                        regf(p[0], p[1], p[2]) = 1.0;
                        ropU(p[0], p[1], p[2]) = (dt > 0.0) ? u.h_line / dt : 0.0;
                        ev_s(p[0], p[1], p[2]) = 1.0;
                        thk_s(p[0], p[1], p[2]) = u.h_line;
                        rop_s(p[0], p[1], p[2]) = (dt > 0.0) ? u.h_line / dt : 0.0;
                        D(p[0], p[1], p[2]) = 0.0;
                    }
                    else
                    {
                        ff(p[0], p[1], p[2]) = fn;
                        H(p[0], p[1], p[2]) += u.dH;
                    }
                    any_write = true;
                }
                // the witnesses on every surface cell, fired or not
                amrex::Array4<Set::Scalar> const& fsc = (*fr_score_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& fmu = (*fr_mu_mf[lev]).array(mfi);
                for (std::size_t ci = 0; ci < cells.size(); ++ci)
                {
                    const FrCell& c = cells[ci];
                    if (!bx.contains(amrex::IntVect(AMREX_D_DECL(c.i, c.j, c.k)))) continue;
                    const int l = win[ci];
                    fsc(c.i, c.j, c.k) = (l >= 0) ? out[static_cast<std::size_t>(l)].score : 0.0;
                    fmu(c.i, c.j, c.k) = (l >= 0) ? lines[static_cast<std::size_t>(l)].mu : 0.0;
                }
            }
            {
                int aw = any_write ? 1 : 0;
                amrex::ParallelDescriptor::ReduceIntMax(aw);
                any_write = (aw != 0);
            }

            // ---- the event log: one row per firing ----------------------
            if (sp_emit_removal_events_csv && amrex::ParallelDescriptor::IOProcessor())
            {
                const std::string csv_path = plot_file + "_removal_events.csv";
                const bool need_header = !sp_removal_events_csv_header_written;
                std::ofstream ofs(csv_path, (need_header ? std::ios::trunc : std::ios::app)
                                            | std::ios::out);
                if (ofs.is_open())
                {
                    if (need_header)
                        ofs << "time,col_i,col_j,regime,k_top,T_top,a_f,Sp_top,"
                               "h_scan,h_applied,n_voided"
                            << ",T_surf,ell_eff,h_f,T_bar_flake,"
                               "rc_flag,d_ell_rel,ell_raw,delta,q_ratio"
                            // P2-9 deliverable 6, appended: the line
                            << ",dir,mu,k,run,h_line,h_n,f_after,score,ts_fb,sheet,"
                               "through,q_net,k_c"
                            // P2-14a: under d > 0 only, the line's reference
                            << (cref_fl ? ",T_cref" : "")
                            << "\n";
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    for (const Ev& ev : evs)
                    {
                        const FrLine& L = lines[ev.l];
                        const LineOut& o = out[ev.l];
                        const FrCell& c = cells[static_cast<std::size_t>(L.cell)];
                        const Set::Scalar* Gl = G.data() + L.off * NF;
                        const int flag = (o.rec.flat ? 1 : 0) + (o.rec.ell_clamped ? 2 : 0)
                                       + (o.rec.T_clamped ? 4 : 0) + (o.rec.h_f_capped ? 8 : 0)
                                       + (o.rec.h_f_bounded ? 16 : 0);
                        ofs << time << "," << c.i << "," << c.j << ",1,"
                            << k_top_global[static_cast<std::size_t>(c.col)] << ","
                            << Gl[0] << "," << Gc[2 * ev.l] << "," << 0.0 << ","
                            << o.h_col << "," << ev.h_line << "," << ev.n_v
                            << "," << o.rec.T_surf << ","
                            << (std::isfinite(o.rec.ell_eff) ? o.rec.ell_eff : 0.0)
                            << "," << o.rec.h_f << "," << o.rec.T_bar << ","
                            << Set::Scalar(flag) << ",0,"
                            << (std::isfinite(o.rec.ell_eff) ? o.rec.ell_eff : 0.0)
                            << ",0,0";
                        ofs << "," << L.dir << "," << L.mu << "," << c.k << "," << o.run
                            << "," << ev.h_line << "," << L.mu * ev.h_line
                            << "," << ev.f_after << "," << o.score << "," << o.fb
                            << "," << c.sheet << "," << o.through << "," << o.q_net
                            << "," << Gc[2 * ev.l + 1];
                        if (cref_fl) ofs << "," << o.T_cref;
                        ofs << "\n";
                    }
                }
            }
            if (sp_emit_removal_events_csv) sp_removal_events_csv_header_written = true;
            (void) nx_dom;
        }

        // the tail every removal step runs (Removal.H's own): T back from H,
        // ghosts refilled, the derived masks rebuilt
        if (any_write)
        {
            InvertHtoT(lev);
            Util::RealFillBoundary(*temp_mf[lev], geom[lev]);
            // P2-9 rev 2 item 2: the periodic + physical ghosts, through the
            // solver's own BC (RefillTempGhosts, the column path's call too).
            RefillTempGhosts(lev, time);
        }
        UpdateSurfaceMaskFromPhi(lev);
        RefreshKappaEff(lev);
        {
            double tr = amrex::second() - fr_t0;
            amrex::ParallelDescriptor::ReduceRealMax(tr);
            fr_t_removal = static_cast<Set::Scalar>(tr);
        }
    }

#endif
