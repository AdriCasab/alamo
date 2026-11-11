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
        int Lg = 0;                // gathered length (Dmax + the shift's reach)
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
    Set::Scalar fr_misaligned = 0.0, fr_nodepth_c[3] = {0.0, 0.0, 0.0};
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
    // P2-15: THE BUCKLING RATIO, one helper, one caller (RecedeAlongLines).
    // A compressed flake of thickness x buckles over a span L only if
    // L / x >= C_L, C_L = pi / sqrt(6 (1 + nu) eps_th) (Rauenzahn & Tester 1989
    // Eqs. 15-16, after Preston and Dey & Kranz). eps_th is the flake's thermal
    // strain over the substrate: under the rock law the ONE strain table
    // (Grimsel::strain_of, which takes K), else beta (T_s - T0_s). nu is
    // ElasticConstants().nu (the law's under RockLaw(), P2-14a). Returns
    // bk = C_L x / L_span; +inf for a non-finite x or eps_th <= 0, so the
    // caller's `blocked = !(bk <= 1)` withdraws them.
    // ---------------------------------------------------------------------
    static Set::Scalar BuckleRatio(const Numeric::Material::Grimsel* law, Set::Scalar nu,
                                   Set::Scalar beta, Set::Scalar T_s, Set::Scalar T0_s,
                                   Set::Scalar x, Set::Scalar L_span)
    {
        const Set::Scalar inf = std::numeric_limits<Set::Scalar>::infinity();
        if (!std::isfinite(x)) return inf;
        const Set::Scalar eps_th = law ? law->strain_of(T_s) - law->strain_of(T0_s)
                                       : beta * (T_s - T0_s);
        if (!(eps_th > 0.0) || !(1.0 + nu > 0.0) || !(L_span > 0.0)) return inf;
        const Set::Scalar C_L = M_PI / std::sqrt(6.0 * (1.0 + nu) * eps_th);
        return C_L * x / L_span;
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
        fr_nodepth_c[0] = fr_nodepth_c[1] = fr_nodepth_c[2] = 0.0;
        fr_n_molten = fr_n_lid = fr_n_vapor = fr_n_vapor_int = 0.0;   // P2-24, per step
        // P2-15: the buckling witnesses, per step (registered only under the key)
        bk_span_out = sp_buckle_span;
        bk_n_eval = bk_n_block = bk_n_inf = 0.0;
        bk_n_block_floor = bk_n_block_wall = bk_n_block_roof = 0.0;
        bk_max_pass = 0.0;
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
                L.off = off;
                off += static_cast<std::size_t>(L.Lg);
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
                    dhi = std::min(dhi, L.Lg - 1);
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

            // ---- P2-24: THE VAPOUR RULE, detected before the scan. A LISTED cell (touching the void) whose
            // enthalpy has passed its phase table's H_vap_hi -- it has received its whole vaporisation enthalpy,
            // through deposition and conduction, ONCE -- leaves in the ONE WRITE below as regime 2, booking H V
            // (ledger_E_vapor). Its lines are not scanned. An INTERIOR cell past H_vap_hi (buried under a solid
            // lid: the beam's peak is interior) is NOT voided -- a sealed bubble has no physics here -- and is
            // counted (fr_n_vapor_int). ONE rock phase (SR 13; Parse aborts otherwise).
            const bool vap_on = vapor_enabled && !phase_tables.empty();
            const Set::Scalar H_vap = vap_on ? phase_tables[0].H_vap_hi() : 0.0;
            std::vector<char> vap_cell(cells.size(), 0);
            std::vector<Set::Scalar> vap_H(cells.size(), 0.0), vap_T(cells.size(), 0.0), vap_f(cells.size(), 0.0);
            if (vap_on)
                for (std::size_t l = 0; l < nl; ++l)
                {
                    const std::size_t ci = static_cast<std::size_t>(lines[l].cell);
                    if (vap_cell[ci]) continue;
                    const Set::Scalar* Gl = G.data() + lines[l].off * NF;
                    if (Gl[3] > 0.5 && Gl[2] >= H_vap)
                    { vap_cell[ci] = 1; vap_H[ci] = Gl[2]; vap_T[ci] = Gl[0]; vap_f[ci] = Gl[4]; }
                }
            if (vap_on && fl_index)
            {
                // the INTERIOR count: solid cells past H_vap_hi that touch no void (not listed)
                Set::Scalar n_int = 0.0;
                for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    amrex::Array4<const Set::Scalar> const& rm = (*removed_mf[lev]).const_array(mfi);
                    amrex::Array4<const Set::Scalar> const& Hc = (*H_mf[lev]).const_array(mfi);
                    amrex::Array4<const int> const& fi = fl_index->const_array(mfi);
                    const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                    for (int k = lo.z; k <= hi.z; ++k)
                    for (int j = lo.y; j <= hi.y; ++j)
                    for (int i = lo.x; i <= hi.x; ++i)
                        if (rm(i, j, k) <= 0.5 && Hc(i, j, k) >= H_vap && fi(i, j, k, 1) == 0) n_int += 1.0;
                }
                amrex::ParallelDescriptor::ReduceRealSum(n_int);
                fr_n_vapor_int = n_int;
                bc_n_vapor_int += n_int;
            }
            long n_molten = 0, n_lid = 0;

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
                int run = 0, fire = 0, through = 0, fb = 0, nodepth = 0;
                Set::Scalar q_net = 0.0;   // the fallback's face balance (logged)
                Set::Scalar score = 0.0, h_col = 0.0;
                Set::Scalar bk = 0.0;      // P2-15: C_L x / L_span on a line that fired (logged)
                // P2-19: logged on a fired line, read by nothing (pure reads; see below)
                Set::Scalar ell_p = std::numeric_limits<Set::Scalar>::quiet_NaN();
                Set::Scalar z_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();
                Set::Scalar ell_pk = std::numeric_limits<Set::Scalar>::quiet_NaN();
                Set::Scalar h_E = std::numeric_limits<Set::Scalar>::quiet_NaN();
                int bk_blocked = 0;
                SurfaceRecon rec;
            };
            std::vector<LineOut> out(nl);
            // P2-15: the buckling condition's witnesses on the lines. Not rank-summed for the
            // same reason (every rank evaluates every line): [0] evaluated, [1] blocked,
            // [2] blocked with non-finite x, [3..5] blocked by class floor/wall/roof.
            const bool bk_on = sp_buckle_span > 0.0;
            long bk_cnt[6] = {0, 0, 0, 0, 0, 0};
            Set::Scalar bk_pass_max = 0.0;
            // P2-23: spall.buckle_log_csv -- one row per line B evaluates, on the
            // steps whose thermo row is written (the row this step's values reach,
            // amr.thermo.plot_int: thermo is written before the NEXT step's advance,
            // so this step is step istep + 1). Every rank evaluates every line; the
            // I/O rank writes.
            struct BkRow { std::size_t l = 0; int through = 0, blocked = 0;
                           Set::Scalar x = 0.0, bk = 0.0, bk_lb = 0.0; };
            std::vector<BkRow> bk_rows;
            const bool bk_log = bk_on && sp_buckle_log_csv && thermo.plot_int > 0
                                && ((istep[lev] + 1) % thermo.plot_int == 0);
            for (std::size_t l = 0; l < nl; ++l)
            {
                const FrLine& L = lines[l];
                const Set::Scalar* Gl = G.data() + L.off * NF;
                LineOut& o = out[l];
                if (vap_cell[static_cast<std::size_t>(L.cell)]) continue;   // P2-24: it leaves as vapour
                int run = 0;
                while (run < L.Lg && Gl[run * NF + 3] > 0.5) ++run;
                o.run = run;
                if (run < 1) continue;   // cannot happen: the cell is solid
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
                // P2-19: the line's own temperature rise over the scan window (the
                // fallback's depth and the logged diagnostics read it; nothing else)
                auto rise_at_d = [&] (int d) -> Set::Scalar { return T_at_d(d) - T0_at_d(d); };
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
                    // on e*: T_s = the face Newton's root (sign-agnostic: it is the
                    // surface temperature whether the face gains or loses heat).
                    // ⚠⚠ P2-19: the flake depth is the E-FOLD DEPTH OF THE LINE'S OWN
                    // RISE from T_s - T0 (ProfileEFoldDepth, Shared.H), defined on any
                    // profile -- no longer ell_q = k_c (T_s - T0) / q_net, which is
                    // undefined when the face loses heat. The flake is c_f ell_p under
                    // ReconFromX's two caps (/ mu). A rise that never e-folds inside
                    // the window has no depth: no recession (fr_nodepth). q_net and the
                    // old tangent depth ell_q are kept for the log only. Nothing defaults
                    // to one cell -- ReconFromX's flat h_f = dz is NOT read. (o.fb is set
                    // on every fallback evaluation, as before, so fr_fb_* is unchanged.)
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
                    o.q_net = q_net;
                    if (q_net > 0.0 && Tf > T0r && k_c > 0.0)
                        o.rec.ell_q = k_c * (Tf - T0r) / q_net;   // the old depth (log only)
                    const Set::Scalar ell_p = ProfileEFoldDepth(rise_at_d, Tf - T0r, dx, nprof);
                    if (std::isfinite(ell_p))
                    {
                        Set::Scalar hf = sp_flake_c_f * ell_p;
                        const Set::Scalar h_max = spall_h_spall_max / L.mu;
                        const Set::Scalar h_bnd = sp_flake_h_max / L.mu;
                        if (h_max > 0.0 && hf > h_max) { hf = h_max; o.rec.h_f_capped = true; }
                        if (h_bnd > 0.0 && hf > h_bnd) { hf = h_bnd; o.rec.h_f_bounded = true; }
                        o.rec.h_f = hf;
                        o.rec.ell_eff = ell_p;
                    }
                    else
                    {
                        o.rec.h_f = 0.0;
                        o.nodepth = 1;
                    }
                }
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
                        const Set::Scalar T0_z =
                            T0_at_d(static_cast<int>(std::floor(z_v / dx)));
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
                // ⚠⚠ P2-24: NO FRACTURE IN OR OVER A LIQUID. A line is not scanned when its surface cell is fully
                // liquid (fr_n_molten) or ANY cell inside its window is (a solid lid over a melt: fr_n_lid) -- the
                // K_I integral must not run through a liquid. f_L is the beam's own (MeltFraction). A partly
                // molten cell (0 < f_L < 1) still carries K_Ic (declared). The surface reconstruction above has
                // run, so the line's T_s still carries to the next step.
                if (MeltFraction(Tc) >= 1.0) { ++n_molten; continue; }
                {
                    bool lid = false;
                    for (int d = 1; d <= nprof && !lid; ++d) lid = (MeltFraction(T_at_d(d)) >= 1.0);
                    if (lid) { ++n_lid; continue; }
                }
                int s_last = -1;
                bool failed = false;
                Set::Scalar score = 0.0;
                if (std::isfinite(a_f) && a_f > 0.0)
                    o.h_col = DepthScanH(nprof, 0, -(L.Lg + 1), dx, 1.0, 1.0, a_f, true,
                                         k_i_at_zcrack, kic_pa_at_z, kic_pa_at_cell,
                                         [] (int) {}, s_last, failed, score);
                o.score = score;
                o.fire = (o.h_col > 0.0) ? 1 : 0;
                // the crack reaches the void behind the line (a sheet, or any
                // line whose run ends inside the window): no rock to stop it
                if (!o.fire && score >= 1.0 && !failed && s_last == nprof
                    && nprof == run - 1 && run - 1 < L.Dmax)
                { o.fire = 1; o.through = 1; }
                // P2-19: DIAGNOSTICS on every fired line, main branch included -- pure
                // reads, consumed by nothing (the event log only): ell_p, the e-fold
                // depth of the rise from this line's T_s (on a fallback line it is the
                // depth just used; on a main-branch line it shows the thickness a
                // switch of branch would jump to); the peak's depth and the depth of
                // its own e-fold below it; and the ENERGY THIN LIMIT
                // h_E = K_Ic(T_s)^2 (1 + nu) / sigma_s^2 [m], the CIRCULAR-blister form
                // (a straight-sided blister gives 2 K^2 / sigma^2, 1.5x larger),
                // sigma_s the relieved criterion stress at the face [Pa] (ScanSigmaT
                // applies sp_relief), K_Ic in Pa sqrt(m). Logged, never applied.
                if (o.fire)
                {
                    o.ell_p = ProfileEFoldDepth(rise_at_d, Ts - T0r, dx, nprof);
                    ProfilePeak(rise_at_d, Ts - T0r, dx, nprof, o.z_peak, o.ell_pk);
                    const Set::Scalar sig_s = ScanSigmaT(scan_g, T_at_zv(0.0), T0_at_d(0), 1.0,
                                                         conf_off, pE, pB, pNu);
                    const Set::Scalar kic_s = ScanKIcPa(Ts);
                    o.h_E = (sig_s != 0.0) ? kic_s * kic_s * (1.0 + pNu) / (sig_s * sig_s)
                                           : std::numeric_limits<Set::Scalar>::infinity();
                }
                // ⚠⚠ P2-15: THE SECOND CONDITION. A line that fired is withdrawn unless its
                // flake can buckle over the span: bk = C_L x / L_span <= 1. x is the flake's
                // NORMAL thickness: a through line's whole solid run (L_av, the same
                // fractions R3 reads), else c_f ell_eff -- both measured ALONG the line, so
                // x mu (ell_eff is the along-line length ell/mu, R1). Withdrawn = no fire, no
                // through; the onset criterion above is unchanged.
                if (bk_on && o.fire)
                {
                    Set::Scalar x;
                    if (o.through)
                    {
                        Set::Scalar L_av = 0.0;
                        for (int d = 0; d < o.run; ++d) L_av += Gl[d * NF + 4] * dx;
                        x = L.mu * L_av;
                    }
                    else
                    {
                        const Set::Scalar ell = o.rec.ell_eff;
                        x = (std::isfinite(ell) && ell > 0.0)
                            ? L.mu * sp_flake_c_f * ell
                            : std::numeric_limits<Set::Scalar>::infinity();
                    }
                    o.bk = BuckleRatio(scan_g, pNu, pB, o.rec.T_surf, T0_at_d(0), x, sp_buckle_span);
                    // P2-23: the record, taken BEFORE the withdrawal below. bk_lb is the
                    // LOWER BOUND B would read if the layer were the scan window itself
                    // (ell = nprof dx, along the line) -- the same helper, T_s, T0, nu, mu
                    // and c_f, x replaced; LOGGED only (x, bk and the withdrawal untouched).
                    if (bk_log)
                    {
                        BkRow br;
                        br.l = l; br.through = o.through; br.x = x; br.bk = o.bk;
                        br.bk_lb = BuckleRatio(scan_g, pNu, pB, o.rec.T_surf, T0_at_d(0),
                                               L.mu * sp_flake_c_f * (Set::Scalar(nprof) * dx),
                                               sp_buckle_span);
                        br.blocked = !(o.bk <= 1.0) ? 1 : 0;
                        bk_rows.push_back(br);
                    }
                    ++bk_cnt[0];
                    if (!(o.bk <= 1.0))
                    {
                        o.fire = 0;
                        o.through = 0;
                        o.bk_blocked = 1;
                        ++bk_cnt[1];
                        if (!std::isfinite(x)) ++bk_cnt[2];
                        ++bk_cnt[3 + FrClass(L.dir)];
                    }
                    else if (o.bk > bk_pass_max) bk_pass_max = o.bk;
                }
            }

            // P2-15: the lines' share, added once (not rank-summed)
            if (bk_on)
            {
                bk_n_eval        += static_cast<Set::Scalar>(bk_cnt[0]);
                bk_n_block       += static_cast<Set::Scalar>(bk_cnt[1]);
                bk_n_inf         += static_cast<Set::Scalar>(bk_cnt[2]);
                bk_n_block_floor += static_cast<Set::Scalar>(bk_cnt[3]);
                bk_n_block_wall  += static_cast<Set::Scalar>(bk_cnt[4]);
                bk_n_block_roof  += static_cast<Set::Scalar>(bk_cnt[5]);
                bk_max_pass       = std::max(bk_max_pass, bk_pass_max);
            }
            // P2-23: the blocked-line record. `time` is the step's START (as the event
            // log); its thermo row is stamped time + dt. class 0 floor / 1 wall / 2 roof;
            // fire 1 on every row (B is evaluated only on a fired line); through, fb,
            // nodepth the line's own flags before the withdrawal; lengths in m along the
            // line except x (the NORMAL flake thickness B reads, inf on a no-depth line);
            // bk = C_L x / L_span; bk_lb its lower bound at the window; blocked = !(bk <= 1).
            if (bk_log && amrex::ParallelDescriptor::IOProcessor())
            {
                const std::string csv_path = plot_file + "_buckle_log.csv";
                const bool need_header = !sp_buckle_log_header_written;
                std::ofstream ofs(csv_path, (need_header ? std::ios::trunc : std::ios::app)
                                            | std::ios::out);
                if (ofs.is_open())
                {
                    if (need_header)
                        ofs << "time,col_i,col_j,class,fire,through,fb,nodepth,score,h_col,"
                               "T_surf,ell_p,z_peak,ell_pk,ell_eff,x,bk,bk_lb,blocked\n";
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    for (const BkRow& br : bk_rows)
                    {
                        const FrLine& L = lines[br.l];
                        const LineOut& o = out[br.l];
                        const FrCell& c = cells[static_cast<std::size_t>(L.cell)];
                        ofs << time << "," << c.i << "," << c.j << "," << FrClass(L.dir)
                            << ",1," << br.through << "," << o.fb << "," << o.nodepth
                            << "," << o.score << "," << o.h_col << "," << o.rec.T_surf
                            << "," << o.ell_p << "," << o.z_peak << "," << o.ell_pk
                            << "," << o.rec.ell_eff << "," << br.x << "," << br.bk
                            << "," << br.bk_lb << "," << br.blocked << "\n";
                    }
                }
            }
            if (bk_log) sp_buckle_log_header_written = true;
            // P2-24: the exclusions' counts (every rank evaluates every line: not rank-summed)
            fr_n_molten = static_cast<Set::Scalar>(n_molten);
            fr_n_lid    = static_cast<Set::Scalar>(n_lid);

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
                if (o.nodepth) fr_nodepth_c[cls] += 1.0;
                // no depth (the rise never e-folds) and no crack to the void:
                // nothing to recede (P2-19: was "no heated skin", q_net <= 0)
                if (o.fire && o.nodepth && !o.through)
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

            // P2-24: the vapour rule's books -- the cells that leave, the enthalpy they carry (H V, booked ONCE:
            // it came in through deposition and conduction), and whether each is a column top (a +z face listed:
            // the column's phi moves down by the cell's remaining solid thickness, as a +z line's recession does)
            std::vector<std::size_t> vap_list;
            std::vector<char> vap_top(cells.size(), 0);
            for (std::size_t ci = 0; ci < cells.size(); ++ci)
            {
                if (!vap_cell[ci]) continue;
                vap_list.push_back(ci);
                ledger_E_vapor += vap_H[ci] * Vcell;
                for (int q = cells[ci].q0; q < cells[ci].q0 + cells[ci].nq; ++q)
                    if (fl_faces[static_cast<std::size_t>(q)].dir == 5) vap_top[ci] = 1;
            }
            fr_n_vapor = static_cast<Set::Scalar>(vap_list.size());
            bc_n_vapor += fr_n_vapor;

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
                // P2-24: the vapour voids (regime 2); the thermal state stays frozen, as a spalled void's
                if (!vap_list.empty())
                {
                    amrex::Array4<Set::Scalar> const& vev  = (*vapor_event_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& vthk = (*vapor_thickness_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& vrop = (*RoP_vap_mf[lev]).array(mfi);
                    for (const std::size_t ci : vap_list)
                    {
                        const FrCell& c = cells[ci];
                        const Set::Scalar h_v = vap_f[ci] * dx;
                        if (vap_top[ci] && c.i >= bx.smallEnd(0) && c.i <= bx.bigEnd(0)
                            && c.j >= bx.smallEnd(1) && c.j <= bx.bigEnd(1))
                            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
                                phi(c.i, c.j, k) -= h_v;
                        if (!bx.contains(amrex::IntVect(AMREX_D_DECL(c.i, c.j, c.k)))) continue;
                        ff(c.i, c.j, c.k) = 0.0;
                        rm(c.i, c.j, c.k) = 1.0;
                        reg(c.i, c.j, c.k) = 2.0;
                        regf(c.i, c.j, c.k) = 2.0;
                        ropU(c.i, c.j, c.k) = (dt > 0.0) ? h_v / dt : 0.0;
                        D(c.i, c.j, c.k) = 0.0;
                        vev(c.i, c.j, c.k) = 1.0;
                        vthk(c.i, c.j, c.k) = h_v;
                        vrop(c.i, c.j, c.k) = (dt > 0.0) ? h_v / dt : 0.0;
                        any_write = true;
                    }
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
                            // P2-15: under spall.buckle_span > 0 only, the line's ratio
                            << (bk_on ? ",bk" : "")
                            // P2-19, appended: the line's profile depths and the energy
                            // thin limit [m] (diagnostics; nothing reads them)
                            << ",ell_p,z_peak,ell_pk,h_E"
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
                        if (bk_on) ofs << "," << o.bk;
                        ofs << "," << o.ell_p << "," << o.z_peak << "," << o.ell_pk
                            << "," << o.h_E;
                        ofs << "\n";
                    }
                    // P2-24: one row per vaporised cell, regime 2; h_applied = the cell's remaining solid
                    // thickness; the line columns not defined for a vapour cell read 0 / NaN
                    const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
                    for (const std::size_t ci : vap_list)
                    {
                        const FrCell& c = cells[ci];
                        const Set::Scalar h_v = vap_f[ci] * dx;
                        ofs << time << "," << c.i << "," << c.j << ",2,"
                            << k_top_global[static_cast<std::size_t>(c.col)] << ","
                            << vap_T[ci] << ",0,0,0," << h_v << ",1"
                            << ",0,0,0,0,0,0,0,0,0";
                        ofs << "," << fl_faces[static_cast<std::size_t>(c.qline)].dir << ",1," << c.k << ",0,"
                            << h_v << "," << h_v << ",0,0,0," << c.sheet << ",0,0,0";
                        if (bk_on) ofs << ",0";
                        ofs << "," << nan << "," << nan << "," << nan << "," << nan << "\n";
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
