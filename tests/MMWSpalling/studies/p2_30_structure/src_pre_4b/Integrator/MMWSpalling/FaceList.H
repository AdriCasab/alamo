    // The solid fraction at first use: 1 on solid, 0 on void (after any seed).
    void FrInitFraction(int lev)
    {
        if (fr_f_init) return;
        for (amrex::MFIter mfi(*fr_f_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar> const& ff = (*fr_f_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& rm = (*removed_mf[lev]).const_array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                ff(i, j, k) = (rm(i, j, k) > 0.5) ? 0.0 : 1.0; });
        }
        fr_f_init = true;
    }

    // The fraction's ghost rule: box to box and periodic by FillBoundary; out
    // of the domain, ABOVE the top is gas (0); a SYMMETRY PLANE (the patch
    // centre lies on that side) is mirrored; every other side takes the
    // nearest in-domain cell.
    void FillFractionGhosts(int lev)
    {
        amrex::MultiFab& Fm = *fr_f_mf[lev];
        Fm.FillBoundary(geom[lev].periodicity());
        const amrex::Box& dom = geom[lev].Domain();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar* PHI = geom[lev].ProbHi();
        const Set::Scalar tol = 1.0e-9 * (PHI[0] - PLO[0]);
        const bool mir_lo[3] = {std::abs(surface_patch.x0 - PLO[0]) <= tol,
                                std::abs(surface_patch.y0 - PLO[1]) <= tol, false};
        const bool mir_hi[3] = {std::abs(surface_patch.x0 - PHI[0]) <= tol,
                                std::abs(surface_patch.y0 - PHI[1]) <= tol, false};
        const bool per[3] = {geom[lev].isPeriodic(0), geom[lev].isPeriodic(1),
                             geom[lev].isPeriodic(2)};
        for (amrex::MFIter mfi(Fm, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box gb = mfi.fabbox();
            amrex::Array4<Set::Scalar> const& a = Fm.array(mfi);
            const amrex::Dim3 lo = amrex::lbound(gb), hi = amrex::ubound(gb);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int p[3] = {i, j, k};
                bool out = false, gas = false;
                for (int d = 0; d < 3; ++d)
                {
                    if (per[d]) continue;
                    const int dl = dom.smallEnd(d), dh = dom.bigEnd(d);
                    if (p[d] < dl) { out = true; p[d] = mir_lo[d] ? 2 * dl - 1 - p[d] : dl; }
                    else if (p[d] > dh)
                    {
                        out = true;
                        if (d == 2) gas = true;
                        p[d] = mir_hi[d] ? 2 * dh + 1 - p[d] : dh;
                    }
                }
                if (!out) continue;
                if (gas) { a(i, j, k) = 0.0; continue; }
                if (!gb.contains(amrex::IntVect(AMREX_D_DECL(p[0], p[1], p[2]))))
                    Util::Abort(INFO, "FillFractionGhosts: a ghost's source cell (", p[0], ", ",
                                p[1], ", ", p[2], ") is outside its own fab -- a box narrower "
                                "than 3 cells at a domain side");
                a(i, j, k) = a(p[0], p[1], p[2]);
            }
        }
    }

    // F3 -- THE AREA RULE. The one weight function.
    static Set::Scalar FaceWeight(const Set::Vector& g, Set::Scalar& nx,
                                  Set::Scalar& ny, Set::Scalar& nz, int& fb)
    {
        const Set::Scalar gn = std::sqrt(g(0) * g(0) + g(1) * g(1) + g(2) * g(2));
        if (!(gn > 0.0)) { nx = ny = nz = 0.0; fb = 1; return 1.0; }
        nx = -g(0) / gn; ny = -g(1) / gn; nz = -g(2) / gn;
        fb = 0;
        return 1.0 / (std::abs(nx) + std::abs(ny) + std::abs(nz));
    }

    // F2 -- THE HEAT RULE. The one heat resolver: h at a face centre.
    Set::Scalar FaceH(Set::Scalar x, Set::Scalar y, Set::Scalar s, Set::Scalar time)
    {
        const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
        const Set::Scalar r2 = dxp * dxp + dyp * dyp;
        const Set::Scalar r = std::sqrt(r2);
        return FaceHr(x, y, r, s, time);
    }
    // P2-11b: the same with r given -- face_heat = path passes r := s_path,
    // s := H.
    Set::Scalar FaceHr(Set::Scalar x, Set::Scalar y, Set::Scalar r, Set::Scalar s,
                       Set::Scalar time)
    {
        Set::Scalar h = 0.0;
        if (surface_patch.jet_enthalpy)
            h = (s > 0.0) ? surface_patch.h_f(x, y, r, s, time) : Set::Scalar(0.0);
        else
            h = surface_patch.h_f(x, y, r, s, time);
        if (!(std::isfinite(h) && h >= 0.0))
            Util::Abort(INFO, "surface_patch.face_list: h_expr must be finite and >= 0 "
                        "at every face centre; got h = ", h, " at x = ", x, ", y = ", y,
                        ", s = ", s);
        return h;
    }

    // F1 + F3 -- THE LIST BUILDER. Once per step, before BuildFlameColumns.
    void BuildFaceList(int lev, const amrex::Vector<int>& k_top_col_in)
    {
        fl_ready = false;
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1), klo = domain.smallEnd(2);
        const int ihi = domain.bigEnd(0),   jhi = domain.bigEnd(1),   khi = domain.bigEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        // cubic cells: checked once at Parse (CheckInterlocks, P2-27)
        const amrex::Vector<int> k_top_col =
            k_top_col_in.empty() ? ColumnTopSolid(lev) : k_top_col_in;
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;

        // P2-9 rev 2 item 4: F3's field is the SOLID FRACTION (3 ghosts, the
        // fraction's own ghost rule), and its normal is the smoothed one.
        const amrex::BoxArray& ba = removed_mf[lev]->boxArray();
        const amrex::DistributionMapping& dm = removed_mf[lev]->DistributionMap();
        FrInitFraction(lev);
        FillFractionGhosts(lev);

        // the local faces, packed FL_NF doubles each
        std::vector<Set::Scalar> loc;
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& Sa  = (*fr_f_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp  = (*kappa_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (rem(i, j, k) > 0.5) continue;
                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (dxp * dxp + dyp * dyp > r2max) continue;
                bool open[6];
                open[0] = (i > ilo) && rem(i - 1, j, k) > 0.5;
                open[1] = (i < ihi) && rem(i + 1, j, k) > 0.5;
                open[2] = (j > jlo) && rem(i, j - 1, k) > 0.5;
                open[3] = (j < jhi) && rem(i, j + 1, k) > 0.5;
                open[4] = (k > klo) && rem(i, j, k - 1) > 0.5;
                open[5] = (k == khi) || rem(i, j, k + 1) > 0.5;
                if (!(open[0] || open[1] || open[2] || open[3] || open[4] || open[5]))
                    continue;
                const int col = (j - jlo) * nx + (i - ilo);
                // P2-9 rev 2: THE normal (FrSmoothedGradient), the one function
                // every consumer's n and w come from
                const Set::Vector g = FrSmoothedGradient(Sa, i, j, k, DX[0]);
                Set::Scalar fnx = 0.0, fny = 0.0, fnz = 0.0;
                int fb = 0;
                const Set::Scalar w = FaceWeight(g, fnx, fny, fnz, fb);
                // the conductivity exactly as the kernel forms it (kappa_phase;
                // the damage factor exp(-a_D D) is retired, P2-27)
                const Set::Scalar K = kp(i, j, k);
                for (int d = 0; d < 6; ++d)
                {
                    if (!open[d]) continue;
                    const int cls = (d == 5) ? ((k == k_top_col[static_cast<std::size_t>(col)])
                                                ? FL_TOP : FL_UP)
                                  : (d == 4) ? FL_DOWN : FL_SIDE;
                    const Set::Scalar pk[FL_NF] = {
                        Set::Scalar(i), Set::Scalar(j), Set::Scalar(k), Set::Scalar(d),
                        Set::Scalar(cls), Set::Scalar(col), Set::Scalar(fb), w,
                        fnx, fny, fnz, T_old(i, j, k), K};
                    loc.insert(loc.end(), pk, pk + FL_NF);
                }
            }
        }
        // gather: disjoint rank segments, summed -- identical on every rank
        const int np = amrex::ParallelDescriptor::NProcs();
        const int me = amrex::ParallelDescriptor::MyProc();
        std::vector<int> cnt(static_cast<std::size_t>(np), 0);
        cnt[static_cast<std::size_t>(me)] = static_cast<int>(loc.size() / FL_NF);
        amrex::ParallelDescriptor::ReduceIntSum(cnt.data(), np);
        std::size_t off = 0, tot = 0;
        for (int p = 0; p < np; ++p)
        {
            if (p == me) off = tot;
            tot += static_cast<std::size_t>(cnt[static_cast<std::size_t>(p)]);
        }
        std::vector<Set::Scalar> all(tot * FL_NF, 0.0);
        std::copy(loc.begin(), loc.end(), all.begin() + static_cast<std::ptrdiff_t>(off * FL_NF));
        if (tot > 0)
            amrex::ParallelDescriptor::ReduceRealSum(all.data(), static_cast<int>(tot * FL_NF));
        fl_faces.assign(tot, FaceEntry());
        for (std::size_t q = 0; q < tot; ++q)
        {
            const Set::Scalar* p = all.data() + q * FL_NF;
            FaceEntry& e = fl_faces[q];
            e.i = static_cast<int>(p[0]); e.j = static_cast<int>(p[1]);
            e.k = static_cast<int>(p[2]); e.dir = static_cast<int>(p[3]);
            e.cls = static_cast<int>(p[4]); e.col = static_cast<int>(p[5]);
            e.fb = static_cast<int>(p[6]); e.w = p[7];
            e.nx = p[8]; e.ny = p[9]; e.nz = p[10]; e.T = p[11]; e.K = p[12];
        }
        std::sort(fl_faces.begin(), fl_faces.end(),
                  [](const FaceEntry& a, const FaceEntry& b) {
                      if (a.col != b.col) return a.col < b.col;
                      if (a.k != b.k) return a.k < b.k;
                      return a.dir < b.dir; });
        // per column: the entry range
        fl_col_first.assign(ncols, 0);
        fl_col_count.assign(ncols, 0);
        Set::Scalar n_cls[4] = {0.0, 0.0, 0.0, 0.0}, n_fb = 0.0;
        for (std::size_t q = 0; q < fl_faces.size(); ++q)
        {
            const FaceEntry& e = fl_faces[q];
            const std::size_t c = static_cast<std::size_t>(e.col);
            if (fl_col_count[c] == 0) fl_col_first[c] = static_cast<int>(q);
            ++fl_col_count[c];
            n_cls[e.cls] += 1.0;
            if (e.fb) n_fb += 1.0;
        }
        fl_n_top = n_cls[FL_TOP]; fl_n_up = n_cls[FL_UP];
        fl_n_side = n_cls[FL_SIDE]; fl_n_down = n_cls[FL_DOWN];
        fl_w_fallback = n_fb;
        // per cell: the entry range (a cell's faces are contiguous by the sort)
        if (!fl_index || fl_index->boxArray() != ba || fl_index->DistributionMap() != dm)
            fl_index = std::make_unique<amrex::iMultiFab>(ba, dm, 2, 0);
        fl_index->setVal(0);
        if (fl_w_mf[lev]) fl_w_mf[lev]->setVal(0.0);
        for (amrex::MFIter mfi(*fl_index, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<int> const& fi = fl_index->array(mfi);
            amrex::Array4<Set::Scalar> fw;
            if (fl_w_mf[lev]) fw = (*fl_w_mf[lev]).array(mfi);
            for (std::size_t q = 0; q < fl_faces.size(); ++q)
            {
                const FaceEntry& e = fl_faces[q];
                if (!bx.contains(amrex::IntVect(AMREX_D_DECL(e.i, e.j, e.k)))) continue;
                if (fi(e.i, e.j, e.k, 1) == 0) fi(e.i, e.j, e.k, 0) = static_cast<int>(q);
                fi(e.i, e.j, e.k, 1) += 1;
                if (fw) fw(e.i, e.j, e.k) += e.w;
            }
        }
        fl_ready = true;
    }

    // F2 over the list. Called by BuildFlameColumns once z_n is known and
    // BEFORE the march, which bills every face at this h.
    void FaceListHeat(int lev, Set::Scalar z_n, Set::Scalar time,
                      const amrex::Vector<int>& k_top_col)
    {
        if (!fl_ready) return;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        long mism = 0;
        // P2-23: no flame (the patch off): every listed face is billed at h = 0
        // -- its radiative and convective losses only. No FaceH, and no
        // flame_h to compare against (it is empty).
        if (!surface_patch.enabled)
        {
            for (auto& e : fl_faces) e.h = 0.0;
        }
        // P2-28: a PRESCRIBED patch -- no flame: every face h = 0, and a TOP face (the
        // column top's +z face) whose column centre lies inside the patch radius holds
        // T_expr(t), the column term's own test (x, y at the cell centre, r^2 <= R^2).
        // Side and roof faces are never prescribed: the patch is a flat top boundary.
        else if (surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature)
        {
            const Set::Scalar Tp = surface_patch.T_at(time);
            const Set::Scalar r2 = surface_patch.radius * surface_patch.radius;
            for (auto& e : fl_faces)
            {
                e.h = 0.0;
                e.T_presc = -1.0;
                if (e.cls != FL_TOP) continue;
                const Set::Scalar x = PLO[0] + (Set::Scalar(e.i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(e.j) + 0.5) * DX[1];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (dxp * dxp + dyp * dyp <= r2) e.T_presc = Tp;
            }
            for (const auto& e : fl_faces)
                if (e.cls == FL_TOP) flame_h[static_cast<std::size_t>(e.col)] = e.h;
        }
        // P2-11b: face_heat = path -- every face from the surface path; the
        // column's flame_h (read only by column diagnostics under the list)
        // takes its top face's value, so fl_h_top_mismatch stays an identity.
        else if (surface_patch.face_heat_path)
        {
            FaceListPath(lev, time, k_top_col);
            for (const auto& e : fl_faces)
                if (e.cls == FL_TOP) flame_h[static_cast<std::size_t>(e.col)] = e.h;
        }
        else
        for (auto& e : fl_faces)
        {
            Set::Scalar ox = 0.5, oy = 0.5;
            if (e.dir == 0) ox = 0.0; else if (e.dir == 1) ox = 1.0;
            if (e.dir == 2) oy = 0.0; else if (e.dir == 3) oy = 1.0;
            // ⚠ the column-centre and top-face expressions are BuildFlameColumns'
            // own, token for token, so the top face reproduces flame_h exactly.
            const Set::Scalar x = (ox == 0.5) ? PLO[0] + (Set::Scalar(e.i) + 0.5) * DX[0]
                                              : PLO[0] + (Set::Scalar(e.i) + ox) * DX[0];
            const Set::Scalar y = (oy == 0.5) ? PLO[1] + (Set::Scalar(e.j) + 0.5) * DX[1]
                                              : PLO[1] + (Set::Scalar(e.j) + oy) * DX[1];
            Set::Scalar z_f;
            if (e.dir == 5)      z_f = PLO[2] + Set::Scalar(e.k + 1) * DX[2];
            else if (e.dir == 4) z_f = PLO[2] + Set::Scalar(e.k) * DX[2];
            else                 z_f = PLO[2] + (Set::Scalar(e.k) + 0.5) * DX[2];
            e.h = FaceH(x, y, z_n - z_f, time);
            if (e.cls == FL_TOP
                && !(e.h == flame_h[static_cast<std::size_t>(e.col)])) ++mism;
        }
        fl_h_top_mismatch = static_cast<Set::Scalar>(mism);
        if (fl_h_mf[lev])
        {
            fl_h_mf[lev]->setVal(0.0);
            for (amrex::MFIter mfi(*fl_index, false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const int> const& fi = fl_index->const_array(mfi);
                amrex::Array4<Set::Scalar> const& fh = (*fl_h_mf[lev]).array(mfi);
                const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                    for (int q = fi(i, j, k, 0); q < fi(i, j, k, 0) + fi(i, j, k, 1); ++q)
                        fh(i, j, k) = std::max(fh(i, j, k), fl_faces[static_cast<std::size_t>(q)].h);
            }
        }
    }

    // The per-class robin power, from the replicated list with the kernel's
    // own inputs (T_old, kc, each face's h and gas temperature) -- the kernel's
    // deposit term by term, summed by class. Every rank computes the same.
    void FaceListDiagnostics(int lev)
    {
        fl_P_top = fl_P_up = fl_P_side = fl_P_down = 0.0;
        if (!fl_ready) return;
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar A = DX[0] * DX[1];
        Set::Scalar P[4] = {0.0, 0.0, 0.0, 0.0};
        for (const auto& e : fl_faces)
        {
            Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
            // P2-28: a PRESCRIBED face -- the kernel's half-cell Dirichlet flux 2 K (T_presc - T) / dz
            if (e.T_presc > 0.0)
            {
                P[e.cls] += e.w * (2.0 * e.K * (e.T_presc - e.T) / DX[2]) * A;
                continue;
            }
            // P2-9: the kernel's own term, FrFaceFlux, for every class
            // (P2-11b: at the face's own gas temperature)
            FrFaceFlux(e, e.Tg, e.K, e.T, DX[0], q_in, q_out, q_rob, q_los);
            P[e.cls] += e.w * q_rob * A;
        }
        fl_P_top = P[FL_TOP]; fl_P_up = P[FL_UP]; fl_P_side = P[FL_SIDE]; fl_P_down = P[FL_DOWN];
        const Set::Scalar Ps = fl_P_top + fl_P_up + fl_P_side + fl_P_down;
        // P2-23: the march's identity exists only with a march (the patch on);
        // flameless, jet_P_face is never set and the identity is vacuous.
        // P2-28: a prescribed patch has no flame and no march either
        const bool presc = surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        fl_id_march = (surface_patch.enabled && !presc && std::abs(Ps) > 0.0)
                    ? std::abs(jet_P_face - Ps) / std::abs(Ps) : 0.0;
    }

    // S1: per-cell surface heat flux [W/m^2] for the microstructure path,
    // shared by the H update and EnergyLedger so both use identical
    // arithmetic. Called for every non-void cell.
    //   q_in    added to the cell as q_in/dz   (nonzero only when in_patch)
    //   q_out   removed from the cell as q_out/dz
    //   q_robin, q_loss   the ledger split, with q_robin - q_loss == q_in - q_out
    //
    // Cell form (default, and all cells outside the patch):
    //   q_in = q_robin = h*(T_gas - T_c)   [in_patch]
    //   q_out = q_loss = surf*(eps(T_c)*sigma*(T_c^4 - T_a^4) + h_c*(T_c - T_a))
    // (byte-identical to the pre-S1 inline expressions).
    //
    // Face form (in_patch && face): the top solid cell centre is dz/2 below
    // its top face. Solve
    //   (2k_c/dz)(T_f - T_c) = h(T_gas - T_f) - eps(T_f)*sigma*(T_f^4 - T_a^4) - h_c(T_f - T_a)
    // for T_f with safeguarded Newton (start T_c, <= 10 iterations, stop when
    // |dT| <= 1e-9*T_c). The residual is monotone increasing in T_f, so the
    // root is bracketed by [min, max](T_c, T_gas, T_a); T_a is needed in the
    // bracket once losses can pull T_f below T_c. A Newton step leaving the
    // bracket is replaced by bisection. The cell receives the net face flux
    // q_in = (2k_c/dz)(T_f - T_c) and no cell-centred loss (q_out = 0);
    // q_robin = h(T_gas - T_f), and q_loss = q_robin - q_in (= the face
    // losses at T_f to Newton tolerance, but the ledger closes to round-off).
    //
    static void SurfaceCellFlux(bool in_patch, bool face_form, Set::Scalar surf_c,
                                Set::Scalar h, Set::Scalar T_gas,
                                Set::Scalar kc, Set::Scalar dz, Set::Scalar Tc,
                                Set::Scalar eps_lo, Set::Scalar eps_hi, Set::Scalar T_eps_b,
                                Set::Scalar sig_SB, Set::Scalar h_c, Set::Scalar T_a,
                                Set::Scalar& q_in, Set::Scalar& q_out,
                                Set::Scalar& q_robin, Set::Scalar& q_loss,
                                Set::Scalar T_surf_in = 0.0,
                                Set::Scalar* T_face_out = nullptr)
    {
        // T_face_out (P1f deliverable 0): the temperature of the face the Robin
        // flux was evaluated at, read by the face list (FaceRemoval.H). `T_surf_in`
        // is an INPUT and the face form's `Tf` is otherwise a local, so this is
        // the only way out for it.
        //
        // ⚠⚠ EACH BRANCH WRITES THE TEMPERATURE ITS OWN q_robin USED, and the CELL
        // FORM WRITES NOTHING. There q_robin = h*(T_gas - Tc) is built on the CELL
        // AVERAGE, which is not a face temperature; writing Tc would hand the
        // reader a plausible number that is the wrong quantity.
        //
        // ⚠ NOT the flux inversion T_gas - q_robin/h. That is algebraically exact
        // and is (k3)'s cross-check, but it divides by h, which h_expr drives to
        // 1e-3 W/m^2K in the 60->76 mm taper, and its meaning changes with branch.
        // P0b deliverable 1: robin_form = surface. The flux AND the losses are
        // evaluated at the reconstructed surface temperature, including the
        // eps_lo/eps_hi selection at T_eps_b -- the packet is explicit that the
        // losses must not be left on T_c.
        //   q_robin = h (T_gas - T_surf)
        //   q_loss  = eps(T_surf) sig (T_surf^4 - T_a^4) + h_c (T_surf - T_a)
        // ✅ No Newton: the half-cell solve exists to INFER a face temperature,
        // and with the reconstruction in hand there is nothing to infer. The
        // `face` form's G = 2k/dz is a LINEAR half-cell profile, which is
        // exactly where it disagrees with the exponential reconstruction --
        // at these rates ell_eff sits below dz/2.
        // ⚠ A caller that passes nothing gets T_surf_in = 0, so `cell` and
        // `face` stay byte-identical (the face list passes its cell's T_s).
        if (in_patch && T_surf_in > 0.0)
        {
            const Set::Scalar Ts    = T_surf_in;
            const Set::Scalar eps_s = (Ts <= T_eps_b) ? eps_lo : eps_hi;
            q_robin = h * (T_gas - Ts);
            q_loss  = eps_s * sig_SB * (Ts * Ts * Ts * Ts - T_a * T_a * T_a * T_a)
                    + h_c * (Ts - T_a);
            q_in    = q_robin - q_loss;
            q_out   = 0.0;
            if (T_face_out) *T_face_out = Ts;   // P1f: the reconstructed surface T
            return;
        }
        if (in_patch && face_form)
        {
            const Set::Scalar G = 2.0 * kc / dz;
            Set::Scalar lo = std::min(std::min(Tc, T_gas), T_a);
            Set::Scalar hi = std::max(std::max(Tc, T_gas), T_a);
            Set::Scalar Tf = Tc;
            for (int it = 0; it < 10; ++it)
            {
                const Set::Scalar eps = (Tf <= T_eps_b) ? eps_lo : eps_hi;
                const Set::Scalar T3  = Tf * Tf * Tf;
                const Set::Scalar F   = G * (Tf - Tc) - h * (T_gas - Tf)
                                      + eps * sig_SB * (T3 * Tf - T_a * T_a * T_a * T_a)
                                      + h_c * (Tf - T_a);
                if (F == 0.0) break;
                if (F > 0.0) hi = Tf; else lo = Tf;
                const Set::Scalar dF = G + h + 4.0 * eps * sig_SB * T3 + h_c;
                Set::Scalar Tn = Tf - F / dF;
                if (!(Tn >= lo && Tn <= hi)) Tn = 0.5 * (lo + hi);
                const Set::Scalar step = Tn - Tf;
                Tf = Tn;
                if (std::abs(step) <= 1.0e-9 * Tc) break;
            }
            q_in    = G * (Tf - Tc);
            q_out   = 0.0;
            q_robin = h * (T_gas - Tf);
            q_loss  = q_robin - q_in;
            // P1f: the Newton's own root -- the temperature this q_robin used.
            // ⚠ Inferred on a LINEAR half-cell (G = 2k/dz) while the criterion
            // reads an exponential; the two disagree once ell_lat < dz/2, which
            // criterion (k2) scores rather than assumes away.
            if (T_face_out) *T_face_out = Tf;
            return;
        }
        // ⚠⚠ P1f: THE CELL FORM DELIBERATELY WRITES NOTHING TO T_face_out.
        // q_robin below is h*(T_gas - Tc) on the CELL AVERAGE; Tc is not a face
        // temperature and never was. A caller that needs a face temperature must
        // see its own initial value survive, not a plausible substitute. (a1).
        const Set::Scalar eps = (Tc <= T_eps_b) ? eps_lo : eps_hi;
        q_in    = in_patch ? h * (T_gas - Tc) : Set::Scalar(0.0);
        q_out   = surf_c * (eps * sig_SB * (Tc*Tc*Tc*Tc - T_a*T_a*T_a*T_a)
                            + h_c * (Tc - T_a));
        q_robin = q_in;
        q_loss  = q_out;
    }
