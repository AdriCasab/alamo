    // P2-11b item 5: every listed face's gas temperature, after the march.
    // position: the column's value, exactly. path: the march wrote e.Tg per
    // face; a face it did not visit (blocked, unreached, s <= 0) takes T_ent,
    // and the column's flame_Tg (column diagnostics only under the list) takes
    // its top face's value. Then the path-mode plotfile fields.
    void FaceListTg(int lev)
    {
        if (!fl_ready) return;
        // P2-23: no flame (the patch off): no gas -- and flame_Tg is EMPTY, so
        // nothing is indexed
        if (!surface_patch.enabled)
        {
            for (auto& e : fl_faces) e.Tg = 0.0;
            return;
        }
        if (!surface_patch.face_heat_path)
        {
            for (auto& e : fl_faces) e.Tg = flame_Tg[static_cast<std::size_t>(e.col)];
            return;
        }
        for (const auto& e : fl_faces)
            if (e.cls == FL_TOP) flame_Tg[static_cast<std::size_t>(e.col)] = e.Tg;
        if (fl_s_mf[lev] && fl_Tg_mf[lev])
        {
            fl_s_mf[lev]->setVal(-1.0);
            fl_Tg_mf[lev]->setVal(0.0);
            for (amrex::MFIter mfi(*fl_index, false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const int> const& fi = fl_index->const_array(mfi);
                amrex::Array4<Set::Scalar> const& fs = (*fl_s_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& ft = (*fl_Tg_mf[lev]).array(mfi);
                const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                    for (int q = fi(i, j, k, 0); q < fi(i, j, k, 0) + fi(i, j, k, 1); ++q)
                    {
                        const FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                        if (e.blk || !(e.sp >= 0.0)) continue;
                        fs(i, j, k) = e.sp;
                        ft(i, j, k) = std::max(ft(i, j, k), e.Tg);
                    }
            }
        }
    }

    // D2b: burner-on-feet descent (surface_patch.nozzle_descent = feet), lev 0,
    // once per step, before the stand-offs are computed. The foot annulus is
    // the live columns with r in [foot_r_inner, foot_r_outer] about (x0, y0);
    // z_face = PLO_z + (k_top + 1)*dz.
    //   pads: the annulus is split into foot_npads equal sectors of its
    //     observed azimuth range [a_min, a_max] (atan2 about (x0, y0)); pad
    //     p = min(npads - 1, floor((a - a_min)/(a_max - a_min)*npads)) (all in
    //     pad 0 if a_max == a_min). A pad's height is the nearest-rank
    //     quantile of its z_face values: sorted ascending, index
    //     ceil(q*n) - 1 (q = foot_pad_quantile; q = 1 is the max). z_foot =
    //     the mean of the pad heights.
    // z_target = z_foot + foot_standoff; z_n = min(z_n_prev, max(z_target,
    // z_n_prev - nozzle_feed_max*dt)) (never rises; first step z_n =
    // z_target). foot_carry_cols counts columns with z_face >= their pad's
    // height; foot_stall_time is the simulated time since z_n last decreased
    // (0 on a step that descends). A stall is logged, never aborted. k_top_col is the domain-wide reduced
    // mask top, so every rank computes the same z_n. Aborts if the annulus or
    // a pad has no live column. nozzle_z is NOT checkpointed.
    // D2c (foot_body_clearance = 1, pads only): each annulus column whose
    // k_top exceeds its pad's nearest-rank k (the same pick on the integer
    // tops) is armed in foot_clear_k and cut to the pad height at the next
    // removal pass (the column path (deleted, P2-28)). Clipping values above a nearest-rank
    // quantile down to it leaves the quantile unchanged (no ratchet).
    Set::Scalar FeetDescent(int lev, Set::Scalar dt, const amrex::Vector<int>& k_top_col)
    {
        if (!SurfaceOwner(lev))
            Util::Abort(INFO, "surface_patch.nozzle_descent = feet runs on the surface-owning level only");
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const auto& sp = surface_patch;
        const Set::Scalar ri2 = sp.foot_r_inner * sp.foot_r_inner;
        const Set::Scalar ro2 = sp.foot_r_outer * sp.foot_r_outer;
        std::vector<Set::Scalar> zf, az;
        std::vector<int> kf;                // D2c clearance: k_top and column
        std::vector<std::size_t> cf;
        const bool clear = sp.foot_body_clearance;
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t col = static_cast<std::size_t>(j) * nx + i;
            const int k_top = k_top_col[col];
            if (k_top == std::numeric_limits<int>::min()) continue;
            const Set::Scalar x = PLO[0] + (Set::Scalar(i + domain.smallEnd(0)) + 0.5) * DX[0];
            const Set::Scalar y = PLO[1] + (Set::Scalar(j + domain.smallEnd(1)) + 0.5) * DX[1];
            const Set::Scalar dxp = x - sp.x0, dyp = y - sp.y0;
            const Set::Scalar r2 = dxp * dxp + dyp * dyp;
            if (r2 < ri2 || r2 > ro2) continue;
            zf.push_back(PLO[2] + Set::Scalar(k_top + 1) * DX[2]);
            az.push_back(std::atan2(dyp, dxp));
            if (clear) { kf.push_back(k_top); cf.push_back(col); }
        }
        if (zf.empty())
            Util::Abort(INFO, "surface_patch feet: the foot annulus [", sp.foot_r_inner, ", ",
                        sp.foot_r_outer, "] m has no live column");
        const int npads = sp.foot_npads;
        std::vector<int> pad(zf.size(), 0);
        if (npads > 1)
        {
            const Set::Scalar a_min = *std::min_element(az.begin(), az.end());
            const Set::Scalar a_max = *std::max_element(az.begin(), az.end());
            if (a_max > a_min)
                for (std::size_t c = 0; c < zf.size(); ++c)
                    pad[c] = std::min(npads - 1, static_cast<int>(
                                 std::floor((az[c] - a_min) / (a_max - a_min) * npads)));
        }
        std::vector<Set::Scalar> pad_h(npads, 0.0);
        std::vector<int> pad_k(npads, 0);
        Set::Scalar z_foot = 0.0;
        for (int p = 0; p < npads; ++p)
        {
            std::vector<Set::Scalar> v;
            for (std::size_t c = 0; c < zf.size(); ++c) if (pad[c] == p) v.push_back(zf[c]);
            if (clear)
            {
                // the same nearest-rank pick on the integer tops (z_face is
                // monotone in k, so this is the pad height's cell index)
                std::vector<int> vk;
                for (std::size_t c = 0; c < kf.size(); ++c) if (pad[c] == p) vk.push_back(kf[c]);
                std::sort(vk.begin(), vk.end());
                if (!vk.empty())
                {
                    long ik = static_cast<long>(std::ceil(sp.foot_pad_quantile * Set::Scalar(vk.size()))) - 1;
                    ik = std::max(0L, std::min(ik, static_cast<long>(vk.size()) - 1));
                    pad_k[p] = vk[static_cast<std::size_t>(ik)];
                }
            }
            if (v.empty())
                Util::Abort(INFO, "surface_patch feet: pad ", p, " of ", npads,
                            " has no live column in the foot annulus");
            std::sort(v.begin(), v.end());
            Set::Scalar hp = 0.0;
            {
                long idx = static_cast<long>(std::ceil(sp.foot_pad_quantile * Set::Scalar(v.size()))) - 1;
                idx = std::max(0L, std::min(idx, static_cast<long>(v.size()) - 1));
                hp = v[static_cast<std::size_t>(idx)];
            }
            pad_h[p] = hp;
            z_foot += hp;
        }
        z_foot /= Set::Scalar(npads);
        Set::Scalar carry = 0.0;
        for (std::size_t c = 0; c < zf.size(); ++c) if (zf[c] >= pad_h[pad[c]]) carry += 1.0;
        const Set::Scalar z_target = z_foot + sp.foot_standoff;
        if (!feet_init)
        {
            nozzle_z  = z_target;
            feet_init = true;
        }
        else
        {
            const Set::Scalar z_prev = nozzle_z;
            nozzle_z = std::min(nozzle_z, std::max(z_target, nozzle_z - sp.nozzle_feed_max * dt));
            // pad heights move in whole cells, so a per-step "no descent"
            // count grows on almost every step; time since the last descent
            // is the stall signal
            foot_stall_time = (nozzle_z < z_prev) ? 0.0 : foot_stall_time + dt;
        }
        foot_z          = z_foot;
        foot_cols       = Set::Scalar(zf.size());
        foot_carry_cols = carry;
        if (clear)
        {
            // D2c: annulus columns standing proud of their pad are cut to the
            // pad height at the next removal pass (the column path (deleted, P2-28), before PASS 1).
            foot_clear_k.assign(static_cast<std::size_t>(nx) * ny, std::numeric_limits<int>::min());
            for (std::size_t c = 0; c < kf.size(); ++c)
                if (kf[c] > pad_k[pad[c]]) foot_clear_k[cf[c]] = pad_k[pad[c]];
        }
        return nozzle_z;
    }

    // D2a2: wall-jet enthalpy closure (surface_patch.jet_closure = enthalpy).
    // Live in-patch columns with s > 0 are binned by r into bins of width
    // jet_dr (default max(dx, dy)). Stagnation (free-jet centreline decay):
    //   s_c    = mean s over the innermost non-empty bin
    //   phi    = min(1, jet_core_length*jet_D/s_c)   (1 for jet_stagnation =
    //            nozzle, or when no column marches: s_c = 0)
    //   T_stag = T_mix + (jet_T_nozzle - T_mix)*phi, T_mix = T_ent, or (D2b,
    //            jet_T_ent_mode = exhaust) the previous step's T_exhaust
    //   m      = jet_mdot, or jet_mdot/phi with jet_entrained_mass = 1
    // then marching outward from the innermost bin with T_gas = T_stag
    // (explicit in r):
    //   P_bin   = sum_bin max(0, q_robin(c; T_gas)) * dx*dy
    //   flame_Tg[c] = T_gas for the bin's columns
    //   share   = min(P_bin, m*cp*(T_gas - T_ent))
    //   T_gas  <- T_gas - share/(m*cp)             (>= T_ent)
    // q_robin is SurfaceCellFlux's h*(T_gas - T_s) with the same inputs as the
    // H update (top-cell T_old and k_c, robin_form, and, for pinned, the pin
    // from the pin rule (deleted, P2-28) at this h and T_gas), so the surface temperature and branch
    // match the kernel: the pin rule (deleted, P2-28), which runs next, reproduces the pin
    // exactly from the final flame vectors. Columns hotter than the gas
    // (q_robin < 0) neither absorb nor return enthalpy. jet_P_face = sum share
    // is the enthalpy taken from the nozzle stream; if a bin would overdraw
    // (share < P_bin), the excess is supplied by the T_ent floor (the
    // recirculation reservoir) and is not counted. With jet_P_cap =
    // m*cp*(T_stag - T_ent) and jet_P_decay = jet_mdot*cp*(jet_T_nozzle - T_ent)
    // - jet_P_cap (the excess the decay line discards; 0 for nozzle or
    // entrained mass), the invariant (abort) is jet_P_face <= jet_P_cap,
    // jet_P_face + jet_P_exhaust == jet_P_cap, jet_P_decay >= 0 (== 0 unless
    // decay without entrained mass), all to 1e-9 of the nozzle budget. Under
    // exhaust the budget is jet_mdot*cp*(jet_T_nozzle - T_ent) + jet_P_recirc,
    // jet_P_recirc = (m - jet_mdot)*cp*(T_rec - T_ent), and the march floor
    // stays T_ent (a T_rec floor would ratchet T_exhaust upward).
    // jet_r_reach = max column r over bins with P_bin > 0. Every rank gathers the same top-cell data and runs the same
    // serial march, so the vectors are rank-identical. With
    // jet_profile_interval > 0, rank 0 appends one row per non-empty bin
    // (time,bin,r_lo,n_cols,T_gas,P_bin,P_share) to <plot_file>_jet_profile.csv
    // at the first step of each interval (T_gas = the bin's inlet value, the
    // one its columns receive). Not checkpointed.
    // D2c options (all default-off; details at each branch):
    //   jet_decay_diameter = momentum: phi uses D_e = jet_De_ref*sqrt(T_mix/T_ent);
    //   jet_bin_update = exponential: columns see T_eff (exact affine bin);
    //     profile column T_eff;
    //   jet_negative_flux = count: P_bin sums the signed q_robin (jet_P_neg);
    //     profile column P_neg;
    //   jet_free_surface = 1: ambient dilution after each bin holding
    //     free-surface columns; T_rec (exhaust) = the inlet T of the first such
    //     bin (jet_T_mouth; T_exhaust if none); profile columns n_fs, m_ratio
    //     (m after the bin / m_stag). Empty bins do not dilute.
    void JetEnthalpyMarch(int lev, Set::Scalar time, const amrex::Vector<int>& k_top_col)
    {
        if (!SurfaceOwner(lev))
            Util::Abort(INFO, "surface_patch.jet_closure = enthalpy runs on the surface-owning level only");
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;
        // P2-28: the march bills the face list's faces (its column flux, and the top-cell gather it
        // read, are deleted with the column path; Parse requires the list under jet_closure = enthalpy)
        if (!(surface_patch.face_list && fl_ready))
            Util::Abort(INFO, "the jet march bills the face list's faces: the list is not built (P2-28)");

        // Bin the marching columns by r.
        const Set::Scalar dr = surface_patch.jet_dr > 0.0 ? surface_patch.jet_dr
                                                          : std::max(DX[0], DX[1]);
        std::vector<std::vector<std::size_t>> bins;
        std::vector<Set::Scalar> col_r(ncols, 0.0);
        // D2c free surface: z0 = each column's face height at the first march
        // call (host state); a column is free-surface when r > foot_r_outer
        // and z0 - z_face <= jet_fs_aspect*(r - foot_r_outer).
        std::vector<unsigned char> col_fs;
        if (surface_patch.jet_free_surface)
        {
            if (jet_fs_z0.size() != ncols)
            {
                jet_fs_z0.assign(ncols, 0.0);
                for (std::size_t c = 0; c < ncols; ++c)
                    if (k_top_col[c] != std::numeric_limits<int>::min())
                        jet_fs_z0[c] = PLO[2] + Set::Scalar(k_top_col[c] + 1) * DX[2];
            }
            col_fs.assign(ncols, 0);
        }
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t col = static_cast<std::size_t>(j) * nx + i;
            if (k_top_col[col] == std::numeric_limits<int>::min()) continue;
            if (!(flame_s[col] > 0.0)) continue;
            const Set::Scalar x = PLO[0] + (Set::Scalar(i + ilo_dom) + 0.5) * DX[0];
            const Set::Scalar y = PLO[1] + (Set::Scalar(j + jlo_dom) + 0.5) * DX[1];
            const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
            const Set::Scalar r2 = dxp * dxp + dyp * dyp;
            if (r2 > r2max) continue;
            const Set::Scalar r = std::sqrt(r2);
            const std::size_t b = static_cast<std::size_t>(std::floor(r / dr));
            if (b >= bins.size()) bins.resize(b + 1);
            bins[b].push_back(col);
            col_r[col] = r;
            if (surface_patch.jet_free_surface && r > surface_patch.foot_r_outer)
            {
                const Set::Scalar z_face = PLO[2] + Set::Scalar(k_top_col[col] + 1) * DX[2];
                if (jet_fs_z0[col] - z_face <= surface_patch.jet_fs_aspect * (r - surface_patch.foot_r_outer))
                    col_fs[col] = 1;
            }
        }

        // Stagnation: s_c = mean s over the innermost non-empty bin (bin 0
        // unless the axis columns are excluded); phi = min(1, L_core*D/s_c)
        // (1 under jet_stagnation = nozzle, or when no column marches);
        // T_stag = T_ent + (T_nozzle - T_ent)*phi; march mass m = mdot, or
        // mdot/phi under jet_entrained_mass = 1.
        const Set::Scalar T_ent = surface_patch.jet_T_ent;
        Set::Scalar s_c = 0.0;
        for (const auto& bin : bins)
        {
            if (bin.empty()) continue;
            for (std::size_t c : bin) s_c += flame_s[c];
            s_c /= Set::Scalar(bin.size());
            break;
        }
        // D2b: exhaust recirculation mixes the nozzle flow with the previous
        // step's exhaust (T_rec) instead of T_ent; the floor stays T_ent.
        const bool exhaust = surface_patch.jet_T_ent_exhaust;
        if (exhaust && !jet_rec_init) { jet_T_rec_next = T_ent; jet_rec_init = true; }
        const Set::Scalar T_mix = exhaust ? jet_T_rec_next : T_ent;
        // D2c: decay diameter, jet_D (nozzle) or the momentum diameter in
        // surroundings at T_mix: D_e = jet_De_ref*sqrt(T_mix/T_ent) (D_e =
        // 2 mdot/sqrt(pi rho J) at fixed mdot and J, rho ~ 1/T).
        const Set::Scalar D_dec = surface_patch.jet_De_momentum
                                ? surface_patch.jet_De_ref * std::sqrt(T_mix / T_ent)
                                : surface_patch.jet_D;
        Set::Scalar phi = 1.0;
        if (surface_patch.jet_stag_decay && s_c > 0.0)
            phi = std::min(Set::Scalar(1.0),
                           surface_patch.jet_core_length * D_dec / s_c);
        const Set::Scalar T_stag = T_mix + (surface_patch.jet_T_nozzle - T_mix) * phi;
        const Set::Scalar m_jet  = surface_patch.jet_entrained_mass
                                 ? surface_patch.jet_mdot / phi : surface_patch.jet_mdot;
        // D2c: m (and mcp) change only under free-surface dilution.
        Set::Scalar m_cur = m_jet;
        Set::Scalar mcp   = m_jet * surface_patch.jet_cp;
        const Set::Scalar mcp0  = mcp;
        const Set::Scalar dA    = DX[0] * DX[1];
        const bool expo  = surface_patch.jet_bin_exponential;
        const bool count = surface_patch.jet_neg_count;
        const bool fsurf = surface_patch.jet_free_surface;
        Set::Scalar T_gas = T_stag;
        Set::Scalar P_face = 0.0, r_reach = 0.0;
        Set::Scalar fl_Pm = 0.0;   // P2-8: the list faces billed, before any cap
        Set::Scalar P_neg = 0.0, fs_cols = 0.0, T_mouth = 0.0;
        bool have_mouth = false;
        // the robin power of column c's listed faces at gas temperature Tg, billed to its bin (P2-9 R1:
        // every listed face of the column -- the top included -- at its cell's T_s, FrFaceFlux, x w; the
        // face area dx*dz against the bin's dA = dx*dy, cubic)
        auto q_col = [&](std::size_t c, Set::Scalar Tg) -> Set::Scalar
        {
            Set::Scalar q_tot = 0.0;
            for (int q = fl_col_first[c]; q < fl_col_first[c] + fl_col_count[c]; ++q)
            {
                const FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                FrFaceFlux(e, Tg, e.K, e.T, DX[0], a_in, a_out, a_rob, a_los);
                q_tot += e.w * a_rob * (DX[2] / DX[0]);
            }
            return q_tot;
        };
        // P2-11b: face_heat = path -- the march bins FACES by floor(s_path/dr):
        // every unblocked, reached face whose column marches (flame_s > 0) and
        // whose centre lies within the patch radius (its TRUE r). A face is
        // billed alone, exactly as q_col bills it inside its column. Faces the
        // march does not visit keep T_ent. Every other rule is the column
        // march's (s_c from the column bins, decay, entrainment, exponential
        // update, negative flux, invariants); free surface and the pin are
        // refused under path (Parse).
        const bool path = surface_patch.face_heat_path;
        std::vector<std::vector<std::size_t>> fbins;
        std::vector<Set::Scalar> face_r;
        if (path)
        {
            face_r.assign(fl_faces.size(), 0.0);
            for (std::size_t q = 0; q < fl_faces.size(); ++q)
            {
                FaceEntry& e = fl_faces[q];
                e.Tg = T_ent;
                if (e.blk || !(e.sp >= 0.0)) continue;
                if (!(flame_s[static_cast<std::size_t>(e.col)] > 0.0)) continue;
                Set::Scalar x, y, z;
                FaceCentre(e, PLO, DX, x, y, z);
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                const Set::Scalar r2 = dxp * dxp + dyp * dyp;
                if (r2 > r2max) continue;
                const std::size_t b = static_cast<std::size_t>(std::floor(e.sp / dr));
                if (b >= fbins.size()) fbins.resize(b + 1);
                fbins[b].push_back(q);
                face_r[q] = std::sqrt(r2);
            }
        }
        auto q_face = [&](std::size_t q, Set::Scalar Tg) -> Set::Scalar
        {
            const FaceEntry& e = fl_faces[q];
            Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
            FrFaceFlux(e, Tg, e.K, e.T, DX[0], a_in, a_out, a_rob, a_los);
            return e.w * a_rob * (DX[2] / DX[0]);
        };
        // one march over items: columns (position) or faces (path)
        auto q_item = [&](std::size_t c, Set::Scalar Tg) -> Set::Scalar
        { return path ? q_face(c, Tg) : q_col(c, Tg); };
        const std::vector<std::vector<std::size_t>>& mbins = path ? fbins : bins;
        // Optional profile log: one row per non-empty bin (inlet T_gas, P_bin,
        // the jet's share), rank 0, at the first step of each interval.
        bool prof = false;
        std::ostringstream prow;
        if (surface_patch.jet_profile_interval > 0.0)
        {
            const long idx = static_cast<long>(std::floor(time / surface_patch.jet_profile_interval));
            if (idx > jet_profile_last) { prof = true; jet_profile_last = idx; }
            prow.setf(std::ios::scientific);
            prow.precision(16);
        }
        for (std::size_t ib = 0; ib < mbins.size(); ++ib)
        {
            const auto& bin = mbins[ib];
            if (bin.empty()) continue;
            const Set::Scalar T_in = T_gas;
            // D2c exponential: the bin's columns see the bin-mean gas
            // temperature of an affine q(T) = G(T - T_eq) exchanger,
            //   T_eff = T_eq + (T_in - T_eq)(1 - e^-N)/N,  N = G/(m cp),
            // G and P_in from q at T_in and T_in - 1 K. P_bin is re-evaluated
            // at T_eff, so T_out = T_in - P_bin/(m cp) stays exact. G <= 0
            // falls back to explicit (T_eff = T_in).
            Set::Scalar T_col = T_gas;
            if (expo)
            {
                Set::Scalar P_in = 0.0, G = 0.0;
                for (std::size_t c : bin)
                {
                    Set::Scalar q1 = q_item(c, T_in), q0 = q_item(c, T_in - 1.0);
                    if (!count) { q1 = std::max(Set::Scalar(0.0), q1); q0 = std::max(Set::Scalar(0.0), q0); }
                    P_in += q1;
                    G    += q1 - q0;
                }
                P_in *= dA;
                G    *= dA;
                if (G > 0.0)
                {
                    const Set::Scalar T_eq = T_in - P_in / G;
                    const Set::Scalar N    = G / mcp;
                    T_col = T_eq + (T_in - T_eq) * (-std::expm1(-N) / N);
                }
            }
            Set::Scalar P_bin = 0.0, P_nb = 0.0;
            for (std::size_t c : bin)
            {
                if (path) fl_faces[c].Tg = T_col;
                else      flame_Tg[c] = T_col;
                const Set::Scalar q_rob = q_item(c, T_col);
                fl_Pm += q_rob * dA;
                if (count)
                {
                    // D2c: rock hotter than the gas heats the gas.
                    P_bin += q_rob;
                    if (q_rob < 0.0) P_nb += q_rob;
                }
                else if (q_rob > 0.0) P_bin += q_rob;
            }
            P_bin *= dA;
            P_neg += P_nb * dA;
            if (P_bin > 0.0)
                for (std::size_t c : bin) r_reach = std::max(r_reach, path ? face_r[c] : col_r[c]);
            // avail >= 0, so a negative P_bin (count) always takes the second
            // branch and raises T_gas; the clamp acts only when P_bin > 0.
            const Set::Scalar avail = mcp * (T_gas - T_ent);
            const Set::Scalar share = (P_bin >= avail) ? avail : P_bin;
            if (P_bin >= avail) { P_face += avail; T_gas = T_ent; }
            else                { P_face += P_bin; T_gas = T_gas - P_bin / mcp; }
            // D2c free surface: after the bin's heat, dilute with ambient air
            //   m_out = m (r_hi/r_lo)^(e n_fs/n),  T_gas <- T_ent + (T_gas - T_ent) m/m_out
            // (the excess enthalpy flux m cp (T - T_ent) is conserved).
            std::size_t n_fs = 0;
            if (fsurf)
            {
                for (std::size_t c : bin) n_fs += col_fs[c];
                if (n_fs > 0)
                {
                    if (!have_mouth) { T_mouth = T_in; have_mouth = true; }
                    const Set::Scalar r_lo = Set::Scalar(ib) * dr;
                    const Set::Scalar r_hi = Set::Scalar(ib + 1) * dr;
                    if (!(r_lo > 0.0))
                        Util::Abort(INFO, "surface_patch.jet_free_surface: a free-surface column "
                                          "is in the axis bin (r_lo = 0)");
                    const Set::Scalar m_out = m_cur * std::pow(r_hi / r_lo,
                        surface_patch.jet_fs_exponent * Set::Scalar(n_fs) / Set::Scalar(bin.size()));
                    T_gas = T_ent + (T_gas - T_ent) * (m_cur / m_out);
                    m_cur = m_out;
                    mcp   = m_cur * surface_patch.jet_cp;
                    fs_cols += Set::Scalar(n_fs);
                }
            }
            if (prof)
            {
                std::ostringstream row;
                row.copyfmt(prow);
                row << time << "," << ib << "," << Set::Scalar(ib) * dr << "," << bin.size()
                     << "," << T_in << "," << P_bin << "," << share;
                if (expo)  row << "," << T_col;
                if (fsurf) row << "," << n_fs << "," << m_cur / m_jet;
                if (count) row << "," << P_nb * dA;
                row << "\n";
                CsvCheckRow(row.str(), JetProfileColumns(expo, fsurf, count).size(), "_jet_profile.csv");
                prow << row.str();
            }
        }
        if (prof && amrex::ParallelDescriptor::IOProcessor())
        {
            const std::string path = plot_file + "_jet_profile.csv";
            const bool first = !jet_profile_header_written;
            std::ofstream ofs(path, (first ? std::ios::trunc : std::ios::app) | std::ios::out);
            if (ofs.is_open())
            {
                if (first) ofs << CsvHeader(JetProfileColumns(expo, fsurf, count));
                ofs << prow.str();
            }
        }
        if (prof) jet_profile_header_written = true;
        const Set::Scalar P_nozzle = surface_patch.jet_mdot * surface_patch.jet_cp
                                   * (surface_patch.jet_T_nozzle - T_ent);
        jet_P_cap     = mcp0 * (T_stag - T_ent);
        jet_T_stag    = T_stag;
        jet_s_c       = s_c;
        Set::Scalar P_budget = P_nozzle;
        if (exhaust)
        {
            // P_face + P_exhaust = mdot*cp*(T_nozzle - T_ent) + P_recirc.
            // With entrained mass m = mdot/phi, P_cap = m cp (T_stag - T_ent)
            // = mdot cp (T_nozzle - T_mix) + m cp (T_mix - T_ent) = P_budget
            // for any T_mix, so P_decay = 0 holds for exhaust (T_mix = the
            // previous step's T_rec).
            jet_T_rec    = T_mix;
            jet_P_recirc = (m_jet - surface_patch.jet_mdot) * surface_patch.jet_cp * (T_mix - T_ent);
            P_budget     = P_nozzle + jet_P_recirc;
            // D2c: under free-surface dilution the recirculating gas is the
            // gas in the hole: the inlet of the first free-surface bin.
            jet_T_rec_next = have_mouth ? T_mouth : T_gas;
        }

        jet_P_decay   = P_budget - jet_P_cap;
        jet_P_face    = P_face;
        fl_P_march = fl_Pm;
        // the gas leaving the floor, which is the gas leaving the domain (P2-10:
        // the wall march that carried it on up the bore is deleted).
        jet_T_exhaust = T_gas;
        jet_P_exhaust = mcp * (T_gas - T_ent);
        jet_r_reach   = r_reach;
        jet_D_e       = D_dec;
        jet_fs_cols   = fs_cols;
        jet_m_ratio   = m_cur / m_jet;
        jet_T_mouth   = have_mouth ? T_mouth : T_gas;
        jet_P_neg     = P_neg;
        // Invariants (abort), tol = 1e-9 of the budget:
        //   P_face <= P_cap + tol;  P_face + P_exhaust == P_cap
        //     (P_exhaust = m_exit cp (T_exhaust - T_ent); dilution conserves it);
        //   P_decay >= -tol, and == 0 unless decay without entrained mass;
        //   drop:  P_face >= 0;
        //   count: P_face is the signed net, P_face - P_neg >= -tol (P_neg <= 0).
        const Set::Scalar tol = 1.0e-9 * std::max(P_budget, Set::Scalar(1.0));
        if (!(std::isfinite(P_face) &&
              (count ? P_face - P_neg >= -tol : P_face >= 0.0) &&
              P_face <= jet_P_cap + tol &&
              std::abs(P_face + jet_P_exhaust - jet_P_cap) <= tol &&
              jet_P_decay >= -tol &&
              (surface_patch.jet_stag_decay && !surface_patch.jet_entrained_mass
                   ? true : std::abs(jet_P_decay) <= tol)))
            Util::Abort(INFO, "jet enthalpy invariant violated at time ", time,
                        ": jet_P_face = ", P_face, ", jet_P_exhaust = ", jet_P_exhaust,
                        ", jet_P_cap = ", jet_P_cap, ", jet_P_decay = ", jet_P_decay,
                        ", jet_P_neg = ", P_neg,
                        ", nozzle budget = ", P_nozzle);
    }

    // =====================================================================
    // P2-8 -- THE SURFACE FACE LIST: one list, one heat rule, one area rule.
    // =====================================================================
    // F1 (BuildFaceList)  every face between a SOLID cell and a VOID cell in
    //     +-x, +-y, +-z, inside the patch, from removed_mf. The domain's +z
    //     boundary is the rock's open surface (gas above); the other domain
    //     boundaries are NOT faces (symmetry planes, base, outer walls).
    //     Built ONCE per step, after the column tops and before the march,
    //     gathered to every rank and sorted canonically (column, k, direction),
    //     so every consumer reads the SAME list in the SAME order whatever the
    //     decomposition.
    // F2 (FaceH, via FaceListHeat)  h = h_expr(x_f, y_f, r_f, s_f, t) at the
    //     FACE CENTRE, s_f = z_n - z_f; 0 where s_f <= 0 (enthalpy closure).
    //     T_gas = flame_Tg[c].
    //     At a column-top +z face the arguments are BuildFlameColumns' own, so
    //     it reproduces flame_h[c] bit for bit (fl_h_top_mismatch == 0).
    // F3 (FaceWeight)  w = 1/(|nx|+|ny|+|nz|), n the unit outward normal of the
    //     SMOOTHED solid fraction (FrSmoothedGradient, P2-9 rev 2). Flat floor:
    //     w = 1 exactly. 45 deg staircase: 1/sqrt2 per face, summing to the
    //     plane. Zero gradient => w = 1, counted (fl_w_fallback). One w per
    //     CELL, applied to each of its faces.
    // F4  every listed face -- the column top included -- is heated at its
    //     cell's one T_s, x w (FrFaceFlux, P2-9).
    // F5  the kernel, the ledger, the march (billing each face to its column's
    //     bin) and the line rule read F1-F4 and nothing else.
    struct FaceEntry
    {
        int i = 0, j = 0, k = 0;   // the solid cell (global index)
        int dir = 0;               // 0 -x, 1 +x, 2 -y, 3 +y, 4 -z, 5 +z
        int cls = 0;               // FL_TOP / FL_UP / FL_SIDE / FL_DOWN
        int col = 0;               // column index (j - jlo)*nx + (i - ilo)
        int fb = 0;                // 1 = w by zero-gradient fallback
        Set::Scalar w = 1.0;       // F3
        Set::Scalar nx = 0.0, ny = 0.0, nz = 0.0;   // F3 unit outward normal
        Set::Scalar h = 0.0;       // F2
        Set::Scalar T = 0.0, K = 0.0;   // the cell's T_old and conductivity
        // P2-9 R1 (FrAssignTs; not in the gather): the CELL's surface
        // temperature, or tsfb = 1 for the face Newton on fl_faces[qline]
        Set::Scalar Ts = 0.0;
        int tsfb = 1, qline = 0;
        // P2-11b item 5: the face's gas temperature, set after the march and
        // read by every per-face consumer (kernel, march, diagnostics, ledger,
        // the fallback flake). Under face_heat = position it IS flame_Tg[col].
        Set::Scalar Tg = 0.0;
        // P2-11b item 4 (face_heat = path only): the surface-path distance of
        // the face's cell from the stagnation sources [m] (-1: blocked or
        // unreached), and 1 = blocked by surface_patch.barrier_expr.
        Set::Scalar sp = -1.0;
        int blk = 0;
        // P2-28: a PRESCRIBED top face -- the patch's held temperature T_expr(t) [K]
        // (FaceListHeat; -1 = none). The kernel, the ledger and the diagnostics bill it
        // by the half-cell Dirichlet term 2 kc (T_presc - T_c) / dz^2, and the criterion
        // reads it as the cell's surface temperature (FrAssignTs, RecedeAlongLines R1).
        Set::Scalar T_presc = -1.0;
    };
    static constexpr int FL_TOP = 0, FL_UP = 1, FL_SIDE = 2, FL_DOWN = 3;
    static constexpr int FL_NF = 13;     // doubles per packed entry
    std::vector<FaceEntry>   fl_faces;
    std::vector<int>         fl_col_first, fl_col_count;
    std::unique_ptr<amrex::iMultiFab> fl_index;   // per cell: first entry, count
    bool fl_ready = false;
    Set::Field<Set::Scalar> fl_w_mf, fl_h_mf;
    Set::Scalar fl_n_top = 0.0, fl_n_up = 0.0, fl_n_side = 0.0, fl_n_down = 0.0;
    Set::Scalar fl_P_top = 0.0, fl_P_up = 0.0, fl_P_side = 0.0, fl_P_down = 0.0;
    Set::Scalar fl_P_march = 0.0, fl_w_fallback = 0.0;
    Set::Scalar fl_h_top_mismatch = 0.0;
    // b2's identities, formed in DOUBLE (thermo prints %g, six digits, which
    // cannot resolve a 1e-9 bar): |jet_P_face - sum fl_P| / |sum fl_P| (the
    // march vs the kernel's terms) and |ledger_P_robin - sum fl_P| / |sum fl_P|
    // (the ledger's per-cell mirror vs the same).
    Set::Scalar fl_id_march = 0.0, fl_id_ledger = 0.0;
    // P2-11b item 6 (face_heat = path only; registered last): the stagnation
    // standoff H = s_c [m], the source / blocked / unreached counts, the
    // largest reached s_path [m], the path's wall time per step [s]; plotfile
    // fl_s (m; -1 blocked or unreached) and fl_Tg (K).
    Set::Scalar fl_H_stag = 0.0, fl_n_src = 0.0, fl_n_blocked = 0.0;
    Set::Scalar fl_n_unreached = 0.0, fl_s_max = 0.0, fl_t_path = 0.0;
    Set::Field<Set::Scalar> fl_s_mf, fl_Tg_mf;

    // =====================================================================
    // P2-11b -- THE SURFACE PATH (surface_patch.face_heat = path). Heat follows
    // the jet along the rock: every listed face is heated by the deck's h_expr
    // evaluated with r := s_path (the distance ALONG THE SURFACE from the
    // stagnation point) and s := H (the stagnation standoff), and the march
    // bins FACES by s_path. One routine, ported from
    // studies/p2_11b_path/offline/path.py (a0):
    //   sources  the cell under the +z face (FL_TOP) of each column in the
    //            march's innermost non-empty r-bin -- s_c's own population --
    //            at d0 = that face centre's radius from (x0, y0); H = s_c
    //   barrier  a face is BLOCKED when barrier_expr(x, y, z, r) > 0 at the
    //            point 0.5 dx outside its centre: h = 0, out of the graph and
    //            the march
    //   graph    the cells holding at least one unblocked listed face; an edge
    //            joins two nodes whose index offset o has |o| <= 4 cells,
    //            weight dx |o| (Euclidean). Plain binary-heap Dijkstra.
    //   unreached faces get h = 0 and are counted.
    // Every rank holds the same replicated list, so every rank computes the
    // same distances. Returns H (= s_c; 0 when no column marches).
    // =====================================================================
    static const std::vector<std::array<int, 4>>& FrPathOffsets()
    {
        // (di, dj, dk, index into FrPathLen); |o|^2 <= 16
        static const std::vector<std::array<int, 4>> off = [] {
            std::vector<std::array<int, 4>> o;
            for (int a = -4; a <= 4; ++a)
            for (int b = -4; b <= 4; ++b)
            for (int c = -4; c <= 4; ++c)
            {
                const int n2 = a * a + b * b + c * c;
                if (n2 == 0 || n2 > 16) continue;
                o.push_back({a, b, c, n2});
            }
            return o;
        }();
        return off;
    }

    Set::Scalar FaceListPath(int lev, Set::Scalar time,
                             const amrex::Vector<int>& k_top_col)
    {
        const double t0 = amrex::second();
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dx = DX[0];
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;

        // 1. the march's innermost non-empty r-bin, built exactly as
        //    JetEnthalpyMarch builds it (column centres, dr, flame_s > 0, the
        //    radius cut), and s_c = the mean s over it in the same order
        const Set::Scalar dr = surface_patch.jet_dr > 0.0 ? surface_patch.jet_dr
                                                          : std::max(DX[0], DX[1]);
        std::size_t b_in = std::numeric_limits<std::size_t>::max();
        std::vector<std::pair<std::size_t, Set::Scalar>> src_cols;   // (col, r)
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t col = static_cast<std::size_t>(j) * nx + i;
            if (k_top_col[col] == std::numeric_limits<int>::min()) continue;
            if (!(flame_s[col] > 0.0)) continue;
            const Set::Scalar x = PLO[0] + (Set::Scalar(i + ilo_dom) + 0.5) * DX[0];
            const Set::Scalar y = PLO[1] + (Set::Scalar(j + jlo_dom) + 0.5) * DX[1];
            const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
            const Set::Scalar r2 = dxp * dxp + dyp * dyp;
            if (r2 > r2max) continue;
            const Set::Scalar r = std::sqrt(r2);
            const std::size_t b = static_cast<std::size_t>(std::floor(r / dr));
            if (b < b_in) { b_in = b; src_cols.clear(); }
            if (b == b_in) src_cols.emplace_back(col, r);
        }
        Set::Scalar H = 0.0;
        if (!src_cols.empty())
        {
            for (const auto& sc : src_cols) H += flame_s[sc.first];
            H /= Set::Scalar(src_cols.size());
        }

        // 2. the cells (a cell's faces are contiguous: the list is sorted by
        //    column, k, dir) and the barrier
        const std::size_t nf = fl_faces.size();
        std::vector<int> cell_q0, cell_nq, cell_i, cell_j, cell_k;
        std::vector<unsigned char> cell_node;
        long n_blk = 0;
        for (std::size_t q = 0; q < nf; ++q)
        {
            FaceEntry& e = fl_faces[q];
            if (q == 0 || e.i != fl_faces[q - 1].i || e.j != fl_faces[q - 1].j
                       || e.k != fl_faces[q - 1].k)
            {
                cell_q0.push_back(static_cast<int>(q));
                cell_nq.push_back(0);
                cell_i.push_back(e.i); cell_j.push_back(e.j); cell_k.push_back(e.k);
                cell_node.push_back(0);
            }
            ++cell_nq.back();
            e.blk = 0;
            e.sp  = -1.0;
            if (surface_patch.has_barrier)
            {
                Set::Scalar x, y, z;
                FaceCentre(e, PLO, DX, x, y, z);
                const int ax = e.dir / 2;
                const Set::Scalar sg = (e.dir % 2 == 1) ? 0.5 : -0.5;
                if (ax == 0) x += sg * DX[0];
                else if (ax == 1) y += sg * DX[1];
                else z += sg * DX[2];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (surface_patch.barrier_f(x, y, z, std::sqrt(dxp * dxp + dyp * dyp)) > 0.0)
                { e.blk = 1; ++n_blk; }
            }
            if (!e.blk) cell_node.back() = 1;
        }
        const std::size_t nc = cell_q0.size();

        // 3. a dense index over the nodes' bounding box
        int i0 = std::numeric_limits<int>::max(), j0 = i0, k0 = i0;
        int i1 = std::numeric_limits<int>::min(), j1 = i1, k1 = i1;
        for (std::size_t n = 0; n < nc; ++n)
        {
            if (!cell_node[n]) continue;
            i0 = std::min(i0, cell_i[n]); i1 = std::max(i1, cell_i[n]);
            j0 = std::min(j0, cell_j[n]); j1 = std::max(j1, cell_j[n]);
            k0 = std::min(k0, cell_k[n]); k1 = std::max(k1, cell_k[n]);
        }
        std::vector<Set::Scalar> dist(nc, std::numeric_limits<Set::Scalar>::infinity());
        long n_src = 0;
        if (i1 >= i0)
        {
            const long bx = i1 - i0 + 1, by = j1 - j0 + 1, bz = k1 - k0 + 1;
            std::vector<int> idx(static_cast<std::size_t>(bx * by * bz), -1);
            auto at = [&] (int i, int j, int k) -> long {
                return (static_cast<long>(k - k0) * by + (j - j0)) * bx + (i - i0); };
            for (std::size_t n = 0; n < nc; ++n)
                if (cell_node[n]) idx[static_cast<std::size_t>(at(cell_i[n], cell_j[n], cell_k[n]))]
                    = static_cast<int>(n);

            // 4. the sources: the cell under each source column's top face
            using QE = std::pair<Set::Scalar, int>;
            std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
            for (const auto& sc : src_cols)
            {
                const std::size_t c = sc.first;
                for (int q = fl_col_first[c]; q < fl_col_first[c] + fl_col_count[c]; ++q)
                {
                    const FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                    if (e.cls != FL_TOP || e.blk) continue;
                    const long a = at(e.i, e.j, e.k);
                    const int n = idx[static_cast<std::size_t>(a)];
                    if (n < 0) continue;
                    ++n_src;
                    if (sc.second < dist[static_cast<std::size_t>(n)])
                    {
                        dist[static_cast<std::size_t>(n)] = sc.second;
                        pq.emplace(sc.second, n);
                    }
                }
            }
            // 5. Dijkstra
            static const std::array<Set::Scalar, 17> len = [] {
                std::array<Set::Scalar, 17> l{};
                for (int m = 0; m < 17; ++m) l[static_cast<std::size_t>(m)] = std::sqrt(Set::Scalar(m));
                return l;
            }();
            const auto& OFF = FrPathOffsets();
            while (!pq.empty())
            {
                const QE top = pq.top();
                pq.pop();
                const std::size_t n = static_cast<std::size_t>(top.second);
                if (top.first > dist[n]) continue;
                const int ci = cell_i[n], cj = cell_j[n], ck = cell_k[n];
                for (const auto& o : OFF)
                {
                    const int p = ci + o[0], q = cj + o[1], r = ck + o[2];
                    if (p < i0 || p > i1 || q < j0 || q > j1 || r < k0 || r > k1) continue;
                    const int m = idx[static_cast<std::size_t>(at(p, q, r))];
                    if (m < 0) continue;
                    const Set::Scalar nd = top.first + dx * len[static_cast<std::size_t>(o[3])];
                    if (nd < dist[static_cast<std::size_t>(m)])
                    {
                        dist[static_cast<std::size_t>(m)] = nd;
                        pq.emplace(nd, m);
                    }
                }
            }
        }

        // 6. the heat: FaceH with r := s_path, s := H; blocked / unreached 0
        long n_unr = 0;
        Set::Scalar s_max = 0.0;
        for (std::size_t n = 0; n < nc; ++n)
        {
            for (int q = cell_q0[n]; q < cell_q0[n] + cell_nq[n]; ++q)
            {
                FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                if (e.blk) { e.h = 0.0; continue; }
                if (!std::isfinite(dist[n])) { e.h = 0.0; ++n_unr; continue; }
                e.sp = dist[n];
                s_max = std::max(s_max, e.sp);
                Set::Scalar x, y, z;
                FaceCentre(e, PLO, DX, x, y, z);
                e.h = FaceHr(x, y, e.sp, H, time);
            }
        }
        fl_H_stag      = H;
        fl_n_src       = static_cast<Set::Scalar>(n_src);
        fl_n_blocked   = static_cast<Set::Scalar>(n_blk);
        fl_n_unreached = static_cast<Set::Scalar>(n_unr);
        fl_s_max       = s_max;
        {
            double tp = amrex::second() - t0;
            amrex::ParallelDescriptor::ReduceRealMax(tp);
            fl_t_path = static_cast<Set::Scalar>(tp);
        }
        return H;
    }

    // a listed face's centre (the expressions FaceListHeat uses, token for token)
    static void FaceCentre(const FaceEntry& e, const Set::Scalar* PLO, const Set::Scalar* DX,
                           Set::Scalar& x, Set::Scalar& y, Set::Scalar& z_f)
    {
        Set::Scalar ox = 0.5, oy = 0.5;
        if (e.dir == 0) ox = 0.0; else if (e.dir == 1) ox = 1.0;
        if (e.dir == 2) oy = 0.0; else if (e.dir == 3) oy = 1.0;
        x = (ox == 0.5) ? PLO[0] + (Set::Scalar(e.i) + 0.5) * DX[0]
                        : PLO[0] + (Set::Scalar(e.i) + ox) * DX[0];
        y = (oy == 0.5) ? PLO[1] + (Set::Scalar(e.j) + 0.5) * DX[1]
                        : PLO[1] + (Set::Scalar(e.j) + oy) * DX[1];
        if (e.dir == 5)      z_f = PLO[2] + Set::Scalar(e.k + 1) * DX[2];
        else if (e.dir == 4) z_f = PLO[2] + Set::Scalar(e.k) * DX[2];
        else                 z_f = PLO[2] + (Set::Scalar(e.k) + 0.5) * DX[2];
    }

    // =====================================================================
    // P2-9 rev 2 item 4 -- THE NORMAL. n = -grad F_s / |grad F_s|: F = the
    // solid FRACTION (fr_f), F_s = F smoothed by a normalised Gaussian,
    // sigma = 1 cell truncated at +-2 cells (a fixed 5x5x5 stencil, separable
    // weights), grad by central differences of F_s. P2-8's +-1 indicator
    // gradient read every step edge as a 45 deg corner: +5 to +10 % area at
    // intermediate tilts (rev 1 a2 vi). ONE function: BuildFaceList calls it for
    // every listed cell, and the kernel, ledger, march, R0 and diagnostics read
    // the n and w it leaves in the list.
    // =====================================================================
    static Set::Scalar FrGauss(int a)
    {
        static const Set::Scalar e1 = std::exp(-0.5), e2 = std::exp(-2.0);
        static const Set::Scalar sum = 1.0 + 2.0 * e1 + 2.0 * e2;
        const int b = a < 0 ? -a : a;
        return (b == 0 ? 1.0 : (b == 1 ? e1 : e2)) / sum;
    }
    static Set::Scalar FrSmoothAt(const amrex::Array4<const Set::Scalar>& F, int i, int j, int k)
    {
        Set::Scalar s = 0.0;
        for (int c = -2; c <= 2; ++c)
        for (int b = -2; b <= 2; ++b)
        for (int a = -2; a <= 2; ++a)
            s += FrGauss(a) * FrGauss(b) * FrGauss(c) * F(i + a, j + b, k + c);
        return s;
    }
    static Set::Vector FrSmoothedGradient(const amrex::Array4<const Set::Scalar>& F,
                                         int i, int j, int k, Set::Scalar dx)
    {
        Set::Vector g;
        g(0) = (FrSmoothAt(F, i + 1, j, k) - FrSmoothAt(F, i - 1, j, k)) / (2.0 * dx);
        g(1) = (FrSmoothAt(F, i, j + 1, k) - FrSmoothAt(F, i, j - 1, k)) / (2.0 * dx);
        g(2) = (FrSmoothAt(F, i, j, k + 1) - FrSmoothAt(F, i, j, k - 1)) / (2.0 * dx);
        return g;
    }
