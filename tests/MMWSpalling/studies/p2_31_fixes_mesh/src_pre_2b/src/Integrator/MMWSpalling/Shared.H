#ifndef INTEGRATOR_MMWSPALLING_SHARED_H
#define INTEGRATOR_MMWSPALLING_SHARED_H

// Fragment of class Integrator::MMWSpalling (P2-18, A3). This file is
// textually #included INSIDE the class body of src/Integrator/MMWSpalling.H,
// immediately before FaceRemoval.H. It is a class-body fragment, NOT a standalone
// header: it must not be compiled or #included anywhere else, and it has no
// #include / namespace / class of its own (everything it needs is already
// included by MMWSpalling.H before the class).
    // =====================================================================
    // P2-9: THE DEPTH SCAN'S PIECES AND THE REMAP, ONE COPY EACH.
    // =====================================================================

    // P2-9 rev 2 item 3: ONE copy of the elastic constants every T-derived
    // stress site reads -- phase 0's E, beta, nu = E/(2 mu) - 1, and the
    // confining offset p nu/(1 - nu) (Kant Eq. 15), 0 when confining.p is 0 or
    // nu is outside (-1, 1). nu is 0 unless phase 0 exists and either the rock
    // law is selected or E > 0 and mu > 0.
    struct ElasticConsts { Set::Scalar E = 0.0, beta = 0.0, nu = 0.0, conf = 0.0; };
    ElasticConsts ElasticConstants() const
    {
        ElasticConsts c;
        if (phases.empty()) return c;
        c.E    = phases[0].E;
        c.beta = phases[0].beta;
        // P2-14a (DR-1): ONE Poisson's ratio. Under the rock law it is the law's
        // own (`grimsel.nu`, the nu LawSigma divides by); phase0.E/mu are then
        // inert. Without the law, the scalar pair E/(2 mu) - 1.
        // P2-27 stage 4: under the rock law E / beta / mu are REFUSED when the
        // phase's thermal law is grimsel too (Parse), so they read 0 here and
        // the test below does not ask for them: every consumer takes the law then
        // (LawSigma in ScanSigmaT, the strain table in BuckleRatio).
        const Numeric::Material::Grimsel *law = RockLaw();
        if (!law && !(phases[0].E > 0.0 && phases[0].mu > 0.0)) return c;
        c.nu = law ? law->nu : phases[0].E / (2.0 * phases[0].mu) - 1.0;
        if (confining_p != 0.0 && c.nu > -1.0 && c.nu < 1.0)
            c.conf = confining_p * c.nu / (1.0 - c.nu);
        return c;
    }

    // P2-9 rev 2 item 2: EVERY temperature ghost after a removal write,
    // through the solver's own BC (HeatConduction::bc). The integrator refills
    // ghosts only at the START of a step (Integrator.cpp:1204), and
    // Util::RealFillBoundary is a bare FillBoundary() -- box to box, neither the
    // periodic copy nor the physical BC -- so without this the firing step's
    // kernel would read the PRE-remap temperature across every domain side. This is
    // what FillPatch does at level 0 with source == destination
    // (FillBoundary(periodicity) then the physical BC); FillPatch itself is
    // private to Integrator, the BC object is protected.
    void RefillTempGhosts(int lev, Set::Scalar time)
    {
        HeatConduction::bc->define(geom[lev]);
        HeatConduction::bc->FillBoundary(*temp_mf[lev], 0, temp_mf[lev]->nComp(), time, 0);
    }

    // P2-19: THE DEPTH OF A LINE'S OWN TEMPERATURE RISE, defined on ANY profile
    // (no fit, no face flux, so also under an interior peak). The samples are the
    // face (z = 0, value rise_s) and the cell centres z_d = (d + 1/2) dz,
    // d = 0..nmax, value rise_at_d(d) [K]; z and dz run ALONG the line.
    // ProfileFallDepth: walking down from the sample (z_a, r_a) through
    // d = d0..nmax, the first depth where the rise falls to `target`: ln(rise)
    // linear between the bracketing samples (exact for an exponential at any
    // dz), the rise itself linear only when the lower sample's rise is <= 0.
    // NaN if it never falls that far. Pure; unit/profile_depth compiles and runs
    // this text.
    template <class RISE>
    static Set::Scalar ProfileFallDepth(RISE&& rise_at_d, Set::Scalar z_a, Set::Scalar r_a,
                                        Set::Scalar target, Set::Scalar dz, int d0, int nmax)
    {
        for (int d = d0; d <= nmax; ++d)
        {
            const Set::Scalar z_b = (static_cast<Set::Scalar>(d) + 0.5) * dz;
            const Set::Scalar r_b = rise_at_d(d);
            if (r_b <= target)
            {
                if (r_b > 0.0)
                    return z_a + (z_b - z_a) * (std::log(r_a / target) / std::log(r_a / r_b));
                return z_a + (z_b - z_a) * ((r_a - target) / (r_a - r_b));
            }
            z_a = z_b;
            r_a = r_b;
        }
        return std::numeric_limits<Set::Scalar>::quiet_NaN();
    }

    // The e-fold depth from the FACE: the first z at which the rise falls to
    // rise_s / e. NaN when rise_s <= 0 or it never does for d <= nmax.
    template <class RISE>
    static Set::Scalar ProfileEFoldDepth(RISE&& rise_at_d, Set::Scalar rise_s,
                                         Set::Scalar dz, int nmax)
    {
        if (!(rise_s > 0.0)) return std::numeric_limits<Set::Scalar>::quiet_NaN();
        return ProfileFallDepth(rise_at_d, 0.0, rise_s, rise_s / std::exp(1.0), dz, 0, nmax);
    }

    // The profile's peak: z_peak = the depth of the largest sample (the first, on
    // a tie) and ell_pk = the DEPTH (from the face) at which the rise has fallen
    // to the peak value / e below it. Both NaN when the largest sample is <= 0;
    // ell_pk NaN (z_peak 0) when the peak is the face itself.
    template <class RISE>
    static void ProfilePeak(RISE&& rise_at_d, Set::Scalar rise_s, Set::Scalar dz, int nmax,
                            Set::Scalar& z_peak, Set::Scalar& ell_pk)
    {
        const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
        int d_pk = -1;
        Set::Scalar r_pk = rise_s;
        for (int d = 0; d <= nmax; ++d)
        {
            const Set::Scalar r = rise_at_d(d);
            if (r > r_pk) { r_pk = r; d_pk = d; }
        }
        if (!(r_pk > 0.0)) { z_peak = nan; ell_pk = nan; return; }
        if (d_pk < 0) { z_peak = 0.0; ell_pk = nan; return; }
        z_peak = (static_cast<Set::Scalar>(d_pk) + 0.5) * dz;
        ell_pk = ProfileFallDepth(rise_at_d, z_peak, r_pk, r_pk / std::exp(1.0), dz, d_pk + 1, nmax);
    }

    // T at depth z_v along a line of cells of length dz. Inside the top
    // half-cell: the lerp from the surface temperature (T_face_of(T_top)) to
    // the first cell centre; deeper: linear between cell centres. T_at_didx
    // reads the line's profile by index and clamps.
    template <class TD, class TF>
    static Set::Scalar ScanTAtZv(TD&& T_at_didx, TF&& T_face_of,
                                 Set::Scalar dz, Set::Scalar z_v)
    {
        const Set::Scalar half_dz =
            0.5 * dz;
        Set::Scalar T_z;
        if (z_v <= half_dz)
        {
            const Set::Scalar T_top =
                T_at_didx(0);
            const Set::Scalar T_face = T_face_of(T_top);
            Set::Scalar alpha =
                z_v / half_dz;
            if (alpha < 0.0) alpha = 0.0;
            if (alpha > 1.0) alpha = 1.0;
            T_z = (1.0 - alpha) * T_face
                +        alpha  * T_top;
        }
        else
        {
            const Set::Scalar u =
                z_v / dz - 0.5;
            const Set::Scalar uf =
                std::floor(u);
            int n = static_cast<int>(uf);
            if (n < 0) n = 0;
            const Set::Scalar alpha =
                u - uf;  // in [0, 1)
            T_z = (1.0 - alpha)
                  * T_at_didx(n)
                +        alpha
                  * T_at_didx(n + 1);
        }
        return T_z;
    }

    // K_Ic(T) under the configured law, in Pa*m^0.5. ⚠ P2-2d d7: the depth
    // scan's two reads (at a cell, at a depth) are this ONE function, so they
    // cannot disagree.
    Set::Scalar ScanKIcPa(Set::Scalar T_at) const
    {
        const Numeric::Material::Grimsel *kic_g_loc =
            phases.empty() ? nullptr : &phases[0].grimsel;
        Set::Scalar K_Ic_MPa;
        if (kic_law == KIcLaw::GrimselTable && kic_g_loc)
            K_Ic_MPa = kic_scale
                     * kic_g_loc->K_Ic_of(T_at);
        else if (kic_law == KIcLaw::NasseriTable)
            K_Ic_MPa =
                Numeric::SpCriterion::KIcNasseri(
                    T_at, kic_scale);
        else
            K_Ic_MPa =
                Numeric::SpCriterion::KIcLinear(
                    T_at, kic_K_Ic0,
                    kic_T_ref, kic_T_melt);
        return Numeric::SpCriterion
            ::MPa_sqrtm_to_Pa_sqrtm(K_Ic_MPa);
    }

    // The T-derived driving stress at one sample (1-D confinement,
    // compression-positive): the rock law through LawSigma when one is
    // selected, else the scalar thermoelastic form.
    Set::Scalar ScanSigmaT(const Numeric::Material::Grimsel *scan_g,
                           Set::Scalar T_z, Set::Scalar T0_z,
                           Set::Scalar sp_conf_offset, Set::Scalar phase_E_loc,
                           Set::Scalar phase_beta_loc, Set::Scalar phase_nu_loc) const
    {
        if (scan_g)
            return sp_relief
                   * LawSigma(scan_g, T_z, T0_z)
                   + sp_conf_offset;
        // ⚠ pre-P2-2b expression, VERBATIM.
        return sp_relief * phase_E_loc
               * phase_beta_loc
               * (T_z - T0_z)
               / (1.0 - phase_nu_loc)
               + sp_conf_offset;
    }

    // The depth scan proper: z_crack walks the line in 1-cell strides and
    // stops at the first K_I < K_Ic, after max_steps, or below scan_kbot;
    // the bracketing cell is then bisected 8 times. Returns h_col (0 = no
    // crack below the surface). `on_pass(k_tip)` sees every passing stride;
    // `s_last_pass` / `found_failure` say how the walk ended; `score0` is
    // K_I/K_Ic at the surface (s = 0) -- the line's score (R2). ⚠ K_Ic is
    // evaluated BEFORE K_I at every sample.
    template <class KI, class KICZ, class KICC, class PASS>
    static Set::Scalar DepthScanH(int max_steps, int k_top, int scan_kbot,
                                  Set::Scalar dz, Set::Scalar a_f, bool rec_ok,
                                  KI&& k_i_at_zcrack, KICZ&& kic_pa_at_z,
                                  KICC&& kic_pa_at_cell, PASS&& on_pass,
                                  int& s_last_pass, bool& found_failure,
                                  Set::Scalar& score0)
    {
        Set::Scalar h_col = 0.0;
        s_last_pass = -1;
        found_failure = false;
        score0 = 0.0;
        for (int s = 0; s <= max_steps; ++s)
        {
            const Set::Scalar z_crack =
                Set::Scalar(s) * dz;
            const int k_tip = k_top - s;
            if (k_tip < scan_kbot) break;

            Set::Scalar K_Ic_Pa = rec_ok
                ? kic_pa_at_z(z_crack + a_f)
                : kic_pa_at_cell(k_tip);
            Set::Scalar K_I_val =
                k_i_at_zcrack(z_crack);
            if (s == 0 && K_Ic_Pa > 0.0) score0 = K_I_val / K_Ic_Pa;

            if (K_I_val >= K_Ic_Pa)
            {
                s_last_pass = s;
                h_col       = z_crack;
                on_pass(k_tip);
            }
            else
            {
                found_failure = true;
                break;
            }
        }
        if (found_failure && s_last_pass >= 0 &&
            s_last_pass < max_steps &&
            (k_top - (s_last_pass + 1)) >= scan_kbot)
        {
            Set::Scalar z_lo = Set::Scalar(s_last_pass)
                               * dz; // pass
            Set::Scalar z_hi = Set::Scalar(s_last_pass + 1)
                               * dz; // fail
            for (int b = 0; b < 8; ++b)
            {
                Set::Scalar z_mid =
                    0.5 * (z_lo + z_hi);
                int k_mid = k_top -
                    static_cast<int>(
                        std::floor(z_mid / dz));
                if (k_mid < scan_kbot) k_mid = scan_kbot;
                if (k_mid > k_top)      k_mid = k_top;
                Set::Scalar K_Ic_mid = rec_ok
                    ? kic_pa_at_z(z_mid + a_f)
                    : kic_pa_at_cell(k_mid);
                Set::Scalar K_I_mid =
                    k_i_at_zcrack(z_mid);
                if (K_I_mid >= K_Ic_mid)
                    z_lo = z_mid;
                else
                    z_hi = z_mid;
            }
            h_col = z_lo;
        }
        return h_col;
    }

    // P0b's conservative remap of one line's enthalpy profile by a flake h_f
    // along the line (called by RecedeAlongLine, FaceRemoval.H):
    //     new[d] = (1 - b) old[d + m] + b old[d + m + 1],
    //     b = (1 - e^{-r/ell}) / (1 - e^{-dz/ell}),  h_f = m dz + r,
    // indices past `nmax` clamped to it (ambient at depth on a long line).
    // ell = infinity (a flat line) gives b = r/dz. Returns the telescoped
    // sum of (old - new) over 0..nmax, per unit volume.
    static Set::Scalar RemapProfile(const Set::Scalar* Ho, int nmax,
                                    Set::Scalar h_f, Set::Scalar ell,
                                    Set::Scalar dz, Set::Scalar* Hn)
    {
        const Set::Scalar x = dz / ell;
        const int m = static_cast<int>(std::floor(h_f / dz));
        const Set::Scalar r = h_f - Set::Scalar(m) * dz;
        // b in [0, 1): the fraction of the NEXT cell that the shift
        // pulls in, weighted by the exponential sub-cell profile.
        const Set::Scalar den = 1.0 - std::exp(-x);
        const Set::Scalar b = (den > 0.0)
            ? (1.0 - std::exp(-r / ell)) / den : (r / dz);
        auto Hat = [&] (int d) -> Set::Scalar {
            if (d < 0)     d = 0;
            if (d > nmax)  d = nmax;
            return Ho[d];
        };
        for (int d = 0; d <= nmax; ++d)
            Hn[d] = (1.0 - b) * Hat(d + m) + b * Hat(d + m + 1);
        // The energy the flake carries away: the telescoped difference,
        // which IS rho*c*h_f*(T_bar_f - T0) for an exponential profile.
        Set::Scalar dE = 0.0;
        for (int d = 0; d <= nmax; ++d) dE += Ho[d] - Hn[d];
        return dE;
    }

    // =====================================================================
    // P2-30: ONE COLUMN LIST PER LOG FILE. The ordered column NAMES of each CSV this class writes live here; the header is
    // built from the list, and every row writer counts its own fields against it before writing (it assembles the row in
    // a local stream carrying the writer's own format, counts, aborts on a mismatch, then writes). A row writer keeps its
    // own << chain and formats (a column can print as `0` in one row type and `0.0` in another). This catches a missing
    // or an extra field, NOT two fields swapped.
    // =====================================================================
    static std::string CsvHeader(const std::vector<std::string>& cols)
    {
        std::string s;
        for (std::size_t i = 0; i < cols.size(); ++i)
        {
            if (i) s += ",";
            s += cols[i];
        }
        return s + "\n";
    }

    // `row` is one finished row (its trailing newline is not a field); every field is a number, so a comma is a separator
    static void CsvCheckRow(const std::string& row, std::size_t n_cols, const char* file)
    {
        std::size_t n = 1;
        for (char ch : row) if (ch == ',') ++n;
        if (n != n_cols)
            Util::Abort(INFO, file, ": a row has ", n, " fields, its header has ", n_cols, ": ", row);
    }

    // <plot_file>_removal_events.csv: every row writer under the face list (firing, vapour, feet) shares this one list
    static std::vector<std::string> FrEventsColumns(bool bk_on)
    {
        std::vector<std::string> c = {"time", "col_i", "col_j", "regime", "k_top", "T_top", "a_f", "Sp_top",
                                      "h_scan", "h_applied", "n_voided",
                                      "T_surf", "ell_eff", "h_f", "T_bar_flake",
                                      "rc_flag", "d_ell_rel", "ell_raw", "delta", "q_ratio",
                                      // P2-9 deliverable 6, appended: the line
                                      "dir", "mu", "k", "run", "h_line", "h_n", "f_after", "score", "ts_fb", "sheet",
                                      "through", "q_net", "k_c"};
        if (bk_on) c.push_back("bk");                      // P2-15: under spall.buckle_span > 0 only, the line's ratio
        // P2-19, appended: the line's profile depths and the energy thin limit [m] (diagnostics; nothing reads them)
        for (const char* s : {"ell_p", "z_peak", "ell_pk", "h_E"}) c.push_back(s);
        return c;
    }

    // <plot_file>_buckle_log.csv
    static std::vector<std::string> FrBuckleColumns()
    {
        return {"time", "col_i", "col_j", "class", "fire", "through", "fb", "nodepth", "score", "h_col",
                "T_surf", "ell_p", "z_peak", "ell_pk", "ell_eff", "x", "bk", "bk_lb", "blocked"};
    }

    // <plot_file>_jet_profile.csv (the optional columns follow the march's own switches)
    static std::vector<std::string> JetProfileColumns(bool expo, bool fsurf, bool count)
    {
        std::vector<std::string> c = {"time", "bin", "r_lo", "n_cols", "T_gas", "P_bin", "P_share"};
        if (expo) c.push_back("T_eff");
        if (fsurf) { c.push_back("n_fs"); c.push_back("m_ratio"); }
        if (count) c.push_back("P_neg");
        return c;
    }

    // <plot_file>_side_power.csv
    static std::vector<std::string> SidePowerColumns()
    {
        return {"time", "col_i", "col_j", "r", "P_side"};
    }

    // P2-30 (SR 14): the surface patch's three constants, ONE definition, read by the kernel and by EnergyLedger.
    // ⚠ `sp_on` is NOT here: the kernel's carries `owner` and the ledger's does not (equal on one level only).
    struct PatchConst { Set::Scalar inv_dz2; bool presc; Set::Scalar r2; };
    PatchConst PatchConstants(const Set::Scalar* DX) const
    {
        PatchConst p;
        p.inv_dz2 = 1.0 / (DX[2] * DX[2]);
        p.presc   = surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        p.r2      = surface_patch.radius * surface_patch.radius;
        return p;
    }

    // =====================================================================
    // P2-18 A3: THE STEP PROLOGUE of the removal step: RemovalStep
    // (FaceRemoval.H) calls both functions before the lines are read.
    // =====================================================================
    // The guards: removal on (P2-28: or the face list's evaluate-only mode, whose
    // criterion runs here and removes nothing) and the surface-owning level only
    // (AMR: coarser levels receive phi/removed by average-down). Returns false
    // when the step has nothing to do on this level.
    bool RemovalStepBegins(int lev)
    {
        if (!list_state()) return false;
        if (!SurfaceOwner(lev)) return false;
        return true;
    }

    // Per-step diagnostic fields are reset here. `phi` and `removed` carry
    // state across steps and are not reset. The face list writes these on the
    // cells it voids (FaceRemoval.H) and relies on this reset.
    void ResetRemovalStepFields(int lev)
    {
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
        // Step 17: per-step reset for the per-cell regime field under
        // sp_weibull. The field is registered only when
        // spallation_model == SpWeibull && (spall_enabled || vapor_enabled).
        if (spallation_model == SpallationModel::SpWeibull &&
            (spall_enabled || vapor_enabled))
        {
            regime_field_mf[lev]->setVal(0.0);
        }
    }

#endif // INTEGRATOR_MMWSPALLING_SHARED_H
