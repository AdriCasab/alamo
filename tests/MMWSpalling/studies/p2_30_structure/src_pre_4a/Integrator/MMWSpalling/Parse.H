    static void Parse(MMWSpalling& value, IO::ParmParse& pp)
    {
        // P2-10: THE RETIRED KEYS. The flame wall machinery the face list made
        // redundant -- the side-face heat family, the lateral criterion (regime 6,
        // banking, the thermal shift), per-cell removal and the wall march with
        // its cold body -- is deleted. A deck that still names one of its keys
        // aborts here, naming the key, instead of running with it ignored
        // (`spall.*` sits outside any strict ParmParse scope, so this list is
        // the only guard those keys have).
        for (const char* k : {
                                "spall.lateral_drive",
                                "spall.lateral_removal",
                                "spall.lateral_subcell",
                                "spall.lateral_thermal",
                                "spall.per_cell_removal",
                                "surface_patch.side_face_flux",
                                "surface_patch.side_face_h",
                                "surface_patch.side_face_T",
                                "surface_patch.side_face_h_r_min",
                                "surface_patch.side_face_h_law",
                                "surface_patch.side_face_h_datum",
                                "surface_patch.side_face_factor_above",
                                "surface_patch.side_face_area_law",
                                "surface_patch.side_face_step_height",
                                "surface_patch.side_face_staircase",
                                "surface_patch.step_diag_radius",
                                "surface_patch.jet_wall_march",
                                "surface_patch.jet_wall_mdot",
                                "surface_patch.jet_body",
                                "surface_patch.jet_body_T",
                                "surface_patch.jet_body_h_scale"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-10 (the flame wall machinery: "
                            "side-face heat, the lateral criterion, per-cell removal and "
                            "the wall march were deleted); remove it from the input");
        // P2-16: THE CONSTRAINT DEPTH (P2-2e, carried onto the face list in P2-14a) is
        // retired: its onset climbs toward the minimum flux for every d, which Kant 2016
        // Fig. 11 does not show, and spall.buckle_span (P2-15) is the minimum-flux
        // mechanism kept. Same guard as above: `spallation.*` is outside any strict
        // scope, so an unread key would otherwise run silently ignored.
        for (const char* k : {"spallation.sp.constraint_depth"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-16 (the constraint depth was "
                            "deleted; spall.buckle_span carries the minimum flux); remove "
                            "it from the input");
        // P2-20: THE TYPED BEAM ABSORPTION AND REFLECTION are retired: one dielectric
        // description per rock (beam.eps_r, beam.tan_delta) now gives both, alpha =
        // 2 pi sqrt(eps_r) tan_delta / lambda and R = ((sqrt(eps_r) - 1)/(sqrt(eps_r) + 1))^2
        // (MMWBeam.H). Checked HERE, before the beam parse, so the message names the
        // replacement rather than the strict scope's generic "specified but not used",
        // and with the beam off (P0 = 0) as well.
        for (const char* k : {"beam.alpha_expr", "beam.R"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-20 (the typed beam absorption and "
                            "reflection): declare beam.eps_r and beam.tan_delta instead, from "
                            "which alpha and R are computed; remove it from the input");
        // P2-27 (census C1a): THE DAMAGE LINEAGE is retired -- the Drucker-Prager
        // damage law, the grain-boundary cohesive zone and the damage-law
        // removal (the C_h sqrt(alpha t) thickness fired by D >= threshold and
        // the principal-stress ratio; decision 1 of 2026-11-05: LEFM is THE
        // model). Every key of the two families aborts, whatever its value, and
        // so do the damage-law removal's spall.* keys; damage.enabled,
        // conductivity.kappa_damage_alpha and spall.sample are retired AT A
        // DEFAULT (the table below), so a kept deck naming them at the old
        // default still runs. Checked HERE, before any key is queried, so a
        // family entry can only come from the input.
        for (const char* k : {"spall.alpha", "spall.C_h", "spall.t_spall",
                              "spall.damage_threshold", "spall.xi_s", "spall.shallow_depth",
                              "spall.prescribed_ratio.enabled", "spall.prescribed_ratio.expr"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-27 (the damage-law removal: the "
                            "Drucker-Prager damage lineage is deleted, the sp_weibull LEFM "
                            "criterion is the one removal model); remove it from the input");
        for (const char* family : {"damage", "cohesive"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                if (k != "damage.enabled")
                    Util::Abort(INFO, k, " was retired in P2-27 (the Drucker-Prager damage "
                                "law and the grain-boundary cohesive zone are deleted: LEFM "
                                "is the model); remove it from the input");
        // P2-28 (census H): THE COLUMN REMOVAL PATH is retired -- one removal path, the face list.
        // Its own keys abort by PRESENCE, whatever the value: the bit (regime 3), the column CSVs,
        // the cluster labeller's area, the AMR surface band and site A's verification hooks.
        for (const char* k : {"spall.clusters_csv", "spall.h_col_events_csv", "weibull.A_crit",
                              "surface.refine_depth", "surface.refine_height"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-28 (the column removal path is deleted: "
                            "the face list is the one removal path); remove it from the input");
        for (const char* family : {"bit", "spallation.sp.prescribed_stress",
                                   "spallation.sp.prescribed_T_surface"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                Util::Abort(INFO, k, " was retired in P2-28 (the column removal path is deleted: "
                            "the face list is the one removal path); remove it from the input");
        // P2-29 (census D): THE SEPARATE HEAT KERNELS are retired -- one heat path, AdvanceMicrostructure. The
        // global material (`material.*`), its implicit solver (`implicit.*`) and the const-alpha closed form's
        // phase scalars abort by PRESENCE; a deck states its rock as one phase.
        for (const char* k : {"phase.T_ref", "phase.T_m", "phase.T_v", "phase.lh_m", "phase.lh_v"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-29 (the const-alpha path is deleted: one heat path); state "
                            "the rock as one phase -- microstructure.nphases = 1, microstructure.phase0.{T_m, "
                            "fus_width, L_m, ...}; remove it from the input");
        for (const char* family : {"material", "implicit"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                Util::Abort(INFO, k, " was retired in P2-29 (the separate heat kernels and the implicit solver are "
                            "deleted: one heat path); state the rock as one phase -- microstructure.nphases = 1, "
                            "microstructure.phase0.law = constant | grimsel | zhang; remove it from the input");
        // P2-29 (census D): a deck without microstructure.nphases ran a separate kernel; there is none now
        if (!pp.contains("microstructure.nphases"))
            Util::Abort(INFO, "microstructure.nphases is required since P2-29 (one heat path: the const-alpha and "
                        "material kernels are deleted): state the rock as one phase -- microstructure.mode = "
                        "expression, microstructure.phase_expr = 0, microstructure.nphases = 1, "
                        "microstructure.phase0.law = constant | grimsel | zhang");
        // P2-18 A2: KEYS RETIRED AT A DEFAULT. A deleted parse is not a guard: a
        // top-level key (`spall.*`, `microstructure.*`, ...) is then silently
        // ignored, and a strict-scope key (`surface_patch.*`, `beam.*`, the
        // per-phase `grimsel.*`) aborts as "specified but not used" whatever its
        // value. Kept decks name retiring keys at their old defaults far more
        // often than at live values (census 2026-11-05 §4 A2), so each entry
        // here is read ONCE, which marks it used for the strict scopes, and
        // compared with the value the deleted code accepted: equal (or "*" = the
        // key was never read, any value is correct) -> the deck runs unchanged;
        // anything else -> abort naming the key, the value and the packet that
        // retired it. Deleting the parse of a key with a single old default ADDS
        // A ROW HERE; it does not go in the abort list above.
        struct RetiredAtDefault { const char* key; const char* value; const char* note; };
        static const RetiredAtDefault retired_at_default[] = {
            {"surface_patch.face",  "zhi",
             "P2-18: consumed, checked == zhi, read by nothing (only the +z patch exists)"},
            {"surface_patch.T_ext", "*",
             "P2-18: parsed into a member nothing read"},
            {"surface_patch.type",  "prescribed_temperature",
             "P2-18: the legacy Step 6d key; surface_patch.mode decides"},
            // P2-26 (census B3): closed study switches; their strict-scope
            // siblings (surface_patch.*) are simply unparsed and abort there
            {"spall.flake_length",   "twocell",
             "P2-26: the flux length closed at P1g-3 (rejected on provenance); the two-cell length is the only one"},
            {"spall.flake_fraction", "cf",
             "P2-26: the cone-median fraction closed at P1g-3; spall.flake_c_f is the only fraction"},
            {"weibull.size_volume",  "cell",
             "P2-26: the flux stressed volume closed at P2-2d (no target); the flaw is the per-cell draw"},
            {"weibull.V_ref",        "*",
             "P2-26: read only under the deleted size_volume = flux"},
            // P2-26 (census C2, C3): closed, and deleted with one declared re-pin.
            // Numeric rows compare as doubles (P2-27): 0.0 and 0 are the same default.
            {"spall.surface_normal",          "0",
             "P2-26: the Step 20 surface-normal frame; its recession half was superseded at D2z-1 (user Q7: retire)"},
            {"spall.flake_coherence_length",  "0",
             "P2-26: the column path's coherence cap, a needle regulariser of the vertical-only law"},
            {"spall.melt_aware",              "0",
             "P2-26: the P1 melt skeleton f(Lambda_L) (user Q6: delete); the face list's melt rule is P2-24's"},
            {"melt.lambda_crit",              "*",
             "P2-26: read only by the deleted melt skeleton"},
            // P2-27 (census C1a): the damage lineage, retired at its defaults
            {"damage.enabled",                  "0",
             "P2-27: the Drucker-Prager damage law is deleted (LEFM is the model)"},
            {"conductivity.kappa_damage_alpha", "0.0",
             "P2-27: the damage factor exp(-a D) on the conductivity is deleted with D's evolution"},
            {"spall.sample",                    "top_cell",
             "P2-27: the damage-law removal's D-sampling window is deleted with it"},
            // P2-28 (census H): the column path's selectors, each retired at the ONE value the face list runs
            {"spall.surface_T",      "reconstructed",
             "P2-28: the face list's R1 IS the reconstruction (site A's cell_average / diagnostic forms are deleted)"},
            {"spall.removal_energy", "flake",
             "P2-28: the face list's R3 IS the flake remap (the column path's whole-cell removal is deleted)"},
            {"weibull.detach_mode",  "per_face",
             "P2-28: the other mode was the cluster labeller, deleted with site A"},
            {"surface.follow_mask",  "1",
             "P2-28: the face list's `top` class IS the mask column top, and removal requires the list"},
            // P2-29 (census D): the separate kernels' selectors, each retired at the value the one heat path runs
            {"solver",                  "explicit",
             "P2-29: the implicit solver is deleted; the one heat path is forward Euler"},
            {"alpha",                   "0 0",
             "P2-29: the homogeneous eigenstrain is deleted; a phase's beta is the thermal expansion"},
            {"surface_patch.foot_rule", "pads",
             "P2-29: the max / mean feet rules were closed (no live deck); pads is the feet rule"},
        };
        // P2-27: a ONE-token value that parses fully as a double is compared as
        // a double, so `0.0` equals the table's "0" (P2-26 DR-3: the string
        // compare refused 622 record argvs that wrote a default as 0.0).
        // P2-29: TOKEN-WISE -- several tokens compare as doubles too, the same
        // count and each equal, so `0.0 0.0` equals the table's "0 0"; a word,
        // or a different count, is compared as a string (and differs).
        auto doubles = [] (const std::string& s, std::vector<double>& d) {
            std::istringstream in(s);
            std::string t;
            d.clear();
            while (in >> t)
            {
                char* end = nullptr;
                const double x = std::strtod(t.c_str(), &end);
                if (end != t.c_str() + t.size()) return false;
                d.push_back(x);
            }
            return !d.empty();
        };
        for (const RetiredAtDefault& r : retired_at_default)
        {
            if (!pp.contains(r.key)) continue;
            std::vector<std::string> toks;
            static_cast<amrex::ParmParse&>(pp).queryarr(r.key, toks);   // marks the entry used
            std::string got;
            for (const std::string& t : toks) got += (got.empty() ? "" : " ") + t;
            if (std::string(r.value) == "*" || got == r.value) continue;
            std::vector<double> d_got, d_old;
            if (doubles(got, d_got) && doubles(r.value, d_old) && d_got == d_old) continue;
            Util::Abort(INFO, r.key, " = ", got, " was retired at its default (", r.value,
                        ") -- ", r.note, "; remove it from the input or restore the code");
        }
        // Reuse HeatConduction's parse: hc.heat.alpha (diffusivity κ/(ρ c_p)),
        // hc.heat.refinement_threshold, hc.ic.*, hc.bc.temp.* (BC on T).
        // The "Temp" / "Temp_old" fields it registers serve as our T field.
        pp.queryclass<HeatConduction>("hc", value);
        pp.queryclass<Mechanics<Model::Solid::Affine::Isotropic>>("el", value);
        pp.queryclass<Numeric::MMWBeam>("beam", value.beam);

        // E1: beam-deposition closure counter (see BeamEnergyBalance).
        // Non-extensive: the Integrator neither zeroes nor MPI-sums these.
        value.RegisterIntegratedVariable(&value.beam_P_inc,          "beam_P_inc",          false);
        value.RegisterIntegratedVariable(&value.beam_P_dep,          "beam_P_dep",          false);
        value.RegisterIntegratedVariable(&value.beam_P_trans,        "beam_P_trans",        false);
        value.RegisterIntegratedVariable(&value.beam_E_inc,          "beam_E_inc",          false);
        value.RegisterIntegratedVariable(&value.beam_E_dep,          "beam_E_dep",          false);
        value.RegisterIntegratedVariable(&value.beam_E_trans,        "beam_E_trans",        false);
        value.RegisterIntegratedVariable(&value.beam_closure_err,    "beam_closure_err",    false);
        value.RegisterIntegratedVariable(&value.beam_invariant_viol, "beam_invariant_viol", false);
        // P2-21: the optical-depth deciders, only on a deck whose loss table
        // varies AND whose beam can deposit (a beam-OFF deck that includes the
        // rock file, e.g. the face gate's HOS_S and CAL_S, keeps its header).
        // P2-24: the melt keys reach the deposition only through the column optical-depth sum
        if (value.beam.melt_active() && value.beam.can_deposit() && !value.beam.alpha_varies())
            Util::Abort(INFO, "beam.melt_factor / beam.R_melt (P2-24) act only through the column optical-depth "
                        "sum: they need a varying beam.tan_delta_T (two equal pairs give a flat one); got "
                        "melt_factor ", value.beam.melt_factor, ", R_melt ", value.beam.R_melt,
                        " with a flat table");
        if (value.beam.alpha_varies() && value.beam.can_deposit())
        {
            value.RegisterIntegratedVariable(&value.beam_tau_max,     "beam_tau_max",     false);
            value.RegisterIntegratedVariable(&value.beam_alpha_max,   "beam_alpha_max",   false);
            value.RegisterIntegratedVariable(&value.beam_col_err_max, "beam_col_err_max", false);
            value.RegisterIntegratedVariable(&value.bc_T_peak,        "bc_T_peak",        false);
            value.RegisterIntegratedVariable(&value.bc_z_peak,        "bc_z_peak",        false);
            value.RegisterIntegratedVariable(&value.bc_z_melt,        "bc_z_melt",        false);
            value.RegisterIntegratedVariable(&value.bc_tau,           "bc_tau",           false);
            value.RegisterIntegratedVariable(&value.bc_ell_p,         "bc_ell_p",         false);
            value.RegisterIntegratedVariable(&value.bc_ell_pk,        "bc_ell_pk",        false);
            // P2-24: the melt and vapour record, appended
            value.RegisterIntegratedVariable(&value.beam_P_beam,    "beam_P_beam",    false);
            value.RegisterIntegratedVariable(&value.bc_T_rad,       "bc_T_rad",       false);
            value.RegisterIntegratedVariable(&value.bc_T_view,      "bc_T_view",      false);
            value.RegisterIntegratedVariable(&value.bc_A_melt,      "bc_A_melt",      false);
            value.RegisterIntegratedVariable(&value.bc_f_L_peak,    "bc_f_L_peak",    false);
            value.RegisterIntegratedVariable(&value.bc_z_liq,       "bc_z_liq",       false);
            value.RegisterIntegratedVariable(&value.bc_R_col,       "bc_R_col",       false);
            value.RegisterIntegratedVariable(&value.bc_n_vapor,     "bc_n_vapor",     false);
            value.RegisterIntegratedVariable(&value.bc_n_vapor_int, "bc_n_vapor_int", false);
        }

        // Surface-loss parameters (eq. 16). Defaults disable both losses, so
        // step 3 regression is unchanged.
        pp.query_default("losses.epsilon",      value.eps_rad,      0.0);
        pp.query_default("losses.sigma_SB",     value.sigma_SB,     0.0);
        pp.query_default("losses.h_conv",       value.h_conv,       0.0);
        pp.query_default("losses.T_amb",        value.T_amb,        0.0);
        // Optional T-dependent emissivity: `losses.epsilon` is used for
        // T <= T_eps_break, `losses.epsilon_high` for T > T_eps_break. With
        // T_eps_break left at +∞ both default to a single emissivity.
        pp.query_default("losses.epsilon_high", value.eps_rad_hi,   value.eps_rad);
        pp.query_default("losses.T_eps_break",  value.T_eps_break,  1.0e30);

        // Step 6d: optional prescribed-temperature surface patch (Hu LRST
        // conduction validation). When `surface_patch.enabled = 1`, a
        // circular patch on a single domain face injects a Dirichlet-like
        // flux  kappa * (T_patch(t) - T_cell) / dz^2  on in-patch top-row
        // cells. Default disabled keeps every prior MMWSpalling test
        // unchanged. Currently only the explicit material path consumes
        // this term; the implicit path hardcodes Neumann domain BCs.
        pp.queryclass<SurfacePatch>("surface_patch", value.surface_patch);

        // Step 10: optional surface advancement and spall detachment.
        // Disabled by default so existing static-surface behavior is
        // preserved. When enabled, registers a level-set `phi`, a moving
        // `removed` mask, and per-step event diagnostics, and reuses the
        // microstructure path for resetting detached cells.
        ParseSpallSettings(value, pp);

        // Step 11: optional vaporisation removal. Default-off and shares
        // phi / removed with the spall path. When enabled, registers
        // RoP_vap, vapor_event, vapor_thickness, regime, and unified RoP
        // diagnostics.
        ParseVaporSettings(value, pp);

        // Step 12: spallation-model toggle. `spallation.model = sp_weibull`
        // loads the LEFM Sp = K_I/K_Ic(T) onset criterion; the default `none`
        // loads no criterion (a thermal-only deck). P2-27: `damage_law` is
        // retired with the damage lineage and aborts by value.
        ParseSpallationModel(value, pp);

        // Register conserved enthalpy H and phase fractions. No BC needed
        // (these have zero ghost cells); RegisterNewFab passes nullptr to
        // AddField, which auto-substitutes the framework's bcnothing.
        value.RegisterNewFab(value.H_mf,        nullptr, 1, 0, "H",        true);
        value.RegisterNewFab(value.H_old_mf,    nullptr, 1, 0, "H_old",    false);
        value.RegisterNewFab(value.Lambda_S_mf, nullptr, 1, 0, "Lambda_S", true);
        value.RegisterNewFab(value.Lambda_L_mf, nullptr, 1, 0, "Lambda_L", true);
        value.RegisterNewFab(value.Lambda_V_mf, nullptr, 1, 0, "Lambda_V", true);

        // Surface mask: 1.0 on the topmost z-row of each level, 0.0 elsewhere.
        // Step 11 will repurpose this as a level-set-driven moving surface.
        value.RegisterNewFab(value.surface_mf, nullptr, 1, 0, "surface", true);

        // Step 5/6: mineral microstructure. Activated when
        // `microstructure.nphases` is present, so existing tests (Stefan, beam,
        // equilibrium, Zhang) keep their homogeneous behaviour. Each
        // `microstructure.phaseN.*` block defines one mineral phase with a
        // target volume fraction and the per-phase thermophysical / mechanical
        // properties used by later steps (heterogeneous κ, ρ, Cp, E for steps 6
        // and 7; β, A_D for step 8; α_attenuation for the Beer-Lambert source).
        if (pp.contains("microstructure.nphases"))
        {
            // Phase-id source: random Voronoi tessellation (default — Step 5)
            // or a user-supplied amrex::Parser expression of (x,y,z) (added
            // in Step 6 so the heterogeneous-kappa test can prescribe a
            // deterministic laminate without random seeding).
            std::string mode_str;
            pp.query_default("microstructure.mode", mode_str, std::string("voronoi"));
            if      (mode_str == "voronoi")    value.micro_mode = MicroMode::Voronoi;
            else if (mode_str == "expression") value.micro_mode = MicroMode::Expression;
            else Util::Abort(INFO, "microstructure.mode must be 'voronoi' or 'expression', got: ", mode_str);

            if (value.micro_mode == MicroMode::Voronoi)
            {
                pp.query_required("microstructure.number_of_grains", value.num_grains);
                pp.query_default ("microstructure.seed",             value.micro_seed, 1);
            }
            else // Expression
            {
                std::string expr;
                pp.query_required("microstructure.phase_expr", expr);
                value.phase_parser = amrex::Parser(expr);
                value.phase_parser.registerVariables({"x", "y", "z"});
                value.phase_f = value.phase_parser.compile<3>();

                // Step 6c: optional independent grain-id expression so tests
                // can prescribe same-mineral grain neighbours. Absent ⇒
                // grain_id = phase (preserves Step 6/6b expression inputs).
                std::string gexpr;
                if (pp.contains("microstructure.grain_expr"))
                {
                    pp.query_required("microstructure.grain_expr", gexpr);
                    value.grain_parser = amrex::Parser(gexpr);
                    value.grain_parser.registerVariables({"x", "y", "z"});
                    value.grain_f = value.grain_parser.compile<3>();
                    value.has_grain_expr = true;
                }
            }

            int nphases = 0;
            pp.query_required("microstructure.nphases", nphases);
            Util::Assert(INFO, TEST(nphases >= 1));
            value.phases.resize(nphases);
            for (int i = 0; i < nphases; ++i)
            {
                std::string pfx = "microstructure.phase" + std::to_string(i);
                pp.query_required(pfx + ".fraction", value.phases[i].fraction);
                // P2-27 (census C1a): A_D was the damage law's growth rate; the
                // law is deleted, so the key is optional and refused unless 0.0
                {
                    Set::Scalar A_D = 0.0;
                    pp.query_default(pfx + ".A_D", A_D, 0.0);
                    if (A_D != 0.0)
                        Util::Abort(INFO, pfx, ".A_D = ", A_D, " was retired in P2-27 (the "
                                    "Drucker-Prager damage law it fed is deleted); remove it or set 0.0");
                }
                Util::Assert(INFO, TEST(value.phases[i].fraction >= 0.0));

                // P2-2 change 2. Default `constant` == the shipped behaviour,
                // so an input that does not name this key is byte-identical.
                {
                    std::string plaw;
                    pp.query_default(pfx + ".law", plaw, std::string("constant"));
                    // ⚠⚠ PARSED UNCONDITIONALLY, AND THAT IS THE POINT.
                    // The law's coefficients are read whether or not THIS key
                    // turns the thermal half on, because P2-2b's mechanical
                    // half (`spall.stress_law = grimsel`) consults the same
                    // object. When this was inside the `grimsel` branch, a leg
                    // running the mechanical set alone -- which is exactly how
                    // the Kant ladder separates conduction from mechanics --
                    // got a DEFAULT-CONSTRUCTED law, so
                    // `<pfx>.grimsel.alpha_cap_theta` was silently ignored and
                    // (c7)'s two alpha fits would have come back IDENTICAL,
                    // reporting no sensitivity where the key had simply never
                    // been read. Every coefficient is defaulted to the
                    // reference's recommendation, so parsing it always is
                    // free, and `queryclass` consuming the keys is what keeps
                    // a mistyped one from passing unnoticed.
                    pp.queryclass<Numeric::Material::Grimsel>(
                        pfx + ".grimsel", value.phases[i].grimsel);
                    // P2-29: Zhang's granite (Zhang Table 1, Numeric::Material::Zhang) as a phase law,
                    // parsed unconditionally as grimsel's is. Its latent scalars are the PHASE's (the
                    // table copy below carries them, the grimsel rule), so the law's own copies under
                    // `<pfx>.zhang.` would be read and then overwritten -- refused instead.
                    for (const char* lk : {".T_ref", ".T_m", ".fus_width", ".L_m",
                                           ".T_vap_lo", ".T_vap_hi", ".L_v"})
                        if (pp.contains(pfx + ".zhang" + lk))
                            Util::Abort(INFO, pfx, ".zhang", lk, " is not read: the phase's own ", pfx, lk,
                                        " carries the latent scalars of every law (P2-29)");
                    pp.queryclass<Numeric::Material::Zhang>(
                        pfx + ".zhang", value.phases[i].zhang);
                    if (plaw == "grimsel")
                        value.phases[i].grimsel_law = true;
                    else if (plaw == "zhang")
                        value.phases[i].zhang_law = true;
                    else if (plaw != "constant")
                        Util::Abort(INFO, "microstructure.phase", i,
                                    ".law must be 'constant', 'grimsel' or 'zhang', got: ",
                                    plaw);
                }

                // P2-27 stage 4: THE SIX PROPERTY SCALARS, where they are read.
                // Under `law = grimsel` the law fills kappa / rho / Cp (and the
                // phase's enthalpy table) at every cell and step, so the three
                // scalars are INERT and refused; E / beta / mu are inert too when
                // the criterion also takes its stress from the law
                // (spall.stress_law = grimsel, RockLaw()) -- with the scalar stress
                // law they are the criterion's E beta (T - T0)/(1 - nu) and stay
                // required. A constant-law phase needs all six (the parse sites
                // are one statement per key: lateral_flux_lock's P2-2b (a1) rule).
                {
                    // P2-29: a zhang phase's law fills the THERMAL scalars too; its mechanics stay
                    // the scalars (law_s is grimsel-only: Zhang has no stress law)
                    const bool law_T = value.phases[i].grimsel_law || value.phases[i].zhang_law;
                    const bool law_s = value.phases[i].grimsel_law && value.sp_stress_law_grimsel;
                    const char* lname = value.phases[i].zhang_law ? "zhang" : "grimsel";
                    auto refuse = [&pp, &pfx, lname] (const char* key, const char* with) {
                        if (pp.contains(pfx + key))
                            Util::Abort(INFO, pfx, key, " is inert under ", pfx, ".law = ", lname, with,
                                        " (the law fills it at every cell and step); remove it (P2-27)");
                    };
                    const char* W_S = " with spall.stress_law = grimsel";
                    if (law_T) refuse(".kappa", "");
                    else pp.query_required(pfx + ".kappa",    value.phases[i].kappa);
                    if (law_s) refuse(".beta", W_S);
                    else pp.query_required(pfx + ".beta",     value.phases[i].beta);
                    if (law_s) refuse(".E", W_S);
                    else pp.query_required(pfx + ".E",        value.phases[i].E);
                    if (law_s) refuse(".mu", W_S);
                    else pp.query_required(pfx + ".mu",       value.phases[i].mu);
                    if (law_T) refuse(".rho", "");
                    else pp.query_required(pfx + ".rho",      value.phases[i].rho);
                    if (law_T) refuse(".Cp", "");
                    else pp.query_required(pfx + ".Cp",       value.phases[i].Cp);
                }

                // Step 6b: optional phase-change scalars per mineral phase.
                // Defaults disable phase change so existing Step 6 input
                // files (e.g. heterogeneous_kappa) remain pure linear-H
                // materials. Latent heats are physical specific heats
                // (J/kg), consistent with Material::effective_rhocp.
                pp.query_default(pfx + ".T_ref",     value.phases[i].T_ref,     273.15);
                pp.query_default(pfx + ".T_m",       value.phases[i].T_m,       1.0e30);
                pp.query_default(pfx + ".fus_width", value.phases[i].fus_width, 0.0);
                pp.query_default(pfx + ".L_m",       value.phases[i].L_m,       0.0);
                pp.query_default(pfx + ".T_vap_lo",  value.phases[i].T_vap_lo,  1.0e30);
                pp.query_default(pfx + ".T_vap_hi",  value.phases[i].T_vap_hi,  1.0e30);
                pp.query_default(pfx + ".L_v",       value.phases[i].L_v,       0.0);

                // Optional quartz alpha-beta transformation strain. Defaults
                // are the atmospheric quartz transition cited in main.tex:
                // T0 = 573 C, dV/V = 0.63 %. The linear strain is derived from
                // the volumetric jump unless the input overrides it directly.
                int ab_enabled = 0;
                pp.query_default(pfx + ".alpha_beta.enabled", ab_enabled, 0);
                pp.query_default(pfx + ".alpha_beta.T0",
                                 value.phases[i].alpha_beta_T0, 846.15);
                pp.query_default(pfx + ".alpha_beta.width",
                                 value.phases[i].alpha_beta_width, 10.0);
                pp.query_default(pfx + ".alpha_beta.dV",
                                 value.phases[i].alpha_beta_dV, 0.0063);
                Util::Assert(INFO, TEST(value.phases[i].alpha_beta_width >= 0.0));
                Util::Assert(INFO, TEST(value.phases[i].alpha_beta_dV > -1.0));
                Set::Scalar ab_linear =
                    std::cbrt(1.0 + value.phases[i].alpha_beta_dV) - 1.0;
                pp.query_default(pfx + ".alpha_beta.linear_strain", ab_linear, ab_linear);
                value.phases[i].alpha_beta_enabled = (ab_enabled != 0);
                value.phases[i].alpha_beta_linear_strain =
                    value.phases[i].alpha_beta_enabled ? ab_linear : 0.0;
            }
            // Normalise fractions so the per-grain sampler is robust to small
            // input rounding (3 phases at 0.30 / 0.50 / 0.20 stays exact).
            Set::Scalar fsum = 0.0;
            for (auto &p : value.phases) fsum += p.fraction;
            Util::Assert(INFO, TEST(fsum > 0.0));
            for (auto &p : value.phases) p.fraction /= fsum;

            value.microstructure_enabled = true;

            // 1 ghost on phase_mf so neighbour comparisons for is_gb need only
            // a Util::RealFillBoundary call; properties are pure cell-local so
            // they are written with 0 ghosts.
            value.RegisterNewFab(value.phase_mf, nullptr, 1, 1, "phase",        true);
            value.RegisterNewFab(value.is_gb_mf, nullptr, 1, 0, "is_gb",        true);
            // Step 6c: explicit grain-topology fields. `grain_id` carries one
            // ghost cell (matches `phase`) so neighbour comparisons can use
            // Util::RealFillBoundary. `is_phase_boundary` is the renamed
            // semantics of the old `is_gb` (kept as an alias);
            // `is_grain_boundary` is the new true-grain neighbour flag.
            value.RegisterNewFab(value.grain_id_mf,           nullptr, 1, 1, "grain_id",          true);
            value.RegisterNewFab(value.is_grain_boundary_mf,  nullptr, 1, 0, "is_grain_boundary", true);
            value.RegisterNewFab(value.is_phase_boundary_mf,  nullptr, 1, 0, "is_phase_boundary", true);
            value.RegisterNewFab(value.kappa_mf, nullptr, 1, 1, "kappa_phase",  true);
            // Step 7: beta_phase, E_phase, mu_phase, T_ref_phase carry 1 ghost
            // so cell-to-node averaging at box boundaries can read neighbour
            // cells via Util::RealFillBoundary.
            value.RegisterNewFab(value.beta_mf,    nullptr, 1, 1, "beta_phase",  true);
            value.RegisterNewFab(value.E_phase_mf, nullptr, 1, 1, "E_phase",     true);
            value.RegisterNewFab(value.mu_phase_mf,nullptr, 1, 1, "mu_phase",    true);
            value.RegisterNewFab(value.T_ref_phase_mf,
                                                   nullptr, 1, 1, "T_ref_phase", true);
            // ⚠⚠ P2-2a change A: THE STRESS-FREE REFERENCE.
            //
            // The STORED INITIAL TEMPERATURE, copied once in Initialize and
            // never updated. It replaces `phases[0].T_ref` -- a DECLARED input
            // constant that can DISAGREE with the deck's own initial condition
            // (`unit/gradient_blind` types 293.15 against an 820 K IC, so the
            // model booked 527 K of superheat that no heating produced).
            //
            // ⚠⚠ REGISTERED UNPLOTTED (`false`) and SEPARATE from
            // `T_ref_phase_mf`, for two independent reasons:
            //   * `T_ref_phase_mf` is read by the EIGENSTRAIN (:1298), so
            //     reusing it would couple the criterion into the elastic solve;
            //   * it is registered `true`, i.e. it lands in the 34-component
            //     hashed cell plotfiles, so reusing it re-pins the whole tree.
            //
            // ⚠ ON REGRID it is filled by INTERPOLATING ITSELF -- that is what
            // RegisterNewFab gives, and it is why nothing re-writes it in
            // Regrid. Re-deriving it from the current temperature would reset
            // the reference to the hot state and return the defect as a mesh
            // artefact.
            value.RegisterNewFab(value.T_init_mf,
                                                   nullptr, 1, 1, "T_init", false);
            // Per-phase alpha-beta transformation parameters and the
            // resulting scalar free strain. Parameter fields are hidden
            // diagnostics; thermal_free_strain is plotted for verification.
            value.RegisterNewFab(value.alpha_beta_linear_strain_mf,
                                                   nullptr, 1, 1, "alpha_beta_linear_strain", false);
            value.RegisterNewFab(value.alpha_beta_T0_mf,
                                                   nullptr, 1, 1, "alpha_beta_T0", false);
            value.RegisterNewFab(value.alpha_beta_width_mf,
                                                   nullptr, 1, 1, "alpha_beta_width", false);
            value.RegisterNewFab(value.thermal_free_strain_mf,
                                                   nullptr, 1, 1, "thermal_free_strain", true);
            value.RegisterNewFab(value.rho_mf,   nullptr, 1, 0, "rho_phase",    true);
            value.RegisterNewFab(value.Cp_mf,    nullptr, 1, 0, "Cp_phase",     true);

            // Step 6c: the grain-boundary contact conductance (MICRO; 0 = off).
            pp.query_default("conductivity.h_gb0",
                             value.h_gb_conductance,   0.0);

            // P2-28: and under the face list's evaluate-only mode (list_state)
            if (value.list_state())
            {
                // phi is a level-set state: phi = z_surface(x,y) - z_cell.
                // Material cells have phi >= 0; void cells (already removed)
                // have phi < 0. Standard AMReX FillPatch preserves it across
                // regrids. phi and removed are shared by spall and vapor
                // removal modes; per-mode diagnostics are registered only
                // when their owning mode is enabled.
                value.RegisterNewFab(value.phi_mf,             nullptr, 1, 0, "phi",             true);
                // 1 ghost: the void-aware conduction face stencil reads
                // `removed` at neighbour cells (i+-1, ...) to zero the flux
                // across void faces. With 0 ghosts those neighbour reads hit
                // unallocated memory at box seams -> decomposition-dependent
                // temperature once cells void (catastrophic under the no-reset
                // void behaviour, where void cells retain their hot enthalpy).
                value.RegisterNewFab(value.removed_mf,         nullptr, 1, 1, "removed",         true);
                value.RegisterNewFab(value.regime_mf,          nullptr, 1, 0, "regime",          true);
                value.RegisterNewFab(value.RoP_mf,             nullptr, 1, 0, "RoP",             true);
            }
            if (value.spall_enabled)
            {
                value.RegisterNewFab(value.spall_event_mf,     nullptr, 1, 0, "spall_event",     true);
                value.RegisterNewFab(value.spall_thickness_mf, nullptr, 1, 0, "spall_thickness", true);
                value.RegisterNewFab(value.RoP_spall_mf,       nullptr, 1, 0, "RoP_spall",       true);
            }
            if (value.vapor_enabled)
            {
                value.RegisterNewFab(value.vapor_event_mf,     nullptr, 1, 0, "vapor_event",     true);
                value.RegisterNewFab(value.vapor_thickness_mf, nullptr, 1, 0, "vapor_thickness", true);
                value.RegisterNewFab(value.RoP_vap_mf,         nullptr, 1, 0, "RoP_vap",         true);
            }
            // Step 17: cell-resolved regime field under sp_weibull. The
            // scalar `regime_mf` remains the per-cell single-regime
            // indicator (registered before regime_field; kept); `regime_field_mf` is
            // the canonical per-cell diagnostic for the sp_weibull
            // per-cell v_n logic — distinct so plotfile consumers can
            // tell the two diagnostics apart.
            if (value.spallation_model == SpallationModel::SpWeibull &&
                (value.spall_enabled || value.vapor_enabled))
            {
                value.RegisterNewFab(value.regime_field_mf, nullptr, 1, 0,
                                     "regime_field", true);
            }
            // Step 15b: the per-cell Weibull flaw length (the face list's R2 reads it at each line's
            // surface cell). P2-28: site A's fields (Sp_field, Sp_cluster_id, h_spall_field,
            // crack_event_count, T_surf_field, ell_eff_field) and the recon_* columns are deleted with it.
            if (value.spallation_model == SpallationModel::SpWeibull && value.weibull_enabled)
                value.RegisterNewFab(value.flaw_a_mf, nullptr, 1, 0, "flaw_a", true);
        }

        // P2-27 (Oct census B14): every cross-key ABORT of this block now lives
        // in CheckInterlocks, called as the LAST statement of Parse. What stays
        // here is what is not an abort: the keys this block parses and the
        // registrations, in their order (CODE_CENSUS F1: registration order is
        // load-bearing for the 2639 gate).

        // S1: surface-resolution study switches (default-off; microstructure
        // path, single level).
        {
            int le = 0;
            pp.query_default("energy_ledger.enabled", le, 0);
            value.surface_follow_mask = value.surface_patch.face_list;
            value.ledger_enabled      = (le != 0);
            // P0b deliverable 1: robin_form = surface. Same prerequisites as
            // pinned -- it is the same closure with the surface temperature
            // computed rather than remembered -- PLUS the reconstruction
            // itself, without which the form has no surface to sit on.
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
            {
                value.RegisterIntegratedVariable(&value.surf_side_share, "surf_side_share", false);
            }
            if (value.ledger_enabled)
            {
                value.RegisterIntegratedVariable(&value.ledger_E_H,     "ledger_E_H",     false);
                value.RegisterIntegratedVariable(&value.ledger_E_src,   "ledger_E_src",   false);
                // P2-28 (P2-27 DR-1): RENAMED -- the all-cell residual omits the flake credit by
                // design, so it is NOT a closure on a spalling deck; ledger_err_solid closes
                value.RegisterIntegratedVariable(&value.ledger_err_raw, "ledger_err_raw", false);
                value.RegisterIntegratedVariable(&value.ledger_P_beam,  "ledger_P_beam",  false);
                value.RegisterIntegratedVariable(&value.ledger_P_robin, "ledger_P_robin", false);
                value.RegisterIntegratedVariable(&value.ledger_P_loss,  "ledger_P_loss",  false);
                value.RegisterIntegratedVariable(&value.ledger_surface_missing_cols,
                                                 "ledger_surface_missing_cols", false);
                // P0b: appended last, and ALSO gated on the reconstruction --
                // a tree input with the ledger on but no P0b key must not gain
                // a thermo.dat column, or the 2639-file identity breaks on the
                // header alone whatever the numbers do. P2-28: the reconstruction's one home is the face list.
                if (value.surface_patch.face_list)
                {
                    value.RegisterIntegratedVariable(&value.ledger_E_flake,
                                                     "ledger_E_flake", false);
                    value.RegisterIntegratedVariable(&value.ledger_E_H_solid,
                                                     "ledger_E_H_solid", false);
                    value.RegisterIntegratedVariable(&value.ledger_err_solid,
                                                     "ledger_err_solid", false);
                }
            }
            // D2a: the spatial-flame diagnostic, after every pre-D2a variable (P2-28: the pinned
            // branch counters and patch_P_robin -- the column patch diagnostics -- are deleted).
            if (value.surface_patch.enabled && value.surface_patch.use_expr)
            {
                value.RegisterIntegratedVariable(&value.patch_min_standoff, "patch_min_standoff", false);
            }
            // D2a2: wall-jet enthalpy closure diagnostics, after every
            // pre-D2a2 variable.
            if (value.surface_patch.enabled && value.surface_patch.jet_enthalpy)
            {
                value.RegisterIntegratedVariable(&value.jet_P_face,    "jet_P_face",    false);
                value.RegisterIntegratedVariable(&value.jet_P_cap,     "jet_P_cap",     false);
                value.RegisterIntegratedVariable(&value.jet_T_stag,    "jet_T_stag",    false);
                value.RegisterIntegratedVariable(&value.jet_s_c,       "jet_s_c",       false);
                value.RegisterIntegratedVariable(&value.jet_P_decay,   "jet_P_decay",   false);
                value.RegisterIntegratedVariable(&value.jet_T_exhaust, "jet_T_exhaust", false);
                value.RegisterIntegratedVariable(&value.jet_P_exhaust, "jet_P_exhaust", false);
                value.RegisterIntegratedVariable(&value.jet_r_reach,   "jet_r_reach",   false);
            }
            // D2b: feet descent and exhaust recirculation, after every pre-D2b
            // variable.
            if (value.surface_patch.enabled && value.surface_patch.feet)
            {
                value.RegisterIntegratedVariable(&value.nozzle_z,           "nozzle_z",           false);
                value.RegisterIntegratedVariable(&value.foot_z,             "foot_z",             false);
                value.RegisterIntegratedVariable(&value.foot_cols,          "foot_cols",          false);
                value.RegisterIntegratedVariable(&value.foot_carry_cols,    "foot_carry_cols",    false);
                value.RegisterIntegratedVariable(&value.foot_stall_time,    "foot_stall_time",    false);
            }
            if (value.surface_patch.enabled && value.surface_patch.jet_enthalpy &&
                value.surface_patch.jet_T_ent_exhaust)
            {
                value.RegisterIntegratedVariable(&value.jet_T_rec,    "jet_T_rec",    false);
                value.RegisterIntegratedVariable(&value.jet_P_recirc, "jet_P_recirc", false);
            }
            // D2c: each block only when its key is non-default, after every
            // pre-D2c variable.
            if (value.surface_patch.enabled && value.surface_patch.jet_enthalpy)
            {
                const auto& sp = value.surface_patch;
                if (sp.jet_De_momentum)
                    value.RegisterIntegratedVariable(&value.jet_D_e, "jet_D_e", false);
                if (sp.jet_free_surface)
                {
                    value.RegisterIntegratedVariable(&value.jet_fs_cols, "jet_fs_cols", false);
                    value.RegisterIntegratedVariable(&value.jet_m_ratio, "jet_m_ratio", false);
                    value.RegisterIntegratedVariable(&value.jet_T_mouth, "jet_T_mouth", false);
                }
                if (sp.jet_neg_count)
                    value.RegisterIntegratedVariable(&value.jet_P_neg, "jet_P_neg", false);
            }
            if (value.surface_patch.enabled && value.surface_patch.feet &&
                value.surface_patch.foot_body_clearance)
            {
                value.RegisterIntegratedVariable(&value.foot_mech_vol,  "foot_mech_vol",  false);
                value.RegisterIntegratedVariable(&value.foot_mech_cols, "foot_mech_cols", false);
                // P2-2f: appended after foot_mech_cols and gated on the SAME
                // key, so no input without foot_body_clearance = 1 gains a
                // thermo.dat column and the 2639-file identity is untouched.
                value.RegisterIntegratedVariable(&value.foot_mech_E,    "foot_mech_E",    false);
            }
            // P2-2 changes 1 and 4: the witnesses (b5) and (b3) read.
            // ⚠ Registered ONLY when a phase actually uses the grimsel law or
            // the relief factor is set, so a default leg gains NO column and
            // the 2639-file identity gate is untouched. This is the "keyed off
            // by default" half of the rule, NOT a re-pin request.
            {
                bool any_law = false;
                for (const auto& ph : value.phases) any_law = any_law || ph.grimsel_law;
                if (any_law)
                {
                    value.RegisterIntegratedVariable(&value.prop_refresh_n,
                                                     "prop_refresh_n", false);
                    value.RegisterIntegratedVariable(&value.prop_kappa_mean,
                                                     "prop_kappa_mean", false);
                    value.RegisterIntegratedVariable(&value.prop_beta_mean,
                                                     "prop_beta_mean", false);
                    value.RegisterIntegratedVariable(&value.prop_E_mean,
                                                     "prop_E_mean", false);
                    value.RegisterIntegratedVariable(&value.prop_T_mean,
                                                     "prop_T_mean", false);
                }
                if (value.sp_relief < 1.0)
                    value.RegisterIntegratedVariable(&value.sp_relief,
                                                     "sp_relief", false);
            }

        }

        // P2-15: the buckling condition's witnesses, registered ONLY under the key
        // (its requirements -- the face list, the reconstruction -- are checked in
        // CheckInterlocks).
        if (value.sp_buckle_span > 0.0)
        {
            value.RegisterIntegratedVariable(&value.bk_span_out,      "bk_span",          false);
            value.RegisterIntegratedVariable(&value.bk_n_eval,        "bk_n_eval",        false);
            value.RegisterIntegratedVariable(&value.bk_n_block,       "bk_n_block",       false);
            value.RegisterIntegratedVariable(&value.bk_n_inf,         "bk_n_inf",         false);
            value.RegisterIntegratedVariable(&value.bk_n_block_floor, "bk_n_block_floor", false);
            value.RegisterIntegratedVariable(&value.bk_n_block_wall,  "bk_n_block_wall",  false);
            value.RegisterIntegratedVariable(&value.bk_n_block_roof,  "bk_n_block_roof",  false);
            value.RegisterIntegratedVariable(&value.bk_max_pass,      "bk_max_pass",      false);
        }

        // P2-8: THE FACE LIST. Its thermo columns and plotfile fields are
        // registered LAST and ONLY under the key, so a key-off run gains no
        // column and no field (CODE_CENSUS F1: registration order is
        // load-bearing for the 2639 gate). Every feature it does not carry is
        // refused in CheckInterlocks.
        if (value.surface_patch.face_list)
        {
            const SurfacePatch& sp = value.surface_patch;
            value.RegisterIntegratedVariable(&value.fl_n_top,   "fl_n_top",   false);
            value.RegisterIntegratedVariable(&value.fl_n_up,    "fl_n_up",    false);
            value.RegisterIntegratedVariable(&value.fl_n_side,  "fl_n_side",  false);
            value.RegisterIntegratedVariable(&value.fl_n_down,  "fl_n_down",  false);
            value.RegisterIntegratedVariable(&value.fl_P_top,   "fl_P_top",   false);
            value.RegisterIntegratedVariable(&value.fl_P_up,    "fl_P_up",    false);
            value.RegisterIntegratedVariable(&value.fl_P_side,  "fl_P_side",  false);
            value.RegisterIntegratedVariable(&value.fl_P_down,  "fl_P_down",  false);
            value.RegisterIntegratedVariable(&value.fl_P_march, "fl_P_march", false);
            value.RegisterIntegratedVariable(&value.fl_w_fallback, "fl_w_fallback", false);
            value.RegisterIntegratedVariable(&value.fl_h_top_mismatch,
                                             "fl_h_top_mismatch", false);
            value.RegisterIntegratedVariable(&value.fl_id_march,  "fl_id_march",  false);
            value.RegisterIntegratedVariable(&value.fl_id_ledger, "fl_id_ledger", false);
            value.RegisterNewFab(value.fl_w_mf, nullptr, 1, 0, "fl_w", true);
            value.RegisterNewFab(value.fl_h_mf, nullptr, 1, 0, "fl_h", true);

            // P2-9 deliverable 5: ONE breaking rule under the same key.
            value.RegisterIntegratedVariable(&value.fr_n_c[0],    "fr_n_floor",    false);
            value.RegisterIntegratedVariable(&value.fr_n_c[1],    "fr_n_wall",     false);
            value.RegisterIntegratedVariable(&value.fr_n_c[2],    "fr_n_roof",     false);
            value.RegisterIntegratedVariable(&value.fr_v_c[0],    "fr_v_floor",    false);
            value.RegisterIntegratedVariable(&value.fr_v_c[1],    "fr_v_wall",     false);
            value.RegisterIntegratedVariable(&value.fr_v_c[2],    "fr_v_roof",     false);
            value.RegisterIntegratedVariable(&value.fr_fb_c[0],   "fr_fb_floor",   false);
            value.RegisterIntegratedVariable(&value.fr_fb_c[1],   "fr_fb_wall",    false);
            value.RegisterIntegratedVariable(&value.fr_fb_c[2],   "fr_fb_roof",    false);
            value.RegisterIntegratedVariable(&value.fr_smax_c[0], "fr_smax_floor", false);
            value.RegisterIntegratedVariable(&value.fr_smax_c[1], "fr_smax_wall",  false);
            value.RegisterIntegratedVariable(&value.fr_smax_c[2], "fr_smax_roof",  false);
            value.RegisterIntegratedVariable(&value.fr_clip,      "fr_clip",       false);
            value.RegisterIntegratedVariable(&value.fr_cap,       "fr_cap",        false);
            value.RegisterIntegratedVariable(&value.fr_sheet,     "fr_sheet",      false);
            value.RegisterIntegratedVariable(&value.fr_multi,     "fr_multi",      false);
            value.RegisterIntegratedVariable(&value.fr_overlap,   "fr_overlap",    false);
            value.RegisterIntegratedVariable(&value.fr_id_remap,  "fr_id_remap",   false);
            value.RegisterIntegratedVariable(&value.fr_ncell,     "fr_ncell",      false);
            value.RegisterIntegratedVariable(&value.fr_through,   "fr_through",    false);
            value.RegisterIntegratedVariable(&value.fr_E_drop,    "fr_E_drop",     false);
            value.RegisterIntegratedVariable(&value.fr_fresh,     "fr_fresh",      false);
            // P2-9 rev 2 item 7
            value.RegisterIntegratedVariable(&value.fr_misaligned,  "fr_misaligned",  false);
            value.RegisterIntegratedVariable(&value.fr_nodepth_c[0], "fr_nodepth_floor", false);
            value.RegisterIntegratedVariable(&value.fr_nodepth_c[1], "fr_nodepth_wall",  false);
            value.RegisterIntegratedVariable(&value.fr_nodepth_c[2], "fr_nodepth_roof",  false);
            value.RegisterIntegratedVariable(&value.fr_t_removal,   "fr_t_removal",   false);
            value.RegisterIntegratedVariable(&value.fr_t_step,      "fr_t_step",      false);
            // P2-24: lines NOT scanned because their surface cell is fully liquid, or a fully liquid cell lies
            // inside their window (a solid lid over a melt); per step, appended LAST on every face-list deck
            value.RegisterIntegratedVariable(&value.fr_n_molten,    "fr_n_molten",    false);
            value.RegisterIntegratedVariable(&value.fr_n_lid,       "fr_n_lid",       false);
            if (value.vapor_enabled)
            {
                // P2-24: the vapour rule's counts (per step) and the enthalpy that left with it (cumulative, J)
                value.RegisterIntegratedVariable(&value.fr_n_vapor,     "fr_n_vapor",     false);
                value.RegisterIntegratedVariable(&value.fr_n_vapor_int, "fr_n_vapor_int", false);
                value.RegisterIntegratedVariable(&value.ledger_E_vapor, "ledger_E_vapor", false);
            }
            // 3 ghosts: the smoothed normal reads +-3 cells (5x5x5 kernel, then
            // a central difference)
            value.RegisterNewFab(value.fr_f_mf,     nullptr, 1, 3, "fr_f",     true);
            value.RegisterNewFab(value.fr_score_mf, nullptr, 1, 0, "fr_score", true);
            value.RegisterNewFab(value.fr_mu_mf,    nullptr, 1, 0, "fr_mu",    true);

            // P2-11b: face_heat = path -- its outputs are registered LAST and
            // only under path (a position run gains none).
            if (sp.face_heat_path)
            {
                value.RegisterIntegratedVariable(&value.fl_H_stag,      "fl_H_stag",      false);
                value.RegisterIntegratedVariable(&value.fl_n_src,       "fl_n_src",       false);
                value.RegisterIntegratedVariable(&value.fl_n_blocked,   "fl_n_blocked",   false);
                value.RegisterIntegratedVariable(&value.fl_n_unreached, "fl_n_unreached", false);
                value.RegisterIntegratedVariable(&value.fl_s_max,       "fl_s_max",       false);
                value.RegisterIntegratedVariable(&value.fl_t_path,      "fl_t_path",      false);
                value.RegisterNewFab(value.fl_s_mf,  nullptr, 1, 0, "fl_s",  true);
                value.RegisterNewFab(value.fl_Tg_mf, nullptr, 1, 0, "fl_Tg", true);
            }
        }

        if (value.vapor_enabled && !value.vapor_A_beam_override)
        {
            // Default beam-area convention: Gaussian effective area
            //   A_beam = 0.5 * pi * omega0^2
            // (matches the (1-R)*P0/A_eff intensity scaling used by the
            // Gaussian source, where peak intensity is 2P/(pi*omega0^2)).
            constexpr Set::Scalar pi = 3.141592653589793238462643383279502884;
            value.vapor_A_beam = 0.5 * pi * value.beam.omega0 * value.beam.omega0;
            if (!(std::isfinite(value.vapor_A_beam) && value.vapor_A_beam > 0.0))
                Util::Abort(INFO, "vapor: default A_beam is non-positive (",
                            value.vapor_A_beam, "); set vapor.A_beam explicitly");
        }
        if (value.vapor_enabled && !value.vapor_T0_override)
            value.vapor_T0 = value.T_amb;

        // P2-27: LAST, so every member it reads is parsed
        CheckInterlocks(value);
    }

    // ⚠⚠ P2-27 (Oct census B14): THE CROSS-KEY INTERLOCKS, IN ONE PLACE. Called
    // ONCE, as the last statement of Parse, so every member a condition reads has
    // already been parsed. Each condition and message is VERBATIM from where it
    // used to sit, and they run in the order they used to fire: the sub-parse
    // checks (ParseSpallSettings, ParseBitSettings, in Parse's call order; the
    // two ParseSpallationModel ones died with the damage lineage, P2-27 stage 2),
    // then Parse's own end-of-parse block, then the face
    // list's cubic-cell requirement (a per-call runtime check in BuildFaceList
    // until P2-27; geom is known at Parse). ONLY ABORTS live here: every key
    // parse and every registration stays in Parse, in its order (registration
    // order is load-bearing, CODE_CENSUS F1). SurfacePatch::Parse keeps its own
    // pp.contains-based checks.
    // ⚠ A deck that breaks TWO rules can now meet a different one first: a
    // sub-parse that used to abort early now runs on to its later key queries
    // (bit.enabled with spall and vapour off still reads bit.feed_rate, which
    // is required). One broken rule fires exactly as before.
    static void CheckInterlocks(MMWSpalling& value)
    {
        // ---- from ParseSpallSettings --------------------------------------
        if (value.sp_flake_h_max > 0.0 &&
            value.sp_flake_h_max >= value.spall_h_spall_max)
            Util::Abort(INFO, "spall.flake_h_max (", value.sp_flake_h_max,
                        ") must be SMALLER than spall.h_spall_max (",
                        value.spall_h_spall_max, "): the first is a physical "
                        "flake thickness, the second is the depth the scan "
                        "searches. A bound at or above the search depth is "
                        "inert and hides that the search depth is doing the "
                        "bounding -- which is the fault this key exists to fix.");

        // ---- Parse's end-of-parse block -----------------------------------
        // P2-27 (census C1a): with the damage-law removal deleted, the only
        // removal criterion is sp_weibull's. Under the default `none` a spall
        // or vapour key would otherwise run and remove nothing, silently.
        if ((value.spall_enabled || value.vapor_enabled) &&
            value.spallation_model != SpallationModel::SpWeibull)
            Util::Abort(INFO, "spall.enabled / vapor.enabled need spallation.model = sp_weibull: "
                              "the damage-law removal was retired in P2-27, so under the default "
                              "spallation.model = none nothing would fire");
        // P2-28 (census H): ONE REMOVAL PATH. The column path is deleted, so removal runs on the face
        // list or not at all (the list's evaluate-only mode, spall and vapour off, is the onset decks').
        if ((value.spall_enabled || value.vapor_enabled) && !value.surface_patch.face_list)
            Util::Abort(INFO, "spall.enabled / vapor.enabled need surface_patch.face_list = 1: the column "
                              "removal path was retired in P2-28 (one removal path, the face list)");
        // P2-28: the jet march bills the face list's faces; its column flux is deleted
        if (value.surface_patch.enabled && value.surface_patch.jet_enthalpy &&
            !value.surface_patch.face_list)
            Util::Abort(INFO, "surface_patch.jet_closure = enthalpy needs surface_patch.face_list = 1: the "
                              "march's column flux was retired in P2-28");

        if (value.spall_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spall.enabled requires microstructure.nphases");

        // ⚠⚠ P2-2d d7: CROSS-KEY, AND IT HAS TO BE HERE. `grimsel_table` reads
        // K_Ic0 and the at-temperature ratio off phase 0's law object, so it
        // needs a phase to read. `ParseSpallationModel` runs BEFORE
        // `phases` is resized, so the check belongs with the other
        // cross-key aborts, not beside the key's own parse.
        if (value.kic_law == KIcLaw::GrimselTable && value.phases.empty())
            Util::Abort(INFO, "spallation.sp.kic.law = grimsel_table needs "
                              "phase 0 to carry the Grimsel law: K_Ic0 and the "
                              "at-temperature ratio come from it, never from a "
                              "new scalar. Set microstructure.nphases >= 1.");

        if (value.vapor_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: vapor.enabled requires microstructure.nphases");

        if (value.spallation_model == SpallationModel::SpWeibull &&
            !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spallation.model = sp_weibull requires "
                              "microstructure.nphases (Step 15 Weibull samples will key "
                              "on grain-boundary faces)");

        if (value.ledger_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "energy_ledger.enabled requires "
                              "the microstructure path (microstructure.nphases)");
        if (value.ledger_enabled && value.max_level > 0)
            Util::Abort(INFO, "energy_ledger.enabled is single-level only "
                              "(amr.max_level = 0)");
        if (value.surface_patch.enabled &&
            value.surface_patch.robin_form == SurfacePatch::RobinForm::Face &&
            !value.microstructure_enabled)
            Util::Abort(INFO, "surface_patch.robin_form = face is implemented on "
                              "the microstructure path only");
        // A1: the face-form Newton (SurfaceCellFlux) relies on a residual
        // that is monotone increasing in T_f. With epsilon_high < epsilon
        // the radiation term steps DOWN at T_eps_break, so the face
        // balance can have more than one root.
        if (value.surface_patch.enabled &&
            value.surface_patch.robin_form == SurfacePatch::RobinForm::Face &&
            value.eps_rad_hi < value.eps_rad)
            Util::Abort(INFO, "surface_patch.robin_form = face requires "
                              "losses.epsilon_high >= losses.epsilon (got ", value.eps_rad_hi,
                              " < ", value.eps_rad, "): the face-temperature Newton assumes a "
                              "monotone residual with a unique root");
        // P0b deliverable 1: robin_form = surface, the same prerequisites as
        // pinned PLUS the reconstruction.
        if (value.surface_patch.enabled &&
            value.surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
        {
            if (!value.microstructure_enabled)
                Util::Abort(INFO, "surface_patch.robin_form = surface is implemented "
                                  "on the microstructure path only");
            // P2-28: the form's surface temperature is the face list's (R1): the column
            // reconstruction it used without the list is deleted
            if (!value.surface_patch.face_list)
                Util::Abort(INFO, "surface_patch.robin_form = surface requires surface_patch.face_list = 1: "
                                  "the column reconstruction it read without the list was retired in P2-28");
            // P2-28: or the face list's evaluate-only mode (the criterion runs, nothing is removed)
            if (!((value.spall_enabled || value.evaluate_only()) &&
                  value.spallation_model == SpallationModel::SpWeibull))
                Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                  "spall.enabled = 1 and spallation.model = sp_weibull");
            if (value.eps_rad_hi < value.eps_rad)
                Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                  "losses.epsilon_high >= losses.epsilon (got ",
                                  value.eps_rad_hi, " < ", value.eps_rad, ")");
        }
        // A1: removal-event log validity.
        if (value.sp_emit_removal_events_csv)
        {
            // P2-28: evaluate-only writes the rows production would (nothing removed)
            if (!((value.spall_enabled || value.evaluate_only()) &&
                  value.spallation_model == SpallationModel::SpWeibull))
                Util::Abort(INFO, "spall.removal_events_csv = 1 requires spall.enabled = 1 "
                                  "and spallation.model = sp_weibull");
        }
        // P2-28: the flame map is evaluated per FACE (F2): the face list
        if (value.surface_patch.enabled && value.surface_patch.use_expr && !value.surface_patch.face_list)
            Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr need surface_patch.face_list = 1 (F2 "
                              "evaluates the map at every face centre; the column path was retired in P2-28)");
        if (value.surface_patch.enabled && value.surface_patch.feet &&
            value.surface_patch.foot_body_clearance)
        {
            if (!value.removal_enabled())
                Util::Abort(INFO, "surface_patch.foot_body_clearance = 1 requires removal "
                                  "(spall.enabled or vapor.enabled)");
        }

        // P2-15: the buckling condition lives on the face list only (the column path never
        // reads the key -- retire, don't stack) and reads the reconstructed skin length, so
        // either missing would leave the key SILENTLY INERT: both abort, naming the key.
        // P2-23: the record of B's evaluations needs B, or it would be silently empty
        if (value.sp_buckle_log_csv && !(value.sp_buckle_span > 0.0))
            Util::Abort(INFO, "spall.buckle_log_csv = 1 (P2-23) logs the buckling condition's "
                        "evaluations: it needs spall.buckle_span > 0");
        if (value.sp_buckle_span > 0.0)
        {
            if (!value.surface_patch.face_list)
                Util::Abort(INFO, "spall.buckle_span > 0 (P2-15, the buckling thickness "
                            "condition) is evaluated only on the face list: it needs "
                            "surface_patch.face_list = 1; got spall.buckle_span = ",
                            value.sp_buckle_span, " m");
        }

        // P2-8: THE FACE LIST, validated where every key it depends on has been parsed.
        if (value.surface_patch.face_list)
        {
            const SurfacePatch& sp = value.surface_patch;
            // P2-11b: face_heat = path. Every feature it does not carry ABORTS,
            // naming the key and the step that owns it -- FIRST in this block,
            // so a robin_form the list itself refuses still names its key here,
            // and in an order where each abort is reachable on its own (the
            // free surface needs feet).
            if (sp.face_heat_path)
            {
                auto ph_abort = [] (const char* key) {
                    Util::Abort(INFO, "surface_patch.face_heat = path (P2-11b, heat along "
                                "the surface path) does not carry ", key, "; that is P2-12's");
                };
                if (sp.jet_free_surface) ph_abort("surface_patch.jet_free_surface = 1");
                if (sp.feet) ph_abort("surface_patch.nozzle_descent = feet");
                if (sp.robin_form != SurfacePatch::RobinForm::Surface)
                    ph_abort("surface_patch.robin_form other than surface");
                if (!sp.jet_enthalpy) ph_abort("surface_patch.jet_closure other than enthalpy");
            }
            // P2-23: with the patch OFF there is no flame to map -- every listed face
            // is billed at h = 0 (its losses only); the map is required only when on.
            // P2-28: a PRESCRIBED patch is carried: its TOP faces take the half-cell
            // Dirichlet term at T_expr(t) (FaceListHeat), every other face h = 0.
            if (sp.enabled && sp.mode == SurfacePatch::Mode::ConvectiveFlame && !sp.use_expr)
                Util::Abort(INFO, "surface_patch.face_list = 1 requires mode = "
                            "convective_flame with an h_expr map: F2 evaluates h_expr "
                            "at every face centre");
            if (sp.robin_form != SurfacePatch::RobinForm::Surface)
                Util::Abort(INFO, "surface_patch.face_list = 1 requires robin_form = "
                            "surface: F4 keeps the column-top face on the surface form");
            if (!value.microstructure_enabled)
                Util::Abort(INFO, "surface_patch.face_list = 1 requires the "
                            "microstructure path (its consumers live in "
                            "AdvanceMicrostructure)");

            // P2-9 deliverable 5: ONE breaking rule under the same key. Every
            // path it does not carry ABORTS, naming the step that owns it.
            auto fr_abort = [] (const char* what, const char* owner) {
                Util::Abort(INFO, "surface_patch.face_list = 1 (P2-9, one breaking "
                            "rule along one line per surface cell) does not carry ",
                            what, "; that is ", owner, "'s");
            };
            // P2-24: regime 2 IS carried now -- the face list's vapour rule (a LISTED cell whose enthalpy has
            // passed its phase table's H_vap_hi leaves, once; FaceRemoval.H). It reads phase 0's table: one rock
            // phase (SR 13), with a declared vaporisation band.
            if (value.vapor_enabled)
            {
                if (value.phases.size() != 1)
                    Util::Abort(INFO, "vapor.enabled under the face list (P2-24) reads ONE rock phase's enthalpy "
                                "table; got microstructure.nphases = ", value.phases.size());
                for (const auto& ph : value.phases)
                    if (!(ph.T_vap_lo < 1.0e29 && ph.L_v > 0.0))
                        Util::Abort(INFO, "vapor.enabled under the face list (P2-24) needs the rock phase's "
                                    "vaporisation band: microstructure.phase0.T_vap_lo < 1e30 and L_v > 0; got "
                                    "T_vap_lo ", ph.T_vap_lo, ", L_v ", ph.L_v);
            }
            // P2-27 stage 5 (census G): regime 4 IS carried now -- the pad clearance
            // runs on the face list (FrFootClearance, FaceRemoval.H)
            if (!value.sp_local_thermoelastic)
                fr_abort("a stress source that is not temperature-derived "
                         "(spall.stress_source != local_thermoelastic)", "the MMW step");
            if (value.max_level > 0) fr_abort("amr.max_level > 0", "a later step");
            // ⚠ beyond the packet's list (declared): the line rule has no
            // route for these, so each would be silently dropped
            // P2-28: spall AND vapour off is the evaluate-only mode (the onset decks); vapour without spall stays
            // refused (the vapour rule rides on the spalling list)
            if (value.spallation_model != SpallationModel::SpWeibull ||
                (!value.spall_enabled && value.vapor_enabled))
                fr_abort("a removal model other than sp_weibull spallation", "the MMW step");

            // F1: cubic cells (the face list's areas and lines assume dx = dy = dz).
            // Until P2-27 this was checked in BuildFaceList on every call; the
            // face list is single-level (above), so level 0's cell size decides.
            const Set::Scalar* DX = value.geom[0].CellSize();
            if (!(std::abs(DX[0] - DX[1]) <= 1.0e-12 * DX[0] &&
                  std::abs(DX[0] - DX[2]) <= 1.0e-12 * DX[0]))
                Util::Abort(INFO, "surface_patch.face_list = 1 requires cubic cells; got ",
                            DX[0], ", ", DX[1], ", ", DX[2]);
        }
    }
