    // Step 6d / Step 15a: surface patch driver. A circular patch on the +z
    // face applies one of two heat-load modes inside the patch footprint:
    //   mode = prescribed_T (default, Step 6d) — holds T = T_expr(t) on the top
    //     FACE and injects the half-cell Dirichlet flux
    //     q = 2*kappa*(T_patch - T_cell)/dz as a divergence contribution
    //     2*kappa*(T_patch - T_cell)/dz^2 on the in-patch top row.
    //     Under the face list the top faces carry it.
    //   mode = convective_flame (Step 15a) — applies a Robin flame-jet BC
    //     q = h_fl*(T_flame - T_surface), injected (robin_form = cell, the
    //     default) as the divergence
    //     contribution h_fl*(T_flame - T_cell)/dz on the in-patch top row.
    // Outside the patch the existing zhi BC continues to drive the face flux
    // (typically Neumann/adiabatic). Disabled by default; the default mode
    // is prescribed_T.
    struct SurfacePatch
    {
        static constexpr const char *name = "surface_patch";

        enum class Mode { PrescribedTemperature, ConvectiveFlame };
        enum class RobinForm { Cell, Face, Surface };

        bool enabled = false;
        Mode mode = Mode::PrescribedTemperature;
        RobinForm robin_form = RobinForm::Cell;
        // P2-8: THE SURFACE FACE LIST (`surface_patch.face_list`, default 0).
        // With it on, every solid face touching void (all six directions) is
        // listed once per step and heated, weighted and judged from that one
        // list (BuildFaceList / FaceListHeat / FaceWeight). 0 => no list: the
        // patch acts on the column-top cells only, and removal, the h_expr
        // map, jet_closure = enthalpy and robin_form = surface are refused
        // (CheckInterlocks).
        bool face_list = false;
        // P2-11b: surface_patch.face_heat = position (default, today) | path
        // (heat along the surface path; FaceListPath), and the gas barrier
        // surface_patch.barrier_expr(x, y, z, r), read only under path.
        bool face_heat_path = false;
        bool has_barrier = false;
        amrex::Parser            barrier_parser;
        amrex::ParserExecutor<4> barrier_f;
        Set::Scalar x0 = 0.0, y0 = 0.0, radius = 0.0;
        Set::Scalar h_conv  = 0.0;  // convective_flame: flame coefficient h_fl
        Set::Scalar T_flame = 0.0;  // convective_flame: flame temperature
        amrex::Parser            T_parser;
        amrex::ParserExecutor<1> T_f;
        // D2a: spatial flame expressions (x, y, r, s, t); use_expr = h_expr given
        // (with T_flame_expr, or alone under jet_closure = enthalpy).
        bool use_expr = false;
        Set::Scalar nozzle_z0 = 0.0, nozzle_feed = 0.0, nozzle_collision_radius = 0.0;
        amrex::Parser            h_parser, Tg_parser;
        amrex::ParserExecutor<5> h_f, Tg_f;
        // D2a2: jet_closure = enthalpy (T_flame_expr not read, Tg_f uncompiled).
        bool jet_enthalpy = false;
        Set::Scalar jet_mdot = 0.0, jet_T_nozzle = 0.0, jet_cp = 0.0;
        Set::Scalar jet_T_ent = 293.15, jet_dr = 0.0, jet_profile_interval = 0.0;
        Set::Scalar side_power_interval = 0.0;   // rev7 d2, 0 = off
        Set::Scalar jet_D = 0.0, jet_core_length = 5.0;
        bool jet_stag_decay = true, jet_entrained_mass = false;
        bool jet_T_ent_exhaust = false;     // D2b: jet_T_ent_mode = exhaust:
                                            // recirculation on
        // D2c: march options (see JetEnthalpyMarch); all default-off.
        bool jet_De_momentum = false;       // jet_decay_diameter = momentum
        Set::Scalar jet_De_ref = 0.0;       // full-jet D_e at jet_T_ent [m]
        bool jet_free_surface = false;
        Set::Scalar jet_fs_aspect = 1.0, jet_fs_exponent = 1.0;
        bool jet_bin_exponential = false;   // jet_bin_update = exponential
        bool jet_neg_count = false;         // jet_negative_flux = count
        // D2b: burner-on-feet descent (nozzle_descent = feet).
        bool feet = false;
        Set::Scalar foot_r_inner = 0.0, foot_r_outer = 0.0, foot_standoff = 0.0;
        int foot_npads = 3;
        Set::Scalar foot_pad_quantile = 0.9;
        Set::Scalar nozzle_feed_max = std::numeric_limits<Set::Scalar>::infinity();
        bool foot_body_clearance = false;   // D2c

        // prescribed_T only: patch surface temperature at time t. Not valid
        // under convective_flame (T_f is left uncompiled in that mode).
        Set::Scalar T_at(Set::Scalar t) const { return T_f(t); }

        static void Parse(SurfacePatch &v, IO::ParmParse &pp)
        {
            int en = 0;
            pp.query_default("enabled", en, 0);
            v.enabled = (en != 0);

            // P2-23: THE FACE LIST WITHOUT A FLAME. The keys the face list itself
            // needs -- the switch, the surface form, and the listing mask (centre
            // and radius) -- are read BEFORE the `enabled` return, so a flameless
            // deck that includes the shared model form (face_list = 1, robin_form =
            // surface) reads them instead of aborting on unread keys. With the
            // patch ON the mask (x0, y0, radius) is required; with it OFF the mask
            // defaults to nothing listed (radius 0). Every flame-only key and
            // every robin_form check stays behind the return.
            {
                int fl = 0;
                pp.query_default("face_list", fl, 0);
                if (fl != 0 && fl != 1)
                    Util::Abort(INFO, "surface_patch.face_list must be 0 or 1, got: ", fl);
                v.face_list = (fl == 1);
            }
            std::string robin_str;
            pp.query_default("robin_form", robin_str, std::string("cell"));
            if      (robin_str == "cell")   v.robin_form = RobinForm::Cell;
            else if (robin_str == "face")   v.robin_form = RobinForm::Face;
            else if (robin_str == "surface") v.robin_form = RobinForm::Surface;
            else if (robin_str == "pinned")
                Util::Abort(INFO, "surface_patch.robin_form = pinned was retired in P2-28 (the pin rule "
                                  "and the column path are deleted): use surface (the face list) or "
                                  "cell / face (a deck without removal)");
            else Util::Abort(INFO, "surface_patch.robin_form must be 'cell', 'face' "
                                   "or 'surface', got: ", robin_str);
            if (v.enabled)
            {
                pp.query_required("x0",     v.x0);
                pp.query_required("y0",     v.y0);
                pp.query_required("radius", v.radius);
            }
            else
            {
                pp.query_default("x0",     v.x0,     0.0);
                pp.query_default("y0",     v.y0,     0.0);
                pp.query_default("radius", v.radius, 0.0);
            }
            if (!v.enabled) return;
            // P2-11b: read only under the face list; position = today.
            if (pp.contains("face_heat"))
            {
                if (!v.face_list)
                    Util::Abort(INFO, "surface_patch.face_heat is read only under "
                                      "surface_patch.face_list = 1");
                std::string fh;
                pp.query_required("face_heat", fh);
                if      (fh == "position") v.face_heat_path = false;
                else if (fh == "path")     v.face_heat_path = true;
                else Util::Abort(INFO, "surface_patch.face_heat must be 'position' or "
                                       "'path', got: ", fh);
            }
            if (pp.contains("barrier_expr"))
            {
                if (!v.face_heat_path)
                    Util::Abort(INFO, "surface_patch.barrier_expr is read only under "
                                      "surface_patch.face_heat = path");
                std::string be;
                pp.query_required("barrier_expr", be);
                v.barrier_parser = amrex::Parser(be);
                v.barrier_parser.registerVariables({"x", "y", "z", "r"});
                v.barrier_f = v.barrier_parser.compile<4>();
                v.has_barrier = true;
            }

            // `mode` is the Step 15a authority; default prescribed_T keeps
            // legacy inputs unchanged.
            std::string mode_str;
            pp.query_default("mode", mode_str, std::string("prescribed_T"));
            if (mode_str == "prescribed_T")
                v.mode = Mode::PrescribedTemperature;
            else if (mode_str == "convective_flame")
                v.mode = Mode::ConvectiveFlame;
            else
                Util::Abort(INFO, "surface_patch.mode must be 'prescribed_T' "
                                  "or 'convective_flame', got: ", mode_str);

            // S1: Robin discretisation for convective_flame
            // (see SurfaceCellFlux). `cell` (default) applies
            // h*(T_flame - T_cell) and cell-centred losses; `face` solves the
            // top-face balance for T_f and injects (2k/dz)*(T_f - T_cell).
            // (robin_form is read above, before the `enabled` return: P2-23)
            if (v.robin_form == RobinForm::Face && v.mode != Mode::ConvectiveFlame)
                Util::Abort(INFO, "surface_patch.robin_form = face requires "
                                  "surface_patch.mode = convective_flame");
            if (v.robin_form == RobinForm::Surface)
            {
                // P2-28: a PRESCRIBED patch on the face list keeps the list's form (its non-
                // prescribed faces are billed by FrFaceFlux at h = 0: their losses only)
                if (v.mode != Mode::ConvectiveFlame && !v.face_list)
                    Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                      "surface_patch.mode = convective_flame");
                // P2-28: pinned_idle_cycles is INERT (there is no pin rule: robin_form = pinned aborts above).
                // Read here, where the strict scope would otherwise abort, and reported inert.
                if (pp.contains("pinned_idle_cycles"))
                {
                    Set::Scalar pic = 0.0;
                    pp.query_required("pinned_idle_cycles", pic);
                    Util::Message(INFO, "surface_patch.pinned_idle_cycles = ", pic, " is INERT: the pin "
                                  "rule it governed was retired in P2-28 (robin_form = pinned)");
                }
            }

            // (x0, y0, radius are read above, before the `enabled` return: P2-23)

            // T_expr: required under prescribed_T, optional (consumed but
            // unused) under convective_flame.
            std::string T_expr;
            if (v.mode == Mode::PrescribedTemperature)
            {
                pp.query_required("T_expr", T_expr);
                v.T_parser = amrex::Parser(T_expr);
                v.T_parser.registerVariables({"t"});
                v.T_f = v.T_parser.compile<1>();
            }
            else
            {
                pp.query_default("T_expr", T_expr, std::string(""));
            }

            // h_conv / T_flame. Under convective_flame, h_conv is the
            // flame heat-transfer coefficient h_fl and T_flame the flame
            // temperature: q = h_fl*(T_flame - T_surface). Under prescribed_T
            // they are reserved/unused.

            // D2a: spatial flame, h_expr / T_flame_expr in (x, y, r, s, t) [SI],
            // r = distance of the column centre from (x0, y0), s = stand-off of
            // the nozzle z_n(t) = nozzle_z0 - nozzle_feed*t above the column's
            // mask top face. convective_flame only; both expressions or both
            // constants, never a mix. (That the expressions need
            // surface_patch.face_list = 1 is checked in CheckInterlocks.)
            const bool has_hx = pp.contains("h_expr");
            const bool has_Tx = pp.contains("T_flame_expr");
            // D2a2: wall-jet enthalpy closure. none (default) = D2a; enthalpy =
            // T_gas(r) marched outward from jet_T_nozzle, taking the face
            // absorption out of the jet (see BuildFlameColumns). Under enthalpy
            // h_expr is required and T_flame_expr is forbidden.
            std::string closure_str;
            pp.query_default("jet_closure", closure_str, std::string("none"));
            if      (closure_str == "none")     v.jet_enthalpy = false;
            else if (closure_str == "enthalpy") v.jet_enthalpy = true;
            else Util::Abort(INFO, "surface_patch.jet_closure must be 'none' or "
                                   "'enthalpy', got: ", closure_str);
            if (v.jet_enthalpy)
            {
                if (!has_hx)
                    Util::Abort(INFO, "surface_patch.jet_closure = enthalpy requires "
                                      "surface_patch.h_expr");
                if (has_Tx)
                    Util::Abort(INFO, "surface_patch.T_flame_expr is not read under "
                                      "jet_closure = enthalpy (T_gas is marched); remove it");
                pp.query_required("jet_mdot",     v.jet_mdot);
                pp.query_required("jet_T_nozzle", v.jet_T_nozzle);
                pp.query_required("jet_cp",       v.jet_cp);
                pp.query_required("jet_D",        v.jet_D);
                pp.query_default("jet_T_ent",     v.jet_T_ent, 293.15);
                // Free-jet centreline decay at the stagnation bin (default) or
                // the nozzle temperature there; optionally carry the diluted mass.
                // (query_default adds the key, so read "given" first; D2c)
                const bool core_given = pp.contains("jet_core_length");
                pp.query_default("jet_core_length", v.jet_core_length, 5.0);
                std::string stag_str;
                pp.query_default("jet_stagnation", stag_str, std::string("decay"));
                if      (stag_str == "decay")  v.jet_stag_decay = true;
                else if (stag_str == "nozzle") v.jet_stag_decay = false;
                else Util::Abort(INFO, "surface_patch.jet_stagnation must be 'decay' or "
                                       "'nozzle', got: ", stag_str);
                int ent = 0;
                pp.query_default("jet_entrained_mass", ent, 0);
                if (ent != 0 && ent != 1)
                    Util::Abort(INFO, "surface_patch.jet_entrained_mass must be 0 or 1, got: ", ent);
                v.jet_entrained_mass = (ent == 1);
                // D2b: temperature of the entrained gas at the stagnation mix.
                // fixed (default) = jet_T_ent; exhaust = the previous step's
                // jet_T_exhaust (blind-hole recirculation; needs entrained mass).
                std::string tem;
                pp.query_default("jet_T_ent_mode", tem, std::string("fixed"));
                if      (tem == "fixed")     v.jet_T_ent_exhaust = false;
                else if (tem == "exhaust")   v.jet_T_ent_exhaust = true;
                else Util::Abort(INFO, "surface_patch.jet_T_ent_mode must be 'fixed' or "
                                       "'exhaust', got: ", tem);
                if (v.jet_T_ent_exhaust && !v.jet_entrained_mass)
                    Util::Abort(INFO, "surface_patch.jet_T_ent_mode = ", tem, " requires "
                                      "surface_patch.jet_entrained_mass = 1");
                // D2c: decay diameter. nozzle (default) = jet_D; momentum =
                // jet_De_ref*sqrt(T_mix/jet_T_ent), jet_De_ref being the FULL
                // jet's momentum diameter 2*mdot/sqrt(pi*rho*J) in surroundings
                // at jet_T_ent (not scaled with the sector's jet_mdot share).
                // Ricou-Spalding entrainment = momentum + jet_core_length =
                // 3.125 + jet_entrained_mass = 1 (no separate key).
                std::string dd;
                pp.query_default("jet_decay_diameter", dd, std::string("nozzle"));
                if      (dd == "nozzle")   v.jet_De_momentum = false;
                else if (dd == "momentum") v.jet_De_momentum = true;
                else Util::Abort(INFO, "surface_patch.jet_decay_diameter must be 'nozzle' or "
                                       "'momentum', got: ", dd);
                if (v.jet_De_momentum)
                {
                    if (!v.jet_stag_decay)
                        Util::Abort(INFO, "surface_patch.jet_decay_diameter = momentum requires "
                                          "surface_patch.jet_stagnation = decay");
                    if (!core_given)
                        Util::Abort(INFO, "surface_patch.jet_decay_diameter = momentum requires "
                                          "an explicit surface_patch.jet_core_length (units of D_e)");
                    pp.query_required("jet_De_ref", v.jet_De_ref);
                    if (!(std::isfinite(v.jet_De_ref) && v.jet_De_ref > 0.0))
                        Util::Abort(INFO, "surface_patch.jet_De_ref must be finite and > 0, got: ",
                                    v.jet_De_ref);
                }
                else if (pp.contains("jet_De_ref"))
                    Util::Abort(INFO, "surface_patch.jet_De_ref is only read under "
                                      "jet_decay_diameter = momentum");
                // D2c: free-surface dilution (needs nozzle_descent = feet,
                // checked below once the descent mode is parsed).
                int fs = 0;
                pp.query_default("jet_free_surface", fs, 0);
                if (fs != 0 && fs != 1)
                    Util::Abort(INFO, "surface_patch.jet_free_surface must be 0 or 1, got: ", fs);
                v.jet_free_surface = (fs == 1);
                if (v.jet_free_surface)
                {
                    pp.query_default("jet_fs_aspect",   v.jet_fs_aspect,   1.0);
                    pp.query_default("jet_fs_exponent", v.jet_fs_exponent, 1.0);
                    if (!(std::isfinite(v.jet_fs_aspect) && v.jet_fs_aspect >= 0.0 &&
                          std::isfinite(v.jet_fs_exponent) && v.jet_fs_exponent >= 0.0))
                        Util::Abort(INFO, "surface_patch.jet_fs_aspect / jet_fs_exponent must be "
                                          "finite and >= 0");
                }
                else if (pp.contains("jet_fs_aspect") || pp.contains("jet_fs_exponent"))
                    Util::Abort(INFO, "surface_patch.jet_fs_aspect / jet_fs_exponent are only "
                                      "read under jet_free_surface = 1");
                std::string bu;
                pp.query_default("jet_bin_update", bu, std::string("explicit"));
                if      (bu == "explicit")    v.jet_bin_exponential = false;
                else if (bu == "exponential") v.jet_bin_exponential = true;
                else Util::Abort(INFO, "surface_patch.jet_bin_update must be 'explicit' or "
                                       "'exponential', got: ", bu);
                std::string nf;
                pp.query_default("jet_negative_flux", nf, std::string("drop"));
                if      (nf == "drop")  v.jet_neg_count = false;
                else if (nf == "count") v.jet_neg_count = true;
                else Util::Abort(INFO, "surface_patch.jet_negative_flux must be 'drop' or "
                                       "'count', got: ", nf);
                pp.query_default("jet_dr",        v.jet_dr,    0.0); // 0 -> max(dx, dy)
                // Optional per-bin profile log <plot_file>_jet_profile.csv every
                // jet_profile_interval seconds of simulated time (0 = off).
                pp.query_default("jet_profile_interval", v.jet_profile_interval, 0.0);
                if (!(std::isfinite(v.jet_mdot) && v.jet_mdot > 0.0 &&
                      std::isfinite(v.jet_cp) && v.jet_cp > 0.0 &&
                      std::isfinite(v.jet_T_ent) && v.jet_T_ent > 0.0 &&
                      std::isfinite(v.jet_T_nozzle) && v.jet_T_nozzle >= v.jet_T_ent &&
                      std::isfinite(v.jet_D) && v.jet_D > 0.0 &&
                      std::isfinite(v.jet_core_length) && v.jet_core_length > 0.0 &&
                      std::isfinite(v.jet_dr) && v.jet_dr >= 0.0 &&
                      std::isfinite(v.jet_profile_interval) && v.jet_profile_interval >= 0.0))
                    Util::Abort(INFO, "surface_patch.jet_*: need jet_mdot > 0, jet_cp > 0, "
                                      "jet_T_nozzle >= jet_T_ent > 0, jet_D > 0, "
                                      "jet_core_length > 0, jet_dr >= 0, "
                                      "jet_profile_interval >= 0 (all finite)");
            }
            else
            {
                for (const char* k : {"jet_decay_diameter", "jet_De_ref", "jet_free_surface",
                                      "jet_fs_aspect", "jet_fs_exponent", "jet_bin_update",
                                      "jet_negative_flux"})
                    if (pp.contains(k))
                        Util::Abort(INFO, "surface_patch.", k, " requires "
                                          "surface_patch.jet_closure = enthalpy");
            }
            // rev7 d2: per-column side-power log (EnergyLedger), seconds of
            // simulated time between rows. 0 = off (no file).
            // P2-26 (census B1): parsed under EVERY jet_closure -- its
            // consumer is closure-independent, so a jet_closure = none deck
            // can set it too.
            pp.query_default("side_power_interval", v.side_power_interval, 0.0);
            if (!(std::isfinite(v.side_power_interval) &&
                  v.side_power_interval >= 0.0))
                Util::Abort(INFO, "surface_patch.side_power_interval must be "
                                  "finite and >= 0, got ", v.side_power_interval);
            // D2b: nozzle path. prescribed (default) = D2a nozzle_z0 -
            // nozzle_feed*t; feet = the burner rests on a foot annulus and
            // descends with the rock under it (see FeetDescent).
            std::string desc;
            pp.query_default("nozzle_descent", desc, std::string("prescribed"));
            if      (desc == "prescribed") v.feet = false;
            else if (desc == "feet")       v.feet = true;
            else Util::Abort(INFO, "surface_patch.nozzle_descent must be 'prescribed' or "
                                   "'feet', got: ", desc);
            if (v.feet && !(has_hx || has_Tx))
                Util::Abort(INFO, "surface_patch.nozzle_descent = feet requires the flame "
                                  "expressions (surface_patch.h_expr)");
            if (v.jet_free_surface && !v.feet)
                Util::Abort(INFO, "surface_patch.jet_free_surface = 1 requires "
                                  "surface_patch.nozzle_descent = feet");
            if (!v.feet && pp.contains("foot_body_clearance"))
                Util::Abort(INFO, "surface_patch.foot_body_clearance requires "
                                  "surface_patch.nozzle_descent = feet");
            if (has_hx || has_Tx)
            {
                if (v.mode != Mode::ConvectiveFlame)
                    Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr require "
                                      "surface_patch.mode = convective_flame");
                if (!(has_hx && (has_Tx || v.jet_enthalpy)))
                    Util::Abort(INFO, "surface_patch.h_expr and surface_patch.T_flame_expr "
                                      "must be given together");
                if (pp.contains("h_conv") || pp.contains("T_flame"))
                    Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr cannot be combined "
                                      "with surface_patch.h_conv / T_flame");
                std::string h_expr, Tg_expr;
                pp.query_required("h_expr",       h_expr);
                v.h_parser = amrex::Parser(h_expr);
                v.h_parser.registerVariables({"x", "y", "r", "s", "t"});
                v.h_f = v.h_parser.compile<5>();
                if (!v.jet_enthalpy)
                {
                    pp.query_required("T_flame_expr", Tg_expr);
                    v.Tg_parser = amrex::Parser(Tg_expr);
                    v.Tg_parser.registerVariables({"x", "y", "r", "s", "t"});
                    v.Tg_f = v.Tg_parser.compile<5>();
                }
                if (v.feet)
                {
                    if (pp.contains("nozzle_z0") || pp.contains("nozzle_feed"))
                        Util::Abort(INFO, "surface_patch.nozzle_z0 / nozzle_feed cannot be "
                                          "given with nozzle_descent = feet");
                    pp.query_required("foot_r_inner",  v.foot_r_inner);
                    pp.query_required("foot_r_outer",  v.foot_r_outer);
                    pp.query_required("foot_standoff", v.foot_standoff);
                    if (!(std::isfinite(v.foot_r_inner) && std::isfinite(v.foot_r_outer) &&
                          v.foot_r_inner >= 0.0 && v.foot_r_inner < v.foot_r_outer))
                        Util::Abort(INFO, "surface_patch.foot_r_inner / foot_r_outer: need "
                                          "0 <= foot_r_inner < foot_r_outer (finite)");
                    if (!(std::isfinite(v.foot_standoff) && v.foot_standoff > 0.0))
                        Util::Abort(INFO, "surface_patch.foot_standoff must be finite and > 0");
                    pp.query_default("foot_npads", v.foot_npads, 3);
                    if (v.foot_npads < 1)
                        Util::Abort(INFO, "surface_patch.foot_npads must be >= 1");
                    pp.query_default("foot_pad_quantile", v.foot_pad_quantile, 0.9);
                    if (!(v.foot_pad_quantile > 0.0 && v.foot_pad_quantile <= 1.0))
                        Util::Abort(INFO, "surface_patch.foot_pad_quantile must be in (0, 1]");
                    pp.query_default("nozzle_feed_max", v.nozzle_feed_max,
                                     std::numeric_limits<Set::Scalar>::infinity());
                    if (!(v.nozzle_feed_max > 0.0))
                        Util::Abort(INFO, "surface_patch.nozzle_feed_max must be > 0");
                    // D2c: burner-body clearance (rock above a pad is removed
                    // mechanically at the next removal pass; see FeetDescent).
                    int cl = 0;
                    pp.query_default("foot_body_clearance", cl, 0);
                    if (cl != 0 && cl != 1)
                        Util::Abort(INFO, "surface_patch.foot_body_clearance must be 0 or 1, got: ", cl);
                    v.foot_body_clearance = (cl == 1);
                }
                else
                {
                    pp.query_required("nozzle_z0", v.nozzle_z0);
                    pp.query_default("nozzle_feed", v.nozzle_feed, 0.0);
                    if (!(std::isfinite(v.nozzle_z0) && std::isfinite(v.nozzle_feed)))
                        Util::Abort(INFO, "surface_patch.nozzle_z0 / nozzle_feed must be finite");
                }
                // Columns with r <= nozzle_collision_radius abort on s <= 0 (the
                // nozzle is in the rock). Default = the patch radius. Outside it a
                // column may sit above the descending nozzle (unrecessed rock
                // beside the hole); s <= 0 is then passed to the expressions,
                // which must stay finite and > 0 (under jet_closure = enthalpy
                // h is 0 there and h_expr is not evaluated).
                // D2b: under feet the default is foot_r_inner (inside the collar).
                pp.query_default("nozzle_collision_radius", v.nozzle_collision_radius,
                                 v.feet ? v.foot_r_inner : v.radius);
                if (!(std::isfinite(v.nozzle_collision_radius) && v.nozzle_collision_radius >= 0.0))
                    Util::Abort(INFO, "surface_patch.nozzle_collision_radius must be finite and >= 0");
                v.use_expr = true;
                return;
            }

            pp.query_default("h_conv",  v.h_conv,  0.0);
            pp.query_default("T_flame", v.T_flame, 0.0);
            if (v.mode == Mode::ConvectiveFlame)
            {
                if (!(v.h_conv > 0.0))
                    Util::Abort(INFO, "surface_patch.h_conv (h_fl) must be > 0 "
                                      "for convective_flame mode");
                if (!(v.T_flame > 0.0))
                    Util::Abort(INFO, "surface_patch.T_flame must be > 0 for "
                                      "convective_flame mode");
            }
        }
    };
    SurfacePatch surface_patch;
