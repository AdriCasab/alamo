    static void ParseSpallSettings(MMWSpalling& value, IO::ParmParse& pp)
    {
        int on = 0;
        pp.query_default("spall.enabled", on, 0);
        value.spall_enabled = (on != 0);

        // Stress source for the Sp = K_I/K_Ic criterion. Parsed BEFORE the
        // spall.enabled early-return because it also governs the Sp_field
        // ONSET diagnostic (site A (deleted, P2-28)), which runs under
        // spallation.model = sp_weibull even when material removal is off
        // (e.g. the Kant onset verification runs). Default `mechanics`
        // (stress_mf from the elastic MLMG solve) preserves existing
        // behaviour. `local_thermoelastic` reconstructs sigma_xx from the
        // temperature field (Kant 1-D confinement), letting el.type = disable.
        std::string stress_src;
        pp.query_default("spall.stress_source", stress_src,
                         std::string("mechanics"));
        if (stress_src == "mechanics")
            value.sp_local_thermoelastic = false;
        else if (stress_src == "local_thermoelastic")
            value.sp_local_thermoelastic = true;
        else
            Util::Abort(INFO, "spall.stress_source must be 'mechanics' or "
                              "'local_thermoelastic', got ", stress_src);

        // ⚠⚠ P2-2b: the MECHANICAL half of the measured property set.
        //   scalar  (default) -- the criterion reads the per-phase E and beta
        //                        scalars, i.e. the pre-P2-2b arithmetic.
        //   grimsel           -- the criterion takes E(T) and the accumulated
        //                        thermal strain from the law, through the one
        //                        `LawSigma` entry point.
        // ⚠ This is SEPARATE from `microstructure.phaseN.law`, which is the
        // THERMAL half (lambda, rho, c_p into conduction). One key for both
        // would make the Kant ladder's K1 and K2 rungs identical and destroy
        // (c1)'s attribution -- the reason this key exists.
        // ⚠ `grimsel` here without the phase law is legal and meaningful: it
        // is the mechanical set ALONE, which is what separates the two halves.
        std::string stress_law;
        pp.query_default("spall.stress_law", stress_law,
                         std::string("scalar"));
        if (stress_law == "scalar")
            value.sp_stress_law_grimsel = false;
        else if (stress_law == "grimsel")
            value.sp_stress_law_grimsel = true;
        else
            Util::Abort(INFO, "spall.stress_law must be 'scalar' or "
                              "'grimsel', got ", stress_law);

        // A1: every-removal event log (validated at the end of Parse, once
        // spallation.model is known). P2-28: the clusters and h_col CSVs are retired with the column path.
        int csv_on = 0;
        pp.query_default("spall.removal_events_csv", csv_on, 0);
        value.sp_emit_removal_events_csv = (csv_on != 0);

        // P0b deliverable 3: the flake thickness is c_f * ell_eff, with NO dz
        // clamp. ⚠ The min(., dz) clamp P0 used is the measured culprit for the
        // 20.6 % residual in its mesh ladder: at 0.5 mm, ell_eff (0.992 mm) is
        // already larger than dz, so the "flake" was the whole cell again and
        // T_bar collapsed back onto the cell average.
        // ⚠ c_f IS THE ONE EMPIRICAL CONSTANT THIS PACKET INTRODUCES. Under
        // `flake` the thickness is IMPOSED (the column path (deleted, P2-28) bypasses the depth scan's
        // h_col), so it is not derived from fracture here. 0.3-0.5 is what
        // Wilkinson's 100-1000 um spall thicknesses and Rauenzahn's size
        // distributions imply at ell_eff = 0.8-1.6 mm. Its lever on the level is
        // weak (about +-5 % across 0.3-0.5), which is the real argument that it
        // is not a tuning knob; the successor is to run the K_I(z) scan on the
        // reconstructed exponential so c_f becomes an OUTPUT (P1).
        pp.query_default("spall.flake_c_f", value.sp_flake_c_f, 0.4);
        if (!(std::isfinite(value.sp_flake_c_f) && value.sp_flake_c_f > 0.0))
            Util::Abort(INFO, "spall.flake_c_f must be finite and > 0, got ",
                        value.sp_flake_c_f);
        // Step 16: the depth the Sp-vs-depth scan searches (10 mm safety cap).
        // P2-26 (census B1): parsed HERE, before the flake_h_max check below and
        // before the spall-off return, so the check reads the DECK's value -- it
        // used to read the member's initialiser, because this parse ran after it.
        // the column reconstruction (deleted, P2-28) reads it on a spall-off onset deck as well.
        pp.query_default("spall.h_spall_max", value.spall_h_spall_max, 0.010);
        if (!(std::isfinite(value.spall_h_spall_max) &&
              value.spall_h_spall_max > 0.0))
            Util::Abort(INFO, "spall.h_spall_max must be finite and > 0, got ",
                        value.spall_h_spall_max);
        // Revision 7 deliverable 1. Default 0 = INERT, so every frozen leg and
        // every control is byte-identical without it. See ReconFromX for what
        // this bound is and, more importantly, what it is not.
        pp.query_default("spall.flake_h_max", value.sp_flake_h_max, 0.0);
        if (!(std::isfinite(value.sp_flake_h_max) && value.sp_flake_h_max >= 0.0))
            Util::Abort(INFO, "spall.flake_h_max must be finite and >= 0 "
                              "(0 = inert), got ", value.sp_flake_h_max);
        // (its cross-key check against spall.h_spall_max is in CheckInterlocks, P2-27)

        if (!value.spall_enabled) return;

        // Step 16: the Sp-vs-depth scan's segment count (default sp_n_segments);
        // only consumed under spallation.model = sp_weibull. Its depth,
        // spall.h_spall_max, is parsed above the flake_h_max check (P2-26).
        pp.query_default("spall.h_spall_n_segments",
                         value.spall_h_spall_n_segments, 0);
        if (value.spall_h_spall_n_segments < 0)
            Util::Abort(INFO, "spall.h_spall_n_segments must be >= 0, got ",
                        value.spall_h_spall_n_segments);

        // D2q-2b item 0: DEBUG-ONLY seeded interior void. Six ints in INDEX
        // space, `ilo jlo klo ihi jhi khi`, flipped to void once after
        // Initialize. Absent = off, and nothing in the code path runs.
        //
        // It exists so the paths that must tolerate SOLID ABOVE VOID -- the
        // depth-scan clamp, the three k_top scans, ColumnSolidAboveCount and
        // beam invariant (a) -- can be exercised at their real call sites
        // before any mechanism creates such a state. D2q-1 shipped that
        // handling with nothing able to produce an interior void, so it was
        // proved inert but never proved correct.
        //
        // Deliberately writes NO diagnostic column: a new thermo.dat or
        // removal_events.csv field breaks the 2639-file identity on headers
        // alone, whatever the numbers do.

        // NOTE the `spall.` prefix: ParseSpallSettings receives the TOP-level
        // ParmParse, so every key in this function carries its own prefix.
        // Writing "seed_void" here parses nothing and AMReX only *warns*
        // ("Unused ParmParse Variables"), so the run proceeds with the feature
        // silently off.
        std::vector<int> sv;
        pp.queryarr("spall.seed_void", sv);
        if (!sv.empty())
        {
            // P1b: SEVERAL boxes, 6 ints each. Six ints is exactly what it
            // was, so every existing input is unchanged. ⚠ The extension is
            // not cosmetic: a standing column with void on two sides CANNOT be
            // built from one box -- if (i+1, j) and (i, j+1) are both inside a
            // box then so is (i, j) -- so a single box can express a pit, a
            // ridge or a step, but never the tower the support rule exists to
            // remove. The unit test needs that geometry prescribed, not grown.
            if (sv.size() % 6 != 0)
                Util::Abort(INFO, "spall.seed_void takes a MULTIPLE of 6 ints "
                            "(ilo jlo klo ihi jhi khi per box), got ",
                            sv.size());
            for (std::size_t b = 0; b < sv.size(); b += 6)
                for (int d = 0; d < 3; ++d)
                    if (sv[b + d] > sv[b + d + 3])
                        Util::Abort(INFO, "spall.seed_void box ", b / 6,
                                    ": lo > hi in direction ", d,
                                    ": ", sv[b + d], " > ", sv[b + d + 3]);
            value.spall_seed_void = sv;
        }
    }

    static void ParseVaporSettings(MMWSpalling& value, IO::ParmParse& pp)
    {
        int on = 0;
        pp.query_default("vapor.enabled", on, 0);
        value.vapor_enabled = (on != 0);
        if (!value.vapor_enabled) return;

        if (pp.contains("vapor.A_beam"))
        {
            value.vapor_A_beam_override = true;
            pp.query_required("vapor.A_beam", value.vapor_A_beam);
            if (!(std::isfinite(value.vapor_A_beam) && value.vapor_A_beam > 0.0))
                Util::Abort(INFO, "vapor.A_beam must be finite and > 0, got ",
                            value.vapor_A_beam);
        }
        if (pp.contains("vapor.T0"))
        {
            value.vapor_T0_override = true;
            pp.query_required("vapor.T0", value.vapor_T0);
            if (!std::isfinite(value.vapor_T0))
                Util::Abort(INFO, "vapor.T0 must be finite, got ", value.vapor_T0);
        }

        int presc = 0;
        pp.query_default("vapor.prescribed_rop.enabled", presc, 0);
        value.vapor_prescribed_rop_enabled = (presc != 0);
        if (value.vapor_prescribed_rop_enabled)
        {
            std::string expr;
            pp.query_required("vapor.prescribed_rop.expr", expr);
            value.vapor_rop_parser = amrex::Parser(expr);
            value.vapor_rop_parser.registerVariables({"x", "y", "z", "t"});
            value.vapor_rop_f = value.vapor_rop_parser.compile<4>();
        }
    }

    // ====================================================================
    // P2-2b deliverable 4: THE ONE THERMOELASTIC STRESS ENTRY POINT.
    // ====================================================================
    // ⚠⚠ EVERY law-path stress in the rock goes through here, so
    // `Grimsel::sigma_thermal` is called from exactly ONE place in src/ and
    // `unit/lateral_flux_lock` can assert that as a fact rather than a habit.
    // Grimsel.H's header records three independent WRONG answers for this one
    // stress, each from a site that wrote its own formula:
    //   * `alpha(T_s)*dT`            OVERSTATES  (alpha triples 300->500 degC)
    //   * accumulated strain, cold E UNDERSTATES (E halves by 600 degC)
    //   * `INTEGRAL(E*alpha dtheta)` wrong outright -- elasticity has no path
    //     memory; stress is the CURRENT stiffness times the CURRENT strain.
    // ⇒ accumulate the STRAIN, multiply by the stiffness AT T_hot, which is
    // what `sigma_thermal` does and what `(b2)` scores to 1e-12 against an
    // independent twin BETWEEN table nodes as well as on them.
    //
    // ⚠ `static` and taking the law by POINTER on purpose. Every caller --
    // the `Sp_field` diagnostic site (an `amrex::LoopOnCpu` with a `[&]`
    // capture today; it was a `ParallelFor` under `[=, this]` when this was
    // written) and the depth scans -- hoists one pointer before its loop; a
    // 15.7 KB by-value copy of the two property tables would be carried into
    // every closure.
    //
    // ⚠⚠ THE SIGN. Returns COMPRESSION-POSITIVE for T_hot > T_cold, matching
    // `SpCriterion::K_I`'s convention and the double-negation the pre-P2-2b
    // expressions perform. The column scan's stress_mf branch (the column path (deleted, P2-28))
    // returns `-sg(0,0)` for the same
    // reason. A reversed bracket here would discard every hot flake and stop
    // all drilling -- it is what rev 4 of this packet would have shipped.
    //
    // ⚠ Under the law, `nu` is the LAW's (parsed, 0.30 by default), not the
    // deck's `E/(2 mu) - 1`. On both production decks those agree at 0.30; if
    // a deck ever disagrees, the law wins on the law path by construction.
    static Set::Scalar LawSigma(const Numeric::Material::Grimsel *g,
                                Set::Scalar T_hot, Set::Scalar T_cold)
    {
        return g->sigma_thermal(T_hot, T_cold);
    }

    // The rock law for the criterion, or nullptr under the shipped `constant`
    // default. ⚠ ONE place decides this, so "is the law on?" cannot be
    // answered differently at two of the three stress sites -- which is the
    // failure mode `MMWSpalling.H`'s two-firing-site invariant exists to stop.
    // ⚠ Phase 0 by construction: the T-derived sigma_xx path is already
    // phase-0-only (it aborts on an empty `phases`), and widening that is not
    // in this packet.
    const Numeric::Material::Grimsel *RockLaw() const
    {
        // ⚠⚠ GATED ON `spall.stress_law`, NOT on the phase's thermal law.
        // Gating this on `phases[0].grimsel_law` -- which is what the first
        // implementation did -- ties the criterion to conduction, so the Kant
        // ladder's K1 (thermal only) and K2 (+ mechanical) become the same
        // run and (c1) cannot attribute the onset shift to either half.
        // ⚠ The law OBJECT still comes from phase 0, because that is where the
        // parsed constants live; what this key decides is only whether the
        // criterion consults it.
        return (sp_stress_law_grimsel && !phases.empty())
             ? &phases[0].grimsel : nullptr;
    }
