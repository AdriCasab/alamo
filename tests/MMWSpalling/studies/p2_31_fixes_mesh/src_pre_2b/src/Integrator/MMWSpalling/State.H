    // Step 5: Voronoi mineral microstructure. `phase_mf` is a cell-centred
    // per-cell mineral phase ID (Set::Scalar storage; integers cast to float).
    // The property fields cache the heterogeneous κ / β / E / ρ / Cp lookup so
    // later steps can read them as ordinary AMReX MultiFabs without
    // re-querying the phase table on every kernel.
    struct MineralPhase
    {
        Set::Scalar fraction = 0.0;
        Set::Scalar kappa = 0.0;
        Set::Scalar beta  = 0.0;
        Set::Scalar E     = 0.0;
        Set::Scalar mu    = 0.0;
        Set::Scalar rho   = 0.0;
        Set::Scalar Cp    = 0.0;
        // Step 6b: phase-change scalars used to build the per-phase
        // H<->T table. Defaults disable phase change.
        Set::Scalar T_ref     = 273.15;
        Set::Scalar T_m       = 1.0e30;
        Set::Scalar fus_width = 0.0;
        Set::Scalar L_m       = 0.0;
        Set::Scalar T_vap_lo  = 1.0e30;
        Set::Scalar T_vap_hi  = 1.0e30;
        Set::Scalar L_v       = 0.0;
        // P2-2 change 2: the per-phase property LAW.
        //   constant (DEFAULT) -- the scalar kappa/beta/E/mu/rho/Cp above, all
        //                         six required.
        //   grimsel            -- the MEASURED Grimsel / Central Aare laws,
        //                         Numeric::Material::Grimsel, including the
        //                         tabulated thermal strain (change 5).
        // ⚠ Under a law the scalars it fills are REFUSED at the parse and stay
        // 0 here: kappa/rho/Cp under grimsel or zhang, and beta/E/mu as well
        // under grimsel with spall.stress_law = grimsel (Parse, P2-27 stage 4).
        bool grimsel_law = false;
        Numeric::Material::Grimsel grimsel;
        //   zhang (P2-29)      -- Zhang's granite (Zhang Table 1, Numeric::Material::Zhang): the
        //                         THERMAL scalars only (kappa a conductivity); the mechanics stay scalars.
        bool zhang_law = false;
        Numeric::Material::Zhang zhang;
        bool alpha_beta_enabled = false;
        Set::Scalar alpha_beta_T0 = 846.15;
        Set::Scalar alpha_beta_width = 10.0;
        Set::Scalar alpha_beta_dV = 0.0063;
        Set::Scalar alpha_beta_linear_strain = 0.0;
    };
    enum class MicroMode { Voronoi, Expression };
    std::vector<int> dt_limit_checked;   // P2-27: per level, CheckTimeStepLimit ran
    int  num_grains  = 0;
    int  micro_seed  = 1;
    MicroMode micro_mode = MicroMode::Voronoi;
    amrex::Parser            phase_parser;
    amrex::ParserExecutor<3> phase_f;
    // Step 6c: optional independent grain-id parser (expression mode only).
    amrex::Parser            grain_parser;
    amrex::ParserExecutor<3> grain_f;
    bool                     has_grain_expr = false;
    std::vector<MineralPhase> phases;
    // Step 6b: per-phase Constant material + H<->T table. One per phase,
    // index-aligned with `phases`. Built once, on the first
    // InitializeMicrostructure call.
    std::vector<Numeric::Material::Constant> phase_materials;
    std::vector<Numeric::Material::Table>    phase_tables;
    std::vector<int>          grain_phase_ids;
    std::vector<Set::Scalar>  grain_centres; // packed (x,y,z) per grain

    Set::Field<Set::Scalar> phase_mf;
    // Per-level flag: phase/grain_id/kappa ghost layers filled since the last
    // InitializeMicrostructure (see AdvanceMicrostructure).
    std::vector<char>       micro_ghost_filled;
    Set::Field<Set::Scalar> is_gb_mf;
    // Step 6c: explicit grain topology. grain_id distinguishes individual
    // grains (including same-mineral neighbours); is_grain_boundary marks
    // cells with a face neighbour of a different grain_id;
    // is_phase_boundary marks cells with a face neighbour of a different
    // mineral phase. is_gb is preserved as an alias for is_phase_boundary.
    Set::Field<Set::Scalar> grain_id_mf;
    Set::Field<Set::Scalar> is_grain_boundary_mf;
    Set::Field<Set::Scalar> is_phase_boundary_mf;
    Set::Field<Set::Scalar> kappa_mf, beta_mf, E_phase_mf, rho_mf, Cp_mf;
    // Step 7: per-phase shear modulus and reference temperature for the
    // heterogeneous thermoelastic update. Cell-centred, 1-ghost; mu is
    // averaged onto the nodes in UpdateModel (CellToNodeAverageInDomain),
    // T_ref enters the cell-centred free strain (ThermalFreeStrain).
    Set::Field<Set::Scalar> mu_phase_mf, T_ref_phase_mf;
    // P2-2a change A: the stored initial temperature -- the stress-free state.
    // Written at exactly ONE site (Initialize) and never again; (b9) locks it.
    Set::Field<Set::Scalar> T_init_mf;
    // Optional quartz alpha-beta transformation strain and the total scalar
    // free strain beta*(T-T_ref) + transition strain. Cell-centred, 1 ghost.
    Set::Field<Set::Scalar> alpha_beta_linear_strain_mf;
    Set::Field<Set::Scalar> alpha_beta_T0_mf;
    Set::Field<Set::Scalar> alpha_beta_width_mf;
    Set::Field<Set::Scalar> thermal_free_strain_mf;

    // Step 6c: the grain-boundary contact conductance (MICRO; 0 = off).
    Set::Scalar h_gb_conductance   = 0.0;

    // Step 10: surface advancement and spall detachment. Default-off.
    // `phi = z_surface(x,y) - z_cell` is the column's level-set state:
    // z_surface starts at the domain top and drops by each recession of the
    // column's +z line, vaporised top cell or pad cut (FaceRemoval.H).
    // ⚠ It is NOT the void test: `removed` is the binary void mask, the one
    // every consumer reads (a top cell more than half consumed is solid with
    // phi < 0; a wall line's void and a seeded void do not move phi).
    // `spall_event`, `spall_thickness`, and `RoP_spall` are per-step
    // diagnostics reset at the start of each removal step.
    bool spall_enabled = false;

    // D2q-2b item 0: debug-only seeded interior void, `ilo jlo klo ihi jhi khi`
    // per box, in index space. Empty = off (the default), and ApplySeedVoid returns
    // immediately. See the parse block for why it carries no diagnostic.
    std::vector<int> spall_seed_void;

    // Step 16: Sp-vs-depth h_spall scan parameters (sp_weibull branch only).
    // `spall.h_spall_max` is the depth the scan searches, in metres normal to
    // the surface (a line's window is the length h_spall_max / mu along it,
    // FaceRemoval.H), and the cap on the flake thickness (ReconFromX). The
    // segment count `spall.h_spall_n_segments` falls back to `sp_n_segments`
    // when set to 0 (default): the Simpson resolution of the K_I integral.
    Set::Scalar spall_h_spall_max         = 0.010;
    int         spall_h_spall_n_segments  = 0;
    // Stress source for the Sp = K_I/K_Ic criterion (spall.stress_source).
    // `local_thermoelastic` forms the 1-D-confined stress algebraically from
    // the temperature field, compression-positive E*beta*(T - T0)/(1 - nu) or
    // the rock law's (ScanSigmaT; T0 the stored initial temperature), so the
    // criterion runs with mechanics disabled (el.type = disable) — no MLMG
    // solve, no convergence cliff. Valid
    // for the brittle-spallation regime where the constrained heated skin is
    // the dominant stress driver (Kant/von Rohr formulation). Loses grain-scale
    // stress heterogeneity and cavity stress redistribution; see
    // mlmg-shelved-local-thermoelastic-notes.md.
    // ⚠ The face list -- the one criterion site -- carries ONLY this source:
    // the parse default `mechanics` (the elastic solve's stress_mf) aborts
    // under surface_patch.face_list (CheckInterlocks).
    bool sp_local_thermoelastic = false;
    // ⚠⚠ P2-2b: THE MECHANICAL HALF OF THE GRIMSEL SET, KEYED SEPARATELY FROM
    // THE THERMAL HALF. `microstructure.phaseN.law = grimsel` routes lambda,
    // rho and c_p into conduction; this routes E and the thermal strain into
    // the CRITERION. They must be independently switchable or the Kant ladder
    // collapses -- `K1` (thermal only) and `K2` (+ mechanical) would be the
    // same run, and `(c1)` could not say which property moved the onset.
    // Default false ⇒ the criterion keeps its pre-P2-2b scalar arithmetic bit
    // for bit even when the thermal law is on.
    bool sp_stress_law_grimsel = false;

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> removed_mf;
    Set::Field<Set::Scalar> spall_event_mf;
    Set::Field<Set::Scalar> spall_thickness_mf;
    Set::Field<Set::Scalar> RoP_spall_mf;

    // Step 11: vaporisation removal. Default-off and shares phi / removed
    // with the spall path. The rule is the face list's (P2-24, FaceRemoval.H):
    // a LISTED cell whose enthalpy has passed its phase table's H_vap_hi
    // leaves whole, as regime 2; its lines are not scanned that step, and
    // RoP_vap is the cell's remaining solid thickness / dt.
    bool vapor_enabled = false;

    Set::Field<Set::Scalar> vapor_event_mf;
    Set::Field<Set::Scalar> vapor_thickness_mf;
    Set::Field<Set::Scalar> RoP_vap_mf;
    Set::Field<Set::Scalar> regime_mf;
    Set::Field<Set::Scalar> RoP_mf;
    // Step 17: cell-resolved per-firing-cell regime diagnostic under
    // sp_weibull. {0 none, 1 spall, 2 vap, 4 pad clearance}. Registered only when
    // spallation_model == SpWeibull && (spall_enabled || vapor_enabled).
    // Distinct from the scalar `regime_mf` (registered earlier; kept).
    Set::Field<Set::Scalar> regime_field_mf;

    bool removal_enabled() const { return spall_enabled || vapor_enabled; }
    // ⚠ P2-28: THE FACE LIST'S EVALUATE-ONLY MODE -- the list on with spall and vapour off (the onset decks: Kant, the
    // pressure sweep, Hu). The state the list reads (phi, removed, regime, RoP) exists and the criterion runs on every
    // listed line; nothing is removed. Removal itself requires the list (CheckInterlocks), so list_state() is true
    // exactly under surface_patch.face_list.
    bool list_state() const { return removal_enabled() || surface_patch.face_list; }
    bool evaluate_only() const { return surface_patch.face_list && !removal_enabled(); }

    // Step 12: spallation-model toggle. `none` (default) loads no criterion;
    // `sp_weibull` activates the LEFM Sp = K_I/K_Ic(T) onset criterion.
    enum class SpallationModel { None, SpWeibull };
    SpallationModel spallation_model = SpallationModel::None;

    // Sp criterion parameters (sp_weibull branch only). Deterministic
    // surface flaw length a0 (m); under weibull.enabled the per-cell flaw
    // field (flaw_a_mf, Step 15b) is read instead.
    Set::Scalar sp_a0 = 0.0;
    int         sp_n_segments = 64;     // composite Simpson segments

    // Lithostatic confining pressure (Pa, Kant sign convention: positive =
    // compressive). Applied as the +p·nu/(1-nu) offset that ElasticConstants
    // forms and ScanSigmaT adds to the criterion stress (Shared.H), NOT as an
    // elastic-solver BC. Default 0 short-circuits to byte-identical.
    Set::Scalar confining_p = 0.0;

    // ⚠⚠ P2-2d deliverable 7: `GrimselTable` routes K_Ic to the MEASURED
    // at-temperature ratio for the rock the experiments were run on
    // (Grimsel.H's T3_KIc, ten rows). `NasseriTable` is WESTERLY granite
    // (SpCriterion.H:86) applied to Grimsel rock, and stays the DEFAULT so the
    // shipped tree is byte-identical.
    // ⚠ K_Ic0 comes from the LAW (`Grimsel::K_Ic0`), never a new scalar.
    // ⚠⚠ The Grimsel ratio FLOORS at 600 degC (0.60 repeated at 600 and 650),
    // so K_Ic >= 0.900 always -- that floor is the cliff `(a5)` maps, and it is
    // why nothing is scored above 873 K.
    enum class KIcLaw { NasseriTable, Linear, GrimselTable };
    KIcLaw      kic_law         = KIcLaw::NasseriTable;
    Set::Scalar kic_scale       = 1.0;  // multiplicative scale on Nasseri table
    Set::Scalar kic_K_Ic0       = 0.0;  // linear law: K_Ic(T_ref)
    Set::Scalar kic_T_ref       = 0.0;  // linear law: reference temperature (K)
    Set::Scalar kic_T_melt      = 0.0;  // linear law: melting span (K)

    // P0b deliverable 3.
    Set::Scalar sp_flake_c_f = 0.4;
    // Revision 7 deliverable 1: the EMPIRICAL flake-thickness bound, in metres.
    // 0 = inert. ⚠ A STOPGAP for a gradient-blind criterion; see ReconFromX.
    Set::Scalar sp_flake_h_max = 0.0;
    // The share of P_robin the side faces carry (the ledger; registered only
    // under robin_form = surface).
    Set::Scalar surf_side_share = 0.0;
    bool side_power_header_written = false;
    long side_power_last = -1;
    // P0b (c''): cumulative enthalpy carried out of the domain by flakes, the
    // sink the remap creates. ⚠ The existing ledger_E_H sums H over ALL cells
    // INCLUDING void ones (whose H is deliberately frozen rather than reset),
    // so removal has never been a sink in its eyes. This is what closes it.
    Set::Scalar ledger_E_flake = 0.0;
    // P2-24: the face list's vapour rule and the melt exclusions (FaceRemoval.H)
    Set::Scalar ledger_E_vapor = 0.0;   // cumulative H V of the vaporised cells [J]
    Set::Scalar fr_n_vapor     = 0.0;   // per step: surface (listed) cells vaporised
    Set::Scalar fr_n_vapor_int = 0.0;   // per step: INTERIOR solid cells past H_vap_hi (not voided)
    Set::Scalar fr_n_molten    = 0.0;   // per step: lines whose surface cell is fully liquid (not scanned)
    Set::Scalar fr_n_lid       = 0.0;   // per step: lines with a fully liquid cell in their window (not scanned)
    // P0b (c''): the SOLID-ONLY enthalpy sum and its residual.
    // ⚠⚠ THE EXISTING LEDGER CANNOT SEE REMOVAL AT ALL. `ledger_E_H` sums H
    // over every cell including void ones, whose H is deliberately frozen
    // rather than reset (resetting to ambient produced a checkerboard in the
    // temperature plotfiles). So a voided cell's enthalpy stays in E_H forever
    // and material leaving the domain is invisible to it -- measured here:
    // dE_H = 1363.9 J against E_src = 2264.3 J and E_flake = 2193.0 J, which
    // cannot be reconciled in any combination. Worse, on a VOIDING step the
    // remap's index shift slides the profile down one cell while the voided
    // cell keeps its own H, so all-cell E_H gains a spurious cell's worth of
    // enthalpy per void. Summing solid cells only removes both effects, and
    // then the flake energy IS the sink exactly.
    Set::Scalar ledger_E_H_solid   = 0.0;
    Set::Scalar ledger_err_solid   = 0.0;
    Set::Scalar ledger_E_H0_sol_s  = 0.0;
    Set::Scalar ledger_E_H0_sol_c  = 0.0;


    // Step 15b: Weibull per-cell flaw-length distribution (sp_weibull branch
    // only). Default-off; `weibull.enabled = 1` requires sp_weibull (every
    // deck has the microstructure: microstructure.nphases is required).
    // `flaw_a_mf` is cell-centred: one
    // Weibull sample per cell, characteristic length keyed on
    // is_grain_boundary_mf. The long plan §15b says "face-centred", but
    // ALAMO's microstructure topology fields are all cell-centred and Sp is
    // evaluated per surface cell, so a cell-centred flaw field is the
    // faithful translation of "per grain-boundary face".
    bool        weibull_enabled = false;
    Set::Scalar weibull_a0_gb   = 0.0;  // Weibull scale, grain-boundary cells (m)
    Set::Scalar weibull_a0_ig   = 0.0;  // Weibull scale, intragranular cells (m)
    Set::Scalar weibull_m       = 15.0; // Weibull shape (modulus)
    // The Weibull SIZE CORRECTION's reference volume: the per-cell flaw scale is
    // multiplied by (V_cell/V_0)^(1/m), so a smaller cell draws a smaller flaw
    // and a finer mesh does not break earlier just by drawing more often --
    // mesh objectivity, IF V_0 is held fixed across meshes. ⚠ What the code does
    // in practice (P2-27): every production deck sets V_0 = its own cell volume
    // (validation/shared/check_shared.py rule iv), so the factor is exactly 1 and
    // the correction is OFF -- a0 is then a per-cell flaw scale and the mesh
    // dependence is accepted. Whether to switch it on is a physics decision left
    // to the mesh-check runs (MODEL_CARD §6 row 11, §11). Default (1 mm)^3.
    Set::Scalar weibull_V0      = 1.0e-9; // m^3
    int         weibull_seed    = 0;
    Set::Field<Set::Scalar> flaw_a_mf;


    // A1: every-removal event log `<plot_file>_removal_events.csv`
    // (spall.removal_events_csv; sp_weibull + spall.enabled, or the face
    // list's evaluate-only mode). One row per fired surface cell (regime 1),
    // per vaporised cell (2) and per pad-cleared column (4); the columns are
    // FrEventsColumns (Shared.H): time, col_i, col_j, regime, k_top, T_top,
    // a_f, Sp_top, h_scan, h_applied, n_voided, then the reconstruction and
    // line columns. In evaluate-only mode a row is a firing with n_voided 0.
    bool sp_emit_removal_events_csv = false;
    bool sp_removal_events_csv_header_written = false;

    // ------------------------------------------------------------------
    // P0: a real surface temperature for the firing criterion.
    //
    // For a steady ablation front T(z) = T0 + (T_s - T0) e^{-z/l}, l = kappa/v,
    // the cell average and the true surface are related by
    //     T_s = T0 + (Tbar - T0) * (dz/l) / (1 - e^{-dz/l}),
    // so the criterion's "surface" is mesh- AND rate-dependent: at Tbar = 823 K
    // the true surface is 1215 K at dz = 2 mm and 1003 K at 1 mm. D2d's
    // measured +43 % ROP from 2 mm to 1 mm is that ratio.
    //
    // THE FIX (the reconstruction, the face list's R1). For an exponential profile
    // the top two cell AVERAGES have the exact ratio e^{dz/l}, so l follows
    // from the profile alone -- no recession rate, no flame coupling, no
    // one-step lag, no persistent state. Everything below is a function of
    //     x = dz / l_eff = log((T_top - T0) / (T_top_minus_1 - T0)).
    //
    // ⚠ T0 IS THE LINE'S START CELL'S INITIAL TEMPERATURE (T_init_mf, the
    // stress's own reference, with the stress's fallback phases[0].T_ref when
    // the field is absent or non-finite; FaceRemoval.H). Until P2-31 the
    // reconstruction read phases[0].T_ref, which a deck may type away from its
    // initial condition (unit/gradient_blind: 293.15 K against an 820 K start).
    // It is NOT a new constant and it is NOT tunable here.


    // The ONE place the reconstruction is written: RawSurfaceX, then
    // ReconFromX. Both callers (the face list's line criterion and its
    // surface temperature for a fresh cell, FaceRemoval.H) go through them, so
    // they cannot drift apart. Pure functions of their arguments --
    // no state, no fields, no MPI -- which is also what makes them unit-testable
    // and decomposition-independent.
    //
    //   x       = log((T_top - T0)/(T_1 - T0))          = dz / ell_eff
    //   T_surf  = T0 + (T_top - T0) * x/(1 - e^{-x})    -> T_top as x -> 0
    //   h_f     = c_f * ell_eff (NO dz clamp), then the two caps
    //   T_bar   = T0 + (T_surf - T0) * (1 - e^{-u})/u,  u = h_f/ell_eff
    //
    // Guards, in order: a non-decaying or non-finite pair degrades to the cell
    // average (flat = true, factor exactly 1.0 -- no correction exactly where
    // none is available); x is capped at 4 (ell_eff >= dz/4) so a near-tie in
    // the two averages cannot produce an unbounded surface; and T_surf is
    // capped at the configuration's gas temperature when one exists.
    // The raw two-cell estimate x = dz/ell_eff. Returns 0 when the pair carries
    // no usable decaying gradient ("flat"), which the CALLER resolves -- the
    // face list then takes T_s from the face heat balance (the face Newton)
    // and the flake depth from the line's own e-fold (FaceRemoval.H), and
    // does not call ReconFromX.
    static Set::Scalar RawSurfaceX(Set::Scalar T_top, Set::Scalar T_1,
                                   Set::Scalar T0)
    {
        const Set::Scalar dT0 = T_top - T0;
        const Set::Scalar dT1 = T_1  - T0;
        if (!(std::isfinite(dT0) && std::isfinite(dT1) &&
              dT0 > 0.0 && dT1 > 0.0 && dT0 > dT1))
            return 0.0;
        const Set::Scalar x = std::log(dT0 / dT1);
        return (std::isfinite(x) && x > 0.0) ? x : 0.0;
    }

    // Everything from a given x = dz/ell_eff. Both callers pass RawSurfaceX's
    // value and call only when it is > 0; x <= 0 gives the flat result
    // (T_surf = T_top, h_f = dz), which no caller reads.
    static SurfaceRecon ReconFromX(Set::Scalar T_top, Set::Scalar T0,
                                   Set::Scalar dz, Set::Scalar x,
                                   Set::Scalar T_gas_cap, Set::Scalar c_f,
                                   Set::Scalar h_f_max,
                                   Set::Scalar h_f_bound = 0.0)
    {
        SurfaceRecon r;
        const Set::Scalar dT0 = T_top - T0;
        if (!(x > 0.0)) { x = 0.0; r.flat = true; }
        if (x > 4.0)    { x = 4.0; r.ell_clamped = true; }

        const Set::Scalar fac = (x > 0.0) ? x / (1.0 - std::exp(-x)) : 1.0;
        r.T_surf  = T0 + dT0 * fac;
        r.ell_eff = (x > 0.0) ? dz / x
                              : std::numeric_limits<Set::Scalar>::infinity();
        if (T_gas_cap > 0.0 && r.T_surf > T_gas_cap)
        { r.T_surf = T_gas_cap; r.T_clamped = true; }

        // P0b deliverable 3: h_f = c_f * ell_eff, NO dz clamp, so
        //   T_bar - T0 = (T_surf - T0) * (1 - e^{-c_f}) / c_f
        // is a MESH-INDEPENDENT multiple of the surface rise (P6). A
        // min(ell_eff, dz) clamp would make the multiple mesh-dependent and
        // collapse it to the cell average wherever ell_eff exceeds dz.
        r.h_f   = (x > 0.0) ? c_f * r.ell_eff : dz;
        // ⚠⚠ ell_eff IS CLAMPED FROM BELOW (>= dz/4) BUT NOT FROM ABOVE, so a
        // nearly-flat column -- x small but not zero -- gives an UNBOUNDED
        // h_f. This is not hypothetical: the Meier rerun aborted at t = 129 s
        // with a column trying to void 287 cells (574 mm, more than the whole
        // block) off a two-cell ratio of 1.0014.
        // The cap is spall.h_spall_max, the key that already defines how deep
        // this model looks for a flake at all. ⚠ THAT IS NOT THE dz CLAMP THE
        // PACKET FORBIDS: dz is a MESH quantity and clamping there is what
        // made T_bar mesh-dependent in P0; h_spall_max is PHYSICAL and
        // mesh-independent, so the cap cannot reintroduce mesh dependence.
        // It binds on nearly-flat profiles only, and a bind is reported
        // (rc_flag bit 8 in the event log).
        if (h_f_max > 0.0 && r.h_f > h_f_max)
        { r.h_f = h_f_max; r.h_f_capped = true; }
        // ⚠⚠ REVISION 7 DELIVERABLE 1 -- THE FLAKE BOUND, AND IT IS A STOPGAP.
        //
        // `spall.flake_h_max` bounds the FLAKE THICKNESS and nothing else.
        // ⚠ It is NOT spall.h_spall_max. That key is a SEARCH DEPTH -- how far
        // down this model looks for a flake at all -- and revision 5 measured
        // it acting as the thickness bound on 13.5 % of Meier's firings, which
        // is a search parameter setting a physical answer. The two are now
        // separate quantities with separate keys.
        //
        // WHY IT IS NEEDED. On a side-heated Meier column the two-cell ratio
        // goes nearly flat -- 1.3 where a front receding at 5.9 m/h demands
        // 117 -- so ell_eff reads 7.5 mm against a physical kappa/v of
        // 0.42 mm, 18x too long, and h_f = c_f*ell_eff follows it. A flat
        // ratio means THE GRADIENT IS UNRESOLVED ("deeper than two cells can
        // resolve"), not that the flake is thick.
        //
        // ⚠⚠ WHAT IT IS NOT. It does not fix the physics. The firing criterion
        // is GRADIENT-BLIND by construction -- a Tada edge crack over a
        // 4-20 um flaw sees an effectively uniform stress, so a flat slab at
        // the firing temperature fires exactly as a steep skin does. Nothing
        // in the K_I test says a preheated slab should spall less. The
        // gradient-aware part of the model is the buckling thickness
        // condition (spall.buckle_span, BuckleRatio in FaceRemoval.H; off at
        // its default 0), which withdraws a fired line whose flake is too
        // thick to buckle. This key buys the Meier comparison; it does not
        // earn it.
        //
        // ⚠ It is applied to h_f ONLY. x and T_surf stay raw, because on a
        // genuinely flat profile the correct surface temperature IS the
        // top-cell average and the raw reconstruction returns exactly that
        // (fac -> 1 as x -> 0). Holding a steeper x would fire the column
        // hundreds of degrees early -- that is the h500 clamp failure.
        //
        // EMPIRICAL, and declared so: Wilkinson's measured 0.1-1 mm flake band
        // and Rauenzahn Fig 4.12. Default 0 = inert.
        if (h_f_bound > 0.0 && r.h_f > h_f_bound)
        { r.h_f = h_f_bound; r.h_f_bounded = true; }
        const Set::Scalar u = (x > 0.0) ? r.h_f / r.ell_eff : 0.0;
        r.T_bar = T0 + (r.T_surf - T0)
                  * ((u > 0.0) ? (1.0 - std::exp(-u)) / u : 1.0);
        return r;
    }

    // Upper bound on the gas temperature this configuration can present to the
    // rock, used only as a ceiling on the reconstruction. Time-independent and
    // read from Parse'd values, so it introduces NO dependence on
    // BuildFlameColumns -- which runs AFTER the criterion (Advance: RemovalStep,
    // then AdvanceMicrostructure). Returns 0 when no gas temperature is defined, which
    // means "no ceiling"; the x <= 4 cap still bounds the result.
    Set::Scalar SurfaceGasTempCap() const
    {
        if (!surface_patch.enabled) return 0.0;
        Set::Scalar cap = 0.0;
        if (surface_patch.T_flame      > cap) cap = surface_patch.T_flame;
        if (surface_patch.jet_T_nozzle > cap) cap = surface_patch.jet_T_nozzle;
        return cap;
    }

    // splitmix64 finalizer — deterministic 64-bit hash used to seed the
    // per-cell Weibull RNG from the GLOBAL cell index so flaw_a_mf is
    // bit-identical regardless of MPI rank count / BoxArray partition.
    static inline std::uint64_t SplitMix64(std::uint64_t z)
    {
        z += 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    static void ParseSpallationModel(MMWSpalling& value, IO::ParmParse& pp)
    {
        std::string model_str;
        // P2-27: the default `none` loads no criterion (spallation.sp.a0 is not
        // required: the early return below).
        pp.query_default("spallation.model", model_str, std::string("none"));
        if      (model_str == "none")       value.spallation_model = SpallationModel::None;
        else if (model_str == "sp_weibull") value.spallation_model = SpallationModel::SpWeibull;
        else if (model_str == "damage_law")
            Util::Abort(INFO, "spallation.model = damage_law was retired in P2-27 (the Drucker-Prager "
                        "damage lineage is deleted; LEFM is the model): use sp_weibull, or none for a "
                        "deck without a criterion");
        else Util::Abort(INFO,
            "spallation.model must be 'none' or 'sp_weibull', got: ", model_str);

        // Step 18: confining pressure (Pa, default 0). Parsed regardless
        // of spallation_model so a criterion-free input that sets the key
        // runs; the offset is only applied in the sp_weibull criterion
        // stress (ElasticConstants -> ScanSigmaT, Shared.H).
        pp.query_default("confining.p", value.confining_p, 0.0);
        if (!std::isfinite(value.confining_p))
            Util::Abort(INFO, "confining.p must be finite, got ",
                        value.confining_p);

        // P2-15: the buckling span, read on EVERY deck, BEFORE the sp_weibull early return:
        // ALAMO does not abort on an unused top-level key (standing rule 1), so a key read
        // only on one branch would be silently inert on the other. Its cross-key requirement
        // (the face list) is checked in CheckInterlocks.
        pp.query_default("spall.buckle_span", value.sp_buckle_span, 0.0);
        if (!(std::isfinite(value.sp_buckle_span) && value.sp_buckle_span >= 0.0))
            Util::Abort(INFO, "spall.buckle_span must be finite and >= 0 m, got: ",
                        value.sp_buckle_span);
        // P2-23: the blocked-line record (default off; a LOG, read by nothing).
        // Checked against the span in CheckInterlocks, beside the span's own face-list check.
        {
            int bl = 0;
            pp.query_default("spall.buckle_log_csv", bl, 0);
            if (bl != 0 && bl != 1)
                Util::Abort(INFO, "spall.buckle_log_csv must be 0 or 1, got: ", bl);
            value.sp_buckle_log_csv = (bl == 1);
        }
        if (value.spallation_model != SpallationModel::SpWeibull) return;

        pp.query_required("spallation.sp.a0", value.sp_a0);
        if (!(std::isfinite(value.sp_a0) && value.sp_a0 > 0.0))
            Util::Abort(INFO, "spallation.sp.a0 must be finite and > 0, got ", value.sp_a0);

        pp.query_default("spallation.sp.n_segments", value.sp_n_segments, 64);
        // P2-2 change 4. Default 1.0 == identity, so an input that does not
        // name this key is byte-identical. 0 < R <= 1: it RELIEVES.
        pp.query_default("spallation.sp.stress_relief", value.sp_relief, 1.0);
        if (!(std::isfinite(value.sp_relief) && value.sp_relief > 0.0
              && value.sp_relief <= 1.0))
            Util::Abort(INFO, "spallation.sp.stress_relief must be finite and in "
                              "(0, 1] -- it RELIEVES the thermoelastic stress, so "
                              "a value above 1 would amplify it. Got: ",
                        value.sp_relief);
        if (value.sp_n_segments < 2)
            Util::Abort(INFO, "spallation.sp.n_segments must be >= 2, got ",
                        value.sp_n_segments);
        if (value.sp_n_segments % 2 != 0) ++value.sp_n_segments;

        std::string law_str;
        pp.query_default("spallation.sp.kic.law", law_str, std::string("nasseri_table"));
        if      (law_str == "nasseri_table") value.kic_law = KIcLaw::NasseriTable;
        else if (law_str == "linear")        value.kic_law = KIcLaw::Linear;
        else if (law_str == "grimsel_table") value.kic_law = KIcLaw::GrimselTable;
        else Util::Abort(INFO,
            "spallation.sp.kic.law must be 'nasseri_table', 'linear' or "
            "'grimsel_table', got: ", law_str);
        // ⚠⚠ `grimsel_table` reads K_Ic0 and the at-temperature ratio off the
        // PHASE-0 LAW OBJECT, which P2-2b parses unconditionally -- so this
        // works whether or not the thermal law is on, and a leg can select the
        // toughness route alone.
        // ⚠ ITS "needs a phase" CHECK CANNOT LIVE HERE. `ParseSpallationModel`
        // runs before Parse's microstructure block resizes `value.phases`, so a
        // `phases.empty()` test at this point is ALWAYS true and aborts every
        // legal leg. The check sits with the other cross-key aborts, in
        // CheckInterlocks, instead.

        // ⚠⚠ `grimsel_table` TAKES `scale` AND NOT `K_Ic0`. The else-branch
        // below `query_required`s `spallation.sp.kic.K_Ic0`, which is right for
        // `linear` (it has no table to read a level off) and WRONG for
        // `grimsel_table`: its K_Ic0 is the LAW's, and requiring a scalar here
        // would reintroduce exactly the hand-typed per-phase constant the
        // packet forbids -- and would let a leg set a level that silently
        // disagrees with the law it is supposed to be reading.
        if (value.kic_law == KIcLaw::NasseriTable
            || value.kic_law == KIcLaw::GrimselTable)
        {
            pp.query_default("spallation.sp.kic.scale", value.kic_scale, 1.0);
            if (!(std::isfinite(value.kic_scale) && value.kic_scale > 0.0))
                Util::Abort(INFO, "spallation.sp.kic.scale must be finite and > 0, got ",
                            value.kic_scale);
        }
        else
        {
            pp.query_required("spallation.sp.kic.K_Ic0",  value.kic_K_Ic0);
            pp.query_required("spallation.sp.kic.T_ref", value.kic_T_ref);
            pp.query_required("spallation.sp.kic.T_melt", value.kic_T_melt);
            if (!(std::isfinite(value.kic_K_Ic0) && value.kic_K_Ic0 > 0.0))
                Util::Abort(INFO, "spallation.sp.kic.K_Ic0 must be finite and > 0, got ",
                            value.kic_K_Ic0);
            if (!(std::isfinite(value.kic_T_melt) && value.kic_T_melt > 0.0))
                Util::Abort(INFO, "spallation.sp.kic.T_melt must be finite and > 0, got ",
                            value.kic_T_melt);
        }

        // Step 15b: Weibull per-cell flaw-length distribution. Default-off.
        // Only parsed under sp_weibull (we already returned for none);
        // every deck has the microstructure (microstructure.nphases is
        // required). With
        // weibull.enabled = 0 the Sp evaluator uses the scalar
        // spallation.sp.a0 exactly as Step 12, which keeps the §12c
        // deterministic-a0 verification path intact.
        int weibull_on = 0;
        pp.query_default("weibull.enabled", weibull_on, 0);
        value.weibull_enabled = (weibull_on != 0);
        // ⚠⚠ P1g-3 deliverable 4: `weibull.m` IS A MATERIAL CONSTANT, NOT A
        // WEIBULL-MODE KEY, and it is parsed OUTSIDE the mode guard for that
        // reason: inside it, a run with `weibull.enabled = 0` would never
        // execute the parse and would silently take the MEMBER default 15.0
        // whatever the deck says. Its one reader today is InitializeFlaws,
        // which runs only under weibull.enabled.
        // ⚠ IDENTITY-SAFE BY CONSTRUCTION: the query_default's fallback (15.0)
        // equals the member default, so any run that does not set the
        // key is byte-identical either way. (b1) is the witness, not the claim.
        pp.query_default("weibull.m", value.weibull_m, 15.0);
        if (!(std::isfinite(value.weibull_m) && value.weibull_m > 0.0))
            Util::Abort(INFO, "weibull.m must be finite and > 0, got ",
                        value.weibull_m);
        if (value.weibull_enabled)
        {
            pp.query_required("weibull.a0_gb", value.weibull_a0_gb);
            if (!(std::isfinite(value.weibull_a0_gb) && value.weibull_a0_gb > 0.0))
                Util::Abort(INFO, "weibull.a0_gb must be finite and > 0, got ",
                            value.weibull_a0_gb);

            // Default intragranular scale is a0_gb / 5 (Rossi 2018
            // rho_B/rho_I ~ 5.8 partition; revised-model-approach.md §3).
            value.weibull_a0_ig = value.weibull_a0_gb / 5.0;
            pp.query_default("weibull.a0_ig", value.weibull_a0_ig,
                             value.weibull_a0_ig);
            if (!(std::isfinite(value.weibull_a0_ig) && value.weibull_a0_ig > 0.0))
                Util::Abort(INFO, "weibull.a0_ig must be finite and > 0, got ",
                            value.weibull_a0_ig);



            // Mesh-objectivity reference volume (weakest-link size effect): the
            // per-cell flaw scale is weighted by (V_cell/V_0)^(1/m). Default
            // (1 mm)^3 = 1e-9 m^3 is the design resolution, so the weighting is
            // exactly 1 there and existing 1 mm runs are byte-identical; set
            // this to the cell volume at which a0_gb/a0_ig were calibrated.
            pp.query_default("weibull.V0", value.weibull_V0, 1.0e-9);
            if (!(std::isfinite(value.weibull_V0) && value.weibull_V0 > 0.0))
                Util::Abort(INFO, "weibull.V0 must be finite and > 0, got ",
                            value.weibull_V0);

            pp.query_default("weibull.seed", value.weibull_seed, 0);

        }
    }
