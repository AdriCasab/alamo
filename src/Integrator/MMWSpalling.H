#ifndef INTEGRATOR_MMWSPALLING_H
#define INTEGRATOR_MMWSPALLING_H

// Step 4b: + material-driven T-dependent ρ(T), κ(T), c_p(T) and time-varying
//           beam power. When a `material` block is supplied (Zhang granite or
//           any future subclass of Numeric::Material::Material), Advance uses
//           a conservative flux-divergence form `∂H/∂t = ∇·(κ(T)∇T) + Q − q_loss`
//           and inverts H↔T through a one-time tabulated map that bakes in
//           the smoothed latent-heat plateaux. Without a `material` block the
//           integrator falls back to step 4's constant-α path.
// Step 4: + surface losses (radiation + convection).
// Step 3: + Beer-Lambert + Gaussian MMW volumetric source.
// Step 2: enthalpy formulation of the heat equation.
//
// Inherits HeatConduction (so we reuse its temp_mf/temp_old_mf with BC + ghost
// cells as our T field) and Mechanics (for the elastic solve). The base
// HeatConduction::Advance evolves T directly via ∂T/∂t = α ∇²T; we OVERRIDE
// Advance to instead evolve enthalpy H via the conservative
//
//   ∂H/∂t = ∇·(κ ∇T) + q
//
// (algebraically equivalent to main.tex eq. 14 ∂H/∂t = ∇·[(κ/(ρ c_p)) ∇H] when
// the coefficient is interpreted as κ · dT/dH, which vanishes on phase
// plateaux — see CLAUDE.md). H is the conserved variable; T is recovered by
// inverting H each step. T is held at T_m or T_v across phase plateaux, so
// adjacent cells in a plateau have ∇T = 0 and latent heat is correctly
// absorbed without conduction "leaking" through the front.
//
// See in-main-tex-you-will-quizzical-treasure.md.

#include "Model/Solid/Affine/Isotropic.H"
#include "Integrator/Mechanics.H"
#include "Integrator/HeatConduction.H"
#include "Numeric/Stencil.H"
#include "Numeric/MMWBeam.H"
#include "Numeric/DruckerPrager.H"
#include "Numeric/CohesiveZone.H"
#include "Numeric/SpCriterion.H"
#include "Numeric/Material/Material.H"
#include "Numeric/Material/Zhang.H"
#include "Numeric/Material/Constant.H"
#include "Numeric/Material/Table.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "BC/Constant.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>

#include <eigen3/Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <ios>
#include <limits>
#include <random>

namespace Integrator
{
class MMWSpalling :
    virtual public HeatConduction,
    virtual public Mechanics<Model::Solid::Affine::Isotropic>
{
public:
    static constexpr const char* name = "mmwspalling";

    MMWSpalling() :
        HeatConduction(3),
        Mechanics<Model::Solid::Affine::Isotropic>()
    {}

    ~MMWSpalling() = default;

    MMWSpalling(IO::ParmParse& pp) : MMWSpalling() { Parse(*this, pp); }

    static void Parse(MMWSpalling& value, IO::ParmParse& pp)
    {
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

        // Cache the diffusivity for use in Advance (HeatConduction's alpha is
        // private; cheaper than another queryclass).
        pp.query_required("hc.heat.alpha", value.therm_alpha, Unit::ThermalDiffusivity());

        // Phase-change parameters (single material in step 2; per-mineral in
        // step 5). Defaults disable phase change.
        pp.query_default("phase.T_ref", value.T_ref, 0.0);
        pp.query_default("phase.T_m",   value.T_m,   1.0e30);
        pp.query_default("phase.T_v",   value.T_v,   1.0e30);
        pp.query_default("phase.lh_m",  value.lh_m,  0.0);  // L_m / c_p [K]
        pp.query_default("phase.lh_v",  value.lh_v,  0.0);  // L_v / c_p [K]

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

        // Optional material model. When omitted, the constant-α path of
        // step 2 is used (driven by hc.heat.alpha). When present, Advance
        // takes the variable-property flux-divergence path and the
        // enthalpy table is built from the material's ρ(T)·c_p(T).
        if (pp.contains("material.type"))
            pp.select_default<Numeric::Material::Zhang,
                              Numeric::Material::Constant>("material", value.material);

        // Step 6d: optional prescribed-temperature surface patch (Hu LRST
        // conduction validation). When `surface_patch.enabled = 1`, a
        // circular patch on a single domain face injects a Dirichlet-like
        // flux  kappa * (T_patch(t) - T_cell) / dz^2  on in-patch top-row
        // cells. Default disabled keeps every prior MMWSpalling test
        // unchanged. Currently only the explicit material path consumes
        // this term; the implicit path hardcodes Neumann domain BCs.
        pp.queryclass<SurfacePatch>("surface_patch", value.surface_patch);

        // Step 8: optional Drucker-Prager damage evolution. Disabled by
        // default so existing passive `damage.ic.*` inputs keep their Step 6
        // conductivity-only behavior.
        ParseDamageSettings(value, pp);

        // Step 9: optional grain-boundary cohesive-zone law. Disabled by
        // default; when enabled, registers gb_* diagnostic fields and runs
        // an irreversible bilinear traction-separation update keyed on
        // `is_grain_boundary`. Step 9 only computes law/history/diagnostics;
        // full feedback into the elastic operator is deferred.
        ParseCohesiveSettings(value, pp);

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

        // Step 12: spallation-model toggle. `spallation.model = damage_law`
        // (default) preserves every Step 8-11b behavior; `sp_weibull`
        // activates the LEFM Sp = K_I/K_Ic(T) onset criterion (deterministic
        // a_0 in this step; Weibull randomization arrives in Step 15) and
        // bypasses CZM. Mixing branches is a hard parser error: parsing the
        // toggle here, after the damage/cohesive blocks, lets us reject
        // `sp_weibull` + `damage.enabled = 1` and `sp_weibull` +
        // `cohesive.enabled = 1` immediately.
        ParseSpallationModel(value, pp);

        // Step 19: optional mechanical bit feed (drawworks). Default-off.
        // Parsed after spall/vapor so it can require one of them.
        ParseBitSettings(value, pp);

        // Solver choice for the diffusion term: "explicit" forward Euler
        // (default; matches steps 2-4) or "implicit" backward Euler (built
        // on amrex::MLABecLaplacian + MLMG). The source Q and surface loss
        // q_loss stay explicit in both — they are not stiff once diffusion
        // is. Currently only the material path supports "implicit"; the
        // constant-α path is always explicit.
        std::string solver_str;
        pp.query_default("solver", solver_str, std::string("explicit"));
        if      (solver_str == "explicit") value.solver_kind = SolverKind::Explicit;
        else if (solver_str == "implicit") value.solver_kind = SolverKind::Implicit;
        else Util::Abort(INFO, "solver must be 'explicit' or 'implicit', got: ", solver_str);
        pp.query_default("implicit.tol_rel",  value.implicit_tol_rel,  1.0e-10);
        pp.query_default("implicit.tol_abs",  value.implicit_tol_abs,  0.0);
        pp.query_default("implicit.max_iter", value.implicit_max_iter, 100);
        pp.query_default("implicit.verbose",  value.implicit_verbose,  0);

        // Thermal-expansion eigenstrain coefficient per model component.
        pp_queryarr("alpha", value.alpha);

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
                pp.query_default (pfx + ".name",     value.phases[i].name,
                                  std::string("phase") + std::to_string(i));
                pp.query_required(pfx + ".fraction", value.phases[i].fraction);
                pp.query_required(pfx + ".kappa",    value.phases[i].kappa);
                pp.query_required(pfx + ".beta",     value.phases[i].beta);
                pp.query_required(pfx + ".E",        value.phases[i].E);
                pp.query_required(pfx + ".mu",       value.phases[i].mu);
                pp.query_required(pfx + ".rho",      value.phases[i].rho);
                pp.query_required(pfx + ".Cp",       value.phases[i].Cp);
                pp.query_required(pfx + ".alpha_attenuation",
                                  value.phases[i].alpha_atten);
                pp.query_required(pfx + ".A_D",      value.phases[i].A_D);
                Util::Assert(INFO, TEST(value.phases[i].fraction >= 0.0));

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

            // Step 6/8: damage field, damage-modified conductivity, GB
            // conductance hook, and bulk effective-conductivity diagnostic.
            // With `damage.enabled = 0`, D remains passive and is either zero
            // everywhere or preset via `damage.ic.*`. With
            // `damage.enabled = 1`, Step 8 evolves D irreversibly after the
            // mechanics solve.
            pp.query_default("conductivity.kappa_damage_alpha",
                             value.kappa_damage_alpha, 0.0);
            pp.query_default("conductivity.h_gb0",
                             value.h_gb_conductance,   0.0);

            // Optional damage IC. Omitted ⇒ D ≡ 0 everywhere.
            if (pp.contains("damage.ic.type"))
                pp.select<IC::Constant, IC::Expression>("damage.ic",
                                                        value.ic_damage,
                                                        value.geom);

            // 1 ghost on D for face-centred damage averaging in the heat
            // kernel (matches phase_mf's ghost layer).
            value.RegisterNewFab(value.D_mf,         nullptr, 1, 1, "D",          true);
            value.RegisterNewFab(value.kappa_eff_mf, nullptr, 1, 0, "kappa_eff",  true);
            value.RegisterNewFab(value.k_eff_mf,     nullptr, 1, 0, "k_eff",      true);
            if (value.damage_enabled)
            {
                value.RegisterNewFab(value.damage_F_mf,    nullptr, 1, 0, "damage_F_DP", true);
                value.RegisterNewFab(value.damage_Psi_mf,  nullptr, 1, 0, "damage_Psi",  true);
                value.RegisterNewFab(value.damage_rate_mf, nullptr, 1, 0, "damage_rate", true);
            }
            if (value.cohesive_enabled)
            {
                value.RegisterNewFab(value.gb_delta_mf,           nullptr, 1, 0, "gb_delta",           true);
                value.RegisterNewFab(value.gb_delta_max_mf,       nullptr, 1, 0, "gb_delta_max",       true);
                value.RegisterNewFab(value.gb_traction_mf,        nullptr, 1, 0, "gb_traction",        true);
                value.RegisterNewFab(value.gb_cohesive_damage_mf, nullptr, 1, 0, "gb_cohesive_damage", true);
            }
            if (value.removal_enabled())
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
            // indicator (damage_law back-compat); `regime_field_mf` is
            // the canonical per-cell diagnostic for the sp_weibull
            // per-cell v_n logic — distinct so plotfile consumers can
            // tell the two diagnostics apart.
            if (value.spallation_model == SpallationModel::SpWeibull &&
                (value.spall_enabled || value.vapor_enabled))
            {
                value.RegisterNewFab(value.regime_field_mf, nullptr, 1, 0,
                                     "regime_field", true);
            }
            // Step 12: surface Sp diagnostic for the sp_weibull branch only.
            // Written into the top z-cell of each (i, j) column after every
            // static mechanics solve; zero on cells below the surface.
            if (value.spallation_model == SpallationModel::SpWeibull)
            {
                value.RegisterNewFab(value.Sp_field_mf, nullptr, 1, 0, "Sp_field", true);
                // Step 15b: per-cell Weibull flaw length, sampled at
                // Initialize (cell-centred, 0 ghost — read only at the
                // owning cell in UpdateSpAfterMechanics).
                if (value.weibull_enabled)
                {
                    value.RegisterNewFab(value.flaw_a_mf, nullptr, 1, 0, "flaw_a", true);
                    // Step 15c: connected-component cluster ID for the
                    // Sp >= 1 surface zone. 0 = not firing; positive int =
                    // cluster ID at the top-z plane; 0 everywhere else.
                    value.RegisterNewFab(value.Sp_cluster_id_mf, nullptr, 1, 0,
                                         "Sp_cluster_id", true);
                }
                // Step 16: per-column h_spall diagnostic, written at the
                // top-z cell of each firing column under sp_weibull. Gated
                // on spall.enabled so damage_law / Kant onset runs (no
                // removal) don't carry the field.
                if (value.spall_enabled)
                {
                    value.RegisterNewFab(value.h_spall_field_mf, nullptr, 1, 0,
                                         "h_spall_field", true);
                    // Step 16b: time-integrated 3-D crack-event count
                    // field. Edge-triggered (first crossing only); NOT
                    // reset per step. See rossi-validation-diagnostic-
                    // design.md §2.
                    value.RegisterNewFab(value.crack_event_count_mf, nullptr,
                                         1, 0, "crack_event_count", true);
                }
            }
        }

        if (value.damage_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: damage.enabled requires microstructure.nphases");

        if (value.cohesive_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: cohesive.enabled requires microstructure.nphases");

        if (value.spall_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spall.enabled requires microstructure.nphases");

        if (value.vapor_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: vapor.enabled requires microstructure.nphases");

        if (value.spallation_model == SpallationModel::SpWeibull &&
            !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spallation.model = sp_weibull requires "
                              "microstructure.nphases (Step 15 Weibull samples will key "
                              "on grain-boundary faces)");

        // S1: surface-resolution study switches (default-off; microstructure
        // path, single level).
        {
            int fm = 0, le = 0;
            pp.query_default("surface.follow_mask",   fm, 0);
            pp.query_default("energy_ledger.enabled", le, 0);
            value.surface_follow_mask = (fm != 0);
            value.ledger_enabled      = (le != 0);
            if (value.surface_follow_mask)
            {
                if (!value.removal_enabled())
                    Util::Abort(INFO, "surface.follow_mask = 1 requires spall.enabled "
                                      "or vapor.enabled");
                if (value.surface_patch.enabled &&
                    value.surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature)
                    Util::Abort(INFO, "surface.follow_mask = 1 is not supported with "
                                      "surface_patch.mode = prescribed_T");
            }
            if ((value.surface_follow_mask || value.ledger_enabled) &&
                !value.microstructure_enabled)
                Util::Abort(INFO, "surface.follow_mask / energy_ledger.enabled require "
                                  "the microstructure path (microstructure.nphases)");
            if ((value.surface_follow_mask || value.ledger_enabled) && value.max_level > 0)
                Util::Abort(INFO, "surface.follow_mask / energy_ledger.enabled are "
                                  "single-level only (amr.max_level = 0)");
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
            // C1: pinned-surface closure prerequisites (the pin comes from the
            // removal pipeline's per-column firing T at the mask top).
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Pinned)
            {
                if (!value.microstructure_enabled)
                    Util::Abort(INFO, "surface_patch.robin_form = pinned is implemented on "
                                      "the microstructure path only");
                if (!value.surface_follow_mask)
                    Util::Abort(INFO, "surface_patch.robin_form = pinned requires "
                                      "surface.follow_mask = 1");
                if (!(value.spall_enabled &&
                      value.spallation_model == SpallationModel::SpWeibull))
                    Util::Abort(INFO, "surface_patch.robin_form = pinned requires spall.enabled = 1 "
                                      "and spallation.model = sp_weibull");
                if (value.max_level > 0)
                    Util::Abort(INFO, "surface_patch.robin_form = pinned is single-level only "
                                      "(amr.max_level = 0)");
                if (value.eps_rad_hi < value.eps_rad)
                    Util::Abort(INFO, "surface_patch.robin_form = pinned requires "
                                      "losses.epsilon_high >= losses.epsilon (got ", value.eps_rad_hi,
                                      " < ", value.eps_rad, "): the face-temperature Newton assumes a "
                                      "monotone residual with a unique root");
                value.RegisterIntegratedVariable(&value.pinned_cols_pinned, "pinned_cols_pinned", false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face,   "pinned_cols_face",   false);
                value.RegisterIntegratedVariable(&value.pinned_q_mean,      "pinned_q_mean",      false);
            }
            // A1: removal-event log validity.
            if (value.sp_emit_removal_events_csv)
            {
                if (!(value.spall_enabled &&
                      value.spallation_model == SpallationModel::SpWeibull))
                    Util::Abort(INFO, "spall.removal_events_csv = 1 requires spall.enabled = 1 "
                                      "and spallation.model = sp_weibull");
                if (value.max_level > 0)
                    Util::Abort(INFO, "spall.removal_events_csv is single-level only "
                                      "(amr.max_level = 0)");
            }
            if (value.ledger_enabled)
            {
                value.RegisterIntegratedVariable(&value.ledger_E_H,     "ledger_E_H",     false);
                value.RegisterIntegratedVariable(&value.ledger_E_src,   "ledger_E_src",   false);
                value.RegisterIntegratedVariable(&value.ledger_err,     "ledger_err",     false);
                value.RegisterIntegratedVariable(&value.ledger_P_beam,  "ledger_P_beam",  false);
                value.RegisterIntegratedVariable(&value.ledger_P_robin, "ledger_P_robin", false);
                value.RegisterIntegratedVariable(&value.ledger_P_loss,  "ledger_P_loss",  false);
                value.RegisterIntegratedVariable(&value.ledger_surface_missing_cols,
                                                 "ledger_surface_missing_cols", false);
            }
            // D2a: pinned branch-reason counters (sum = pinned_cols_face) and the
            // spatial-flame diagnostics. Registered after every pre-D2a variable
            // so existing thermo.dat columns keep their positions.
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Pinned)
            {
                value.RegisterIntegratedVariable(&value.pinned_cols_face_unset,   "pinned_cols_face_unset",   false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face_idle,    "pinned_cols_face_idle",    false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face_qneg,    "pinned_cols_face_qneg",    false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face_minrule, "pinned_cols_face_minrule", false);
            }
            if (value.surface_patch.enabled && value.surface_patch.use_expr)
            {
                if (!value.surface_follow_mask)
                    Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr require "
                                      "surface.follow_mask = 1");
                value.RegisterIntegratedVariable(&value.patch_min_standoff, "patch_min_standoff", false);
                value.RegisterIntegratedVariable(&value.patch_P_robin,      "patch_P_robin",      false);
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
    }

protected:
    void Initialize(int lev) override
    {
        // Step 6: simultaneous `material` + `microstructure` is not yet
        // composed (the heterogeneous path uses constant per-phase ρ, c_p, κ;
        // the material path uses T-dependent polynomial ρ(T), c_p(T), κ(T)).
        // Reject the combination until a composition rule is defined.
        if (material && microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: simultaneous `material` and `microstructure` blocks "
                              "are not supported (Step 6 has no composition rule). Drop one.");

        if (material && !material_table.built()) BuildEnthalpyTable();   // build once
        HeatConduction::Initialize(lev);                      // sets temp_mf, temp_old_mf
        Mechanics<Model::Solid::Affine::Isotropic>::Initialize(lev);
        if (removal_enabled()) InitPhi(lev);
        InitSurfaceMask(lev);
        // InitializeMicrostructure must run before ConvertTtoH on the
        // microstructure path because the H(T) inversion is per-phase.
        if (microstructure_enabled) InitializeMicrostructure(lev, /*announce=*/true,
                                                             /*reset_damage=*/true);
        ConvertTtoH(lev);
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        HeatConduction::TimeStepBegin(a_time, a_step);
        if (microstructure_enabled && max_level > 0)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);
    }

    void TimeStepComplete(Set::Scalar a_time, int a_step) override
    {
        HeatConduction::TimeStepComplete(a_time, a_step);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepComplete(a_time, a_step);
        (void)a_time;
        (void)a_step;
        if (microstructure_enabled && max_level > 0)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::Regrid(lev, time);
        InitSurfaceMask(lev);
        if (microstructure_enabled)
            InitializeMicrostructure(lev, /*announce=*/false,
                                     /*reset_damage=*/!damage_enabled);
        RepairMechanicsAfterRegrid(time);
    }

    void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::UpdateModel(a_step, a_time);

        // Step 8b: prescribed-temperature surface_patch acts as a flux source
        // on the cell-centred heat conduction (AdvanceMaterial), but the
        // elastic eigenstrain below reads cell-centred T through
        // CellToNodeAverage. With Neumann zhi on temp_mf, the ghost row equals
        // the topmost real cell, so the surface node sees a half-cell-lagged
        // temperature instead of the prescribed Dirichlet sp_T(t). Override
        // the top-node temperature inside the patch with sp_T(a_time) so that
        // the eigenstrain at z=z_top reflects the BC Hu prescribes.
        // Step 15a: the top-node eigenstrain override below is meaningful
        // only under mode = prescribed_T. Under convective_flame the surface
        // temperature is whatever the heat solve produces, so the top node
        // reads solved temp_mf through CellToNodeAverage with no override.
        const bool        sp_on  = surface_patch.enabled &&
            surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        const Set::Scalar sp_x0c = surface_patch.x0;
        const Set::Scalar sp_y0c = surface_patch.y0;
        const Set::Scalar sp_r2c = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_Tnow = sp_on ? surface_patch.T_at(a_time) : 0.0;

        if (microstructure_enabled)
        {
            // Step 7/alpha-beta: heterogeneous mechanics. Build each nodal
            // model from cell-centred phase stiffness and cell-centred free
            // strain, both averaged onto the node. The free strain is
            // beta*(T-T_ref) plus an optional phase-local alpha-beta
            // transition strain, so nonlinear mineral transitions are
            // evaluated before nodal averaging.
            // Stiffness uses the (E, mu) → lambda relation from
            // Affine::Isotropic::Parse:
            //     lambda = mu * (E - 2*mu) / (3*mu - E).
            // F0 = eps_free_node * I.
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                Util::RealFillBoundary(*beta_mf[lev],        geom[lev]);
                Util::RealFillBoundary(*E_phase_mf[lev],     geom[lev]);
                Util::RealFillBoundary(*mu_phase_mf[lev],    geom[lev]);
                Util::RealFillBoundary(*T_ref_phase_mf[lev], geom[lev]);
                Util::RealFillBoundary(*alpha_beta_linear_strain_mf[lev], geom[lev]);
                Util::RealFillBoundary(*alpha_beta_T0_mf[lev],            geom[lev]);
                Util::RealFillBoundary(*alpha_beta_width_mf[lev],         geom[lev]);
                UpdateMicrostructureFreeStrain(lev);
                if (damage_enabled)
                    Util::RealFillBoundary(*D_mf[lev], geom[lev]);

                const Set::Scalar* DXl  = geom[lev].CellSize();
                const Set::Scalar* PLOl = geom[lev].ProbLo();
                const amrex::Box domain = geom[lev].Domain();
                const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
                const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
                const int dom_klo = domain.smallEnd(2), dom_khi = domain.bigEnd(2);
                const int sp_topk_cell = domain.bigEnd(2);
                const int sp_topk_node = geom[lev].Domain().bigEnd(2) + 1;
                for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.grownnodaltilebox();
                    amrex::Array4<Model::Solid::Affine::Isotropic> const& model =
                        model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& epsc =
                        thermal_free_strain_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& Ec   = E_phase_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& muc  = mu_phase_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& betc = beta_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& Trc  = T_ref_phase_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& abec = alpha_beta_linear_strain_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& abTc = alpha_beta_T0_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& abWc = alpha_beta_width_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> Dc;
                    if (damage_enabled) Dc = D_mf[lev]->array(mfi);
                    const bool damage_on = damage_enabled;
                    const Set::Scalar stiff_floor = damage_stiffness_floor;
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Set::Scalar epsn = CellToNodeAverageInDomain(
                            epsc, i, j, k, 0,
                            dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                            0.1);
                        if (sp_on && k == sp_topk_node) {
                            Set::Scalar xn = PLOl[0] + Set::Scalar(i) * DXl[0];
                            Set::Scalar yn = PLOl[1] + Set::Scalar(j) * DXl[1];
                            Set::Scalar dxp = xn - sp_x0c, dyp = yn - sp_y0c;
                            if (dxp*dxp + dyp*dyp <= sp_r2c) {
                                Set::Scalar sum = 0.0;
                                int count = 0;
                                for (int di = 0; di <= 1; ++di) {
                                    int ci = i - di;
                                    if (ci < dom_ilo || ci > dom_ihi) continue;
                                    for (int dj = 0; dj <= 1; ++dj) {
                                        int cj = j - dj;
                                        if (cj < dom_jlo || cj > dom_jhi) continue;
                                        sum += ThermalFreeStrain(
                                            sp_Tnow,
                                            betc(ci, cj, sp_topk_cell, 0),
                                            Trc (ci, cj, sp_topk_cell, 0),
                                            abec(ci, cj, sp_topk_cell, 0),
                                            abTc(ci, cj, sp_topk_cell, 0),
                                            abWc(ci, cj, sp_topk_cell, 0));
                                        ++count;
                                    }
                                }
                                if (count > 0) epsn = sum / Set::Scalar(count);
                            }
                        }
                        Set::Scalar En = CellToNodeAverageInDomain(
                            Ec, i, j, k, 0,
                            dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                            1.0e13);
                        Set::Scalar mun = CellToNodeAverageInDomain(
                            muc, i, j, k, 0,
                            dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                            1.0e13);
                        if (damage_on) {
                            Set::Scalar Dn = CellToNodeAverageInDomain(
                                Dc, i, j, k, 0,
                                dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                                1.0);
                            if (!std::isfinite(Dn)) Dn = 0.0;
                            Dn = std::max(Set::Scalar(0.0), std::min(Set::Scalar(1.0), Dn));
                            Set::Scalar factor = std::max(stiff_floor, 1.0 - Dn);
                            En  *= factor;
                            mun *= factor;
                        }
                        Set::Scalar denom = 3.0 * mun - En;
                        // Guard against the degenerate (lambda denominator
                        // small / non-finite moduli) input combinations.
                        Set::Scalar lam = (std::abs(denom) > 1.0e-30)
                                        ? mun * (En - 2.0 * mun) / denom
                                        : 0.0;
                        if (!std::isfinite(lam) || !std::isfinite(mun)) {
                            lam = 0.0; mun = 0.0;
                        }
                        Set::Matrix F0 = epsn * Set::Matrix::Identity();
                        model(i, j, k).Define(mun, lam, F0);
                    });
                }
                Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            }
            return;
        }

        // Homogeneous path: existing alpha[]/eta[] thermal eigenstrain.
        Set::Scalar Tref = T_ref;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            const Set::Scalar* DXl  = geom[lev].CellSize();
            const Set::Scalar* PLOl = geom[lev].ProbLo();
            const amrex::Box domain = geom[lev].Domain();
            const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
            const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
            const int dom_klo = domain.smallEnd(2), dom_khi = domain.bigEnd(2);
            const int sp_topk_node = geom[lev].Domain().bigEnd(2) + 1;
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<Model::Solid::Affine::Isotropic> const& model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& eta  = eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& temp = temp_old_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix F0 = Set::Matrix::Zero();
                    Set::Scalar Tavg = CellToNodeAverageInDomain(
                        temp, i, j, k, 0,
                        dom_ilo, dom_ihi, dom_jlo, dom_jhi, dom_klo, dom_khi,
                        1.0e5);
                    if (sp_on && k == sp_topk_node) {
                        Set::Scalar xn = PLOl[0] + Set::Scalar(i) * DXl[0];
                        Set::Scalar yn = PLOl[1] + Set::Scalar(j) * DXl[1];
                        Set::Scalar dxp = xn - sp_x0c, dyp = yn - sp_y0c;
                        if (dxp*dxp + dyp*dyp <= sp_r2c) Tavg = sp_Tnow;
                    }
                    for (int n = 0; n < eta.nComp(); n++)
                    {
                        F0 += (eta(i, j, k, n) * alpha[n]) * (Tavg - Tref) * Set::Matrix::Identity();
                    }
                    model(i, j, k).F0 = F0;
                });
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

    // Override HeatConduction::Advance with the enthalpy method.
    // Branch: material + implicit → AdvanceMaterialImplicit (BE on diffusion);
    //         material + explicit → AdvanceMaterial (forward Euler);
    //         no material         → AdvanceConstAlpha (step 4 path).
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::Advance(lev, time, dt);
        UpdateDamageAfterMechanics(lev, time, dt);
        UpdateCohesiveAfterMechanics(lev, time, dt);
        // Step 12: Sp = K_I/K_Ic(T) surface diagnostic for the sp_weibull
        // branch. No-op under damage_law. Reads the same nodal stress field
        // produced by the static elastic solve (Step 7) used by the
        // hu_breakage_index f_b pipeline.
        UpdateSpAfterMechanics(lev, time, dt);
        // Step 15c: 4-connected component labeller over the Sp >= 1 surface
        // zone. Default-off everywhere except sp_weibull + weibull.enabled;
        // writes Sp_cluster_id_mf and per-step cluster diagnostics CSV.
        UpdateSpClusters(lev, time, dt);
        // Step 10/11: removal (spall and/or vaporisation) must happen
        // BEFORE the H/Temp swaps below so the reset H(T_amb)/T_amb values
        // become the "old" state consumed by the thermal advance step.
        UpdateRemovalAfterCohesive(lev, time, dt);
        std::swap(*H_mf[lev], *H_old_mf[lev]);
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        if      (microstructure_enabled)                          AdvanceMicrostructure(lev, time, dt);
        else if (material && solver_kind == SolverKind::Implicit) AdvanceMaterialImplicit(lev, time, dt);
        else if (material)                                        AdvanceMaterial(lev, time, dt);
        else                                                      AdvanceConstAlpha(lev, time, dt);
        // Step 6b: every path evolves H_mf as the conserved state, so
        // dispatch a single InvertHtoT here. The microstructure variant
        // selects the per-phase Material::Table by cell phase ID.
        InvertHtoT(lev);
    }

    void AdvanceConstAlpha(int lev, Set::Scalar time, Set::Scalar dt)
    {
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar a_th = therm_alpha;
        const Set::Scalar P_now = beam.power_at(time);
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
        const int         bprof = beam.profile;
        const Set::Scalar brad  = beam.radius;
        auto alpha_f = beam.alpha_f;
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];

        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> H_old = H_old_mf.Patch(lev, mfi);
            Set::Patch<Set::Scalar>       H     = H_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> T_old = temp_old_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> surf  = surface_mf.Patch(lev, mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar lap = a_th * Numeric::Laplacian(T_old, i, j, k, 0, DX);
                Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                Set::Scalar Tv = T_old(i, j, k);
                Set::Scalar Q  = BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                            bR, bprof, brad, alpha_f, x, y, z, Tv,
                                            DX[2]);
                Set::Scalar eps = (Tv <= T_eps_b) ? eps_lo : eps_hi;
                Set::Scalar q_loss = surf(i, j, k)
                                   * (eps * sig_SB * (Tv*Tv*Tv*Tv - T_a*T_a*T_a*T_a)
                                      + h_c * (Tv - T_a));
                H(i, j, k) = H_old(i, j, k) + dt * (lap + Q - q_loss * inv_dz);
            });
        }
    }

    void AdvanceMaterial(int lev, Set::Scalar time, Set::Scalar dt)
    {
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Numeric::Material::Material *mat = material;
        const Set::Scalar P_now = beam.power_at(time);
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
        const int         bprof = beam.profile;
        const Set::Scalar brad  = beam.radius;
        auto alpha_f = beam.alpha_f;
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];
        const Set::Scalar inv_dx2 = 1.0 / (DX[0] * DX[0]);
        const Set::Scalar inv_dy2 = 1.0 / (DX[1] * DX[1]);
        const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);

        // Removed cells are void and must not exchange heat with surviving
        // material (mirrors AdvanceMicrostructure). Without this, a voided cell
        // — reset to T_amb on removal — would act as a full-conductivity cold
        // sink for its neighbours. Fill the `removed` ghost layer so the
        // stencil can detect void neighbours.
        const bool void_aware = removal_enabled();
        if (void_aware) Util::RealFillBoundary(*removed_mf[lev], geom[lev]);
        const amrex::Box& vdomain = geom[lev].Domain();
        const int vdlo[3] = {vdomain.smallEnd(0), vdomain.smallEnd(1), vdomain.smallEnd(2)};
        const int vdhi[3] = {vdomain.bigEnd(0),   vdomain.bigEnd(1),   vdomain.bigEnd(2)};

        // Step 6d / Step 15a: surface patch on the +z face. prescribed_T
        // injects the Dirichlet face flux; convective_flame injects the Robin
        // flame-jet flux q = h_fl*(T_flame - T_cell).
        const bool        sp_on    = surface_patch.enabled;
        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
        const int         sp_topk  = geom[lev].Domain().bigEnd(2);
        const Set::Scalar sp_x0    = surface_patch.x0;
        const Set::Scalar sp_y0    = surface_patch.y0;
        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_T     = (sp_on && sp_presc)
                                   ? surface_patch.T_at(time) : 0.0;
        const Set::Scalar sp_hfl   = surface_patch.h_conv;
        const Set::Scalar sp_Tfl   = surface_patch.T_flame;

        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> H_old = H_old_mf.Patch(lev, mfi);
            Set::Patch<Set::Scalar>       H     = H_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> T_old = temp_old_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> surf  = surface_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> rem   =
                void_aware ? removed_mf.Patch(lev, mfi) : surface_mf.Patch(lev, mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tc = T_old(i, j, k);
                // Void cells (already removed) hold no material: freeze them at
                // their reset state and exchange no heat with neighbours.
                const bool self_void = void_aware && rem(i, j, k) > 0.5;
                if (self_void) { H(i, j, k) = H_old(i, j, k); return; }
                Set::Scalar Txm = T_old(i-1, j, k), Txp = T_old(i+1, j, k);
                Set::Scalar Tym = T_old(i, j-1, k), Typ = T_old(i, j+1, k);
                Set::Scalar Tzm = T_old(i, j, k-1), Tzp = T_old(i, j, k+1);
                Set::Scalar kc  = mat->kappa(Tc);
                // Harmonic-average κ on faces (handles discontinuous κ across
                // interfaces correctly; degenerates to arithmetic mean for
                // smooth κ(T)). A face touching a void neighbour carries zero
                // flux (insulating cavity wall); domain-boundary faces keep the
                // normal stencil (rem ghosts there are not meaningful).
                auto kface = [](Set::Scalar a, Set::Scalar b) {
                    return (a + b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
                };
                auto nvoid = [&](int di, int dj, int dk) -> bool {
                    const int ii = i+di, jj = j+dj, kk = k+dk;
                    if (!void_aware) return false;
                    if (ii < vdlo[0] || ii > vdhi[0] || jj < vdlo[1] ||
                        jj > vdhi[1] || kk < vdlo[2] || kk > vdhi[2]) return false;
                    return rem(ii, jj, kk) > 0.5;
                };
                Set::Scalar fxp = nvoid( 1,0,0) ? 0.0 : kface(kc, mat->kappa(Txp)) * (Txp - Tc) * inv_dx2;
                Set::Scalar fxm = nvoid(-1,0,0) ? 0.0 : kface(kc, mat->kappa(Txm)) * (Tc - Txm) * inv_dx2;
                Set::Scalar fyp = nvoid(0, 1,0) ? 0.0 : kface(kc, mat->kappa(Typ)) * (Typ - Tc) * inv_dy2;
                Set::Scalar fym = nvoid(0,-1,0) ? 0.0 : kface(kc, mat->kappa(Tym)) * (Tc - Tym) * inv_dy2;
                Set::Scalar fzp = nvoid(0,0, 1) ? 0.0 : kface(kc, mat->kappa(Tzp)) * (Tzp - Tc) * inv_dz2;
                Set::Scalar fzm = nvoid(0,0,-1) ? 0.0 : kface(kc, mat->kappa(Tzm)) * (Tc - Tzm) * inv_dz2;
                Set::Scalar div_kgrad = (fxp - fxm) + (fyp - fym) + (fzp - fzm);

                Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];

                // Step 6d: top-row patch contribution. The patch acts as a
                // Dirichlet face value at distance dz/2 from the cell centre,
                // matching the existing FV ghost convention used by the rest
                // of the stencil (ghost = T_ghost; flux = k*(T_ghost-T)/dx).
                // With hc.bc.temp.zhi = neumann the existing fzp is zero on
                // top cells, so this term is the entire +z face flux there;
                // outside the patch the Neumann BC continues to apply.
                if (sp_on && k == sp_topk) {
                    Set::Scalar dxp = x - sp_x0, dyp = y - sp_y0;
                    if (dxp * dxp + dyp * dyp <= sp_r2) {
                        if (sp_presc)
                            div_kgrad += kc * (sp_T - Tc) * inv_dz2;
                        else
                            div_kgrad += sp_hfl * (sp_Tfl - Tc) * inv_dz;
                    }
                }

                Set::Scalar Q = BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                           bR, bprof, brad, alpha_f, x, y, z, Tc,
                                           DX[2]);

                Set::Scalar eps = (Tc <= T_eps_b) ? eps_lo : eps_hi;
                Set::Scalar q_loss = surf(i, j, k)
                                   * (eps * sig_SB * (Tc*Tc*Tc*Tc - T_a*T_a*T_a*T_a)
                                      + h_c * (Tc - T_a));

                H(i, j, k) = H_old(i, j, k) + dt * (div_kgrad + Q - q_loss * inv_dz);
            });
        }
    }

    // Backward-Euler diffusion via amrex::MLABecLaplacian:
    //   (ρcp/dt) · T_new − ∇·(κ ∇T_new) = (ρcp/dt) · T_old + Q − q_loss/dz
    // ρ, c_p, κ are evaluated at T_old (one-step linearisation). The source
    // term and the radiative loss stay explicit in T_old. Only single-level
    // (max_level = 0) is currently supported. Domain BCs are read from the
    // BC handler attached to temp_mf so the existing input syntax continues
    // to work.
    void AdvanceMaterialImplicit(int lev, Set::Scalar time, Set::Scalar dt)
    {
        Util::Assert(INFO, TEST(material != nullptr));
        Util::Assert(INFO, TEST(finest_level == 0));

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Numeric::Material::Material *mat = material;
        const Set::Scalar P_now = beam.power_at(time);
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
        const int         bprof = beam.profile;
        const Set::Scalar brad  = beam.radius;
        auto alpha_f = beam.alpha_f;
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];

        const amrex::BoxArray ba = grids[lev];
        const amrex::DistributionMapping &dm = dmap[lev];

        amrex::MultiFab A_mf  (ba, dm, 1, 0);
        amrex::MultiFab rhs_mf(ba, dm, 1, 0);
        amrex::MultiFab T_new (ba, dm, 1, 1);
        std::array<amrex::MultiFab, AMREX_SPACEDIM> B_mf;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray ba_face = ba;
            ba_face.surroundingNodes(d);
            B_mf[d].define(ba_face, dm, 1, 0);
        }

        // Build A = ρcp/dt and rhs = (ρcp/dt)·T_old + Q − q_loss/dz.
        for (amrex::MFIter mfi(A_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& surf  = (*surface_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> A   = A_mf.array(mfi);
            amrex::Array4<Set::Scalar> rhs = rhs_mf.array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tc = T_old(i, j, k);
                // Use ρ·cp_eff (with smoothed latent-heat pulses) as the
                // implicit thermal-mass coefficient — this is the slope
                // dH/dT used by the H↔T table, so the linear solver "sees"
                // the same effective heat capacity that the enthalpy
                // formulation would. Without it, cells crossing the
                // fusion / vaporisation plateau in one step skip the
                // latent-heat absorption.
                Set::Scalar rhocp = mat->effective_rhocp(Tc);
                Set::Scalar adiag = rhocp / dt;
                Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                Set::Scalar Q = BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                           bR, bprof, brad, alpha_f, x, y, z, Tc,
                                           DX[2]);
                Set::Scalar eps = (Tc <= T_eps_b) ? eps_lo : eps_hi;
                Set::Scalar q_loss = surf(i, j, k)
                                   * (eps * sig_SB * (Tc*Tc*Tc*Tc - T_a*T_a*T_a*T_a)
                                      + h_c * (Tc - T_a));
                A(i, j, k)   = adiag;
                rhs(i, j, k) = adiag * Tc + Q - q_loss * inv_dz;
            });
        }

        // Build face-centred B = harmonic-averaged κ from cell T_old.
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            const int sx = (d == 0), sy = (d == 1), sz = (d == 2);
            for (amrex::MFIter mfi(B_mf[d], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> B = B_mf[d].array(mfi);
                amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                    Set::Scalar kl = mat->kappa(T_old(i - sx, j - sy, k - sz));
                    Set::Scalar kr = mat->kappa(T_old(i, j, k));
                    B(i, j, k) = (kl + kr > 0.0) ? 2.0 * kl * kr / (kl + kr) : 0.0;
                });
            }
        }

        amrex::LPInfo info;
        info.setMetricTerm(false);
        amrex::MLABecLaplacian linop({geom[lev]}, {ba}, {dm}, info);
        linop.setMaxOrder(2);

        // Domain BCs: hardcoded Neumann (zero-flux) on every face. This
        // matches the adiabatic side/bottom + radiation/convection-as-source
        // top setup used by the Zhang test and any other purely thermal-loss
        // configuration. (HeatConduction's BC pointer is private, so we
        // can't introspect the user's BC handler from here without touching
        // that file. Step-2 inputs that need Dirichlet BCs use the const-α
        // explicit path anyway.)
        std::array<amrex::LinOpBCType, AMREX_SPACEDIM> lobc, hibc;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            lobc[d] = hibc[d] = amrex::LinOpBCType::Neumann;
        linop.setDomainBC(lobc, hibc);

        // setLevelBC needs ghost cells filled (used for Dirichlet values).
        linop.setLevelBC(0, &(*temp_old_mf[lev]));
        linop.setScalars(1.0, 1.0);
        linop.setACoeffs(0, A_mf);
        amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> B_ptrs{
            AMREX_D_DECL(&B_mf[0], &B_mf[1], &B_mf[2])};
        linop.setBCoeffs(0, B_ptrs);

        amrex::MLMG mlmg(linop);
        mlmg.setMaxIter(implicit_max_iter);
        mlmg.setVerbose(implicit_verbose);
        mlmg.setBottomVerbose(0);

        // Initial guess T_new = T_old.
        amrex::MultiFab::Copy(T_new, *temp_old_mf[lev], 0, 0, 1, 0);
        mlmg.solve({&T_new}, {&rhs_mf}, implicit_tol_rel, implicit_tol_abs);

        // Push T_new back into the H field via the table; InvertHtoT (called
        // by Advance) then re-derives temp_mf and the phase fractions.
        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& Tn = T_new.array(mfi);
            amrex::Array4<Set::Scalar> H = (*H_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=, this] (int i, int j, int k) {
                H(i, j, k) = material_table.H_of_T(Tn(i, j, k));
            });
        }
    }

    // Step 6b: heterogeneous-property explicit enthalpy update for the
    // microstructure path. H is the conserved variable:
    //
    //   ∂H/∂t = ∇·(κ_eff(x,D) ∇T) + Q − q_loss/dz
    //
    // where κ_eff(x,D) = κ_phase(x) · exp(−α_D · D). Face conductivities are
    // built from the harmonic mean of the cell-centred κ_eff values, with an
    // optional grain-boundary contact-conductance term `h_gb(D) = h_gb0·(1−D)`
    // added in series whenever the two cells across the face have different
    // phase IDs.
    //
    // T is recovered cell-by-cell in InvertHtoT using the per-phase
    // Numeric::Material::Table (built from the phase's Constant material).
    // For phases without phase change this reduces to T_new = T_old +
    // dt·(div+Q-loss)/(rho*Cp), matching the Step 6 behaviour exactly.
    void AdvanceMicrostructure(int lev, Set::Scalar time, Set::Scalar dt)
    {
        // Refresh damage ghost cells (D may be heterogeneous; D itself does
        // not evolve in Step 6 but its ghost layer still has to be filled
        // before the harmonic-average face stencil reads neighbour values).
        Util::RealFillBoundary(*D_mf[lev],         geom[lev]);
        Util::RealFillBoundary(*phase_mf[lev],     geom[lev]);
        Util::RealFillBoundary(*grain_id_mf[lev],  geom[lev]);
        Util::RealFillBoundary(*kappa_mf[lev],     geom[lev]);
        // Step 10: if spall is active, removed cells are void and must not
        // exchange heat with surviving material. Filling the ghost layer of
        // `removed` lets the face stencil below detect a void neighbour
        // through ParmParse-style ghost reads.
        const bool void_aware = removal_enabled();
        if (void_aware) Util::RealFillBoundary(*removed_mf[lev], geom[lev]);

        // Domain extents for the face stencil's boundary-face fallback below.
        const amrex::Box& domain = geom[lev].Domain();
        const int dlo[3] = {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
        const int dhi[3] = {domain.bigEnd(0),   domain.bigEnd(1),   domain.bigEnd(2)};

        // E1: global top solid cell per (i, j) column for the beam path (see
        // the BeamSource call below). Empty unless void-aware with the beam on,
        // or (S1) with surface.follow_mask / the energy ledger on.
        const Set::Scalar P_now = beam.power_at(time);
        const bool follow_mask = surface_follow_mask;
        const amrex::Vector<int> k_top_col =
            (void_aware && (P_now > 0.0 || follow_mask || ledger_enabled))
            ? ColumnTopSolid(lev) : amrex::Vector<int>();
        const int* ktc = k_top_col.data();
        const int ncx = dhi[0] - dlo[0] + 1;

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
        const int         bprof = beam.profile;
        const Set::Scalar brad  = beam.radius;
        auto alpha_f = beam.alpha_f;
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];
        const Set::Scalar inv_dx2 = 1.0 / (DX[0] * DX[0]);
        const Set::Scalar inv_dy2 = 1.0 / (DX[1] * DX[1]);
        const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);
        const Set::Scalar a_D     = kappa_damage_alpha;
        const Set::Scalar h_gb0   = h_gb_conductance;
        const Set::Scalar dx_ax[3] = {DX[0], DX[1], DX[2]};

        // Step 11b: the Hu onset validation uses a homogeneous one-phase
        // microstructure so the DP damage fields exist. Apply the same
        // surface patch here as in AdvanceMaterial(): prescribed_T injects the
        // Dirichlet face flux; convective_flame (Step 15a) injects the Robin
        // flame-jet flux q = h_fl*(T_flame - T_cell).
        const bool        sp_on    = surface_patch.enabled;
        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
        const int         sp_topk  = geom[lev].Domain().bigEnd(2);
        const Set::Scalar sp_x0    = surface_patch.x0;
        const Set::Scalar sp_y0    = surface_patch.y0;
        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_T     = (sp_on && sp_presc)
                                   ? surface_patch.T_at(time) : 0.0;
        // D2a: one per-column (h, T_gas) evaluation per call, shared by the
        // kernel, the pin rule, the ledger and the diagnostics.
        if (sp_on) BuildFlameColumns(lev, time, dt, k_top_col);
        const Set::Scalar* fh  = sp_on ? flame_h.data()  : nullptr;
        const Set::Scalar* fTg = sp_on ? flame_Tg.data() : nullptr;
        // C1: pinned runs the face Newton (face_form = true) plus the pin.
        const bool        sp_pinned = sp_on && surface_patch.robin_form ==
                                      SurfacePatch::RobinForm::Pinned;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face || sp_pinned;
        if (sp_pinned && lev == 0) BuildPinEffective(lev, time);
        const Set::Scalar* tpe = (sp_pinned && lev == 0) ? pin_T_eff.data() : nullptr;

        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> T_old = temp_old_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> H_old = H_old_mf   .Patch(lev, mfi);
            Set::Patch<Set::Scalar>       H     = H_mf       .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> ph    = phase_mf   .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> gid   = grain_id_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> kp    = kappa_mf   .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> Df    = D_mf       .Patch(lev, mfi);
            Set::Patch<Set::Scalar>       keff_out = kappa_eff_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> surf  = surface_mf.Patch(lev, mfi);
            // When void_aware is false, rem and phi are bound to surface_mf
            // as harmless placeholders that the kernel never reads (phi_mf is
            // only registered when removal is enabled).
            Set::Patch<const Set::Scalar> rem   =
                void_aware ? removed_mf.Patch(lev, mfi) : surface_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> phi   =
                void_aware ? phi_mf.Patch(lev, mfi) : surface_mf.Patch(lev, mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tc = T_old(i, j, k);
                Set::Scalar Dc = Df(i, j, k);
                Set::Scalar kc = kp(i, j, k) * std::exp(-a_D * Dc);
                const bool self_void = void_aware && rem(i, j, k) > 0.5;
                if (self_void) kc = 0.0;
                int pc = static_cast<int>(ph(i, j, k));
                int gc = static_cast<int>(gid(i, j, k));
                (void)pc;

                auto k_face = [&](int di, int dj, int dk, int axis) {
                    int ii = i + di, jj = j + dj, kk = k + dk;
                    // Domain-boundary face: kappa_mf and D_mf have no BC
                    // handler so their ghost cells are uninitialised.
                    // Use the cell-centred κ_eff directly so the Dirichlet
                    // ghost in T_old still drives a flux of
                    //   q ≈ kc · (T_BC − T_cell) / dx
                    // (the standard FV boundary treatment when the BC value
                    // sits in the ghost cell at distance dx).
                    if (ii < dlo[0] || ii > dhi[0] ||
                        jj < dlo[1] || jj > dhi[1] ||
                        kk < dlo[2] || kk > dhi[2]) return kc;
                    if (void_aware && rem(ii, jj, kk) > 0.5) return Set::Scalar(0.0);
                    Set::Scalar k_n = kp(ii, jj, kk)
                                    * std::exp(-a_D * Df(ii, jj, kk));
                    int p_n = static_cast<int>(ph(ii, jj, kk));
                    int g_n = static_cast<int>(gid(ii, jj, kk));
                    (void)p_n; // Step 6c: GB hook now keys on grain_id, not phase.
                    Set::Scalar inv = 0.0;
                    if (kc > 0.0 && k_n > 0.0)
                        inv = 0.5 / kc + 0.5 / k_n;
                    else
                        return Set::Scalar(0.0);
                    // Step 6c: physical grain-boundary contact resistance is
                    // keyed on true grain identity, not on mineral phase. Two
                    // adjacent same-mineral grains still see h_gb0.
                    if (h_gb0 > 0.0 && g_n != gc)
                    {
                        Set::Scalar D_face = 0.5 * (Dc + Df(ii, jj, kk));
                        Set::Scalar h_eff  = h_gb0 * (1.0 - D_face);
                        if (h_eff > 0.0) inv += 1.0 / (h_eff * dx_ax[axis]);
                        else             return Set::Scalar(0.0);
                    }
                    return 1.0 / inv;
                };

                Set::Scalar k_xp = k_face( 1, 0, 0, 0);
                Set::Scalar k_xm = k_face(-1, 0, 0, 0);
                Set::Scalar k_yp = k_face( 0, 1, 0, 1);
                Set::Scalar k_ym = k_face( 0,-1, 0, 1);
                Set::Scalar k_zp = k_face( 0, 0, 1, 2);
                Set::Scalar k_zm = k_face( 0, 0,-1, 2);

                Set::Scalar div_kgrad =
                    (k_xp * (T_old(i+1,j,k) - Tc) - k_xm * (Tc - T_old(i-1,j,k))) * inv_dx2 +
                    (k_yp * (T_old(i,j+1,k) - Tc) - k_ym * (Tc - T_old(i,j-1,k))) * inv_dy2 +
                    (k_zp * (T_old(i,j,k+1) - Tc) - k_zm * (Tc - T_old(i,j,k-1))) * inv_dz2;

                Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];

                // S1: with surface.follow_mask the patch and the losses sit on
                // the mask column top; otherwise the fixed domain top and the
                // phi-window surface_mf (unchanged).
                const int k_top_c = (follow_mask && ktc)
                                  ? ktc[(j - dlo[1]) * ncx + (i - dlo[0])] : sp_topk;
                const Set::Scalar surf_c = follow_mask
                                         ? Set::Scalar(k == k_top_c ? 1.0 : 0.0)
                                         : surf(i, j, k);
                bool in_patch = false;
                if (!self_void && sp_on && k == k_top_c) {
                    Set::Scalar dxp = x - sp_x0, dyp = y - sp_y0;
                    if (dxp * dxp + dyp * dyp <= sp_r2) {
                        if (sp_presc)
                            div_kgrad += kc * (sp_T - Tc) * inv_dz2;
                        else
                            in_patch = true;
                    }
                }
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                const int col_c = (j - dlo[1]) * ncx + (i - dlo[0]);
                if (!self_void)
                    SurfaceCellFlux(in_patch, sp_face, surf_c,
                                    fh  ? fh[col_c]  : Set::Scalar(0.0),
                                    fTg ? fTg[col_c] : Set::Scalar(0.0),
                                    kc, DX[2], Tc,
                                    eps_lo, eps_hi, T_eps_b, sig_SB, h_c, T_a,
                                    q_in, q_out, q_rob, q_los,
                                    tpe ? tpe[col_c] : Set::Scalar(0.0));
                if (in_patch) div_kgrad += q_in * inv_dz;

                // E1: Beer-Lambert path is measured through material from the
                // TOP FACE OF THE FIRST SOLID CELL, so it matches the `removed`
                // mask. Removal keeps solid <=> phi >= 0 with phi linear down
                // the column, so the top solid cell has phi in [0, dz) and a
                // cell m cells below it has n_above = floor(phi/dz) == m; its
                // path interval [m*dz, (m+1)*dz] tiles the column with no gap
                // and no overlap for ANY phi_top. Passing raw phi instead left
                // the slice between the level set and the top solid face inside
                // a void cell (Q = 0): attenuated but never deposited. Before any
                // removal phi = (n + 1/2)*dz, so path == phi to round-off.
                // n_above is taken from the mask (k_top - k), which equals
                // floor(phi/dz) except at round-off ties phi_top ~ 0 / ~ dz,
                // where the floor picks the wrong cell (see SolidFacePath).
                Set::Scalar Q = self_void ? Set::Scalar(0.0)
                                          : BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                                       bR, bprof, brad, alpha_f, x, y, z, Tc,
                                                       DX[2],
                                                       // ktc is null unless void-aware (Q = 0 when P_now == 0).
                                                       (void_aware && ktc)
                                                         ? SolidFacePath(ktc[(j - dlo[1]) * ncx + (i - dlo[0])] - k, DX[2])
                                                         : Set::Scalar(-1.0));

                (void)q_rob; (void)q_los;
                const Set::Scalar q_loss = q_out;

                // Conservative H update; T is recovered in InvertHtoT
                // using the per-phase Numeric::Material::Table. Detached
                // (void) cells exchange no heat and get no source, so their
                // H stays frozen at its value when voided (removal never
                // writes H).
                H(i, j, k) = H_old(i, j, k) + dt * (div_kgrad + Q - q_loss * inv_dz);
                keff_out(i, j, k) = kc; // diagnostic, refreshed each step
            });
        }

        // E1: standing beam-deposition closure counter (thermo.dat).
        if (P_now > 0.0) BeamEnergyBalance(lev, time, dt, k_top_col);
        // S1: whole-domain energy ledger (thermo.dat), after the H update.
        if (ledger_enabled) EnergyLedger(lev, time, dt, k_top_col);
        // C1/D2a: pinned-closure branch counts and spatial-flame diagnostics.
        if (lev == 0 && (sp_pinned || (sp_on && surface_patch.use_expr)))
            PatchDiagnostics(lev, k_top_col);
    }

    // D2a: per-column flame for this AdvanceMicrostructure call.
    //   constants: flame_h = h_conv, flame_Tg = T_flame (bit-identical path).
    //   expressions: at the column centre (x, y), r = |(x, y) - (x0, y0)|,
    //     z_face = PLO_z + (k_top + 1)*dz (mask top face), z_n = FeetDescent
    //     (D2b, nozzle_descent = feet) or nozzle_z0 -
    //     nozzle_feed*time, s = z_n - z_face; h = h_expr(x,y,r,s,t),
    //     T_gas = T_flame_expr(...). Columns with no solid cell get 0.
    // Aborts if a live column with r <= nozzle_collision_radius (default: the
    // patch radius) has s <= 0 (nozzle in the rock), or if an in-patch column
    // gets a non-finite / non-positive h or T_gas. patch_min_standoff is the
    // min s over live columns with r <= nozzle_collision_radius (in the patch).
    // k_top_col is the reduced
    // ColumnTopSolid, so every rank builds identical vectors.
    void BuildFlameColumns(int lev, Set::Scalar time, Set::Scalar dt,
                           const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        if (!surface_patch.use_expr)
        {
            flame_h.assign(ncols, surface_patch.h_conv);
            flame_Tg.assign(ncols, surface_patch.T_flame);
            flame_s.assign(ncols, 0.0);
            return;
        }
        flame_h.assign(ncols, 0.0);
        flame_Tg.assign(ncols, 0.0);
        flame_s.assign(ncols, 0.0);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;
        const Set::Scalar rc2   = surface_patch.nozzle_collision_radius *
                                  surface_patch.nozzle_collision_radius;
        const Set::Scalar z_n = surface_patch.feet
                              ? FeetDescent(lev, dt, k_top_col)
                              : surface_patch.nozzle_z0 - surface_patch.nozzle_feed * time;
        Set::Scalar s_min = std::numeric_limits<Set::Scalar>::infinity();
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t col = static_cast<std::size_t>(j) * nx + i;
            const int k_top = k_top_col[col];
            if (k_top == std::numeric_limits<int>::min()) continue;
            const Set::Scalar x = PLO[0] + (Set::Scalar(i + domain.smallEnd(0)) + 0.5) * DX[0];
            const Set::Scalar y = PLO[1] + (Set::Scalar(j + domain.smallEnd(1)) + 0.5) * DX[1];
            const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
            const Set::Scalar r2 = dxp * dxp + dyp * dyp;
            const Set::Scalar r = std::sqrt(r2);
            const Set::Scalar z_face = PLO[2] + Set::Scalar(k_top + 1) * DX[2];
            const Set::Scalar s = z_n - z_face;
            flame_s[col] = s;
            if (r2 > r2max) continue;
            const bool in_coll = r2 <= rc2;
            if (in_coll && !(s > 0.0))
                Util::Abort(INFO, "surface_patch nozzle collision: stand-off s = ", s,
                            " <= 0 at column (i, j) = (", i + domain.smallEnd(0), ", ",
                            j + domain.smallEnd(1), "), time = ", time,
                            " (z_nozzle = ", z_n, ", z_face = ", z_face, ")");
            if (surface_patch.jet_enthalpy)
            {
                // D2a2: no jet flux above the nozzle plane (s <= 0); T_gas is
                // set by JetEnthalpyMarch (T_ent here, a harmless finite value).
                const Set::Scalar hj = (s > 0.0) ? surface_patch.h_f(x, y, r, s, time)
                                                 : Set::Scalar(0.0);
                if (!(std::isfinite(hj) && hj >= 0.0))
                    Util::Abort(INFO, "surface_patch.h_expr must be finite and >= 0 in the "
                                "patch; got h = ", hj, " at column (", i + domain.smallEnd(0),
                                ", ", j + domain.smallEnd(1), "), r = ", r, ", s = ", s,
                                ", time = ", time);
                flame_h[col]  = hj;
                flame_Tg[col] = surface_patch.jet_T_ent;
                if (in_coll) s_min = std::min(s_min, s);
                continue;
            }
            const Set::Scalar h  = surface_patch.h_f(x, y, r, s, time);
            const Set::Scalar Tg = surface_patch.Tg_f(x, y, r, s, time);
            if (!(std::isfinite(h) && h > 0.0 && std::isfinite(Tg) && Tg > 0.0))
                Util::Abort(INFO, "surface_patch.h_expr / T_flame_expr must be finite and > 0 "
                            "in the patch; got h = ", h, ", T_gas = ", Tg, " at column (",
                            i + domain.smallEnd(0), ", ", j + domain.smallEnd(1), "), r = ", r,
                            ", s = ", s, ", time = ", time);
            flame_h[col]  = h;
            flame_Tg[col] = Tg;
            if (in_coll) s_min = std::min(s_min, s);
        }
        patch_min_standoff = std::isfinite(s_min) ? s_min : 0.0;
        if (surface_patch.jet_enthalpy) JetEnthalpyMarch(lev, time, k_top_col);
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
    //   max / mean: the whole-annulus reduction (one "pad").
    // z_target = z_foot + foot_standoff; z_n = min(z_n_prev, max(z_target,
    // z_n_prev - nozzle_feed_max*dt)) (never rises; first step z_n =
    // z_target). foot_carry_cols counts columns with z_face >= their pad's
    // height; foot_stall_time is the simulated time since z_n last decreased
    // (0 on a step that descends). A stall is logged, never aborted. k_top_col is the domain-wide reduced
    // mask top, so every rank computes the same z_n. Aborts if the annulus or
    // a pad has no live column. nozzle_z is NOT checkpointed.
    Set::Scalar FeetDescent(int lev, Set::Scalar dt, const amrex::Vector<int>& k_top_col)
    {
        if (lev != 0)
            Util::Abort(INFO, "surface_patch.nozzle_descent = feet is single-level only");
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const auto& sp = surface_patch;
        const Set::Scalar ri2 = sp.foot_r_inner * sp.foot_r_inner;
        const Set::Scalar ro2 = sp.foot_r_outer * sp.foot_r_outer;
        std::vector<Set::Scalar> zf, az;
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
        }
        if (zf.empty())
            Util::Abort(INFO, "surface_patch feet: the foot annulus [", sp.foot_r_inner, ", ",
                        sp.foot_r_outer, "] m has no live column");
        const bool pads = sp.foot_rule == SurfacePatch::FootRule::Pads;
        const int npads = pads ? sp.foot_npads : 1;
        std::vector<int> pad(zf.size(), 0);
        if (pads && npads > 1)
        {
            const Set::Scalar a_min = *std::min_element(az.begin(), az.end());
            const Set::Scalar a_max = *std::max_element(az.begin(), az.end());
            if (a_max > a_min)
                for (std::size_t c = 0; c < zf.size(); ++c)
                    pad[c] = std::min(npads - 1, static_cast<int>(
                                 std::floor((az[c] - a_min) / (a_max - a_min) * npads)));
        }
        std::vector<Set::Scalar> pad_h(npads, 0.0);
        Set::Scalar z_foot = 0.0;
        for (int p = 0; p < npads; ++p)
        {
            std::vector<Set::Scalar> v;
            for (std::size_t c = 0; c < zf.size(); ++c) if (pad[c] == p) v.push_back(zf[c]);
            if (v.empty())
                Util::Abort(INFO, "surface_patch feet: pad ", p, " of ", npads,
                            " has no live column in the foot annulus");
            std::sort(v.begin(), v.end());
            Set::Scalar hp = 0.0;
            if (sp.foot_rule == SurfacePatch::FootRule::Max)
                hp = v.back();
            else if (sp.foot_rule == SurfacePatch::FootRule::Mean)
            {
                for (Set::Scalar z : v) hp += z;
                hp /= Set::Scalar(v.size());
            }
            else
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
    // from PinRule at this h and T_gas), so the surface temperature and branch
    // match the kernel: BuildPinEffective, which runs next, reproduces the pin
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
    void JetEnthalpyMarch(int lev, Set::Scalar time, const amrex::Vector<int>& k_top_col)
    {
        if (lev != 0)
            Util::Abort(INFO, "surface_patch.jet_closure = enthalpy is single-level only");
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0), ny = domain.length(1);
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;
        const bool pinned = surface_patch.robin_form == SurfacePatch::RobinForm::Pinned;
        const bool face   = surface_patch.robin_form == SurfacePatch::RobinForm::Face || pinned;

        // Top-cell T_old and k_c per column (one owner per column; exact sum).
        std::vector<Set::Scalar> jTc(ncols, 0.0), jkc(ncols, 0.0);
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp    = (*kappa_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& Df    = (*D_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t col = static_cast<std::size_t>(j - jlo_dom) * nx + (i - ilo_dom);
                const int k = k_top_col[col];
                if (k < lo.z || k > hi.z) continue;
                jTc[col] = T_old(i, j, k);
                jkc[col] = kp(i, j, k) * std::exp(-kappa_damage_alpha * Df(i, j, k));
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(jTc.data(), static_cast<int>(ncols));
        amrex::ParallelDescriptor::ReduceRealSum(jkc.data(), static_cast<int>(ncols));

        // Bin the marching columns by r.
        const Set::Scalar dr = surface_patch.jet_dr > 0.0 ? surface_patch.jet_dr
                                                          : std::max(DX[0], DX[1]);
        std::vector<std::vector<std::size_t>> bins;
        std::vector<Set::Scalar> col_r(ncols, 0.0);
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
        Set::Scalar phi = 1.0;
        if (surface_patch.jet_stag_decay && s_c > 0.0)
            phi = std::min(Set::Scalar(1.0),
                           surface_patch.jet_core_length * surface_patch.jet_D / s_c);
        // D2b: exhaust recirculation mixes the nozzle flow with the previous
        // step's exhaust (T_rec) instead of T_ent; the floor stays T_ent.
        const bool exhaust = surface_patch.jet_T_ent_exhaust;
        if (exhaust && !jet_rec_init) { jet_T_rec_next = T_ent; jet_rec_init = true; }
        const Set::Scalar T_mix = exhaust ? jet_T_rec_next : T_ent;
        const Set::Scalar T_stag = T_mix + (surface_patch.jet_T_nozzle - T_mix) * phi;
        const Set::Scalar m_jet  = surface_patch.jet_entrained_mass
                                 ? surface_patch.jet_mdot / phi : surface_patch.jet_mdot;
        const bool have_pin = pinned && pin_T.size() == ncols;
        const Set::Scalar mcp   = m_jet * surface_patch.jet_cp;
        const Set::Scalar dA    = DX[0] * DX[1];
        Set::Scalar T_gas = T_stag;
        Set::Scalar P_face = 0.0, r_reach = 0.0;
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
        for (std::size_t ib = 0; ib < bins.size(); ++ib)
        {
            const auto& bin = bins[ib];
            if (bin.empty()) continue;
            const Set::Scalar T_in = T_gas;
            Set::Scalar P_bin = 0.0;
            for (std::size_t c : bin)
            {
                flame_Tg[c] = T_gas;
                int reason = 0;
                const Set::Scalar Tpin = have_pin
                    ? PinRule(c, flame_h[c], T_gas, time, DX[2], reason) : Set::Scalar(0.0);
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                SurfaceCellFlux(true, face, 1.0, flame_h[c], T_gas,
                                jkc[c], DX[2], jTc[c], eps_rad, eps_rad_hi, T_eps_break,
                                sigma_SB, h_conv, T_amb, q_in, q_out, q_rob, q_los, Tpin);
                if (q_rob > 0.0) P_bin += q_rob;
            }
            P_bin *= dA;
            if (P_bin > 0.0)
                for (std::size_t c : bin) r_reach = std::max(r_reach, col_r[c]);
            const Set::Scalar avail = mcp * (T_gas - T_ent);
            const Set::Scalar share = (P_bin >= avail) ? avail : P_bin;
            if (P_bin >= avail) { P_face += avail; T_gas = T_ent; }
            else                { P_face += P_bin; T_gas = T_gas - P_bin / mcp; }
            if (prof)
                prow << time << "," << ib << "," << Set::Scalar(ib) * dr << "," << bin.size()
                     << "," << T_in << "," << P_bin << "," << share << "\n";
        }
        if (prof && amrex::ParallelDescriptor::IOProcessor())
        {
            const std::string path = plot_file + "_jet_profile.csv";
            const bool first = !jet_profile_header_written;
            std::ofstream ofs(path, (first ? std::ios::trunc : std::ios::app) | std::ios::out);
            if (ofs.is_open())
            {
                if (first) ofs << "time,bin,r_lo,n_cols,T_gas,P_bin,P_share\n";
                ofs << prow.str();
            }
        }
        if (prof) jet_profile_header_written = true;
        const Set::Scalar P_nozzle = surface_patch.jet_mdot * surface_patch.jet_cp
                                   * (surface_patch.jet_T_nozzle - T_ent);
        jet_P_cap     = mcp * (T_stag - T_ent);
        jet_T_stag    = T_stag;
        jet_s_c       = s_c;
        Set::Scalar P_budget = P_nozzle;
        if (exhaust)
        {
            // P_face + P_exhaust = mdot*cp*(T_nozzle - T_ent) + P_recirc
            jet_T_rec    = T_mix;
            jet_P_recirc = (m_jet - surface_patch.jet_mdot) * surface_patch.jet_cp * (T_mix - T_ent);
            P_budget     = P_nozzle + jet_P_recirc;
            jet_T_rec_next = T_gas;
        }
        jet_P_decay   = P_budget - jet_P_cap;
        jet_P_face    = P_face;
        jet_T_exhaust = T_gas;
        jet_P_exhaust = mcp * (T_gas - T_ent);
        jet_r_reach   = r_reach;
        const Set::Scalar tol = 1.0e-9 * std::max(P_budget, Set::Scalar(1.0));
        if (!(std::isfinite(P_face) && P_face >= 0.0 && P_face <= jet_P_cap + tol &&
              std::abs(P_face + jet_P_exhaust - jet_P_cap) <= tol &&
              jet_P_decay >= -tol &&
              (surface_patch.jet_stag_decay && !surface_patch.jet_entrained_mass
                   ? true : std::abs(jet_P_decay) <= tol)))
            Util::Abort(INFO, "jet enthalpy invariant violated at time ", time,
                        ": jet_P_face = ", P_face, ", jet_P_exhaust = ", jet_P_exhaust,
                        ", jet_P_cap = ", jet_P_cap, ", jet_P_decay = ", jet_P_decay,
                        ", nozzle budget = ", P_nozzle);
    }

    // C1: per-column pin used by this step: pin_T_eff[col] = pin_T[col] unless
    //   the column never fired (pin_T == 0), the flame cannot heat it past
    //   the pin (q_pin = h(T_gas - pin_T) <= 0), or it has been idle for
    //   time - pin_t >= pinned_idle_cycles * t_cell, with the cell heating
    //   time t_cell = pin_rhoCp*dz*(pin_T - T_amb)/q_pin.
    // Removal (which updates pin_T/pin_t) runs before the thermal advance in
    // Advance(), so a column that fired this step has time - pin_t = 0 and
    // pinned_idle_cycles = 0 never pins (face-identical test hook).
    void BuildPinEffective(int lev, Set::Scalar time)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const std::size_t ncols = static_cast<std::size_t>(domain.length(0)) * domain.length(1);
        if (pin_T.size() != ncols)
        {
            pin_T.assign(ncols, 0.0);
            pin_t.assign(ncols, 0.0);
            pin_rhoCp.assign(ncols, 0.0);
        }
        pin_T_eff.assign(ncols, 0.0);
        pin_reason.assign(ncols, 0);
        const Set::Scalar dz = geom[lev].CellSize()[2];
        for (std::size_t c = 0; c < ncols; ++c)
            // D2a: this column's flame (h, T_gas), from BuildFlameColumns.
            pin_T_eff[c] = PinRule(c, flame_h[c], flame_Tg[c], time, dz, pin_reason[c]);
    }

    // C1 pin rule for one column at flame (h, T_gas); returns the pin (0 = face
    // form) and sets reason (see pin_reason). Shared by BuildPinEffective and
    // the D2a2 JetEnthalpyMarch. pin_T/pin_t/pin_rhoCp must be sized.
    Set::Scalar PinRule(std::size_t c, Set::Scalar h, Set::Scalar T_gas,
                        Set::Scalar time, Set::Scalar dz, int& reason) const
    {
        reason = 0;
        if (pin_T[c] <= 0.0) return 0.0;                                // reason 0: unset
        const Set::Scalar q_pin = h * (T_gas - pin_T[c]);
        if (q_pin <= 0.0) { reason = 2; return 0.0; }                   // q_pin <= 0
        const Set::Scalar t_cell = pin_rhoCp[c] * dz * (pin_T[c] - T_amb) / q_pin;
        if (time - pin_t[c] >= surface_patch.pinned_idle_cycles * t_cell)
        { reason = 1; return 0.0; }                                     // idle
        reason = 3;
        return pin_T[c];
    }

    // C1: thermo diagnostics for the pinned closure. Re-evaluates
    // SurfaceCellFlux on the in-patch top cells with the H update's inputs
    // (T_old, k_c, flame_h/Tg, pin_T_eff), so the branch taken matches the H
    // update. D2a: face-branch reason counters (pin_reason + branch flag) and,
    // with expressions, patch_P_robin = sum q_robin*dx*dy (any robin_form).
    void PatchDiagnostics(int lev, const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const int ncx = domain.length(0);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;
        const bool pinned = surface_patch.robin_form == SurfacePatch::RobinForm::Pinned;
        const bool face   = surface_patch.robin_form == SurfacePatch::RobinForm::Face || pinned;
        Set::Scalar n_pin = 0.0, n_face = 0.0, q_sum = 0.0;
        Set::Scalar n_unset = 0.0, n_idle = 0.0, n_qneg = 0.0, n_min = 0.0;
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp    = (*kappa_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& Df    = (*D_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const int col = (j - jlo_dom) * ncx + (i - ilo_dom);
                const int k = k_top_col[col];
                if (k < lo.z || k > hi.z) continue;   // top cell not owned here
                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (dxp * dxp + dyp * dyp > r2max) continue;
                const Set::Scalar Tc = T_old(i, j, k);
                const Set::Scalar kc = kp(i, j, k) * std::exp(-kappa_damage_alpha * Df(i, j, k));
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                int branch = 0;
                SurfaceCellFlux(true, face, 1.0, flame_h[col], flame_Tg[col],
                                kc, DX[2], Tc, eps_rad, eps_rad_hi, T_eps_break, sigma_SB,
                                h_conv, T_amb, q_in, q_out, q_rob, q_los,
                                pinned ? pin_T_eff[col] : Set::Scalar(0.0), &branch);
                if (branch == 1) n_pin += 1.0;
                else
                {
                    n_face += 1.0;
                    if (pinned)
                    {
                        switch (pin_reason[col])
                        {
                        case 0:  n_unset += 1.0; break;
                        case 1:  n_idle  += 1.0; break;
                        case 2:  n_qneg  += 1.0; break;
                        default: n_min   += 1.0; break;
                        }
                    }
                }
                q_sum += q_rob;
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(n_pin);
        amrex::ParallelDescriptor::ReduceRealSum(n_face);
        amrex::ParallelDescriptor::ReduceRealSum(q_sum);
        if (pinned)
        {
            amrex::ParallelDescriptor::ReduceRealSum(n_unset);
            amrex::ParallelDescriptor::ReduceRealSum(n_idle);
            amrex::ParallelDescriptor::ReduceRealSum(n_qneg);
            amrex::ParallelDescriptor::ReduceRealSum(n_min);
            pinned_cols_pinned = n_pin;
            pinned_cols_face   = n_face;
            pinned_q_mean      = (n_pin + n_face > 0.0) ? q_sum / (n_pin + n_face) : 0.0;
            pinned_cols_face_unset   = n_unset;
            pinned_cols_face_idle    = n_idle;
            pinned_cols_face_qneg    = n_qneg;
            pinned_cols_face_minrule = n_min;
        }
        patch_P_robin = q_sum * DX[0] * DX[1];
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
    // C1 pinned form (face_form = true and T_pin > 0, robin_form = pinned):
    // the same face Newton runs first; T_s = min(T_f, T_pin). If T_f <= T_pin
    // (weak flame, or T_pin = 0: column never fired / idle) the face-form
    // outputs are returned unchanged (*branch = 0). Otherwise the surface is
    // pinned at the column's last firing temperature (*branch = 1): the flame
    // flux AND the losses are both evaluated at T_s = T_pin,
    //   q_robin = h(T_gas - T_s),  q_loss = eps(T_s)*sigma*(T_s^4 - T_a^4) + h_c(T_s - T_a),
    // and the cell absorbs the net q_in = q_robin - q_loss (q_out = 0). The
    // conduction consistency G(T_s - T_c) = q_in is deliberately NOT enforced:
    // that is the closure (the resolved surface rides at T_fire while the
    // coarse cell average lags below it). Non-pinned callers pass T_pin = 0,
    // so the cell and face arithmetic are untouched.
    AMREX_GPU_HOST_DEVICE
    static void SurfaceCellFlux(bool in_patch, bool face_form, Set::Scalar surf_c,
                                Set::Scalar h, Set::Scalar T_gas,
                                Set::Scalar kc, Set::Scalar dz, Set::Scalar Tc,
                                Set::Scalar eps_lo, Set::Scalar eps_hi, Set::Scalar T_eps_b,
                                Set::Scalar sig_SB, Set::Scalar h_c, Set::Scalar T_a,
                                Set::Scalar& q_in, Set::Scalar& q_out,
                                Set::Scalar& q_robin, Set::Scalar& q_loss,
                                Set::Scalar T_pin = 0.0, int* branch = nullptr)
    {
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
            if (!(T_pin > 0.0 && Tf > T_pin))
            {
                q_in    = G * (Tf - Tc);
                q_out   = 0.0;
                q_robin = h * (T_gas - Tf);
                q_loss  = q_robin - q_in;
                if (branch) *branch = 0;
                return;
            }
            // C1 pinned branch: T_s = min(T_f, T_pin) = T_pin.
            const Set::Scalar Ts   = T_pin;
            const Set::Scalar eps_s = (Ts <= T_eps_b) ? eps_lo : eps_hi;
            q_robin = h * (T_gas - Ts);
            q_loss  = eps_s * sig_SB * (Ts * Ts * Ts * Ts - T_a * T_a * T_a * T_a)
                    + h_c * (Ts - T_a);
            q_in    = q_robin - q_loss;
            q_out   = 0.0;
            if (branch) *branch = 1;
            return;
        }
        const Set::Scalar eps = (Tc <= T_eps_b) ? eps_lo : eps_hi;
        q_in    = in_patch ? h * (T_gas - Tc) : Set::Scalar(0.0);
        q_out   = surf_c * (eps * sig_SB * (Tc*Tc*Tc*Tc - T_a*T_a*T_a*T_a)
                            + h_c * (Tc - T_a));
        q_robin = q_in;
        q_loss  = q_out;
    }

    // S1: whole-domain energy ledger at level 0, called after the H update.
    //   E_H   = sum_all H*V (void cells included; their H is frozen)
    //   E_src += dt*(P_beam + P_robin - P_loss)
    //   ledger_err = (E_H - E_H0 - E_src)/max(|E_src|, 1e-30)
    // P_beam reuses the E1 BeamEnergyBalance sum (same BeamSource calls);
    // P_robin/P_loss re-evaluate SurfaceCellFlux with the H update's inputs
    // (prescribed_T patch flux, if any, is counted in P_robin).
    // EXACT only with Neumann-0 outer temperature BCs on every face: a
    // Dirichlet/Robin domain BC puts a flux through the boundary faces that
    // the ledger does not count. The ledger state (ledger_started, E_H0,
    // E_src) is NOT checkpointed: a restarted run re-baselines at its first
    // step. Likewise the C1 pin state (pin_T/pin_t/pin_rhoCp) is not
    // checkpointed: after a restart every column uses the face form until it
    // fires again. E_H is a compensated (Neumaier) sum so the
    // difference E_H - E_H0 keeps the small per-step increments.
    // surface_missing_cols counts columns that have a mask-top solid cell but
    // no surface_mf cell (the phi-window dropout at phi_top ~ dz ties).
    void EnergyLedger(int lev, Set::Scalar time, Set::Scalar dt,
                      const amrex::Vector<int>& k_top_col)
    {
        if (lev != 0) return;
        const bool void_aware  = removal_enabled();
        const bool follow_mask = surface_follow_mask;
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const int khi_dom = domain.bigEnd(2);
        const int ncx = domain.length(0), ncy = domain.length(1);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz = DX[2];
        const Set::Scalar A  = DX[0] * DX[1];
        const Set::Scalar V  = A * dz;
        const Set::Scalar inv_dz2 = 1.0 / (dz * dz);
        const Set::Scalar a_D = kappa_damage_alpha;

        const bool        sp_on    = surface_patch.enabled;
        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
        const bool        sp_pinned = sp_on && surface_patch.robin_form ==
                                      SurfacePatch::RobinForm::Pinned && lev == 0;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face || sp_pinned;
        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_T     = (sp_on && sp_presc) ? surface_patch.T_at(time) : 0.0;

        auto neumaier = [](Set::Scalar& s, Set::Scalar& c, Set::Scalar x) {
            const Set::Scalar t = s + x;
            if (std::abs(s) >= std::abs(x)) c += (s - t) + x;
            else                            c += (x - t) + s;
            s = t;
        };

        Set::Scalar EH_s = 0.0, EH_c = 0.0, EH0_s = 0.0, EH0_c = 0.0;
        Set::Scalar P_robin = 0.0, P_loss = 0.0;
        std::vector<int> has_surf(static_cast<std::size_t>(ncx) * ncy, 0);
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& H     = (*H_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& H_old = (*H_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp    = (*kappa_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& Df    = (*D_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& surf  = (*surface_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& rem   =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                neumaier(EH_s, EH_c, H(i, j, k) * V);
                if (!ledger_started) neumaier(EH0_s, EH0_c, H_old(i, j, k) * V);
                const int col = (j - jlo_dom) * ncx + (i - ilo_dom);
                if (surf(i, j, k) > 0.5) has_surf[col] = 1;
                const bool self_void = void_aware && rem(i, j, k) > 0.5;
                if (self_void) continue;

                const Set::Scalar Tc = T_old(i, j, k);
                const Set::Scalar kc = kp(i, j, k) * std::exp(-a_D * Df(i, j, k));
                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const int k_top_c = (follow_mask && void_aware) ? k_top_col[col] : khi_dom;
                const Set::Scalar surf_c = follow_mask
                                         ? Set::Scalar(k == k_top_c ? 1.0 : 0.0)
                                         : surf(i, j, k);
                bool in_patch = false;
                if (sp_on && k == k_top_c)
                {
                    const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                    if (dxp * dxp + dyp * dyp <= sp_r2)
                    {
                        if (sp_presc) P_robin += kc * (sp_T - Tc) * inv_dz2 * V;
                        else          in_patch = true;
                    }
                }
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                SurfaceCellFlux(in_patch, sp_face, surf_c,
                                sp_on ? flame_h[col]  : Set::Scalar(0.0),
                                sp_on ? flame_Tg[col] : Set::Scalar(0.0), kc, dz, Tc,
                                eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                                q_in, q_out, q_rob, q_los,
                                sp_pinned ? pin_T_eff[col] : Set::Scalar(0.0));
                P_robin += q_rob * A;
                P_loss  += q_los * A;
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(EH_s);
        amrex::ParallelDescriptor::ReduceRealSum(EH_c);
        amrex::ParallelDescriptor::ReduceRealSum(P_robin);
        amrex::ParallelDescriptor::ReduceRealSum(P_loss);
        amrex::ParallelDescriptor::ReduceIntMax(has_surf.data(), static_cast<int>(has_surf.size()));
        if (!ledger_started)
        {
            amrex::ParallelDescriptor::ReduceRealSum(EH0_s);
            amrex::ParallelDescriptor::ReduceRealSum(EH0_c);
            ledger_E_H0_s  = EH0_s;
            ledger_E_H0_c  = EH0_c;
            ledger_started = true;
        }

        int n_missing = 0;
        for (std::size_t c = 0; c < has_surf.size(); ++c)
        {
            const bool has_solid = void_aware ? (k_top_col[c] != std::numeric_limits<int>::min())
                                              : true;
            if (has_solid && !has_surf[c]) ++n_missing;
        }

        const Set::Scalar P_beam = (beam.power_at(time) > 0.0) ? beam_P_dep : 0.0;
        ledger_P_beam  = P_beam;
        ledger_P_robin = P_robin;
        ledger_P_loss  = P_loss;
        ledger_E_src  += dt * (P_beam + P_robin - P_loss);
        ledger_E_H     = EH_s + EH_c;
        const Set::Scalar dE = (EH_s - ledger_E_H0_s) + (EH_c - ledger_E_H0_c);
        ledger_err = (dE - ledger_E_src) / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));
        ledger_surface_missing_cols = Set::Scalar(n_missing);
    }

    // E1: optical path from the top face of the first solid cell to the
    // centre of a cell with n_above solid cells above it in its column (see
    // the BeamSource call in AdvanceMicrostructure). With the removal
    // invariant n_above == floor(phi/dz), but the packet formula
    // (floor(phi/dz) + 1/2)*dz is NOT evaluated from phi: at a round-off tie
    // the column shift leaves e.g. phi = -5e-18 (voided) above phi = dz
    // exactly (top solid), floor gives 1, and the whole top slice (98% of the
    // column power at alpha*dz = 4) is lost -- seen on the Meier dev harness,
    // where it stalled one column for 22 s. Counting solid cells from the
    // `removed` mask is exact by construction. Do NOT use
    // floor((phi - dz/2)/dz) either: for phi_top in [0, dz/2) the top and
    // second cells both get path dz/2 and the surface slice is deposited twice.
    AMREX_GPU_HOST_DEVICE
    static Set::Scalar SolidFacePath(int n_above, Set::Scalar dz)
    {
        return (Set::Scalar(n_above) + 0.5) * dz;
    }

    // E1: global topmost solid (removed <= 0.5) cell index per (i, j)
    // column at level lev, INT_MIN for fully removed columns. Index is
    // (j - dlo_y)*nx + (i - dlo_x). MPI-reduced, so identical on every rank
    // and independent of the z box decomposition (cf. Removal.H PASS 1,
    // which additionally requires phi >= 0).
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

    // E1: per-step closure of the microstructure-path beam deposition.
    //
    //   P_inc   = sum over lit columns of (1-R)*Pi*dx*dy at the top solid cell
    //   P_dep   = sum over solid cells of Q*dx*dy*dz  (same BeamSource call)
    //   P_trans = (1-R)*Pi*exp(-a*s_bot)*dx*dy at k == domain bottom
    //   closure : P_dep + P_trans == P_inc
    //
    // Exact (to round-off) only for CONSTANT alpha and the UNIFORM profile.
    // With alpha(T) each cell uses its own a_k, the column sum no longer
    // telescopes and P_dep can exceed P_inc; a correct alpha(T) deposition
    // needs the cumulative optical depth tau = sum_k a_k*dz (follow-on step,
    // required before real alpha(T) MMW work). A Gaussian profile adds an
    // O((dz/z_R)^2) deviation because Pi evolves with the geometric depth, so
    // the P_dep <= P_inc double-count guard is applied to the uniform profile
    // only.
    //
    // Invariant (void-aware runs), per solid cell with n = k_top - k:
    //   (a) contiguity: "cell above is void or k is the domain top" <=> n == 0
    //       (a mismatch is a hole or floating solid in the column);
    //   (b) phi/mask sync: phi - n*dz in [0, dz), with a 1e-6*dz round-off
    //       band, i.e. floor(phi/dz) == n up to ties. A violation means phi
    //       and the `removed` mask disagree by more than round-off.
    // Single level: lev == 0 only (multi-level unsupported; composite sums
    // would double-count covered cells).
    void BeamEnergyBalance(int lev, Set::Scalar time, Set::Scalar dt,
                           const amrex::Vector<int>& k_top_col)
    {
        if (lev != 0) return;
        const bool void_aware = removal_enabled();
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
        auto alpha_f = beam.alpha_f;

        Set::Scalar P_inc = 0.0, P_dep = 0.0, P_trans = 0.0, n_viol = 0.0;
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).array(mfi);
            // Placeholders as in AdvanceMicrostructure when not void-aware.
            amrex::Array4<const Set::Scalar> const& rem =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& phi =
                void_aware ? (*phi_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (void_aware && rem(i, j, k) > 0.5) continue;
                const int n_above = void_aware
                    ? k_top_col[(j - jlo_dom) * ncx + (i - ilo_dom)] - k : 0;
                const Set::Scalar path = void_aware ? SolidFacePath(n_above, dz)
                                                    : Set::Scalar(-1.0);
                // removed_mf ghosts were filled at the top of
                // AdvanceMicrostructure, so rem(k+1) is valid at box seams.
                const bool is_top_solid =
                    (k == khi_dom) || (void_aware && rem(i, j, k + 1) > 0.5);
                if (void_aware)
                {
                    if (is_top_solid != (n_above == 0)) n_viol += 1.0;        // (a)
                    const Set::Scalar d = phi(i, j, k) - Set::Scalar(n_above) * dz;
                    if (d < -1.0e-6 * dz || d >= (1.0 + 1.0e-6) * dz) n_viol += 1.0; // (b)
                }

                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                const Set::Scalar Pi = BeamIncidentIntensity(P_now, beam.omega0, beam.lambda,
                                                             beam.x0, beam.y0, beam.z0,
                                                             beam.profile, beam.radius, x, y, z);
                if (Pi == 0.0) continue;

                const Set::Scalar Tc = T_old(i, j, k);
                const Set::Scalar Q = BeamSource(P_now, beam.omega0, beam.lambda,
                                                 beam.x0, beam.y0, beam.z0, beam.R,
                                                 beam.profile, beam.radius, alpha_f,
                                                 x, y, z, Tc, dz, path);
                P_dep += Q * V;
                if (is_top_solid) P_inc += (1.0 - beam.R) * Pi * A;
                if (k == klo_dom)
                {
                    const Set::Scalar a = alpha_f(Tc);
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
    // with n_above = solid cells above in the column (== floor(phi/dz), since
    // removal keeps the top solid cell at phi in [0, dz)).
    // Raw phi leaks the slice between the level set and that face.
    // dz_abs < 0 (the default) falls back to the geometric depth z_0 - z,
    // which is correct whenever no material above this cell has been removed
    // (then z_surface == z_0 and phi == z_0 - z exactly). The Gaussian waist
    // still evolves with the geometric distance z_0 - z regardless, since the
    // beam propagates geometrically through the void cavity.
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
                                  amrex::ParserExecutor<1> alpha_f,
                                  Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                  Set::Scalar T,
                                  Set::Scalar dz_cell = 0.0,
                                  Set::Scalar dz_abs = -1.0)
    {
        Set::Scalar dz = z0 - z;
        Set::Scalar Pi = BeamIncidentIntensity(P, w0, lam, x0, y0, z0,
                                               profile, radius, x, y, z);
        if (Pi == 0.0) return 0.0;
        Set::Scalar a    = alpha_f(T);
        Set::Scalar path = (dz_abs >= 0.0) ? dz_abs : dz;
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
    // uniform disk. The Gaussian waist evolves with the geometric z0 - z.
    AMREX_GPU_HOST_DEVICE
    static Set::Scalar BeamIncidentIntensity(Set::Scalar P, Set::Scalar w0, Set::Scalar lam,
                                             Set::Scalar x0, Set::Scalar y0, Set::Scalar z0,
                                             int profile, Set::Scalar radius,
                                             Set::Scalar x, Set::Scalar y, Set::Scalar z)
    {
        Set::Scalar dz = z0 - z;
        if (dz < 0.0 || P == 0.0) return 0.0;
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
            Set::Scalar w  = w0 * std::sqrt(1.0 + (dz / zR) * (dz / zR));
            Pi = (2.0 * P / (M_PI * w0 * w0)) * (w0 / w) * (w0 / w)
               * std::exp(-2.0 * r2 / (w * w));
        }
        return Pi;
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow) override
    {
        // Reuses HeatConduction's tagging (∇T) and Mechanics' tagging.
        HeatConduction::TagCellsForRefinement(lev, a_tags, time, ngrow);
        Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev, a_tags, time, ngrow);
    }

private:
    static void ParseDamageSettings(MMWSpalling& value, IO::ParmParse& pp)
    {
        int damage_on = 0;
        pp.query_default("damage.enabled", damage_on, 0);
        value.damage_enabled = (damage_on != 0);

        if (value.damage_enabled)
        {
            pp.query_required("damage.phi_deg",  value.damage_phi_deg);
            pp.query_required("damage.cohesion", value.damage_cohesion);
        }
        else
        {
            pp.query_default("damage.phi_deg",  value.damage_phi_deg,  0.0);
            pp.query_default("damage.cohesion", value.damage_cohesion, 1.0);
        }
        pp.query_default("damage.psi0",            value.damage_psi0,            0.0);
        pp.query_default("damage.n",               value.damage_n,               1.0);
        pp.query_default("damage.m",               value.damage_m,               1.0);
        pp.query_default("damage.gb_multiplier",   value.damage_gb_multiplier,   1.0);
        pp.query_default("damage.stiffness_floor", value.damage_stiffness_floor, 1.0e-6);

        int prescribed_on = 0;
        pp.query_default("damage.prescribed_stress.enabled", prescribed_on, 0);
        value.damage_prescribed_stress_enabled = (prescribed_on != 0);
        if (value.damage_prescribed_stress_enabled)
        {
            auto compile = [&pp](const std::string& key,
                                 amrex::Parser& parser,
                                 amrex::ParserExecutor<4>& exec,
                                 const std::string& def) {
                std::string expr;
                pp.query_default(key, expr, def);
                parser = amrex::Parser(expr);
                parser.registerVariables({"x", "y", "z", "t"});
                exec = parser.compile<4>();
            };
            compile("damage.prescribed_stress.xx_expr",
                    value.damage_sxx_parser, value.damage_sxx_f, "0.0");
            compile("damage.prescribed_stress.yy_expr",
                    value.damage_syy_parser, value.damage_syy_f, "0.0");
            compile("damage.prescribed_stress.zz_expr",
                    value.damage_szz_parser, value.damage_szz_f, "0.0");
            compile("damage.prescribed_stress.xy_expr",
                    value.damage_sxy_parser, value.damage_sxy_f, "0.0");
            compile("damage.prescribed_stress.yz_expr",
                    value.damage_syz_parser, value.damage_syz_f, "0.0");
            compile("damage.prescribed_stress.zx_expr",
                    value.damage_szx_parser, value.damage_szx_f, "0.0");
        }

        if (!value.damage_enabled) return;

        using MechanicsBase = Base::Mechanics<Model::Solid::Affine::Isotropic>;
        if (!value.damage_prescribed_stress_enabled &&
            value.m_type == MechanicsBase::Type::Disable)
            Util::Abort(INFO, "damage.enabled requires active mechanics unless "
                              "damage.prescribed_stress.enabled = 1");

        if (!(std::isfinite(value.damage_phi_deg) &&
              value.damage_phi_deg >= 0.0 && value.damage_phi_deg < 90.0))
            Util::Abort(INFO, "damage.phi_deg must be finite and in [0, 90), got ",
                        value.damage_phi_deg);
        if (!(std::isfinite(value.damage_cohesion) && value.damage_cohesion > 0.0))
            Util::Abort(INFO, "damage.cohesion must be finite and > 0, got ",
                        value.damage_cohesion);
        if (!(std::isfinite(value.damage_psi0) && value.damage_psi0 >= 0.0))
            Util::Abort(INFO, "damage.psi0 must be finite and >= 0, got ",
                        value.damage_psi0);
        if (!(std::isfinite(value.damage_n) && value.damage_n >= 0.0))
            Util::Abort(INFO, "damage.n must be finite and >= 0, got ",
                        value.damage_n);
        if (!(std::isfinite(value.damage_m) && value.damage_m >= 0.0))
            Util::Abort(INFO, "damage.m must be finite and >= 0, got ",
                        value.damage_m);
        if (!(std::isfinite(value.damage_gb_multiplier) &&
              value.damage_gb_multiplier >= 0.0))
            Util::Abort(INFO, "damage.gb_multiplier must be finite and >= 0, got ",
                        value.damage_gb_multiplier);
        if (!(std::isfinite(value.damage_stiffness_floor) &&
              value.damage_stiffness_floor > 0.0 &&
              value.damage_stiffness_floor <= 1.0))
            Util::Abort(INFO, "damage.stiffness_floor must be finite and in (0, 1], got ",
                        value.damage_stiffness_floor);

        constexpr Set::Scalar pi = 3.141592653589793238462643383279502884;
        value.damage_phi_rad = value.damage_phi_deg * pi / 180.0;
        value.damage_alpha_dp =
            Numeric::DruckerPrager::AlphaFromPhiRadians(value.damage_phi_rad);
        value.damage_k_dp =
            Numeric::DruckerPrager::KFromCohesionPhiRadians(value.damage_cohesion,
                                                            value.damage_phi_rad);
        if (!(std::isfinite(value.damage_alpha_dp) &&
              std::isfinite(value.damage_k_dp) &&
              value.damage_k_dp > 0.0))
            Util::Abort(INFO, "invalid Drucker-Prager parameters: alpha=",
                        value.damage_alpha_dp, ", k=", value.damage_k_dp);
    }

    void RefreshKappaEff(int lev)
    {
        const Set::Scalar a_D = kappa_damage_alpha;
        for (amrex::MFIter mfi(*kappa_eff_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& kp = (*kappa_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& Df = (*D_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& ke = (*kappa_eff_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Dv = Df(i, j, k);
                if (!std::isfinite(Dv)) Dv = 0.0;
                Dv = std::max(Set::Scalar(0.0), std::min(Set::Scalar(1.0), Dv));
                ke(i, j, k) = kp(i, j, k) * std::exp(-a_D * Dv);
            });
        }
    }

    static Set::Scalar DamageGrowthRHS(Set::Scalar D,
                                       Set::Scalar psi,
                                       Set::Scalar A_D,
                                       Set::Scalar psi0,
                                       Set::Scalar n,
                                       Set::Scalar m)
    {
        if (psi < psi0 || A_D <= 0.0) return 0.0;
        Set::Scalar intact = std::max(Set::Scalar(0.0), 1.0 - D);
        Set::Scalar drive = psi / (1.0 + psi0);
        if (drive < 0.0) return 0.0;
        return A_D * std::pow(drive, n) * std::pow(intact, m);
    }

    static Set::Scalar IntegrateDamageRK4(Set::Scalar D,
                                          Set::Scalar dt,
                                          Set::Scalar psi,
                                          Set::Scalar A_D,
                                          Set::Scalar psi0,
                                          Set::Scalar n,
                                          Set::Scalar m)
    {
        auto clamp01 = [](Set::Scalar x) {
            return std::max(Set::Scalar(0.0), std::min(Set::Scalar(1.0), x));
        };
        Set::Scalar k1 = DamageGrowthRHS(clamp01(D), psi, A_D, psi0, n, m);
        Set::Scalar k2 = DamageGrowthRHS(clamp01(D + 0.5 * dt * k1), psi, A_D, psi0, n, m);
        Set::Scalar k3 = DamageGrowthRHS(clamp01(D + 0.5 * dt * k2), psi, A_D, psi0, n, m);
        Set::Scalar k4 = DamageGrowthRHS(clamp01(D +       dt * k3), psi, A_D, psi0, n, m);
        return D + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    void UpdateDamageAfterMechanics(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (!damage_enabled) return;
        Util::Assert(INFO, TEST(microstructure_enabled));

        Util::RealFillBoundary(*D_mf[lev], geom[lev]);
        Util::RealFillBoundary(*phase_mf[lev], geom[lev]);
        Util::RealFillBoundary(*is_grain_boundary_mf[lev], geom[lev]);
        if (!damage_prescribed_stress_enabled)
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        (void)DX;
        const int n_phases = static_cast<int>(phases.size());
        const Set::Scalar alpha_dp = damage_alpha_dp;
        const Set::Scalar k_dp     = damage_k_dp;
        const Set::Scalar psi0     = damage_psi0;
        const Set::Scalar n_exp    = damage_n;
        const Set::Scalar m_exp    = damage_m;
        const Set::Scalar gb_mult  = damage_gb_multiplier;
        const bool prescribed = damage_prescribed_stress_enabled;
        auto sxx_f = damage_sxx_f;
        auto syy_f = damage_syy_f;
        auto szz_f = damage_szz_f;
        auto sxy_f = damage_sxy_f;
        auto syz_f = damage_syz_f;
        auto szx_f = damage_szx_f;

        for (amrex::MFIter mfi(*D_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& D = (*D_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gb = (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Fout = (*damage_F_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Psiout = (*damage_Psi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rateout = (*damage_rate_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> sig_node;
            if (!prescribed) sig_node = (*stress_mf[lev]).array(mfi);

            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                Set::Matrix sigma = Set::Matrix::Zero();
                if (prescribed)
                {
                    Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                    sigma(0, 0) = sxx_f(x, y, z, time);
#if AMREX_SPACEDIM > 1
                    sigma(1, 1) = syy_f(x, y, z, time);
                    Set::Scalar sxy = sxy_f(x, y, z, time);
                    sigma(0, 1) = sxy;
                    sigma(1, 0) = sxy;
#endif
#if AMREX_SPACEDIM > 2
                    sigma(2, 2) = szz_f(x, y, z, time);
                    Set::Scalar syz = syz_f(x, y, z, time);
                    Set::Scalar szx = szx_f(x, y, z, time);
                    sigma(1, 2) = syz;
                    sigma(2, 1) = syz;
                    sigma(2, 0) = szx;
                    sigma(0, 2) = szx;
#endif
                }
                else
                {
                    sigma = Numeric::Interpolate::NodeToCellAverage(sig_node, i, j, k, 0);
                }

                if (!FiniteMatrix(sigma))
                    Util::Abort(INFO, "non-finite averaged stress before Drucker-Prager at cell (",
                                i, ",", j, ",", k, ")");

                Set::Scalar Fv = Numeric::DruckerPrager::YieldFunction(sigma, alpha_dp, k_dp);
                Set::Scalar psiv = Fv / k_dp;
                if (!std::isfinite(Fv) || !std::isfinite(psiv))
                    Util::Abort(INFO, "non-finite Drucker-Prager value at cell (",
                                i, ",", j, ",", k, "): F=", Fv, ", Psi=", psiv);

                int p = static_cast<int>(std::round(ph(i, j, k)));
                if (p < 0 || p >= n_phases) p = 0;
                Set::Scalar A_D = phases[p].A_D;
                if (gb(i, j, k) > 0.5) A_D *= gb_mult;

                Set::Scalar D0 = D(i, j, k);
                if (!std::isfinite(D0))
                    Util::Abort(INFO, "non-finite damage D at cell (",
                                i, ",", j, ",", k, "): D=", D0);
                D0 = std::max(Set::Scalar(0.0), std::min(Set::Scalar(1.0), D0));

                Set::Scalar rate = DamageGrowthRHS(D0, psiv, A_D, psi0, n_exp, m_exp);
                Set::Scalar Dcand = IntegrateDamageRK4(D0, dt, psiv, A_D, psi0, n_exp, m_exp);
                if (!std::isfinite(rate) || !std::isfinite(Dcand))
                    Util::Abort(INFO, "non-finite damage update at cell (",
                                i, ",", j, ",", k, "): rate=", rate,
                                ", candidate=", Dcand);
                Set::Scalar Dnew = std::max(D0,
                    std::max(Set::Scalar(0.0), std::min(Set::Scalar(1.0), Dcand)));

                D(i, j, k) = Dnew;
                Fout(i, j, k) = Fv;
                Psiout(i, j, k) = psiv;
                rateout(i, j, k) = rate;
            });
        }

        Util::RealFillBoundary(*D_mf[lev], geom[lev]);
        RefreshKappaEff(lev);
    }

    static void ParseCohesiveSettings(MMWSpalling& value, IO::ParmParse& pp)
    {
        int on = 0;
        pp.query_default("cohesive.enabled", on, 0);
        value.cohesive_enabled = (on != 0);
        if (!value.cohesive_enabled) return;

        pp.query_required("cohesive.K_n", value.cohesive_K_n);
        if (!(std::isfinite(value.cohesive_K_n) && value.cohesive_K_n > 0.0))
            Util::Abort(INFO, "cohesive.K_n must be finite and > 0, got ", value.cohesive_K_n);

        // Calibration: prefer explicit (delta_c, delta_max) when provided;
        // otherwise derive from (G_c, f_t) using the bilinear identities
        //   delta_c   = f_t / K_n
        //   delta_max = 2 * G_c / f_t
        Set::Scalar delta_c = 0.0, delta_max = 0.0;
        bool have_delta_c   = pp.contains("cohesive.delta_c");
        bool have_delta_max = pp.contains("cohesive.delta_max");
        if (have_delta_c)   pp.query_required("cohesive.delta_c",   delta_c);
        if (have_delta_max) pp.query_required("cohesive.delta_max", delta_max);
        if (!have_delta_c || !have_delta_max)
        {
            Set::Scalar G_c = 0.0, f_t = 0.0;
            pp.query_required("cohesive.G_c", G_c);
            pp.query_required("cohesive.f_t", f_t);
            if (!(std::isfinite(G_c) && G_c > 0.0))
                Util::Abort(INFO, "cohesive.G_c must be finite and > 0, got ", G_c);
            if (!(std::isfinite(f_t) && f_t > 0.0))
                Util::Abort(INFO, "cohesive.f_t must be finite and > 0, got ", f_t);
            if (!have_delta_c)   delta_c   = f_t / value.cohesive_K_n;
            if (!have_delta_max) delta_max = 2.0 * G_c / f_t;
        }
        if (!(std::isfinite(delta_c) && delta_c > 0.0))
            Util::Abort(INFO, "cohesive delta_c must be finite and > 0, got ", delta_c);
        if (!(std::isfinite(delta_max) && delta_max > delta_c))
            Util::Abort(INFO, "cohesive delta_max must be finite and > delta_c, got delta_c=",
                        delta_c, ", delta_max=", delta_max);
        value.cohesive_delta_c   = delta_c;
        value.cohesive_delta_max = delta_max;

        int presc_on = 0;
        pp.query_default("cohesive.prescribed_delta.enabled", presc_on, 0);
        value.cohesive_prescribed_enabled = (presc_on != 0);
        if (value.cohesive_prescribed_enabled)
        {
            std::string expr;
            pp.query_required("cohesive.prescribed_delta.expr", expr);
            value.cohesive_delta_parser = amrex::Parser(expr);
            value.cohesive_delta_parser.registerVariables({"x", "y", "z", "t"});
            value.cohesive_delta_f = value.cohesive_delta_parser.compile<4>();
        }
    }

    void UpdateCohesiveAfterMechanics(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (!cohesive_enabled) return;
        Util::Assert(INFO, TEST(microstructure_enabled));
        (void)dt;

        Util::RealFillBoundary(*is_grain_boundary_mf[lev], geom[lev]);

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar K_n       = cohesive_K_n;
        const Set::Scalar delta_c   = cohesive_delta_c;
        const Set::Scalar delta_max = cohesive_delta_max;
        const bool prescribed = cohesive_prescribed_enabled;
        auto delta_f = cohesive_delta_f;

        for (amrex::MFIter mfi(*gb_delta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& gb = (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& d_now    = (*gb_delta_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& d_hist   = (*gb_delta_max_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& trac     = (*gb_traction_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& dam      = (*gb_cohesive_damage_mf[lev]).array(mfi);

            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (!(gb(i, j, k) > 0.5))
                {
                    d_now (i, j, k) = 0.0;
                    d_hist(i, j, k) = 0.0;
                    trac  (i, j, k) = 0.0;
                    dam   (i, j, k) = 0.0;
                    return;
                }

                Set::Scalar delta_input = 0.0;
                if (prescribed)
                {
                    Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                    delta_input = delta_f(x, y, z, time);
                    if (!std::isfinite(delta_input)) delta_input = 0.0;
                }
                if (delta_input < 0.0) delta_input = 0.0;

                Set::Scalar prev_hist = d_hist(i, j, k);
                if (!std::isfinite(prev_hist) || prev_hist < 0.0) prev_hist = 0.0;
                Set::Scalar new_hist = std::max(prev_hist, delta_input);

                Set::Scalar t_now = Numeric::CohesiveZone::SecantTraction(
                    delta_input, new_hist, K_n, delta_c, delta_max);
                Set::Scalar D_co  = Numeric::CohesiveZone::EnvelopeDamage(
                    new_hist, K_n, delta_c, delta_max);
                if (!std::isfinite(t_now)) t_now = 0.0;
                if (!std::isfinite(D_co))  D_co  = 0.0;

                d_now (i, j, k) = delta_input;
                d_hist(i, j, k) = new_hist;
                trac  (i, j, k) = t_now;
                dam   (i, j, k) = D_co;
            });
        }
    }

    void ResetCohesiveDiagnostics(int lev, bool reset_history)
    {
        if (!cohesive_enabled) return;
        gb_delta_mf[lev]->setVal(0.0);
        gb_traction_mf[lev]->setVal(0.0);
        gb_cohesive_damage_mf[lev]->setVal(0.0);
        if (reset_history) gb_delta_max_mf[lev]->setVal(0.0);
    }

    static void ParseSpallSettings(MMWSpalling& value, IO::ParmParse& pp)
    {
        int on = 0;
        pp.query_default("spall.enabled", on, 0);
        value.spall_enabled = (on != 0);

        // Stress source for the Sp = K_I/K_Ic criterion. Parsed BEFORE the
        // spall.enabled early-return because it also governs the Sp_field
        // ONSET diagnostic (UpdateSpAfterMechanics), which runs under
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

        // Step 20: surface-normal (slope-aware) spall criterion + removal.
        // Default off (existing behaviour: criterion and recession strictly
        // vertical). When on, each column's local cavity-wall inclination
        // theta is computed from the gathered column tops; the Sp/K_I depth
        // scan is evaluated along the surface NORMAL (the vertical T profile
        // is compressed by cos(theta) -- the physically correct crack frame,
        // Kant's 1-D model applies normal to the surface), and a flake of
        // normal thickness h recedes the column vertically by h/cos(theta).
        // This is the graph form of level-set normal motion
        // (dh/dt = -v_n*sqrt(1+|grad h|^2)) and lets inclined cavity walls
        // keep spalling instead of freezing into a cone (the rim-baking
        // artifact). Parsed before the early-return because it also governs
        // the Sp_field gate diagnostic (UpdateSpAfterMechanics). Under
        // local_thermoelastic the tangential stress magnitude is orientation-
        // independent (isotropic biaxial), so no tensor rotation is needed --
        // only the scan direction changes.
        int sn_on = 0;
        pp.query_default("spall.surface_normal", sn_on, 0);
        value.spall_surface_normal = (sn_on != 0);

        // P1: melt-aware spallation (default off, byte-identical when off).
        // A melt has no load-bearing skeleton, so under sp_weibull the K_I
        // driving stress sigma_xx is scaled by the skeleton function
        // f(Lambda_L) = max(0, 1 - Lambda_L/lambda_crit) (1 in solid, 0 at the
        // rheologically critical melt fraction lambda_crit ~ 0.4, rigidity
        // percolation). The lever is the DRIVING STRESS, not K_Ic. This alone
        // makes the depth scan self-anchor at the solid-liquid interface: the
        // scan breaks at the first K_I < K_Ic, and f -> 0 zeros K_I in melt, so
        // a molten surface cell stops the scan (h_col = 0, no spall = stall)
        // while a melt cell *below* solid stops it there (h_col = the solid
        // thickness). No separate "anchor at the topmost solid cell" pass is
        // needed (an explicit one would evacuate the melt cap with the flake --
        // that is P2's job, and it would destroy the P1 melt-blanket stall).
        // Parsed before the spall.enabled early-return so the keys are always
        // consumed if present; the feature only acts on the sp_weibull removal
        // path.
        int melt_on = 0;
        pp.query_default("spall.melt_aware", melt_on, 0);
        value.melt_aware = (melt_on != 0);
        pp.query_default("melt.lambda_crit", value.melt_lambda_crit, 0.4);
        if (!(std::isfinite(value.melt_lambda_crit) &&
              value.melt_lambda_crit > 0.0 &&
              value.melt_lambda_crit <= 1.0))
            Util::Abort(INFO, "melt.lambda_crit must be in (0, 1], got ",
                        value.melt_lambda_crit);

        // CSV diagnostic toggles (default off). Parsed before the early-return
        // because the clusters CSV is emitted by UpdateSpClusters under
        // sp_weibull even when material removal is off.
        int csv_on = 0;
        pp.query_default("spall.clusters_csv", csv_on, 0);
        value.sp_emit_clusters_csv = (csv_on != 0);
        csv_on = 0;
        pp.query_default("spall.h_col_events_csv", csv_on, 0);
        value.sp_emit_h_col_events_csv = (csv_on != 0);
        // A1: every-removal event log (validated at the end of Parse, once
        // spallation.model is known).
        csv_on = 0;
        pp.query_default("spall.removal_events_csv", csv_on, 0);
        value.sp_emit_removal_events_csv = (csv_on != 0);

        if (!value.spall_enabled) return;

        pp.query_default("spall.damage_threshold", value.spall_damage_threshold, 0.99);
        pp.query_default("spall.xi_s",             value.spall_xi_s,             8.0);
        pp.query_default("spall.C_h",              value.spall_C_h,              1.0);

        if (pp.contains("spall.alpha"))
        {
            value.spall_alpha_override = true;
            pp.query_required("spall.alpha", value.spall_alpha);
        }
        if (pp.contains("spall.t_spall"))
        {
            value.spall_t_spall_override = true;
            pp.query_required("spall.t_spall", value.spall_t_spall);
        }

        if (!(std::isfinite(value.spall_damage_threshold) &&
              value.spall_damage_threshold > 0.0 &&
              value.spall_damage_threshold <= 1.0))
            Util::Abort(INFO, "spall.damage_threshold must be in (0, 1], got ",
                        value.spall_damage_threshold);
        if (!(std::isfinite(value.spall_xi_s) && value.spall_xi_s > 0.0))
            Util::Abort(INFO, "spall.xi_s must be finite and > 0, got ", value.spall_xi_s);
        if (!(std::isfinite(value.spall_C_h) && value.spall_C_h > 0.0))
            Util::Abort(INFO, "spall.C_h must be finite and > 0, got ", value.spall_C_h);
        if (value.spall_alpha_override && !(std::isfinite(value.spall_alpha) &&
                                            value.spall_alpha > 0.0))
            Util::Abort(INFO, "spall.alpha must be finite and > 0, got ", value.spall_alpha);
        if (value.spall_t_spall_override && !(std::isfinite(value.spall_t_spall) &&
                                              value.spall_t_spall > 0.0))
            Util::Abort(INFO, "spall.t_spall must be finite and > 0, got ", value.spall_t_spall);

        int presc = 0;
        pp.query_default("spall.prescribed_ratio.enabled", presc, 0);
        value.spall_prescribed_ratio_enabled = (presc != 0);
        if (value.spall_prescribed_ratio_enabled)
        {
            std::string expr;
            pp.query_required("spall.prescribed_ratio.expr", expr);
            value.spall_ratio_parser = amrex::Parser(expr);
            value.spall_ratio_parser.registerVariables({"x", "y", "z", "t"});
            value.spall_ratio_f = value.spall_ratio_parser.compile<4>();
        }

        // Step 14b: spall D-sampling mode. Default `top_cell` preserves
        // Step 10/11 behaviour byte-identical. `column_max` and
        // `shallow_max` reduce D over a per-column scan window so the
        // pipeline can fire when DP damage concentrates a few cells
        // below the surface (e.g. Hu's prescribed-T patch geometry).
        std::string sample_str;
        pp.query_default("spall.sample", sample_str, std::string("top_cell"));
        if      (sample_str == "top_cell")
            value.spall_sample_mode = SpallSampleMode::TopCell;
        else if (sample_str == "column_max")
            value.spall_sample_mode = SpallSampleMode::ColumnMax;
        else if (sample_str == "shallow_max")
            value.spall_sample_mode = SpallSampleMode::ShallowMax;
        else
            Util::Abort(INFO, "spall.sample must be 'top_cell', 'column_max', "
                              "or 'shallow_max', got: ", sample_str);

        pp.query_default("spall.shallow_depth", value.spall_shallow_depth, 0.006);
        if (!(std::isfinite(value.spall_shallow_depth) &&
              value.spall_shallow_depth > 0.0))
            Util::Abort(INFO, "spall.shallow_depth must be finite and > 0, got ",
                        value.spall_shallow_depth);

        // Step 16: Sp-vs-depth h_spall scan controls. Only consumed when
        // spallation.model = sp_weibull; the damage_law branch ignores
        // both. Defaults: 10 mm safety cap, segments = sp_n_segments.
        pp.query_default("spall.h_spall_max", value.spall_h_spall_max, 0.010);
        if (!(std::isfinite(value.spall_h_spall_max) &&
              value.spall_h_spall_max > 0.0))
            Util::Abort(INFO, "spall.h_spall_max must be finite and > 0, got ",
                        value.spall_h_spall_max);
        pp.query_default("spall.h_spall_n_segments",
                         value.spall_h_spall_n_segments, 0);
        if (value.spall_h_spall_n_segments < 0)
            Util::Abort(INFO, "spall.h_spall_n_segments must be >= 0, got ",
                        value.spall_h_spall_n_segments);

        pp.query_default("spall.flake_coherence_length", value.spall_flake_coh_len, 0.0);
        if (!(std::isfinite(value.spall_flake_coh_len) && value.spall_flake_coh_len >= 0.0))
            Util::Abort(INFO, "spall.flake_coherence_length must be finite and >= 0, got ",
                        value.spall_flake_coh_len);
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

    // Step 19: mechanical bit feed (drawworks). A rigid flat bit face descends
    // along the axis at `bit.feed_rate` and forces every column within
    // `bit.radius` down to the bit floor z_bit(t) = bit.z0 - feed*(t - start).
    // Requires the spall/removal machinery (phi/removed); abort if neither
    // spall nor vapor is enabled. Axis/floor default to the beam's.
    static void ParseBitSettings(MMWSpalling& value, IO::ParmParse& pp)
    {
        int on = 0;
        pp.query_default("bit.enabled", on, 0);
        value.bit_enabled = (on != 0);
        if (!value.bit_enabled) return;

        if (!value.spall_enabled && !value.vapor_enabled)
            Util::Abort(INFO, "bit.enabled requires spall.enabled or vapor.enabled "
                              "(the bit reuses the phi/removed surface machinery)");

        pp.query_required("bit.feed_rate", value.bit_feed_rate);
        if (!(std::isfinite(value.bit_feed_rate) && value.bit_feed_rate >= 0.0))
            Util::Abort(INFO, "bit.feed_rate must be finite and >= 0, got ",
                        value.bit_feed_rate);
        pp.query_required("bit.radius", value.bit_radius);
        if (!(std::isfinite(value.bit_radius) && value.bit_radius > 0.0))
            Util::Abort(INFO, "bit.radius must be finite and > 0, got ",
                        value.bit_radius);
        pp.query_default("bit.x0",         value.bit_x0,         value.beam.x0);
        pp.query_default("bit.y0",         value.bit_y0,         value.beam.y0);
        pp.query_default("bit.z0",         value.bit_z0,         value.beam.z0);
        pp.query_default("bit.start_time", value.bit_start_time, 0.0);
    }

    // Hoek-Brown anisotropy ratio sigma_1 / sigma_3 (compression-positive).
    // Returns TENSILE_RATIO_SENTINEL (1e30) on mixed/tensile states so the
    // call-site `ratio >= xi_s` gate passes unconditionally there — tensile
    // sigma_zz at a free surface is more spall-prone, not less.
    static AMREX_FORCE_INLINE
    Set::Scalar PrincipalStressRatio(const Set::Matrix& sigma)
    {
        constexpr Set::Scalar TENSILE_RATIO_SENTINEL = 1.0e30;
        Eigen::SelfAdjointEigenSolver<Set::Matrix> es(sigma);
        if (es.info() != Eigen::Success) return 0.0;
        const auto& e = es.eigenvalues();
        Set::Scalar e_min = e(0);
        Set::Scalar e_max = e(e.size() - 1);
        Set::Scalar sigma_1 = -e_min;
        Set::Scalar sigma_3 = -e_max;
        if (!(sigma_3 > 1.0e-30)) return TENSILE_RATIO_SENTINEL;
        Set::Scalar r = sigma_1 / sigma_3;
        return std::isfinite(r) ? r : TENSILE_RATIO_SENTINEL;
    }

    // Fill k_top_col[ncols] with the topmost live (non-void, phi >= 0) cell of
    // each (i, j) column, gathered across MPI z-decomposition via ReduceIntMax
    // (mirrors UpdateRemovalAfterCohesive PASS 1). Without removal active every
    // entry is the fixed domain top `top_k`; with removal it tracks the receded
    // surface so the surface diagnostics (Sp_field, Sp clusters) follow the
    // drilled face instead of reporting 0 over every spalled column. A column
    // that is fully removed gets the INT_MIN sentinel (no live surface ->
    // nothing written there).
    void GatherColumnSurface(int lev, amrex::Vector<int>& k_top_col,
                             int dom_ilo, int dom_jlo, int nx_dom, int top_k)
    {
        const int ncols = static_cast<int>(k_top_col.size());
        if (!removal_enabled())
        {
            for (int c = 0; c < ncols; ++c) k_top_col[c] = top_k;
            return;
        }
        for (int c = 0; c < ncols; ++c) k_top_col[c] = std::numeric_limits<int>::min();
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                int ktl = std::numeric_limits<int>::min();
                for (int k = khi; k >= klo; --k)
                    if (rem(i, j, k) <= 0.5 && phi(i, j, k) >= 0.0) { ktl = k; break; }
                if (ktl > std::numeric_limits<int>::min())
                {
                    const int idx = (j - dom_jlo) * nx_dom + (i - dom_ilo);
                    if (ktl > k_top_col[idx]) k_top_col[idx] = ktl;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(k_top_col.data(), ncols);
    }

    // Sp = K_I/K_Ic(T) surface diagnostic for `spallation.model = sp_weibull`.
    // Writes Sp into the topmost LIVE cell of each (i, j) column, 0 elsewhere.
    // Under active removal the surface recedes, so the write tracks the drilled
    // face (GatherColumnSurface); with no removal it is the fixed domain top.
    // sigma_xx source: prescribed parser σ(x,y,z,t) if sp_prescribed_sigma is
    // enabled (used by the closed-form verification test); otherwise the
    // depth-resolved physical-path NodeToCellAverage of stress_mf. T_surface
    // for K_Ic: prescribed T(t) parser if enabled, otherwise temp_mf top cell.
    // Step 20: per-column cos(theta) of the local cavity-wall inclination,
    // from a domain-gathered column-top array (so the result is identical on
    // every rank regardless of decomposition). Surface height h = (k_top+1)*dz
    // (the constant offset cancels in the gradient); central differences where
    // both neighbours are valid, one-sided at edges/voided columns, slope 0
    // when no valid neighbour. cos = 1/sqrt(1+|grad h|^2), clamped to
    // >= SP_SN_COS_MIN. Columns with k_top < k_valid_min (fully removed /
    // sentinel) get cos = 1 (unused -- they have no surface to scan).
    static void ColumnSurfaceCos(const int* k_top, int nx, int ny,
                                 Set::Scalar dx, Set::Scalar dy, Set::Scalar dz,
                                 int k_valid_min,
                                 amrex::Vector<Set::Scalar>& cos_out)
    {
        cos_out.assign(static_cast<std::size_t>(nx) * ny, 1.0);
        auto h_at = [&](int ci, int cj) -> Set::Scalar {
            return (Set::Scalar(k_top[cj * nx + ci]) + 1.0) * dz;
        };
        auto valid = [&](int ci, int cj) -> bool {
            if (ci < 0 || ci >= nx || cj < 0 || cj >= ny) return false;
            return k_top[cj * nx + ci] >= k_valid_min;
        };
        for (int cj = 0; cj < ny; ++cj)
        for (int ci = 0; ci < nx; ++ci)
        {
            if (!valid(ci, cj)) continue;
            const Set::Scalar h0 = h_at(ci, cj);
            auto slope1d = [&](int im, int jm, int ip, int jp,
                               Set::Scalar dl) -> Set::Scalar {
                const bool vm = valid(im, jm), vp = valid(ip, jp);
                if (vm && vp) return (h_at(ip, jp) - h_at(im, jm)) / (2.0 * dl);
                if (vp)       return (h_at(ip, jp) - h0) / dl;
                if (vm)       return (h0 - h_at(im, jm)) / dl;
                return 0.0;
            };
            const Set::Scalar gx = slope1d(ci - 1, cj, ci + 1, cj, dx);
            const Set::Scalar gy = slope1d(ci, cj - 1, ci, cj + 1, dy);
            Set::Scalar c = 1.0 / std::sqrt(1.0 + gx * gx + gy * gy);
            if (c < SP_SN_COS_MIN) c = SP_SN_COS_MIN;
            cos_out[static_cast<std::size_t>(cj) * nx + ci] = c;
        }
    }

    void UpdateSpAfterMechanics(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (spallation_model != SpallationModel::SpWeibull) return;
        (void)dt;

        const amrex::Box& domain = geom[lev].Domain();
        const int top_k = domain.bigEnd(2);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar* PHI = geom[lev].ProbHi();
        const Set::Scalar z_top_face = PHI[2];
        const Set::Scalar dz = DX[2];

        const bool        use_presc_s = sp_prescribed_sigma_enabled;
        const bool        use_presc_T = sp_prescribed_T_surface_enabled;
        auto              presc_s_f   = sp_prescribed_sigma_f;
        auto              presc_T_f   = sp_prescribed_T_surface_f;
        const bool        use_weibull = weibull_enabled;

        // local_thermoelastic stress source: reconstruct sigma_xx from T(z)
        // instead of sampling stress_mf. Phase-0 params validated here at
        // first use. Mutually exclusive with the prescribed-sigma verification
        // path (use_presc_s wins, preserving that test byte-for-byte).
        const bool        use_local_TE = sp_local_thermoelastic && !use_presc_s;
        Set::Scalar te_E = 0.0, te_beta = 0.0, te_nu = 0.0, te_Tref = 0.0;
        if (use_local_TE)
        {
            if (phases.empty() || !(phases[0].E > 0.0 && phases[0].mu > 0.0))
                Util::Abort(INFO, "spall.stress_source = local_thermoelastic "
                                  "needs phase 0 with E > 0 and mu > 0");
            te_E    = phases[0].E;
            te_beta = phases[0].beta;
            te_nu   = phases[0].E / (2.0 * phases[0].mu) - 1.0;
            te_Tref = phases[0].T_ref;
            if (!(te_nu > -1.0 && te_nu < 0.5))
                Util::Abort(INFO, "local_thermoelastic: phase 0 Poisson's "
                                  "ratio out of physical range, got ", te_nu);
        }
        const Set::Scalar a0_loc      = sp_a0;
        const int         n_seg       = sp_n_segments;
        const KIcLaw      law         = kic_law;
        const Set::Scalar kic_scale_  = kic_scale;
        const Set::Scalar K_Ic0_      = kic_K_Ic0;
        const Set::Scalar T_ref_kic_  = kic_T_ref;
        const Set::Scalar T_melt_     = kic_T_melt;

        // Confining-pressure contribution +p·nu/(1-nu) to compression-positive
        // sigma_xx (Kant Eq. 15, 1-D-confined geometry). Exact-zero short-
        // circuit preserves byte-identity when confining.p is unset.
        // Single-phase scope (uses phases[0]); multi-phase not supported.
        Set::Scalar sp_conf_offset = 0.0;
        if (confining_p != 0.0 && !phases.empty() &&
            phases[0].E > 0.0 && phases[0].mu > 0.0)
        {
            const Set::Scalar nu_local =
                phases[0].E / (2.0 * phases[0].mu) - 1.0;
            if (nu_local > -1.0 && nu_local < 1.0)
                sp_conf_offset = confining_p * nu_local / (1.0 - nu_local);
        }

        if (!use_presc_s && !use_local_TE)
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

        // Per-column live surface: with removal active the surface recedes from
        // the fixed domain top, so the Sp diagnostic must be written at the
        // topmost surviving cell of each column (else every spalled column
        // reads 0 and looks inactive). Gathered z-decomposition-safe.
        const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
        const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
        const int nx_dom = dom_ihi - dom_ilo + 1;
        const int ny_dom = dom_jhi - dom_jlo + 1;
        const int ncols  = nx_dom * ny_dom;
        amrex::Vector<int> k_top_col(ncols, top_k);
        GatherColumnSurface(lev, k_top_col, dom_ilo, dom_jlo, nx_dom, top_k);
        const int dom_ilo_l = dom_ilo, dom_jlo_l = dom_jlo, nx_dom_l = nx_dom;
        const int* k_top_col_p = k_top_col.data();
        auto col_index_sp = [=] (int i, int j) {
            return (j - dom_jlo_l) * nx_dom_l + (i - dom_ilo_l);
        };

        // Step 20: surface-normal frame -- per-column cos(theta) so the K_I
        // scan below samples the vertical T/stress profile at depth/cos(theta)
        // (the crack runs along the local surface normal). nullptr when off:
        // the conversion multiplies by exactly 1.0 and stays bit-identical.
        amrex::Vector<Set::Scalar> sn_cos_col;
        if (spall_surface_normal)
            ColumnSurfaceCos(k_top_col.data(), nx_dom, ny_dom,
                             DX[0], DX[1], dz, domain.smallEnd(2), sn_cos_col);
        const Set::Scalar* sn_cos_p =
            sn_cos_col.empty() ? nullptr : sn_cos_col.data();

        // Sp is a per-column diagnostic; clear everywhere first so cells
        // below the top z-slab stay zero on plotfiles.
        Sp_field_mf[lev]->setVal(0.0);

        for (amrex::MFIter mfi(*Sp_field_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const int vbx_lo_k = mfi.validbox().smallEnd(2);
            amrex::Array4<Set::Scalar>       const& Sp_a  = (*Sp_field_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T_a   = (*temp_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> sig_node;
            if (!use_presc_s && !use_local_TE) sig_node = (*stress_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> flaw_a;
            if (use_weibull) flaw_a = (*flaw_a_mf[lev]).array(mfi);

            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const int ksurf = k_top_col_p[col_index_sp(i, j)];
                if (k != ksurf) return;
                // Receded surface top face (== domain top z when nothing above
                // this column has been removed, since then ksurf == top_k).
                const Set::Scalar z_surf_face =
                    z_top_face - Set::Scalar(top_k - ksurf) * dz;

                Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];

                Set::Scalar Tsurf = use_presc_T ? presc_T_f(time) : T_a(i, j, k);

                // Step 20: under surface_normal the crack frame is the local
                // surface normal; `depth` below is a NORMAL distance and the
                // vertical profile is sampled at depth/cos(theta). icth = 1
                // exactly when the flag is off (bit-identical).
                const Set::Scalar icth = sn_cos_p
                    ? 1.0 / sn_cos_p[col_index_sp(i, j)] : 1.0;

                // K_I integrand consumes -sigma_xx (compression-positive).
                // stress_mf and the prescribed parser are tension-positive,
                // so a compressive heated surface (sigma_xx < 0) yields
                // K_I > 0. The physical path samples stress_mf at the cell
                // containing the requested depth (k_d), clamped into the
                // owned box. z-decomposition-robust column gather is deferred.
                auto sigma_at_depth = [&] (Set::Scalar depth) -> Set::Scalar {
                    const Set::Scalar dv = depth * icth;  // vertical sample depth
                    if (use_presc_s)
                    {
                        Set::Scalar z_phys = z_surf_face - dv;
                        return -presc_s_f(xc, yc, z_phys, time) + sp_conf_offset;
                    }
                    if (use_local_TE)
                    {
                        // 1-D-confined thermoelastic sigma_xx from T(z), same
                        // form and compression-positive convention as the
                        // h_spall scan (Removal.H). Sample temp_mf at the cell
                        // containing the requested depth below the (receded)
                        // surface, clamped into the box.
                        int k_d = ksurf - static_cast<int>(std::floor(dv / dz));
                        if (k_d > ksurf)    k_d = ksurf;
                        if (k_d < vbx_lo_k) k_d = vbx_lo_k;
                        return te_E * te_beta * (T_a(i, j, k_d) - te_Tref)
                               / (1.0 - te_nu) + sp_conf_offset;
                    }
                    int k_d = ksurf - static_cast<int>(std::floor(dv / dz));
                    if (k_d > ksurf)    k_d = ksurf;
                    if (k_d < vbx_lo_k) k_d = vbx_lo_k;
                    Set::Matrix sigma =
                        Numeric::Interpolate::NodeToCellAverage(sig_node, i, j, k_d, 0);
                    return -sigma(0, 0) + sp_conf_offset;
                };

                Set::Scalar a_loc = use_weibull ? flaw_a(i, j, k) : a0_loc;

                Set::Scalar K_I_val =
                    Numeric::SpCriterion::K_I(sigma_at_depth, a_loc, n_seg);

                Set::Scalar K_Ic_MPa;
                if (law == KIcLaw::NasseriTable)
                    K_Ic_MPa = Numeric::SpCriterion::KIcNasseri(Tsurf, kic_scale_);
                else
                    K_Ic_MPa = Numeric::SpCriterion::KIcLinear(
                        Tsurf, K_Ic0_, T_ref_kic_, T_melt_);

                Set::Scalar K_Ic_Pa = Numeric::SpCriterion::MPa_sqrtm_to_Pa_sqrtm(K_Ic_MPa);
                Sp_a(i, j, k) = (K_Ic_Pa > 0.0) ? K_I_val / K_Ic_Pa : 0.0;
            });
        }
    }

    // Step 15c: connected-component labeller on the Sp >= 1 surface zone.
    // Single-level, 4-connected on (i, j) at the top-z plane. Gathers the
    // domain-wide firing mask to the I/O processor via ReduceIntMax, runs a
    // sequential union-find, and broadcasts the merged label image back. For
    // a 24^3 run the top plane is 576 cells — bandwidth and cost negligible.
    // Sp_field_mf is never modified here. Cross-level connectivity is
    // deferred (the active sp_kant_onset / sp_weibull runs are max_level = 0
    // because elastic MLMG on anisotropic cells does not converge).
    void UpdateSpClusters(int lev, Set::Scalar time, Set::Scalar dt)
    {
        (void)dt;
        if (spallation_model != SpallationModel::SpWeibull) return;
        if (!weibull_enabled) return;

        const amrex::Box& domain = geom[lev].Domain();
        const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
        const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
        const int nx_dom = dom_ihi - dom_ilo + 1;
        const int ny_dom = dom_jhi - dom_jlo + 1;
        const int ncols  = nx_dom * ny_dom;
        const int top_k  = domain.bigEnd(2);
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar  dx = DX[0];
        const Set::Scalar  dy = DX[1];

        auto col_idx = [&] (int i, int j) {
            return (j - dom_jlo) * nx_dom + (i - dom_ilo);
        };

        // Per-column live surface (receded under active removal); used to write
        // the cluster-id field at the drilled face so it aligns with Sp_field.
        amrex::Vector<int> k_top_col(ncols, top_k);
        GatherColumnSurface(lev, k_top_col, dom_ilo, dom_jlo, nx_dom, top_k);

        // Domain-wide firing buffer: 1 on firing surface-cells, 0 elsewhere.
        // ReduceIntMax across ranks (the rank that owns a firing cell writes 1,
        // others stay 0 — max collapses cleanly). Sp_field is nonzero only at
        // each column's live surface cell, so scanning the full z-extent finds
        // exactly the firing surfaces (pre-removal that is the domain-top row).
        amrex::Vector<int> fire(ncols, 0);
        for (amrex::MFIter mfi(*Sp_field_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
            amrex::Array4<const Set::Scalar> const& Sp_a =
                (*Sp_field_mf[lev]).array(mfi);
            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            for (int k = klo; k <= khi; ++k)
            {
                if (Sp_a(i, j, k) >= 1.0)
                    fire[col_idx(i, j)] = 1;
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(fire.data(), ncols);

        // I/O processor: 4-connected union-find over the (i_dom, j_dom)
        // plane, optional A_crit gate, then build the final label image.
        // Labels in `labels` are 0 (not firing) or positive cluster IDs
        // assigned densely from 1 in scan order of the first cell of each
        // cluster (so IDs are reproducible across runs).
        amrex::Vector<int> labels(ncols, 0);
        int cluster_count = 0;
        Set::Scalar max_cluster_area = 0.0;
        Set::Scalar total_firing_area = 0.0;
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            // Max possible provisional labels: every firing cell its own.
            amrex::Vector<int> parent(static_cast<std::size_t>(ncols + 1), 0);
            int next_label = 1;
            auto find = [&] (int x) {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            };
            auto unite = [&] (int a, int b) {
                int ra = find(a), rb = find(b);
                if (ra == rb) return;
                if (ra < rb) parent[rb] = ra; else parent[ra] = rb;
            };
            // First pass: assign provisional labels using the two-pass
            // connected-component classical algorithm with union-find on
            // already-labelled left/down neighbours.
            for (int jj = 0; jj < ny_dom; ++jj)
            for (int ii = 0; ii < nx_dom; ++ii)
            {
                int idx = jj * nx_dom + ii;
                if (!fire[idx]) continue;
                int left  = (ii > 0) ? labels[idx - 1]      : 0;
                int down  = (jj > 0) ? labels[idx - nx_dom] : 0;
                if (left == 0 && down == 0)
                {
                    int nl = next_label++;
                    parent[nl] = nl;
                    labels[idx] = nl;
                }
                else if (left != 0 && down == 0) labels[idx] = left;
                else if (left == 0 && down != 0) labels[idx] = down;
                else { labels[idx] = left; unite(left, down); }
            }
            // Second pass: resolve aliases and renumber densely.
            amrex::Vector<int> remap(static_cast<std::size_t>(next_label), 0);
            int dense = 0;
            for (int idx = 0; idx < ncols; ++idx)
            {
                if (!labels[idx]) continue;
                int r = find(labels[idx]);
                if (remap[r] == 0) remap[r] = ++dense;
                labels[idx] = remap[r];
            }
            cluster_count = dense;

            // Per-cluster area (number of firing cells × dx × dy).
            amrex::Vector<int> cell_count(static_cast<std::size_t>(dense + 1), 0);
            for (int idx = 0; idx < ncols; ++idx)
                if (labels[idx]) cell_count[labels[idx]]++;
            const Set::Scalar A_cell = dx * dy;

            if (weibull_detach_mode == DetachMode::ConnectedCluster &&
                weibull_A_crit > 0.0)
            {
                amrex::Vector<int> keep(static_cast<std::size_t>(dense + 1), 0);
                int new_dense = 0;
                for (int c = 1; c <= dense; ++c)
                {
                    Set::Scalar area = cell_count[c] * A_cell;
                    if (area >= weibull_A_crit) keep[c] = ++new_dense;
                }
                amrex::Vector<int> new_count(static_cast<std::size_t>(new_dense + 1), 0);
                for (int idx = 0; idx < ncols; ++idx)
                {
                    int old = labels[idx];
                    int kp  = old ? keep[old] : 0;
                    labels[idx] = kp;
                    if (kp) new_count[kp]++;
                }
                cluster_count = new_dense;
                cell_count.swap(new_count);
                dense = new_dense;
            }

            int total_cells = 0;
            int max_count   = 0;
            for (int c = 1; c <= dense; ++c)
            {
                total_cells += cell_count[c];
                if (cell_count[c] > max_count) max_count = cell_count[c];
            }
            max_cluster_area  = max_count   * A_cell;
            total_firing_area = total_cells * A_cell;
        }

        // Broadcast the label image and cluster diagnostics so every rank
        // can write Sp_cluster_id_mf locally.
        amrex::ParallelDescriptor::Bcast(labels.data(), ncols,
                                         amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&cluster_count, 1,
                                         amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&max_cluster_area, 1,
                                         amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&total_firing_area, 1,
                                         amrex::ParallelDescriptor::IOProcessorNumber());

        // Write Sp_cluster_id_mf everywhere on this level: each column's live
        // (possibly receded) surface cell gets the cluster ID (or 0 for
        // filtered-out / non-firing cells); 0 elsewhere. Tracks the drilled
        // face so the field aligns with Sp_field.
        const int* k_top_col_p = k_top_col.data();
        Sp_cluster_id_mf[lev]->setVal(0.0);
        for (amrex::MFIter mfi(*Sp_cluster_id_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& cl =
                (*Sp_cluster_id_mf[lev]).array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (k != k_top_col_p[col_idx(i, j)]) return;
                cl(i, j, k) = static_cast<Set::Scalar>(labels[col_idx(i, j)]);
            });
        }

        // Per-step CSV diagnostics (I/O processor only). Only emit a row
        // when at least one cluster fires; mirrors the active step's "log
        // when cluster_count > 0" cheapest-path spec. The CSV lives next to
        // the plotfile directory (sibling file `<plot_file>_clusters.csv`)
        // so it persists across runs of one binary invocation but is
        // truncated and re-written on the first labeller call of a fresh
        // run (sp_clusters_csv_header_written is per-MMWSpalling-instance).
        if (lev == 0 && cluster_count > 0 && sp_emit_clusters_csv &&
            amrex::ParallelDescriptor::IOProcessor())
        {
            std::string csv_path = plot_file + "_clusters.csv";
            std::ios::openmode mode = sp_clusters_csv_header_written
                                          ? std::ios::app
                                          : std::ios::trunc;
            std::ofstream ofs(csv_path, mode | std::ios::out);
            if (ofs.is_open())
            {
                if (!sp_clusters_csv_header_written)
                {
                    ofs << "time,cluster_count,max_cluster_area_m2,"
                        << "total_firing_area_m2\n";
                    sp_clusters_csv_header_written = true;
                }
                ofs.setf(std::ios::scientific);
                ofs.precision(8);
                ofs << time << "," << cluster_count << ","
                    << max_cluster_area << "," << total_firing_area << "\n";
            }
        }
    }

    // UpdateRemovalAfterCohesive — extracted to a partial header to keep this
    // file navigable. The definition is still an inline member of MMWSpalling,
    // textually included here at its original position (R1 refactor).
    #include "Integrator/MMWSpalling/Removal.H"

    // Build the global-material H<->T table once. Microstructure tables
    // are built per phase in InitializeMicrostructure(); see phase_tables.
    void BuildEnthalpyTable()
    {
        if (!material) return;
        material_table.Build(*material);
    }

    // Const-α (no material) phase-fraction tagging used by ConvertTtoH and
    // InvertHtoT for the closed-form path. The material and microstructure
    // paths use Numeric::Material::Table::Invert directly.
    void SetPhaseFractionsConstAlpha(Set::Scalar T,
                                     Set::Scalar &LS, Set::Scalar &LL, Set::Scalar &LV) const
    {
        const Set::Scalar Tm = T_m, Tv = T_v;
        if (T < Tm)       { LS = 1.0; LL = 0.0; LV = 0.0; }
        else if (T <= Tv) { LS = 0.0; LL = 1.0; LV = 0.0; }
        else              { LS = 0.0; LL = 0.0; LV = 1.0; }
    }

    void ConvertTtoH(int lev)
    {
        // Step 6b dispatch:
        //   microstructure → per-cell phase table (this->phase_tables[pid])
        //   global material → material_table
        //   else            → const-α closed form
        const Numeric::Material::Table *gtab = material ? &material_table : nullptr;
        const Numeric::Material::Table *ptabs =
            phase_tables.empty() ? nullptr : phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        const bool micro = microstructure_enabled;
        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& H    = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Hold = (*H_old_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph =
                micro ? (*phase_mf[lev]).array(mfi) : amrex::Array4<const Set::Scalar>{};
            amrex::ParallelFor(bx, [=, this] (int i, int j, int k) {
                Set::Scalar Tv = T(i, j, k);
                Set::Scalar Hv = 0.0;
                Set::Scalar ls = 1.0, ll = 0.0, lv = 0.0;
                if (micro) {
                    int p = static_cast<int>(ph(i, j, k));
                    if (p < 0 || p >= n_phases) p = 0;
                    Hv = ptabs[p].H_of_T(Tv);
                    Set::Scalar dummyT;
                    ptabs[p].Invert(Hv, dummyT, ls, ll, lv);
                } else if (gtab) {
                    Hv = gtab->H_of_T(Tv);
                    Set::Scalar dummyT;
                    gtab->Invert(Hv, dummyT, ls, ll, lv);
                } else {
                    const Set::Scalar Tref = T_ref, Tm = T_m, Tvm = T_v;
                    const Set::Scalar HS  = Tm - Tref;
                    const Set::Scalar HL  = HS + lh_m;
                    const Set::Scalar HLV = HL + (Tvm - Tm);
                    const Set::Scalar HV  = HLV + lh_v;
                    if (Tv <= Tm)       Hv = Tv - Tref;
                    else if (Tv <= Tvm) Hv = HL + (Tv - Tm);
                    else                Hv = HV + (Tv - Tvm);
                    SetPhaseFractionsConstAlpha(Tv, ls, ll, lv);
                }
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
                H(i,j,k)    = Hv;
                Hold(i,j,k) = Hv;
            });
        }
    }

    void InvertHtoT(int lev)
    {
        const Numeric::Material::Table *gtab = material ? &material_table : nullptr;
        const Numeric::Material::Table *ptabs =
            phase_tables.empty() ? nullptr : phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        const bool micro = microstructure_enabled;
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& H = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& T  = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph =
                micro ? (*phase_mf[lev]).array(mfi) : amrex::Array4<const Set::Scalar>{};
            amrex::ParallelFor(bx, [=, this] (int i, int j, int k) {
                Set::Scalar Hv = H(i, j, k);
                if (micro) {
                    int p = static_cast<int>(ph(i, j, k));
                    if (p < 0 || p >= n_phases) p = 0;
                    Set::Scalar Tv, ls, ll, lv;
                    ptabs[p].Invert(Hv, Tv, ls, ll, lv);
                    T(i,j,k) = Tv;
                    LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
                } else if (gtab) {
                    Set::Scalar Tv, ls, ll, lv;
                    gtab->Invert(Hv, Tv, ls, ll, lv);
                    T(i,j,k) = Tv;
                    LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
                } else {
                    // Const-α path: phase fraction inferred from H during a
                    // sharp plateau (T sits at T_m or T_v across the latent
                    // heat absorption).
                    const Set::Scalar Tref = T_ref, Tm = T_m, Tvm = T_v;
                    const Set::Scalar Lhm = lh_m, Lhv = lh_v;
                    const Set::Scalar HS  = Tm - Tref;
                    const Set::Scalar HL  = HS + Lhm;
                    const Set::Scalar HLV = HL + (Tvm - Tm);
                    const Set::Scalar HV  = HLV + Lhv;
                    if (Hv < HS) {
                        T(i,j,k)  = Tref + Hv;
                        LS(i,j,k) = 1.0; LL(i,j,k) = 0.0; LV(i,j,k) = 0.0;
                    } else if (Hv <= HL) {
                        T(i,j,k)  = Tm;
                        Set::Scalar fL = (Lhm > 0.0) ? (Hv - HS) / Lhm : 1.0;
                        LL(i,j,k) = fL; LS(i,j,k) = 1.0 - fL; LV(i,j,k) = 0.0;
                    } else if (Hv < HLV) {
                        T(i,j,k)  = Tm + (Hv - HL);
                        LS(i,j,k) = 0.0; LL(i,j,k) = 1.0; LV(i,j,k) = 0.0;
                    } else if (Hv <= HV) {
                        T(i,j,k)  = Tvm;
                        Set::Scalar fV = (Lhv > 0.0) ? (Hv - HLV) / Lhv : 1.0;
                        LV(i,j,k) = fV; LL(i,j,k) = 1.0 - fV; LS(i,j,k) = 0.0;
                    } else {
                        T(i,j,k)  = Tvm + (Hv - HV);
                        LS(i,j,k) = 0.0; LL(i,j,k) = 0.0; LV(i,j,k) = 1.0;
                    }
                }
            });
        }
    }

    void InitSurfaceMask(int lev)
    {
        if (removal_enabled())
        {
            UpdateSurfaceMaskFromPhi(lev);
            return;
        }
        const int top_k = geom[lev].Domain().bigEnd(2);
        for (amrex::MFIter mfi(*surface_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& mask = (*surface_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                mask(i, j, k) = (k == top_k) ? 1.0 : 0.0;
            });
        }
    }

    // Surface cells are material cells (phi >= 0) that are the highest
    // material cell in their column, equivalently 0 <= phi < dz. Implement
    // the simple per-column scan; AMR is single-level for Step 10.
    void UpdateSurfaceMaskFromPhi(int lev)
    {
        const Set::Scalar dz = geom[lev].CellSize()[2];
        for (amrex::MFIter mfi(*surface_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& mask = (*surface_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar p = phi(i, j, k);
                mask(i, j, k) = (p >= 0.0 && p < dz) ? 1.0 : 0.0;
            });
        }
    }

    // Initialize phi to a flat top surface at z = prob_hi[2]. The cell-
    // centred value is z_top - z_center; with cell index k counting from
    // the bottom, that is (Nz - 0.5 - k) * dz.
    void InitPhi(int lev)
    {
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PHI = geom[lev].ProbHi();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar z_top = PHI[2];
        const Set::Scalar z_lo  = PLO[2];
        (void)z_lo;
        for (amrex::MFIter mfi(*phi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar zc = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                phi(i, j, k) = z_top - zc;
                rem(i, j, k) = 0.0;
            });
        }
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
    }

    void RefreshMicrostructureDerivedFields(bool announce)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            InitSurfaceMask(lev);
            InitializeMicrostructure(lev, announce, !damage_enabled);
        }
    }

    void RepairMechanicsAfterRegrid(Set::Scalar time)
    {
        if (m_type == Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            return;

        if (microstructure_enabled && max_level > 0)
            RefreshMicrostructureDerivedFields(/*announce=*/false);

        UpdateModel(/*a_step=*/0, time);

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (m_zero_out_displacement)
                disp_mf[lev]->setVal(Set::Vector::Zero());
            stress_mf[lev]->setVal(Set::Matrix::Zero());
            strain_mf[lev]->setVal(Set::Matrix::Zero());

            Util::RealFillBoundary(*disp_mf[lev], geom[lev]);
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);
            Util::RealFillBoundary(*strain_mf[lev], geom[lev]);
        }
    }

    static bool FiniteMatrix(const Set::Matrix& m)
    {
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            for (int j = 0; j < AMREX_SPACEDIM; ++j)
                if (!std::isfinite(m(i, j))) return false;
        return true;
    }


    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar AlphaBetaFraction(Set::Scalar T, Set::Scalar T0,
                                         Set::Scalar width)
    {
        if (width <= 0.0) return (T >= T0) ? 1.0 : 0.0;
        Set::Scalar lo = T0 - 0.5 * width;
        Set::Scalar s = (T - lo) / width;
        if (s <= 0.0) return 0.0;
        if (s >= 1.0) return 1.0;
        return s * s * (3.0 - 2.0 * s);
    }

    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ThermalFreeStrain(Set::Scalar T,
                                         Set::Scalar beta,
                                         Set::Scalar T_ref,
                                         Set::Scalar alpha_beta_eps,
                                         Set::Scalar alpha_beta_T0,
                                         Set::Scalar alpha_beta_width)
    {
        return beta * (T - T_ref)
             + alpha_beta_eps * AlphaBetaFraction(T, alpha_beta_T0,
                                                  alpha_beta_width);
    }

    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar CellToNodeAverageInDomain(const amrex::Array4<const Set::Scalar>& f,
                                          int i, int j, int k, int comp,
                                          int ilo, int ihi,
                                          int jlo, int jhi,
                                          int klo, int khi,
                                          Set::Scalar max_abs)
    {
        Set::Scalar sum = 0.0;
        int count = 0;
        for (int dk = 0; dk <= 1; ++dk)
        {
            int kk = k - dk;
            if (kk < klo || kk > khi) continue;
            for (int dj = 0; dj <= 1; ++dj)
            {
                int jj = j - dj;
                if (jj < jlo || jj > jhi) continue;
                for (int di = 0; di <= 1; ++di)
                {
                    int ii = i - di;
                    if (ii < ilo || ii > ihi) continue;
                    Set::Scalar v = f(ii, jj, kk, comp);
                    if (!std::isfinite(v) || std::abs(v) > max_abs) continue;
                    sum += v;
                    ++count;
                }
            }
        }
        return (count > 0) ? sum / Set::Scalar(count) : 0.0;
    }

    void UpdateMicrostructureFreeStrain(int lev)
    {
        for (MFIter mfi(*thermal_free_strain_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& temp = temp_old_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& beta = beta_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& Tref = T_ref_phase_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& abe  = alpha_beta_linear_strain_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& abT0 = alpha_beta_T0_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& abW  = alpha_beta_width_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& eps = thermal_free_strain_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                eps(i, j, k) = ThermalFreeStrain(
                    temp(i, j, k, 0),
                    beta(i, j, k, 0),
                    Tref(i, j, k, 0),
                    abe (i, j, k, 0),
                    abT0(i, j, k, 0),
                    abW (i, j, k, 0));
            });
        }
        Util::RealFillBoundary(*thermal_free_strain_mf[lev], geom[lev]);
    }

    // Step 5: build the Voronoi grain field, mineral assignment, heterogeneous
    // property fields, and grain-boundary flag for level `lev`. Grain centres
    // and the grain→phase mapping are generated once (on the first call) using
    // a deterministic mt19937 seeded from `microstructure.seed`, so the result
    // is reproducible and consistent across AMR levels.
    void InitializeMicrostructure(int lev, bool announce = true,
                                  bool reset_damage = true)
    {
        // Step 6b: build per-phase Constant material and H<->T table once
        // (lazy on first level call). Each phase's Constant captures the
        // user-specified rho, Cp, kappa plus the optional phase-change
        // scalars; the Table integrates rho*cp_eff to give H(T) and the
        // inverse, with phase fractions inferred from H.
        if (phase_tables.size() != phases.size())
        {
            phase_materials.clear();
            phase_tables.clear();
            phase_materials.resize(phases.size());
            phase_tables.resize(phases.size());
            for (size_t p = 0; p < phases.size(); ++p)
            {
                phase_materials[p].m_rho   = phases[p].rho;
                phase_materials[p].m_cp    = phases[p].Cp;
                phase_materials[p].m_kappa = phases[p].kappa;
                phase_materials[p].T_ref     = phases[p].T_ref;
                phase_materials[p].T_m       = phases[p].T_m;
                phase_materials[p].fus_width = phases[p].fus_width;
                phase_materials[p].L_m       = phases[p].L_m;
                phase_materials[p].T_vap_lo  = phases[p].T_vap_lo;
                phase_materials[p].T_vap_hi  = phases[p].T_vap_hi;
                phase_materials[p].L_v       = phases[p].L_v;
                phase_tables[p].Build(phase_materials[p]);
            }
        }

        // Lazy one-shot generation of grain centres + per-grain mineral
        // assignment. We do the tessellation locally rather than going through
        // IC::Voronoi because IC::Voronoi::Parse reads its `seed` member only
        // after Define() runs, so its grain placement isn't actually under
        // user-seed control. Doing it inline keeps the microstructure fully
        // reproducible from `microstructure.seed`. Skipped entirely in
        // expression mode — Step 6 prescribes the phase IDs from a
        // user-supplied amrex::Parser.
        if (micro_mode == MicroMode::Voronoi && grain_centres.empty())
        {
            std::mt19937 rng(static_cast<unsigned int>(micro_seed));
            std::uniform_real_distribution<Set::Scalar> uni(0.0, 1.0);
            const Set::Scalar* lo = geom[0].ProbLo();
            const Set::Scalar* hi = geom[0].ProbHi();
            grain_centres.assign(3 * num_grains, 0.0);
            for (int g = 0; g < num_grains; ++g)
            {
                grain_centres[3*g + 0] = lo[0] + (hi[0] - lo[0]) * uni(rng);
                grain_centres[3*g + 1] = lo[1] + (hi[1] - lo[1]) * uni(rng);
                grain_centres[3*g + 2] = lo[2] + (hi[2] - lo[2]) * uni(rng);
            }
            std::vector<Set::Scalar> cdf(phases.size(), 0.0);
            Set::Scalar acc = 0.0;
            for (size_t p = 0; p < phases.size(); ++p)
            {
                acc += phases[p].fraction;
                cdf[p] = acc;
            }
            grain_phase_ids.assign(num_grains, 0);
            for (int g = 0; g < num_grains; ++g)
            {
                Set::Scalar u = uni(rng);
                int pid = (int)phases.size() - 1;
                for (size_t p = 0; p < phases.size(); ++p)
                    if (u <= cdf[p]) { pid = (int)p; break; }
                grain_phase_ids[g] = pid;
            }
        }

        // Assign each cell to a phase ID. Ghost cells get a -1 sentinel so
        // out-of-domain ghosts are identifiable in the is_gb pass. The kernel
        // branches on `micro_mode`: Voronoi looks up the nearest grain centre;
        // Expression evaluates the user parser at the cell centre and rounds
        // to the nearest integer phase index.
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const int n_phases = static_cast<int>(phases.size());
        phase_mf[lev]->setVal(-1.0, 0, 1, phase_mf[lev]->nGrow());
        grain_id_mf[lev]->setVal(-1.0, 0, 1, grain_id_mf[lev]->nGrow());
        if (micro_mode == MicroMode::Voronoi)
        {
            const int n_grains = num_grains;
            const Set::Scalar* gx = grain_centres.data(); // packed [num_grains*3]
            const int* gpid = grain_phase_ids.data();
            for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& gid = (*grain_id_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                    Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                    Set::Scalar best_d2 = std::numeric_limits<Set::Scalar>::infinity();
                    int best_g = 0;
                    for (int g = 0; g < n_grains; ++g)
                    {
                        Set::Scalar dx_ = x - gx[3*g + 0];
                        Set::Scalar dy_ = y - gx[3*g + 1];
                        Set::Scalar dz_ = z - gx[3*g + 2];
                        Set::Scalar d2 = dx_*dx_ + dy_*dy_ + dz_*dz_;
                        if (d2 < best_d2) { best_d2 = d2; best_g = g; }
                    }
                    // Step 6c: grain_id is the unique grain index; phase is
                    // the per-grain mineral assignment, so two adjacent
                    // same-mineral grains have grain_id != grain_id but
                    // identical phase.
                    gid(i, j, k) = static_cast<Set::Scalar>(best_g);
                    ph (i, j, k) = static_cast<Set::Scalar>(gpid[best_g]);
                });
            }
        }
        else // Expression
        {
            auto pf = phase_f;
            const bool has_gexpr = has_grain_expr;
            auto gf = grain_f;
            for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& gid = (*grain_id_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                    Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                    Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                    Set::Scalar z = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                    int p = static_cast<int>(std::round(pf(x, y, z)));
                    if (p < 0)         p = 0;
                    if (p >= n_phases) p = n_phases - 1;
                    ph(i, j, k) = static_cast<Set::Scalar>(p);
                    // Step 6c: grain_id either follows a user expression or
                    // defaults to the phase id. Expression values are
                    // rounded to a non-negative integer label and never used
                    // as a vector index.
                    int g = p;
                    if (has_gexpr)
                    {
                        int gv = static_cast<int>(std::round(gf(x, y, z)));
                        if (gv < 0) gv = 0;
                        g = gv;
                    }
                    gid(i, j, k) = static_cast<Set::Scalar>(g);
                });
            }
        }

        // Property lookup tables snapshotted onto raw pointers so the
        // ParallelFor capture works with std::vector storage outside the
        // lambda. (Lambdas here intentionally run host-side; properties are
        // a per-cell scalar lookup with no fast-path benefit from GPU.)
        std::vector<Set::Scalar> kap(n_phases), bet(n_phases), Ev(n_phases),
                                 muv(n_phases), Trv(n_phases),
                                 ab_eps(n_phases), ab_T0(n_phases),
                                 ab_width(n_phases), rh(n_phases),
                                 Cpv(n_phases);
        for (int p = 0; p < n_phases; ++p)
        {
            kap[p] = phases[p].kappa;
            bet[p] = phases[p].beta;
            Ev[p]  = phases[p].E;
            muv[p] = phases[p].mu;
            Trv[p] = phases[p].T_ref;
            ab_eps[p]   = phases[p].alpha_beta_linear_strain;
            ab_T0[p]    = phases[p].alpha_beta_T0;
            ab_width[p] = phases[p].alpha_beta_width;
            rh[p]  = phases[p].rho;
            Cpv[p] = phases[p].Cp;
        }
        const Set::Scalar* d_kap = kap.data();
        const Set::Scalar* d_bet = bet.data();
        const Set::Scalar* d_Ev  = Ev.data();
        const Set::Scalar* d_mu  = muv.data();
        const Set::Scalar* d_Tr  = Trv.data();
        const Set::Scalar* d_ab_eps = ab_eps.data();
        const Set::Scalar* d_ab_T0  = ab_T0.data();
        const Set::Scalar* d_ab_w   = ab_width.data();
        const Set::Scalar* d_rh  = rh.data();
        const Set::Scalar* d_Cpv = Cpv.data();

        for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& kapa = (*kappa_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& beta = (*beta_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& EE   = (*E_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& muu  = (*mu_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Trf  = (*T_ref_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& abe  = (*alpha_beta_linear_strain_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& abT0 = (*alpha_beta_T0_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& abW  = (*alpha_beta_width_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rho_a= (*rho_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Cp_a = (*Cp_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                kapa (i, j, k) = d_kap[p];
                beta (i, j, k) = d_bet[p];
                EE   (i, j, k) = d_Ev [p];
                muu  (i, j, k) = d_mu [p];
                Trf  (i, j, k) = d_Tr [p];
                abe  (i, j, k) = d_ab_eps[p];
                abT0 (i, j, k) = d_ab_T0 [p];
                abW  (i, j, k) = d_ab_w  [p];
                rho_a(i, j, k) = d_rh [p];
                Cp_a (i, j, k) = d_Cpv[p];
            });
        }
        // Step 7: fill ghost cells of mechanical-property fields so cell-to-node
        // averaging at box boundaries reads coherent neighbour values.
        Util::RealFillBoundary(*beta_mf[lev],        geom[lev]);
        Util::RealFillBoundary(*E_phase_mf[lev],     geom[lev]);
        Util::RealFillBoundary(*mu_phase_mf[lev],    geom[lev]);
        Util::RealFillBoundary(*T_ref_phase_mf[lev], geom[lev]);
        Util::RealFillBoundary(*alpha_beta_linear_strain_mf[lev], geom[lev]);
        Util::RealFillBoundary(*alpha_beta_T0_mf[lev],            geom[lev]);
        Util::RealFillBoundary(*alpha_beta_width_mf[lev],         geom[lev]);
        UpdateMicrostructureFreeStrain(lev);

        // Fill phase_mf and grain_id_mf ghost cells from neighbouring boxes;
        // out-of-domain ghosts retain the -1 sentinel set above.
        Util::RealFillBoundary(*phase_mf[lev],    geom[lev]);
        Util::RealFillBoundary(*grain_id_mf[lev], geom[lev]);

        const amrex::Box& domain = geom[lev].Domain();
        const int dlo[3] = {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
        const int dhi[3] = {domain.bigEnd(0),   domain.bigEnd(1),   domain.bigEnd(2)};
        // Step 6c: build phase- and grain-boundary flags from the same face
        // sweep. is_gb is preserved as a phase-boundary alias so existing
        // tests (Step 5/6/6b) keep their behaviour; later CZM/damage steps
        // should consume is_grain_boundary instead.
        for (amrex::MFIter mfi(*is_gb_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gid = (*grain_id_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& gb  = (*is_gb_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& pbf = (*is_phase_boundary_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& gbf = (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                int p_self = static_cast<int>(ph(i, j, k));
                int g_self = static_cast<int>(gid(i, j, k));
                int p_flag = 0, g_flag = 0;
                const int di[6] = {-1, 1, 0, 0, 0, 0};
                const int dj[6] = { 0, 0,-1, 1, 0, 0};
                const int dk[6] = { 0, 0, 0, 0,-1, 1};
                for (int n = 0; n < 6; ++n)
                {
                    int ii = i + di[n], jj = j + dj[n], kk = k + dk[n];
                    // Skip neighbours that fall outside the physical domain.
                    if (ii < dlo[0] || ii > dhi[0] ||
                        jj < dlo[1] || jj > dhi[1] ||
                        kk < dlo[2] || kk > dhi[2]) continue;
                    int p_n = static_cast<int>(ph (ii, jj, kk));
                    int g_n = static_cast<int>(gid(ii, jj, kk));
                    if (p_n >= 0 && p_n != p_self) p_flag = 1;
                    if (g_n >= 0 && g_n != g_self) g_flag = 1;
                    if (p_flag && g_flag) break;
                }
                pbf(i, j, k) = static_cast<Set::Scalar>(p_flag);
                gbf(i, j, k) = static_cast<Set::Scalar>(g_flag);
                gb (i, j, k) = static_cast<Set::Scalar>(p_flag);
            });
        }

        // Initialise or preserve the damage field. With damage evolution
        // disabled, Step 6/7 AMR repairs intentionally reapply damage.ic.* so
        // passive-D regressions stay exact. With damage evolution enabled,
        // regrid/repair calls preserve the transferred irreversible D field
        // and refresh only derived quantities below.
        if (reset_damage)
        {
            D_mf[lev]->setVal(0.0);
            if (ic_damage) ic_damage->Initialize(lev, D_mf);
        }
        Util::RealFillBoundary(*D_mf[lev], geom[lev]);

        // Damage-modified per-cell conductivity diagnostic (κ_eff = κ_phase ·
        // exp(−α_D · D)). Refreshed each step in AdvanceMicrostructure too,
        // but populated here so the t=0 plotfile is correct.
        RefreshKappaEff(lev);
        if (damage_enabled && reset_damage)
        {
            damage_F_mf[lev]->setVal(0.0);
            damage_Psi_mf[lev]->setVal(0.0);
            damage_rate_mf[lev]->setVal(0.0);
        }

        // Step 9: cohesive-zone diagnostic fields. History (gb_delta_max) is
        // irreversible across timesteps; on AMR repair / regrid we preserve it
        // (similar to evolved damage). On a fresh init we zero everything.
        if (cohesive_enabled)
            ResetCohesiveDiagnostics(lev, /*reset_history=*/reset_damage);

        // Bulk effective-conductivity diagnostic (Step 6):
        //   k_eff = 0.5 · ( Σ_k φ_k κ_k  +  1 / Σ_k φ_k/κ_k )
        // where φ_k is the phase-k cell-volume fraction at this level. The
        // value is broadcast as a uniform per-cell field so it is trivially
        // exposed via yt without registering an integrated variable.
        std::vector<long> count(n_phases, 0);
        for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            // Host-side reduction on a single rank is fine for Step 6 (no
            // MPI tests). Skip the GPU reducer to keep this simple.
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                int p = static_cast<int>(ph(i, j, k));
                if (p >= 0 && p < n_phases) count[p]++;
            });
        }
        long total = 0;
        for (long c : count) total += c;
        Set::Scalar k_arith = 0.0, k_harm_inv = 0.0;
        if (total > 0)
        {
            for (int p = 0; p < n_phases; ++p)
            {
                Set::Scalar phi = static_cast<Set::Scalar>(count[p]) / total;
                k_arith    += phi * phases[p].kappa;
                if (phases[p].kappa > 0.0)
                    k_harm_inv += phi / phases[p].kappa;
            }
        }
        Set::Scalar k_harm = (k_harm_inv > 0.0) ? 1.0 / k_harm_inv : 0.0;
        Set::Scalar k_eff_value = 0.5 * (k_arith + k_harm);
        k_eff_mf[lev]->setVal(k_eff_value);
        if (announce)
            Util::Message(INFO, "microstructure k_eff (level ", lev, ") = ", k_eff_value,
                          "  (arith=", k_arith, ", harm=", k_harm, ")");

        // Step 15b: sample the per-cell Weibull flaw-length field. Runs
        // after is_grain_boundary_mf is populated above; default-off, only
        // active under spallation.model = sp_weibull + weibull.enabled = 1.
        InitializeFlaws(lev);
    }

    // Step 15b: deterministic per-cell Weibull flaw-length sampling.
    // a_f = a_scale * (-ln U)^(1/m), a 2-parameter Weibull with scale
    // a_scale and shape m; a_scale = a0_gb on is_grain_boundary cells,
    // a0_ig otherwise. The RNG is seeded from a splitmix64 hash of the
    // GLOBAL cell index (i, j, k) and weibull.seed, so flaw_a_mf is
    // bit-identical regardless of MPI rank count or BoxArray partition.
    // AMR regrid-repair of flaw_a_mf is a §15c deliverable; this resamples
    // deterministically on any call, which is single-level-safe.
    void InitializeFlaws(int lev)
    {
        if (!weibull_enabled) return;

        const Set::Scalar a0_gb = weibull_a0_gb;
        const Set::Scalar a0_ig = weibull_a0_ig;
        const Set::Scalar inv_m = 1.0 / weibull_m;
        const std::uint64_t seed =
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(weibull_seed));

        // Mesh-objective flaw assignment, split by physical population:
        //
        //   GB cells -> per-FACET. One Weibull draw per adjacent grain pair
        //     (gA, gB), keyed by a hash of the SORTED pair + seed. Grains are
        //     seed-fixed, so a facet draws the same flaw at any resolution and
        //     the facet population (hence the max surface flaw that sets onset,
        //     ΔT_onset ∝ a^(-1/2)) is mesh-invariant. A cell touching several
        //     grains takes the MAX (weakest-link) facet flaw. No vol_factor:
        //     the facet, not the cell, is the physical unit.
        //   IG cells -> per-CELL Weibull (hash of i,j,k + seed) scaled by the
        //     weakest-link size factor (V_cell/V_0)^(1/m); intragranular cracks
        //     are a volumetric population. V_0 (default (1 mm)^3) makes this
        //     exactly 1 at the calibration mesh, so 1 mm IG draws are identical.
        //
        // This removes the per-cell GB tag's mesh-dependence (tagged GB
        // fraction ~ dx/grain -> 0 as dx->0, which drove onset toward the wrong
        // all-IG limit). grain_id carries 1 ghost (ngrow=1) and is ghost-filled
        // below, so the 6-face neighbour grain lookup is valid across ranks;
        // the facet hash needs no MPI communication.
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar V_cell = DX[0] * DX[1] * DX[2];
        const Set::Scalar vol_factor = std::pow(V_cell / weibull_V0, inv_m);

        // Resolvability check for the connected_cluster lateral gate. A_crit is
        // a fixed physical area (m^2), so the gate itself is mesh-objective; but
        // when A_crit < one cell's lateral area dx*dy, a single firing cell
        // already clears it and the "patch vs. speck" filter is a silent no-op.
        // Warn once so a coarse run doesn't quietly bypass the gate.
        if (weibull_detach_mode == DetachMode::ConnectedCluster &&
            weibull_A_crit > 0.0 && !weibull_A_crit_warned)
        {
            const Set::Scalar A_cell = DX[0] * DX[1];
            if (weibull_A_crit < A_cell)
                Util::Warning(INFO, "weibull.A_crit (", weibull_A_crit,
                    " m^2) < cell lateral area dx*dy (", A_cell,
                    " m^2): the connected_cluster gate is a no-op at this "
                    "resolution (one firing cell clears it). Refine the mesh "
                    "or raise weibull.A_crit to filter specks.");
            weibull_A_crit_warned = true;
        }

        // grain_id ghosts must be current for the 6-face neighbour lookup.
        Util::RealFillBoundary(*grain_id_mf[lev], geom[lev]);
        const amrex::Box& domain = geom[lev].Domain();
        const auto dlo = domain.smallEnd();
        const auto dhi = domain.bigEnd();

        // Uniform U in (0, 1) from a 64-bit hash (53-bit mantissa).
        auto u_of = [] (std::uint64_t key) -> Set::Scalar {
            return (static_cast<Set::Scalar>(key >> 11) + 0.5)
                   * (1.0 / 9007199254740992.0);
        };
        auto u32 = [] (int v) -> std::uint64_t {
            return static_cast<std::uint64_t>(static_cast<std::uint32_t>(v));
        };

        for (amrex::MFIter mfi(*flaw_a_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& fa = (*flaw_a_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gb =
                (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gid =
                (*grain_id_mf[lev]).array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (gb(i, j, k) > 0.5)
                {
                    // GRAIN-BOUNDARY cell -> per-facet, weakest adjacent facet.
                    const int g_self = static_cast<int>(gid(i, j, k));
                    const int di[6] = {-1, 1, 0, 0, 0, 0};
                    const int dj[6] = { 0, 0,-1, 1, 0, 0};
                    const int dk[6] = { 0, 0, 0, 0,-1, 1};
                    Set::Scalar a_max = 0.0;
                    for (int n = 0; n < 6; ++n)
                    {
                        int ii = i + di[n], jj = j + dj[n], kk = k + dk[n];
                        if (ii < dlo[0] || ii > dhi[0] ||
                            jj < dlo[1] || jj > dhi[1] ||
                            kk < dlo[2] || kk > dhi[2]) continue;
                        int g_n = static_cast<int>(gid(ii, jj, kk));
                        if (g_n < 0 || g_n == g_self) continue;
                        // Sorted pair -> both sides of the facet draw alike.
                        int gA = g_self < g_n ? g_self : g_n;
                        int gB = g_self < g_n ? g_n : g_self;
                        std::uint64_t key = SplitMix64(seed);
                        key = SplitMix64(key ^ u32(gA));
                        key = SplitMix64(key ^ u32(gB));
                        Set::Scalar a_f =
                            a0_gb * std::pow(-std::log(u_of(key)), inv_m);
                        if (a_f > a_max) a_max = a_f;
                    }
                    if (a_max <= 0.0)
                    {
                        // A GB-tagged cell should always have a different-grain
                        // face neighbour; if not (shouldn't happen), fall back
                        // to a per-cell GB draw so flaw_a is never left at 0.
                        std::uint64_t key = SplitMix64(seed);
                        key = SplitMix64(key ^ u32(i));
                        key = SplitMix64(key ^ u32(j));
                        key = SplitMix64(key ^ u32(k));
                        a_max = a0_gb * std::pow(-std::log(u_of(key)), inv_m);
                    }
                    fa(i, j, k) = a_max;
                }
                else
                {
                    // INTRAGRANULAR cell -> per-cell volumetric (size effect).
                    std::uint64_t key = SplitMix64(seed);
                    key = SplitMix64(key ^ u32(i));
                    key = SplitMix64(key ^ u32(j));
                    key = SplitMix64(key ^ u32(k));
                    fa(i, j, k) =
                        a0_ig * vol_factor * std::pow(-std::log(u_of(key)), inv_m);
                }
            });
        }
    }

    // Conserved enthalpy and phase fractions (own to MMWSpalling; T comes from
    // HeatConduction's temp_mf / temp_old_mf).
    Set::Field<Set::Scalar> H_mf, H_old_mf;
    Set::Field<Set::Scalar> Lambda_S_mf, Lambda_L_mf, Lambda_V_mf;
    Set::Field<Set::Scalar> surface_mf; // 1.0 on surface cells, 0.0 otherwise

    Set::Scalar therm_alpha = NAN; // κ/(ρ c_p), cached from hc.heat.alpha
    Set::Scalar T_ref = 0.0;
    Set::Scalar T_m   = 1.0e30;
    Set::Scalar T_v   = 1.0e30;
    Set::Scalar lh_m  = 0.0;
    Set::Scalar lh_v  = 0.0;

    Set::Scalar eps_rad      = 0.0;     // emissivity below T_eps_break  (-)
    Set::Scalar eps_rad_hi   = 0.0;     // emissivity above T_eps_break  (-)
    Set::Scalar T_eps_break  = 1.0e30;  // breakpoint for ε(T)           (K)
    Set::Scalar sigma_SB = 0.0; // Stefan-Boltzmann constant       (W/m²K⁴)
    Set::Scalar h_conv   = 0.0; // convective heat-transfer coeff. (W/m²K)
    Set::Scalar T_amb    = 0.0; // ambient temperature             (K)

    Numeric::Material::Material *material = nullptr; // optional T-dependent props

    // H↔T lookup table for the global `material` path. Step 6b moved the
    // implementation into a reusable Numeric::Material::Table; the
    // microstructure path uses one Table per phase (see phase_tables).
    Numeric::Material::Table material_table;

    enum class SolverKind { Explicit, Implicit };
    SolverKind  solver_kind        = SolverKind::Explicit;
    Set::Scalar implicit_tol_rel   = 1.0e-10;
    Set::Scalar implicit_tol_abs   = 0.0;
    int         implicit_max_iter  = 100;
    int         implicit_verbose   = 0;

    std::vector<Set::Scalar> alpha; // thermoelastic eigenstrain coefficient per component

    Numeric::MMWBeam beam; // Beer-Lambert + Gaussian volumetric source (step 3)

    // E1: beam-deposition closure counter (BeamEnergyBalance), written to
    // thermo.dat as non-extensive vars (reduced by BeamEnergyBalance itself).
    // Powers [W] are per step; energies [J] and the violation count are
    // cumulative. Microstructure path, level 0 only.
    Set::Scalar beam_P_inc          = 0.0;
    Set::Scalar beam_P_dep          = 0.0;
    Set::Scalar beam_P_trans        = 0.0;
    Set::Scalar beam_E_inc          = 0.0;
    Set::Scalar beam_E_dep          = 0.0;
    Set::Scalar beam_E_trans        = 0.0;
    Set::Scalar beam_closure_err    = 0.0;
    Set::Scalar beam_invariant_viol = 0.0;
    long        beam_viol_steps     = 0;

    // S1: surface.follow_mask applies the convective_flame patch and the
    // surface losses at the `removed`-mask column top (ColumnTopSolid)
    // instead of the fixed domain top / the phi-window surface_mf.
    bool surface_follow_mask = false;

    // S1: whole-domain energy ledger (EnergyLedger), thermo.dat vars.
    // E_H [J] is sum(H*V) over all level-0 cells (void included), E_src [J]
    // the cumulative applied energy, P_* [W] per step. The compensated-sum
    // state (ledger_E_H0_s/c) is internal.
    bool        ledger_enabled               = false;
    // C1: pinned-surface closure state (lev 0, ncols-sized host vectors,
    // identical on every rank). pin_T = top-cell T at the column's last spall
    // firing (0 = never fired), pin_t = that time, pin_rhoCp = rho*Cp of the
    // top-cell phase then. pin_T_eff = the pin used by this step's H update
    // and ledger (0 = face form), built at the start of AdvanceMicrostructure.
    // NOT checkpointed: a restart falls back to face until each column fires.
    std::vector<Set::Scalar> pin_T, pin_t, pin_rhoCp, pin_T_eff;
    Set::Scalar pinned_cols_pinned = 0.0; // thermo: in-patch top cells pinned
    Set::Scalar pinned_cols_face   = 0.0; // thermo: in-patch top cells on face
    Set::Scalar pinned_q_mean      = 0.0; // thermo: mean in-patch q_robin [W/m^2]
    // D2a: why a pinned-run in-patch top cell took the face branch. pin_reason
    // per column from BuildPinEffective: 0 never fired, 1 idle-expired,
    // 2 q_pin <= 0, 3 pin active (face then means T_f <= T_pin, the min rule).
    std::vector<int> pin_reason;
    Set::Scalar pinned_cols_face_unset   = 0.0;
    Set::Scalar pinned_cols_face_idle    = 0.0;
    Set::Scalar pinned_cols_face_qneg    = 0.0;
    Set::Scalar pinned_cols_face_minrule = 0.0;
    // D2a: per-column flame (h, T_gas) and nozzle stand-off s for the current
    // AdvanceMicrostructure call (BuildFlameColumns). Constants path: h_conv /
    // T_flame everywhere. Read by the kernel, BuildPinEffective, EnergyLedger
    // and PatchDiagnostics.
    std::vector<Set::Scalar> flame_h, flame_Tg, flame_s;
    Set::Scalar patch_min_standoff = 0.0; // thermo: min s over in-patch live columns [m]
    Set::Scalar patch_P_robin      = 0.0; // thermo: sum q_robin*dx*dy over patch top cells [W]
    // D2a2: jet_closure = enthalpy diagnostics (JetEnthalpyMarch).
    Set::Scalar jet_P_face    = 0.0; // thermo: enthalpy taken from the jet stream [W]
    Set::Scalar jet_P_cap     = 0.0; // thermo: m*cp*(T_stag - T_ent) [W]
    Set::Scalar jet_T_stag    = 0.0; // thermo: gas temperature entering the axis bin [K]
    Set::Scalar jet_s_c       = 0.0; // thermo: mean stand-off of the axis bin [m]
    Set::Scalar jet_P_decay   = 0.0; // thermo: nozzle excess discarded by the decay line [W]
    Set::Scalar jet_T_exhaust = 0.0; // thermo: T_gas after the outermost bin [K]
    Set::Scalar jet_P_exhaust = 0.0; // thermo: m*cp*(T_exhaust - T_ent) [W]
    Set::Scalar jet_r_reach   = 0.0; // thermo: max column r in a bin with P_bin > 0 [m]
    // D2b: exhaust recirculation (jet_T_ent_mode = exhaust). jet_T_rec = the
    // entrained-gas temperature used this step (= the previous step's
    // jet_T_exhaust; jet_T_ent on the first step). Host state, NOT checkpointed.
    Set::Scalar jet_T_rec     = 0.0; // thermo [K]
    Set::Scalar jet_P_recirc  = 0.0; // thermo: (m - mdot)*cp*(T_rec - T_ent) [W]
    bool        jet_rec_init  = false;
    Set::Scalar jet_T_rec_next = 0.0;
    // D2b: burner-on-feet state (nozzle_descent = feet). nozzle_z is host
    // state, NOT checkpointed (a restart re-lands the burner on the feet).
    bool        feet_init           = false;
    Set::Scalar nozzle_z            = 0.0; // thermo: nozzle exit height [m]
    Set::Scalar foot_z              = 0.0; // thermo: foot (tripod) height [m]
    Set::Scalar foot_cols           = 0.0; // thermo: live annulus columns
    Set::Scalar foot_carry_cols     = 0.0; // thermo: columns at or above their pad height
    Set::Scalar foot_stall_time     = 0.0; // thermo: simulated time since z_n last decreased [s]
    long        jet_profile_last = -1; // last jet_profile_interval index written
    bool        jet_profile_header_written = false;

    bool        ledger_started               = false;
    Set::Scalar ledger_E_H0_s                = 0.0;
    Set::Scalar ledger_E_H0_c                = 0.0;
    Set::Scalar ledger_E_H                   = 0.0;
    Set::Scalar ledger_E_src                 = 0.0;
    Set::Scalar ledger_err                   = 0.0;
    Set::Scalar ledger_P_beam                = 0.0;
    Set::Scalar ledger_P_robin               = 0.0;
    Set::Scalar ledger_P_loss                = 0.0;
    Set::Scalar ledger_surface_missing_cols  = 0.0;

    // Step 6d / Step 15a: surface patch driver. A circular patch on the +z
    // face applies one of two heat-load modes inside the patch footprint:
    //   mode = prescribed_T (default, Step 6d) — holds T = T_expr(t) and
    //     injects the Dirichlet face flux q = kappa*(T_patch - T_cell)/dz as a
    //     divergence contribution kappa*(T_patch - T_cell)/dz^2 on the
    //     in-patch top row.
    //   mode = convective_flame (Step 15a) — applies a Robin flame-jet BC
    //     q = h_fl*(T_flame - T_surface), injected as the divergence
    //     contribution h_fl*(T_flame - T_cell)/dz on the in-patch top row.
    // Outside the patch the existing zhi BC continues to drive the face flux
    // (typically Neumann/adiabatic). Disabled by default; the default mode
    // keeps every Step 6d/7b/8b/11b/14 input byte-for-byte.
    struct SurfacePatch
    {
        static constexpr const char *name = "surface_patch";

        enum class Mode { PrescribedTemperature, ConvectiveFlame };
        enum class RobinForm { Cell, Face, Pinned };

        bool enabled = false;
        Mode mode = Mode::PrescribedTemperature;
        RobinForm robin_form = RobinForm::Cell;
        Set::Scalar pinned_idle_cycles = 2.0; // C1: robin_form = pinned only
        Set::Scalar x0 = 0.0, y0 = 0.0, radius = 0.0;
        Set::Scalar T_ext   = 0.0;
        Set::Scalar h_conv  = 0.0;  // convective_flame: flame coefficient h_fl
        Set::Scalar T_flame = 0.0;  // convective_flame: flame temperature
        amrex::Parser            T_parser;
        amrex::ParserExecutor<1> T_f;
        // D2a: spatial flame expressions (x, y, r, s, t); use_expr = both given.
        bool use_expr = false;
        Set::Scalar nozzle_z0 = 0.0, nozzle_feed = 0.0, nozzle_collision_radius = 0.0;
        amrex::Parser            h_parser, Tg_parser;
        amrex::ParserExecutor<5> h_f, Tg_f;
        // D2a2: jet_closure = enthalpy (T_flame_expr not read, Tg_f uncompiled).
        bool jet_enthalpy = false;
        Set::Scalar jet_mdot = 0.0, jet_T_nozzle = 0.0, jet_cp = 0.0;
        Set::Scalar jet_T_ent = 293.15, jet_dr = 0.0, jet_profile_interval = 0.0;
        Set::Scalar jet_D = 0.0, jet_core_length = 5.0;
        bool jet_stag_decay = true, jet_entrained_mass = false;
        bool jet_T_ent_exhaust = false;     // D2b: jet_T_ent_mode = exhaust
        // D2b: burner-on-feet descent (nozzle_descent = feet).
        enum class FootRule { Pads, Max, Mean };
        bool feet = false;
        FootRule foot_rule = FootRule::Pads;
        Set::Scalar foot_r_inner = 0.0, foot_r_outer = 0.0, foot_standoff = 0.0;
        int foot_npads = 3;
        Set::Scalar foot_pad_quantile = 0.9;
        Set::Scalar nozzle_feed_max = std::numeric_limits<Set::Scalar>::infinity();

        // prescribed_T only: patch surface temperature at time t. Not valid
        // under convective_flame (T_f is left uncompiled in that mode).
        Set::Scalar T_at(Set::Scalar t) const { return T_f(t); }

        static void Parse(SurfacePatch &v, IO::ParmParse &pp)
        {
            int en = 0;
            pp.query_default("enabled", en, 0);
            v.enabled = (en != 0);
            if (!v.enabled) return;

            // Legacy `type` key — always consumed so the strict parser does
            // not choke on existing Step 6d-style inputs.
            std::string type;
            pp.query_default("type", type, std::string("prescribed_temperature"));
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
            if (v.mode == Mode::PrescribedTemperature &&
                type != "prescribed_temperature")
                Util::Abort(INFO, "surface_patch.type must be "
                                  "'prescribed_temperature', got: ", type);

            // S1: Robin discretisation for convective_flame (microstructure
            // path only; see SurfaceCellFlux). `cell` (default) applies
            // h*(T_flame - T_cell) and cell-centred losses; `face` solves the
            // top-face balance for T_f and injects (2k/dz)*(T_f - T_cell).
            std::string robin_str;
            pp.query_default("robin_form", robin_str, std::string("cell"));
            if      (robin_str == "cell")   v.robin_form = RobinForm::Cell;
            else if (robin_str == "face")   v.robin_form = RobinForm::Face;
            else if (robin_str == "pinned") v.robin_form = RobinForm::Pinned;
            else Util::Abort(INFO, "surface_patch.robin_form must be 'cell', 'face' "
                                   "or 'pinned', got: ", robin_str);
            if (v.robin_form == RobinForm::Face && v.mode != Mode::ConvectiveFlame)
                Util::Abort(INFO, "surface_patch.robin_form = face requires "
                                  "surface_patch.mode = convective_flame");
            // C1: pinned-surface closure (see SurfaceCellFlux). A column that
            // has not fired for pinned_idle_cycles x (its own cell heating time)
            // falls back to the face form.
            if (v.robin_form == RobinForm::Pinned)
            {
                if (v.mode != Mode::ConvectiveFlame)
                    Util::Abort(INFO, "surface_patch.robin_form = pinned requires "
                                      "surface_patch.mode = convective_flame");
                pp.query_default("pinned_idle_cycles", v.pinned_idle_cycles, 2.0);
                if (!(std::isfinite(v.pinned_idle_cycles) && v.pinned_idle_cycles >= 0.0))
                    Util::Abort(INFO, "surface_patch.pinned_idle_cycles must be finite "
                                      "and >= 0, got: ", v.pinned_idle_cycles);
            }

            std::string face;
            pp.query_default("face", face, std::string("zhi"));
            if (face != "zhi")
                Util::Abort(INFO, "surface_patch.face must be 'zhi' (only +z "
                                  "patch implemented), got: ", face);
            pp.query_required("x0",     v.x0);
            pp.query_required("y0",     v.y0);
            pp.query_required("radius", v.radius);

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

            // T_ext / h_conv / T_flame. Under convective_flame, h_conv is the
            // flame heat-transfer coefficient h_fl and T_flame the flame
            // temperature: q = h_fl*(T_flame - T_surface). Under prescribed_T
            // they are reserved/unused.
            pp.query_default("T_ext",   v.T_ext,   0.0);

            // D2a: spatial flame, h_expr / T_flame_expr in (x, y, r, s, t) [SI],
            // r = distance of the column centre from (x0, y0), s = stand-off of
            // the nozzle z_n(t) = nozzle_z0 - nozzle_feed*t above the column's
            // mask top face. convective_flame only; both expressions or both
            // constants, never a mix. (surface.follow_mask = 1 is checked in
            // MMWSpalling::Parse.)
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
                if      (tem == "fixed")   v.jet_T_ent_exhaust = false;
                else if (tem == "exhaust") v.jet_T_ent_exhaust = true;
                else Util::Abort(INFO, "surface_patch.jet_T_ent_mode must be 'fixed' or "
                                       "'exhaust', got: ", tem);
                if (v.jet_T_ent_exhaust && !v.jet_entrained_mass)
                    Util::Abort(INFO, "surface_patch.jet_T_ent_mode = exhaust requires "
                                      "surface_patch.jet_entrained_mass = 1");
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
                    std::string rule;
                    pp.query_default("foot_rule", rule, std::string("pads"));
                    if      (rule == "pads") v.foot_rule = FootRule::Pads;
                    else if (rule == "max")  v.foot_rule = FootRule::Max;
                    else if (rule == "mean") v.foot_rule = FootRule::Mean;
                    else Util::Abort(INFO, "surface_patch.foot_rule must be 'pads', 'max' or "
                                           "'mean', got: ", rule);
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
                // which must stay finite and > 0.
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

    // Step 5: Voronoi mineral microstructure. `phase_mf` is a cell-centred
    // per-cell mineral phase ID (Set::Scalar storage; integers cast to float).
    // The property fields cache the heterogeneous κ / β / E / ρ / Cp lookup so
    // later steps can read them as ordinary AMReX MultiFabs without
    // re-querying the phase table on every kernel.
    struct MineralPhase
    {
        std::string name;
        Set::Scalar fraction = 0.0;
        Set::Scalar kappa = 0.0;
        Set::Scalar beta  = 0.0;
        Set::Scalar E     = 0.0;
        Set::Scalar mu    = 0.0;
        Set::Scalar rho   = 0.0;
        Set::Scalar Cp    = 0.0;
        Set::Scalar alpha_atten = 0.0;
        Set::Scalar A_D   = 0.0;
        // Step 6b: phase-change scalars used to build the per-phase
        // H<->T table. Defaults disable phase change.
        Set::Scalar T_ref     = 273.15;
        Set::Scalar T_m       = 1.0e30;
        Set::Scalar fus_width = 0.0;
        Set::Scalar L_m       = 0.0;
        Set::Scalar T_vap_lo  = 1.0e30;
        Set::Scalar T_vap_hi  = 1.0e30;
        Set::Scalar L_v       = 0.0;
        bool alpha_beta_enabled = false;
        Set::Scalar alpha_beta_T0 = 846.15;
        Set::Scalar alpha_beta_width = 10.0;
        Set::Scalar alpha_beta_dV = 0.0063;
        Set::Scalar alpha_beta_linear_strain = 0.0;
    };
    enum class MicroMode { Voronoi, Expression };
    bool microstructure_enabled = false;
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
    // index-aligned with `phases`. Built once in Initialize().
    std::vector<Numeric::Material::Constant> phase_materials;
    std::vector<Numeric::Material::Table>    phase_tables;
    std::vector<int>          grain_phase_ids;
    std::vector<Set::Scalar>  grain_centres; // packed (x,y,z) per grain

    Set::Field<Set::Scalar> phase_mf;
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
    // heterogeneous thermoelastic update. Cell-centred, 1-ghost; sampled by
    // CellToNodeAverage in UpdateModel.
    Set::Field<Set::Scalar> mu_phase_mf, T_ref_phase_mf;
    // Optional quartz alpha-beta transformation strain and the total scalar
    // free strain beta*(T-T_ref) + transition strain. Cell-centred, 1 ghost.
    Set::Field<Set::Scalar> alpha_beta_linear_strain_mf;
    Set::Field<Set::Scalar> alpha_beta_T0_mf;
    Set::Field<Set::Scalar> alpha_beta_width_mf;
    Set::Field<Set::Scalar> thermal_free_strain_mf;

    // Step 6: passive damage, damage-modified κ, GB conductance, k_eff.
    Set::Field<Set::Scalar> D_mf;
    Set::Field<Set::Scalar> kappa_eff_mf; // κ_phase · exp(−α_D · D)
    Set::Field<Set::Scalar> k_eff_mf;     // bulk arithmetic-harmonic average
    Set::Scalar kappa_damage_alpha = 0.0;
    Set::Scalar h_gb_conductance   = 0.0;
    IC::IC<Set::Scalar> *ic_damage = nullptr; // optional, defaults to D ≡ 0

    // Step 8: Drucker-Prager damage evolution. ALAMO mechanics uses positive
    // tensile stress; compression is negative in I1.
    bool damage_enabled = false;
    Set::Scalar damage_phi_deg = 0.0;
    Set::Scalar damage_phi_rad = 0.0;
    Set::Scalar damage_cohesion = 1.0;
    Set::Scalar damage_alpha_dp = 0.0;
    Set::Scalar damage_k_dp = 1.0;
    Set::Scalar damage_psi0 = 0.0;
    Set::Scalar damage_n = 1.0;
    Set::Scalar damage_m = 1.0;
    Set::Scalar damage_gb_multiplier = 1.0;
    Set::Scalar damage_stiffness_floor = 1.0e-6;
    Set::Field<Set::Scalar> damage_F_mf;
    Set::Field<Set::Scalar> damage_Psi_mf;
    Set::Field<Set::Scalar> damage_rate_mf;

    // Default-off deterministic driver for the Step 8 regression. Production
    // runs leave this disabled and use the nodal mechanics stress field.
    bool damage_prescribed_stress_enabled = false;
    amrex::Parser damage_sxx_parser, damage_syy_parser, damage_szz_parser;
    amrex::Parser damage_sxy_parser, damage_syz_parser, damage_szx_parser;
    amrex::ParserExecutor<4> damage_sxx_f, damage_syy_f, damage_szz_f;
    amrex::ParserExecutor<4> damage_sxy_f, damage_syz_f, damage_szx_f;

    // Step 9: bilinear grain-boundary cohesive-zone law. Default-off; only
    // active when cohesive.enabled = 1 with a valid microstructure.
    // Calibrated from any of {K_n, delta_c, delta_max} or {K_n, G_c, f_t}
    // pairs. State fields are populated only on cells flagged as
    // is_grain_boundary; everywhere else they are forced to zero. Step 9
    // computes the law and tracks irreversible history; coupling the
    // resulting traction back into the elastic operator is deferred.
    bool cohesive_enabled = false;
    Set::Scalar cohesive_K_n       = 0.0;
    Set::Scalar cohesive_delta_c   = 0.0;
    Set::Scalar cohesive_delta_max = 0.0;
    bool cohesive_prescribed_enabled = false;
    amrex::Parser cohesive_delta_parser;
    amrex::ParserExecutor<4> cohesive_delta_f;
    Set::Field<Set::Scalar> gb_delta_mf;            // current opening
    Set::Field<Set::Scalar> gb_delta_max_mf;        // irreversible history
    Set::Field<Set::Scalar> gb_traction_mf;         // current normal traction
    Set::Field<Set::Scalar> gb_cohesive_damage_mf;  // 1 - secant/initial

    // Step 10: surface advancement and spall detachment. Default-off.
    // `phi = z_surface(x,y) - z_cell` is a level-set state; cells with
    // phi >= 0 are material, phi < 0 are void. `removed` is the binary
    // void mask. `spall_event`, `spall_thickness`, and `RoP_spall` are
    // per-step diagnostics reset at the start of each spall update.
    bool spall_enabled = false;
    Set::Scalar spall_damage_threshold = 0.99;
    Set::Scalar spall_xi_s = 8.0;
    Set::Scalar spall_C_h  = 1.0;
    bool spall_alpha_override = false;
    Set::Scalar spall_alpha = 0.0;
    bool spall_t_spall_override = false;
    Set::Scalar spall_t_spall = 0.0;
    bool spall_prescribed_ratio_enabled = false;
    amrex::Parser spall_ratio_parser;
    amrex::ParserExecutor<4> spall_ratio_f;
    // Step 14b: how `D` is sampled per column for the
    // `spall_damage_threshold` check. Default `TopCell` keeps Step 10/11
    // behaviour byte-identical. `ColumnMax` scans the entire column for
    // max D; `ShallowMax` scans cells with depth <= `spall_shallow_depth`
    // below the current top (`k_top` of the column). Stress is still
    // sampled at `k_top` regardless of mode — release geometry is
    // unchanged.
    enum class SpallSampleMode { TopCell, ColumnMax, ShallowMax };
    SpallSampleMode spall_sample_mode  = SpallSampleMode::TopCell;
    Set::Scalar     spall_shallow_depth = 0.006;

    // Flake coherence length [m] (rock-mechanics cavity-wall coherence). Default
    // 0 = off. When > 0, a firing column may recede at most this far below its
    // shallowest live 4-neighbour, so the spall cavity floor has a bounded wall
    // slope: a smooth bowl can descend gradually (<= L_coh per lateral cell) but
    // a free-standing high-aspect-ratio needle cannot form (the runaway column
    // is capped until its neighbours catch up). Decomposition-independent
    // (operates on the gathered per-column surface heights).
    Set::Scalar spall_flake_coh_len = 0.0;

    // Step 16: Sp-vs-depth h_spall scan parameters (sp_weibull branch only).
    // `spall.h_spall_max` caps the downward column scan in metres (safety
    // stop; in practice the K_I < K_Ic crossing always fires first). The
    // segment count `spall.h_spall_n_segments` falls back to `sp_n_segments`
    // when set to 0 (default), so the depth scan uses the same Simpson
    // resolution as `UpdateSpAfterMechanics`'s surface Sp.
    Set::Scalar spall_h_spall_max         = 0.010;
    int         spall_h_spall_n_segments  = 0;
    // Stress source for the Sp = K_I/K_Ic criterion (both the Sp_field surface
    // diagnostic/gate AND the h_spall depth scan). Default `mechanics` samples
    // the elastic MLMG solve's stress_mf. `local_thermoelastic` reconstructs
    // sigma_xx algebraically from the temperature field via the 1-D-confined
    // relation sigma_xx = -E*beta*(T - T_ref)/(1 - nu) (phases[0]), so the
    // criterion runs with mechanics disabled (el.type = disable) — no MLMG
    // solve, no convergence cliff, no max_level = 0 resolution ceiling. Valid
    // for the brittle-spallation regime where the constrained heated skin is
    // the dominant stress driver (Kant/von Rohr formulation). Loses grain-scale
    // stress heterogeneity and cavity stress redistribution; see
    // mlmg-shelved-local-thermoelastic-notes.md.
    bool sp_local_thermoelastic = false;
    // Step 20: surface-normal (slope-aware) spall criterion + removal. When
    // true, the K_I depth scan runs in the local surface-normal frame and a
    // flake of normal thickness h recedes its column by h/cos(theta) -- the
    // graph form of level-set normal motion. Default false (vertical-only,
    // bit-identical to all pre-Step-20 behaviour).
    bool spall_surface_normal = false;
    // P1: melt-aware spallation. When melt_aware is true, the sp_weibull
    // driving stress is scaled by MeltSkeleton(Lambda_L) in k_i_at_zcrack so a
    // molten cell carries no crack-driving stress (the depth scan then breaks
    // at the solid-liquid interface on its own). Default false => f == 1
    // everywhere (byte-identical).
    bool        melt_aware       = false;
    Set::Scalar melt_lambda_crit = 0.4;
    // Load-bearing skeleton fraction f(Lambda_L) = max(0, 1 - Lambda_L/
    // lambda_crit): 1 in solid, 0 at the rheologically critical melt fraction
    // (rigidity percolation) and above. Multiplies the elastic driving stress
    // (equivalently E_eff = f*E) so a melt cannot store the stress that drives
    // a crack. lambda_crit is validated in (0, 1] at Parse.
    AMREX_FORCE_INLINE
    Set::Scalar MeltSkeleton(Set::Scalar LL) const
    {
        Set::Scalar f = 1.0 - LL / melt_lambda_crit;
        if (f < 0.0) f = 0.0;
        if (f > 1.0) f = 1.0;
        return f;
    }
    // Floor on cos(theta) when surface_normal is on: walls steeper than
    // ~81 deg keep a finite normal frame instead of diverging. The column
    // (height-function) representation cannot represent overhangs anyway.
    static constexpr Set::Scalar SP_SN_COS_MIN = 0.15;
    // Cell-centred per-column h_spall diagnostic, written at (i, j, k_top)
    // on the step the column fires under sp_weibull. Registered only when
    // spallation_model == SpWeibull && spall_enabled.
    Set::Field<Set::Scalar> h_spall_field_mf;
    // Step 16b: time-integrated 3-D crack-event count field. Cell-centred,
    // binary 0/1, edge-triggered (option (ii) in
    // rossi-validation-diagnostic-design.md §3 — freeze re-firing without
    // touching sigma). Each (i, j, k) cell records whether K_I(z=depth_of_k;
    // a_f at column top) >= K_Ic(T at cell) ever held at any timestep,
    // with `a_f` fixed at the surface cell's flaw_a per design-note §5/Q2
    // option (a). Once a cell flips 0 -> 1, it stays at 1 (freeze guard).
    // NOT reset per step — the field accumulates across the whole run.
    // Registered only when spallation_model == SpWeibull && spall_enabled.
    Set::Field<Set::Scalar> crack_event_count_mf;

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> removed_mf;
    Set::Field<Set::Scalar> spall_event_mf;
    Set::Field<Set::Scalar> spall_thickness_mf;
    Set::Field<Set::Scalar> RoP_spall_mf;

    // Step 11: vaporisation removal. Default-off and shares phi / removed
    // with the spall path. RoP_vap follows main.tex:
    //   RoP_vap = (1 - R) * P0(t) / (rho * L_tot * A_beam)
    //   L_tot   = Cp * (T_vap_lo - T0) + L_m + L_v
    // The trigger is H >= H_vap = phase_table.H_of_T(T_vap_lo). If both
    // spall and vapor want the same column in the same step, the larger
    // candidate thickness wins and `regime` records the winning mode.
    bool vapor_enabled = false;
    bool vapor_A_beam_override = false;
    Set::Scalar vapor_A_beam = 0.0;
    bool vapor_T0_override = false;
    Set::Scalar vapor_T0 = 0.0;
    bool vapor_prescribed_rop_enabled = false;
    amrex::Parser vapor_rop_parser;
    amrex::ParserExecutor<4> vapor_rop_f;

    // Step 19: mechanical bit feed (drawworks). When bit.enabled, a rigid flat
    // bit face descends along the beam axis at bit.feed_rate and removes any
    // material above it within bit.radius -- forcing a flat, advancing,
    // cylindrical floor that keeps refreshing the thin thermal skin so
    // spallation does not stall into a self-limiting cone. Spallation/vapor may
    // still remove material AHEAD of (below) the bit; the surface recession per
    // column is max(spall, vapor, bit). The bit cut is NOT coherence-capped.
    bool        bit_enabled    = false;
    Set::Scalar bit_feed_rate  = 0.0;   // [m/s] downward feed (ROP)
    Set::Scalar bit_radius     = 0.0;   // [m] bit footprint radius
    Set::Scalar bit_x0         = 0.0;   // [m] bit axis x (default beam.x0)
    Set::Scalar bit_y0         = 0.0;   // [m] bit axis y (default beam.y0)
    Set::Scalar bit_z0         = 0.0;   // [m] initial bit floor (default beam.z0)
    Set::Scalar bit_start_time = 0.0;   // [s] time the feed begins
    Set::Field<Set::Scalar> vapor_event_mf;
    Set::Field<Set::Scalar> vapor_thickness_mf;
    Set::Field<Set::Scalar> RoP_vap_mf;
    Set::Field<Set::Scalar> regime_mf;
    Set::Field<Set::Scalar> RoP_mf;
    // Step 17: cell-resolved per-firing-cell regime diagnostic under
    // sp_weibull. {0 none, 1 spall, 2 vap}. Registered only when
    // spallation_model == SpWeibull && (spall_enabled || vapor_enabled).
    // Distinct from the scalar `regime_mf` (which stays for damage_law
    // backward compat).
    Set::Field<Set::Scalar> regime_field_mf;

    bool removal_enabled() const { return spall_enabled || vapor_enabled; }

    // Step 12: spallation-model toggle. `damage_law` (default) preserves
    // every Step 8-11b behavior; `sp_weibull` activates the LEFM
    // Sp = K_I/K_Ic(T) onset criterion. Mixing is enforced as a hard parser
    // error in ParseSpallationModel.
    enum class SpallationModel { DamageLaw, SpWeibull };
    SpallationModel spallation_model = SpallationModel::DamageLaw;

    // Sp criterion parameters (sp_weibull branch only). Deterministic
    // surface flaw length a0 (m); Weibull randomization arrives in Step 15.
    Set::Scalar sp_a0 = 0.0;
    int         sp_n_segments = 64;     // composite Simpson segments

    // Lithostatic confining pressure (Pa, Kant sign convention: positive =
    // compressive). Applied as +p·nu/(1-nu) offset in both sigma_xx lambdas
    // (UpdateSpAfterMechanics and the k_i_at_zcrack depth scan), NOT as an
    // elastic-solver BC. Default 0 short-circuits to byte-identical.
    Set::Scalar confining_p = 0.0;

    enum class KIcLaw { NasseriTable, Linear };
    KIcLaw      kic_law         = KIcLaw::NasseriTable;
    Set::Scalar kic_scale       = 1.0;  // multiplicative scale on Nasseri table
    Set::Scalar kic_K_Ic0       = 0.0;  // linear law: K_Ic(T_ref)
    Set::Scalar kic_T_ref       = 0.0;  // linear law: reference temperature (K)
    Set::Scalar kic_T_melt      = 0.0;  // linear law: melting span (K)

    // Optional prescribed sigma_xx(x, y, z, t) override for verification.
    // Mirrors the damage.prescribed_stress.* pattern used in Step 8. When
    // enabled, K_I integration uses the parser expression instead of the
    // mechanics-produced stress field, isolating the K_I integrator and
    // K_Ic(T) interpolator from FEM-stress-field noise.
    bool sp_prescribed_sigma_enabled = false;
    amrex::Parser            sp_prescribed_sigma_parser;
    amrex::ParserExecutor<4> sp_prescribed_sigma_f;

    // Optional prescribed T_surface(t) override for K_Ic look-up under
    // verification. When unset (production path), K_Ic uses `temp_mf` at
    // the top cell. Used together with sp_prescribed_sigma to give a clean
    // K_I + K_Ic closed-form cross-check without depending on the timing
    // of the H/T swap relative to plotfile output.
    bool sp_prescribed_T_surface_enabled = false;
    amrex::Parser            sp_prescribed_T_surface_parser;
    amrex::ParserExecutor<1> sp_prescribed_T_surface_f;

    Set::Field<Set::Scalar> Sp_field_mf;

    // Step 15b: Weibull per-cell flaw-length distribution (sp_weibull branch
    // only). Default-off; `weibull.enabled = 1` requires sp_weibull +
    // microstructure (the latter is already enforced by the sp_weibull +
    // !microstructure_enabled abort). `flaw_a_mf` is cell-centred: one
    // Weibull sample per cell, characteristic length keyed on
    // is_grain_boundary_mf. The long plan §15b says "face-centred", but
    // ALAMO's microstructure topology fields are all cell-centred and Sp is
    // evaluated per surface cell, so a cell-centred flaw field is the
    // faithful translation of "per grain-boundary face".
    bool        weibull_enabled = false;
    Set::Scalar weibull_a0_gb   = 0.0;  // Weibull scale, grain-boundary cells (m)
    Set::Scalar weibull_a0_ig   = 0.0;  // Weibull scale, intragranular cells (m)
    Set::Scalar weibull_m       = 15.0; // Weibull shape (modulus)
    // Mesh-objectivity reference volume for the weakest-link size effect: the
    // per-cell flaw scale is multiplied by (V_cell/V_0)^(1/m). FIXED physical
    // constant (NOT the run's cell size) so runs at different resolutions share
    // one calibration. Default (1 mm)^3 = the design mesh, where the factor = 1.
    Set::Scalar weibull_V0      = 1.0e-9; // m^3
    int         weibull_seed    = 0;
    Set::Field<Set::Scalar> flaw_a_mf;

    // Step 15c: connected-component spall labeller + detach mode gate.
    // `DetachMode::PerFace` (default) preserves §15b behaviour byte-identical:
    // every Sp >= 1 cell shows up as firing, the labeller writes a positive
    // cluster ID into Sp_cluster_id_mf, and no cluster is zeroed out.
    // `ConnectedCluster` additionally zeroes out clusters whose area is
    // strictly below `weibull_A_crit` (m^2), in Sp_cluster_id_mf only —
    // Sp_field_mf is never modified by the labeller.
    enum class DetachMode { PerFace, ConnectedCluster };
    DetachMode  weibull_detach_mode = DetachMode::PerFace;
    Set::Scalar weibull_A_crit      = 0.0;  // m^2, only used in ConnectedCluster
    bool        weibull_A_crit_warned = false; // one-shot resolvability warning
    Set::Field<Set::Scalar> Sp_cluster_id_mf;
    // Per-step cluster diagnostics CSV: written by the I/O processor only
    // when weibull is enabled and a cluster fires. Tracks the path so a
    // single trailing-newline header is written exactly once per output
    // directory (re-opening with std::ios::app keeps existing data on
    // restart-style runs).
    bool sp_clusters_csv_header_written = false;
    // CSV emission toggles (default OFF — most runs are behaviour studies that
    // only need plotfiles). Enable per-input where a test consumes the file:
    // `spall.clusters_csv` (kant A/B cluster onset) and `spall.h_col_events_csv`
    // (rossi depth-profile). Avoids littering output dirs with unused CSVs.
    bool sp_emit_clusters_csv = false;
    bool sp_emit_h_col_events_csv = false;

    // Step 16c: per-column h_col first-firing event log. Same I/O-proc
    // header-latch pattern as the clusters CSV above; written to
    // `<plot_file>_h_col_events.csv`. One row per column per run (the
    // first step the column's cec(i, j, k_top) flips 0 -> 1). Subsequent
    // firings of the same column are freeze-guarded (cec already 1, no
    // re-emit). Columns: time, col_i, col_j, a_f, h_col, is_gb_at_top,
    // Sp_field_at_top.
    bool sp_h_col_events_csv_header_written = false;

    // A1: every-removal event log `<plot_file>_removal_events.csv`
    // (spall.removal_events_csv; sp_weibull + spall.enabled, level 0). One
    // row per (step, column) whose surface moved: time, col_i, col_j, regime,
    // k_top, T_top, a_f, Sp_top, h_scan, h_applied, n_voided. Unlike the
    // edge-triggered h_col CSV above, it is a complete removal record.
    bool sp_emit_removal_events_csv = false;
    bool sp_removal_events_csv_header_written = false;

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
        pp.query_default("spallation.model", model_str, std::string("damage_law"));
        if      (model_str == "damage_law") value.spallation_model = SpallationModel::DamageLaw;
        else if (model_str == "sp_weibull") value.spallation_model = SpallationModel::SpWeibull;
        else Util::Abort(INFO,
            "spallation.model must be 'damage_law' or 'sp_weibull', got: ", model_str);

        // Step 18: confining pressure (Pa, default 0). Parsed regardless
        // of spallation_model so damage_law inputs that set the key
        // don't trip the strict-key checker; the offset is only applied
        // in the sp_weibull sigma_xx lambdas (UpdateSpAfterMechanics +
        // Step 16d's k_i_at_zcrack), so damage_law sees no effect.
        pp.query_default("confining.p", value.confining_p, 0.0);
        if (!std::isfinite(value.confining_p))
            Util::Abort(INFO, "confining.p must be finite, got ",
                        value.confining_p);

        if (value.spallation_model != SpallationModel::SpWeibull) return;

        // Hard mutual-exclusion: sp_weibull bypasses CZM (Step 15 Weibull
        // statistics carry the GB-bias work instead) and the eq.27 damage
        // evolution (Sp criterion is the onset law on this branch).
        if (value.cohesive_enabled)
            Util::Abort(INFO, "spallation.model = sp_weibull is incompatible with "
                              "cohesive.enabled = 1: CZM is bypassed on the sp_weibull branch");
        if (value.damage_enabled)
            Util::Abort(INFO, "spallation.model = sp_weibull is incompatible with "
                              "damage.enabled = 1: eq.27 damage evolution belongs to "
                              "the damage_law branch");

        pp.query_required("spallation.sp.a0", value.sp_a0);
        if (!(std::isfinite(value.sp_a0) && value.sp_a0 > 0.0))
            Util::Abort(INFO, "spallation.sp.a0 must be finite and > 0, got ", value.sp_a0);

        pp.query_default("spallation.sp.n_segments", value.sp_n_segments, 64);
        if (value.sp_n_segments < 2)
            Util::Abort(INFO, "spallation.sp.n_segments must be >= 2, got ",
                        value.sp_n_segments);
        if (value.sp_n_segments % 2 != 0) ++value.sp_n_segments;

        std::string law_str;
        pp.query_default("spallation.sp.kic.law", law_str, std::string("nasseri_table"));
        if      (law_str == "nasseri_table") value.kic_law = KIcLaw::NasseriTable;
        else if (law_str == "linear")        value.kic_law = KIcLaw::Linear;
        else Util::Abort(INFO,
            "spallation.sp.kic.law must be 'nasseri_table' or 'linear', got: ", law_str);

        if (value.kic_law == KIcLaw::NasseriTable)
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

        int presc_on = 0;
        pp.query_default("spallation.sp.prescribed_stress.enabled", presc_on, 0);
        value.sp_prescribed_sigma_enabled = (presc_on != 0);
        if (value.sp_prescribed_sigma_enabled)
        {
            std::string expr;
            pp.query_required("spallation.sp.prescribed_stress.xx_expr", expr);
            value.sp_prescribed_sigma_parser = amrex::Parser(expr);
            value.sp_prescribed_sigma_parser.registerVariables({"x", "y", "z", "t"});
            value.sp_prescribed_sigma_f = value.sp_prescribed_sigma_parser.compile<4>();
        }

        int presc_T_on = 0;
        pp.query_default("spallation.sp.prescribed_T_surface.enabled", presc_T_on, 0);
        value.sp_prescribed_T_surface_enabled = (presc_T_on != 0);
        if (value.sp_prescribed_T_surface_enabled)
        {
            std::string expr;
            pp.query_required("spallation.sp.prescribed_T_surface.expr", expr);
            value.sp_prescribed_T_surface_parser = amrex::Parser(expr);
            value.sp_prescribed_T_surface_parser.registerVariables({"t"});
            value.sp_prescribed_T_surface_f =
                value.sp_prescribed_T_surface_parser.compile<1>();
        }

        // Step 15b: Weibull per-cell flaw-length distribution. Default-off.
        // Only parsed under sp_weibull (we already returned for damage_law);
        // the microstructure requirement is enforced by the existing
        // sp_weibull + !microstructure_enabled abort. With
        // weibull.enabled = 0 the Sp evaluator uses the scalar
        // spallation.sp.a0 exactly as Step 12, which keeps the §12c
        // deterministic-a0 verification path intact.
        int weibull_on = 0;
        pp.query_default("weibull.enabled", weibull_on, 0);
        value.weibull_enabled = (weibull_on != 0);
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

            pp.query_default("weibull.m", value.weibull_m, 15.0);
            if (!(std::isfinite(value.weibull_m) && value.weibull_m > 0.0))
                Util::Abort(INFO, "weibull.m must be finite and > 0, got ",
                            value.weibull_m);

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

            // Step 15c: detach-mode gate + cluster-area threshold. Default
            // `per_face` preserves §15b behaviour byte-identical. Under
            // `connected_cluster`, the labeller zeroes Sp_cluster_id_mf on
            // any cluster with area strictly below weibull.A_crit (m^2);
            // Sp_field_mf is never modified by this gate.
            std::string detach_str;
            pp.query_default("weibull.detach_mode", detach_str,
                             std::string("per_face"));
            if      (detach_str == "per_face")
                value.weibull_detach_mode = DetachMode::PerFace;
            else if (detach_str == "connected_cluster")
                value.weibull_detach_mode = DetachMode::ConnectedCluster;
            else Util::Abort(INFO,
                "weibull.detach_mode must be 'per_face' or 'connected_cluster', "
                "got: ", detach_str);

            if (value.weibull_detach_mode == DetachMode::ConnectedCluster)
            {
                pp.query_required("weibull.A_crit", value.weibull_A_crit);
                if (!(std::isfinite(value.weibull_A_crit) &&
                      value.weibull_A_crit > 0.0))
                    Util::Abort(INFO,
                        "weibull.A_crit must be finite and > 0 under "
                        "detach_mode = connected_cluster, got ",
                        value.weibull_A_crit);
            }
        }
    }
};
} // namespace Integrator
#endif
