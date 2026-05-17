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

        int debug_amr_mechanics = 0;
        pp.query_default("mmw.debug_amr_mechanics", debug_amr_mechanics, 0);
        value.debug_amr_mechanics = (debug_amr_mechanics != 0);

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
                value.RegisterNewFab(value.removed_mf,         nullptr, 1, 0, "removed",         true);
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
        if (debug_amr_mechanics && finest_level > 0)
            CheckAMRMechanicsFieldsFinite("before mechanics solve", /*include_stress=*/false);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);
        if (debug_amr_mechanics && finest_level > 0)
            CheckAMRMechanicsFieldsFinite("after mechanics solve", /*include_stress=*/true);
        DiagnoseAMRResidual(a_time, a_step);
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
        if (debug_amr_mechanics && finest_level > 0)
            CheckAMRMechanicsFieldsFinite("after regrid repair", /*include_stress=*/true);
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
            ApplyAMRF0Consistency();
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
        ApplyAMRF0Consistency();
    }

    // AMR multi-AMR-level consistency for thermo-elastic eigenstrain F0.
    //
    // Per-level UpdateModel builds F0 from each level's local cell-centered
    // state via CellToNodeAverageInDomain. Even when the underlying cell field
    // is averaged-down across levels, the cell-to-node stencil produces
    // different nodal F0 values at coarse-fine overlap (the coarse node
    // averages 8 coarse cells while the fine node at the same physical
    // location averages 8 fine cells). The composite MLMG static-elastic
    // solve cannot reconcile the per-level rhs mismatch (rhs_lev = b_lev -
    // div(C_lev*F0_lev)) at the AMR overlap, leading to an irreducible
    // residual floor. See tests/MMWSpalling/hu_spall_onset/AMR_FAILURE_ANALYSIS.md.
    //
    // Modes:
    //   probe      (env ALAMO_LEV1_F0_FROM_COARSE): overwrite lev>=1's F0 with
    //              InterpFromCoarseLevel of lev-1's F0. Diagnostic only.
    //   production (default): average-down nodal F0 fine->coarse so each
    //              coarse level's F0 at AMR overlap matches the fine-level
    //              restricted view. Standard "state lives on the finest
    //              covering level" convention.
    void ApplyAMRF0Consistency()
    {
        if (finest_level <= 0) return;
        const int nF0 = AMREX_SPACEDIM * AMREX_SPACEDIM;
        const bool probe_mode =
            (std::getenv("ALAMO_LEV1_F0_FROM_COARSE") != nullptr);

        amrex::Vector<amrex::MultiFab> F0_mf(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            F0_mf[lev].define(model_mf[lev]->boxArray(),
                              model_mf[lev]->DistributionMap(),
                              nF0, model_mf[lev]->nGrow());
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<Set::Scalar> const& f0a = F0_mf[lev].array(mfi);
                amrex::Array4<Model::Solid::Affine::Isotropic const> const& m =
                    model_mf[lev]->const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    const Set::Matrix& F0 = m(i, j, k).F0;
                    for (int p = 0; p < AMREX_SPACEDIM; ++p)
                        for (int q = 0; q < AMREX_SPACEDIM; ++q)
                            f0a(i, j, k, p * AMREX_SPACEDIM + q) = F0(p, q);
                });
            }
            F0_mf[lev].FillBoundary(geom[lev].periodicity());
        }

        if (probe_mode)
        {
            amrex::PhysBCFunctNoOp cphysbc, fphysbc;
            amrex::Vector<amrex::BCRec> bcs(nF0, amrex::BCRec());
            for (int lev = 1; lev <= finest_level; ++lev)
            {
                amrex::InterpFromCoarseLevel(
                    F0_mf[lev], 0.0, F0_mf[lev - 1],
                    0, 0, nF0, geom[lev - 1], geom[lev],
                    cphysbc, 0, fphysbc, 0,
                    refRatio(lev - 1),
                    &amrex::node_bilinear_interp, bcs, 0);
            }
        }
        else
        {
            for (int lev = finest_level; lev >= 1; --lev)
            {
                amrex::average_down_nodal(F0_mf[lev], F0_mf[lev - 1],
                                          refRatio(lev - 1));
            }
        }

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<Model::Solid::Affine::Isotropic> const& m =
                    model_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar const> const& f0a =
                    F0_mf[lev].const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix F0;
                    for (int p = 0; p < AMREX_SPACEDIM; ++p)
                        for (int q = 0; q < AMREX_SPACEDIM; ++q)
                            F0(p, q) = f0a(i, j, k, p * AMREX_SPACEDIM + q);
                    m(i, j, k).F0 = F0;
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
                                            bR, alpha_f, x, y, z, Tv);
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
        auto alpha_f = beam.alpha_f;
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];
        const Set::Scalar inv_dx2 = 1.0 / (DX[0] * DX[0]);
        const Set::Scalar inv_dy2 = 1.0 / (DX[1] * DX[1]);
        const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);

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
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tc = T_old(i, j, k);
                Set::Scalar Txm = T_old(i-1, j, k), Txp = T_old(i+1, j, k);
                Set::Scalar Tym = T_old(i, j-1, k), Typ = T_old(i, j+1, k);
                Set::Scalar Tzm = T_old(i, j, k-1), Tzp = T_old(i, j, k+1);
                Set::Scalar kc  = mat->kappa(Tc);
                // Harmonic-average κ on faces (handles discontinuous κ across
                // interfaces correctly; degenerates to arithmetic mean for
                // smooth κ(T)).
                auto kface = [](Set::Scalar a, Set::Scalar b) {
                    return (a + b > 0.0) ? 2.0 * a * b / (a + b) : 0.0;
                };
                Set::Scalar fxp = kface(kc, mat->kappa(Txp)) * (Txp - Tc) * inv_dx2;
                Set::Scalar fxm = kface(kc, mat->kappa(Txm)) * (Tc - Txm) * inv_dx2;
                Set::Scalar fyp = kface(kc, mat->kappa(Typ)) * (Typ - Tc) * inv_dy2;
                Set::Scalar fym = kface(kc, mat->kappa(Tym)) * (Tc - Tym) * inv_dy2;
                Set::Scalar fzp = kface(kc, mat->kappa(Tzp)) * (Tzp - Tc) * inv_dz2;
                Set::Scalar fzm = kface(kc, mat->kappa(Tzm)) * (Tc - Tzm) * inv_dz2;
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
                                           bR, alpha_f, x, y, z, Tc);

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
                                           bR, alpha_f, x, y, z, Tc);
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

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar P_now = beam.power_at(time);
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
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
        const Set::Scalar sp_hfl   = surface_patch.h_conv;
        const Set::Scalar sp_Tfl   = surface_patch.T_flame;

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

                if (!self_void && sp_on && k == sp_topk) {
                    Set::Scalar dxp = x - sp_x0, dyp = y - sp_y0;
                    if (dxp * dxp + dyp * dyp <= sp_r2) {
                        if (sp_presc)
                            div_kgrad += kc * (sp_T - Tc) * inv_dz2;
                        else
                            div_kgrad += sp_hfl * (sp_Tfl - Tc) * inv_dz;
                    }
                }

                // Beer-Lambert path is measured through material from the
                // current (receding) surface: phi(i,j,k) = z_surface - z.
                // Before any removal in this column phi == z0 - z, so this is
                // a no-op until the surface actually recedes.
                Set::Scalar Q = self_void ? Set::Scalar(0.0)
                                          : BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                                       bR, alpha_f, x, y, z, Tc,
                                                       void_aware
                                                         ? std::max(Set::Scalar(0.0), phi(i, j, k))
                                                         : Set::Scalar(-1.0));

                Set::Scalar eps = (Tc <= T_eps_b) ? eps_lo : eps_hi;
                Set::Scalar q_loss = self_void ? Set::Scalar(0.0)
                                               : surf(i, j, k)
                                                 * (eps * sig_SB * (Tc*Tc*Tc*Tc - T_a*T_a*T_a*T_a)
                                                    + h_c * (Tc - T_a));

                // Conservative H update; T is recovered in InvertHtoT
                // using the per-phase Numeric::Material::Table. Detached
                // (void) cells exchange no heat with neighbours, so their
                // H stays at the ambient value written by the spall reset.
                H(i, j, k) = H_old(i, j, k) + dt * (div_kgrad + Q - q_loss * inv_dz);
                keff_out(i, j, k) = kc; // diagnostic, refreshed each step
            });
        }
    }

    // Cell-centred Beer-Lambert + Gaussian volumetric source (eqs. 1-3 of
    // main.tex), scaled by (1 − R) for surface reflectivity. Returns 0
    // above the surface (z > z_0).
    //
    // `dz_abs` is the optical path the beam travels through *material* to
    // reach this cell. When removal is active and the surface has receded,
    // the beam crosses the void cavity unattenuated, so the Beer-Lambert
    // attenuation must be measured from the current surface — pass
    // dz_abs = z_surface(x,y) - z there (equivalently phi at this cell).
    // dz_abs < 0 (the default) falls back to the geometric depth z_0 - z,
    // which is correct whenever no material above this cell has been removed
    // (then z_surface == z_0 and phi == z_0 - z exactly). The Gaussian waist
    // still evolves with the geometric distance z_0 - z regardless, since the
    // beam propagates geometrically through the void cavity.
    AMREX_GPU_HOST_DEVICE
    static Set::Scalar BeamSource(Set::Scalar P, Set::Scalar w0, Set::Scalar lam,
                                  Set::Scalar x0, Set::Scalar y0, Set::Scalar z0,
                                  Set::Scalar R,
                                  amrex::ParserExecutor<1> alpha_f,
                                  Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                  Set::Scalar T,
                                  Set::Scalar dz_abs = -1.0)
    {
        Set::Scalar dz = z0 - z;
        if (dz < 0.0 || P == 0.0) return 0.0;
        Set::Scalar zR = M_PI * w0 * w0 / lam;
        Set::Scalar w  = w0 * std::sqrt(1.0 + (dz / zR) * (dz / zR));
        Set::Scalar r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        Set::Scalar Pi = (2.0 * P / (M_PI * w0 * w0)) * (w0 / w) * (w0 / w)
                       * std::exp(-2.0 * r2 / (w * w));
        Set::Scalar a    = alpha_f(T);
        Set::Scalar path = (dz_abs >= 0.0) ? dz_abs : dz;
        return (1.0 - R) * Pi * a * std::exp(-a * path);
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
        if (debug_amr_mechanics && !damage_prescribed_stress_enabled)
            CheckMatrixFieldFinite(stress_mf, "stress_mf", lev,
                                   "before damage update");

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
            amrex::Box valid_node = mfi.validbox();
            valid_node.convert(amrex::IntVect::TheNodeVector());
            amrex::Box domain_node = geom[lev].Domain();
            domain_node.convert(amrex::IntVect::TheNodeVector());

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
                    if (debug_amr_mechanics)
                    {
                        for (int dk = 0; dk <= 1; ++dk)
                            for (int dj = 0; dj <= 1; ++dj)
                                for (int di = 0; di <= 1; ++di)
                                {
                                    int ii = i + di;
                                    int jj = j + dj;
                                    int kk = k + dk;
                                    if (!FiniteMatrix(sig_node(ii, jj, kk)))
                                    {
                                        amrex::IntVect iv(AMREX_D_DECL(ii, jj, kk));
                                        AbortBadAMRValue("damage stress sampling",
                                                         "stress_mf", lev, -1,
                                                         iv, valid_node, domain_node);
                                    }
                                }
                    }
                    sigma = Numeric::Interpolate::NodeToCellAverage(sig_node, i, j, k, 0);
                }

                if (!FiniteMatrix(sigma))
                {
                    if (!prescribed)
                    {
                        for (int dk = 0; dk <= 1; ++dk)
                            for (int dj = 0; dj <= 1; ++dj)
                                for (int di = 0; di <= 1; ++di)
                                {
                                    int ii = i + di;
                                    int jj = j + dj;
                                    int kk = k + dk;
                                    if (!FiniteMatrix(sig_node(ii, jj, kk)))
                                    {
                                        amrex::IntVect iv(AMREX_D_DECL(ii, jj, kk));
                                        AbortBadAMRValue("damage stress sampling",
                                                         "stress_mf", lev, -1,
                                                         iv, valid_node, domain_node);
                                    }
                                }
                    }
                    Util::Abort(INFO, "non-finite averaged stress before Drucker-Prager at cell (",
                                i, ",", j, ",", k, ")");
                }

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

    // Step 10: principal-stress ratio for the production spall path. ALAMO
    // stores tensile stress as positive, so the most compressive principal
    // stress is the smallest (most negative) eigenvalue. We recast to
    // compressive-positive (sigma_1 = -e_min, sigma_3 = -e_max). When all
    // three principals are compressive (sigma_3 > 0), the ratio
    // sigma_1 / sigma_3 is the canonical Hoek-Brown anisotropy measure.
    //
    // Step 14b refinement: when any principal stress is tensile or the
    // state is degenerate (sigma_3 <= 0, or a non-finite ratio results),
    // the function returns `TENSILE_RATIO_SENTINEL` instead of 0. Mixed/
    // tensile states are unconditionally above the Hoek-Brown threshold —
    // a free surface with tensile sigma_zz (e.g. Hu's prescribed-T patch
    // generating sigma_zz ~ +55 MPa at the top cell) is *more* susceptible
    // to mode-I spall, not less, so the spall gate should pass. The
    // sentinel (1e30) is unconditionally above any plausible `xi_s` and
    // preserves the existing `if (ratio >= xi_s)` call-site semantics
    // without a signature change.
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

    // Step 10/11: per-column unified surface advancement. Computes a
    // candidate spallation thickness (when spall is enabled) and a
    // candidate vaporisation thickness (when vapor is enabled) for the
    // top material cell of each column, advances phi by the larger of
    // the two, and tags `regime` with the winning mode (1 = spall,
    // 2 = vapor). Detached cells reset H/T/D/phase fractions exactly as
    // in Step 10.
    //
    // Parallel correctness: domain-wide per-column buffers indexed by
    // (j*nx_dom + i) carry the global topmost material cell, the
    // per-column candidate thicknesses, and the chosen regime. A first
    // Step 12: Sp = K_I/K_Ic(T) surface diagnostic. Runs only under
    // `spallation.model = sp_weibull`. Writes Sp into the topmost cell of
    // each (i, j) column and 0 elsewhere; the rest of the integrator's
    // physics is unchanged (DP/CZM evolution stays off on this branch).
    //
    // sigma_xx source per column:
    //   - sp_prescribed_sigma_enabled: parser sigma_xx(x, y, z, t),
    //     evaluated at integration depths z = z_top_face - xi * a0. Used by
    //     the Step 12c verification.
    //   - else: nodal-stress NodeToCellAverage at the top cell, treated as
    //     constant over [0, a0]. Valid while a0 << dz; revisit when Step 15
    //     introduces a-distributions with a0 comparable to dz.
    //
    // T_surface source for K_Ic:
    //   - sp_prescribed_T_surface_enabled: parser T(t).
    //   - else: temp_mf at the top cell.
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
        const Set::Scalar a0_loc      = sp_a0;
        const int         n_seg       = sp_n_segments;
        const KIcLaw      law         = kic_law;
        const Set::Scalar kic_scale_  = kic_scale;
        const Set::Scalar K_Ic0_      = kic_K_Ic0;
        const Set::Scalar T_ref_kic_  = kic_T_ref;
        const Set::Scalar T_melt_     = kic_T_melt;

        // Step 18: confining-pressure contribution to compression-positive
        // sigma_xx. Per Kant Eq. 15, p enters as +p·nu/(1-nu) under 1-D-
        // confined geometry. Phase-0 nu = E/(2·mu) - 1 (single-phase
        // scope; multi-phase is out of scope per the active step). Exact-
        // zero short-circuit guarantees byte-identity when confining.p
        // is unset (the default).
        Set::Scalar sp_conf_offset = 0.0;
        if (confining_p != 0.0 && !phases.empty() &&
            phases[0].E > 0.0 && phases[0].mu > 0.0)
        {
            const Set::Scalar nu_local =
                phases[0].E / (2.0 * phases[0].mu) - 1.0;
            if (nu_local > -1.0 && nu_local < 1.0)
                sp_conf_offset = confining_p * nu_local / (1.0 - nu_local);
        }

        if (!use_presc_s)
            Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

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
            if (!use_presc_s) sig_node = (*stress_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> flaw_a;
            if (use_weibull) flaw_a = (*flaw_a_mf[lev]).array(mfi);

            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (k != top_k) return;

                Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];

                Set::Scalar Tsurf = use_presc_T ? presc_T_f(time) : T_a(i, j, k);

                // Step 15b: the K_I integrand consumes -sigma_xx
                // (compression-positive). ALAMO's stress_mf and the
                // prescribed expression are both tension-positive sigma_xx,
                // so a compressive heated surface (sigma_xx < 0) gives
                // K_I > 0. The physical path samples stress_mf at the cell
                // containing the requested depth (k_d), clamped into the
                // owned box (the integration window is assumed to lie
                // within the top box's z-extent for this packet;
                // z-decomposition-robust column gather is deferred).
                auto sigma_at_depth = [&] (Set::Scalar depth) -> Set::Scalar {
                    // Step 18: +p·nu/(1-nu) confining-pressure contribution
                    // added to compression-positive sigma_xx. Default-zero
                    // short-circuit at the offset capture above keeps every
                    // pre-Step-18 input byte-identical.
                    if (use_presc_s)
                    {
                        Set::Scalar z_phys = z_top_face - depth;
                        return -presc_s_f(xc, yc, z_phys, time) + sp_conf_offset;
                    }
                    int k_d = top_k - static_cast<int>(std::floor(depth / dz));
                    if (k_d > top_k)    k_d = top_k;
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

        // Domain-wide firing buffer: 1 on firing top-cells, 0 elsewhere.
        // ReduceIntMax across ranks (the rank that owns the (i, j, top_k)
        // cell writes the value, others stay at 0 — max collapses cleanly).
        amrex::Vector<int> fire(ncols, 0);
        for (amrex::MFIter mfi(*Sp_field_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            if (top_k < klo || top_k > khi) continue;
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
            amrex::Array4<const Set::Scalar> const& Sp_a =
                (*Sp_field_mf[lev]).array(mfi);
            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                if (Sp_a(i, j, top_k) >= 1.0)
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

        // Write Sp_cluster_id_mf everywhere on this level: top-z plane gets
        // the cluster ID (or 0 for filtered-out / non-firing cells); 0
        // elsewhere.
        Sp_cluster_id_mf[lev]->setVal(0.0);
        for (amrex::MFIter mfi(*Sp_cluster_id_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& cl =
                (*Sp_cluster_id_mf[lev]).array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (k != top_k) return;
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
        if (lev == 0 && cluster_count > 0 &&
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

    // pass populates each rank's local maximum k_top via ReduceIntMax;
    // a second pass lets the rank that owns (i, j, k_top_global)
    // compute the spall and vapor candidates and write into the
    // per-column buffers, which are summed across ranks. The final
    // pass applies the column-wide phi shift and per-cell detach
    // bookkeeping locally on every rank using the reduced data. This
    // makes the update correct for arbitrary AMReX z-decompositions,
    // not just column-coherent partitions.
    void UpdateRemovalAfterCohesive(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (!removal_enabled()) return;
        Util::Assert(INFO, TEST(microstructure_enabled));

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz   = DX[2];

        const amrex::Box& domain = geom[lev].Domain();
        const int dom_ilo = domain.smallEnd(0), dom_ihi = domain.bigEnd(0);
        const int dom_jlo = domain.smallEnd(1), dom_jhi = domain.bigEnd(1);
        const int nx_dom = dom_ihi - dom_ilo + 1;
        const int ny_dom = dom_jhi - dom_jlo + 1;
        const int ncols  = nx_dom * ny_dom;

        // Spall calibration ---------------------------------------------
        const bool spall_on = spall_enabled;
        const Set::Scalar D_th = spall_damage_threshold;
        const Set::Scalar xi_s = spall_xi_s;
        Set::Scalar h_spall = 0.0;
        if (spall_on)
        {
            const Set::Scalar C_h  = spall_C_h;
            const Set::Scalar alpha_eff   = spall_alpha_override ? spall_alpha   : therm_alpha;
            const Set::Scalar t_spall_eff = spall_t_spall_override ? spall_t_spall : (time + dt);
            if (!(std::isfinite(alpha_eff) && alpha_eff > 0.0))
                Util::Abort(INFO, "spall: thermal diffusivity for h_spall is non-positive (",
                            alpha_eff, "); set spall.alpha or hc.heat.alpha");
            if (!(std::isfinite(t_spall_eff) && t_spall_eff > 0.0))
                Util::Abort(INFO, "spall: t_spall for h_spall is non-positive (", t_spall_eff,
                            "); set spall.t_spall");
            h_spall = C_h * std::sqrt(alpha_eff * t_spall_eff);
            if (!(std::isfinite(h_spall) && h_spall >= 0.0))
                Util::Abort(INFO, "spall: computed h_spall is non-finite (", h_spall, ")");
        }
        const bool spall_prescribed = spall_prescribed_ratio_enabled;
        auto ratio_f = spall_ratio_f;
        const bool have_stress = spall_on && !spall_prescribed;
        if (have_stress) Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

        // Vapor calibration ---------------------------------------------
        const bool vapor_on = vapor_enabled;
        const Set::Scalar A_beam = vapor_A_beam;
        const Set::Scalar T0_v   = vapor_T0;
        const Set::Scalar P_now  = vapor_on ? beam.power_at(time) : 0.0;
        const Set::Scalar oneR   = 1.0 - beam.R;
        // Step 17: local Gaussian-beam parameters for per-cell v_vap
        // under sp_weibull. Identical to the values used by the heat
        // source path; `beam.power_at(time)` is already captured in
        // `P_now`. At the surface (z = z0) the beam waist is w0 with
        // no spreading; per-area surface power density is then
        // (2·P / (π·w0²)) · exp(-2·r² / w0²). Captured outside the
        // MFIter so the inner block can read them.
        const Set::Scalar bw0_v  = beam.omega0;
        const Set::Scalar bx0_v  = beam.x0;
        const Set::Scalar by0_v  = beam.y0;
        const bool vapor_prescribed = vapor_prescribed_rop_enabled;
        auto vapor_rop_eval = vapor_rop_f;
        if (vapor_on && !(std::isfinite(A_beam) && A_beam > 0.0))
            Util::Abort(INFO, "vapor: A_beam is non-positive (", A_beam,
                        "); set vapor.A_beam");

        // Per-step diagnostic fields are reset here. `phi` and `removed`
        // carry state across steps and are not reset.
        regime_mf[lev]->setVal(0.0);
        RoP_mf[lev]->setVal(0.0);
        if (spall_on)
        {
            spall_event_mf[lev]->setVal(0.0);
            spall_thickness_mf[lev]->setVal(0.0);
            RoP_spall_mf[lev]->setVal(0.0);
        }
        if (vapor_on)
        {
            vapor_event_mf[lev]->setVal(0.0);
            vapor_thickness_mf[lev]->setVal(0.0);
            RoP_vap_mf[lev]->setVal(0.0);
        }
        // Step 17: per-step reset for the per-cell regime field under
        // sp_weibull. The field is registered only when
        // spallation_model == SpWeibull && (spall_on || vapor_on).
        if (spallation_model == SpallationModel::SpWeibull &&
            (spall_on || vapor_on))
        {
            regime_field_mf[lev]->setVal(0.0);
        }
        // Step 16: sp_weibull spall path flag. The h_spall_field_mf
        // diagnostic accumulates as a running max per column (written
        // at the INITIAL top-z slice z = top_k_dom, which doesn't move
        // as the surface erodes) — see PASS 2 write. The field is NOT
        // reset per step, so the last plotfile records the deepest h_col
        // every column has ever produced.
        const bool sp_weibull_spall =
            spall_on && spallation_model == SpallationModel::SpWeibull;

        const Numeric::Material::Table *ptabs =
            phase_tables.empty() ? nullptr : phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        const Set::Scalar T_amb_local = T_amb;
        const Set::Scalar dt_local = dt;

        auto col_index = [=](int i, int j) {
            return (j - dom_jlo) * nx_dom + (i - dom_ilo);
        };

        // PASS 1: each rank fills its local maximum k_top per column,
        // sentinel INT_MIN for columns it does not own at all. ReduceIntMax
        // collapses to the global topmost material cell per (i, j). After
        // the reduce, every rank knows the same k_top_global.
        amrex::Vector<int> k_top_global(ncols, std::numeric_limits<int>::min());
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
                int k_top_local = std::numeric_limits<int>::min();
                for (int k = khi; k >= klo; --k)
                {
                    if (rem(i, j, k) <= 0.5 && phi(i, j, k) >= 0.0) { k_top_local = k; break; }
                }
                if (k_top_local > std::numeric_limits<int>::min())
                {
                    int idx = col_index(i, j);
                    if (k_top_local > k_top_global[idx]) k_top_global[idx] = k_top_local;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(k_top_global.data(), ncols);

        // PASS 1.5 (Step 14b): per-column max D over the configured scan
        // window. The window depends on `spall.sample`:
        //   top_cell    -> single cell at k_top (preserves Step 10/11
        //                  behaviour byte-identical when this mode is
        //                  selected by default).
        //   column_max  -> [box.smallEnd(2), k_top] within each rank's
        //                  local box; ReduceRealMax collapses across
        //                  the MPI z-decomposition.
        //   shallow_max -> [k_top - floor(shallow_depth/dz), k_top]
        //                  within each rank's local box; same reduction.
        // Buffers initialise to 0.0 so columns with no contributor land
        // at 0 (the same effective "no fire" result as the old "skip"
        // path).
        amrex::Vector<Set::Scalar> D_max_col(ncols, 0.0);
        if (spall_on)
        {
            const SpallSampleMode sample_mode = spall_sample_mode;
            const Set::Scalar shallow_d       = spall_shallow_depth;
            const Set::Scalar dz_local        = DX[2];
            const int k_band_shallow =
                static_cast<int>(std::floor(shallow_d / dz_local));

            for (amrex::MFIter mfi(*D_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& Dfab =
                    (*D_mf[lev]).array(mfi);
                const int klo_box = bx.smallEnd(2), khi_box = bx.bigEnd(2);
                const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
                const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);

                for (int j = jlo; j <= jhi; ++j)
                for (int i = ilo; i <= ihi; ++i)
                {
                    int idx = col_index(i, j);
                    int k_top = k_top_global[idx];
                    if (k_top == std::numeric_limits<int>::min()) continue;

                    int k_hi_scan, k_lo_scan;
                    switch (sample_mode)
                    {
                        case SpallSampleMode::TopCell:
                            // Only the rank that owns k_top contributes;
                            // every other rank skips so the reduction
                            // picks up exactly D(i, j, k_top).
                            if (k_top < klo_box || k_top > khi_box) continue;
                            k_hi_scan = k_top;
                            k_lo_scan = k_top;
                            break;
                        case SpallSampleMode::ColumnMax:
                            k_hi_scan = std::min(k_top, khi_box);
                            if (k_hi_scan < klo_box) continue;
                            k_lo_scan = klo_box;
                            break;
                        case SpallSampleMode::ShallowMax:
                            k_hi_scan = std::min(k_top, khi_box);
                            if (k_hi_scan < klo_box) continue;
                            k_lo_scan = std::max(klo_box,
                                                 k_top - k_band_shallow);
                            if (k_lo_scan > k_hi_scan) continue;
                            break;
                        default:
                            continue;  // unreachable
                    }

                    Set::Scalar local_max = 0.0;
                    for (int k = k_lo_scan; k <= k_hi_scan; ++k)
                    {
                        Set::Scalar v = Dfab(i, j, k);
                        if (!std::isfinite(v)) continue;
                        if (v > local_max) local_max = v;
                    }
                    if (local_max > D_max_col[idx]) D_max_col[idx] = local_max;
                }
            }
            amrex::ParallelDescriptor::ReduceRealMax(D_max_col.data(), ncols);
        }

        // PASS 2: only the rank that owns (i, j, k_top_global) writes
        // into the per-column candidate buffers; others contribute 0.
        // ReduceRealSum / ReduceIntSum collapse the per-rank
        // contributions onto every rank.
        amrex::Vector<int> spall_fires_col(ncols, 0);
        amrex::Vector<int> vapor_fires_col(ncols, 0);
        amrex::Vector<Set::Scalar> h_s_cand_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> h_v_cand_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> rop_s_cand_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> rop_v_cand_col(ncols, 0.0);

        // Step 16: capture sp_weibull-specific K_Ic + scan parameters
        // outside the MFIter so the lambda in the depth scan can close over
        // them. These values are unused when spallation_model == DamageLaw.
        const KIcLaw      law_loc          = kic_law;
        const Set::Scalar kic_scale_loc    = kic_scale;
        const Set::Scalar kic_K_Ic0_loc    = kic_K_Ic0;
        const Set::Scalar kic_T_ref_loc    = kic_T_ref;
        const Set::Scalar kic_T_melt_loc   = kic_T_melt;
        const Set::Scalar h_spall_max_loc  = spall_h_spall_max;
        const int         h_spall_n_seg_loc =
            (spall_h_spall_n_segments > 0) ? spall_h_spall_n_segments
                                           : sp_n_segments;

        // Step 16d: capture phase 0 material params + the prescribed-T
        // face callable for the T-derived sigma_xx reconstruction
        // (option (a) in rossi-validation-diagnostic-design.md §8).
        // The bisection sample uses
        //   sigma_xx(z) = -E * beta * (T(z) - T_ref) / (1 - nu)
        // with T(z) linearly interpolated between sp_T_face_loc(time) at
        // z = 0 and cell-centred temp_mf values at z = (n + 0.5) * dz.
        // Single-phase Voronoi only — multi-phase requires per-cell phase
        // lookup which is out of scope per the active step's Guardrails.
        // Fail fast if surface_patch isn't in prescribed_T mode (the face
        // callable would be invalid).
        Set::Scalar phase_E_loc     = 0.0;
        Set::Scalar phase_beta_loc  = 0.0;
        Set::Scalar phase_nu_loc    = 0.0;
        Set::Scalar phase_T_ref_loc = 0.0;
        auto sp_T_face_loc          = surface_patch.T_f;  // ParserExecutor<1>
        if (sp_weibull_spall)
        {
            if (phases.empty())
                Util::Abort(INFO, "Step 16d (sp_weibull + spall.enabled): "
                                  "phases vector is empty; need at least "
                                  "phase 0 for T-derived sigma_xx");
            phase_E_loc    = phases[0].E;
            phase_beta_loc = phases[0].beta;
            if (!(phases[0].E > 0.0 && phases[0].mu > 0.0))
                Util::Abort(INFO, "Step 16d: phase 0 must have E > 0 and "
                                  "mu > 0 for the T-derived sigma_xx path; "
                                  "got E=", phases[0].E, ", mu=",
                                  phases[0].mu);
            phase_nu_loc    = phases[0].E / (2.0 * phases[0].mu) - 1.0;
            phase_T_ref_loc = phases[0].T_ref;
            if (!(phase_nu_loc > -1.0 && phase_nu_loc < 0.5))
                Util::Abort(INFO, "Step 16d: phase 0 Poisson's ratio out "
                                  "of physical range; got nu=", phase_nu_loc);
            if (surface_patch.mode !=
                SurfacePatch::Mode::PrescribedTemperature)
                Util::Abort(INFO, "Step 16d (sp_weibull + spall.enabled): "
                                  "surface_patch.mode = prescribed_T is "
                                  "required for the T-derived sigma_xx path");
        }

        // Step 18: confining-pressure offset to compression-positive
        // sigma_xx inside the Step 16d k_i_at_zcrack lambda. Same
        // formula as the UpdateSpAfterMechanics offset (Kant Eq. 15
        // +p·nu/(1-nu)). Exact-zero short-circuit preserves byte-
        // identity when confining.p is unset.
        Set::Scalar sp_conf_offset = 0.0;
        if (sp_weibull_spall && confining_p != 0.0 &&
            phase_nu_loc > -1.0 && phase_nu_loc < 1.0)
        {
            sp_conf_offset = confining_p * phase_nu_loc / (1.0 - phase_nu_loc);
        }

        // Step 16c: per-rank buffer of first-firing CSV events. A column's
        // event is emitted ONCE per run (on the step its cec(i, j, k_top)
        // flips 0 -> 1). After PASS 2's MFIter loop ends, each rank in
        // turn appends its rows to <plot_file>_h_col_events.csv with the
        // I/O proc writing the header on first call. The buffer holds
        // (time, col_i, col_j, a_f, h_col, is_gb_at_top, Sp_at_top).
        struct HColEvent {
            Set::Scalar t, a_f, h_col, sp;
            int col_i, col_j, is_gb;
        };
        std::vector<HColEvent> local_h_col_events;

        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> sig_node;
            if (have_stress) sig_node = (*stress_mf[lev]).array(mfi);
            // Step 16: per-cell h_spall diagnostic write target + the
            // Sp_cluster_id detach gate / flaw_a / temperature reads for
            // the depth scan. All four fields exist iff sp_weibull_spall
            // is true (Register block in Parse()).
            // Step 16b: crack_event_count_mf write target (edge-triggered
            // 3-D first-crossing count field, design note §2).
            amrex::Array4<Set::Scalar>      hsf;
            amrex::Array4<Set::Scalar>      cec;
            amrex::Array4<const Set::Scalar> clu, fa_a, T_a;
            // Step 16c: bind Sp_field + is_grain_boundary on the top
            // slice to populate the per-column h_col CSV event row.
            amrex::Array4<const Set::Scalar> Sp_top, is_gb_top;
            if (sp_weibull_spall)
            {
                hsf       = (*h_spall_field_mf[lev]).array(mfi);
                cec       = (*crack_event_count_mf[lev]).array(mfi);
                clu       = (*Sp_cluster_id_mf[lev]).array(mfi);
                fa_a      = (*flaw_a_mf[lev]).array(mfi);
                T_a       = (*temp_mf[lev]).array(mfi);
                Sp_top    = (*Sp_field_mf[lev]).array(mfi);
                is_gb_top = (*is_grain_boundary_mf[lev]).array(mfi);
            }

            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);

            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                int idx = col_index(i, j);
                int k_top = k_top_global[idx];
                if (k_top == std::numeric_limits<int>::min()) continue;
                if (k_top < klo || k_top > khi) continue; // not owned here

                int p_top = static_cast<int>(std::round(ph(i, j, k_top)));
                if (p_top < 0 || p_top >= n_phases) p_top = 0;

                Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar zc = PLO[2] + (Set::Scalar(k_top) + 0.5) * DX[2];

                if (spall_on)
                {
                    if (sp_weibull_spall)
                    {
                        // Step 16: Sp-vs-depth zero-crossing h_spall.
                        // Detach gate is Sp_cluster_id_mf > 0 (consumes
                        // the §15c labeller — propagates per_face vs
                        // connected_cluster + A_crit cleanly to material
                        // removal). h_spall is the deepest crack-tip depth
                        // z_crack at which K_I(z_crack; a_f) >= K_Ic(T(z))
                        // for the hypothetical edge crack of fixed length
                        // a_f starting at z_crack below the surface.
                        Set::Scalar cluster_id = clu(i, j, k_top);
                        if (cluster_id > 0.0)
                        {
                            Set::Scalar a_f = fa_a(i, j, k_top);
                            // Step 16c: capture cec(i, j, k_top) BEFORE the
                            // scan modifies it, so the post-scan check can
                            // detect a 0 -> 1 transition (i.e. this column's
                            // first firing this run).
                            const Set::Scalar cec_at_top_pre =
                                cec(i, j, k_top);
                            if (std::isfinite(a_f) && a_f > 0.0 && have_stress)
                            {
                                // Walk z_crack downward from the surface
                                // in 1-cell strides. Stop at the first
                                // K_I < K_Ic crossing OR h_spall_max OR
                                // the bottom of the local box. The
                                // sigma_at_depth callable samples
                                // stress_mf at the absolute depth (z_crack
                                // + xi*a_f) below the surface; the
                                // compression-positive (-σ_xx) convention
                                // lives here, matching
                                // UpdateSpAfterMechanics.
                                //
                                // Step 16b: also drives the edge-triggered
                                // crack_event_count_mf update (a passing
                                // K_I >= K_Ic at any (i, j, k_tip) cell
                                // flips that cell's count from 0 to 1,
                                // never re-tested afterwards — option (ii)
                                // freeze in rossi-validation-diagnostic-
                                // design.md §3) and the bisection that
                                // sharpens h_col below cell resolution.
                                Set::Scalar h_col = 0.0;
                                const int max_steps_cap = static_cast<int>(
                                    std::floor(h_spall_max_loc / dz));
                                const int max_steps_box = k_top - klo;
                                const int max_steps =
                                    std::min(max_steps_cap, max_steps_box);
                                const int klo_loc = klo;
                                const int khi_loc = khi;

                                // Local helper: K_Ic(T) under the
                                // configured law, in Pa·m^0.5.
                                auto kic_pa_at_cell =
                                    [&] (int k_cell) -> Set::Scalar {
                                    Set::Scalar T_at = T_a(i, j, k_cell);
                                    Set::Scalar K_Ic_MPa;
                                    if (law_loc == KIcLaw::NasseriTable)
                                        K_Ic_MPa =
                                            Numeric::SpCriterion::KIcNasseri(
                                                T_at, kic_scale_loc);
                                    else
                                        K_Ic_MPa =
                                            Numeric::SpCriterion::KIcLinear(
                                                T_at, kic_K_Ic0_loc,
                                                kic_T_ref_loc, kic_T_melt_loc);
                                    return Numeric::SpCriterion
                                        ::MPa_sqrtm_to_Pa_sqrtm(K_Ic_MPa);
                                };

                                // Local helper: K_I(z_crack; a_f) using
                                // the same Tada integrator as
                                // UpdateSpAfterMechanics, with the
                                // compression-positive (-sigma_xx)
                                // convention in the caller.
                                //
                                // Step 16d: T-derived sigma_xx with the
                                // prescribed-T face BC (option (a) per
                                // rossi-validation-diagnostic-design.md
                                // §8). T(z) is linearly interpolated
                                // between bracketing T values:
                                //   z = 0:           sp_T_face_loc(time)
                                //   z = (n+0.5)*dz:  T_a(i, j, k_top - n)
                                // sigma_xx(z) from the thermoelastic
                                // relation under 1-D confinement:
                                //   sigma_xx = -E * beta * (T - T_ref) / (1 - nu)
                                // Single-phase Voronoi only — uses
                                // phases[0] material params captured
                                // outside the MFIter loop. The surface
                                // BC recovery (z ∈ [0, dz/2]) is the
                                // root-cause fix for the Step 16c
                                // option (L) flat-extrapolation artefact
                                // documented in design note §7.
                                auto k_i_at_zcrack =
                                    [&] (Set::Scalar z_crack) -> Set::Scalar {
                                    auto sigma_at_depth =
                                        [&] (Set::Scalar d) -> Set::Scalar {
                                        const Set::Scalar z_abs =
                                            z_crack + d;
                                        const Set::Scalar half_dz =
                                            0.5 * dz;
                                        Set::Scalar T_z;
                                        if (z_abs <= half_dz)
                                        {
                                            // Surface region: lerp
                                            // between T_face (prescribed)
                                            // and top-cell-centre T.
                                            const Set::Scalar T_face =
                                                sp_T_face_loc(time);
                                            const Set::Scalar T_top =
                                                T_a(i, j, k_top);
                                            Set::Scalar alpha =
                                                z_abs / half_dz;
                                            if (alpha < 0.0) alpha = 0.0;
                                            if (alpha > 1.0) alpha = 1.0;
                                            T_z = (1.0 - alpha) * T_face
                                                +        alpha  * T_top;
                                        }
                                        else
                                        {
                                            // Interior: linear T interp
                                            // between adjacent cell
                                            // centres at z = (k + 0.5)*dz.
                                            const Set::Scalar u =
                                                z_abs / dz - 0.5;
                                            const Set::Scalar uf =
                                                std::floor(u);
                                            int n = static_cast<int>(uf);
                                            int k_above = k_top - n;
                                            int k_below = k_above - 1;
                                            if (k_above > k_top)   k_above = k_top;
                                            if (k_above < klo_loc) k_above = klo_loc;
                                            if (k_above > khi_loc) k_above = khi_loc;
                                            if (k_below > k_top)   k_below = k_top;
                                            if (k_below < klo_loc) k_below = klo_loc;
                                            if (k_below > khi_loc) k_below = khi_loc;
                                            const Set::Scalar alpha =
                                                u - uf;  // in [0, 1)
                                            T_z = (1.0 - alpha)
                                                  * T_a(i, j, k_above)
                                                +        alpha
                                                  * T_a(i, j, k_below);
                                        }
                                        // Thermoelastic relation under
                                        // 1-D confinement: the actual
                                        // sigma_xx is -E*beta*dT/(1-nu)
                                        // (negative = compressive when
                                        // T > T_ref). The K_I integrator
                                        // expects the compression-positive
                                        // value (= -sigma_xx_mechanics),
                                        // matching UpdateSpAfterMechanics's
                                        // "return -sigma(0, 0)" pattern.
                                        // So we negate twice → return the
                                        // positive E*beta*dT/(1-nu).
                                        //
                                        // Step 18: + p·nu/(1-nu)
                                        // confining-pressure contribution
                                        // (Kant Eq. 15). Captured as
                                        // sp_conf_offset outside the
                                        // MFIter; default 0 short-circuit
                                        // preserves byte-identity.
                                        return phase_E_loc * phase_beta_loc
                                               * (T_z - phase_T_ref_loc)
                                               / (1.0 - phase_nu_loc)
                                               + sp_conf_offset;
                                    };
                                    return Numeric::SpCriterion::K_I(
                                        sigma_at_depth, a_f,
                                        h_spall_n_seg_loc);
                                };

                                int  s_last_pass = -1;
                                bool found_failure = false;
                                for (int s = 0; s <= max_steps; ++s)
                                {
                                    const Set::Scalar z_crack =
                                        Set::Scalar(s) * dz;
                                    const int k_tip = k_top - s;
                                    if (k_tip < klo) break;

                                    Set::Scalar K_Ic_Pa =
                                        kic_pa_at_cell(k_tip);
                                    Set::Scalar K_I_val =
                                        k_i_at_zcrack(z_crack);

                                    if (K_I_val >= K_Ic_Pa)
                                    {
                                        s_last_pass = s;
                                        h_col       = z_crack;
                                        // Step 16b: edge-triggered event
                                        // count at (i, j, k_tip).
                                        if (cec(i, j, k_tip) < 0.5)
                                            cec(i, j, k_tip) = 1.0;
                                    }
                                    else
                                    {
                                        found_failure = true;
                                        break;
                                    }
                                }

                                // Step 16b: sub-cell binary search if the
                                // integer scan brackets the crossing
                                // between s_last_pass and s_last_pass + 1
                                // AND we have an in-box deeper cell to
                                // evaluate. Refines h_col to dz / 256
                                // (~6 um at dz = 1.5625 mm) — the GB/IG
                                // quantization fix from design note §4.
                                if (found_failure && s_last_pass >= 0 &&
                                    s_last_pass < max_steps &&
                                    (k_top - (s_last_pass + 1)) >= klo)
                                {
                                    Set::Scalar z_lo = Set::Scalar(s_last_pass)
                                                       * dz;       // pass
                                    Set::Scalar z_hi = Set::Scalar(s_last_pass + 1)
                                                       * dz;       // fail
                                    for (int b = 0; b < 8; ++b)
                                    {
                                        Set::Scalar z_mid =
                                            0.5 * (z_lo + z_hi);
                                        int k_mid = k_top -
                                            static_cast<int>(
                                                std::floor(z_mid / dz));
                                        if (k_mid < klo) k_mid = klo;
                                        if (k_mid > khi) k_mid = khi;
                                        Set::Scalar K_Ic_mid =
                                            kic_pa_at_cell(k_mid);
                                        Set::Scalar K_I_mid =
                                            k_i_at_zcrack(z_mid);
                                        if (K_I_mid >= K_Ic_mid)
                                            z_lo = z_mid;
                                        else
                                            z_hi = z_mid;
                                    }
                                    h_col = z_lo;
                                }

                                if (h_col > 0.0)
                                {
                                    spall_fires_col[idx] = 1;
                                    h_s_cand_col[idx]    = h_col;
                                    rop_s_cand_col[idx]  =
                                        (dt_local > 0.0) ? h_col / dt_local : 0.0;
                                    // Running max at the INITIAL top-z
                                    // slice (immovable; column anchor for
                                    // the Rossi depth-profile post-process).
                                    // Each column records the deepest h_col
                                    // ever observed across the simulation.
                                    const int top_k_dom = domain.bigEnd(2);
                                    if (top_k_dom >= klo && top_k_dom <= khi)
                                    {
                                        Set::Scalar prev =
                                            hsf(i, j, top_k_dom);
                                        if (h_col > prev) hsf(i, j, top_k_dom) = h_col;
                                    }
                                    // Step 16c: first-firing CSV event.
                                    // Detect the cec(i, j, k_top) 0 -> 1
                                    // transition from this step's integer
                                    // scan. The freeze guard inside the
                                    // scan already set cec to 1 if any
                                    // crossing passed; we use the pre-loop
                                    // snapshot to decide whether THIS step
                                    // is the column's first-firing step.
                                    if (cec_at_top_pre < 0.5)
                                    {
                                        HColEvent ev;
                                        ev.t      = time;
                                        ev.col_i  = i;
                                        ev.col_j  = j;
                                        ev.a_f    = a_f;
                                        ev.h_col  = h_col;
                                        ev.is_gb  = (is_gb_top(i, j, k_top)
                                                     >= 0.5) ? 1 : 0;
                                        ev.sp     = Sp_top(i, j, k_top);
                                        local_h_col_events.push_back(ev);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        // Step 14b: D sampled per the configured spall.sample
                        // mode (PASS 1.5 above). top_cell mode preserves the
                        // pre-Step-14b behaviour bit-for-bit because the scan
                        // window collapses to k = k_top.
                        Set::Scalar Dv = D_max_col[idx];
                        if (!std::isfinite(Dv)) Dv = 0.0;
                        if (Dv >= D_th)
                        {
                            Set::Scalar ratio = 0.0;
                            if (spall_prescribed)
                            {
                                ratio = ratio_f(xc, yc, zc, time);
                                if (!std::isfinite(ratio)) ratio = 0.0;
                            }
                            else
                            {
                                Set::Matrix sigma = Numeric::Interpolate::NodeToCellAverage(
                                    sig_node, i, j, k_top, 0);
                                ratio = PrincipalStressRatio(sigma);
                            }
                            if (ratio >= xi_s)
                            {
                                spall_fires_col[idx] = 1;
                                h_s_cand_col[idx]    = h_spall;
                                rop_s_cand_col[idx]  = (dt_local > 0.0) ? h_spall / dt_local : 0.0;
                            }
                        }
                    }
                }

                if (vapor_on && ptabs)
                {
                    Set::Scalar H_top = H(i, j, k_top);
                    Set::Scalar H_vap = ptabs[p_top].H_of_T(phases[p_top].T_vap_lo);
                    Set::Scalar Cp_p  = phases[p_top].Cp;
                    Set::Scalar rho_p = phases[p_top].rho;
                    Set::Scalar L_tot = Cp_p * (phases[p_top].T_vap_lo - T0_v)
                                      + phases[p_top].L_m + phases[p_top].L_v;
                    Set::Scalar rop_v = 0.0;
                    if (vapor_prescribed)
                    {
                        // Prescribed-RoP path stays as-is for synthetic
                        // unit tests (e.g. regime_low_high_power).
                        rop_v = vapor_rop_eval(xc, yc, zc, time);
                        if (!std::isfinite(rop_v) || rop_v < 0.0) rop_v = 0.0;
                    }
                    else if (sp_weibull_spall)
                    {
                        // Step 17: LOCAL v_vap from per-cell Gaussian
                        // beam power flux at the surface (W/m²). At the
                        // top face z = z0, the beam waist is w0 with no
                        // spreading:
                        //   Q_local(x, y, t) = (2·P(t) / (π·w0²))
                        //                       · exp(-2·r² / w0²)
                        // r² = (x - x0)² + (y - y0²). Compared to the
                        // bulk averaged form below (used under
                        // damage_law), the local v_vap at the beam
                        // centre equals the bulk value, and decays
                        // exponentially off-centre — matching the
                        // revised-approach §5 "interpreted locally
                        // rather than averaged" prescription.
                        if (std::isfinite(L_tot) && L_tot > 0.0 &&
                            rho_p > 0.0 && bw0_v > 0.0)
                        {
                            const Set::Scalar pi  = 3.141592653589793238;
                            const Set::Scalar r2  = (xc - bx0_v) * (xc - bx0_v)
                                                  + (yc - by0_v) * (yc - by0_v);
                            const Set::Scalar Q_local =
                                (2.0 * P_now / (pi * bw0_v * bw0_v))
                                * std::exp(-2.0 * r2 / (bw0_v * bw0_v));
                            rop_v = oneR * Q_local / (rho_p * L_tot);
                        }
                        if (!std::isfinite(rop_v) || rop_v < 0.0) rop_v = 0.0;
                    }
                    else
                    {
                        // damage_law: BULK-averaged RoP_vap path
                        // (byte-identical to Step 11). At the beam
                        // centre this matches the local Q above; off
                        // centre it incorrectly applies the centre
                        // value uniformly. damage_law tests
                        // (hu_end_to_end etc.) lock this behaviour in.
                        if (std::isfinite(L_tot) && L_tot > 0.0 && rho_p > 0.0)
                            rop_v = oneR * P_now / (rho_p * L_tot * A_beam);
                        if (!std::isfinite(rop_v) || rop_v < 0.0) rop_v = 0.0;
                    }
                    if (H_top >= H_vap && rop_v > 0.0)
                    {
                        vapor_fires_col[idx] = 1;
                        h_v_cand_col[idx]    = rop_v * dt_local;
                        rop_v_cand_col[idx]  = rop_v;
                    }
                }
            }
        }
        if (ncols > 0)
        {
            amrex::ParallelDescriptor::ReduceIntSum(spall_fires_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceIntSum(vapor_fires_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(h_s_cand_col.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(h_v_cand_col.data(),   ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rop_s_cand_col.data(), ncols);
            amrex::ParallelDescriptor::ReduceRealSum(rop_v_cand_col.data(), ncols);
        }

        // Step 16c: rank-ordered append of per-column first-firing events
        // to <plot_file>_h_col_events.csv. Each rank writes its own
        // subset; barriers serialize the file writes so the header (from
        // the I/O proc on the first emission) and subsequent rows don't
        // interleave. Append mode preserves cross-step accumulation; the
        // I/O proc latch ensures the header is emitted exactly once per
        // run (re-runs of the same plot_file truncate-and-rewrite).
        if (sp_weibull_spall && lev == 0)
        {
            const std::string csv_path =
                plot_file + "_h_col_events.csv";
            const int io_rank =
                amrex::ParallelDescriptor::IOProcessorNumber();
            const int my_rank = amrex::ParallelDescriptor::MyProc();
            const int n_ranks = amrex::ParallelDescriptor::NProcs();
            for (int r = 0; r < n_ranks; ++r)
            {
                if (my_rank == r)
                {
                    const bool need_header =
                        (r == io_rank) &&
                        !sp_h_col_events_csv_header_written;
                    std::ios::openmode mode = need_header
                                                  ? std::ios::trunc
                                                  : std::ios::app;
                    std::ofstream ofs(csv_path,
                                      mode | std::ios::out);
                    if (ofs.is_open())
                    {
                        if (need_header)
                        {
                            ofs << "time,col_i,col_j,a_f,h_col,"
                                   "is_gb_at_top,Sp_field_at_top\n";
                            sp_h_col_events_csv_header_written = true;
                        }
                        ofs.setf(std::ios::scientific);
                        ofs.precision(8);
                        for (auto const& ev : local_h_col_events)
                        {
                            ofs << ev.t << "," << ev.col_i << ","
                                << ev.col_j << "," << ev.a_f << ","
                                << ev.h_col << "," << ev.is_gb << ","
                                << ev.sp << "\n";
                        }
                    }
                }
                amrex::ParallelDescriptor::Barrier();
            }
        }

        // Per-column winner from reduced data (deterministic on every rank).
        // Step 17: under sp_weibull, spall has priority over vap when both
        // gates fire at the same column. Under damage_law, the
        // winner-takes-all on h (Step 11 behaviour) is preserved
        // byte-identical via the else-branch.
        amrex::Vector<int> regime_col(ncols, 0);
        amrex::Vector<Set::Scalar> h_chosen_col(ncols, 0.0);
        amrex::Vector<Set::Scalar> rop_chosen_col(ncols, 0.0);
        for (int idx = 0; idx < ncols; ++idx)
        {
            const bool sf = spall_fires_col[idx] != 0;
            const bool vf = vapor_fires_col[idx] != 0;
            const Set::Scalar hs = h_s_cand_col[idx];
            const Set::Scalar hv = h_v_cand_col[idx];
            int reg_choice = 0;
            Set::Scalar h_c = 0.0;
            if (sp_weibull_spall)
            {
                // Spall priority at the same (x, y) per revised-approach
                // §5: v_n = v_spall if Sp fires, else v_vap if H >= H_vap.
                if      (sf) { reg_choice = 1; h_c = hs; }
                else if (vf) { reg_choice = 2; h_c = hv; }
            }
            else
            {
                if (sf && vf)
                {
                    if (hv > hs) { reg_choice = 2; h_c = hv; }
                    else         { reg_choice = 1; h_c = hs; }
                }
                else if (sf) { reg_choice = 1; h_c = hs; }
                else if (vf) { reg_choice = 2; h_c = hv; }
            }
            regime_col[idx]     = reg_choice;
            h_chosen_col[idx]   = h_c;
            rop_chosen_col[idx] = (dt_local > 0.0) ? h_c / dt_local : 0.0;
        }

        // PASS 3: every rank applies the per-column phi shift to its
        // local cells and writes per-cell detach diagnostics. Both
        // tasks read from the reduced per-column buffers, so the
        // result is identical regardless of how the BoxArray splits z.
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar>       const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& rem = (*removed_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& reg = (*regime_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& ropU= (*RoP_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& T   = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LS  = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LL  = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LV  = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);

            amrex::Array4<Set::Scalar> ev_s, thk_s, rop_s;
            if (spall_on)
            {
                ev_s  = (*spall_event_mf[lev]).array(mfi);
                thk_s = (*spall_thickness_mf[lev]).array(mfi);
                rop_s = (*RoP_spall_mf[lev]).array(mfi);
            }
            amrex::Array4<Set::Scalar> ev_v, thk_v, rop_v;
            if (vapor_on)
            {
                ev_v  = (*vapor_event_mf[lev]).array(mfi);
                thk_v = (*vapor_thickness_mf[lev]).array(mfi);
                rop_v = (*RoP_vap_mf[lev]).array(mfi);
            }
            // Step 17: per-cell regime field under sp_weibull. Registered
            // only when sp_weibull_spall (gated above in Parse()), so
            // binding is conditional too.
            amrex::Array4<Set::Scalar> regf;
            if (sp_weibull_spall)
                regf = (*regime_field_mf[lev]).array(mfi);

            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);

            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                const int idx = col_index(i, j);
                const int reg_choice = regime_col[idx];
                if (reg_choice == 0) continue;

                const Set::Scalar h_c     = h_chosen_col[idx];
                const Set::Scalar rop_c   = rop_chosen_col[idx];
                const bool sf             = spall_fires_col[idx] != 0;
                const bool vf             = vapor_fires_col[idx] != 0;
                const Set::Scalar h_s_c   = h_s_cand_col[idx];
                const Set::Scalar h_v_c   = h_v_cand_col[idx];
                const Set::Scalar rop_s_c = rop_s_cand_col[idx];
                const Set::Scalar rop_v_c = rop_v_cand_col[idx];

                for (int k = klo; k <= khi; ++k)
                {
                    Set::Scalar phi_old = phi(i, j, k);
                    Set::Scalar phi_new = phi_old - h_c;
                    phi(i, j, k) = phi_new;
                    if (rem(i, j, k) > 0.5) continue;
                    if (phi_old >= 0.0 && phi_new < 0.0)
                    {
                        rem(i, j, k)  = 1.0;
                        reg(i, j, k)  = static_cast<Set::Scalar>(reg_choice);
                        ropU(i, j, k) = rop_c;
                        // Step 17: also write the cell-resolved
                        // regime_field under sp_weibull. Same value as
                        // the scalar `reg(i, j, k)` here — the
                        // distinction is the registered field name
                        // (`regime` vs `regime_field`) so plotfile
                        // consumers can pick the canonical per-cell
                        // diagnostic under sp_weibull.
                        if (sp_weibull_spall)
                            regf(i, j, k) =
                                static_cast<Set::Scalar>(reg_choice);
                        if (spall_on)
                        {
                            ev_s (i, j, k) = sf ? 1.0 : 0.0;
                            thk_s(i, j, k) = sf ? h_s_c : 0.0;
                            rop_s(i, j, k) = rop_s_c;
                        }
                        if (vapor_on)
                        {
                            ev_v (i, j, k) = vf ? 1.0 : 0.0;
                            thk_v(i, j, k) = vf ? h_v_c : 0.0;
                            rop_v(i, j, k) = rop_v_c;
                        }
                        D(i, j, k) = 0.0;
                        T(i, j, k) = T_amb_local;
                        if (ptabs)
                        {
                            int p = static_cast<int>(std::round(ph(i, j, k)));
                            if (p < 0 || p >= n_phases) p = 0;
                            Set::Scalar Hv = ptabs[p].H_of_T(T_amb_local);
                            Set::Scalar Tdummy, ls, ll, lv;
                            ptabs[p].Invert(Hv, Tdummy, ls, ll, lv);
                            H(i, j, k)  = Hv;
                            LS(i, j, k) = ls;
                            LL(i, j, k) = ll;
                            LV(i, j, k) = lv;
                        }
                        else
                        {
                            H(i, j, k)  = 0.0;
                            LS(i, j, k) = 1.0;
                            LL(i, j, k) = 0.0;
                            LV(i, j, k) = 0.0;
                        }
                    }
                }
            }
        }

        // Rebuild diagnostics that depend on the new mask/D field.
        UpdateSurfaceMaskFromPhi(lev);
        RefreshKappaEff(lev);
        (void)dz;
    }

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

    static bool FiniteVector(const Set::Vector& v)
    {
        for (int n = 0; n < AMREX_SPACEDIM; ++n)
            if (!std::isfinite(v(n))) return false;
        return true;
    }

    static bool FiniteMatrix(const Set::Matrix& m)
    {
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            for (int j = 0; j < AMREX_SPACEDIM; ++j)
                if (!std::isfinite(m(i, j))) return false;
        return true;
    }

    static bool FiniteModel(const Model::Solid::Affine::Isotropic& model)
    {
        return std::isfinite(model.ddw.Lambda())
            && std::isfinite(model.ddw.Mu())
            && FiniteMatrix(model.F0);
    }

    std::string AMRLocationLabel(int lev,
                                 const amrex::IntVect& iv,
                                 const amrex::Box& valid,
                                 const amrex::Box& domain) const
    {
        if (valid.contains(iv)) return "valid";
        if (!domain.contains(iv)) return "physical-boundary ghost";
        if (lev > 0) return "coarse-fine/box ghost";
        return "box ghost";
    }

    void AbortBadAMRValue(const char* stage,
                          const char* field,
                          int lev,
                          int comp,
                          const amrex::IntVect& iv,
                          const amrex::Box& valid,
                          const amrex::Box& domain) const
    {
        Util::Abort(INFO, "non-finite ", field, " during ", stage,
                    " at level ", lev, ", comp ", comp,
                    ", index (", iv[0], ",", iv[1], ",", iv[2], ") [",
                    AMRLocationLabel(lev, iv, valid, domain), "]");
    }

    void CheckScalarFieldFinite(const Set::Field<Set::Scalar>& field,
                                const char* name,
                                int lev,
                                const char* stage,
                                bool nodal) const
    {
        amrex::Box domain = geom[lev].Domain();
        if (nodal) domain.convert(amrex::IntVect::TheNodeVector());
        for (amrex::MFIter mfi(*field[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box valid = mfi.validbox();
            const amrex::Box scan = (*field[lev])[mfi].box();
            amrex::Array4<const Set::Scalar> const& arr = field[lev]->array(mfi);
            const int ncomp = field[lev]->nComp();
            amrex::LoopOnCpu(scan, [&] (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                for (int n = 0; n < ncomp; ++n)
                    if (!std::isfinite(arr(i, j, k, n)))
                        AbortBadAMRValue(stage, name, lev, n, iv, valid, domain);
            });
        }
    }

    void CheckVectorFieldFinite(const Set::Field<Set::Vector>& field,
                                const char* name,
                                int lev,
                                const char* stage) const
    {
        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        for (amrex::MFIter mfi(*field[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box valid = mfi.validbox();
            const amrex::Box scan = valid;
            amrex::Array4<const Set::Vector> const& arr = field[lev]->array(mfi);
            amrex::LoopOnCpu(scan, [&] (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                    if (!std::isfinite(arr(i, j, k)(n)))
                        AbortBadAMRValue(stage, name, lev, n, iv, valid, domain);
            });
        }
    }

    void CheckMatrixFieldFinite(const Set::Field<Set::Matrix>& field,
                                const char* name,
                                int lev,
                                const char* stage) const
    {
        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        for (amrex::MFIter mfi(*field[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box valid = mfi.validbox();
            const amrex::Box scan = valid;
            amrex::Array4<const Set::Matrix> const& arr = field[lev]->array(mfi);
            amrex::LoopOnCpu(scan, [&] (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const Set::Matrix& m = arr(i, j, k);
                for (int r = 0; r < AMREX_SPACEDIM; ++r)
                    for (int c = 0; c < AMREX_SPACEDIM; ++c)
                        if (!std::isfinite(m(r, c)))
                            AbortBadAMRValue(stage, name, lev,
                                             r * AMREX_SPACEDIM + c,
                                             iv, valid, domain);
            });
        }
    }

    void CheckModelFieldFinite(const char* stage) const
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box valid = mfi.validbox();
                amrex::Box scan = valid;
                scan.grow(1);
                scan &= domain;
                amrex::Array4<const Model::Solid::Affine::Isotropic> const& arr =
                    model_mf[lev]->array(mfi);
                amrex::LoopOnCpu(scan, [&] (int i, int j, int k) {
                    amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                    const auto& model = arr(i, j, k);
                    if (!std::isfinite(model.ddw.Lambda()) ||
                        std::abs(model.ddw.Lambda()) > 1.0e13)
                        AbortBadAMRValue(stage, "model_mf(lambda)", lev, 0,
                                         iv, valid, domain);
                    if (!std::isfinite(model.ddw.Mu()) ||
                        std::abs(model.ddw.Mu()) > 1.0e13)
                        AbortBadAMRValue(stage, "model_mf(mu)", lev, 1,
                                         iv, valid, domain);
                    for (int r = 0; r < AMREX_SPACEDIM; ++r)
                        for (int c = 0; c < AMREX_SPACEDIM; ++c)
                            if (!std::isfinite(model.F0(r, c)) ||
                                std::abs(model.F0(r, c)) > 1.0)
                                AbortBadAMRValue(stage, "model_mf(F0)", lev,
                                                 2 + r * AMREX_SPACEDIM + c,
                                                 iv, valid, domain);
                });
            }
        }
    }

    void CheckAMRMechanicsFieldsFinite(const char* stage, bool include_stress) const
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            CheckVectorFieldFinite(disp_mf, "disp_mf", lev, stage);
            if (microstructure_enabled)
            {
                CheckScalarFieldFinite(D_mf, "D_mf", lev, stage, /*nodal=*/false);
                CheckScalarFieldFinite(thermal_free_strain_mf,
                                       "thermal_free_strain_mf",
                                       lev, stage, /*nodal=*/false);
            }
            if (include_stress)
            {
                CheckMatrixFieldFinite(stress_mf, "stress_mf", lev, stage);
                CheckMatrixFieldFinite(strain_mf, "strain_mf", lev, stage);
            }
        }
        CheckModelFieldFinite(stage);
    }

    // §4.4 diagnostic: post-mechanics-solve residual analysis.
    // Gated on env var ALAMO_AMR_DIAG=1 + el.print_residual=1.
    // Reports per-level L2(res), per-z-slab L2 on finest level, and
    // rigid-body-mode projection of the finest-level residual.
    // See tests/MMWSpalling/hu_spall_onset/AMR_FAILURE_ANALYSIS.md §4.4.
    void DiagnoseAMRResidual(Set::Scalar a_time, int a_step) const
    {
        if (std::getenv("ALAMO_AMR_DIAG") == nullptr) return;
        if (!m_print_residual) return;
        if (finest_level == 0) return;

        std::stringstream ss;
        ss << "=== AMR residual diag STEP " << a_step
           << " (t=" << a_time << ") ===\n";

        // Per-level L2 norm and max
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            Set::Scalar l2sq = 0.0;
            Set::Scalar maxabs = 0.0;
            Set::Scalar maxx = 0.0, maxy = 0.0, maxz = 0.0;
            const auto* dx  = geom[lev].CellSize();
            const auto* plo = geom[lev].ProbLo();
            for (amrex::MFIter mfi(*res_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox();
                amrex::Array4<const Set::Vector> const& res =
                    res_mf[lev]->array(mfi);
                amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                    const Set::Vector& r = res(i, j, k);
                    Set::Scalar n2 = r.dot(r);
                    l2sq += n2;
                    Set::Scalar n = std::sqrt(n2);
                    if (n > maxabs) {
                        maxabs = n;
                        maxx = plo[0] + Set::Scalar(i) * dx[0];
                        maxy = plo[1] + Set::Scalar(j) * dx[1];
                        maxz = plo[2] + Set::Scalar(k) * dx[2];
                    }
                });
            }
            amrex::ParallelDescriptor::ReduceRealSum(l2sq);
            amrex::ParallelDescriptor::ReduceRealMax(maxabs);
            ss << "  L" << lev << ": L2(res)=" << std::sqrt(l2sq)
               << "  max|res|=" << maxabs
               << " at (" << maxx << "," << maxy << "," << maxz << ")\n";
        }

        // Per-z-slab L2 on finest level
        const int lev = finest_level;
        amrex::Box dom_lev = geom[lev].Domain();
        dom_lev.convert(amrex::IntVect::TheNodeVector());
        const int kmin = dom_lev.smallEnd(2);
        const int kmax = dom_lev.bigEnd(2);
        const int nk = kmax - kmin + 1;
        const auto* dx  = geom[lev].CellSize();
        const auto* plo = geom[lev].ProbLo();
        std::vector<Set::Scalar> z_l2sq(nk, 0.0);
        for (amrex::MFIter mfi(*res_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box bx = mfi.validbox();
            amrex::Array4<const Set::Vector> const& res =
                res_mf[lev]->array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                if (k >= kmin && k <= kmax) {
                    const Set::Vector& r = res(i, j, k);
                    z_l2sq[k - kmin] += r.dot(r);
                }
            });
        }
        for (int k = 0; k < nk; ++k)
            amrex::ParallelDescriptor::ReduceRealSum(z_l2sq[k]);
        ss << "  L" << lev << " z-slab L2(res):\n";
        for (int k = 0; k < nk; ++k) {
            Set::Scalar zc = plo[2] + Set::Scalar(kmin + k) * dx[2];
            ss << "    z=" << zc << "  L2=" << std::sqrt(z_l2sq[k]) << "\n";
        }

        // Rigid-body-mode projection on finest-level residual
        Set::Scalar dot_m[6] = {0,0,0,0,0,0};
        Set::Scalar nrm_m[6] = {0,0,0,0,0,0};
        Set::Scalar res_l2sq = 0.0;
        for (amrex::MFIter mfi(*res_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box bx = mfi.validbox();
            amrex::Array4<const Set::Vector> const& res =
                res_mf[lev]->array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                Set::Scalar x = plo[0] + Set::Scalar(i) * dx[0];
                Set::Scalar y = plo[1] + Set::Scalar(j) * dx[1];
                Set::Scalar z = plo[2] + Set::Scalar(k) * dx[2];
                const Set::Vector& r = res(i, j, k);
                res_l2sq += r.dot(r);
                // Translation modes (mode_norm per node = 1)
                dot_m[0] += r(0); nrm_m[0] += 1.0;
                dot_m[1] += r(1); nrm_m[1] += 1.0;
                dot_m[2] += r(2); nrm_m[2] += 1.0;
                // Rotation about x: u = (0, -z, y)
                dot_m[3] += -z * r(1) + y * r(2); nrm_m[3] += z*z + y*y;
                // Rotation about y: u = (z, 0, -x)
                dot_m[4] += z * r(0) - x * r(2);  nrm_m[4] += z*z + x*x;
                // Rotation about z: u = (-y, x, 0)
                dot_m[5] += -y * r(0) + x * r(1); nrm_m[5] += y*y + x*x;
            });
        }
        for (int m = 0; m < 6; ++m) {
            amrex::ParallelDescriptor::ReduceRealSum(dot_m[m]);
            amrex::ParallelDescriptor::ReduceRealSum(nrm_m[m]);
        }
        amrex::ParallelDescriptor::ReduceRealSum(res_l2sq);
        const char* names[6] = {"TX","TY","TZ","RX","RY","RZ"};
        ss << "  L" << lev << " RBM projection (||proj*mode||^2 / ||res||^2):\n";
        Set::Scalar total_rbm = 0.0;
        for (int m = 0; m < 6; ++m) {
            Set::Scalar proj = (nrm_m[m] > 0.0) ? dot_m[m] / nrm_m[m] : 0.0;
            Set::Scalar contrib = proj * proj * nrm_m[m];
            total_rbm += contrib;
            Set::Scalar pct = (res_l2sq > 0.0)
                ? 100.0 * contrib / res_l2sq : 0.0;
            ss << "    " << names[m]
               << ": proj=" << proj
               << "  contrib=" << contrib
               << "  (" << pct << "% of res^2)\n";
        }
        Set::Scalar tot_pct = (res_l2sq > 0.0)
            ? 100.0 * total_rbm / res_l2sq : 0.0;
        ss << "  Total RBM: " << total_rbm
           << " / ||res||^2=" << res_l2sq
           << " = " << tot_pct << "%\n";

        // §9.5 hypothesis test: dump disp_mf at Level_1's zhi physical
        // boundary (k=hi_z, the topmost valid nodes) and the ghost row
        // above (k=hi_z+1). If the ghost row is uninitialized/zero
        // while the valid row is meaningful, alamo's RealFillBoundary
        // is failing to fill physical-BC ghosts on Level_1.
        const int kbnd  = dom_lev.bigEnd(2);
        const int kghst = kbnd + 1;
        Set::Scalar disp_l2_at_kbnd  = 0.0;
        Set::Scalar disp_l2_at_kghst = 0.0;
        Set::Scalar disp_z_max_kbnd  = 0.0;
        Set::Scalar disp_z_max_kghst = 0.0;
        int n_kbnd = 0, n_kghst = 0;
        const int ngrow = disp_mf[lev]->nGrow();
        ss << "  L" << lev << " ghost-fill probe (nGrow=" << ngrow
           << "  k_valid_top=" << kbnd << "  k_ghost=" << kghst << "):\n";
        for (amrex::MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.validbox();
            bx.grow(amrex::IntVect(0, 0, ngrow));  // include ghosts in z
            amrex::Array4<const Set::Vector> const& d =
                disp_mf[lev]->const_array(mfi);
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
            const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
            for (int k = klo; k <= khi; ++k) {
                if (k != kbnd && k != kghst) continue;
                for (int j = jlo; j <= jhi; ++j)
                    for (int i = ilo; i <= ihi; ++i) {
                        const Set::Vector& v = d(i, j, k);
                        Set::Scalar n2 = v.dot(v);
                        Set::Scalar zmag = std::abs(v(2));
                        if (k == kbnd) {
                            disp_l2_at_kbnd += n2;
                            disp_z_max_kbnd = std::max(disp_z_max_kbnd, zmag);
                            ++n_kbnd;
                        } else {
                            disp_l2_at_kghst += n2;
                            disp_z_max_kghst = std::max(disp_z_max_kghst, zmag);
                            ++n_kghst;
                        }
                    }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(disp_l2_at_kbnd);
        amrex::ParallelDescriptor::ReduceRealSum(disp_l2_at_kghst);
        amrex::ParallelDescriptor::ReduceRealMax(disp_z_max_kbnd);
        amrex::ParallelDescriptor::ReduceRealMax(disp_z_max_kghst);
        amrex::ParallelDescriptor::ReduceIntSum(n_kbnd);
        amrex::ParallelDescriptor::ReduceIntSum(n_kghst);
        ss << "    k=" << kbnd << " (valid zhi top): n=" << n_kbnd
           << "  L2(disp)=" << std::sqrt(disp_l2_at_kbnd)
           << "  max|disp_z|=" << disp_z_max_kbnd << "\n"
           << "    k=" << kghst << " (ghost above zhi): n=" << n_kghst
           << "  L2(disp)=" << std::sqrt(disp_l2_at_kghst)
           << "  max|disp_z|=" << disp_z_max_kghst << "\n";

        Util::Message(INFO, ss.str());
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

        for (amrex::MFIter mfi(*flaw_a_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar>       const& fa = (*flaw_a_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& gb =
                (*is_grain_boundary_mf[lev]).array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                std::uint64_t key = SplitMix64(seed);
                key = SplitMix64(key ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)));
                key = SplitMix64(key ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(j)));
                key = SplitMix64(key ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(k)));
                // 53-bit mantissa -> U in the open interval (0, 1).
                Set::Scalar U = (static_cast<Set::Scalar>(key >> 11) + 0.5)
                                * (1.0 / 9007199254740992.0);
                Set::Scalar a_scale = (gb(i, j, k) > 0.5) ? a0_gb : a0_ig;
                fa(i, j, k) = a_scale * std::pow(-std::log(U), inv_m);
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

        bool enabled = false;
        Mode mode = Mode::PrescribedTemperature;
        Set::Scalar x0 = 0.0, y0 = 0.0, radius = 0.0;
        Set::Scalar T_ext   = 0.0;
        Set::Scalar h_conv  = 0.0;  // convective_flame: flame coefficient h_fl
        Set::Scalar T_flame = 0.0;  // convective_flame: flame temperature
        amrex::Parser            T_parser;
        amrex::ParserExecutor<1> T_f;

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
    bool debug_amr_mechanics = false;

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

    // Step 16: Sp-vs-depth h_spall scan parameters (sp_weibull branch only).
    // `spall.h_spall_max` caps the downward column scan in metres (safety
    // stop; in practice the K_I < K_Ic crossing always fires first). The
    // segment count `spall.h_spall_n_segments` falls back to `sp_n_segments`
    // when set to 0 (default), so the depth scan uses the same Simpson
    // resolution as `UpdateSpAfterMechanics`'s surface Sp.
    Set::Scalar spall_h_spall_max         = 0.010;
    int         spall_h_spall_n_segments  = 0;
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

    // Step 18: lithostatic confining pressure (Pa). Positive value =
    // compressive pre-stress on side faces (Kant convention). Enters the
    // K_I integrand via the compression-positive sigma_xx offset
    // +p·nu/(1-nu) per Kant Eq. 15 (kant-2017-central-aare-validation-
    // reference.md line 96). Default 0 = byte-identical to all
    // pre-Step-18 tests. Phase-0 nu is used (single-phase scope, same
    // as Step 16d's T-derived sigma_xx path); multi-phase Voronoi is
    // out of scope. Applied as a POST-SOLVE offset in BOTH sigma_xx
    // lambdas (UpdateSpAfterMechanics's Sp_field path AND Step 16d's
    // k_i_at_zcrack depth-scan path), NOT as a new BC on the elastic
    // solver — avoids reopening the MLMG limitations from Step 16
    // takeaway #1.
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
    Set::Field<Set::Scalar> Sp_cluster_id_mf;
    // Per-step cluster diagnostics CSV: written by the I/O processor only
    // when weibull is enabled and a cluster fires. Tracks the path so a
    // single trailing-newline header is written exactly once per output
    // directory (re-opening with std::ios::app keeps existing data on
    // restart-style runs).
    bool sp_clusters_csv_header_written = false;

    // Step 16c: per-column h_col first-firing event log. Same I/O-proc
    // header-latch pattern as the clusters CSV above; written to
    // `<plot_file>_h_col_events.csv`. One row per column per run (the
    // first step the column's cec(i, j, k_top) flips 0 -> 1). Subsequent
    // firings of the same column are freeze-guarded (cec already 1, no
    // re-emit). Columns: time, col_i, col_j, a_f, h_col, is_gb_at_top,
    // Sp_field_at_top.
    bool sp_h_col_events_csv_header_written = false;

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
