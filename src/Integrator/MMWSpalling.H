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
            if (value.spall_enabled)
            {
                // phi is a level-set state: phi = z_surface(x,y) - z_cell.
                // Material cells have phi >= 0; void cells (already removed)
                // have phi < 0. Standard AMReX FillPatch preserves it across
                // regrids.
                value.RegisterNewFab(value.phi_mf,             nullptr, 1, 0, "phi",             true);
                value.RegisterNewFab(value.removed_mf,         nullptr, 1, 0, "removed",         true);
                value.RegisterNewFab(value.spall_event_mf,     nullptr, 1, 0, "spall_event",     true);
                value.RegisterNewFab(value.spall_thickness_mf, nullptr, 1, 0, "spall_thickness", true);
                value.RegisterNewFab(value.RoP_spall_mf,       nullptr, 1, 0, "RoP_spall",       true);
            }
        }

        if (value.damage_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: damage.enabled requires microstructure.nphases");

        if (value.cohesive_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: cohesive.enabled requires microstructure.nphases");

        if (value.spall_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spall.enabled requires microstructure.nphases");
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
        if (spall_enabled) InitPhi(lev);
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
        const bool        sp_on  = surface_patch.enabled;
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
                        Set::Scalar epsn = Numeric::Interpolate::CellToNodeAverage(epsc, i, j, k, 0);
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
                        Set::Scalar En  = Numeric::Interpolate::CellToNodeAverage(Ec,   i, j, k, 0);
                        Set::Scalar mun = Numeric::Interpolate::CellToNodeAverage(muc,  i, j, k, 0);
                        if (damage_on) {
                            Set::Scalar Dn = Numeric::Interpolate::CellToNodeAverage(Dc, i, j, k, 0);
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
            const int sp_topk_node = geom[lev].Domain().bigEnd(2) + 1;
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<Model::Solid::Affine::Isotropic> const& model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& eta  = eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& temp = temp_old_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix F0 = Set::Matrix::Zero();
                    Set::Scalar Tavg = Numeric::Interpolate::CellToNodeAverage(temp, i, j, k, 0);
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
        // Step 10: spall detection/removal must happen BEFORE the H/Temp
        // swaps below so the reset H(T_amb)/T_amb values become the "old"
        // state consumed by the thermal advance step.
        UpdateSpallAfterCohesive(lev, time, dt);
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

        // Step 6d: prescribed-temperature circular patch on the +z face.
        const bool        sp_on   = surface_patch.enabled;
        const int         sp_topk = geom[lev].Domain().bigEnd(2);
        const Set::Scalar sp_x0   = surface_patch.x0;
        const Set::Scalar sp_y0   = surface_patch.y0;
        const Set::Scalar sp_r2   = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_T    = sp_on ? surface_patch.T_at(time) : 0.0;

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
                    if (dxp * dxp + dyp * dyp <= sp_r2)
                        div_kgrad += kc * (sp_T - Tc) * inv_dz2;
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
        const bool void_aware = spall_enabled;
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
            Set::Patch<const Set::Scalar> rem;
            if (void_aware) rem = removed_mf.Patch(lev, mfi);
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
                Set::Scalar Q = BeamSource(P_now, bw0, blam, bx0, by0, bz0,
                                           bR, alpha_f, x, y, z, Tc);

                Set::Scalar eps = (Tc <= T_eps_b) ? eps_lo : eps_hi;
                Set::Scalar q_loss = surf(i, j, k)
                                   * (eps * sig_SB * (Tc*Tc*Tc*Tc - T_a*T_a*T_a*T_a)
                                      + h_c * (Tc - T_a));

                // Conservative H update; T is recovered in InvertHtoT
                // using the per-phase Numeric::Material::Table.
                H(i, j, k) = H_old(i, j, k) + dt * (div_kgrad + Q - q_loss * inv_dz);
                keff_out(i, j, k) = kc; // diagnostic, refreshed each step
            });
        }
    }

    // Cell-centred Beer-Lambert + Gaussian volumetric source (eqs. 1-3 of
    // main.tex), scaled by (1 − R) for surface reflectivity. Returns 0
    // above the surface (z > z_0).
    AMREX_GPU_HOST_DEVICE
    static Set::Scalar BeamSource(Set::Scalar P, Set::Scalar w0, Set::Scalar lam,
                                  Set::Scalar x0, Set::Scalar y0, Set::Scalar z0,
                                  Set::Scalar R,
                                  amrex::ParserExecutor<1> alpha_f,
                                  Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                  Set::Scalar T)
    {
        Set::Scalar dz = z0 - z;
        if (dz < 0.0 || P == 0.0) return 0.0;
        Set::Scalar zR = M_PI * w0 * w0 / lam;
        Set::Scalar w  = w0 * std::sqrt(1.0 + (dz / zR) * (dz / zR));
        Set::Scalar r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        Set::Scalar Pi = (2.0 * P / (M_PI * w0 * w0)) * (w0 / w) * (w0 / w)
                       * std::exp(-2.0 * r2 / (w * w));
        Set::Scalar a  = alpha_f(T);
        return (1.0 - R) * Pi * a * std::exp(-a * dz);
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
    }

    // Step 10: principal-stress ratio for the production spall path. ALAMO
    // stores tensile stress as positive, so the most compressive principal
    // stress is the smallest (most negative) eigenvalue. We recast to
    // compressive-positive (sigma_1 = -e_min, sigma_3 = -e_max) and require
    // both to be compressive (>0) for a meaningful Hoek-Brown ratio. When
    // the state is tensile, mixed, or numerically degenerate, the ratio is
    // returned as zero so spall is suppressed. Step 10 unit-tests use the
    // prescribed-ratio hook; this code path is documented but not exercised
    // by Step 10 regressions.
    static AMREX_FORCE_INLINE
    Set::Scalar PrincipalStressRatio(const Set::Matrix& sigma)
    {
        Eigen::SelfAdjointEigenSolver<Set::Matrix> es(sigma);
        if (es.info() != Eigen::Success) return 0.0;
        const auto& e = es.eigenvalues();
        Set::Scalar e_min = e(0);
        Set::Scalar e_max = e(e.size() - 1);
        Set::Scalar sigma_1 = -e_min;
        Set::Scalar sigma_3 = -e_max;
        if (!(sigma_3 > 1.0e-30) || !(sigma_1 > 0.0)) return 0.0;
        Set::Scalar r = sigma_1 / sigma_3;
        return std::isfinite(r) ? r : 0.0;
    }

    void UpdateSpallAfterCohesive(int lev, Set::Scalar time, Set::Scalar dt)
    {
        if (!spall_enabled) return;
        Util::Assert(INFO, TEST(microstructure_enabled));

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz   = DX[2];
        const Set::Scalar D_th = spall_damage_threshold;
        const Set::Scalar xi_s = spall_xi_s;
        const Set::Scalar C_h  = spall_C_h;
        const Set::Scalar alpha_eff   = spall_alpha_override ? spall_alpha   : therm_alpha;
        const Set::Scalar t_spall_eff = spall_t_spall_override ? spall_t_spall : (time + dt);
        if (!(std::isfinite(alpha_eff) && alpha_eff > 0.0))
            Util::Abort(INFO, "spall: thermal diffusivity for h_spall is non-positive (",
                        alpha_eff, "); set spall.alpha or hc.heat.alpha");
        if (!(std::isfinite(t_spall_eff) && t_spall_eff > 0.0))
            Util::Abort(INFO, "spall: t_spall for h_spall is non-positive (", t_spall_eff,
                        "); set spall.t_spall");
        const Set::Scalar h_spall = C_h * std::sqrt(alpha_eff * t_spall_eff);
        if (!(std::isfinite(h_spall) && h_spall >= 0.0))
            Util::Abort(INFO, "spall: computed h_spall is non-finite (", h_spall, ")");

        const bool prescribed = spall_prescribed_ratio_enabled;
        auto ratio_f = spall_ratio_f;
        const bool have_stress = !prescribed; // production path reads stress_mf
        if (have_stress) Util::RealFillBoundary(*stress_mf[lev], geom[lev]);

        // Per-step diagnostic fields are reset here. `phi` and `removed`
        // carry state across steps and are not reset.
        spall_event_mf[lev]->setVal(0.0);
        spall_thickness_mf[lev]->setVal(0.0);
        RoP_spall_mf[lev]->setVal(0.0);

        const Numeric::Material::Table *ptabs =
            phase_tables.empty() ? nullptr : phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        const Set::Scalar T_amb_local = T_amb;
        const Set::Scalar dt_local = dt;

        // Single-level, CPU-side per-column scan. Box-decomposition aware:
        // each MFIter box is searched for its top-most material cell in
        // each (i,j) column inside the box, so multi-box decompositions in
        // x/y are fine. The unit test runs single-rank single-box, but the
        // code is correct for column-coherent box partitions.
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar>       const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& rem = (*removed_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& ev  = (*spall_event_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& thk = (*spall_thickness_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& rop = (*RoP_spall_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& H   = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& T   = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LS  = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LL  = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& LV  = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& ph  = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> sig_node;
            if (have_stress) sig_node = (*stress_mf[lev]).array(mfi);

            const int klo = bx.smallEnd(2);
            const int khi = bx.bigEnd(2);
            const int ilo = bx.smallEnd(0);
            const int ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1);
            const int jhi = bx.bigEnd(1);

            for (int j = jlo; j <= jhi; ++j)
            for (int i = ilo; i <= ihi; ++i)
            {
                // Topmost material cell in this column inside this box.
                int k_top = -1;
                for (int k = khi; k >= klo; --k)
                {
                    if (rem(i, j, k) <= 0.5 && phi(i, j, k) >= 0.0) { k_top = k; break; }
                }
                if (k_top < 0) continue;

                Set::Scalar Dv = D(i, j, k_top);
                if (!std::isfinite(Dv)) Dv = 0.0;
                if (Dv < D_th) continue;

                Set::Scalar ratio = 0.0;
                Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                Set::Scalar zc = PLO[2] + (Set::Scalar(k_top) + 0.5) * DX[2];
                if (prescribed)
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
                if (ratio < xi_s) continue;

                // Event! Apply h_spall to the entire column inside this box.
                for (int k = klo; k <= khi; ++k)
                {
                    Set::Scalar phi_old = phi(i, j, k);
                    Set::Scalar phi_new = phi_old - h_spall;
                    phi(i, j, k) = phi_new;
                    if (rem(i, j, k) > 0.5) continue; // already void
                    if (phi_old >= 0.0 && phi_new < 0.0)
                    {
                        // Detach this cell.
                        rem(i, j, k) = 1.0;
                        ev (i, j, k) = 1.0;
                        thk(i, j, k) = h_spall;
                        rop(i, j, k) = (dt_local > 0.0) ? h_spall / dt_local : 0.0;
                        D(i, j, k)   = 0.0;
                        T(i, j, k)   = T_amb_local;
                        if (ptabs)
                        {
                            int p = static_cast<int>(std::round(ph(i, j, k)));
                            if (p < 0 || p >= n_phases) p = 0;
                            Set::Scalar Hv = ptabs[p].H_of_T(T_amb_local);
                            Set::Scalar Tdummy, ls, ll, lv;
                            ptabs[p].Invert(Hv, Tdummy, ls, ll, lv);
                            H(i, j, k) = Hv;
                            LS(i, j, k) = ls;
                            LL(i, j, k) = ll;
                            LV(i, j, k) = lv;
                        }
                        else
                        {
                            H(i, j, k) = 0.0;
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
        if (spall_enabled)
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
            amrex::Array4<Set::Scalar> const& ev  = (*spall_event_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& th  = (*spall_thickness_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rop = (*RoP_spall_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar zc = PLO[2] + (Set::Scalar(k) + 0.5) * DX[2];
                phi(i, j, k) = z_top - zc;
                rem(i, j, k) = 0.0;
                ev (i, j, k) = 0.0;
                th (i, j, k) = 0.0;
                rop(i, j, k) = 0.0;
            });
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

    // Step 6d: prescribed-temperature surface patch driver. A circular patch
    // on the +z face holds T = T_expr(t) and injects the Dirichlet face flux
    //   q = kappa * (T_patch - T_cell) / dz
    // as a divergence contribution kappa * (T_patch - T_cell) / dz^2 on
    // top-row cells whose centre falls inside the patch. Outside the patch,
    // the existing zhi BC continues to drive the face flux (typically
    // Neumann/adiabatic for the Hu validation). Disabled by default.
    struct SurfacePatch
    {
        static constexpr const char *name = "surface_patch";

        bool enabled = false;
        Set::Scalar x0 = 0.0, y0 = 0.0, radius = 0.0;
        Set::Scalar T_ext  = 0.0;
        Set::Scalar h_conv = 0.0;
        amrex::Parser            T_parser;
        amrex::ParserExecutor<1> T_f;

        Set::Scalar T_at(Set::Scalar t) const { return T_f(t); }

        static void Parse(SurfacePatch &v, IO::ParmParse &pp)
        {
            int en = 0;
            pp.query_default("enabled", en, 0);
            v.enabled = (en != 0);
            if (!v.enabled) return;
            std::string type;
            pp.query_default("type", type, std::string("prescribed_temperature"));
            if (type != "prescribed_temperature")
                Util::Abort(INFO, "surface_patch.type must be 'prescribed_temperature', got: ", type);
            std::string face;
            pp.query_default("face", face, std::string("zhi"));
            if (face != "zhi")
                Util::Abort(INFO, "surface_patch.face must be 'zhi' (only +z patch implemented), got: ", face);
            pp.query_required("x0",     v.x0);
            pp.query_required("y0",     v.y0);
            pp.query_required("radius", v.radius);
            std::string T_expr;
            pp.query_required("T_expr", T_expr);
            v.T_parser = amrex::Parser(T_expr);
            v.T_parser.registerVariables({"t"});
            v.T_f = v.T_parser.compile<1>();
            // T_ext / h_conv are reserved for an optional non-patch
            // convective exterior; not consumed by this step.
            pp.query_default("T_ext",  v.T_ext,  0.0);
            pp.query_default("h_conv", v.h_conv, 0.0);
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
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> removed_mf;
    Set::Field<Set::Scalar> spall_event_mf;
    Set::Field<Set::Scalar> spall_thickness_mf;
    Set::Field<Set::Scalar> RoP_spall_mf;
};
} // namespace Integrator
#endif
