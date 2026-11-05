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
#include "Numeric/Material/Grimsel.H"
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
// P1g-2c deliverable 1: the lateral thermal shift's write-back collects
// per-cell deltas in a small ordered map before touching the field, so a cell
// written by two rows is composed instead of overwritten. Removal.H is included
// INSIDE the class body, so this cannot live there.
#include <map>
#include <cstdlib>   // D2ab: std::strtod, for side_face_T = live
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
        // ⚠⚠ P1g-3 deliverable 7: THE TWO LEN FALLBACK COUNTERS, REGISTERED HERE
        // AND THE PLACEMENT TOOK FOUR ATTEMPTS. Registration must be BOTH after
        // the keys are parsed AND outside every mode guard:
        //   - inside `if (sp_lateral_drive)` (within the side-face-flux guard):
        //     unreachable on p0_steady_column, which has no surface patch
        //   - inside `if (weibull_enabled && ...)`: unreachable there too, that
        //     leg sets weibull.enabled = 0
        //   - at the top of Parse beside the beam counters: runs BEFORE
        //     ParseSpallSettings (:158), so the key was still false
        // p0_steady_column is the ONLY leg in this step that fires vertically and
        // therefore the only one where these counters can move, so each wrong
        // placement showed up as "column absent", which is indistinguishable
        // from "counter reads zero" unless the header is checked.
        // ⚠ Gated on the LENGTH key, so a keys-off run's thermo.dat column set is
        // unchanged and (b1)'s 0-of-2639 cannot fail for a reporting reason.
        if (value.sp_flake_len_flux)
        {
            value.RegisterIntegratedVariable(&value.sp_len_fb_h0,
                                             "len_fb_h0", false);
            value.RegisterIntegratedVariable(&value.sp_len_fb_q0,
                                             "len_fb_q0", false);
            value.RegisterIntegratedVariable(&value.lat_len_fb_h0_n,
                                             "lat_len_fb_h0", false);
            value.RegisterIntegratedVariable(&value.lat_len_fb_q0_n,
                                             "lat_len_fb_q0", false);
        }

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
                    if (plaw == "grimsel")
                        value.phases[i].grimsel_law = true;
                    else if (plaw != "constant")
                        Util::Abort(INFO, "microstructure.phase", i,
                                    ".law must be 'constant' or 'grimsel', got: ",
                                    plaw);
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
                // P0: registered LAST so every pre-P0 plotfile variable keeps
                // its position, and only under the key.
                if (value.sp_surface_T_recon)
                {
                    value.RegisterNewFab(value.T_surf_field_mf, nullptr, 1, 0,
                                         "T_surf_field", true);
                    value.RegisterNewFab(value.ell_eff_field_mf, nullptr, 1, 0,
                                         "ell_eff_field", true);
                }
                if (value.sp_surface_T_recon)
                {
                    value.RegisterIntegratedVariable(&value.recon_cols,
                                                     "recon_cols", false);
                    value.RegisterIntegratedVariable(&value.recon_flat,
                                                     "recon_flat", false);
                    value.RegisterIntegratedVariable(&value.recon_ell_clamped,
                                                     "recon_ell_clamped", false);
                    value.RegisterIntegratedVariable(&value.recon_T_clamped,
                                                     "recon_T_clamped", false);
                    value.RegisterIntegratedVariable(&value.recon_fallback,
                                                     "recon_fallback", false);
                    value.RegisterIntegratedVariable(&value.recon_hf_capped,
                                                     "recon_hf_capped", false);
                    value.RegisterIntegratedVariable(&value.recon_hf_bounded,
                                                     "recon_hf_bounded", false);
                    value.RegisterIntegratedVariable(&value.recon_dT_mean,
                                                     "recon_dT_mean", false);
                }
                // P1b: the acceptance test for the support rule is that
            }
        }

        if (value.damage_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: damage.enabled requires microstructure.nphases");

        if (value.cohesive_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: cohesive.enabled requires microstructure.nphases");

        if (value.spall_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spall.enabled requires microstructure.nphases");

        // ⚠⚠ P2-2d d7: CROSS-KEY, AND IT HAS TO BE HERE. `grimsel_table` reads
        // K_Ic0 and the at-temperature ratio off phase 0's law object, so it
        // needs a phase to read. `ParseSpallationModel` (:206) runs BEFORE
        // `phases` is resized (:293), so the check belongs with the other
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
            if (value.ledger_enabled && value.max_level > 0)
                Util::Abort(INFO, "energy_ledger.enabled is single-level only "
                                  "(amr.max_level = 0)");
            // AMR (2026-09-21): surface-band refinement. The finest level owns
            // the surface (SurfaceOwner): it alone applies the surface patch,
            // the losses, the beam and the removal pass; coarser levels run
            // conduction only and take the surface state by average-down.
            // TagCellsForRefinement keeps every column's top solid cell, the
            // refine_height above it and the refine_depth below it on the
            // finest level, so the band follows the receding surface.
            pp.query_default("surface.refine_depth",  value.surface_refine_depth,  0.0);
            pp.query_default("surface.refine_height", value.surface_refine_height, 0.0);
            if (value.max_level > 0 && value.removal_enabled())
            {
                if (!(value.surface_refine_depth > 0.0))
                    Util::Abort(INFO, "amr.max_level > 0 with spall/vapor removal requires "
                                      "surface.refine_depth > 0 (the surface band must live "
                                      "on the finest level)");
                for (int l = 1; l <= value.max_level; ++l)
                    if (value.nsubsteps[l] != 1)
                        Util::Abort(INFO, "amr.max_level > 0 with spall/vapor removal requires "
                                          "amr.nsubsteps = 1 (the surface closures advance "
                                          "once per step on the finest level)");
            }
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
                if (value.eps_rad_hi < value.eps_rad)
                    Util::Abort(INFO, "surface_patch.robin_form = pinned requires "
                                      "losses.epsilon_high >= losses.epsilon (got ", value.eps_rad_hi,
                                      " < ", value.eps_rad, "): the face-temperature Newton assumes a "
                                      "monotone residual with a unique root");
                value.RegisterIntegratedVariable(&value.pinned_cols_pinned, "pinned_cols_pinned", false);
                value.RegisterIntegratedVariable(&value.pinned_cols_face,   "pinned_cols_face",   false);
                value.RegisterIntegratedVariable(&value.pinned_q_mean,      "pinned_q_mean",      false);
            }
            // P0b deliverable 1: robin_form = surface. Same prerequisites as
            // pinned -- it is the same closure with the surface temperature
            // computed rather than remembered -- PLUS the reconstruction
            // itself, without which the form has no surface to sit on.
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
            {
                if (!value.microstructure_enabled)
                    Util::Abort(INFO, "surface_patch.robin_form = surface is implemented "
                                      "on the microstructure path only");
                if (!value.surface_follow_mask)
                    Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                      "surface.follow_mask = 1");
                if (!(value.spall_enabled &&
                      value.spallation_model == SpallationModel::SpWeibull))
                    Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                      "spall.enabled = 1 and spallation.model = sp_weibull");
                if (value.eps_rad_hi < value.eps_rad)
                    Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                      "losses.epsilon_high >= losses.epsilon (got ",
                                      value.eps_rad_hi, " < ", value.eps_rad, ")");
                if (!value.sp_surface_T_apply)
                    Util::Abort(INFO, "surface_patch.robin_form = surface requires "
                                      "spall.surface_T = reconstructed: the form evaluates "
                                      "the flux AND the losses at the reconstructed surface "
                                      "temperature, and without the reconstruction there is "
                                      "no surface temperature to evaluate them at.");
                // Deliverable 1.4: pinned_idle_cycles governs the pin rule,
                // which this form replaces. Say so rather than ignore it.
                if (value.surface_patch.pinned_idle_cycles != 2.0)
                    Util::Message(INFO, "surface_patch.pinned_idle_cycles = ",
                                  value.surface_patch.pinned_idle_cycles,
                                  " is INERT under robin_form = surface: the pin rule it "
                                  "governs is replaced by the per-column reconstruction.");
                value.RegisterIntegratedVariable(&value.surf_cols_robin, "surf_cols_robin", false);
                value.RegisterIntegratedVariable(&value.surf_q_mean,     "surf_q_mean",     false);
                value.RegisterIntegratedVariable(&value.surf_side_share, "surf_side_share", false);
                value.RegisterIntegratedVariable(&value.surf_side_P_top,
                                                 "surf_side_P_top",    false);
                value.RegisterIntegratedVariable(&value.surf_side_P_buried,
                                                 "surf_side_P_buried", false);
                value.RegisterIntegratedVariable(&value.surf_side_P_rin,
                                                 "surf_side_P_rin",  false);
                value.RegisterIntegratedVariable(&value.surf_side_P_rout,
                                                 "surf_side_P_rout", false);
            }
            // A1: removal-event log validity.
            if (value.sp_emit_removal_events_csv)
            {
                if (!(value.spall_enabled &&
                      value.spallation_model == SpallationModel::SpWeibull))
                    Util::Abort(INFO, "spall.removal_events_csv = 1 requires spall.enabled = 1 "
                                      "and spallation.model = sp_weibull");
            }
            // P0 SCOPE, part 2: spallation.model is known here (Parse():174
            // runs after ParseSpallSettings at :158). See the part-1 note.
            if (value.sp_surface_T_recon || value.sp_removal_energy_flake)
            {
                if (value.spallation_model != SpallationModel::SpWeibull)
                    Util::Abort(INFO,
                        "spall.surface_T / spall.removal_energy apply to the "
                        "sp_weibull criterion only; got spallation.model = "
                        "damage_law. The damage-law path computes its removal "
                        "thickness from a calibrated ratio, where a multi-cell "
                        "flake is the intended answer (unit/spall_event "
                        "deliberately removes TWO cells per event).");
            }
            // P0b needs material removal to mean anything at all.
            if (value.sp_removal_energy_flake && !value.spall_enabled)
                Util::Abort(INFO, "spall.removal_energy = flake requires "
                                  "spall.enabled = 1");
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
                // P0b: appended last, and ALSO gated on the reconstruction --
                // a tree input with the ledger on but no P0b key must not gain
                // a thermo.dat column, or the 2639-file identity breaks on the
                // header alone whatever the numbers do.
                if (value.sp_surface_T_recon)
                {
                    value.RegisterIntegratedVariable(&value.ledger_E_flake,
                                                     "ledger_E_flake", false);
                    value.RegisterIntegratedVariable(&value.ledger_E_H_solid,
                                                     "ledger_E_H_solid", false);
                    value.RegisterIntegratedVariable(&value.ledger_err_solid,
                                                     "ledger_err_solid", false);
                }
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
                // P2-3: gated on the march's own key and appended after every
                // pre-P2-3 jet column, so a leg without the key gains nothing
                // and no existing column moves. ⚠ NOT named `jet_T_mouth`:
                // that column already exists under jet_free_surface with a
                // different meaning (D2c dilution).
                if (sp.jet_wall_march)
                {
                    value.RegisterIntegratedVariable(&value.jet_P_wall,
                                                     "jet_P_wall", false);
                    value.RegisterIntegratedVariable(&value.wall_T_mouth,
                                                     "wall_T_mouth", false);
                    value.RegisterIntegratedVariable(&value.wall_bins,
                                                     "wall_bins", false);
                    value.RegisterIntegratedVariable(&value.wall_L_fit,
                                                     "wall_L_fit", false);
                }
            }
            if (value.surface_patch.enabled && value.surface_patch.feet &&
                value.surface_patch.foot_body_clearance)
            {
                if (!value.removal_enabled())
                    Util::Abort(INFO, "surface_patch.foot_body_clearance = 1 requires removal "
                                      "(spall.enabled or vapor.enabled)");
                value.RegisterIntegratedVariable(&value.foot_mech_vol,  "foot_mech_vol",  false);
                value.RegisterIntegratedVariable(&value.foot_mech_cols, "foot_mech_cols", false);
                // P2-2f: appended after foot_mech_cols and gated on the SAME
                // key, so no input without foot_body_clearance = 1 gains a
                // thermo.dat column and the 2639-file identity is untouched.
                value.RegisterIntegratedVariable(&value.foot_mech_E,    "foot_mech_E",    false);
            }
            // D2k: side-face diagnostics, registered LAST and only when the key
            // is on, so every pre-D2k thermo.dat column keeps its position and
            // the reference set stays byte-identical with the key off.
            if (value.surface_patch.enabled && value.surface_patch.side_face_flux)
            {
                if (!value.removal_enabled())
                    Util::Abort(INFO, "surface_patch.side_face_flux = 1 requires removal "
                                      "(spall.enabled or vapor.enabled): with no void cells "
                                      "there are no exposed vertical faces");
                if (!value.surface_follow_mask)
                    Util::Abort(INFO, "surface_patch.side_face_flux = 1 requires "
                                      "surface.follow_mask = 1");
                value.RegisterIntegratedVariable(&value.side_P_face,  "side_P_face",  false);
                // P1e: registered only under the key, so a keys-off run keeps
                // its exact thermo column set and the frozen tree stays identical.
                // P1f: the lateral mode's own counters, registered only under the
                // key. ⚠ lat_cand vs lat_fires separates "nothing qualified" from
                // "nothing fired", and lat_no_faceT is the runtime witness that the
                // face temperature actually arrived (inertness case 1).
                if (value.sp_lateral_drive)
                {
                    value.RegisterIntegratedVariable(&value.lat_cand_n,
                                                     "lat_cand", false);
                    value.RegisterIntegratedVariable(&value.lat_fires,
                                                     "lat_fires", false);
                    value.RegisterIntegratedVariable(&value.lat_voided_n,
                                                     "lat_voided", false);
                    value.RegisterIntegratedVariable(&value.lat_no_faceT_n,
                                                     "lat_no_faceT", false);
                    value.RegisterIntegratedVariable(&value.lat_ell_med,
                                                     "lat_ell_med", false);
                    // deliverable 12. ⚠ Inside the SAME `sp_lateral_drive`
                    // guard, so a keys-off run keeps its exact thermo.dat
                    // column set and no re-pin is owed -- asserted, not assumed.
                    value.RegisterIntegratedVariable(&value.lat_qmax,
                                                     "lat_qmax", false);
                    value.RegisterIntegratedVariable(&value.lat_Tface_max,
                                                     "lat_Tface_max", false);
                }
                // P1g-1: nested INSIDE the drive guard above is not enough --
                // these must appear only when the banking law is actually on,
                // or an OFF leg's thermo.dat column set moves and (a1) fails
                // for a reporting reason rather than a physics one.
                if (value.sp_lateral_thermal)
                {
                    value.RegisterIntegratedVariable(&value.lat_shift_n,
                                                     "lat_shift", false);
                    value.RegisterIntegratedVariable(&value.lat_shift_refuse,
                                                     "lat_shift_refuse", false);
                    value.RegisterIntegratedVariable(&value.lat_shift_E,
                                                     "lat_shift_E", false);
                    value.RegisterIntegratedVariable(&value.lat_shift_E_est,
                                                     "lat_shift_E_est", false);
                    value.RegisterIntegratedVariable(&value.lat_shift_ratio,
                                                     "lat_shift_ratio", false);
                    value.RegisterIntegratedVariable(&value.lat_overlap_n,
                                                     "lat_overlap", false);
                    value.RegisterIntegratedVariable(&value.lat_orphan_n,
                                                     "lat_orphan", false);
                    // P1g-2c: the truth diagnostic, the coverage counter and the
                    // three-way refusal split. ⚠ All inside the SAME
                    // sp_lateral_thermal guard, so a leg with the key off keeps
                    // its exact column set and the identity gate is untouched.
                    value.RegisterIntegratedVariable(&value.lat_shift_real,
                                                     "lat_shift_real", false);
                    value.RegisterIntegratedVariable(&value.lat_shift_realfrac,
                                                     "lat_shift_realfrac", false);
                    value.RegisterIntegratedVariable(&value.lat_touch_multi,
                                                     "lat_touch_multi", false);
                    value.RegisterIntegratedVariable(&value.lat_twoaxis_n,
                                                     "lat_twoaxis", false);
                    value.RegisterIntegratedVariable(&value.lat_ref_ell,
                                                     "lat_ref_ell", false);
                    value.RegisterIntegratedVariable(&value.lat_ref_void,
                                                     "lat_ref_void", false);
                    value.RegisterIntegratedVariable(&value.lat_ref_short,
                                                     "lat_ref_short", false);
                }
                if (value.sp_lateral_subcell)
                {
                    value.RegisterIntegratedVariable(&value.lat_bank_n,
                                                     "lat_bank", false);
                    value.RegisterIntegratedVariable(&value.lat_bank_hf,
                                                     "lat_bank_hf", false);
                    value.RegisterIntegratedVariable(&value.lat_accum_sum,
                                                     "lat_accum_sum", false);
                    value.RegisterIntegratedVariable(&value.lat_void_colmax,
                                                     "lat_void_colmax", false);
                }
                if (value.surface_patch.side_face_staircase)
                {
                    value.RegisterIntegratedVariable(&value.side_sc_mean,
                                                     "side_sc_mean", false);
                    value.RegisterIntegratedVariable(&value.side_sc_min,
                                                     "side_sc_min", false);
                    value.RegisterIntegratedVariable(&value.side_sc_cols,
                                                     "side_sc_cols", false);
                    value.RegisterIntegratedVariable(&value.side_P_in35,
                                                     "side_P_in35", false);
                }
                value.RegisterIntegratedVariable(&value.side_P_loss,  "side_P_loss",  false);
                value.RegisterIntegratedVariable(&value.side_faces,   "side_faces",   false);
                // D2r-1: registered ONLY when side_face_factor_above is set.
                // A new thermo column changes thermo.dat, and the frozen A_f02
                // baseline must stay byte-comparable (Goal item 1), so a run
                // that does not set the key keeps its exact column set.
                if (value.surface_patch.side_face_factor_above >= 0.0)
                {
                    value.RegisterIntegratedVariable(&value.side_P_face_above,
                                                     "side_P_face_above", false);
                    value.RegisterIntegratedVariable(&value.side_cols_above,
                                                     "side_cols_above",   false);
                    value.RegisterIntegratedVariable(&value.wall_step_q, "wall_step_q", false);
                    value.RegisterIntegratedVariable(&value.wall_step_n, "wall_step_n", false);
                }
                // D2x: registered ONLY under side_face_area_exact, for the same
                // reason -- a new thermo column changes thermo.dat and would
                // break byte comparison against the frozen baselines.
                if (value.surface_patch.side_face_area_exact)
                    value.RegisterIntegratedVariable(&value.side_w_mean, "side_w_mean", false);
                // D2z: registered ONLY under the step key, so every other
                // thermo.dat header is untouched and the d2z1.md5 tree
                // reference stays valid for runs without it.
                // ⚠ None of these can read only one value: delta_mean is the
                // measured step height, delta_fallback_cols counts a named
                // population, and riser_area_ratio compares the two area laws.
                if (value.surface_patch.side_face_step_height)
                {
                    value.RegisterIntegratedVariable(&value.delta_mean,
                                                     "delta_mean", false);
                    value.RegisterIntegratedVariable(&value.delta_fallback_cols,
                                                     "delta_fallback_cols", false);
                    value.RegisterIntegratedVariable(&value.riser_area_ratio,
                                                     "riser_area_ratio", false);
                    // The drain witness must report the EFFECTIVE flux under the
                    // same face rule as the kernel, or correct code fails
                    // criterion (b): WallStepDiagnostic uses the RAW cell
                    // conductivity and knows nothing about a factor applied
                    // inside k_face, so scaling the face would move it only
                    // through the temperatures.
                    value.RegisterIntegratedVariable(&value.wall_step_q_eff,
                                                     "wall_step_q_eff", false);
                    // D2z-2. All APPENDED, nothing above moves, and all still
                    // behind side_face_step_height so a keys-off run registers
                    // no new column and witness/d2z1.md5 stays valid.
                    // deliverable 1: the normaliser behind delta_mean, made
                    // auditable. delta_faces / delta_n (the pair count) is the
                    // mean riser height in cells.
                    value.RegisterIntegratedVariable(&value.step_delta_faces,
                                                     "delta_faces", false);
                    // deliverable 3: attribution for SW -- without it a bad
                    // Meier result cannot be pinned on the multi-cell branch.
                    value.RegisterIntegratedVariable(&value.side_P_face_multi,
                                                     "side_P_face_multi", false);
                    // deliverable 2, in the order the packet fixes.
                    value.RegisterIntegratedVariable(&value.wall_step_n_in,
                                                     "wall_step_n_in", false);
                    value.RegisterIntegratedVariable(&value.wall_step_q_in,
                                                     "wall_step_q_in", false);
                    value.RegisterIntegratedVariable(&value.wall_step_q_eff_in,
                                                     "wall_step_q_eff_in", false);
                    value.RegisterIntegratedVariable(&value.step_g_mean,
                                                     "step_g_mean", false);
                    value.RegisterIntegratedVariable(&value.step_g_max,
                                                     "step_g_max", false);
                    value.RegisterIntegratedVariable(&value.wall_step_f_multi,
                                                     "wall_step_f_multi", false);
                    // D2aa deliverables 3 and 4. APPENDED; nothing above moves.
                    // Registered under side_face_step_height, as D2z-2's were:
                    // every D2aa leg sets it and `S0` sets none of it, so the
                    // identity leg and witness/d2z1.md5 stay valid.
                    value.RegisterIntegratedVariable(&value.side_A_bore, "side_A_bore", false);
                    value.RegisterIntegratedVariable(&value.side_A_pads, "side_A_pads", false);
                    value.RegisterIntegratedVariable(&value.side_A_wall, "side_A_wall", false);
                    value.RegisterIntegratedVariable(&value.side_P_bore, "side_P_bore", false);
                    value.RegisterIntegratedVariable(&value.side_P_pads, "side_P_pads", false);
                    value.RegisterIntegratedVariable(&value.side_P_wall, "side_P_wall", false);
                    value.RegisterIntegratedVariable(&value.wall_step_f_multi_pads,
                                                     "wall_step_f_multi_pads", false);
                    // D2ab deliverable 1, the witness. APPENDED; nothing above
                    // moves. Registered only when the explicit pair is GIVEN,
                    // so `S0` (which sets neither) gains no column and the
                    // header-delta identity and witness/d2z1.md5 both hold.
                    // On a constant leg this reads back side_face_T exactly; on
                    // a `live` leg it tracks jet_T_exhaust one step behind.
                    if (value.surface_patch.side_face_h >= 0.0)
                        value.RegisterIntegratedVariable(&value.side_T_eff,
                                                         "side_T_wall", false);
                    // P2-1 deliverables 2 and 6. APPENDED, and registered ONLY
                    // under side_face_h_law = derived, so a `typed` leg gains no
                    // column, the header-delta identity holds and the frozen
                    // tree needs no re-pin. ⚠ This is the "keyed off by default"
                    // half of the rule -- it is NOT also declared as a re-pin.
                    if (value.surface_patch.side_face_h_derived)
                    {
                        value.RegisterIntegratedVariable(&value.side_h_used_min,
                                                         "side_h_used_min", false);
                        value.RegisterIntegratedVariable(&value.side_h_used_p05,
                                                         "side_h_used_p05", false);
                        value.RegisterIntegratedVariable(&value.side_h_used_p50,
                                                         "side_h_used_p50", false);
                        value.RegisterIntegratedVariable(&value.side_h_used_p95,
                                                         "side_h_used_p95", false);
                        value.RegisterIntegratedVariable(&value.side_h_used_max,
                                                         "side_h_used_max", false);
                        value.RegisterIntegratedVariable(&value.side_h_used_n,
                                                         "side_h_used_n", false);
                        // the law's own state, so the branch that decides each
                        // face is an OUTPUT in stated units rather than an
                        // inference from the result (deliverable 2).
                        value.RegisterIntegratedVariable(&value.wall_h_corner,
                                                         "wall_h_corner", false);
                        value.RegisterIntegratedVariable(&value.wall_z_fill,
                                                         "wall_z_fill", false);
                        value.RegisterIntegratedVariable(&value.wall_r_bore,
                                                         "wall_r_bore", false);
                        value.RegisterIntegratedVariable(&value.wall_h_duct,
                                                         "wall_h_duct", false);
                        value.RegisterIntegratedVariable(&value.wall_uDe,
                                                         "wall_uDe", false);
                        value.RegisterIntegratedVariable(&value.wall_sector,
                                                         "wall_sector", false);
                        value.RegisterIntegratedVariable(&value.wall_law_on,
                                                         "wall_law_on", false);
                        // P2-3, APPENDED LAST inside the SAME already-gated
                        // block: no leg without the derived law gains a column,
                        // and every pre-existing column keeps its position.
                        // `wall_datum_r` is the radius the height datum actually
                        // truncated at, and `side_h_used_n_above` is the count
                        // of ABOVE-PLANE faces the law reached -- the witness
                        // P2-3 (a5) needs and which did not exist: `side_h_used_n`
                        // is band-wide and already reads 8968 on the m13 control,
                        // where side_P_face_above is exactly 0.
                        value.RegisterIntegratedVariable(&value.wall_datum_r,
                                                         "wall_datum_r", false);
                        value.RegisterIntegratedVariable(&value.side_h_used_n_above,
                                                         "side_h_used_n_above", false);
                    }
                    // ⚠⚠ P2-4 COLUMNS GO HERE, AFTER `side_h_used_n_above`, and
                    // NOT in the jet_enthalpy block where they were first put.
                    // Registration ORDER is load-bearing for the 2639-file gate
                    // (CODE_CENSUS F1): the jet block is registered BEFORE the
                    // derived-law, property and pcr_* blocks, so appending there
                    // would have shifted every one of those columns on any leg
                    // carrying both. They are gated on `jet_body` /
                    // `jet_wall_mdot = fresh`, both off by default, so no
                    // existing leg gains a column either way -- but "no leg has
                    // it on" is a property of today's inputs, and position is a
                    // property of the code.
                }
            }
            // ⚠⚠ P2-4 COLUMNS: POSITION after the derived-law block, but gated
            // on `jet_body` ALONE and therefore OUTSIDE its braces. The first
            // attempt put them inside, which silently made the body's columns
            // conditional on `side_face_h_law = derived` as well -- the fixture
            // (which does not name the law) then showed `jet_P_body` with ZERO
            // rows and four rows failed on a missing column rather than on the
            // physics. Position and gating are separate choices and were
            // conflated.
            // ⚠ Registration ORDER is load-bearing for the 2639-file gate
            // (CODE_CENSUS F1), which is why the position is kept here rather
            // than moved to the end of Parse.
            if (value.surface_patch.jet_body)
            {
                value.RegisterIntegratedVariable(&value.jet_P_body,
                                                 "jet_P_body", false);
                value.RegisterIntegratedVariable(&value.body_h_used,
                                                 "body_h_used", false);
                value.RegisterIntegratedVariable(&value.body_z_imm,
                                                 "body_z_imm", false);
                value.RegisterIntegratedVariable(&value.body_L,
                                                 "body_L", false);
                value.RegisterIntegratedVariable(&value.body_noLaw_steps,
                                                 "body_noLaw_n", false);
                // the per-REGIME price has been owed since P2-3, which reported
                // ONE price over floor + wall while the wall carried 42 % of the
                // volume. `ledger_E_flake` is a single cumulative total over
                // four += terms; this is its regime-6 (lateral/wall) part,
                // booked at the same site, so floor = total - wall.
                value.RegisterIntegratedVariable(&value.ledger_E_flake_wall,
                                                 "ledger_E_flake_wall", false);
            }
            if (value.surface_patch.jet_wall_mdot
                    == SurfacePatch::WallMdot::Fresh)
                value.RegisterIntegratedVariable(&value.jet_P_retained,
                                                 "jet_P_retained", false);

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

            // D2z-1 Goal 4 (revised): NO firings_per_void column. Under
            // whole-cell it is 1.000 by construction, so it carries no
            // information, and it appended after the eight jet_* columns and
            // broke jet_enthalpy criterion (a). The invariant lives in the
            // n_voided == 1 assertion in Removal.H instead. Registered here is
            // only the conditional Step 20 witness, so the tree audit is
            // absolute: no column appears or disappears anywhere else.
            if (value.spall_surface_normal)
                value.RegisterIntegratedVariable(&value.sn_recession_superseded,
                                                 "sn_recession_superseded", false);
        }

        // D2z: alpha = kappa/(rho*cp) must come FROM THE RUN. Without a
        // material there is no rho*cp to divide by, and the only way to
        // proceed would be a hardcoded diffusivity -- which is the tuning
        // constant the packet forbids. Refuse instead of inventing one.
        if (value.surface_patch.enabled &&
            value.surface_patch.side_face_step_height
            && !value.material && !value.microstructure_enabled)
            Util::Abort(INFO, "surface_patch.side_face_step_height = 1 needs a material: "
                              "the conductance correction uses alpha = kappa/(rho*cp) "
                              "taken from mat->effective_rhocp(T), and hardcoding a "
                              "diffusivity would reintroduce exactly the tuning constant "
                              "this closure exists to remove");

        // D2q-2b item 5: per-cell removal preconditions, cross-checked here
        // because they live in three different parse blocks.
        if (value.spall_per_cell_removal)
        {
            if (!value.removal_enabled())
                Util::Abort(INFO, "spall.per_cell_removal = 1 requires removal "
                                  "(spall.enabled or vapor.enabled)");
            if (!(value.surface_patch.enabled && value.surface_patch.side_face_flux))
                Util::Abort(INFO, "spall.per_cell_removal = 1 requires "
                                  "surface_patch.side_face_flux = 1: without it an exposed "
                                  "lateral face is never heated, so the criterion evaluated "
                                  "there reads bulk rock and the mechanism is silently inert");
            // side_face_flux already implies this, but state it directly: the
            // requirement belongs to per-cell removal on its own terms (phi is
            // still the sole input to UpdateSurfaceMaskFromPhi, and a cell
            // voided off-top does not move its column's phi). See RESULTS.md Q1.
            if (!value.surface_follow_mask)
                Util::Abort(INFO, "spall.per_cell_removal = 1 requires "
                                  "surface.follow_mask = 1: the phi-window surface mask "
                                  "would keep selecting a cell this mechanism has voided");
            // The ONE case where the column gate and the per-cell gate would
            // genuinely differ. Under weibull the column path gates on the
            // 4-connected cluster labeller, which reduces to Sp_top >= 1 only
            // while the area filter is inactive. With the filter ON there is no
            // lateral analogue of a top-z-plane labeller, so refuse rather than
            // silently fall back to Sp >= 1. Generalising the labeller to the
            // exposed-face set is a separate packet.
            if (value.weibull_enabled &&
                value.weibull_detach_mode == DetachMode::ConnectedCluster &&
                value.weibull_A_crit > 0.0)
                Util::Abort(INFO, "spall.per_cell_removal = 1 is incompatible with "
                                  "weibull.detach_mode = connected_cluster and "
                                  "weibull.A_crit > 0: the 4-connected cluster filter is "
                                  "defined over the top-z plane and has no lateral "
                                  "analogue, so the column gate and the per-cell gate "
                                  "would differ. Unset A_crit (the filter is then a no-op "
                                  "and both gates reduce to Sp >= 1), or drop "
                                  "per_cell_removal");
            // Item 8: overhang census. Registered LAST and only with the key
            // on, so every pre-D2q-2b thermo.dat column keeps its position.
            value.RegisterIntegratedVariable(&value.pcr_overhang_cells, "pcr_overhang_cells", false);
            value.RegisterIntegratedVariable(&value.pcr_lateral_cells,  "pcr_lateral_cells",  false);

        }

        // P2-2e: the constraint-depth witnesses. ⚠ Registered LAST and ONLY
        // under d > 0, so a d = 0 leg gains no column and every existing column
        // keeps its position -- the "keyed off by default" half of the rule,
        // NOT also declared as a re-pin (CODE_CENSUS F1: order is load-bearing).
        if (value.sp_cref_d > 0.0)
        {
            value.RegisterIntegratedVariable(&value.cref_d_out,   "cref_d",       false);
            value.RegisterIntegratedVariable(&value.cref_n_vert,  "cref_n_vert",  false);
            value.RegisterIntegratedVariable(&value.cref_n_lat,   "cref_n_lat",   false);
            value.RegisterIntegratedVariable(&value.cref_n_clamp, "cref_n_clamp", false);
            value.RegisterIntegratedVariable(&value.cref_n_neg,   "cref_n_neg",   false);
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
        // ⚠⚠ P2-2a change A: THE ONLY WRITE OF THE STRESS-FREE REFERENCE.
        // Immediately after HeatConduction::Initialize has applied the deck's
        // initial condition and before anything advances, so this is the
        // temperature at which the rock is unstressed BY DEFINITION. Never
        // written again -- not in Advance, not in UpdatePropertyFields, not in
        // Regrid (where RegisterNewFab interpolates it from itself). (b9)
        // asserts the single site; (b2) asserts it does not drift.
        StoreInitialTemperature(lev);
        Mechanics<Model::Solid::Affine::Isotropic>::Initialize(lev);
        if (removal_enabled()) InitPhi(lev);
        if (removal_enabled()) ApplySeedVoid(lev);
        InitSurfaceMask(lev);
        // InitializeMicrostructure must run before ConvertTtoH on the
        // microstructure path because the H(T) inversion is per-phase.
        if (microstructure_enabled) InitializeMicrostructure(lev, /*announce=*/true,
                                                             /*reset_damage=*/true);
        ConvertTtoH(lev);
    }

    // P2-2a change A: copy the initial temperature into the stress-free
    // reference. ⚠ The ONE write site. Deliberately a plain copy of temp_mf
    // rather than a re-evaluation of the IC parser: the reference must be the
    // temperature the solver actually starts from, including whatever clamping
    // or BC fill the IC path applied, or the bracket would not be exactly zero
    // on a column that started uniform.
    void StoreInitialTemperature(int lev)
    {
        const std::size_t l = static_cast<std::size_t>(lev);
        if (T_init_mf.size() <= l || temp_mf.size() <= l) return;
        if (!T_init_mf[lev] || !temp_mf[lev]) return;
        const int ng = std::min(T_init_mf[lev]->nGrow(), temp_mf[lev]->nGrow());
        amrex::MultiFab::Copy(*T_init_mf[lev], *temp_mf[lev], 0, 0, 1, ng);
        Util::RealFillBoundary(*T_init_mf[lev], geom[lev]);
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        HeatConduction::TimeStepBegin(a_time, a_step);
        // Step 7c/7d repair for the mechanics AMR path. With mechanics off the
        // thermal/removal AMR path keeps its fields across steps (regrid
        // repair happens in Regrid), so skip the per-step re-initialisation.
        if (microstructure_enabled && max_level > 0 &&
            m_type != Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);
    }

    void TimeStepComplete(Set::Scalar a_time, int a_step) override
    {
        HeatConduction::TimeStepComplete(a_time, a_step);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepComplete(a_time, a_step);
        (void)a_time;
        (void)a_step;
        if (microstructure_enabled && max_level > 0 &&
            m_type != Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::Regrid(lev, time);
        // AMR: a remade level keeps exact phi/removed where it already
        // existed and gets conservatively interpolated values where it is
        // new; the lateral slope limiter smears the per-column height at
        // steps. Rebuild the column-linear phi from each column's top cell.
        if (removal_enabled() && lev > 0) RepairRemovalStateAfterRegrid(lev);
        InitSurfaceMask(lev);
        if (microstructure_enabled)
        {
            InitializeMicrostructure(lev, /*announce=*/false,
                                     /*reset_damage=*/!damage_enabled);
            // H and T were interpolated independently; re-derive T and the
            // phase fractions from the conserved H on the remade level.
            if (lev > 0) InvertHtoT(lev);
        }
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
    // P0b: declared here, ahead of its first use as a parameter type in
    // SharedSurfaceRecon. An in-class member function BODY is a
    // complete-class context, but its DECLARATION is not.
    struct SurfaceRecon
    {
        Set::Scalar T_surf  = 0.0;   // reconstructed surface temperature [K]
        Set::Scalar ell_eff = 0.0;   // recovered decay length [m]
        // P1g-3: the FLUX-based heated depth actually used, 0 when the shipped
        // two-cell length was used instead. Logged so the event row records the
        // length the thickness was built from, not merely the one available.
        Set::Scalar ell_q   = 0.0;
        Set::Scalar h_f     = 0.0;   // flake thickness min(ell_eff, dz) [m]
        Set::Scalar T_bar   = 0.0;   // flake-average temperature [K]
        bool ell_clamped    = false; // the ell_eff >= dz/4 floor bound
        bool h_f_capped     = false; // h_f hit the spall.h_spall_max SEARCH depth
        bool h_f_bounded    = false; // h_f hit spall.flake_h_max (deliverable 1)
        bool T_clamped      = false; // the gas-temperature ceiling bound
        bool flat           = false; // no usable gradient; T_surf == T_top
    };

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        // P0b deliverable 1.1: ONE surface temperature per column for this
        // whole step -- criterion, BC and removal. Must be first: every
        // consumer below reads it, and temp_mf is swapped into temp_old_mf
        // partway through.
        BuildSurfaceTColumns(lev, dt);
        // P2-1 deliverable 2: the wall law's per-step state and the ONE height
        // datum, built here because the criterion -- the earliest of the four
        // reader roles -- runs below in UpdateRemovalAfterCohesive. Building it
        // once is what makes the four roles agree per face; building it per
        // role would make them disagree by a sub-step. No-op, and no
        // allocation, under the default side_face_h_law = typed.
        UpdateWallLawState(lev);
        // P2-2 change 1: the property fields, refreshed at the cell
        // temperature BEFORE the criterion and the conduction dispatch read
        // them. No-op unless a phase sets law = grimsel.
        UpdatePropertyFields(lev);
        Mechanics<Model::Solid::Affine::Isotropic>::Advance(lev, time, dt);
        UpdateDamageAfterMechanics(lev, time, dt);
        UpdateCohesiveAfterMechanics(lev, time, dt);
        // Step 12: Sp = K_I/K_Ic(T) surface diagnostic for the sp_weibull
        // branch. No-op under damage_law. Reads the same nodal stress field
        // produced by the static elastic solve (Step 7) used by the
        // hu_breakage_index f_b pipeline.
        // ⚠⚠ P2-2e: the constraint reference is gathered HERE, immediately
        // before the FIRST of the two firing sites. UpdateSpAfterMechanics runs
        // at this line and UpdateRemovalAfterCohesive (the depth scan) below it,
        // and only UpdateSpClusters -- which writes neither the temperature nor
        // the void mask -- runs between them, so one gather here IS the state
        // both sites see. Removal.H's own gathers are function-local and built
        // AFTER this call, which is why the reference cannot live there. No-op,
        // no allocation and no reduction at d = 0.
        BuildConstraintRefColumns(lev);
        UpdateSpAfterMechanics(lev, time, dt);
        // Step 15c: 4-connected component labeller over the Sp >= 1 surface
        // zone. Default-off everywhere except sp_weibull + weibull.enabled;
        // writes Sp_cluster_id_mf and per-step cluster diagnostics CSV.
        UpdateSpClusters(lev, time, dt);
        // Step 10/11: removal (spall and/or vaporisation) must happen
        // BEFORE the H/Temp swaps below so the reset H(T_amb)/T_amb values
        // become the "old" state consumed by the thermal advance step.
        UpdateRemovalAfterCohesive(lev, time, dt);
        // D2q-2b item 6: lateral per-cell removal, AFTER the column pass so the
        // top-cell path is untouched and this step's newly exposed faces are
        // the ones scored next step. No-op with the key off.
        PerCellRemoval(lev, time, dt);
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
        // Perf (2026-09-18): phase, grain_id and kappa_phase are written only
        // by InitializeMicrostructure (Initialize / Regrid), so their ghost
        // layers are filled once after each (re)initialisation instead of
        // every step. D is refilled every step (removal and damage write it).
        if (lev >= static_cast<int>(micro_ghost_filled.size()))
            micro_ghost_filled.resize(lev + 1, 0);
        if (!micro_ghost_filled[lev])
        {
            Util::RealFillBoundary(*phase_mf[lev],     geom[lev]);
            Util::RealFillBoundary(*grain_id_mf[lev],  geom[lev]);
            Util::RealFillBoundary(*kappa_mf[lev],     geom[lev]);
            micro_ghost_filled[lev] = 1;
        }
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
        // AMR: only the surface-owning level applies the beam, the surface
        // patch, the side faces and the losses (coarser levels: conduction
        // only, their covered cells are overwritten by average-down).
        const bool        owner = SurfaceOwner(lev);
        const Set::Scalar P_now = owner ? beam.power_at(time) : Set::Scalar(0.0);
        const bool follow_mask = surface_follow_mask;
        const amrex::Vector<int> k_top_col =
            (void_aware && (P_now > 0.0 || follow_mask || ledger_enabled))
            ? ColumnTopSolid(lev) : amrex::Vector<int>();
        const int* ktc = k_top_col.data();
        const int ncx = dhi[0] - dlo[0] + 1;
        // D2z: advance the per-column recession phase off this same mask,
        // before anything reads it this step. Cheap (one pass over columns) and
        // rank-consistent, because k_top_col is already domain-reduced.
        if (surface_patch.side_face_step_height && owner && !k_top_col.empty())
        {
            step_split_nx = static_cast<std::size_t>(ncx);
            UpdateColumnPhase(k_top_col, time);
            AssertStepSplit(k_top_col, time);
        }
        const Set::Scalar* c_tl = col_t_last.empty() ? nullptr : col_t_last.data();
        const Set::Scalar* c_dt = col_dt.empty()     ? nullptr : col_dt.data();
        const int*         c_nv = col_nvoid.empty()  ? nullptr : col_nvoid.data();
        // D2q-1: the beam path counts SOLID cells, not geometric distance.
        // Only needed when the beam is actually on (BeamSource returns 0 at
        // P = 0), which keeps the large beam-free runs free of the gather.
        const amrex::Vector<int> solid_above_col =
            (void_aware && P_now > 0.0 && !k_top_col.empty())
            ? ColumnSolidAboveCount(lev, k_top_col) : amrex::Vector<int>();
        const int* sab    = solid_above_col.data();
        const int  nz_dom = dhi[2] - dlo[2] + 1;
        const int  klo_d  = dlo[2];

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
        // Perf (2026-09-18): exp(-a_D*D) == 1 exactly when a_D == 0, so skip
        // the 7 exp() calls per cell (12 % of step time) on undamaged runs.
        const bool        a_D_on  = (a_D != 0.0);
        const Set::Scalar h_gb0   = h_gb_conductance;
        const Set::Scalar dx_ax[3] = {DX[0], DX[1], DX[2]};
        // D2k: side-face flux. One SurfaceCellFlux call serves both lateral
        // axes and the face-form conductance 2k/dx equals the top face's
        // 2k/dz, so require cubic cells rather than silently using dz for a
        // face whose normal spacing is dx.
        const bool        side_on = surface_patch.enabled && surface_patch.side_face_flux && owner;
        const Set::Scalar side_f  = surface_patch.side_face_factor;
        const Set::Scalar side_f_above = surface_patch.side_face_factor_above;
        // D2x: both new keys default to the previous behaviour, so an unset run
        // takes the identical branch and the identical expression below.
        const Set::Scalar side_r_min2 = surface_patch.side_face_r_min
                                      * surface_patch.side_face_r_min;
        const bool side_area_exact = surface_patch.side_face_area_exact;
        // D2z: the step-height closure's kernel locals.
        //
        // alpha = kappa / (rho*cp) is the thermal DIFFUSIVITY, taken from where
        // both halves are consumed: the per-cell conductivity field kp (== the
        // stencil's own kc) and mat->effective_rhocp(T), evaluated at the SAME
        // cell temperature the conduction stencil uses (T_old(i,j,k)).
        // ⚠ alpha is NOT the codebase's `kappa`, which is CONDUCTIVITY.
        // dz*v/alpha is dimensionless; dz*v/conductivity is not, and the error
        // would be silent -- no abort, no failed assertion.
        // ⚠ The global `material` is NULL under microstructure -- which is
        // every input this closure runs on. Use the CELL'S PHASE material, the
        // same object the phase tables are built from, and fall back to the
        // global only when there is no microstructure.
        const Numeric::Material::Constant* mat_ph =
            phase_materials.empty() ? nullptr : phase_materials.data();
        const int mat_nph = static_cast<int>(phase_materials.size());
        const Numeric::Material::Material* mat_k = material;
        const bool side_step_h = surface_patch.side_face_step_height;
        // P1e: the staircase orientation correction. Default off, and the
        // parser refuses it without side_face_step_height, because all four
        // AreaOfDelta sites live behind that key -- a silent no-op is the
        // failure mode the guard prevents.
        const bool side_staircase = surface_patch.side_face_staircase;
        const bool wetted_law  =
            (surface_patch.side_face_area_law == SurfacePatch::SideAreaLaw::Wetted);
        const Set::Scalar t_now = time;
        // D2q-2b item 7: cap the coefficient a SIDE face sees, ONLY when
        // per-cell removal is on. -1 = no cap, which is every pre-D2q-2b run.
        // Applied identically in the kernel and in the ledger, or the ledger
        // stops reproducing the kernel's flux.
        const Set::Scalar side_hmax =
            (surface_patch.side_face_flux && spall_per_cell_removal)
            ? surface_patch.side_face_h_max : Set::Scalar(-1.0);
        // D2aa: the explicit lateral (h, T) pair, hoisted for the device lambda.
        const Set::Scalar side_h_ex = surface_patch.side_face_h;
        // D2ab deliverable 2: the annulus pair reaches only faces beyond the
        // burner radius; inside it a face keeps the inherited impinging map.
        // 0 (the default) squares to 0 and r2_col >= 0 always holds, so an
        // unset run takes the identical branch -- identity by construction.
        const Set::Scalar side_h_rmin2 = surface_patch.side_face_h_r_min
                                       * surface_patch.side_face_h_r_min;
        const Set::Scalar inv_dx  = 1.0 / DX[0];
        const Set::Scalar inv_dy  = 1.0 / DX[1];
        if (side_on && !(std::abs(DX[0] - DX[1]) <= 1.0e-12 * DX[0] &&
                         std::abs(DX[0] - DX[2]) <= 1.0e-12 * DX[0]))
            Util::Abort(INFO, "surface_patch.side_face_flux = 1 requires cubic cells "
                        "(dx = dy = dz); got ", DX[0], ", ", DX[1], ", ", DX[2]);

        // Step 11b: the Hu onset validation uses a homogeneous one-phase
        // microstructure so the DP damage fields exist. Apply the same
        // surface patch here as in AdvanceMaterial(): prescribed_T injects the
        // Dirichlet face flux; convective_flame (Step 15a) injects the Robin
        // flame-jet flux q = h_fl*(T_flame - T_cell).
        const bool        sp_on    = surface_patch.enabled && owner;
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
        // D2ab deliverable 1: the RESOLVED wall temperature, not
        // surface_patch.side_face_T -- under `live` the key holds the -1
        // sentinel and the number lives in side_T_eff.
        // ⚠ Hoisted HERE, not with the other side_* locals above, because
        // BuildFlameColumns is what resolves it. Reading it earlier would hand
        // the kernel the PREVIOUS step's value while the march and the ledger
        // used this one -- a split the power identity would catch, but only
        // after a 16-minute leg.
        // P2-3 site 1 of 5 (KERNEL): the per-face value is resolved inside the
        // loop by SideFaceTgas; this remains the fallback the resolver returns
        // for a below-plane face or with the march off.
        const Set::Scalar side_T_ex = side_T_eff;
        const Set::Scalar* fh  = sp_on ? flame_h.data()  : nullptr;
        const Set::Scalar* fTg = sp_on ? flame_Tg.data() : nullptr;
        // D2r-1: per-column stand-off, so the side factor can differ above and
        // below the nozzle plane. Same array the march and the pin rule read.
        const Set::Scalar* fs  = sp_on ? flame_s.data()  : nullptr;
        // C1: pinned runs the face Newton (face_form = true) plus the pin.
        const bool        sp_pinned = sp_on && surface_patch.robin_form ==
                                      SurfacePatch::RobinForm::Pinned;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face || sp_pinned;
        // P0b rev5 deliverable 7: the lateral rule fires only under `surface`,
        // and rev7 deliverable 6 lets it be switched off there without
        // touching the top faces.
        const bool        sp_surfform = surface_patch.robin_form ==
                                        SurfacePatch::RobinForm::Surface
                                     && surface_patch.side_T_form ==
                                        SurfacePatch::SideTForm::Follow;
        if (sp_pinned) BuildPinEffective(lev, time);
        const Set::Scalar* tpe = sp_pinned ? pin_T_eff.data() : nullptr;

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
                Set::Scalar kc = a_D_on ? kp(i, j, k) * std::exp(-a_D * Dc) : kp(i, j, k);
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
                    Set::Scalar k_n = a_D_on ? kp(ii, jj, kk) * std::exp(-a_D * Df(ii, jj, kk))
                                             : kp(ii, jj, kk);
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
                    Set::Scalar kf = 1.0 / inv;
                    // D2z Step 3: correct the LATERAL face conductance across a
                    // staircase step. THIS is the site -- the harmonic face
                    // conductance -- not the level-set shift in Removal.H,
                    // which is recession, not conduction.
                    //
                    // The interface at the SHORTER column's top level is solid
                    // on both sides: no void, so the void test above leaves it
                    // alone. The model treats the taller column's cell there as
                    // buried interior when it physically sits only delta*dz
                    // below the surface, and that is the drain -- k*dT/dx into
                    // rock nothing has heated.
                    //
                    // g(delta) scales it to what a cell at its TRUE depth would
                    // pass. ⚠ SYMMETRIC: the test and delta depend only on the
                    // unordered pair, so the face carries the same conductance
                    // evaluated from either side and energy is conserved.
                    if (side_step_h && axis < 2 && dk == 0 && ktc && c_tl
                        && (mat_k || (mat_ph && pc >= 0 && pc < mat_nph)))
                    {
                        const int c1 = (j  - dlo[1]) * ncx + (i  - dlo[0]);
                        const int c2 = (jj - dlo[1]) * ncx + (ii - dlo[0]);
                        const int k1 = ktc[c1], k2 = ktc[c2];
                        const int kS = k1 < k2 ? k1 : k2;      // shorter column's top
                        const int kT = k1 < k2 ? k2 : k1;
                        const int cS = k1 < k2 ? c1 : c2;
                        const int cT = k1 < k2 ? c2 : c1;
                        if (kT > kS && k == kS)
                        {
                            int fb = 0;
                            const Set::Scalar d =
                                PairDelta(ktc, c_tl, c_dt, c_nv, cT, cS, t_now, fb);
                            // v = dz / dt_i of the TALLER column -- it is that
                            // column's cell being treated as buried, so its own
                            // near-surface length alpha/v is the relevant one.
                            // Symmetric, because "taller" is a property of the
                            // unordered pair.
                            const Set::Scalar rc =
                                (mat_ph && pc >= 0 && pc < mat_nph)
                                ? mat_ph[pc].effective_rhocp(Tc)
                                : mat_k->effective_rhocp(Tc);
                            Set::Scalar a_diff = 0.0;
                            if (rc > 0.0 && c_dt[cT] > 0.0 && kc > 0.0)
                            {
                                const Set::Scalar alpha_d = kc / rc;
                                const Set::Scalar v = DX[2] / c_dt[cT];
                                if (alpha_d > 0.0) a_diff = DX[2] * v / alpha_d;
                            }
                            kf *= StepConductance(d, kT - kS, a_diff);
                        }
                    }
                    return kf;
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
                if (!self_void && owner)
                    SurfaceCellFlux(in_patch, sp_face, surf_c,
                                    fh  ? fh[col_c]  : Set::Scalar(0.0),
                                    fTg ? fTg[col_c] : Set::Scalar(0.0),
                                    kc, DX[2], Tc,
                                    eps_lo, eps_hi, T_eps_b, sig_SB, h_c, T_a,
                                    q_in, q_out, q_rob, q_los,
                                    tpe ? tpe[col_c] : Set::Scalar(0.0),
                                    nullptr,
                                    SurfaceRobinT(
                                        static_cast<std::size_t>(col_c)));
                if (in_patch) div_kgrad += q_in * inv_dz;

                // D2k: flux on EXPOSED VERTICAL FACES of a staircase wall.
                // A solid cell whose column lies in the patch and which has a
                // void lateral neighbour is a piece of exposed wall: the gas
                // reaches it, and before D2k the model heated only the column
                // top, so a firing cell beside a taller cold column lost
                // k(T_top - T_wall)/dx through a face that was never heated
                // (~240 kW/m2 at 2 mm, ~480 at 1 mm, against q_pin ~ 520 --
                // the 1/dx scaling that makes refinement worse, D2j-0b).
                //
                // The face uses the SAME SurfaceCellFlux call as the column's
                // top face -- the column's h and T_gas, its pin, the same
                // losses -- evaluated at this cell's own T. Reusing the shared
                // routine is load-bearing: it carries the state rule
                // T_s = min(T_f, T_pin), so a side-heated cell cannot be
                // driven above the firing temperature and become an
                // unremovable heat reservoir.
                //
                // The `s <= 0` exclusion applies identically and for free: a
                // column above the nozzle plane has flame_h = 0 (set in
                // BuildFlameColumns under jet_enthalpy, or by the h_expr
                // itself), so q_side is exactly 0 there. Nothing here weakens
                // the D2j-0b C1/C2/C3 finding.
                //
                // Out-of-domain lateral neighbours are NOT faces: removed_mf's
                // ghosts are undefined outside the domain (the k_face lambda
                // above takes the same care).
                const Set::Scalar r2_col = (x - sp_x0) * (x - sp_x0)
                                         + (y - sp_y0) * (y - sp_y0);
                const bool col_in_patch = side_on && sp_on && !sp_presc && !self_void
                                        && r2_col <= sp_r2 && r2_col >= side_r_min2;
                if (col_in_patch)
                {
                    const bool vxm = (i > dlo[0]) && rem(i-1, j, k) > 0.5;
                    const bool vxp = (i < dhi[0]) && rem(i+1, j, k) > 0.5;
                    const bool vym = (j > dlo[1]) && rem(i, j-1, k) > 0.5;
                    const bool vyp = (j < dhi[1]) && rem(i, j+1, k) > 0.5;
                    const int nfx = (vxm ? 1 : 0) + (vxp ? 1 : 0);
                    const int nfy = (vym ? 1 : 0) + (vyp ? 1 : 0);
                    if (nfx + nfy > 0)
                    {
                        Set::Scalar s_in = 0.0, s_out = 0.0, s_rob = 0.0, s_los = 0.0;
                        Set::Scalar h_side = fh ? fh[col_c] : Set::Scalar(0.0);
                        Set::Scalar T_side = fTg ? fTg[col_c] : Set::Scalar(0.0);
                        // D2aa: the explicit lateral pair overrides the inherited
                        // top-face map. Both or neither, enforced at parse.
                        // D2ab deliverable 2: gated on the COLUMN radius,
                        // which is the same r2_col the patch test above uses.
                        if (side_h_ex >= 0.0 && r2_col >= side_h_rmin2)
                        { h_side = side_h_ex;
                          // P2-3 site 1 of 5 (KERNEL -- what actually heats the
                          // rock): the gas temperature at THIS face's height.
                          T_side = SideFaceTgas(static_cast<int>(col_c), k); }
                        // P2-1 deliverable 4, ROLE 1 of 4 (the KERNEL -- what
                        // actually heats the rock). One resolver, per face.
                        Set::Scalar h_face =
                            SideFaceH(h_side, r2_col, static_cast<int>(col_c), k);
                        if (side_hmax > 0.0 && h_face > side_hmax) h_face = side_hmax;
                        // P0b rev5 deliverable 7 -- ONE RULE, TWO REGIMES for a
                        // lateral face. Revision 2 left every side face on the
                        // cell average, and the Meier rerun then measured the
                        // side faces carrying 47 % of P_robin over the rate
                        // window (62 % by the end of the run), so "a lateral
                        // face has no column profile" stopped being a small
                        // approximation.
                        //   7a  HOLE WALL, r >= side_face_h_r_min: the FACE
                        //       Newton. A non-receding wall conducts inward
                        //       diffusively -- penetration passes dz at 1.45 s
                        //       (2 mm) and is 8.3 cells by 100 s -- so its
                        //       near-face profile is very nearly LINEAR, which
                        //       is the half-cell profile the face form already
                        //       assumes. Not the surface reconstruction: that
                        //       fits an EXPONENTIAL of the column's recession
                        //       rate, and a wall is not receding.
                        //   7b  FLOOR RISER, r < side_face_h_r_min: the owning
                        //       COLUMN's T_surf. A riser is the vertical step
                        //       of a receding floor, so the column's own
                        //       exponential IS its profile -- exact for a
                        //       one-cell step, and the multi-cell riser share
                        //       is logged to bound the rest.
                        // ⚠ Both regimes are reached ONLY under robin_form =
                        // surface. Under cell/face/pinned lat_face == sp_face
                        // and lat_Ts == 0, so the controls and every frozen
                        // study from d2v to d2ad stay byte-identical.
                        const bool lat_wall = sp_surfform
                                            && (r2_col >= side_h_rmin2);
                        const bool lat_face = sp_surfform ? true : sp_face;
                        const Set::Scalar lat_Ts = (sp_surfform && !lat_wall)
                            ? SurfaceRobinT(static_cast<std::size_t>(col_c))
                            : Set::Scalar(0.0);
                        SurfaceCellFlux(true, lat_face, 1.0,
                                        h_face,
                                        T_side,
                                        kc, DX[0], Tc,
                                        eps_lo, eps_hi, T_eps_b, sig_SB, h_c, T_a,
                                        s_in, s_out, s_rob, s_los,
                                        tpe ? tpe[col_c] : Set::Scalar(0.0),
                                        nullptr, lat_Ts);
                        // dx == dy == dz is enforced when the key is on, so one
                        // call serves both axes and inv_dx == inv_dy.
                        // D2r-1: above the nozzle plane (s <= 0) the side
                        // factor may differ. Sentinel < 0 keeps the single
                        // uniform factor, so unset is byte-identical.
                        const bool above_c =
                            (side_f_above >= 0.0 && fs && !(fs[col_c] > 0.0));
                        const Set::Scalar sf_c = above_c ? side_f_above : side_f;
                        // D2x: the derived area applies to BELOW-plane faces
                        // only. Above the plane the D2r-1 sentinel path is
                        // untouched, so a leg setting side_face_factor_above =
                        // 0.0 still puts exactly 0 W there; honest area on the
                        // free surface would blanket it (the D2u fault).
                        // D2z: the step-height closure. delta uses BOTH
                        // columns' phases; the per-cell split is Step 3b.
                        // BELOW-plane only, same as D2x -- above the plane the
                        // D2r-1 sentinel path is untouched.
                        if (side_step_h && !above_c && ktc && c_tl)
                        {
                            const int kT = ktc[col_c];
                            Set::Scalar wsum = 0.0;
                            int fb = 0;
                            const Set::Scalar pT =
                                ColPhi(c_tl, c_dt, c_nv, col_c, t_now, fb);
                            const int nbr[4] = {col_c - 1, col_c + 1,
                                                col_c - ncx, col_c + ncx};
                            const bool has[4] = {vxm, vxp, vym, vyp};
                            // P1e: riser cell counts per face, 0 where the
                            // face is not open. SITE 1 OF 4 -- the flux kernel,
                            // i.e. what actually heats the rock.
                            int rn[4] = {0, 0, 0, 0};
                            for (int d = 0; d < 4; ++d)
                            {
                                if (!has[d]) continue;
                                const int cS = nbr[d];
                                const int kS = ktc[cS];
                                const int n  = kT - kS;
                                if (n < 1) continue;
                                rn[d] = n;
                                int fbS = 0;
                                const Set::Scalar pS =
                                    ColPhi(c_tl, c_dt, c_nv, cS, t_now, fbS);
                                const Set::Scalar w =
                                    StepCellWeight(n, k, kT, kS, pT, pS);
                                wsum += AreaOfDelta(w < 0.0 ? 0.0 : w, wetted_law);
                            }
                            // P1e: the factor is PER COLUMN, so it multiplies the
                            // summed area, never an individual face's share.
                            if (side_staircase)
                                wsum *= StaircaseFactor(rn[0], rn[1], rn[2], rn[3]);
                            div_kgrad += (s_in - s_out) * wsum * inv_dx;
                        }
                        else if (side_area_exact && !above_c && ktc)
                        {
                            const int nc = ktc[col_c];
                            Set::Scalar wsum = 0.0;
                            if (vxm) wsum += RiserWeight(nc - ktc[col_c - 1],   DX[0], DX[2]);
                            if (vxp) wsum += RiserWeight(nc - ktc[col_c + 1],   DX[0], DX[2]);
                            if (vym) wsum += RiserWeight(nc - ktc[col_c - ncx], DX[0], DX[2]);
                            if (vyp) wsum += RiserWeight(nc - ktc[col_c + ncx], DX[0], DX[2]);
                            div_kgrad += (s_in - s_out) * wsum * inv_dx;
                        }
                        else
                        div_kgrad += sf_c * (s_in - s_out)
                                   * (Set::Scalar(nfx) * inv_dx + Set::Scalar(nfy) * inv_dy);
                    }
                }

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
                                                       // D2q-1: `sab` is the solid-cell count above k;
                                                       // it is null when the beam is off, where the
                                                       // path is unused and k_top - k is kept.
                                                       (void_aware && ktc)
                                                         ? SolidFacePath(
                                                             sab ? sab[static_cast<std::size_t>(col_c) * nz_dom
                                                                       + (k - klo_d)]
                                                                 : ktc[col_c] - k, DX[2])
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
        if (owner && (sp_pinned || (sp_on && surface_patch.use_expr)))
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
        // D2ab deliverable 1: resolve the wall gas temperature ONCE per step,
        // here, before the march (called below), the conduction kernel and the
        // ledger read it. All three must see the SAME number: the march charges
        // the side term into jet_P_face, so a kernel/march split breaks the
        // power identity -- the one check that catches it.
        //
        // ⚠ Under `live` this is the march's own jet_T_exhaust from the
        // PREVIOUS step. The march that computes this step's value runs below,
        // so using its own output inside itself would be circular. Same
        // one-step lag as jet_T_rec (D2b, :2793) and at dt = 16 ms the exhaust
        // moves well under 1 K per step.
        //
        // ⚠ First step: no exhaust exists yet, and the ONE value that must not
        // be used as the seed is the 293 K ambient -- an h chosen to heat the
        // wall would REFRIGERATE it. Seed with the nozzle temperature, the
        // gas's upper bound, for exactly one step out of ~31 000.
        if (surface_patch.side_face_T_live)
            side_T_eff = (jet_T_exhaust > 0.0) ? jet_T_exhaust
                                               : surface_patch.jet_T_nozzle;
        else
            side_T_eff = surface_patch.side_face_T;
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
        const Set::Scalar pad_open = surface_patch.pad_shield_open;
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
                // D2v: jet_above_h evaluates the same expression above the
                // plane. The map must be radially confined there by its own
                // `s` branch -- a scalar would blanket the 89 % of columns that
                // sit above the plane late in a feet run (D2o-0's 1400 K
                // blanket, and the fault that withdrew D2u).
                const Set::Scalar hj = (s > 0.0 || surface_patch.jet_above_h)
                                     ? surface_patch.h_f(x, y, r, s, time)
                                     : Set::Scalar(0.0);
                if (!(std::isfinite(hj) && hj >= 0.0))
                    Util::Abort(INFO, "surface_patch.h_expr must be finite and >= 0 in the "
                                "patch; got h = ", hj, " at column (", i + domain.smallEnd(0),
                                ", ", j + domain.smallEnd(1), "), r = ", r, ", s = ", s,
                                ", time = ", time);
                flame_h[col]  = hj;
                // D2aa deliverable 1: PAD-RING SHIELDING. The pads REST ON THE
                // ROCK (Meier drawing: Ø56 bore, three exhaust slots, three arc
                // pads on r = 28-40 mm under 12 mm of steel). Under them there
                // is essentially no gas, so the floor there is ADIABATIC except
                // where the slots open. Applied HERE, in BuildFlameColumns, so
                // the kernel, the ledger AND the enthalpy march all inherit one
                // value -- shielding the kernel alone would break the power
                // identity, criterion (c).
                // ⚠ NOT a sink: "the steel is a sink" needs a contact
                // conductance and a body temperature, neither of which is known
                // (risk 3). Adiabatic-except-at-the-slots is derivable; a sink
                // would be invented.
                // ⚠ pad_shield_open is the ring's OPEN FRACTION and it is a
                // BOUND, not a value: the slot width is on record nowhere
                // (`meier_thesis_notes.md:40` says only "three exhaust slots"),
                // so the leg is run at the ends of the bracket and reported as
                // bounded. Sentinel < 0 = inert = exactly today's behaviour.
                if (pad_open >= 0.0 && r >= surface_patch.foot_r_inner
                                    && r <  surface_patch.foot_r_outer)
                    flame_h[col] *= pad_open;
                // D2v: above the plane the march never visits this column
                // (JetEnthalpyMarch skips !(flame_s > 0) at both its column
                // pass and its bin pass), so this value is final -- and the
                // heat it delivers is NOT deducted from the jet's enthalpy
                // budget. Sentinel < 0 keeps the D2a2 behaviour exactly.
                flame_Tg[col] = (s > 0.0 || !(surface_patch.jet_T_above >= 0.0))
                              ? surface_patch.jet_T_ent : surface_patch.jet_T_above;
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
            // D2aa deliverable 1: PAD-RING SHIELDING. The pads REST ON THE
            // ROCK (Meier drawing: Ø56 bore, three exhaust slots, three arc
            // pads on r = 28-40 mm under 12 mm of steel). Under them there
            // is essentially no gas, so the floor there is ADIABATIC except
            // where the slots open. Applied HERE, in BuildFlameColumns, so
            // the kernel, the ledger AND the enthalpy march all inherit one
            // value -- shielding the kernel alone would break the power
            // identity, criterion (c).
            // ⚠ NOT a sink: "the steel is a sink" needs a contact
            // conductance and a body temperature, neither of which is known
            // (risk 3). Adiabatic-except-at-the-slots is derivable; a sink
            // would be invented.
            // ⚠ pad_shield_open is the ring's OPEN FRACTION and it is a
            // BOUND, not a value: the slot width is on record nowhere
            // (`meier_thesis_notes.md:40` says only "three exhaust slots"),
            // so the leg is run at the ends of the bracket and reported as
            // bounded. Sentinel < 0 = inert = exactly today's behaviour.
            if (pad_open >= 0.0 && r >= surface_patch.foot_r_inner
                            && r <  surface_patch.foot_r_outer)
                flame_h[col] *= pad_open;
            flame_Tg[col] = Tg;
            if (in_coll) s_min = std::min(s_min, s);
        }
        patch_min_standoff = std::isfinite(s_min) ? s_min : 0.0;
        if (surface_patch.jet_enthalpy) JetEnthalpyMarch(lev, time, k_top_col);
        // D2r-1 item 3: gated on the same key as its thermo columns, so an
        // unset run does no extra work and writes the same thermo.dat.
        if (surface_patch.side_face_factor_above >= 0.0)
            WallStepDiagnostic(lev, k_top_col, time);
    }

    // D2z: per-column recession phase, advanced once per step from the
    // domain-reduced k_top_col. A column whose k_top DROPPED since last step
    // voided at least one cell; record the interval since its previous void as
    // the TRAILING interval and stamp t_last.
    //
    // ⚠ TRAILING, not hindsight. The offline estimator can look at the actual
    // next void; the live one cannot, so it extrapolates with the previous
    // interval. Intervals vary 2-10 %, so this is second order, and the clamp
    // at phi = 1 degrades an overdue column gracefully instead of letting it
    // run away.
    //
    // A column with fewer than two voids has no interval and gets phi = 1 at
    // the use site, i.e. delta = the discrete cell count and g = 1 -- exactly
    // today's behaviour. That population is counted, never silently absorbed.
    void UpdateColumnPhase(const amrex::Vector<int>& k_top_col, Set::Scalar time)
    {
        const std::size_t ncols = k_top_col.size();
        if (col_ktop_prev.size() != ncols)
        {
            col_ktop_prev.assign(ncols, std::numeric_limits<int>::min());
            col_t_last.assign(ncols, 0.0);
            col_dt.assign(ncols, 0.0);
            col_nvoid.assign(ncols, 0);
        }
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int k = k_top_col[c];
            if (k == std::numeric_limits<int>::min()) continue;
            const int kp = col_ktop_prev[c];
            if (kp != std::numeric_limits<int>::min() && k < kp)
            {
                if (col_nvoid[c] >= 1)
                {
                    const Set::Scalar dt = time - col_t_last[c];
                    if (dt > 0.0) col_dt[c] = dt;
                }
                col_t_last[c] = time;
                ++col_nvoid[c];
            }
            col_ktop_prev[c] = k;
        }
    }

    // D2z Step 3b: ASSERT the per-cell split, do not trust it.
    //
    // The packet's split ((1-phi_T) top, 1 per middle, phi_S lowest) sums to
    // delta - 1: the phi_S piece belongs to a cell with NO void neighbour, so
    // it must be FOLDED INTO the lowest exposed cell. Without this assertion
    // the four sites can be edited "identically" and still disagree with the
    // ledger, and the disagreement would look like physics.
    //
    // Per pair: sum StepCellWeight over the exposed cells k in [kS+1, kT] and
    // require it to equal delta to round-off.
    void AssertStepSplit(const amrex::Vector<int>& k_top_col, Set::Scalar time)
    {
        if (!amrex::ParallelDescriptor::IOProcessor()) return;
        const std::size_t ncols = k_top_col.size();
        if (ncols == 0 || col_t_last.empty()) return;
        const int nx = static_cast<int>(step_split_nx);
        if (nx <= 0) return;
        const int* ktc = k_top_col.data();
        const int ny = static_cast<int>(ncols) / nx;
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const int cT = j * nx + i;
            if (ktc[cT] == std::numeric_limits<int>::min()) continue;
            const int nb[2] = {(i + 1 < nx) ? cT + 1 : -1,
                               (j + 1 < ny) ? cT + nx : -1};
            for (int d = 0; d < 2; ++d)
            {
                if (nb[d] < 0) continue;
                int ca = cT, cb = nb[d];
                if (ktc[cb] == std::numeric_limits<int>::min()) continue;
                if (ktc[cb] > ktc[ca]) { const int t = ca; ca = cb; cb = t; }
                const int kT = ktc[ca], kS = ktc[cb];
                const int n = kT - kS;
                if (n < 1) continue;
                int fb = 0;
                const Set::Scalar pT = ColPhi(col_t_last.data(), col_dt.data(),
                                              col_nvoid.data(), ca, time, fb);
                const Set::Scalar pS = ColPhi(col_t_last.data(), col_dt.data(),
                                              col_nvoid.data(), cb, time, fb);
                Set::Scalar sum = 0.0;
                for (int k = kS + 1; k <= kT; ++k)
                    sum += StepCellWeight(n, k, kT, kS, pT, pS);
                const Set::Scalar d_pair = Set::Scalar(n) + pS - pT;
                if (std::abs(sum - d_pair) > 1.0e-10 * std::max(1.0, std::abs(d_pair)))
                    Util::Abort(INFO, "D2z Step 3b split is inconsistent at columns (",
                                ca, ",", cb, ") n = ", n, ": per-cell pieces sum to ",
                                sum, " but delta = ", d_pair,
                                ". The four side sites and the ledger would disagree "
                                "and it would look like physics.");
            }
        }
    }

    // D2r-1 item 3 (owed since D2q-2c): THE DRAIN, measured. When a column
    // spalls it drops one cell, so its hot new top cell sits against the buried
    // interior of the taller column beside it and loses k*(T_top - T_wall)/dx
    // sideways. D2o-0b measured 302 kW/m^2 -- 58 % of the ~520 arriving -- and
    // it doubles when dx halves, because nothing heats the receiving rock.
    // Reported as the mean over exactly those faces, with the face count, so a
    // mesh pair can be compared directly against that number.
    // Diagnostic only: reads state, writes no field, changes no trajectory.
    void WallStepDiagnostic(int lev, const amrex::Vector<int>& k_top_col,
                            Set::Scalar time_now = 0.0)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int nx = domain.length(0);
        const int ilo_dom = domain.smallEnd(0), jlo_dom = domain.smallEnd(1);
        const int ihi_dom = domain.bigEnd(0),   jhi_dom = domain.bigEnd(1);
        const Set::Scalar* DX = geom[lev].CellSize();
        const bool void_aware = removal_enabled();
        Set::Scalar q_sum = 0.0, n_face = 0.0, q_eff = 0.0;
        // D2z-2 deliverable 2: the IN-GATE twins. This loop walks EVERY column
        // in the domain, and heating stops at surface_patch.radius, so the
        // patch edge carries permanent cliffs whose outer column has never
        // fired: nvoid < 2 => phi = 1 => delta = n, and col_dt = 0 => a = 0, so
        // StepConductance returns exactly 1. Those faces enter at g = 1 and
        // outnumber the real ones -- D2z's ungated mean g read 0.843 where the
        // in-patch value is g(0.257) = 0.335.
        // ⚠ surface_patch.radius is 0.2 on Meier (the whole domain), so it is
        // NOT usable as this gate; step_diag_radius is its own key.
        // ⚠ The gate is applied to the TALLER column -- the one whose face IS
        // the riser, and the one that sits OUTSIDE the patch at an edge cliff.
        // delta_offline.py applies the radius the same way, which is what
        // criterion (b) means by specifying both sides identically.
        const Set::Scalar sdr  = surface_patch.step_diag_radius;
        const bool gate_on     = (sdr > 0.0);
        const Set::Scalar sdr2 = sdr * sdr;
        const Set::Scalar* PLO_d = geom[lev].ProbLo();
        Set::Scalar q_in = 0.0, qe_in = 0.0, n_in = 0.0;
        Set::Scalar g_sum = 0.0, g_max = 0.0, n_multi = 0.0;
        // D2aa deliverable 4: the same multi-cell count restricted to the pad
        // ring [foot_r_inner, foot_r_outer). Required by the false-pass test --
        // without it an in-band rate cannot be told from the step drain
        // returning on the ring the burner rests on.
        Set::Scalar n_pads = 0.0, n_pads_multi = 0.0;
        const Set::Scalar ri2_p = surface_patch.foot_r_inner * surface_patch.foot_r_inner;
        const Set::Scalar ro2_p = surface_patch.foot_r_outer * surface_patch.foot_r_outer;
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp    = (*kappa_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& Df    = (*D_mf[lev]).const_array(mfi);
            // Lateral neighbours are read through the GHOST layer, exactly as
            // the conduction kernel and the side-face branch already do, so the
            // measurement does not depend on the box decomposition.
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& ph_w = (*phase_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            const int nbr[4][2] = {{-1,0},{1,0},{0,-1},{0,1}};
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(j - jlo_dom) * nx + (i - ilo_dom);
                const int k = k_top_col[c];
                if (k == std::numeric_limits<int>::min()) continue;
                if (k < lo.z || k > hi.z) continue;
                const Set::Scalar Tc = T_old(i, j, k);
                const Set::Scalar kc = kp(i, j, k) * std::exp(-kappa_damage_alpha * Df(i, j, k));
                for (const auto& d : nbr)
                {
                    const int i2 = i + d[0], j2 = j + d[1];
                    if (i2 < ilo_dom || i2 > ihi_dom || j2 < jlo_dom || j2 > jhi_dom) continue;
                    const std::size_t c2 =
                        static_cast<std::size_t>(j2 - jlo_dom) * nx + (i2 - ilo_dom);
                    const int k2 = k_top_col[c2];
                    if (k2 == std::numeric_limits<int>::min()) continue;
                    // The neighbour must be TALLER, so the cell beside this top
                    // is that column's buried interior, not its own top.
                    if (k2 <= k) continue;
                    if (void_aware && rem(i2, j2, k) > 0.5) continue;   // not solid there
                    const Set::Scalar q_raw = kc * (Tc - T_old(i2, j2, k)) / DX[0];
                    q_sum  += q_raw;
                    // D2z criterion (b): the EFFECTIVE flux, under the SAME
                    // face rule the kernel applies. This diagnostic reads the
                    // RAW cell conductivity, so without this it would move only
                    // through the temperatures and a correct conductance
                    // correction would look like it did nothing.
                    Set::Scalar g = 1.0;
                    if (surface_patch.side_face_step_height && !col_t_last.empty()
                        && (material || !phase_materials.empty()))
                    {
                        int fb = 0;
                        const Set::Scalar d =
                            PairDelta(k_top_col.data(), col_t_last.data(), col_dt.data(),
                                      col_nvoid.data(), static_cast<int>(c2),
                                      static_cast<int>(c), time_now, fb);
                        const int pcw = static_cast<int>(ph_w(i, j, k));
                        const Set::Scalar rc =
                            (!phase_materials.empty() && pcw >= 0
                             && pcw < static_cast<int>(phase_materials.size()))
                            ? phase_materials[pcw].effective_rhocp(Tc)
                            : material->effective_rhocp(Tc);
                        Set::Scalar a_d = 0.0;
                        if (rc > 0.0 && col_dt[c2] > 0.0 && kc > 0.0)
                            a_d = DX[2] * (DX[2] / col_dt[c2]) / (kc / rc);
                        g = StepConductance(d, k2 - k, a_d);
                    }
                    q_eff  += q_raw * g;
                    n_face += 1.0;
                    // The in-gate population: taller column (i2, j2) inside
                    // step_diag_radius of the patch centre. Ungated columns
                    // keep their old meaning above, untouched.
                    const Set::Scalar x2 = PLO_d[0] + (i2 + 0.5) * DX[0]
                                         - surface_patch.x0;
                    const Set::Scalar y2 = PLO_d[1] + (j2 + 0.5) * DX[1]
                                         - surface_patch.y0;
                    if (!gate_on || (x2 * x2 + y2 * y2) <= sdr2)
                    {
                        q_in  += q_raw;
                        qe_in += q_raw * g;
                        n_in  += 1.0;
                        g_sum += g;
                        if (g > g_max) g_max = g;
                        if (k2 - k >= 2) n_multi += 1.0;
                    }
                    // Pad ring: independent of step_diag_radius, because the
                    // ring is a fixed piece of hardware, not a scoring choice.
                    {
                        const Set::Scalar xp = PLO_d[0] + (i2 + 0.5) * DX[0]
                                             - surface_patch.x0;
                        const Set::Scalar yp = PLO_d[1] + (j2 + 0.5) * DX[1]
                                             - surface_patch.y0;
                        const Set::Scalar rr = xp * xp + yp * yp;
                        if (rr >= ri2_p && rr < ro2_p)
                        {
                            n_pads += 1.0;
                            if (k2 - k >= 2) n_pads_multi += 1.0;
                        }
                    }
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(q_sum);
        amrex::ParallelDescriptor::ReduceRealSum(q_eff);
        amrex::ParallelDescriptor::ReduceRealSum(n_face);
        amrex::ParallelDescriptor::ReduceRealSum(q_in);
        amrex::ParallelDescriptor::ReduceRealSum(qe_in);
        amrex::ParallelDescriptor::ReduceRealSum(n_in);
        amrex::ParallelDescriptor::ReduceRealSum(g_sum);
        amrex::ParallelDescriptor::ReduceRealSum(n_multi);
        amrex::ParallelDescriptor::ReduceRealMax(g_max);
        amrex::ParallelDescriptor::ReduceRealSum(n_pads);
        amrex::ParallelDescriptor::ReduceRealSum(n_pads_multi);
        wall_step_f_multi_pads = (n_pads > 0.0) ? n_pads_multi / n_pads : 0.0;
        wall_step_q     = (n_face > 0.0) ? q_sum / n_face : 0.0;
        wall_step_q_eff = (n_face > 0.0) ? q_eff / n_face : 0.0;
        wall_step_n_in     = n_in;
        wall_step_q_in     = (n_in > 0.0) ? q_in  / n_in : 0.0;
        wall_step_q_eff_in = (n_in > 0.0) ? qe_in / n_in : 0.0;
        step_g_mean        = (n_in > 0.0) ? g_sum / n_in : 0.0;
        step_g_max         = g_max;
        wall_step_f_multi  = (n_in > 0.0) ? n_multi / n_in : 0.0;
        wall_step_n     = n_face;
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
    // D2c (foot_body_clearance = 1, pads only): each annulus column whose
    // k_top exceeds its pad's nearest-rank k (the same pick on the integer
    // tops) is armed in foot_clear_k and cut to the pad height at the next
    // removal pass (Removal.H). Clipping values above a nearest-rank
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
        if (clear)
        {
            // D2c: annulus columns standing proud of their pad are cut to the
            // pad height at the next removal pass (Removal.H, before PASS 1).
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
    //            jet_T_ent_mode = exhaust) the previous step's T_exhaust, or
    //            (D2e, rec_fixed) the constant jet_T_rec_fixed
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
        const Set::Scalar side_r_min2 = surface_patch.side_face_r_min
                                      * surface_patch.side_face_r_min;
        // D2ab deliverable 2: the annulus pair's radius gate. ⚠ This is NOT
        // side_face_r_min: that one selects which COLUMNS take a side term at
        // all, this one selects which of them take the EXPLICIT pair rather
        // than the inherited map. Two different radii, deliberately separate.
        const Set::Scalar side_h_rmin2 = surface_patch.side_face_h_r_min
                                       * surface_patch.side_face_h_r_min;
        const bool side_area_exact = surface_patch.side_face_area_exact;
        const bool side_step_h_m = surface_patch.side_face_step_height;
        const bool side_staircase_m = surface_patch.side_face_staircase;
        const bool wetted_law_m =
            (surface_patch.side_face_area_law == SurfacePatch::SideAreaLaw::Wetted);
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;
        const bool pinned = surface_patch.robin_form == SurfacePatch::RobinForm::Pinned;
        const bool face   = surface_patch.robin_form == SurfacePatch::RobinForm::Face || pinned;
        // P0b rev5 deliverable 7, march copy of the lateral rule (rev7 d6 switch).
        const bool surf_form_m = surface_patch.robin_form ==
                                 SurfacePatch::RobinForm::Surface
                              && surface_patch.side_T_form ==
                                 SurfacePatch::SideTForm::Follow;

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

        // D2k: exposed vertical faces, so the side power comes OUT OF THE JET
        // BUDGET instead of appearing from nowhere (the D2a failure this march
        // exists to prevent; it is invisible in a fixed-h test).
        //
        // The march is serial and replicated, so every rank needs the wall
        // cells' T and kappa, and a column's wall cells are spread over the
        // z-decomposition. A dense (column x k) gather would be tens of MB per
        // step, so gather sparsely in two passes:
        //   1. per column, nslots = 1 + max(k_top - k) over its exposed cells
        //      (ReduceIntMax over ncols ints); flat columns give 0;
        //   2. prefix-sum the counts into offsets -- identical on every rank --
        //      and ReduceRealSum three arrays of that total length.
        // Only columns the march actually visits are gathered (inside the
        // patch and with s > 0), which is what keeps the total small: the
        // tall rim columns above the nozzle plane are excluded here exactly as
        // they are excluded from the march and from the kernel's side term.
        const bool side_on = surface_patch.side_face_flux;
        std::vector<int> side_nslot;
        std::vector<std::size_t> side_off;
        // P2-3: hoisted to function scope because BOTH marches read it -- the
        // floor side term (tag 1) and the wall march (tag 2). It was block-local
        // to the gather while only the gather used it.
        std::vector<unsigned char> march_col;
        std::vector<Set::Scalar> sT, sK, sN;
        std::size_t side_total = 0;
        if (side_on)
        {
            const int ihi_dom = domain.bigEnd(0), jhi_dom = domain.bigEnd(1);
            const bool void_aware = removal_enabled();
            side_nslot.assign(ncols, 0);
            // Which columns the march will visit (same test as the bin loop).
            // ⚠⚠ P2-3: TWO-VALUED NOW. 1 = below plane, which is the FLOOR
            // march's own set and the only one the bin loop and its side term
            // ever visit; 2 = ABOVE plane, gathered for the WALL march.
            //
            // Before this, the gather skipped above-plane columns outright
            // (`if (!(flame_s[c] > 0.0)) continue;`), so the slot arrays held
            // floor columns only. The wall march's first version reused those
            // slots and silently found NOTHING: `side_P_face_above` was 22.4 W
            // -- the faces plainly existed and were being heated -- while
            // `jet_P_wall` sat at exactly 0 and the fixture's (e2) row read
            // ON == OFF to four decimals. A reuse that quietly returns an empty
            // set is worse than a second enumeration, because it looks like a
            // physical result.
            //
            // Gathering both here keeps ONE enumeration with the ONE wetted-area
            // law (this file already has four side-area call sites; a fifth
            // private copy is how they drift) and leaves the floor march's
            // behaviour bit-for-bit unchanged, because every floor-side consumer
            // tests `== 1`.
            march_col.assign(ncols, 0);
            const bool want_above = surface_patch.jet_wall_march;
            for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(j) * nx + i;
                if (k_top_col[c] == std::numeric_limits<int>::min()) continue;
                const bool below = (flame_s[c] > 0.0);
                if (!below && !want_above) continue;
                const Set::Scalar x = PLO[0] + (Set::Scalar(i + ilo_dom) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j + jlo_dom) + 0.5) * DX[1];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (dxp * dxp + dyp * dyp <= r2max)
                    march_col[c] = below ? 1 : 2;
            }
            // Local pass over cells: record (column, slot, T, kappa, faces).
            struct SideHit { std::size_t col; int slot; Set::Scalar T, K, N; };
            std::vector<SideHit> hits;
            if (void_aware)
            for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& kp    = (*kappa_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& Df    = (*D_mf[lev]).const_array(mfi);
                amrex::Array4<const Set::Scalar> const& rem   = (*removed_mf[lev]).const_array(mfi);
                const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(j - jlo_dom) * nx + (i - ilo_dom);
                    if (!march_col[c]) continue;
                    if (rem(i, j, k) > 0.5) continue;
                    // D2x: the radius gate, mirroring the kernel's col_in_patch.
                    if (side_r_min2 > 0.0)
                    {
                        const Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                        const Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                        const Set::Scalar dxp2 = xc - surface_patch.x0, dyp2 = yc - surface_patch.y0;
                        if (dxp2 * dxp2 + dyp2 * dyp2 < side_r_min2) continue;
                    }
                    const bool mxm = (i > ilo_dom && rem(i-1, j, k) > 0.5);
                    const bool mxp = (i < ihi_dom && rem(i+1, j, k) > 0.5);
                    const bool mym = (j > jlo_dom && rem(i, j-1, k) > 0.5);
                    const bool myp = (j < jhi_dom && rem(i, j+1, k) > 0.5);
                    const int nfi = (mxm ? 1 : 0) + (mxp ? 1 : 0)
                                  + (mym ? 1 : 0) + (myp ? 1 : 0);
                    if (nfi == 0) continue;
                    // D2x: the march bills the jet for what the KERNEL deposits,
                    // so it must carry the same per-edge weights. The packet
                    // names two call sites; this is a third, and criterion (e)'s
                    // identity jet_P_face = tops + sides + side_loss is what
                    // proves it has to change with them.
                    Set::Scalar nf = Set::Scalar(nfi);
                    // D2z: the march must bill the jet for exactly what the
                    // kernel deposits, so it carries the same per-cell step
                    // weights. ⚠ :2507's side_f_m defaults to 1.0, which is
                    // numerically the same as the derived branch, so a missed
                    // site here is INVISIBLE except through criterion (c)'s
                    // power identity.
                    if (side_step_h_m && !col_t_last.empty())
                    {
                        const int kT = k_top_col[c];
                        int fb = 0;
                        const Set::Scalar pT = ColPhi(col_t_last.data(), col_dt.data(),
                                                      col_nvoid.data(), c, time, fb);
                        const std::size_t nb[4] = {c - 1, c + 1, c - nx, c + nx};
                        const bool hs[4] = {mxm, mxp, mym, myp};
                        nf = 0.0;
                        // P1e SITE 2 OF 4 -- the march. It bills the jet for what
                        // the kernel deposits, so it must carry the same factor or
                        // the power identity breaks.
                        int rn[4] = {0, 0, 0, 0};
                        for (int d = 0; d < 4; ++d)
                        {
                            if (!hs[d]) continue;
                            const int kS = k_top_col[nb[d]];
                            const int nn = kT - kS;
                            if (nn < 1) continue;
                            rn[d] = nn;
                            int fbS = 0;
                            const Set::Scalar pS = ColPhi(col_t_last.data(), col_dt.data(),
                                                          col_nvoid.data(), nb[d], time, fbS);
                            const Set::Scalar w = StepCellWeight(nn, k, kT, kS, pT, pS);
                            nf += AreaOfDelta(w < 0.0 ? 0.0 : w, wetted_law_m);
                        }
                        if (side_staircase_m)
                            nf *= StaircaseFactor(rn[0], rn[1], rn[2], rn[3]);
                    }
                    else if (side_area_exact)
                    {
                        const int nc2 = k_top_col[c];
                        nf = 0.0;
                        if (mxm) nf += RiserWeight(nc2 - k_top_col[c - 1],  DX[0], DX[2]);
                        if (mxp) nf += RiserWeight(nc2 - k_top_col[c + 1],  DX[0], DX[2]);
                        if (mym) nf += RiserWeight(nc2 - k_top_col[c - nx], DX[0], DX[2]);
                        if (myp) nf += RiserWeight(nc2 - k_top_col[c + nx], DX[0], DX[2]);
                    }
                    const int slot = k_top_col[c] - k;
                    if (slot < 0) continue;
                    hits.push_back({c, slot, T_old(i, j, k),
                                    kp(i, j, k) * std::exp(-kappa_damage_alpha * Df(i, j, k)),
                                    nf});
                    if (slot + 1 > side_nslot[c]) side_nslot[c] = slot + 1;
                }
            }
            amrex::ParallelDescriptor::ReduceIntMax(side_nslot.data(), static_cast<int>(ncols));
            side_off.assign(ncols, 0);
            for (std::size_t c = 0; c < ncols; ++c)
            {
                side_off[c] = side_total;
                side_total += static_cast<std::size_t>(side_nslot[c]);
            }
            if (side_total > 0)
            {
                sT.assign(side_total, 0.0);
                sK.assign(side_total, 0.0);
                sN.assign(side_total, 0.0);
                for (const auto& h : hits)
                {
                    const std::size_t idx = side_off[h.col] + static_cast<std::size_t>(h.slot);
                    sT[idx] = h.T; sK[idx] = h.K; sN[idx] = h.N;
                }
                amrex::ParallelDescriptor::ReduceRealSum(sT.data(), static_cast<int>(side_total));
                amrex::ParallelDescriptor::ReduceRealSum(sK.data(), static_cast<int>(side_total));
                amrex::ParallelDescriptor::ReduceRealSum(sN.data(), static_cast<int>(side_total));
            }
        }

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
        // D2e: rec_fixed takes the exhaust path with T_mix = jet_T_rec_fixed.
        const bool exhaust = surface_patch.jet_T_ent_exhaust;
        if (exhaust && !jet_rec_init) { jet_T_rec_next = T_ent; jet_rec_init = true; }
        const Set::Scalar T_mix = surface_patch.jet_rec_fixed ? surface_patch.jet_T_rec_fixed
                                : exhaust ? jet_T_rec_next : T_ent;
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
        const bool have_pin = pinned && pin_T.size() == ncols;
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
        Set::Scalar P_neg = 0.0, fs_cols = 0.0, T_mouth = 0.0;
        bool have_mouth = false;
        // q_robin of column c at gas temperature Tg (the kernel's inputs).
        // D2k: with side_face_flux the column draws on the jet through its
        // exposed vertical faces too. The extra is returned in the SAME units
        // as the top face (W/m^2 of dx*dy), so every P_bin += q*dA below, and
        // therefore jet_P_face and the cap, pick it up unchanged. The face
        // areas are dy*dz and dx*dz; cells are cubic when the key is on, so
        // each face contributes exactly one dA-equivalent.
        // D2x: under the derived area the weight already lives in sN, so the
        // factor here is 1 -- and side_face_factor may not be set at all.
        // D2z names the step key EXPLICITLY here. Falling through to
        // side_face_factor would give 1.0 by default -- accidentally the right
        // number -- so the omission would not show up anywhere.
        const Set::Scalar side_f_m = (side_area_exact || side_step_h_m)
            ? Set::Scalar(1.0) : surface_patch.side_face_factor;
        auto q_col = [&](std::size_t c, Set::Scalar Tg) -> Set::Scalar
        {
            int reason = 0;
            const Set::Scalar Tpin = have_pin
                ? PinRule(c, flame_h[c], Tg, time, DX[2], reason) : Set::Scalar(0.0);
            Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
            SurfaceCellFlux(true, face, 1.0, flame_h[c], Tg,
                            jkc[c], DX[2], jTc[c], eps_rad, eps_rad_hi, T_eps_break,
                            sigma_SB, h_conv, T_amb, q_in, q_out, q_rob, q_los, Tpin,
                            nullptr, SurfaceRobinT(static_cast<std::size_t>(c)));
            Set::Scalar q_tot = q_rob;
            // ⚠⚠ P2-3: `== 1` and not a truth test. The gather now also holds
            // ABOVE-plane slots (tag 2) for the wall march, and the floor march
            // must keep billing the jet for below-plane faces ONLY -- a truth
            // test here would silently double-charge every bore-wall face to
            // jet_P_face as well as jet_P_wall and break the budget identity.
            if (side_on && side_total > 0 && side_nslot[c] > 0
                && march_col[c] == 1)
            {
                const Set::Scalar w = DX[2] / DX[0];   // == 1 (cubic cells enforced)
                for (int s = 0; s < side_nslot[c]; ++s)
                {
                    const std::size_t idx = side_off[c] + static_cast<std::size_t>(s);
                    if (!(sN[idx] > 0.0)) continue;
                    Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                    // D2aa deliverable 2: the SIDE faces take the explicit pair
                    // when it is given. ⚠ This must be applied HERE as well as
                    // in the kernel and the ledger: the march accounts the side
                    // term into jet_P_face, so an override the march cannot see
                    // breaks the power identity (criterion (c)) -- which is the
                    // one check that catches a kernel<->march split.
                    // D2ab deliverables 1 and 2: the resolved temperature
                    // (side_T_eff, never the key) and the radius gate, applied
                    // here EXACTLY as in the kernel and the ledger. The column
                    // radius is rebuilt from c because the march indexes
                    // columns flat; nx, PLO, DX and the domain corner are all
                    // in scope and are the same ones BuildFlameColumns used.
                    const int ci = static_cast<int>(c % static_cast<std::size_t>(nx));
                    const int cj = static_cast<int>(c / static_cast<std::size_t>(nx));
                    const Set::Scalar cx = PLO[0]
                        + (Set::Scalar(ci + ilo_dom) + 0.5) * DX[0] - surface_patch.x0;
                    const Set::Scalar cy = PLO[1]
                        + (Set::Scalar(cj + jlo_dom) + 0.5) * DX[1] - surface_patch.y0;
                    const bool ex_c = (surface_patch.side_face_h >= 0.0)
                                   && (cx * cx + cy * cy >= side_h_rmin2);
                    const Set::Scalar h_s = ex_c ? surface_patch.side_face_h : flame_h[c];
                    // P2-3 site 3 of 5 (the FLOOR MARCH's own side term).
                    // ⚠ `slot` IS the height: slot = k_top_col[c] - k, so this
                    // row stands for cell k = k_top_col[c] - s. The march
                    // indexes by slot, not by k, so the inverse is taken here
                    // rather than assuming a `k` that is not in scope.
                    const int k_slot = k_top_col[c] - s;
                    const Set::Scalar T_s =
                        ex_c ? SideFaceTgas(static_cast<int>(c), k_slot) : Tg;
                    // P2-1 deliverable 4, ROLE 2 of 4 (the MARCH). ⚠ It bills
                    // the jet for what the kernel deposits, so a law the kernel
                    // applies and the march does not breaks the power identity
                    // -- the one check that catches a kernel<->march split.
                    // ⚠ `slot` IS the height: slot = k_top_col[c] - k, so the
                    // cell this row stands for is k_top_col[c] - s.
                    const Set::Scalar h_face =
                        SideFaceH(h_s, cx * cx + cy * cy, static_cast<int>(c),
                                  k_top_col[c] - s);
                    // P0b rev5 deliverable 7, MIRRORED. ⚠ The march accounts
                    // the side term into jet_P_face, so a rule the kernel
                    // applies and the march does not breaks the power identity
                    // -- criterion (b), and the one check that catches a
                    // kernel<->march split. Same gate, same radius, same
                    // column T_surf.
                    const bool lat_wall_m = surf_form_m
                        && (cx * cx + cy * cy >= side_h_rmin2);
                    const bool lat_face_m = surf_form_m ? true : face;
                    const Set::Scalar lat_Ts_m = (surf_form_m && !lat_wall_m)
                        ? SurfaceRobinT(c) : Set::Scalar(0.0);
                    SurfaceCellFlux(true, lat_face_m, 1.0, h_face, T_s,
                                    sK[idx], DX[0], sT[idx], eps_rad, eps_rad_hi, T_eps_break,
                                    sigma_SB, h_conv, T_amb, a_in, a_out, a_rob, a_los, Tpin,
                                    nullptr, lat_Ts_m);
                    q_tot += side_f_m * a_rob * sN[idx] * w;
                }
            }
            return q_tot;
        };
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
                    Set::Scalar q1 = q_col(c, T_in), q0 = q_col(c, T_in - 1.0);
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
                flame_Tg[c] = T_col;
                const Set::Scalar q_rob = q_col(c, T_col);
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
                for (std::size_t c : bin) r_reach = std::max(r_reach, col_r[c]);
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
                prow << time << "," << ib << "," << Set::Scalar(ib) * dr << "," << bin.size()
                     << "," << T_in << "," << P_bin << "," << share;
                if (expo)  prow << "," << T_col;
                if (fsurf) prow << "," << n_fs << "," << m_cur / m_jet;
                if (count) prow << "," << P_nb * dA;
                prow << "\n";
            }
        }
        if (prof && amrex::ParallelDescriptor::IOProcessor())
        {
            const std::string path = plot_file + "_jet_profile.csv";
            const bool first = !jet_profile_header_written;
            std::ofstream ofs(path, (first ? std::ios::trunc : std::ios::app) | std::ios::out);
            if (ofs.is_open())
            {
                if (first)
                {
                    ofs << "time,bin,r_lo,n_cols,T_gas,P_bin,P_share";
                    if (expo)  ofs << ",T_eff";
                    if (fsurf) ofs << ",n_fs,m_ratio";
                    if (count) ofs << ",P_neg";
                    ofs << "\n";
                }
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
            // previous step's T_rec) and for rec_fixed (T_mix constant) alike.
            jet_T_rec    = T_mix;
            jet_P_recirc = (m_jet - surface_patch.jet_mdot) * surface_patch.jet_cp * (T_mix - T_ent);
            P_budget     = P_nozzle + jet_P_recirc;
            // D2c: under free-surface dilution the recirculating gas is the
            // gas in the hole: the inlet of the first free-surface bin.
            jet_T_rec_next = have_mouth ? T_mouth : T_gas;
        }

        // =================================================================
        // P2-3 deliverable 2 -- THE WALL MARCH.
        // =================================================================
        // The floor march has just taken the gas across the floor and left it
        // at T_gas. Today every side face above the nozzle plane then sees ONE
        // SCALAR, `side_T_eff` = `jet_T_exhaust`, at every height: 1,504 K on
        // the P2-2f m13 leg, identically for a face 10 mm down the bore and one
        // 400 mm down. That uniform hot gas is the mechanism that ballooned
        // D2v's M2 mouth to Ø 137 mm against Meier's 96.
        //
        // Here the gas is marched UP the bore the way it is marched OUT across
        // the floor: bin the above-plane faces by height, and at each bin take
        // the enthalpy the wall actually removes, so the gas arrives at the
        // mouth cooled by the rock it has heated on the way.
        //
        // ⚠⚠ IT REUSES THE FLOOR MARCH'S OWN GATHERED SLOTS AND ADDS NO FIFTH
        // ENUMERATION. `sN/sK/sT` are already rank-identical per-column side
        // slots with the wetted-area weight in `sN`, and `slot = k_top_col[c]-k`
        // inverts to k = k_top_col[c] - s. The side-face area law has FOUR call
        // sites in this file already; a fifth private copy is how they drift.
        //
        // ⚠⚠ `jet_P_wall` IS THE ENTHALPY THE STREAM ACTUALLY LOST, not the raw
        // h*dT*A sum. When a bin would overdraw the stream past the T_ent floor
        // the excess comes from the recirculation reservoir and is NOT billed --
        // exactly the floor march's own rule. Accounting the raw sum instead
        // would break `face + wall + exhaust == cap` the moment a bin clamps.
        //
        // ⚠ NOT `jet_T_exhaust`. The packet specified that this become the mouth
        // temperature; it must not. `side_T_eff` is `jet_T_exhaust` under
        // `side_face_T = live`, and it is (i) the state the wall law derives its
        // gas properties at (:5461, site 5 of 5) and (ii) the resolver's
        // fallback for BELOW-plane faces, which are floor cone flank. Redefining
        // it would derive the wall law at mouth conditions and hand every
        // floor-flank face the cold mouth gas -- changing the floor side heating
        // D2w measured as a rate doubler. The mouth goes to `wall_T_mouth`.
        // ⚠ `jet_T_mouth` was NOT reused either: it already exists as D2c's
        // free-surface dilution temperature (set from `have_mouth` above), and
        // one column with two meanings depending on which key is on is not a
        // column.
        // ⚠ Entrainment is PARKED (user, 2026-10-13): `jet_T_rec_next` keeps the
        // FLOOR-exit value set above and the wall march does not feed it.
        Set::Scalar P_wall = 0.0;
        wall_T_gas.clear();
        wall_k0 = 0;
        wall_T_mouth = T_gas;          // floor exit, if there is no wall at all
        wall_bins    = 0.0;
        wall_L_fit   = 0.0;
        // P2-4, reset per step like every other per-step state above.
        // ⚠ `body_noLaw_steps` is the ONE exception: it is CUMULATIVE, so the
        // reduction differences it over a window and never sums it.
        jet_P_body     = 0.0;
        body_h_used    = 0.0;
        body_z_imm     = 0.0;
        body_L         = 0.0;
        jet_P_retained = 0.0;
        if (surface_patch.jet_wall_march && side_on && side_total > 0)
        {
            // ⚠ named `side_hmax`, the SAME local name the kernel (:2033) and
            // the ledger (:4514) use. `unit/lateral_flux_lock` whitelists the
            // right-hand sides a cap assignment may take, so a private spelling
            // here would read as a fifth, unvetted cap -- which is precisely
            // the distinction that lock exists to keep.
            const Set::Scalar side_hmax =
                (surface_patch.side_face_flux && spall_per_cell_removal)
                ? surface_patch.side_face_h_max : Set::Scalar(-1.0);
            const Set::Scalar side_f_ab = (surface_patch.side_face_factor_above >= 0.0)
                ? surface_patch.side_face_factor_above
                : surface_patch.side_face_factor;
            // ⚠ the gather's OWN tag, not a second reading of flame_s: the two
            // must agree about which columns have slots, and the tag is what
            // actually decided whether slots were written.
            auto above = [&](std::size_t c) {
                return c < march_col.size() && march_col[c] == 2;
            };
            // pass 1: the height range the above-plane faces span.
            int kmin = std::numeric_limits<int>::max();
            int kmax = std::numeric_limits<int>::min();
            for (std::size_t c = 0; c < ncols; ++c)
            {
                if (!above(c) || side_nslot[c] <= 0) continue;
                for (int sl = 0; sl < side_nslot[c]; ++sl)
                {
                    const std::size_t idx = side_off[c] + static_cast<std::size_t>(sl);
                    if (!(sN[idx] > 0.0)) continue;
                    const int k = k_top_col[c] - sl;
                    if (k < kmin) kmin = k;
                    if (k > kmax) kmax = k;
                }
            }
            if (kmax >= kmin)
            {
                wall_k0 = kmin;
                wall_T_gas.assign(static_cast<std::size_t>(kmax - kmin + 1), 0.0);
                std::vector<Set::Scalar> Pk(wall_T_gas.size(), 0.0);
                std::vector<long>        nk(wall_T_gas.size(), 0);
                std::vector<Set::Scalar> Pb(wall_T_gas.size(), 0.0);
                // ===== P2-4: THE COLD BODY IN THE ANNULUS =====
                // ⚠⚠ THE AREA CARRIES THE SECTOR FRACTION. The deck is a
                // QUARTER annulus, so the body's wetted perimeter per bin is
                // 2*pi*r*wall_sector*dz, not 2*pi*r*dz. Dropping `wall_sector`
                // would read as 4x the cooling -- which is why the fixture's
                // (f) row checks the fitted length against `body_L` rather than
                // only checking that the gas got colder.
                // ⚠ `wall_sector` is DERIVED (PatchAzimuthFraction), not typed,
                // so the quarter assumption is an output and not a comment.
                const Set::Scalar body_P_q =
                    2.0 * 3.14159265358979323846 * WJ_R_SHIELD * wall_sector;
                const Set::Scalar body_A_bin = body_P_q * DX[2];
                // ⚠ ONE READ of the law's duct scalar and NO typed fallback.
                // The body stands in the annulus, so the duct branch is its
                // coefficient by construction -- the wall-jet branch below
                // z_fill belongs to the ROCK face, and P2-1's composite p50
                // (372, then 262->167 on W1) is not the body's number.
                // When the law is not ready the term is zero and the step is
                // COUNTED (body_noLaw_n), because an inert sink that nothing
                // reports is the failure mode three keys here already hit.
                const bool body_on = surface_patch.jet_body;
                const Set::Scalar body_h = (body_on && wall_law_ready)
                    ? surface_patch.jet_body_h_scale * wall_h_duct
                    : Set::Scalar(0.0);
                if (body_on && !wall_law_ready) body_noLaw_steps += 1.0;
                // the mass the WALL section carries. Default `jet` keeps mcp,
                // so a leg that does not name the key is byte-identical.
                const Set::Scalar mcp_wall =
                    (surface_patch.jet_wall_mdot == SurfacePatch::WallMdot::Fresh)
                        ? surface_patch.jet_mdot * surface_patch.jet_cp : mcp;
                Set::Scalar P_body = 0.0;
                // pass 2: march UPWARD, lowest above-plane face first -- the
                // gas enters the bore at the floor and leaves at the mouth.
                for (int k = kmin; k <= kmax; ++k)
                {
                    const std::size_t b = static_cast<std::size_t>(k - kmin);
                    const Set::Scalar T_in = T_gas;
                    // the bin's INLET temperature: the value its faces receive,
                    // which is what SideFaceTgas hands the kernel for this k.
                    wall_T_gas[b] = T_in;
                    Set::Scalar q_sum = 0.0;
                    for (std::size_t c = 0; c < ncols; ++c)
                    {
                        if (!above(c) || side_nslot[c] <= 0) continue;
                        const int sl = k_top_col[c] - k;
                        if (sl < 0 || sl >= side_nslot[c]) continue;
                        const std::size_t idx = side_off[c] + static_cast<std::size_t>(sl);
                        if (!(sN[idx] > 0.0)) continue;
                        const int ci = static_cast<int>(c % static_cast<std::size_t>(nx));
                        const int cj = static_cast<int>(c / static_cast<std::size_t>(nx));
                        const Set::Scalar cx = PLO[0]
                            + (Set::Scalar(ci + ilo_dom) + 0.5) * DX[0] - surface_patch.x0;
                        const Set::Scalar cy = PLO[1]
                            + (Set::Scalar(cj + jlo_dom) + 0.5) * DX[1] - surface_patch.y0;
                        const Set::Scalar r2c = cx * cx + cy * cy;
                        const Set::Scalar h_s = (surface_patch.side_face_h >= 0.0)
                            ? surface_patch.side_face_h : flame_h[c];
                        Set::Scalar h_face =
                            SideFaceH(h_s, r2c, static_cast<int>(c), k);
                        if (side_hmax > 0.0 && h_face > side_hmax)
                            h_face = side_hmax;
                        if (!(h_face > 0.0)) continue;
                        // ⚠ the packet's form: h*(T_gas - T_face), one-sided.
                        // A face hotter than the gas does not return enthalpy to
                        // the stream, exactly as the floor march's non-count
                        // branch refuses q_robin < 0.
                        const Set::Scalar dT = T_in - sT[idx];
                        if (!(dT > 0.0)) { ++nk[b]; continue; }
                        // same units as the top face: W/m^2 of dx*dy, with the
                        // wetted weight already in sN and w = dz/dx == 1.
                        q_sum += side_f_ab * h_face * dT * sN[idx];
                        ++nk[b];
                    }
                    const Set::Scalar P_raw = q_sum * dA;
                    // P2-4: the body's debit for THIS bin, at the bin's INLET
                    // temperature, beside the rock faces and one-sided exactly
                    // as they are -- a body colder than the gas takes heat, a
                    // body hotter than it does not give any back.
                    const Set::Scalar dTb = T_in - surface_patch.jet_body_T;
                    const Set::Scalar P_body_k =
                        (body_h > 0.0 && dTb > 0.0) ? body_h * body_A_bin * dTb
                                                    : Set::Scalar(0.0);
                    // the stream can only be cooled to the T_ent floor; the
                    // remainder is reservoir enthalpy and is not billed.
                    // ⚠ BOTH sinks are clamped TOGETHER against that floor, and
                    // then split in proportion to what each asked for. Clamping
                    // them separately would let the pair overdraw the stream and
                    // break `face + wall + body + exhaust == cap`.
                    const Set::Scalar P_ask = P_raw + P_body_k;
                    const Set::Scalar T_out = (mcp_wall > 0.0)
                        ? std::max(T_ent, T_in - P_ask / mcp_wall) : T_in;
                    const Set::Scalar P_got =
                        (mcp_wall > 0.0) ? mcp_wall * (T_in - T_out) : 0.0;
                    const Set::Scalar frac = (P_ask > 0.0) ? P_got / P_ask : 0.0;
                    Pk[b]   = P_raw * frac;
                    Pb[b]   = P_body_k * frac;
                    P_wall += Pk[b];
                    P_body += Pb[b];
                    T_gas   = T_out;
                }
                wall_T_mouth = T_gas;
                jet_P_body   = P_body;
                body_h_used  = body_h;
                body_z_imm   = Set::Scalar(wall_T_gas.size()) * DX[2];
                body_L       = (body_h > 0.0 && body_P_q > 0.0)
                    ? mcp_wall / (body_h * body_P_q) : Set::Scalar(0.0);
                wall_bins    = Set::Scalar(wall_T_gas.size());
                // the cooling length, by least squares on ln(T - T_ent) vs z
                // over the bins that actually cooled. ⚠ reported as 0 when
                // fewer than two such bins exist rather than fitted to one
                // point, which would return an arbitrary length.
                {
                    Set::Scalar Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
                    long n = 0;
                    for (std::size_t b = 0; b < wall_T_gas.size(); ++b)
                    {
                        const Set::Scalar dTb = wall_T_gas[b] - T_ent;
                        if (!(dTb > 0.0)) continue;
                        const Set::Scalar z = Set::Scalar(b) * DX[2];
                        const Set::Scalar y = std::log(dTb);
                        Sx += z; Sy += y; Sxx += z * z; Sxy += z * y; ++n;
                    }
                    if (n >= 2)
                    {
                        const Set::Scalar den = Set::Scalar(n) * Sxx - Sx * Sx;
                        if (std::abs(den) > 0.0)
                        {
                            const Set::Scalar slope =
                                (Set::Scalar(n) * Sxy - Sx * Sy) / den;
                            if (slope < 0.0) wall_L_fit = -1.0 / slope;
                        }
                    }
                }
                // the per-bin profile, rank 0, on the jet profile interval --
                // the same `prof` latch the floor profile uses, so the two files
                // are written on the same steps and can be joined by time.
                if (prof && amrex::ParallelDescriptor::IOProcessor())
                {
                    std::ostringstream wn;
                    wn << plot_file << "_wall_profile.csv";
                    const bool fresh = !std::ifstream(wn.str().c_str()).good();
                    std::ofstream wf(wn.str().c_str(), std::ios::app);
                    if (fresh) wf << "time,k,z,n_faces,T_gas,P_k,P_body_k\n";
                    wf.setf(std::ios::scientific); wf.precision(16);
                    for (std::size_t b = 0; b < wall_T_gas.size(); ++b)
                        wf << time << "," << (kmin + static_cast<int>(b)) << ","
                           << (Set::Scalar(b) * DX[2]) << "," << nk[b] << ","
                           << wall_T_gas[b] << "," << Pk[b] << ","
                           << Pb[b] << "\n";
                }
            }
        }
        jet_P_wall = P_wall;

        jet_P_decay   = P_budget - jet_P_cap;
        jet_P_face    = P_face;
        // ⚠⚠ `jet_T_exhaust` KEEPS ITS MEANING: the gas state leaving the FLOOR.
        // With the wall march on, `T_gas` has since been marched up the bore, so
        // the floor-exit value is recovered from the stream's own bookkeeping
        // rather than by reordering the assignments -- P_wall is by construction
        // mcp*(T_floor_exit - T_mouth).
        // P2-4: the floor-exit value is recovered from BOTH wall-section sinks.
        // ⚠ mcp_w here must be the same mass the bin loop marched, or the
        // recovered floor exit is wrong under `fresh`.
        const Set::Scalar mcp_w =
            (surface_patch.jet_wall_mdot == SurfacePatch::WallMdot::Fresh)
                ? surface_patch.jet_mdot * surface_patch.jet_cp : mcp;
        jet_T_exhaust = (mcp_w > 0.0)
            ? (T_gas + (P_wall + jet_P_body) / mcp_w) : T_gas;
        // the stream that actually leaves the domain is at the MOUTH, so the
        // exhaust term is taken there and `face + wall + exhaust == cap` holds.
        // ⚠⚠ UNDER `fresh` ONLY THE FRESH MASS LEAVES BY THE MOUTH. The
        // entrained share never climbed the bore, so it is booked as RETAINED
        // and the identity becomes
        //   face + wall + body + exhaust_mouth + retained == cap,
        // which holds because
        //   mcp_jet (T_fe - T_ent) = mcp_fresh (T_fe - T_ent) + retained
        // and the fresh part then pays wall + body + exhaust_mouth.
        jet_P_exhaust = mcp_w * (T_gas - T_ent);
        jet_P_retained = (mcp_w == mcp)
            ? Set::Scalar(0.0)
            : (mcp - mcp_w) * (jet_T_exhaust - T_ent);
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
              // P2-3: the wall is a third sink on the SAME stream.
              // P2-4: the body is a fourth, and under `fresh` the retained
              // entrained share is a fifth term that never leaves.
              P_wall >= -tol && jet_P_body >= -tol && jet_P_retained >= -tol &&
              std::abs(P_face + P_wall + jet_P_body + jet_P_exhaust
                       + jet_P_retained - jet_P_cap) <= tol &&
              jet_P_decay >= -tol &&
              (surface_patch.jet_stag_decay && !surface_patch.jet_entrained_mass
                   ? true : std::abs(jet_P_decay) <= tol)))
            Util::Abort(INFO, "jet enthalpy invariant violated at time ", time,
                        ": jet_P_face = ", P_face, ", jet_P_exhaust = ", jet_P_exhaust,
                        ", jet_P_cap = ", jet_P_cap, ", jet_P_decay = ", jet_P_decay,
                        ", jet_P_neg = ", P_neg, ", jet_P_wall = ", P_wall,
                        ", jet_P_body = ", jet_P_body,
                        ", jet_P_retained = ", jet_P_retained,
                        ", nozzle budget = ", P_nozzle);
    }

    // C1: per-column pin used by this step: pin_T_eff[col] = pin_T[col] unless
    //   the column never fired (pin_T == 0), the flame cannot heat it past
    //   the pin (q_pin = h(T_gas - pin_T) <= 0), or it has been idle for
    //   time - pin_t >= pinned_idle_cycles * t_cell, with the cell heating
    //   time t_cell = pin_rho*pin_Cp*dz*(pin_T - T_amb)/q_pin.
    // Removal (which updates pin_T/pin_t) runs before the thermal advance in
    // Advance(), so a column that fired this step has time - pin_t = 0 and
    // pinned_idle_cycles = 0 never pins (face-identical test hook).
    // P0b deliverable 1.1. Called FIRST in Advance(), so the array it produces
    // is read by UpdateSpAfterMechanics (criterion site A), by
    // UpdateRemovalAfterCohesive (criterion site B and the removal) and by the
    // conduction kernel / jet march / ledger (the BC) -- all within one step,
    // all from the same temperatures. Advance() swaps temp_mf into temp_old_mf
    // AFTER removal, so temp_mf here IS the kernel's T_old.
    //
    // ⚠ The two-cell read goes through an owner-write + ReduceRealSum column
    // gather (the jTc/jkc pattern), never a k-1 neighbour read: the pre-P0 code
    // clamps its sample index into the validbox, so a direct read would make
    // the criterion depend on amr.max_grid_size at every box seam.
    // ⚠⚠ P2-2e: THE GATHERED CONSTRAINT REFERENCE. One live-temperature band per
    // column, from the column top down past h_spall_max + d, gathered across the
    // z-decomposition (the T_col_prof pattern: each cell written by the ONE rank
    // whose validbox holds it, then ReduceRealSum). Read by BOTH vertical firing
    // sites through CrefTAtDepth, so they cannot disagree.
    // ⚠ NEVER a box-local read: the Sp_field site's own deeper branch clamps at
    // the BOX bound (vbx_lo_k) and would hit it d/dz cells down.
    // ⚠ The same column-top rule as both sites (topmost cell with removed_mf <=
    // 0.5, ReduceIntMax -- GatherColumnSurface), and the same contiguous-solid-
    // run guard as the depth scan (ColumnSolidRunBelowTop), so the reference
    // never reads through an interior void as if it were constraining rock.
    void BuildConstraintRefColumns(int lev)
    {
        // ⚠ reset PER STEP, before the early return, so a non-owning level or a
        // d = 0 run never leaves a stale count in thermo.dat.
        cref_n_vert = cref_n_lat = cref_n_clamp = cref_n_neg = 0.0;
        cref_d_out  = sp_cref_d;
        if (!(sp_cref_d > 0.0) || !SurfaceOwner(lev)) return;

        const amrex::Box& domain = geom[lev].Domain();
        const int dom_ilo = domain.smallEnd(0), dom_jlo = domain.smallEnd(1);
        const int nx_dom  = domain.length(0);
        const int ncols   = nx_dom * domain.length(1);
        const int top_k   = domain.bigEnd(2);
        const Set::Scalar dz = geom[lev].CellSize()[2];

        cref_ktop.assign(ncols, top_k);
        GatherColumnSurface(lev, cref_ktop, dom_ilo, dom_jlo, nx_dom, top_k);
        const amrex::Vector<int> run = ColumnSolidRunBelowTop(lev, cref_ktop);

        // The band. The scan's deepest sample sits h_spall_max (+ the flaw's own
        // a_f) below the face and its reference d further: cells
        // 0 .. floor((h_spall_max + d)/dz) are the packet's
        // `floor((h_spall_max + d)/dz) + 1`, and ONE more covers the a_f
        // overhang so a read at the deepest crack interpolates rather than
        // clamping. Clamps are still counted (cref_n_clamp).
        cref_colw = static_cast<int>(std::floor((spall_h_spall_max + sp_cref_d) / dz)) + 2;
        cref_T_col.assign(static_cast<std::size_t>(ncols) * cref_colw, 0.0);
        for (amrex::MFIter mfi(*temp_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).const_array(mfi);
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
            {
                const int c  = (j - dom_jlo) * nx_dom + (i - dom_ilo);
                const int kt = cref_ktop[c];
                if (kt == std::numeric_limits<int>::min()) continue;
                for (int m = 0; m < cref_colw; ++m)
                {
                    const int k = kt - m;
                    if (k < klo || k > khi) continue;
                    cref_T_col[static_cast<std::size_t>(c) * cref_colw + m] = T(i, j, k);
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(
            cref_T_col.data(), static_cast<int>(cref_T_col.size()));

        // Usable cells per column: the band, the domain floor and the solid run.
        // `run` already stops at the domain floor (kt - klo + 1 with no void).
        cref_nuse.assign(ncols, 0);
        for (int c = 0; c < ncols; ++c)
            if (cref_ktop[c] != std::numeric_limits<int>::min())
                cref_nuse[c] = std::min(cref_colw, run[c]);
    }

    // T at VERTICAL depth D [m] below the TOP FACE of the buffer's own column top,
    // linear between cell centres ((m + 1/2) dz) -- METRES, never cells. A read
    // beyond the last usable cell (band, domain floor or the end of the contiguous
    // solid run) returns that cell and is COUNTED. A read above the top cell
    // centre returns the top cell: it is not a floor clamp, and for any d above
    // half a cell it cannot occur.
    // ⚠ Callers pass D measured from THEIR face plus (cref_ktop[c] - k_site)*dz,
    // which is 0 whenever the two tops agree -- and they do by construction (one
    // rule, one state) -- but keeps the read physically exact if they ever did not.
    Set::Scalar CrefTAtDepth(std::size_t c, Set::Scalar D, Set::Scalar dz,
                             long& n_read, long& n_clamp) const
    {
        ++n_read;
        const int nuse = cref_nuse[c];
        if (nuse <= 0)
            Util::Abort(INFO, "P2-2e: the constraint reference was read on a column "
                              "with no gathered top -- a firing site scored a column "
                              "BuildConstraintRefColumns did not see");
        const Set::Scalar* g = cref_T_col.data() + c * static_cast<std::size_t>(cref_colw);
        Set::Scalar mf = D / dz - 0.5;
        if (mf < 0.0) mf = 0.0;
        const Set::Scalar mlast = static_cast<Set::Scalar>(nuse - 1);
        if (mf >= mlast)
        {
            if (mf > mlast) ++n_clamp;
            return g[nuse - 1];
        }
        const int m0 = static_cast<int>(std::floor(mf));
        const Set::Scalar w = mf - static_cast<Set::Scalar>(m0);
        return (1.0 - w) * g[m0] + w * g[m0 + 1];
    }

    void BuildSurfaceTColumns(int lev, Set::Scalar dt)
    {
        if (!sp_surface_T_recon || !SurfaceOwner(lev)) return;
        // ⚠ reset PER STEP: RegisterIntegratedVariable integrates what it is
        // given, so a counter that kept accumulating here would be integrated
        // twice and read as a meaningless growing number.
        sp_len_fb_h0 = 0.0;
        sp_len_fb_q0 = 0.0;
        const amrex::Box& domain = geom[lev].Domain();
        const int dom_ilo = domain.smallEnd(0), dom_jlo = domain.smallEnd(1);
        const int nx_dom = domain.length(0);
        const int ncols  = nx_dom * domain.length(1);
        const int top_k  = domain.bigEnd(2);
        const Set::Scalar dz = geom[lev].CellSize()[2];

        amrex::Vector<int> k_top_col(ncols, top_k);
        GatherColumnSurface(lev, k_top_col, dom_ilo, dom_jlo, nx_dom, top_k);
        const amrex::Vector<int> solid_run = ColumnSolidRunBelowTop(lev, k_top_col);

        std::vector<Set::Scalar> Ttop(ncols, 0.0), T1(ncols, 0.0);
        // ⚠⚠ P2-2d deliverable 6: THE PER-CELL CONDUCTIVITY FOR THE VERTICAL
        // FLAW READ. `ell_flux` below used `phases[0].kappa` -- the deck's
        // per-phase SCALAR -- while the LATERAL read (Removal.H, `kc_cell`)
        // uses the per-cell field. They agree under `law = constant` and
        // DIVERGE under the Grimsel law, where lambda halves by 600 degC. That
        // is the same per-phase-constant class P2-2b's zero-exception rule
        // removed for rho/c_p; it survived here only because `(a1)`'s predicate
        // scans `rho` and `Cp` and not `kappa`.
        //
        // ⚠⚠ THE EVALUATION TEMPERATURE IS THE TOP CELL CENTRE, STATED. The
        // field is cell-centred, so this is lambda(T_cell), not lambda(T_surf).
        // It matters: `ell_q^3` spans 8.36e-9 (scalar k = 1.5) to ~1.34e-8
        // (lambda at an 850 K surface) to ~1.67e-8 (a cell centre nearer
        // 768 K), so `V_ref` moves ~25 % on this choice alone. Cell-centre is
        // chosen to MATCH THE LATERAL READ, which is what makes the two
        // comparable at all -- not because the surface value is wrong.
        std::vector<Set::Scalar> Ktop(ncols, 0.0);
        for (amrex::MFIter mfi(*temp_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& vbx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& KA = (*kappa_mf[lev]).const_array(mfi);
            const int vklo = vbx.smallEnd(2), vkhi = vbx.bigEnd(2);
            for (int j = vbx.smallEnd(1); j <= vbx.bigEnd(1); ++j)
            for (int i = vbx.smallEnd(0); i <= vbx.bigEnd(0); ++i)
            {
                const int c  = (j - dom_jlo) * nx_dom + (i - dom_ilo);
                const int kt = k_top_col[c];
                if (kt >= vklo && kt <= vkhi)  Ttop[c] = T(i, j, kt);
                if (kt >= vklo && kt <= vkhi)  Ktop[c] = KA(i, j, kt);
                if (kt - 1 >= vklo && kt - 1 <= vkhi &&
                    kt - 1 >= domain.smallEnd(2))  T1[c] = T(i, j, kt - 1);
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(Ttop.data(), ncols);
        amrex::ParallelDescriptor::ReduceRealSum(T1.data(),   ncols);
        amrex::ParallelDescriptor::ReduceRealSum(Ktop.data(), ncols);

        if (static_cast<int>(surf_x_hold.size()) != ncols)
            surf_x_hold.assign(ncols, 0.0);
        surf_T_col.assign(ncols, 0.0);
        surf_ell_col.assign(ncols, 0.0);
        surf_ellq_col.assign(ncols, 0.0);
        surf_ell_raw.assign(ncols, 0.0);
        surf_hf_col.assign(ncols, 0.0);
        surf_tbar_col.assign(ncols, 0.0);
        surf_flag_col.assign(ncols, 0.0);

        const Set::Scalar cap = SurfaceGasTempCap();
        const Set::Scalar T0  = phases.empty() ? 0.0 : phases[0].T_ref;
        // P0b rev5 deliverable 3: THE SMOOTHING IS GONE. Revision 2 ran the
        // raw x = dz/ell_eff through a 1 s exponential filter and held the
        // step after a firing. Both are deleted, and the packet's own evidence
        // is why: at h_conv = 500 the onset takes 0.3-0.4 s, so a 1 s filter
        // is longer than the event it smooths. It made the column fire at
        // 0.30-0.81 s instead of P0's 1.44-1.81 s -- i.e. before the two-cell
        // estimate had any content -- which drove 2 of 3 legs onto the
        // ell >= dz/4 floor and turned P0's 0.116 % anchor spread into 2.089 %.
        // A filter whose time constant is not derived from the problem is a
        // tuning constant; the raw estimate is the honest one.
        // ⚠ THE FLAT-HOLD STAYS (1.3): it is a different mechanism -- when the
        // two-cell pair carries no usable gradient the column holds its last
        // VALID x rather than collapsing to the cell average. That is what
        // took Meier's `flat` fallback share from 12.8 % to 0.0 %.
        (void) dt;

        for (int c = 0; c < ncols; ++c)
        {
            const Set::Scalar x_raw = (solid_run[c] >= 2)
                ? RawSurfaceX(Ttop[c], T1[c], T0) : 0.0;
            surf_ell_raw[c] = (x_raw > 0.0) ? dz / x_raw : 0.0;
            // 1.3, the flat-hold: a usable pair overwrites the held value, a
            // flat pair keeps it. With the EMA gone this IS the raw estimate
            // on every step that has one -- which is the witness the packet
            // asks for (surf_ell must now equal the raw column).
            if (x_raw > 0.0) surf_x_hold[c] = x_raw;
            const Set::Scalar x_use = surf_x_hold[c];

            // ⚠⚠ P1g-3 deliverables 5 + 7: THIS SITE'S OWN NET FLUX. The floor
            // sees the flame/jet pair, so q_net is
            //     h (T_g - T_s) - eps sig (T_s^4 - T_amb^4) - h_conv (T_s - T_amb)
            // the SAME net form `Removal.H:1503-1517` already uses for q_real --
            // net of radiation and convection, not the delivered flux. Using the
            // delivered flux here and the net flux there would silently compare
            // two different quantities.
            // ⚠ Two fallbacks to `twocell`, each COUNTED, never a clamp:
            // `flame_h` is 0 outside the patch and above the nozzle plane, and
            // q_net can be <= 0 on a cooling face.
            Set::Scalar ell_flux = 0.0;
            if (sp_flake_len_flux)
            {
                const Set::Scalar Ts_pre = ReconFromX(Ttop[c], T0, dz, x_use,
                                                      cap, sp_flake_c_f,
                                                      spall_h_spall_max,
                                                      sp_flake_h_max).T_surf;
                const Set::Scalar hh = (c < static_cast<int>(flame_h.size()))
                                     ? flame_h[c] : 0.0;
                const Set::Scalar tg = (c < static_cast<int>(flame_Tg.size()))
                                     ? flame_Tg[c] : 0.0;
                if (!(hh > 0.0))
                {
                    sp_len_fb_h0 += 1.0;
                }
                else
                {
                    // ⚠ the loss constants are direct members (:7920-7925),
                    // not a `losses` sub-object -- the input KEYS are
                    // `losses.*` but the C++ names are bare.
                    const Set::Scalar ep =
                        (Ts_pre <= T_eps_break) ? eps_rad : eps_rad_hi;
                    const Set::Scalar q_net =
                        hh * (tg - Ts_pre)
                        - ep * sigma_SB
                          * (Ts_pre * Ts_pre * Ts_pre * Ts_pre
                             - T_amb * T_amb * T_amb * T_amb)
                        - h_conv * (Ts_pre - T_amb);
                    // ⚠ P2-2d d6: the PER-CELL conductivity at the top cell
                    // centre, matching the lateral read. Falls back to the
                    // per-phase scalar only if the field is unavailable (it is
                    // written at init for every phase under every law, so this
                    // is belt-and-braces, not a routine path).
                    const Set::Scalar k_use =
                        (c < static_cast<int>(Ktop.size()) && Ktop[c] > 0.0)
                        ? Ktop[c] : phases[0].kappa;
                    if (q_net > 0.0 && Ts_pre > T0)
                        ell_flux = k_use * (Ts_pre - T0) / q_net;
                    else
                        sp_len_fb_q0 += 1.0;
                }
            }
            const Set::Scalar frac_over = sp_flake_frac_weibull
                ? (W_CONE_MEDIAN / weibull_m) : 0.0;
            const SurfaceRecon r = ReconFromX(Ttop[c], T0, dz, x_use, cap,
                                              sp_flake_c_f, spall_h_spall_max,
                                              sp_flake_h_max,
                                              ell_flux, frac_over);
            surf_T_col[c]    = r.T_surf;
            surf_ell_col[c]  = std::isfinite(r.ell_eff) ? r.ell_eff : 0.0;
            surf_ellq_col[c] = std::isfinite(r.ell_q) ? r.ell_q : 0.0;
            surf_hf_col[c]   = std::isfinite(r.h_f) ? r.h_f : 0.0;
            surf_tbar_col[c] = r.T_bar;
            surf_flag_col[c] = (r.flat ? 1.0 : 0.0)
                             + (r.ell_clamped ? 2.0 : 0.0)
                             + (r.T_clamped ? 4.0 : 0.0)
                             + (r.h_f_capped ? 8.0 : 0.0)
                             + (r.h_f_bounded ? 16.0 : 0.0);   // rev7 d1
        }
    }

    // P0b deliverable 1.1: the per-column surface temperature the Robin BC
    // uses, or 0 ("not the surface form") -- the same channel shape as
    // pin_T_eff/tpe. ⚠ TOP FACES ONLY. Deliverable 1.2: a lateral face has no
    // column profile, so side-face callers pass nothing and keep the cell
    // average, exactly as P0 did for lateral criterion faces.
    // P0b: hand back the per-column reconstruction built at the top of the
    // step. ⚠ This REPLACES the two separate gathers P0 used (one in
    // UpdateSpAfterMechanics, one in Removal.H's PASS 1b). They agreed, but
    // they were two copies of the same estimate that could drift, and neither
    // could carry the smoothing or the flat-hold, which are per-column state
    // across steps.
    bool SharedSurfaceRecon(std::size_t c, SurfaceRecon& r) const
    {
        if (c >= surf_T_col.size() || !(surf_T_col[c] > 0.0)) return false;
        r.T_surf  = surf_T_col[c];
        r.ell_eff = (surf_ell_col[c] > 0.0) ? surf_ell_col[c]
                    : std::numeric_limits<Set::Scalar>::infinity();
        r.h_f     = surf_hf_col[c];
        r.ell_q   = (c < surf_ellq_col.size()) ? surf_ellq_col[c] : 0.0;
        r.T_bar   = surf_tbar_col[c];
        const int f = static_cast<int>(surf_flag_col[c] + 0.5);
        r.flat        = (f & 1) != 0;
        r.ell_clamped = (f & 2) != 0;
        r.T_clamped   = (f & 4) != 0;
        r.h_f_capped  = (f & 8) != 0;
        r.h_f_bounded = (f & 16) != 0;
        return true;
    }

    // ⚠⚠ P1 -- THE REFERENCE TEMPERATURE, READ FROM THE EXPONENTIAL.
    //
    // The reference sits at a PHYSICAL depth below the surface, so it must be
    // read from a profile that does not know about the mesh. The gathered
    // column interpolates linearly between CELL CENTRES, so the same physical
    Set::Scalar SurfaceRobinT(std::size_t c) const
    {
        if (surface_patch.robin_form != SurfacePatch::RobinForm::Surface)
            return 0.0;
        return (c < surf_T_col.size()) ? surf_T_col[c] : 0.0;
    }

    void BuildPinEffective(int lev, Set::Scalar time)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const std::size_t ncols = static_cast<std::size_t>(domain.length(0)) * domain.length(1);
        if (pin_T.size() != ncols)
        {
            pin_T.assign(ncols, 0.0);
            pin_t.assign(ncols, 0.0);
            pin_rho.assign(ncols, 0.0);
            pin_Cp.assign(ncols, 0.0);
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
    // the D2a2 JetEnthalpyMarch. pin_T/pin_t/pin_rho/pin_Cp must be sized.
    Set::Scalar PinRule(std::size_t c, Set::Scalar h, Set::Scalar T_gas,
                        Set::Scalar time, Set::Scalar dz, int& reason) const
    {
        reason = 0;
        if (pin_T[c] <= 0.0) return 0.0;                                // reason 0: unset
        const Set::Scalar q_pin = h * (T_gas - pin_T[c]);
        if (q_pin <= 0.0) { reason = 2; return 0.0; }                   // q_pin <= 0
        // P2-2b site 2: the product is formed HERE, at the point of use, from
        // two separately carried per-cell values. ⚠ Under `law = constant` the
        // fields hold the phase scalars, so this is arithmetically identical to
        // the retired stored product -- which is what (b7) measures.
        const Set::Scalar t_cell =
            pin_rho[c] * pin_Cp[c] * dz * (pin_T[c] - T_amb) / q_pin;
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
                                pinned ? pin_T_eff[col] : Set::Scalar(0.0), &branch,
                                SurfaceRobinT(static_cast<std::size_t>(col)));
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
    // D2x: the per-edge riser weight. A surface descending Dz across one cell
    // of width dx has true area dy*sqrt(dx^2 + Dz^2); the TOP face already
    // carries dx*dy, so the riser carries the remainder,
    // dy*(sqrt(dx^2 + Dz^2) - dx), spread over the n cells of that riser, each
    // of face area dx*dz:   w = (sqrt(dx^2 + (n*dz)^2) - dx) / (n*dz)
    // PER DIRECTION -- each lateral neighbour has its own drop, so the four
    // neighbours contribute four independent weights and must not be lumped.
    // Cubic cells give w(1) = 0.4142, w(2) = 0.6180, w(3) = 0.7208 -> 1.
    // w >= 0.4142 BY CONSTRUCTION: a riser face exists only where the surface
    // drops at least one full cell, so the local slope is at least 1. The
    // f = 0.2 it replaces is BELOW the minimum the geometry permits.
    static Set::Scalar RiserWeight(int n, Set::Scalar dx, Set::Scalar dz)
    {
        const Set::Scalar d = Set::Scalar(n < 1 ? 1 : n) * dz;
        return (std::sqrt(dx * dx + d * d) - dx) / d;
    }

    static void SurfaceCellFlux(bool in_patch, bool face_form, Set::Scalar surf_c,
                                Set::Scalar h, Set::Scalar T_gas,
                                Set::Scalar kc, Set::Scalar dz, Set::Scalar Tc,
                                Set::Scalar eps_lo, Set::Scalar eps_hi, Set::Scalar T_eps_b,
                                Set::Scalar sig_SB, Set::Scalar h_c, Set::Scalar T_a,
                                Set::Scalar& q_in, Set::Scalar& q_out,
                                Set::Scalar& q_robin, Set::Scalar& q_loss,
                                Set::Scalar T_pin = 0.0, int* branch = nullptr,
                                Set::Scalar T_surf_in = 0.0,
                                Set::Scalar* T_face_out = nullptr)
    {
        // ⚠⚠ P1f deliverable 0 -- THE FACE TEMPERATURE, OUT.
        //
        // The lateral criterion needs the temperature of the face the Robin flux
        // was evaluated at, and before this there was NO WAY TO GET IT: `T_surf_in`
        // is an INPUT, the face form's `Tf` was a local used for the fluxes and
        // discarded, `SurfaceRobinT` returns the column's VERTICAL TOP temperature,
        // and `side_T_wall` is the gas side. That gap was P1f revision 3's blocker.
        //
        // ⚠ Defaulted to nullptr and appended LAST, so every existing call site is
        // unchanged and the frozen tree must stay 0/2639 -- criterion (a0).
        //
        // ⚠⚠ EACH BRANCH WRITES THE TEMPERATURE ITS OWN q_robin USED, and the CELL
        // FORM WRITES NOTHING. There q_robin = h*(T_gas - Tc) is built on the CELL
        // AVERAGE, which is not a face temperature; writing Tc would hand the
        // criterion a plausible number that is the wrong quantity -- a silent wrong
        // answer, which criterion (a1) exists to catch.
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
        // ⚠ Every existing caller passes nothing, so T_surf_in defaults to 0
        // and `cell`, `face` and `pinned` stay byte-identical. They are this
        // packet's CONTROLS and every frozen study from d2v to d2ad used
        // `pinned`, so that is not a convenience -- it is scored.
        if (in_patch && T_surf_in > 0.0)
        {
            const Set::Scalar Ts    = T_surf_in;
            const Set::Scalar eps_s = (Ts <= T_eps_b) ? eps_lo : eps_hi;
            q_robin = h * (T_gas - Ts);
            q_loss  = eps_s * sig_SB * (Ts * Ts * Ts * Ts - T_a * T_a * T_a * T_a)
                    + h_c * (Ts - T_a);
            q_in    = q_robin - q_loss;
            q_out   = 0.0;
            if (branch) *branch = 2;
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
            if (!(T_pin > 0.0 && Tf > T_pin))
            {
                q_in    = G * (Tf - Tc);
                q_out   = 0.0;
                q_robin = h * (T_gas - Tf);
                q_loss  = q_robin - q_in;
                if (branch) *branch = 0;
                // P1f: the Newton's own root -- the temperature this q_robin used.
                // ⚠ Inferred on a LINEAR half-cell (G = 2k/dz) while the criterion
                // reads an exponential; the two disagree once ell_lat < dz/2, which
                // criterion (k2) scores rather than assumes away.
                if (T_face_out) *T_face_out = Tf;
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
            if (T_face_out) *T_face_out = Ts;   // P1f: Ts == T_pin on this branch
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
    // step. Likewise the C1 pin state (pin_T/pin_t/pin_rho/pin_Cp) is not
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
        // P0b rev5 deliverable 7, ledger copy of the lateral rule. ⚠ It must
        // mirror the kernel EXACTLY -- same gate, same branch, same order --
        // or ledger_err breaks, which is the check that makes the mirror
        // provable rather than asserted.
        const bool        sp_surfform_l = surface_patch.robin_form ==
                                         SurfacePatch::RobinForm::Surface
                                      && surface_patch.side_T_form ==
                                         SurfacePatch::SideTForm::Follow;
        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;
        const Set::Scalar sp_T     = (sp_on && sp_presc) ? surface_patch.T_at(time) : 0.0;

        auto neumaier = [](Set::Scalar& s, Set::Scalar& c, Set::Scalar x) {
            const Set::Scalar t = s + x;
            if (std::abs(s) >= std::abs(x)) c += (s - t) + x;
            else                            c += (x - t) + s;
            s = t;
        };

        Set::Scalar EH_s = 0.0, EH_c = 0.0, EH0_s = 0.0, EH0_c = 0.0;

        Set::Scalar EHsol_s = 0.0, EHsol_c = 0.0;

        Set::Scalar EH0sol_s = 0.0, EH0sol_c = 0.0;
        Set::Scalar P_robin = 0.0, P_loss = 0.0, P_robin_side = 0.0;
        Set::Scalar side_P_topcell = 0.0, side_P_buried = 0.0;
        // ⚠⚠ REVISION 7 DELIVERABLE 2 -- RE-BIN AT THE LIVE GATE.
        // side_P_{bore,pads,wall} split on foot_r_inner/foot_r_outer, and
        // under surface_patch.nozzle_descent = prescribed BOTH ARE ZERO, so
        // every face lands in `wall` and the "100 % of side power is hole
        // wall" reading of revision 5 was an artefact of a dead gate. These
        // two split on side_face_h_r_min, which is the radius the lateral
        // closure itself switches on (40 mm on Meier), so they measure the
        // thing deliverable 7 actually branches on.
        Set::Scalar side_P_rin = 0.0, side_P_rout = 0.0;
        // Per-column side power, for the wall/riser split that no plotfile can
        // give (38 fields, none of them a side flux). Sized only when the key
        // is on.
        const bool sidecsv = surface_patch.side_power_interval > 0.0;
        std::vector<Set::Scalar> sidecol;
        // D2k
        const bool        side_on = sp_on && surface_patch.side_face_flux;
        const Set::Scalar side_f  = surface_patch.side_face_factor;
        // D2r-1: mirrors the kernel's per-region factor. If these two ever
        // disagree the ledger stops reproducing the flux and ledger_err breaks,
        // which is the whole reason the twin exists.
        const Set::Scalar side_f_above = surface_patch.side_face_factor_above;
        // D2x: the same two keys, read identically. ledger_err staying at 1e-15
        // is the proof this mirror is exact.
        const Set::Scalar side_r_min2 = surface_patch.side_face_r_min
                                      * surface_patch.side_face_r_min;
        // D2ab deliverable 2: mirrors the kernel's gate exactly, or ledger_err
        // breaks -- which is what makes the mirror provable, not asserted.
        const Set::Scalar side_h_rmin2 = surface_patch.side_face_h_r_min
                                       * surface_patch.side_face_h_r_min;
        const bool side_area_exact = surface_patch.side_face_area_exact;
        const bool side_step_h_l = surface_patch.side_face_step_height;
        const bool side_staircase_l = surface_patch.side_face_staircase;
        const bool wetted_law_l =
            (surface_patch.side_face_area_law == SurfacePatch::SideAreaLaw::Wetted);
        Set::Scalar delta_sum = 0.0, delta_n = 0.0, delta_fb = 0.0, delta_faces = 0.0;
        Set::Scalar ratio_num = 0.0, ratio_den = 0.0;
        // D2q-2b item 7: cap the coefficient a SIDE face sees, ONLY when
        // per-cell removal is on. -1 = no cap, which is every pre-D2q-2b run.
        // Applied identically in the kernel and in the ledger, or the ledger
        // stops reproducing the kernel's flux.
        const Set::Scalar side_hmax =
            (surface_patch.side_face_flux && spall_per_cell_removal)
            ? surface_patch.side_face_h_max : Set::Scalar(-1.0);
        const int ihi_dom = domain.bigEnd(0), jhi_dom = domain.bigEnd(1);
        Set::Scalar side_rob = 0.0, side_los = 0.0, side_n = 0.0, side_multi = 0.0;
        // P1e diagnostics. sc_min starts at 1.0 so a run with no contributing
        // face reports the identity rather than a sentinel. sc_seen counts
        // DISTINCT columns, which is what side_sc_cols names -- a per-cell count
        // would multiply every column by its riser height.
        Set::Scalar sc_num = 0.0, sc_den = 0.0, sc_min = 1.0, sc_in35 = 0.0;
        std::vector<char> sc_seen;
        if (surface_patch.side_face_staircase && !k_top_col.empty())
            sc_seen.assign(k_top_col.size(), 0);
        // D2aa deliverable 3: the CHARGED lateral area and the side power split
        // by radial band about (x0, y0) -- bore r < foot_r_inner, pads
        // [foot_r_inner, foot_r_outer), wall beyond. Accumulated from exactly
        // the same `Aside` and `(s_in - s_out) * Aside` that feed side_faces
        // and side_P_face, so criterion (d) closes to round-off by
        // construction rather than by agreement.
        Set::Scalar A_bore = 0.0, A_pads = 0.0, A_wall = 0.0;
        Set::Scalar P_bore = 0.0, P_pads = 0.0, P_wall = 0.0;
        Set::Scalar side_ab_P = 0.0, side_ab_n = 0.0;   // D2r-1 witness
        // P2-1 deliverable 6: the distribution of the coefficient ACTUALLY
        // USED, over the faces the law is applied to. ⚠ A fixed-bin histogram,
        // not a gathered sample: quantiles of a per-rank sample are not
        // rank-independent, and this one reduces with a single ReduceLongSum.
        // ⚠ Per BOX, so it MUST be reduced -- the opposite of a per-COLUMN
        // accumulator. Getting that backwards gave exactly 4x on 4 ranks twice
        // in P1g-2b/2c.
        std::vector<long> h_hist;
        // P2-3 (a5): above-plane faces the derived law reached, this rank.
        long n_above_law = 0;
        if (surface_patch.side_face_h_derived)
            h_hist.assign(SIDE_H_BINS + 1, 0);
        // D2x: the realised weight, accumulated ONLY over the faces the
        // derived area actually applies to. Above-plane faces are governed by
        // side_face_factor_above, a DIFFERENT law; counting them here (they are
        // ~79 % of all faces late in a Meier run, at weight 0) would dilute the
        // mean to 0.147 and hide the derived number this diagnostic exists to
        // report. Both are per-face sums, so the ratio is the weighted mean.
        Set::Scalar side_w_sum = 0.0, side_w_n = 0.0;
        std::vector<int> has_surf(static_cast<std::size_t>(ncx) * ncy, 0);
        // rev7 d2: sized here, where ncx/ncy are known, not at the declaration.
        if (sidecsv) sidecol.assign(static_cast<std::size_t>(ncx) * ncy, 0.0);
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
                // P0b (c''): the solid-only sum. See ledger_E_H_solid.
                if (!self_void)
                {
                    neumaier(EHsol_s, EHsol_c, H(i, j, k) * V);
                    if (!ledger_started)
                        neumaier(EH0sol_s, EH0sol_c, H_old(i, j, k) * V);
                }
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
                                sp_pinned ? pin_T_eff[col] : Set::Scalar(0.0),
                                nullptr,
                                SurfaceRobinT(static_cast<std::size_t>(col)));
                P_robin += q_rob * A;
                P_loss  += q_los * A;

                // D2k: the identical side-face term, or ledger_err breaks.
                // Same call, same inputs, same f; the areas are the lateral
                // faces dy*dz and dx*dz instead of the top face dx*dy.
                const Set::Scalar r2_col =
                    (x - surface_patch.x0) * (x - surface_patch.x0)
                  + (y - surface_patch.y0) * (y - surface_patch.y0);
                if (side_on && !sp_presc && sp_on
                    && r2_col <= sp_r2 && r2_col >= side_r_min2)
                {
                    const bool vxm = (i > ilo_dom)     && rem(i-1, j, k) > 0.5;
                    const bool vxp = (i < ihi_dom)     && rem(i+1, j, k) > 0.5;
                    const bool vym = (j > jlo_dom)     && rem(i, j-1, k) > 0.5;
                    const bool vyp = (j < jhi_dom)     && rem(i, j+1, k) > 0.5;
                    const int nfx = (vxm ? 1 : 0) + (vxp ? 1 : 0);
                    const int nfy = (vym ? 1 : 0) + (vyp ? 1 : 0);
                    if (nfx + nfy > 0)
                    {
                        Set::Scalar s_in = 0.0, s_out = 0.0, s_rob = 0.0, s_los = 0.0;
                        Set::Scalar h_side = sp_on ? flame_h[col] : Set::Scalar(0.0);
                        Set::Scalar T_side = sp_on ? flame_Tg[col] : Set::Scalar(0.0);
                        // D2aa: mirrors the kernel exactly, or ledger_err breaks.
                        // D2ab: same radius gate, same resolved temperature.
                        if (surface_patch.side_face_h >= 0.0
                            && r2_col >= side_h_rmin2)
                        { h_side = surface_patch.side_face_h;
                          // P2-3 site 4 of 5 (the LEDGER). Mirrors the kernel
                          // face for face, or ledger_err breaks -- which is what
                          // makes this mirror provable rather than asserted.
                          T_side = SideFaceTgas(static_cast<int>(col), k); }
                        // P2-1 deliverable 4, ROLE 3 of 4 (the LEDGER). Mirrors
                        // the kernel exactly, or ledger_err breaks.
                        Set::Scalar h_face =
                            SideFaceH(h_side, r2_col, static_cast<int>(col), k);
                        if (side_hmax > 0.0 && h_face > side_hmax) h_face = side_hmax;
                        // ⚠⚠ RESTRICTED TO THE LAW'S OWN SCOPE. The first
                        // smoke filled this over EVERY lateral face, including
                        // r < side_face_h_r_min where the typed path uses the
                        // h_expr map -- which spans 167 to 1358 on its own, so
                        // (b2)'s "p95/p05 >= 2x" would have passed with the law
                        // switched OFF. A span test outside the mechanism's
                        // scope measures the configuration, not the mechanism.
                        // ⚠⚠ P2-3: THE SCOPE TEST IS THE LAW'S OWN REACH, NOT
                        // THE TYPED RADIUS. `SideFaceH` applies the derived law
                        // exactly where `wall_corner[col]` is set, which under
                        // `side_face_h_datum = bore` is a DIFFERENT set from
                        // `r2_col >= side_h_rmin2`. Keeping the typed test here
                        // would have left the diagnostic describing a scope the
                        // kernel no longer uses -- the histogram would miss
                        // every bore-wall face the law had just started
                        // heating, which is the one thing (a5) must see.
                        // ⚠ Under the default datum the two tests are
                        // EQUIVALENT: ColumnWallCorner sets wall_corner exactly
                        // where r2_of >= side_face_h_r_min^2 and k_top is valid,
                        // so a `typed` leg is byte-identical.
                        const bool in_law =
                            (col >= 0
                             && col < static_cast<int>(wall_corner.size())
                             && wall_corner[static_cast<std::size_t>(col)]
                                != std::numeric_limits<int>::min());
                        if (!h_hist.empty()
                            && surface_patch.side_face_h >= 0.0
                            && in_law)
                        {
                            int b = static_cast<int>(h_face / SIDE_H_BIN_W);
                            if (b < 0) b = 0;
                            if (b > SIDE_H_BINS) b = SIDE_H_BINS;
                            // ⚠ weighted by the OPEN FACE COUNT, so the
                            // quantiles are per FACE as deliverable 6 states.
                            // Per cell they would under-weight every corner,
                            // which is where the law differs most from typed.
                            h_hist[static_cast<std::size_t>(b)] += (nfx + nfy);
                            // P2-3 (a5): the ABOVE-PLANE share of that reach.
                            // ⚠ the LEDGER's own above-plane form, identical to
                            // `above_c` just below (:4453). The kernel's `fs`
                            // pointer is not in scope in this block, and writing
                            // a second spelling of the same test is how the
                            // kernel and its mirror drift apart.
                            if (sp_on && !(flame_s[col] > 0.0))
                                n_above_law += (nfx + nfy);
                        }
                        // P0b rev5 deliverable 7, mirrored from the kernel.
                        const bool lat_wall_l = sp_surfform_l
                                              && (r2_col >= side_h_rmin2);
                        const bool lat_face_l = sp_surfform_l ? true : sp_face;
                        const Set::Scalar lat_Ts_l = (sp_surfform_l && !lat_wall_l)
                            ? SurfaceRobinT(static_cast<std::size_t>(col))
                            : Set::Scalar(0.0);
                        SurfaceCellFlux(true, lat_face_l, 1.0,
                                        h_face,
                                        T_side,
                                        kc, DX[0], Tc,
                                        eps_rad, eps_rad_hi, T_eps_break, sigma_SB,
                                        h_conv, T_amb, s_in, s_out, s_rob, s_los,
                                        sp_pinned ? pin_T_eff[col] : Set::Scalar(0.0),
                                        nullptr, lat_Ts_l);
                        const bool above_r = (side_f_above >= 0.0)
                                           && sp_on && !(flame_s[col] > 0.0);
                        const bool above_c = sp_on && !(flame_s[col] > 0.0);
                        const Set::Scalar sf_c = above_r ? side_f_above : side_f;
                        // D2x: mirrors the kernel branch exactly, including
                        // which faces it applies to and in what order.
                        Set::Scalar Aside;
                        // D2z-2 deliverable 3: the n >= 2 part of Aside. Zero
                        // on every branch but the step-height one.
                        Set::Scalar A_multi = 0.0;
                        // D2z: mirrors the kernel's step branch EXACTLY --
                        // same faces, same order, same per-cell split. If it
                        // does not, ledger_err blows up, which is precisely the
                        // check that makes the mirror provable rather than
                        // asserted.
                        if (side_step_h_l && !above_r && !k_top_col.empty()
                            && !col_t_last.empty())
                        {
                            const int kT = k_top_col[col];
                            // P1e: riser cell counts, 0 where not open.
                            int rn_l[4] = {0, 0, 0, 0};
                            int fb = 0;
                            const Set::Scalar pT = ColPhi(col_t_last.data(), col_dt.data(),
                                                          col_nvoid.data(), col, time, fb);
                            const int nb[4] = {col - 1, col + 1,
                                               col - ncx, col + ncx};
                            const bool hs[4] = {vxm, vxp, vym, vyp};
                            Set::Scalar wsum = 0.0, dsum = 0.0, wmulti = 0.0;
                            int nfb = 0, nstep = 0, npair = 0;
                            for (int d = 0; d < 4; ++d)
                            {
                                if (!hs[d]) continue;
                                const int kS = k_top_col[nb[d]];
                                const int nn = kT - kS;
                                if (nn < 1) continue;
                                int fbS = 0;
                                const Set::Scalar pS = ColPhi(col_t_last.data(), col_dt.data(),
                                                              col_nvoid.data(), nb[d], time, fbS);
                                const Set::Scalar w = StepCellWeight(nn, k, kT, kS, pT, pS);
                                rn_l[d] = nn;
                                const Set::Scalar aw =
                                    AreaOfDelta(w < 0.0 ? 0.0 : w, wetted_law_l);
                                wsum += aw;
                                // D2z-2 deliverable 3: the share of the area --
                                // and so of side_P_face -- drawn by faces with
                                // n >= 2. Without it a bad SW result cannot be
                                // attributed to the amplification rather than
                                // to anything else.
                                if (nn >= 2) wmulti += aw;
                                ++nstep;                       // stepped FACES (cell x direction)
                                // D2z-2 deliverable 1: delta_mean must be the
                                // PAIR's full step height, not StepCellWeight's
                                // per-cell share -- the two coincide only at
                                // n = 1, and D2z's live column read -87.5 %
                                // against the offline estimator on the 1 mm rim
                                // (26.6 % multi-cell).
                                //
                                // ⚠ Accumulated ONCE PER PAIR, at the taller
                                // column's top exposed cell (k == kT), NOT once
                                // per stepped face. The packet's wording says
                                // "once per stepped face", but this branch runs
                                // for every exposed cell of the riser, so a
                                // per-face accumulation would weight an n-cell
                                // pair n times -- area weighting -- while
                                // delta_offline.py weights every PAIR equally.
                                // Criterion (b) compares the two directly, so
                                // matching its weighting is what "specify both
                                // sides the same way" requires. The per-face
                                // count is not lost: it is the new `delta_faces`
                                // column, and delta_faces / delta_n is the mean
                                // riser height in cells.
                                if (k == kT)
                                {
                                    const Set::Scalar dp = Set::Scalar(nn) + pS - pT;
                                    dsum += (dp > 0.0 ? dp : 0.0);
                                    ++npair;
                                    nfb += (fb | fbS);
                                }
                            }
                            // P1e SITE 3 OF 4 -- THE LEDGER MIRROR. :3980 states
                            // this branch matches the kernel "EXACTLY -- same
                            // faces, same order, same per-cell split. If it does
                            // not, ledger_err blows up." That is precisely why
                            // applying the factor at three sites of four FAILS
                            // LOUDLY instead of silently: criterion (e).
                            // ⚠ A_multi is a SHARE of Aside, so it takes the same
                            // factor or the share stops being a share.
                            if (side_staircase_l)
                            {
                                const Set::Scalar scf = StaircaseFactor(
                                    rn_l[0], rn_l[1], rn_l[2], rn_l[3]);
                                // Diagnostics, area-weighted by the area the
                                // column WOULD have contributed, which is the
                                // weighting the delivered power actually sees.
                                sc_num += scf * wsum;
                                sc_den += wsum;
                                if (scf < sc_min) sc_min = scf;
                                if (col >= 0 && col < static_cast<int>(sc_seen.size())
                                    && !sc_seen[col])
                                { sc_seen[col] = 1; }
                                wsum   *= scf;
                                wmulti *= scf;
                            }
                            Aside = wsum * DX[0] * dz;
                            A_multi = wmulti * DX[0] * dz;
                            delta_sum   += dsum;
                            // ⚠ was (nfx + nfy) -- EVERY exposed lateral face,
                            // including the nn < 1 directions that contribute
                            // nothing to the sum and so biased the mean DOWN.
                            delta_n     += Set::Scalar(npair);
                            delta_faces += Set::Scalar(nstep);
                            delta_fb    += Set::Scalar(nfb);
                            // riser_area_ratio: assigned area over the OTHER
                            // law's, on the same faces, so the two laws can be
                            // compared without running both.
                            Set::Scalar other = 0.0;
                            for (int d = 0; d < 4; ++d)
                            {
                                if (!hs[d]) continue;
                                const int kS = k_top_col[nb[d]];
                                const int nn = kT - kS;
                                if (nn < 1) continue;
                                int fbS = 0;
                                const Set::Scalar pS = ColPhi(col_t_last.data(), col_dt.data(),
                                                              col_nvoid.data(), nb[d], time, fbS);
                                const Set::Scalar w = StepCellWeight(nn, k, kT, kS, pT, pS);
                                other += AreaOfDelta(w < 0.0 ? 0.0 : w, !wetted_law_l);
                            }
                            // P1e SITE 4 OF 4 -- the opposite-law diagnostic. It
                            // exists so the two area laws can be compared without
                            // running both, so it must be corrected the same way
                            // or the comparison acquires a spurious ~0.81 factor.
                            if (side_staircase_l)
                                other *= StaircaseFactor(rn_l[0], rn_l[1],
                                                         rn_l[2], rn_l[3]);
                            ratio_num += wsum;
                            ratio_den += other;
                        }
                        else if (side_area_exact && !above_r && !k_top_col.empty())
                        {
                            const int nc = k_top_col[col];
                            Set::Scalar wsum = 0.0;
                            if (vxm) wsum += RiserWeight(nc - k_top_col[col - 1],   DX[0], dz);
                            if (vxp) wsum += RiserWeight(nc - k_top_col[col + 1],   DX[0], dz);
                            if (vym) wsum += RiserWeight(nc - k_top_col[col - ncx], DX[0], dz);
                            if (vyp) wsum += RiserWeight(nc - k_top_col[col + ncx], DX[0], dz);
                            Aside = wsum * DX[0] * dz;
                            side_w_sum += wsum;
                            side_w_n   += Set::Scalar(nfx + nfy);
                        }
                        else
                        {
                            Aside = sf_c
                                * (Set::Scalar(nfx) * DX[1] * dz + Set::Scalar(nfy) * DX[0] * dz);
                            if (!side_area_exact)
                            {
                                side_w_sum += sf_c * Set::Scalar(nfx + nfy);
                                side_w_n   += Set::Scalar(nfx + nfy);
                            }
                        }
                        P_robin  += s_rob * Aside;
                        P_loss   += s_los * Aside;
                        // P0b deliverable 1.2: side faces have no column
                        // profile, so they keep the CELL AVERAGE. That is a
                        // real approximation, so the share of P_robin they
                        // carry goes on the record rather than being assumed
                        // small -- without it the Meier rerun is not
                        // interpretable.
                        P_robin_side += s_rob * Aside;
                        // P0b rev5 deliverable 9, bin 2: "does the receiving
                        // cell later fire?". The honest online proxy is
                        // whether the cell IS its own column's top cell -- the
                        // one the criterion scores next -- or is buried below
                        // it and cannot fire until the floor recedes past it.
                        // ⚠ Share of power is not share of error: a buried
                        // cell's heat is not wasted, it is stored and arrives
                        // later. This bounds which side power can act on the
                        // RATE, which is the question deliverable 9 asks.
                        if (!k_top_col.empty() && k == k_top_col[col])
                             side_P_topcell += s_rob * Aside;
                        else side_P_buried  += s_rob * Aside;
                        if (r2_col >= side_h_rmin2) side_P_rout += s_rob * Aside;
                        else                        side_P_rin  += s_rob * Aside;
                        // P1e: the fraction of side power landing inside 35 mm.
                        // ⚠ Accumulated as (s_in - s_out)*Aside, the SAME
                        // quantity side_P_face is built from (side_rob), so the
                        // reported number is a true fraction of it rather than of
                        // a different integral. It is the roughness measure that
                        // says whether side heat went to the WALL or, once the
                        // floor roughens, onto the FLOOR -- which is what decides
                        // how to read any rate change.
                        if (surface_patch.side_face_staircase &&
                            r2_col < Set::Scalar(0.035 * 0.035))
                            sc_in35 += (s_in - s_out) * Aside;
                        if (sidecsv && col >= 0 &&
                            col < static_cast<int>(sidecol.size()))
                            sidecol[col] += s_rob * Aside;
                        side_rob += (s_in - s_out) * Aside;
                        side_multi += (s_in - s_out) * A_multi;
                        side_los += s_los * Aside;
                        {
                            const Set::Scalar ri2 = surface_patch.foot_r_inner
                                                  * surface_patch.foot_r_inner;
                            const Set::Scalar ro2 = surface_patch.foot_r_outer
                                                  * surface_patch.foot_r_outer;
                            const Set::Scalar Pf = (s_in - s_out) * Aside;
                            if (r2_col < ri2)       { A_bore += Aside; P_bore += Pf; }
                            else if (r2_col < ro2)  { A_pads += Aside; P_pads += Pf; }
                            else                    { A_wall += Aside; P_wall += Pf; }
                        }
                        side_n   += Set::Scalar(nfx + nfy);
                        // D2r-1 criterion (c): THE MECHANISM'S OWN WITNESS.
                        // If this is 0 on the test leg, the leg heated nothing
                        // and the result is a plumbing failure, not physics.
                        if (above_c)
                        {
                            side_ab_P += (s_in - s_out) * Aside;
                            side_ab_n += Set::Scalar(nfx + nfy);
                        }
                    }
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(EH_s);
        amrex::ParallelDescriptor::ReduceRealSum(EH_c);
        amrex::ParallelDescriptor::ReduceRealSum(EHsol_s);
        amrex::ParallelDescriptor::ReduceRealSum(EHsol_c);
        amrex::ParallelDescriptor::ReduceRealSum(P_robin);
        amrex::ParallelDescriptor::ReduceRealSum(P_loss);
        if (side_on)
        {
            amrex::ParallelDescriptor::ReduceRealSum(side_rob);
            amrex::ParallelDescriptor::ReduceRealSum(side_multi);
            for (Set::Scalar* q : {&A_bore, &A_pads, &A_wall,
                                   &P_bore, &P_pads, &P_wall})
                amrex::ParallelDescriptor::ReduceRealSum(*q);
            side_A_bore = A_bore; side_A_pads = A_pads; side_A_wall = A_wall;
            side_P_bore = P_bore; side_P_pads = P_pads; side_P_wall = P_wall;
            amrex::ParallelDescriptor::ReduceRealSum(side_los);
            amrex::ParallelDescriptor::ReduceRealSum(side_n);
            side_P_face = side_rob;
            // P1e: reduce and publish. ⚠ sc_min takes ReduceRealMin, not Sum --
            // a summed minimum is meaningless and would read ~nranks.
            if (surface_patch.side_face_staircase)
            {
                amrex::ParallelDescriptor::ReduceRealSum(sc_num);
                amrex::ParallelDescriptor::ReduceRealSum(sc_den);
                amrex::ParallelDescriptor::ReduceRealSum(sc_in35);
                amrex::ParallelDescriptor::ReduceRealMin(sc_min);
                long ncol_sc = 0;
                for (char c : sc_seen) ncol_sc += (c ? 1 : 0);
                amrex::ParallelDescriptor::ReduceLongSum(ncol_sc);
                side_sc_mean = (sc_den > 0.0) ? sc_num / sc_den : 1.0;
                side_sc_min  = sc_min;
                side_sc_cols = static_cast<Set::Scalar>(ncol_sc);
                side_P_in35  = (std::abs(side_rob) > 0.0)
                             ? sc_in35 / side_rob : 0.0;
            }
            side_P_face_multi = side_multi;
            side_P_loss = side_los;
            // D2x: the face-count-weighted mean realised weight. It must land
            // in [0.4142, 1] by construction; 0.2 would mean the key is inert.
            if (side_area_exact)
            {
                amrex::ParallelDescriptor::ReduceRealSum(side_w_sum);
                amrex::ParallelDescriptor::ReduceRealSum(side_w_n);
                side_w_mean = (side_w_n > 0.0) ? side_w_sum / side_w_n : 0.0;
            }
            // D2z: delta_mean is a population mean over the LIVE face set and
            // is therefore TIME-WEIGHTED, not at-birth. ⚠ Masked on its own
            // population (below-plane step faces only): D2x's side_w_mean read
            // 0.147 because 79 % of late-run faces are above-plane under a
            // different law and were counted at weight zero.
            if (side_step_h_l)
            {
                amrex::ParallelDescriptor::ReduceRealSum(delta_sum);
                amrex::ParallelDescriptor::ReduceRealSum(delta_n);
                amrex::ParallelDescriptor::ReduceRealSum(delta_fb);
                amrex::ParallelDescriptor::ReduceRealSum(delta_faces);
                amrex::ParallelDescriptor::ReduceRealSum(ratio_num);
                amrex::ParallelDescriptor::ReduceRealSum(ratio_den);
                delta_mean          = (delta_n > 0.0) ? delta_sum / delta_n : 0.0;
                delta_fallback_cols = delta_fb;
                step_delta_faces    = delta_faces;
                riser_area_ratio    = (ratio_den > 0.0) ? ratio_num / ratio_den : 0.0;
            }
            side_faces  = side_n;
            if (!h_hist.empty())
            {
                amrex::ParallelDescriptor::ReduceLongSum(
                    h_hist.data(), static_cast<int>(h_hist.size()));
                long tot = 0;
                for (long v : h_hist) tot += v;
                side_h_used_n = Set::Scalar(tot);
                // P2-3 (a5): per-rank, each rank owning distinct faces, so a
                // sum is correct -- exactly like h_hist above.
                amrex::ParallelDescriptor::ReduceLongSum(n_above_law);
                side_h_used_n_above = Set::Scalar(n_above_law);
                auto qtl = [&](Set::Scalar q) -> Set::Scalar {
                    if (tot <= 0) return Set::Scalar(0.0);
                    const long want = static_cast<long>(q * Set::Scalar(tot));
                    long run = 0;
                    for (std::size_t b = 0; b < h_hist.size(); ++b)
                    {
                        run += h_hist[b];
                        if (run > want)
                            return (Set::Scalar(b) + 0.5) * SIDE_H_BIN_W;
                    }
                    return (Set::Scalar(h_hist.size() - 1) + 0.5) * SIDE_H_BIN_W;
                };
                side_h_used_min = qtl(0.0);
                side_h_used_p05 = qtl(0.05);
                side_h_used_p50 = qtl(0.50);
                side_h_used_p95 = qtl(0.95);
                side_h_used_max = qtl(0.999999);
            }
            if (surface_patch.side_face_factor_above >= 0.0)
            {
                amrex::ParallelDescriptor::ReduceRealSum(side_ab_P);
                amrex::ParallelDescriptor::ReduceRealSum(side_ab_n);
                side_P_face_above = side_ab_P;
                side_cols_above   = side_ab_n;
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(has_surf.data(), static_cast<int>(has_surf.size()));
        if (!ledger_started)
        {
            amrex::ParallelDescriptor::ReduceRealSum(EH0_s);
            amrex::ParallelDescriptor::ReduceRealSum(EH0_c);
            ledger_E_H0_s  = EH0_s;
            ledger_E_H0_c  = EH0_c;
            amrex::ParallelDescriptor::ReduceRealSum(EH0sol_s);
            amrex::ParallelDescriptor::ReduceRealSum(EH0sol_c);
            ledger_E_H0_sol_s = EH0sol_s;
            ledger_E_H0_sol_c = EH0sol_c;
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
        // P0b witnesses for robin_form = surface.
        if (surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
        {
            amrex::ParallelDescriptor::ReduceRealSum(P_robin_side);
            amrex::ParallelDescriptor::ReduceRealSum(side_P_topcell);
            amrex::ParallelDescriptor::ReduceRealSum(side_P_buried);
            amrex::ParallelDescriptor::ReduceRealSum(side_P_rin);
            amrex::ParallelDescriptor::ReduceRealSum(side_P_rout);
            surf_side_share = (std::abs(P_robin) > 0.0)
                              ? P_robin_side / P_robin : 0.0;
            surf_side_P_top    = side_P_topcell;
            surf_side_P_buried = side_P_buried;
            surf_side_P_rin    = side_P_rin;
            surf_side_P_rout   = side_P_rout;
            // rev7 d2: the per-column side-power log. Deliverable 2 is explicit
            // that the wall/riser split "cannot be a re-analysis" -- the
            // plotfiles carry 38 fields and not one of them is a side flux, and
            // the thermo holds only aggregates. One row per column that carries
            // side power, at side_power_interval seconds of simulated time.
            if (sidecsv)
            {
                amrex::ParallelDescriptor::ReduceRealSum(sidecol.data(),
                                                         static_cast<int>(sidecol.size()));
                const long widx = static_cast<long>(
                    std::floor(time / surface_patch.side_power_interval));
                if (widx > side_power_last)
                {
                    side_power_last = widx;
                    if (amrex::ParallelDescriptor::IOProcessor())
                    {
                        const std::string path = plot_file + "_side_power.csv";
                        const bool first = !side_power_header_written;
                        std::ofstream ofs(path,
                            (first ? std::ios::trunc : std::ios::app) | std::ios::out);
                        if (ofs.is_open())
                        {
                            if (first) ofs << "time,col_i,col_j,r,P_side\n";
                            ofs.setf(std::ios::scientific);
                            ofs.precision(10);
                            for (int cj = 0; cj < ncy; ++cj)
                            for (int ci = 0; ci < ncx; ++ci)
                            {
                                const std::size_t c =
                                    static_cast<std::size_t>(cj) * ncx + ci;
                                if (!(std::abs(sidecol[c]) > 0.0)) continue;
                                const Set::Scalar cx = PLO[0]
                                    + (Set::Scalar(ci + ilo_dom) + 0.5) * DX[0]
                                    - surface_patch.x0;
                                const Set::Scalar cy = PLO[1]
                                    + (Set::Scalar(cj + jlo_dom) + 0.5) * DX[1]
                                    - surface_patch.y0;
                                ofs << time << "," << (ci + ilo_dom) << ","
                                    << (cj + jlo_dom) << ","
                                    << std::sqrt(cx * cx + cy * cy) << ","
                                    << sidecol[c] << "\n";
                            }
                        }
                        side_power_header_written = true;
                    }
                }
            }
            Set::Scalar nsurf = 0.0, qsum = 0.0;
            for (std::size_t c = 0; c < surf_T_col.size(); ++c)
                if (surf_T_col[c] > 0.0)
                { nsurf += 1.0; qsum += flame_h[c] * (flame_Tg[c] - surf_T_col[c]); }
            surf_cols_robin = nsurf;
            surf_q_mean     = (nsurf > 0.0) ? qsum / nsurf : 0.0;
        }
        ledger_P_loss  = P_loss;
        ledger_E_src  += dt * (P_beam + P_robin - P_loss);
        ledger_E_H     = EH_s + EH_c;
        const Set::Scalar dE = (EH_s - ledger_E_H0_s) + (EH_c - ledger_E_H0_c);
        ledger_err = (dE - ledger_E_src) / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));
        // P0b (c''): with the flake sink, the SOLID enthalpy must close against
        // the source minus what left with the rock.
        ledger_E_H_solid = EHsol_s + EHsol_c;
        const Set::Scalar dE_sol = (EHsol_s - ledger_E_H0_sol_s)
                                 + (EHsol_c - ledger_E_H0_sol_c);
        // P2-2f: foot_mech_E is the pad-clearance sink -- enthalpy that left the
        // SOLID set mechanically. Without it this row cannot close on any leg
        // with feet + foot_body_clearance, which is one cause and not two: the
        // gap between ledger_err and ledger_err_solid is the flake credit the
        // TOTAL row omits by design, and the accepted CTRL5 shows a larger gap
        // while closing here to 7.78e-04.
        ledger_err_solid = (dE_sol - (ledger_E_src - ledger_E_flake - foot_mech_E))
                           / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));
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

    // AMR: the finest level owns the surface. Every surface closure (patch,
    // losses, beam, Sp, removal, the per-column state and the CSV logs) runs
    // on this level only; coarser levels advance conduction in their covered
    // cells and are overwritten by average-down. Single level: always true.
    bool SurfaceOwner(int lev) const { return max_level == 0 || lev == finest_level; }

    // AMR: rebuild the column-linear level set on a remade level (lev > 0).
    // Removal.H keeps solid <=> phi >= 0 <=> removed <= 0.5 with phi linear
    // down each column, phi(k) = phi(k_top) + (k_top - k)*dz. RemakeLevel
    // copies cells the level already owned exactly and interpolates cells
    // new to the level from the coarse level (cell_cons_interp), whose
    // lateral slope limiter smears phi across neighbouring columns of
    // different height by up to half a step. The column top is always an
    // already-owned cell (the finest level owns the surface), so it anchors
    // the exact profile. Only cells off the profile by more than round-off
    // are rewritten, so an unchanged grid stays bit-identical.
    void RepairRemovalStateAfterRegrid(int lev)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1), klo = domain.smallEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const int ncols = nx * ny;
        const Set::Scalar dz  = geom[lev].CellSize()[2];
        const Set::Scalar tol = 1.0e-6 * dz;
        auto col = [=](int i, int j) { return (j - jlo) * nx + (i - ilo); };

        amrex::Vector<int> k_top(ncols, std::numeric_limits<int>::min());
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& phi = phi_mf[lev]->const_array(mfi);
            amrex::Array4<const Set::Scalar> const& rem = removed_mf[lev]->const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int& kt = k_top[col(i, j)];
                for (int k = hi.z; k >= lo.z && k > kt; --k)
                    // D2q-1: removed_mf is authoritative. The old `&& phi >= 0`
                    // conjunct was redundant (the two agree cell-for-cell while
                    // void is contiguous) and would mask an interior void.
                    if (rem(i, j, k) <= 0.5) { kt = k; break; }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(k_top.data(), ncols);

        // phi extrapolated to the domain floor index klo: owner-only, summed.
        amrex::Vector<Set::Scalar> phi_floor(ncols, 0.0);
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& phi = phi_mf[lev]->const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const int c  = col(i, j);
                const int kt = k_top[c];
                if (kt == std::numeric_limits<int>::min() || kt < lo.z || kt > hi.z) continue;
                phi_floor[c] = phi(i, j, kt) + Set::Scalar(kt - klo) * dz;
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(phi_floor.data(), ncols);

        int n_fixed = 0;
        for (amrex::MFIter mfi(*phi_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar> const& phi = phi_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& rem = removed_mf[lev]->array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const int c = col(i, j);
                if (k_top[c] == std::numeric_limits<int>::min()) continue;
                const Set::Scalar base = phi_floor[c];
                for (int k = lo.z; k <= hi.z; ++k)
                {
                    // D2q-1: removed_mf is authoritative and is NEVER re-derived
                    // from phi here. Interpolation across a regrid can leave the
                    // 0/1 indicator fractional, so rem is only *cleaned* to an
                    // exact 0/1 by the same > 0.5 test every reader already
                    // applies; phi is re-imposed as the column-linear sub-cell
                    // accumulator, anchored on rem's own top solid cell (k_top
                    // above). While void is contiguous `expect < 0` and
                    // `rem > 0.5` agree cell-for-cell, so this is byte-identical;
                    // once stage 2 creates an interior void, the old form would
                    // have silently filled it in.
                    const Set::Scalar expect = base - Set::Scalar(k - klo) * dz;
                    const Set::Scalar r_auth = (rem(i, j, k) > 0.5) ? 1.0 : 0.0;
                    if (std::abs(phi(i, j, k) - expect) <= tol &&
                        rem(i, j, k) == r_auth) continue;
                    phi(i, j, k) = expect;
                    rem(i, j, k) = r_auth;
                    ++n_fixed;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceIntSum(n_fixed);
        if (n_fixed > 0)
            Util::Message(INFO, "AMR regrid repair on level ", lev, ": rebuilt phi/removed in ",
                          n_fixed, " cells");
    }

    // AMR: after the owner's removal pass, verify the finest level still holds
    // every column's top solid cell and the depth-scan window beneath it
    // (Removal.H reads T(z) from this level only; cells outside it read as
    // 0 K and silently stop the scan). Aborts with the remedy otherwise.
    void CheckSurfaceCovered(int lev, Set::Scalar time)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1), klo = domain.smallEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const int ncols = nx * ny;
        const Set::Scalar dz = geom[lev].CellSize()[2];
        const amrex::Vector<int> k_top = ColumnTopSolid(lev);
        amrex::Vector<int> k_bot(ncols, std::numeric_limits<int>::max());
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int& kb = k_bot[(j - jlo) * nx + (i - ilo)];
                if (lo.z < kb) kb = lo.z;
            }
        }
        amrex::ParallelDescriptor::ReduceIntMin(k_bot.data(), ncols);
        const int n_scan = std::max(1, static_cast<int>(std::floor(spall_h_spall_max / dz)));
        int n_lost = 0, n_shallow = 0;
        for (int c = 0; c < ncols; ++c)
        {
            if (k_top[c] == std::numeric_limits<int>::min()) { ++n_lost; continue; }
            const int k_need = std::max(klo, k_top[c] - n_scan);
            if (k_need < k_bot[c]) ++n_shallow;
        }
        if (n_lost > 0 || n_shallow > 0)
            Util::Abort(INFO, "AMR: the surface left the refined band at t = ", time, " s: ",
                        n_lost, " columns have no solid cell on level ", lev, " and ",
                        n_shallow, " columns have fewer than ", n_scan,
                        " refined cells below their surface. Increase surface.refine_depth "
                        "or lower amr.regrid_int.");
    }

    // E1: global topmost solid (removed <= 0.5) cell index per (i, j)
    // column at level lev, INT_MIN for fully removed columns. Index is
    // (j - dlo_y)*nx + (i - dlo_x). MPI-reduced, so identical on every rank
    // and independent of the z box decomposition. D2q-1: every column-top
    // scan in the code now uses this rule -- `removed` alone -- so Removal.H
    // PASS 1 no longer additionally requires phi >= 0.
    // =====================================================================
    // P2-1: THE COMPOSITE WALL LAW
    // =====================================================================
    // The hole wall was heated through ONE hand-typed scalar applied flat at
    // every depth. These four members replace it with a law, under
    // `surface_patch.side_face_h_law = derived` (default `typed`, so every
    // frozen anchor is untouched).
    //
    //      h_wall(z) = h_wj(z)   for z <  z_fill     (wall jet, thin layer)
    //                = h_duct    for z >= z_fill     (fully-developed annulus)
    //
    // The two earlier derivations are not rivals, they are the same flow at two
    // heights: at the corner the gas has just turned off the floor and runs as
    // a thin fast layer on a free stream; once the layer fills the wall-to-
    // shield gap there is no free stream left to entrain and the flow IS the
    // duct.
    //
    // ⚠⚠ A HARD SWITCH AT z_fill, NEVER max(h_wj, h_duct). Under the `max` the
    // duct is never selected -- h_wj is still ~96 W/m^2K at z = 400 mm against
    // the duct's ~58 -- so the crossover would do nothing and the deep wall
    // would keep a model that no longer applies.
    //
    // ⚠ The step at z_fill is a factor ~2.6 and is ACCEPTED, not blended. It is
    // safe because the datum is the LOCAL CORNER FLOOR, so z_fill is a height
    // above a floor that descends through the run: the step sweeps through the
    // rock rather than sitting at one elevation, and cannot etch a ledge. A
    // blend would need a second length that nothing measures.
    //
    // ⚠⚠ THE HALF IN c_f/2 IS SETTLED BY MEASUREMENT, NOT BY READING A
    // DEFINITION. A six-laboratory compilation of plane wall-jet data reports
    // the friction velocity directly, which is convention-free:
    // u_tau/U_max = 0.050-0.064. The friction law below reproduces that band
    // only with the half; the other reading requires 0.083, which no
    // measurement shows.
    //
    // Offline twin, and the derivation of every constant:
    //   tests/MMWSpalling/studies/p2_1_wallh/wall_law.py  (+ WALL_H.txt)
    static constexpr Set::Scalar WJ_VEL   = 1.05;   // v_max = C u D_e / r
    static constexpr Set::Scalar WJ_DELTA = 0.09;   // delta = C r at the corner
    static constexpr Set::Scalar WJ_GROW  = 0.073;  // d delta / dz up the wall
    static constexpr Set::Scalar WJ_CF_A  = 0.0315; // c_f = A Re_delta^B
    static constexpr Set::Scalar WJ_CF_B  = -0.182;
    static constexpr Set::Scalar WJ_PR    = 0.72;   // gas Prandtl number
    static constexpr Set::Scalar WJ_P_ATM = 101325.0;
    static constexpr Set::Scalar WJ_R_GAS = 287.0;
    // ⚠ The burner's water-cooled shield is Oe67 mm => 33.5 mm radius, from the
    // VT5 burner drawings. This is NOT surface_patch.side_face_h_r_min (40 mm,
    // the Oe80 feet): that gates WHICH faces take wall treatment, this is what
    // the annulus GAP is measured against. It is a declared constant rather
    // than a key because the step is allowed to add exactly ONE key.
    static constexpr Set::Scalar WJ_R_SHIELD = 0.0335;
    // deliverable 6's histogram: 5 W/m^2K bins to 1500, with the last bin an
    // overflow catch-all. The law's whole range on the Meier wall is 55-330, so
    // the resolution is ~2 % of the quantity and nothing saturates.
    static constexpr int         SIDE_H_BINS  = 300;
    static constexpr Set::Scalar SIDE_H_BIN_W = 5.0;

    // Sutherland's law for air (White, Viscous Fluid Flow, 2nd ed. Eq. 1-36).
    static Set::Scalar WallMu(Set::Scalar T)
    {
        const Set::Scalar mu0 = 1.716e-5, T0 = 273.15, S = 110.4;
        return mu0 * std::pow(T / T0, 1.5) * (T0 + S) / (T + S);
    }

    // The height at which the layer fills the wall-to-shield gap [m]. Returns 0
    // when the gap is already thinner than the corner layer, i.e. no wall-jet
    // zone at all.
    static Set::Scalar WallZFill(Set::Scalar r)
    {
        const Set::Scalar gap = r - WJ_R_SHIELD;
        const Set::Scalar d_c = WJ_DELTA * r;
        return (gap <= d_c) ? Set::Scalar(0.0) : (gap - d_c) / WJ_GROW;
    }

    // The wall-jet branch at height z above the corner, at radius r.
    Set::Scalar WallJetBranch(Set::Scalar z, Set::Scalar r) const
    {
        if (!(r > 0.0) || !(wall_nu > 0.0)) return Set::Scalar(0.0);
        const Set::Scalar v_c = WJ_VEL * wall_uDe / r;
        const Set::Scalar d_c = WJ_DELTA * r;
        const Set::Scalar d_z = d_c + WJ_GROW * (z > 0.0 ? z : Set::Scalar(0.0));
        // plane-jet momentum conservation: v^2 * delta = const
        const Set::Scalar v = v_c * std::sqrt(d_c / d_z);
        const Set::Scalar Re = v * d_z / wall_nu;
        if (!(Re > 0.0)) return Set::Scalar(0.0);
        const Set::Scalar cf = WJ_CF_A * std::pow(Re, WJ_CF_B);
        const Set::Scalar St = (0.5 * cf) / std::pow(WJ_PR, 2.0 / 3.0);
        return St * wall_rho * wall_cp * v;
    }

    // The composite law [W/m^2K].
    Set::Scalar WallJetH(Set::Scalar z, Set::Scalar r) const
    {
        if (z >= WallZFill(r)) return wall_h_duct;
        return WallJetBranch(z, r);
    }

    // =====================================================================
    // THE ONE PER-FACE RESOLVER -- every lateral role calls exactly this.
    // =====================================================================
    // ⚠⚠ Removal.H:589-594 pre-registers this bug family: three times now the
    // criterion has read something the flux never used (Pa vs MPa, the Newton
    // length d_n vs DX[0], the coefficient cap). Making the per-height
    // resolution a single function that all four roles call is how the fourth
    // is prevented rather than found by a campaign.
    //
    // With the law off this returns `h_typed` unchanged and bit-for-bit, which
    // is what keeps every frozen anchor byte-identical.
    Set::Scalar SideFaceH(Set::Scalar h_typed, Set::Scalar r2_col,
                          int col, int k) const
    {
        if (!surface_patch.side_face_h_derived || !wall_law_ready)
            return h_typed;
        if (col < 0 || col >= static_cast<int>(wall_corner.size()))
            return h_typed;
        const int kc = wall_corner[static_cast<std::size_t>(col)];
        if (kc == std::numeric_limits<int>::min()) return h_typed;
        Set::Scalar z = (Set::Scalar(k - kc) - 0.5) * wall_dz;
        if (z < 0.5 * wall_dz) z = 0.5 * wall_dz;
        const Set::Scalar h = WallJetH(z, std::sqrt(r2_col));
        // P2-3 deliverable 5: the scale multiplies the DERIVED law's output
        // HERE AND NOWHERE ELSE, so a leg running it cannot also be scaling a
        // typed fallback by accident. 1.0 => byte-identical.
        return (h > 0.0) ? h * surface_patch.side_face_h_scale : h_typed;
    }

    // =====================================================================
    // P2-3 deliverable 3 -- ONE RESOLVER FOR THE WALL-GAS TEMPERATURE.
    // =====================================================================
    // `side_T_eff` is ONE SCALAR for every side face at every height: under
    // `side_face_T = live` it is `jet_T_exhaust`, 1,504 K on the P2-2f m13 leg,
    // seen identically by a face 10 mm down the bore and one 400 mm down. That
    // uniform hot gas is the mechanism that ballooned D2v's M2 mouth to Ø 137
    // against Meier's 96.
    //
    // With the wall march on, an ABOVE-PLANE face instead sees the gas
    // temperature at ITS OWN height, which the wall itself has depleted on the
    // way up. Below-plane faces are unchanged: they are floor cone flank, the
    // floor march already owns their gas, and the wall march starts from the
    // floor march's exit state.
    //
    // ⚠⚠ THERE ARE FIVE READERS AND THE CRITERION IS THE ONE THAT GETS MISSED
    // (P2-0 lesson): Removal.H:605 (lateral criterion), :2047 (kernel), :3419
    // (march), :4374 (ledger), :5350 (wall-law gas properties). All five go
    // through here, and `unit/lateral_flux_lock` asserts it.
    // ⚠ Returns `side_T_eff` unchanged whenever the march is off, the face is
    // below plane, or the bin is missing => byte-identical with the key off.
    Set::Scalar SideFaceTgas(int col, int k) const
    {
        if (!surface_patch.jet_wall_march || wall_T_gas.empty()) return side_T_eff;
        if (col < 0 || col >= static_cast<int>(flame_s.size())) return side_T_eff;
        // above plane == flame_s <= 0, the test used at :2250 / :2329 / :8160.
        if (flame_s[static_cast<std::size_t>(col)] > 0.0) return side_T_eff;
        const int b = k - wall_k0;
        if (b < 0 || b >= static_cast<int>(wall_T_gas.size())) return side_T_eff;
        const Set::Scalar T = wall_T_gas[static_cast<std::size_t>(b)];
        return (T > 0.0) ? T : side_T_eff;
    }

    // =====================================================================
    // THE HEIGHT DATUM -- ONE definition, built once per step.
    // =====================================================================
    //   k_corner[col] = the DEEPEST column top on the straight ray from this
    //   column inward to the axis, truncated at the wall band's own inner edge
    //   r = side_face_h_r_min.   z = (k - k_corner - 1/2) dz, floored at dz/2.
    //
    // ⚠ NOT the nearest open neighbour's top. That reads the one-cell
    // STAIRCASE tread, not the hole's floor: on the Meier wall it gives a
    // median height of 7 mm where the crater floor profile puts the band at
    // 30-65 mm. The wall-jet layer is 4.4 mm at the corner and ~15 mm by
    // z = 150 mm, so it spans several treads -- a 2 mm ledge every 10 mm is
    // ROUGHNESS and does not restart the jet. Resetting z at each tread would
    // also make the datum mesh-dependent, which is the D2j-0b failure mode.
    //
    // ⚠ NOT the axis front either (that is the p50 90 mm prior packets quote):
    // the jet turns at the corner under THIS azimuth, and the crater floor
    // rises ~50 mm from the axis to the band.
    //
    // ⚠ r_gate is reused deliberately: inside it the code is already in regime
    // 7b (floor riser), so it is the EXISTING floor/wall boundary and the datum
    // introduces no new radius.
    //
    // ⚠⚠ P2-3: THE TRUNCATION RADIUS IS NOW A PARAMETER, not a read of
    // `side_face_h_r_min`. Under the default (`side_face_h_datum = typed`) the
    // caller passes exactly that key and nothing changes. Under `bore` it
    // passes the MEASURED bore radius less one cell, which is the only reading
    // under which the datum's inner edge is the hole's actual floor/wall corner
    // rather than a typed guess that the bore had already drifted 0.23 mm
    // inside of. See the key's declaration for why the typed number could not
    // simply be lowered instead.
    //
    // Decomposition-independent by construction: it is a pure function of
    // k_top, which ColumnTopSolid has already ReduceIntMax'd.
    amrex::Vector<int> ColumnWallCorner(int lev,
                                        const amrex::Vector<int>& k_top,
                                        Set::Scalar r_trunc) const
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int nx = domain.length(0), ny = domain.length(1);
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar x0 = surface_patch.x0, y0 = surface_patch.y0;
        const Set::Scalar rg2 = r_trunc * r_trunc;
        const int MISS = std::numeric_limits<int>::min();
        amrex::Vector<int> out(k_top.size(), MISS);
        if (k_top.empty()) return out;
        auto r2_of = [&](int i, int j) {
            const Set::Scalar x = PLO[0] + (Set::Scalar(i + ilo) + 0.5) * DX[0] - x0;
            const Set::Scalar y = PLO[1] + (Set::Scalar(j + jlo) + 0.5) * DX[1] - y0;
            return x * x + y * y;
        };
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(j) * nx + i;
            if (k_top[c] == MISS || r2_of(i, j) < rg2) continue;
            int best = k_top[c];
            // march inward along the ray to the axis, one cell at a time
            const int n = (i > j ? i : j);
            for (int s = 1; s <= n; ++s)
            {
                const Set::Scalar t = 1.0 - Set::Scalar(s) / Set::Scalar(n + 1);
                const int ii = static_cast<int>(
                    std::lround((Set::Scalar(i) + 0.5) * t - 0.5));
                const int jj = static_cast<int>(
                    std::lround((Set::Scalar(j) + 0.5) * t - 0.5));
                if (ii < 0 || jj < 0 || ii >= nx || jj >= ny) break;
                if (r2_of(ii, jj) < rg2) break;
                const int kt = k_top[static_cast<std::size_t>(jj) * nx + ii];
                if (kt != MISS && kt < best) best = kt;
            }
            out[c] = best;
        }
        return out;
    }

    // The azimuthal fraction of a full circle of radius r about the patch axis
    // that lies INSIDE the domain. The Meier input is a quarter annulus with
    // the axis at the domain corner, so this is exactly 1/4; a centred patch
    // gives 1. Derived from geometry already in scope, so the step adds no key
    // for it, and logged, so the quarter-domain assumption is an output rather
    // than a comment.
    Set::Scalar PatchAzimuthFraction(int lev, Set::Scalar r) const
    {
        const amrex::Box& domain = geom[lev].Domain();
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar xlo = PLO[0], ylo = PLO[1];
        const Set::Scalar xhi = PLO[0] + Set::Scalar(domain.length(0)) * DX[0];
        const Set::Scalar yhi = PLO[1] + Set::Scalar(domain.length(1)) * DX[1];
        if (!(r > 0.0)) return Set::Scalar(0.0);
        const int N = 2048;
        int in = 0;
        for (int a = 0; a < N; ++a)
        {
            const Set::Scalar th = 2.0 * 3.14159265358979323846
                                 * (Set::Scalar(a) + 0.5) / Set::Scalar(N);
            const Set::Scalar x = surface_patch.x0 + r * std::cos(th);
            const Set::Scalar y = surface_patch.y0 + r * std::sin(th);
            if (x >= xlo && x <= xhi && y >= ylo && y <= yhi) ++in;
        }
        // ⚠⚠ 0, NOT 1, when the probe circle lies wholly outside the domain.
        // Returning 1 is what the first smoke run did, silently: the probe used
        // surface_patch.radius = 0.2 m against a 0.12 m domain, every sample
        // fell outside, and the fallback handed back a FULL circle -- so
        // mdot_full came out 4x too small and the whole law with it. A
        // degenerate probe is a failure to be reported, not a default.
        return (in > 0) ? Set::Scalar(in) / Set::Scalar(N) : Set::Scalar(0.0);
    }

    // =====================================================================
    // The per-step law state. MUST run before the criterion, i.e. first in
    // Advance: the criterion is the earliest of the four roles.
    // =====================================================================
    void UpdateWallLawState(int lev)
    {
        wall_law_ready = false;
        wall_law_on = 0.0;
        if (!surface_patch.side_face_h_derived) { wall_corner.clear(); return; }
        if (!SurfaceOwner(lev) || !removal_enabled()) { wall_corner.clear(); return; }

        const amrex::Vector<int> k_top = ColumnTopSolid(lev);
        if (k_top.empty()) { wall_corner.clear(); return; }
        // ⚠⚠ P2-3: the datum is built BELOW, after r_bore is measured, because
        // under `side_face_h_datum = bore` it needs that radius. Both use the
        // SAME k_top, so moving the call introduces no lag and no second probe
        // -- it is a reordering inside one function, not a staler input.
        wall_corner.clear();

        const amrex::Box& domain = geom[lev].Domain();
        const Set::Scalar* DX = geom[lev].CellSize();
        const int nz = domain.length(2);
        const int klo = domain.smallEnd(2);
        const int nx_dom = domain.length(0), ny_dom = domain.length(1);
        wall_dz = DX[2];

        // THE BORE, from the void mask -- never from the event log, whose
        // estimator is biased ~0.3 mm (P1g-2c carried risk iv). Area-equivalent
        // rather than max-radius: a staircased wall makes the max radius a
        // one-cell artefact. Taken over the live wall band only, so the collar
        // near the mouth cannot stand in for the bore.
        // ⚠ One pass over columns into a histogram, then a prefix sum: the void
        // count at level k is the number of columns whose top lies below k. The
        // naive nz*ncols double loop is 1.1M compares EVERY step.
        const int MISS = std::numeric_limits<int>::min();
        std::vector<long> below(static_cast<std::size_t>(nz) + 1, 0);
        for (std::size_t c = 0; c < k_top.size(); ++c)
        {
            const int kt = k_top[c];
            if (kt == MISS) continue;
            const int kk = kt + 1 - klo;                 // first void level
            if (kk <= 0) below[0] += 1;
            else if (kk < nz) below[static_cast<std::size_t>(kk)] += 1;
        }
        for (int k = 1; k < nz; ++k) below[static_cast<std::size_t>(k)]
            += below[static_cast<std::size_t>(k - 1)];
        // the live band: drop 20 mm above the advancing front and 100 mm below
        // the mouth. ⚠ a band within 20 mm of either front is not comparable,
        // and the collar near the mouth runs ~10 mm wider than the bore, so
        // including it would let the collar stand in for the hole.
        int k_deep = nz;                   // deepest level carrying any void
        for (int k = 0; k < nz; ++k)
            if (below[static_cast<std::size_t>(k)] > 0) { k_deep = k; break; }
        // ⚠ The band must be defined on the hole that EXISTS, not on a fixed
        // depth. The first smoke run drilled 15 mm in 40 s, the 100 mm mouth pad
        // then emptied the band, and the area-equivalent bore collapsed onto the
        // WJ_R_SHIELD clamp -- a 0.3 mm gap, h_duct 770 and z_fill 0, i.e. every
        // face on a duct branch that does not exist. Pads are applied only when
        // the hole is deep enough to carry them; otherwise the middle half of
        // the void range is used.
        const int pad_lo = static_cast<int>(0.020 / DX[2] + 0.5);
        const int pad_hi = static_cast<int>(0.100 / DX[2] + 0.5);
        int b_lo = k_deep, b_hi = nz - 1;
        if ((b_hi - b_lo) > 2 * (pad_lo + pad_hi))
        { b_lo = k_deep + pad_lo; b_hi = (nz - 1) - pad_hi; }
        else if (b_hi > b_lo)
        { const int q = (b_hi - b_lo) / 4; b_lo += q; b_hi -= q; }
        long nvoid = 0; int nk = 0;
        for (int k = b_lo; k <= b_hi && k < nz; ++k)
        {
            const long nv = below[static_cast<std::size_t>(k)];
            if (nv > 0) { nvoid += nv; ++nk; }
        }
        Set::Scalar r_bore = 0.0;
        if (nk > 0)
        {
            // ⚠ Probed at a quarter of the smaller domain extent, NOT at
            // surface_patch.radius: the patch radius may exceed the domain (the
            // Meier input sets 0.2 m on a 0.12 m box). With the axis on a
            // domain corner the azimuth fraction is identical for every radius
            // that fits, so this is exact rather than a sample of a guess.
            const Set::Scalar Lx = Set::Scalar(nx_dom) * DX[0];
            const Set::Scalar Ly = Set::Scalar(ny_dom) * DX[1];
            const Set::Scalar f = PatchAzimuthFraction(
                lev, 0.25 * ((Lx < Ly) ? Lx : Ly));
            wall_sector = f;
            const Set::Scalar A = (Set::Scalar(nvoid) / Set::Scalar(nk))
                                * DX[0] * DX[1] / (f > 0.0 ? f : Set::Scalar(1.0));
            r_bore = std::sqrt(A / 3.14159265358979323846);
        }
        wall_r_bore = r_bore;
        // P2-3 deliverable 4: the height datum, now that the bore is known.
        // ⚠ one cell of slack: the bore-wall COLUMN CENTRES sit just inside the
        // area-equivalent radius (39.65 vs 39.77 mm measured on m13), so a
        // threshold at r_bore exactly would still skip them.
        const Set::Scalar r_trunc =
            (surface_patch.side_face_h_datum_bore && r_bore > 0.0)
                ? std::max(r_bore - DX[0], Set::Scalar(0.0))
                : surface_patch.side_face_h_r_min;
        wall_datum_r = r_trunc;
        wall_corner = ColumnWallCorner(lev, k_top, r_trunc);
        // ⚠⚠ NO ANNULUS YET. Until the bore clears the burner shield by at
        // least two cells there is no gap to run a duct through and no corner
        // for a wall jet to turn at, so the law is NOT applied and the typed
        // coefficient stands. Reported through wall_law_on, so the inert period
        // is an OUTPUT -- an inert law that nothing shows is the failure mode
        // three keys in this line of packets already hit.
        if (!(r_bore > WJ_R_SHIELD + 2.0 * DX[0]))
        {
            wall_law_on = 0.0;
            wall_h_duct = 0.0; wall_h_corner = 0.0; wall_z_fill = 0.0;
            wall_corner.clear();
            return;
        }

        // THE GAS STATE. side_T_eff is the temperature the wall faces actually
        // see (jet_T_exhaust under side_face_T = live), so the law is derived at
        // the state it is applied at, not at nozzle conditions -- the gas
        // reaching the wall has already impinged.
        //
        // ⚠⚠ P2-3 site 5 of 5, AND THE ONE THAT IS NOT A PER-FACE READ. This
        // derives the law's gas PROPERTIES (rho, nu, cp) and its momentum
        // product -- one law per step, with no column and no height -- so
        // `SideFaceTgas(col, k)` has nothing to be called with here. The
        // meaningful substitution is the state at which the wall section
        // BEGINS, i.e. the gas leaving the floor, which is exactly what
        // `side_T_eff` already is: `jet_T_exhaust` is kept at the FLOOR-EXIT
        // value by the wall march (the mouth goes to the new `jet_T_mouth`
        // column instead). So this site needs no change to stay correct, and
        // the lock asserts that reasoning rather than a call it cannot make.
        //
        // ⚠ Had `jet_T_exhaust` been redefined as the MOUTH temperature, this
        // line would silently derive the wall law at mouth conditions, and -- a
        // worse side effect -- every BELOW-PLANE floor-cone-flank face would
        // have seen the cold mouth gas through the resolver's own fallback,
        // changing the floor side heating that D2w measured as a rate doubler.
        const Set::Scalar T_gas = (side_T_eff > 0.0)
            ? side_T_eff : surface_patch.jet_T_nozzle;
        if (!(T_gas > 0.0)) { wall_corner.clear(); return; }
        wall_cp  = surface_patch.jet_cp;
        wall_rho = WJ_P_ATM / (WJ_R_GAS * T_gas);
        wall_nu  = WallMu(T_gas) / wall_rho;

        // THE MOMENTUM PRODUCT. ⚠⚠ u and D_e are never split: writing u alone
        // would mean quoting a nozzle velocity of ~1500 m/s, which is an
        // artefact of reading a momentum-equivalent diameter against a
        // geometric exit area. The momentum-equivalent pair satisfies
        // rho_inf u^2 D_e^2 = rho_j U_j^2 D^2, so
        //     u D_e = U_j D sqrt(rho_j/rho_inf)
        //           = 4 mdot_full / (pi D sqrt(rho_j rho_inf)),
        // which is convention-free and needs no velocity at all.
        // ⚠ jet_mdot is the SECTOR's share (the parse says so of jet_De_ref);
        // the law needs the FULL burner, hence the azimuth fraction.
        const Set::Scalar f_sec = (wall_sector > 0.0) ? wall_sector : Set::Scalar(1.0);
        const Set::Scalar mdot_full = surface_patch.jet_mdot / f_sec;
        const Set::Scalar rho_j = WJ_P_ATM / (WJ_R_GAS * surface_patch.jet_T_nozzle);
        const Set::Scalar D = surface_patch.jet_D;
        wall_uDe = (D > 0.0 && rho_j > 0.0)
            ? 4.0 * mdot_full / (3.14159265358979323846 * D
                                 * std::sqrt(rho_j * wall_rho))
            : Set::Scalar(0.0);

        // THE DUCT BRANCH. Gnielinski on the concentric annulus between the
        // bore and the shield, with the ENTRAINED mass flow (the recirculating
        // below-nozzle flow m_jet = mdot/phi), floored at the laminar Nusselt
        // number so a correlation cannot undershoot its own floor.
        const Set::Scalar gap = r_bore - WJ_R_SHIELD;
        const Set::Scalar D_h = 2.0 * gap;
        const Set::Scalar A_an = 3.14159265358979323846
                               * (r_bore * r_bore - WJ_R_SHIELD * WJ_R_SHIELD);
        const Set::Scalar phi_e = (jet_T_exhaust > 0.0
                                   && surface_patch.jet_T_nozzle > jet_T_exhaust)
            ? (jet_T_stag - jet_T_exhaust)
              / (surface_patch.jet_T_nozzle - jet_T_exhaust)
            : Set::Scalar(1.0);
        const Set::Scalar m_an = mdot_full
            / ((phi_e > 1e-3 && phi_e <= 1.0) ? phi_e : Set::Scalar(1.0));
        const Set::Scalar k_g = 1.5207e-11 * T_gas * T_gas * T_gas
                              - 4.8574e-8 * T_gas * T_gas
                              + 1.0184e-4 * T_gas - 3.9333e-4;
        Set::Scalar h_duct = 0.0;
        if (A_an > 0.0 && D_h > 0.0 && k_g > 0.0)
        {
            const Set::Scalar Re = (m_an / A_an) * D_h / WallMu(T_gas);
            Set::Scalar Nu = 0.012 * (std::pow(Re, 0.87) - 280.0)
                           * std::pow(WJ_PR, 0.4);
            if (Nu < 5.74) Nu = 5.74;     // laminar annulus floor
            h_duct = Nu * k_g / D_h;
        }
        wall_h_duct = h_duct;

        wall_law_ready = (wall_uDe > 0.0 && wall_rho > 0.0 && wall_nu > 0.0
                          && wall_h_duct > 0.0 && wall_sector > 0.0);
        wall_law_on = wall_law_ready ? 1.0 : 0.0;
        wall_h_corner = wall_law_ready ? WallJetH(0.0, r_bore) : Set::Scalar(0.0);
        wall_z_fill   = WallZFill(r_bore);
        if (!wall_law_ready) wall_corner.clear();
    }

    // =====================================================================
    // P2-2 change 1: UpdatePropertyFields — the per-step property refresh.
    // =====================================================================
    // ⚠⚠ THE PROPERTIES ARE STATIC PER-CELL FIELDS FILLED AT INIT.
    // `kappa_mf` and `beta_mf` are each written at exactly ONE site
    // (`InitializeMicrostructure`), so temperature dependence is NOT a
    // parameter swap -- without a refresh site the criterion fires on the rock
    // the fields were initialised with, whatever the laws say. (a5) verified
    // that premise by counting the writers by DECLARED TYPE.
    //
    // ⚠ `RefreshKappaEff` writes `kappa_eff_mf`, a DIFFERENT field: it applies
    // the damage factor exp(-a_D D) to `kappa_mf`. It is not a temperature
    // refresh and must not be mistaken for one.
    //
    // Called FIRST in Advance, beside P2-1's UpdateWallLawState, so the
    // criterion (`:1436`) and the kernel/march/ledger all see ONE set by
    // construction -- the ordering P2-1 validated.
    //
    // ⚠ `beta_mf` carries the SECANT expansion coefficient, strain/(T - T_ref),
    // NOT alpha(T). That is what makes the seven stress sites compute
    // E(T) * INTEGRAL(alpha) / (1-nu) while still reading one field: every site
    // already forms E * beta * (T - T_ref), so writing the secant here makes
    // the product the accumulated strain identically. ⚠⚠ Writing alpha(T) here
    // instead would give every site alpha(T)*dT -- the planner's wrong reading
    // -- and R would silently absorb the factor. That is gate (b6).
    //
    // No-op, and no traversal, unless some phase sets `law = grimsel`.
    void UpdatePropertyFields(int lev)
    {
        if (!microstructure_enabled) return;
        bool any = false;
        for (const auto& ph : phases) any = any || ph.grimsel_law;
        if (!any) return;
        if (lev >= static_cast<int>(phase_mf.size()) || !phase_mf[lev]) return;

        const int n_phases = static_cast<int>(phases.size());
        std::vector<char> use(n_phases, 0);
        for (int p = 0; p < n_phases; ++p) use[p] = phases[p].grimsel_law ? 1 : 0;
        Set::Scalar acc_k = 0.0, acc_b = 0.0, acc_E = 0.0, acc_T = 0.0, acc_n = 0.0;

        for (amrex::MFIter mfi(*phase_mf[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& ph = (*phase_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T  = (*temp_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& Trf =
                (*T_ref_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& kapa = (*kappa_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& beta = (*beta_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& EE   = (*E_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& muu  = (*mu_phase_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& rho_a = (*rho_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Cp_a  = (*Cp_mf[lev]).array(mfi);
            // ⚠ the phase objects are not device-copyable, so this loop is on
            // the host. It is one pass over the level per step, the same cost
            // class as RefreshKappaEff, which already runs every step.
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                if (!use[p]) continue;
                const Numeric::Material::Grimsel& g = phases[p].grimsel;
                Set::Scalar Tc = T(i, j, k);
                if (!(Tc > 0.0) || !std::isfinite(Tc)) continue;
                const Set::Scalar T_r = Trf(i, j, k);
                kapa (i, j, k) = g.lambda(Tc);
                EE   (i, j, k) = g.E_of(Tc);
                muu  (i, j, k) = g.mu_of(Tc);
                rho_a(i, j, k) = g.rho(Tc);
                Cp_a (i, j, k) = g.cp(Tc);
                // THE SECANT COEFFICIENT -- see the warning above.
                const Set::Scalar dT = Tc - T_r;
                beta (i, j, k) = (std::abs(dT) > 1.0e-9)
                    ? (g.strain_of(Tc) - g.strain_of(T_r)) / dT
                    : g.alpha_of(Tc);
                acc_k += kapa(i, j, k); acc_b += beta(i, j, k);
                acc_E += EE(i, j, k);   acc_T += Tc; acc_n += 1.0;
            }
        }
        Util::RealFillBoundary(*kappa_mf[lev],       geom[lev]);
        Util::RealFillBoundary(*beta_mf[lev],        geom[lev]);
        Util::RealFillBoundary(*E_phase_mf[lev],     geom[lev]);
        Util::RealFillBoundary(*mu_phase_mf[lev],    geom[lev]);
        // the diagnostics the (b5) two-temperature probe reads. ⚠ per BOX, so
        // they MUST be reduced -- the P1g-2c lesson.
        amrex::ParallelDescriptor::ReduceRealSum(acc_k);
        amrex::ParallelDescriptor::ReduceRealSum(acc_b);
        amrex::ParallelDescriptor::ReduceRealSum(acc_E);
        amrex::ParallelDescriptor::ReduceRealSum(acc_T);
        amrex::ParallelDescriptor::ReduceRealSum(acc_n);
        if (acc_n > 0.0)
        {
            prop_kappa_mean = acc_k / acc_n;
            prop_beta_mean  = acc_b / acc_n;
            prop_E_mean     = acc_E / acc_n;
            prop_T_mean     = acc_T / acc_n;
        }
        prop_refresh_n += 1.0;
    }

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

    // D2q-1: number of SOLID cells strictly above cell k, up to and including
    // the column top, indexed [col*nz + (k - klo)] with col as in
    // ColumnTopSolid. The beam's optical path must be measured through
    // MATERIAL, not through geometric distance: with an interior void the beam
    // crosses it unattenuated, so `k_top - k` would overstate the path.
    //
    // While void is contiguous every cell in (k, k_top] is solid and this
    // returns exactly `k_top - k`, which is what keeps the beam runs
    // byte-identical.
    //
    // Cost: one ncols*nz int reduction. Built ONLY when the beam is on (the
    // path argument is unused otherwise -- BeamSource returns 0 at P = 0 and
    // BeamEnergyBalance skips columns with Pi == 0), so the Meier and arbiter
    // runs, which are the large ones, never pay for it.
    amrex::Vector<int> ColumnSolidAboveCount(int lev, const amrex::Vector<int>& k_top)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int klo = domain.smallEnd(2);
        const int nx = domain.length(0), ny = domain.length(1), nz = domain.length(2);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        // Each cell is owned by exactly one box at a level, so a sum over ranks
        // gathers the column without double counting.
        amrex::Vector<int> solid(ncols * static_cast<std::size_t>(nz), 0);
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
                if (rem(i, j, k) <= 0.5)
                    solid[static_cast<std::size_t>((j - jlo) * nx + (i - ilo)) * nz
                          + (k - klo)] = 1;
        }
        amrex::ParallelDescriptor::ReduceIntSum(solid.data(), static_cast<int>(solid.size()));

        amrex::Vector<int> cnt(ncols * static_cast<std::size_t>(nz), 0);
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = k_top[c];
            if (kt == std::numeric_limits<int>::min()) continue;
            int acc = 0;
            for (int k = kt; k >= klo; --k)
            {
                const std::size_t p = c * nz + (k - klo);
                cnt[p] = acc;                 // solid cells in (k, k_top]
                if (solid[p]) ++acc;
            }
        }
        return cnt;
    }

    // D2q-1: length of the CONTIGUOUS solid run below each column's top cell,
    // in cells, capped at the column height. The depth scan walks downward
    // from k_top in 1-cell strides assuming every cell it passes is material;
    // this is how far that assumption holds. Equal to the full distance to the
    // domain floor while void is contiguous, so the scan is unchanged today.
    // Reduces only ncols ints (one per column, like ColumnTopSolid): each rank
    // reports the HIGHEST void cell it owns strictly below that column's top,
    // and ReduceIntMax picks the global one. A full per-(column, k) gather here
    // would cost an ncols*nz reduction on every removal step, which the Meier
    // runs cannot afford.
    amrex::Vector<int> ColumnSolidRunBelowTop(int lev, const amrex::Vector<int>& k_top)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int klo = domain.smallEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        amrex::Vector<int> kvoid(ncols, std::numeric_limits<int>::min());
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t c = static_cast<std::size_t>((j - jlo) * nx + (i - ilo));
                const int kt = k_top[c];
                if (kt == std::numeric_limits<int>::min()) continue;
                int kmax = std::min(hi.z, kt - 1);
                for (int k = kmax; k >= lo.z; --k)
                    if (rem(i, j, k) > 0.5) { if (k > kvoid[c]) kvoid[c] = k; break; }
            }
        }
        amrex::ParallelDescriptor::ReduceIntMax(kvoid.data(), static_cast<int>(kvoid.size()));

        amrex::Vector<int> run(ncols, 0);
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = k_top[c];
            if (kt == std::numeric_limits<int>::min()) continue;
            // Cells kt down to kvoid+1 are solid; with no void below the top
            // the run reaches the domain floor.
            run[c] = (kvoid[c] == std::numeric_limits<int>::min())
                   ? (kt - klo + 1) : (kt - kvoid[c]);
        }
        return run;
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
    // Invariant (void-aware runs), per solid cell, where n = the count of SOLID
    // cells in (k, k_top] (D2q-1: a solid count, no longer k_top - k):
    //   (a) n == 0 <=> k == k_top.
    //       RESTATED in D2q-2b. It used to read "cell above is void or k is the
    //       domain top <=> n == 0", a CONTIGUITY test, which an interior void
    //       breaks: the cell under the void satisfies the left side with solid
    //       still above it. The current form is exact at any topology because
    //       k_top is itself solid, so the count vanishes only at the top cell.
    //   (b) RETIRED in D2q-1. It asserted phi - n*dz in [0, dz), i.e. that phi
    //       and the `removed` mask agree -- that phi is a valid void oracle,
    //       which is exactly the assumption D2q-1 removed.
    // The P_inc tally below is gated on the SAME k == k_top test, for the same
    // reason: the old contiguity form counted incident power twice per column
    // once an interior void existed (measured 14.29 % closure error).
    // Single level: lev == 0 only (multi-level unsupported; composite sums
    // would double-count covered cells).
    void BeamEnergyBalance(int lev, Set::Scalar time, Set::Scalar dt,
                           const amrex::Vector<int>& k_top_col)
    {
        if (!SurfaceOwner(lev)) return;
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
        // D2q-1: solid-cell count above k, as in AdvanceMicrostructure.
        const amrex::Vector<int> solid_above =
            (void_aware && P_now > 0.0 && !k_top_col.empty())
            ? ColumnSolidAboveCount(lev, k_top_col) : amrex::Vector<int>();
        const int* sab    = solid_above.data();
        const int  nz_dom = domain.length(2);

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
                const int col_b = (j - jlo_dom) * ncx + (i - ilo_dom);
                const int n_above = void_aware
                    ? (sab ? sab[static_cast<std::size_t>(col_b) * nz_dom + (k - klo_dom)]
                           : k_top_col[col_b] - k)
                    : 0;
                const Set::Scalar path = void_aware ? SolidFacePath(n_above, dz)
                                                    : Set::Scalar(-1.0);
                // D2q-2b: the column top, for BOTH invariant (a) and the P_inc
                // tally below. This replaces the old `is_top_solid` test --
                // "the cell directly above is void, or this is the domain top"
                // -- which is a contiguity test and counts MORE THAN ONE cell
                // per column the moment an interior void exists: the true top,
                // and again the cell sitting just under the void. That
                // double-counted the incident beam power (measured: a 14.29 %
                // closure error against a predicted 140/120 over-count on the
                // seeded probe). Identical while void is contiguous.
                const bool at_col_top = void_aware ? (k == k_top_col[col_b])
                                                   : (k == khi_dom);
                if (void_aware)
                {
                    // (a) RESTATED (D2q-2b item 0.3): `n_above == 0` iff
                    // `k == k_top`.
                    //
                    // The old form compared against `is_top_solid` -- "the cell
                    // directly above is void" -- which is a CONTIGUITY test: a
                    // cell sitting just under an interior void satisfies it
                    // while still having solid above, so it would be counted as
                    // a violation the moment per-cell removal creates one.
                    //
                    // The restated form is exact at any topology, because
                    // n_above counts SOLID cells in (k, k_top] and k_top is
                    // itself solid: below the top there is always at least one,
                    // so the count vanishes exactly at the top cell. With void
                    // contiguous the two forms agree, which keeps identity.
                    if (at_col_top != (n_above == 0)) n_viol += 1.0;
                    // D2q-1: invariant (b) -- phi - n*dz in [0, dz) -- is
                    // RETIRED. It asserted that phi and the `removed` mask
                    // agree, i.e. that phi is a valid void oracle, which is
                    // exactly the assumption this step removes. It reads 0
                    // violations under both forms today, and n_viol feeds the
                    // `beam_invariant_viol` thermo.dat column, so dropping it
                    // is byte-identical.
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
                if (at_col_top) P_inc += (1.0 - beam.R) * Pi * A;
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
        if (P == 0.0) return 0.0; // perf: beam off (Pi would be 0 anyway)
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
        // AMR surface band: every column's top solid cell (on this level's
        // view of the mask) plus refine_height above and refine_depth below.
        // Tagging on every level below max_level puts the whole surface, and
        // the depth-scan / beam window under it, on the finest level.
        if (!(removal_enabled() && surface_refine_depth > 0.0)) return;
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int nx  = domain.length(0);
        const Set::Scalar dz = geom[lev].CellSize()[2];
        const int n_below = static_cast<int>(std::ceil(surface_refine_depth  / dz));
        const int n_above = static_cast<int>(std::ceil(surface_refine_height / dz));
        const amrex::Vector<int> k_top = ColumnTopSolid(lev);
        const int* ktp = k_top.data();
        for (amrex::MFIter mfi(*removed_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const int kt = ktp[(j - jlo) * nx + (i - ilo)];
                if (kt == std::numeric_limits<int>::min()) return;
                if (k >= kt - n_below && k <= kt + n_above) tags(i, j, k) = amrex::TagBox::SET;
            });
        }
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

        // Step 20: surface-normal (slope-aware) spall CRITERION.
        //
        // ⚠ D2z-1: the RECESSION half of Step 20 no longer applies. Whole-cell
        // removal discards the scan's result, and the h/cos(theta) recession
        // described below lived inside the scan, so a column recedes exactly
        // one cell per firing regardless of slope. The paragraph below is kept
        // because the CRITERION half is unchanged and still governs when a
        // sloped column fires; read "and a flake of normal thickness h recedes
        // the column vertically by h/cos(theta)" as historical.
        //
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
        // D2z-1: whole-cell removal discards the scan's result (h_c = dz), and
        // Step 20's slope scaling lives INSIDE the scan
        // (MMWSpalling/Removal.H:691 and :992), while recession is a plain phi
        // shift with no cos(theta) at its own site (:1444). So the key keeps
        // its slope-aware FIRING DECISION -- K_I is still sampled at
        // depth/cos(theta) in the criterion update (:4668), which never goes
        // through the scan -- but loses its slope-aware RECESSION: a wall
        // column drops one cell per firing whatever its slope.
        //
        // This is a MESSAGE, not an abort. Four inputs set the key and two of
        // them (unit/beam_void_closure and the meier input_2d_dev "dev2d" leg)
        // are inside the 2639-file identity sweep, so aborting would stop the
        // sweep from running them at all -- a capability notice must not break
        // the gate that protects everything else. The on-disk witness is the
        // thermo column sn_recession_superseded, registered under exactly this
        // key: it reports the mean 1/cos(theta) over each step's firing
        // columns, i.e. the factor the discarded scaling would have applied.
        if (value.spall_surface_normal)
            Util::Message(INFO, "spall.surface_normal = 1: the slope-aware "
                                "CRITERION is retained, but the slope-aware "
                                "RECESSION scaling is superseded by whole-cell "
                                "removal -- a column now drops exactly one cell "
                                "per firing regardless of slope. Watch the "
                                "thermo column sn_recession_superseded (1.0 = "
                                "nothing lost on that step).");

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

        // ------------------------------------------------------------------
        // P0. BOTH keys are parsed BEFORE the spall.enabled early-return, for
        // the same reason spall.stress_source is: the surface temperature also
        // governs the Sp_field ONSET diagnostic (UpdateSpAfterMechanics),
        // which runs under spallation.model = sp_weibull with material removal
        // OFF -- that is exactly the configuration of the Kant onset anchor.
        // Parsing after the return would leave the key silently unread and the
        // anchor leg would look like a null (D2z-0 lesson 1).
        std::string surfT;
        pp.query_default("spall.surface_T", surfT, std::string("cell_average"));
        if (surfT == "cell_average")
        {
            value.sp_surface_T_mode  = SurfaceTMode::CellAverage;
            value.sp_surface_T_recon = false;
            value.sp_surface_T_apply = false;
        }
        else if (surfT == "diagnostic")
        {
            value.sp_surface_T_mode  = SurfaceTMode::Diagnostic;
            value.sp_surface_T_recon = true;
            value.sp_surface_T_apply = false;
        }
        else if (surfT == "reconstructed")
        {
            value.sp_surface_T_mode  = SurfaceTMode::Reconstructed;
            value.sp_surface_T_recon = true;
            value.sp_surface_T_apply = true;
        }
        else
            Util::Abort(INFO, "spall.surface_T must be 'cell_average', "
                              "'diagnostic' or 'reconstructed', got: ", surfT);

        // P0b deliverable 3: the flake thickness is c_f * ell_eff, with NO dz
        // clamp. ⚠ The min(., dz) clamp P0 used is the measured culprit for the
        // 20.6 % residual in its mesh ladder: at 0.5 mm, ell_eff (0.992 mm) is
        // already larger than dz, so the "flake" was the whole cell again and
        // T_bar collapsed back onto the cell average.
        // ⚠ c_f IS THE ONE EMPIRICAL CONSTANT THIS PACKET INTRODUCES. Under
        // `flake` the thickness is IMPOSED (Removal.H bypasses the depth scan's
        // h_col), so it is not derived from fracture here. 0.3-0.5 is what
        // Wilkinson's 100-1000 um spall thicknesses and Rauenzahn's size
        // distributions imply at ell_eff = 0.8-1.6 mm. Its lever on the level is
        // weak (about +-5 % across 0.3-0.5), which is the real argument that it
        // is not a tuning knob; the successor is to run the K_I(z) scan on the
        // reconstructed exponential so c_f becomes an OUTPUT (P1).
        // ⚠⚠ P1 -- THE GRADIENT-AWARE DRIVING STRESS.
        //
        // `absolute` (default, byte-identical) drives the crack on the stress
        // measured from T_ref, i.e. on the ABSOLUTE superheat. That criterion
        // is GRADIENT-BLIND: stress is local in temperature (Boley & Weiner,
        // half-space exact), so a uniformly preheated slab at the firing
        // temperature fires exactly as a steep thermal skin does. Measured on
        // the Meier rev-7 run, the superheat at firing reads 469 / 468 / 469 /
        // 465 / 466 K across the five radial bands -- identical from axis to
        // wall -- while T_surf - T_top reads 250 / 202 / 135 / 61 / 31 K. The
        // information is there; the criterion discards it.
        pp.query_default("spall.flake_c_f", value.sp_flake_c_f, 0.4);
        if (!(std::isfinite(value.sp_flake_c_f) && value.sp_flake_c_f > 0.0))
            Util::Abort(INFO, "spall.flake_c_f must be finite and > 0, got ",
                        value.sp_flake_c_f);
        // Revision 7 deliverable 1. Default 0 = INERT, so every frozen leg and
        // every control is byte-identical without it. See ReconFromX for what
        // this bound is and, more importantly, what it is not.
        pp.query_default("spall.flake_h_max", value.sp_flake_h_max, 0.0);
        if (!(std::isfinite(value.sp_flake_h_max) && value.sp_flake_h_max >= 0.0))
            Util::Abort(INFO, "spall.flake_h_max must be finite and >= 0 "
                              "(0 = inert), got ", value.sp_flake_h_max);
        if (value.sp_flake_h_max > 0.0 &&
            value.sp_flake_h_max >= value.spall_h_spall_max)
            Util::Abort(INFO, "spall.flake_h_max (", value.sp_flake_h_max,
                        ") must be SMALLER than spall.h_spall_max (",
                        value.spall_h_spall_max, "): the first is a physical "
                        "flake thickness, the second is the depth the scan "
                        "searches. A bound at or above the search depth is "
                        "inert and hides that the search depth is doing the "
                        "bounding -- which is the fault this key exists to fix.");

        // P1f: the lateral mode. Default 0 = INERT.
        // ⚠ These are `spall.*`, so the `surface_patch.side_face_*` else-chain that
        // refuses that family when side_face_flux is off does NOT apply here; the
        // guards below are this pair's only protection. AMReX merely WARNS about an
        // unused ParmParse variable, so a silent no-op is the failure mode.
        pp.query_default("spall.lateral_drive",   value.sp_lateral_drive,   false);
        pp.query_default("spall.lateral_removal", value.sp_lateral_removal, false);
        if (value.sp_lateral_removal && !value.sp_lateral_drive)
            Util::Abort(INFO, "spall.lateral_removal = 1 requires "
                              "spall.lateral_drive = 1: removal voids the cell the "
                              "lateral criterion fired, so without the drive there "
                              "is nothing to remove.");
        // P1g-1 deliverable 1: SUB-CELL LATERAL RECESSION, default OFF.
        // ⚠⚠ Why a third key rather than changing regime 6 in place. Regime 6
        // today voids a whole 2 mm cell for a flake the criterion earned at
        // 1.5 mm -- a 1.33x over-removal on 99.9 % of firings, because
        // h_f_lat sits at the spall.flake_h_max cap. Banking the thickness
        // instead is a CHANGE OF REMOVAL LAW, not a bug fix, so it is keyed
        // and the old behaviour remains reachable for comparison.
        // ⚠ The conjunction is deliberately three-deep -- subcell requires
        // removal, which already requires drive -- so the banking block below
        // is unreachable unless all three are set, and (c3) scores exactly that.
        pp.query_default("spall.lateral_subcell", value.sp_lateral_subcell, false);
        // P1g-2 deliverable 2: THE THERMAL HALF of lateral sub-cell recession.
        // ⚠⚠ NO `requires` ABORT, and that is a lesson paid for: P1g-1 shipped
        // `lateral_subcell requires lateral_removal` and its own negative-control
        // leg (`drive = 0` + `subcell = 1`) then returned rc 6. A control leg that
        // must RUN AND BE INERT cannot be guarded in Parse. Reachability lives at
        // the use site, as `bank_on && sp_lateral_thermal`, and (c5) scores it.
        pp.query_default("spall.lateral_thermal", value.sp_lateral_thermal, false);
        // ⚠⚠ P1g-3 deliverables 5 + 6: TWO KEYS, NOT ONE, AND THAT IS THE POINT.
        // The step changes the flake LENGTH estimate and the FRACTION of it that
        // leaves. One key would confound them, and (b4)'s 2x2 is what separates
        // their levers (FRAC is the larger: x0.29 against LEN's x0.80).
        {
            std::string fl;
            pp.query_default("spall.flake_length", fl, std::string("twocell"));
            if      (fl == "twocell") value.sp_flake_len_flux = false;
            else if (fl == "flux")    value.sp_flake_len_flux = true;
            else Util::Abort(INFO, "spall.flake_length must be 'twocell' or "
                                   "'flux', got: ", fl);
            std::string ff;
            pp.query_default("spall.flake_fraction", ff, std::string("cf"));
            if      (ff == "cf")      value.sp_flake_frac_weibull = false;
            else if (ff == "weibull") value.sp_flake_frac_weibull = true;
            else Util::Abort(INFO, "spall.flake_fraction must be 'cf' or "
                                   "'weibull', got: ", ff);
        }
        // ⚠⚠ DELIBERATELY NOT an abort-if-unaccompanied, unlike lateral_removal
        // above. The packet's NEG leg is exactly `lateral_drive = 0` with
        // `lateral_subcell = 1`, and its bar is that banking events are ZERO --
        // so this combination must RUN and be inert, not abort. An abort here
        // makes (c3) unrunnable, which is how it was caught. Reachability is
        // enforced at the use site by the three-deep conjunction
        //   lat_on && sp_lateral_removal && sp_lateral_subcell
        // which is what the reachability ledger names as the guard.
        // ⚠⚠ P1e IS A PREREQUISITE, NOT A CARRIED-OVER SETTING. The lateral mode is
        // grid-isotropic only because the staircase correction makes side heat
        // isotropic per unit of wall AREA; with it off the lateral drive inherits
        // the same sqrt(2) over-count the vertical mode has, and criterion (c)
        // fails by construction rather than by physics.
        // ⚠ surface_patch is parsed at :139 and ParseSpallSettings at :158, so this
        // flag IS set by the time this runs -- the conjunction was checked, not
        // assumed (the P1f revision-3 blocker was exactly an unchecked conjunction).
        // ⚠⚠ P1f deliverable 7 -- (a1)'s CELL-FORM HALF, AS A DIRECT CALL.
        // The lateral pass only ever calls SurfaceCellFlux with face_form = true,
        // so the cell-form branch -- the one that must write NOTHING to
        // T_face_out -- never gets exercised in anger, and a wrong write there
        // would be SILENT: the criterion would receive a cell average dressed as
        // a face temperature. A source-level assertion (three write sites, none
        // after the cell-form marker) is not runtime cover, so this is the real
        // thing: one direct call with face_form = false and the out-parameter
        // pre-set to a sentinel, checked on every run that turns the mode on.
        // ⚠ Pure arithmetic, so it is safe at Parse time -- no fields needed.
        if (value.sp_lateral_drive || value.sp_lateral_removal)
        {
            Set::Scalar q_i = 0, q_o = 0, q_r = 0, q_l = 0;
            const Set::Scalar sentinel = -12345.0;
            Set::Scalar probe = sentinel;
            // face_form = false => the CELL form. Representative arguments; the
            // branch is selected by the flags, not by the values.
            SurfaceCellFlux(/*in_patch=*/true, /*face_form=*/false, 1.0,
                            500.0, 1200.0, 1.5, 2.0e-3, 700.0,
                            0.8, 0.9, 1000.0, 5.67e-8, 10.0, 293.15,
                            q_i, q_o, q_r, q_l,
                            /*T_pin=*/0.0, /*branch=*/nullptr,
                            /*T_surf_in=*/0.0, &probe);
            if (probe != sentinel)
                Util::Abort(INFO, "P1f self-test: SurfaceCellFlux wrote ", probe,
                            " to T_face_out on its CELL-FORM branch, which must "
                            "write nothing -- q_robin there is h*(T_gas - Tc) on "
                            "the CELL AVERAGE, and Tc is not a face temperature. "
                            "The lateral criterion would be handed the wrong "
                            "quantity with no symptom. Sentinel was ", sentinel);
            // ⚠ And the face form MUST write, or the lateral mode is silently
            // starved: the companion half of the same contract.
            probe = sentinel;
            SurfaceCellFlux(true, /*face_form=*/true, 1.0,
                            500.0, 1200.0, 1.5, 2.0e-3, 700.0,
                            0.8, 0.9, 1000.0, 5.67e-8, 10.0, 293.15,
                            q_i, q_o, q_r, q_l, 0.0, nullptr, 0.0, &probe);
            if (probe == sentinel || !(probe > 0.0))
                Util::Abort(INFO, "P1f self-test: SurfaceCellFlux did NOT write a "
                            "positive face temperature to T_face_out on its "
                            "face-Newton branch; the lateral criterion has no "
                            "temperature to read and would be inert.");
        }
        if ((value.sp_lateral_drive || value.sp_lateral_removal)
            && !value.surface_patch.side_face_staircase)
            Util::Abort(INFO, "spall.lateral_drive / spall.lateral_removal require "
                              "surface_patch.side_face_staircase = 1: the lateral "
                              "mode is grid-isotropic only because the staircase "
                              "area correction is applied. Without it the lateral "
                              "drive carries the same sqrt(2) azimuthal over-count "
                              "the per-column criterion does, and the mode would be "
                              "measuring the grid.");

        std::string remE;
        pp.query_default("spall.removal_energy", remE, std::string("cell"));
        if      (remE == "cell")  value.sp_removal_energy_flake = false;
        else if (remE == "flake") value.sp_removal_energy_flake = true;
        else
            Util::Abort(INFO, "spall.removal_energy must be 'cell' or 'flake', "
                              "got: ", remE);

        // flake recession recedes phi by h_f = min(ell_eff, dz). ell_eff has
        // exactly one source -- the surface reconstruction -- so the pair is
        // not independent and the combination is rejected rather than silently
        // falling back to whole-cell.
        if (value.sp_removal_energy_flake && !value.sp_surface_T_apply)
            Util::Abort(INFO,
                "spall.removal_energy = flake requires spall.surface_T = "
                "reconstructed: the flake thickness is h_f = min(ell_eff, dz) "
                "and ell_eff is produced only by the surface reconstruction. "
                "Set spall.surface_T = reconstructed, or leave "
                "spall.removal_energy = cell.");

        // ⚠ SCOPE, part 1 (the D2z-1 lesson, restated). D2z-1 found that a
        // default-off key HIDES THE BLAST RADIUS OF ITS OWN CONDITION: its
        // override also caught the damage-law path, where a multi-cell flake is
        // the intended answer, and validation/hu/hu_end_to_end (real removal,
        // never run) would have moved silently. So both keys are admitted ONLY
        // on the sp_weibull computed-T path, and that is enforced at PARSE
        // time, where it aborts, rather than by a runtime condition that would
        // leave a mis-scoped input looking like a null leg.
        // `spall.stress_source` is parsed a few lines above, so it can be
        // checked here; `spallation.model` is NOT -- ParseSpallSettings runs at
        // Parse():158 and ParseSpallationModel at :174 -- so that half of the
        // scope check lives in the late validation block instead.
        if ((value.sp_surface_T_recon || value.sp_removal_energy_flake) &&
            !value.sp_local_thermoelastic)
            Util::Abort(INFO,
                "spall.surface_T = reconstructed requires "
                "spall.stress_source = local_thermoelastic: the "
                "reconstruction rebuilds sigma_xx and K_Ic from the "
                "temperature profile, and the mechanics (stress_mf) path "
                "takes sigma_xx from the elastic solve instead.");


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
        // D2q-2b item 5: per-cell removal. Default 0 = off and byte-identical
        // to every pre-D2q-2b result.
        //
        // Preconditions are enforced here rather than trusted to the input:
        //  * `surface_patch.side_face_flux = 1` -- without it an exposed
        //    lateral face is never heated, so the criterion evaluated there
        //    reads bulk rock and the mechanism is silently inert.
        //  * `surface.follow_mask = 1` -- this one is NOT cosmetic. `phi` is
        //    still the sole input to UpdateSurfaceMaskFromPhi (D2q-1 site 7,
        //    unconverted). Per-cell removal can void a cell without moving the
        //    column's phi, so the window `0 <= phi < dz` would keep pointing at
        //    a now-void cell and the surface flux would be applied to a void.
        //    follow_mask takes the surface from the rem-based k_top instead.
        //    See RESULTS.md Q1.
        int pcr = 0;
        pp.query_default("spall.per_cell_removal", pcr, 0);
        if (pcr != 0 && pcr != 1)
            Util::Abort(INFO, "spall.per_cell_removal must be 0 or 1, got: ", pcr);
        value.spall_per_cell_removal = (pcr == 1);


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
                    // D2q-1: removed_mf is authoritative (see ColumnTopSolid).
                    if (rem(i, j, k) <= 0.5) { ktl = k; break; }
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
        if (!SurfaceOwner(lev)) return;   // AMR: Sp lives on the surface-owning level
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
        // P2-2 change 4, sites 1-3 of 7 (the integrator's three branches).
        const Set::Scalar te_relief = sp_relief;
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
        // P2-2b deliverable 4: the rock law for the Sp_field DIAGNOSTIC site.
        // ⚠⚠ THE TWO-FIRING-SITE INVARIANT. This site and `Removal.H`'s depth
        // scan are the two places the criterion is formed, and P0 deliberately
        // closed the gap between them; both are switched by the SAME
        // `RockLaw()`, hoisted here the way `te_E`/`te_beta` already are so a
        // `[=, this]` closure never reads `phases` itself.
        const Numeric::Material::Grimsel *te_g = use_local_TE ? RockLaw()
                                                              : nullptr;
        const Set::Scalar a0_loc      = sp_a0;
        const int         n_seg       = sp_n_segments;
        const KIcLaw      law         = kic_law;
        const Set::Scalar kic_scale_  = kic_scale;
        // ⚠ hoisted for the same reason te_g is: this site sits in a
        // `[=, this]` ParallelFor and must not read `phases` itself.
        const Numeric::Material::Grimsel *kic_g =
            phases.empty() ? nullptr : &phases[0].grimsel;
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

        // ------------------------------------------------------------------
        // P0b: the reconstruction is now built ONCE per step for the whole
        // level (BuildSurfaceTColumns, first thing in Advance) and read here.
        // P0's own gather at this site is deleted -- see SharedSurfaceRecon.
        const bool recon_on = sp_surface_T_recon && use_local_TE &&
                              !use_presc_T && !use_presc_s;
        if (recon_on)
        {
            T_surf_field_mf[lev]->setVal(0.0);
            ell_eff_field_mf[lev]->setVal(0.0);
        }
        // Witness counters for "the reconstruction is inert" (see the packet's
        // witness table). Reduced and reported on thermo.dat only under the key.
        long rec_n = 0, rec_n_flat = 0, rec_n_ellclamp = 0, rec_n_Tclamp = 0,
             rec_n_fallback = 0, rec_n_hfcap = 0, rec_n_hfbnd = 0;
        Set::Scalar rec_sum_dT = 0.0;

        // P2-2e: the constraint reference at THIS site. Per-box accumulation,
        // hand-reduced once after the loop (the rule P1g-1/P1g-2 paid for).
        // ⚠ Only the local_thermoelastic branch carries a stress reference; the
        // prescribed-sigma and stress_mf branches have none to move.
        const bool        cref_on_sp = (sp_cref_d > 0.0) && use_local_TE;
        const Set::Scalar cref_d_sp  = sp_cref_d;
        long sp_cref_vert = 0, sp_cref_lat = 0, sp_cref_clamp = 0, sp_cref_neg = 0;

        // Sp is a per-column diagnostic; clear everywhere first so cells
        // below the top z-slab stay zero on plotfiles.
        Sp_field_mf[lev]->setVal(0.0);

        // D2q-2b item 6: per-cell removal scores every exposed face, so the
        // six-neighbour void test below needs a current ghost layer on
        // removed_mf. Parse() guarantees removal is enabled whenever the key
        // is on.
        const bool pcr_on = spall_per_cell_removal;
        const int dlo_pc[3] = {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
        const int dhi_pc[3] = {domain.bigEnd(0),   domain.bigEnd(1),   domain.bigEnd(2)};
        if (pcr_on) Util::RealFillBoundary(*removed_mf[lev], geom[lev]);

        for (amrex::MFIter mfi(*Sp_field_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const int vbx_lo_k = mfi.validbox().smallEnd(2);
            amrex::Array4<Set::Scalar>       const& Sp_a  = (*Sp_field_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T_a   = (*temp_mf[lev]).array(mfi);
            // P2-2a change A: the stress-free reference, READ at the sampled
            // cell and never written here. Left unbound (and the lambda falls
            // back to the declared constant) if the field does not exist.
            amrex::Array4<const Set::Scalar> T0_a;
            if (T_init_mf.size() > static_cast<std::size_t>(lev) && T_init_mf[lev])
                T0_a = (*T_init_mf[lev]).const_array(mfi);
            // Bound to Sp_field_mf as a harmless placeholder when the key is
            // off; the lambda never reads it then (mirrors the AdvanceMicro-
            // structure convention for rem/phi).
            amrex::Array4<const Set::Scalar> const& rem_pc = pcr_on
                ? (*removed_mf[lev]).const_array(mfi)
                : (*Sp_field_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Matrix> sig_node;
            if (!use_presc_s && !use_local_TE) sig_node = (*stress_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> flaw_a;
            if (use_weibull) flaw_a = (*flaw_a_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> Ts_a, le_a;
            if (recon_on)
            {
                Ts_a = (*T_surf_field_mf[lev]).array(mfi);
                le_a = (*ell_eff_field_mf[lev]).array(mfi);
            }

            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const int ksurf = k_top_col_p[col_index_sp(i, j)];
                // D2q-2b item 6: with per-cell removal the criterion is
                // evaluated on EVERY exposed face, not only the column top.
                // `nx_in, ny_in, nz_in` is the inward normal of the face being
                // scored -- (0,0,-1) is the top face, i.e. exactly the
                // pre-D2q-2b behaviour, which is why the key-off path below is
                // the original code unchanged.
                int nx_in = 0, ny_in = 0, nz_in = -1;
                if (!pcr_on)
                {
                    if (k != ksurf) return;
                }
                else
                {
                    if (rem_pc(i, j, k) > 0.5) return;       // void: nothing to score
                    // An exposed face has a void (or out-of-domain) neighbour.
                    // removed_mf carries one ghost, refilled before this call,
                    // so the six-neighbour read is valid at box seams; the
                    // domain top counts as exposed, matching `k == ksurf`.
                    const bool ez = (k == top_k) || rem_pc(i, j, k + 1) > 0.5;
                    const bool exm = (i > dlo_pc[0]) && rem_pc(i - 1, j, k) > 0.5;
                    const bool exp_ = (i < dhi_pc[0]) && rem_pc(i + 1, j, k) > 0.5;
                    const bool eym = (j > dlo_pc[1]) && rem_pc(i, j - 1, k) > 0.5;
                    const bool eyp = (j < dhi_pc[1]) && rem_pc(i, j + 1, k) > 0.5;
                    if (!(ez || exm || exp_ || eym || eyp)) return;
                    // Score the face whose inward normal is most favourable.
                    // At a0 << dx every sample lands in this cell, so all
                    // normals give the same Sp and the choice is immaterial;
                    // the loop is written generally so it stays correct if the
                    // flaw ever approaches a cell width.
                    if (!ez)
                    {
                        nz_in = 0;
                        if      (exm) { nx_in =  1; }
                        else if (exp_) { nx_in = -1; }
                        else if (eym) { ny_in =  1; }
                        else           { ny_in = -1; }
                    }
                }
                // Receded surface top face (== domain top z when nothing above
                // this column has been removed, since then ksurf == top_k).
                const Set::Scalar z_surf_face =
                    z_top_face - Set::Scalar(top_k - ksurf) * dz;

                Set::Scalar xc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                Set::Scalar yc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];

                Set::Scalar Tsurf = use_presc_T ? presc_T_f(time) : T_a(i, j, k);

                // P0: reconstruct the SURFACE temperature from the top two
                // cell averages. ⚠ Applied to the TOP FACE only: under
                // spall.per_cell_removal a lateral face has no column profile,
                // and per-cell removal stays off in this programme -- so a
                // lateral face keeps the cell average and that is stated, not
                // silently assumed. Needs a solid second cell; when the
                // contiguous solid run below the top is a single cell there is
                // no profile to read and the cell average stands (counted).
                SurfaceRecon rec;
                bool rec_ok = false;
                if (recon_on && nz_in == -1 && nx_in == 0 && ny_in == 0)
                {
                    const std::size_t c_rec =
                        static_cast<std::size_t>(col_index_sp(i, j));
                    if (SharedSurfaceRecon(c_rec, rec))
                    {
                        rec_ok = sp_surface_T_apply;
                        ++rec_n;
                        if (rec.flat)        ++rec_n_flat;
                        if (rec.ell_clamped) ++rec_n_ellclamp;
                        if (rec.T_clamped)   ++rec_n_Tclamp;
                        if (rec.h_f_capped)  ++rec_n_hfcap;
                        if (rec.h_f_bounded) ++rec_n_hfbnd;
                        rec_sum_dT += rec.T_surf - T_a(i, j, k);
                        Ts_a(i, j, k) = rec.T_surf;
                        le_a(i, j, k) = std::isfinite(rec.ell_eff)
                                        ? rec.ell_eff : 0.0;
                    }
                    else ++rec_n_fallback;
                }
                const Set::Scalar T_top_cell = T_a(i, j, k);
                // T(z) inside the top half-cell: the surface endpoint is the
                // reconstruction, the dz/2 endpoint is the cell average. This
                // is the SAME lerp Removal.H's depth scan uses, with the same
                // endpoints, which is what makes the two criterion sites one
                // formula rather than two.
                auto T_recon_at = [&] (Set::Scalar dv) -> Set::Scalar {
                    const Set::Scalar half_dz = 0.5 * dz;
                    Set::Scalar al = dv / half_dz;
                    if (al < 0.0) al = 0.0;
                    if (al > 1.0) al = 1.0;
                    return (1.0 - al) * rec.T_surf + al * T_top_cell;
                };

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
                // ⚠⚠ P1: BOTH firing sites must change together. P0 closed the
                // gap between this Sp_field site and the depth scan
                // deliberately; changing one would re-open it, and the Meier
                // run fires on the scan while the ONSET diagnostic fires here.
                // ⚠ P2-2a: that invariant is why change A is applied at BOTH,
                // and why `differential` was deleted at both.
                //
                // THE STRESS-FREE REFERENCE at a sampled cell. ⚠ Falls back to
                // the DECLARED constant where the field is absent or
                // non-finite, which is the pre-P2-2a behaviour exactly -- a
                // cell with no stored value must not silently acquire a zero
                // bracket and stop firing.
                auto T0_at = [&] (int ii, int jj, int kk) -> Set::Scalar {
                    if (!T0_a) return te_Tref;
                    const Set::Scalar v = T0_a(ii, jj, kk);
                    return std::isfinite(v) ? v : te_Tref;
                };
                // ⚠⚠ P2-2e: THE REFERENCE OF THE DRIVE at normal depth `depth`.
                // d = 0: the stored initial temperature at the sampled cell --
                // exactly the pre-P2-2e argument, so the expressions below are
                // bit-identical. d > 0:
                //  * TOP face: the LIVE temperature d further behind the face,
                //    read from the gathered column (CrefTAtDepth) -- the SAME
                //    values the depth scan reads. NOT T_a(i_d, j_d, k_d): that is
                //    clamped at this box's bottom and would read the wrong cell at
                //    a seam once the read is d/dz cells deep.
                //  * LATERAL face (reached only under spall.per_cell_removal,
                //    OFF on every production leg): the packet's lateral form
                //    T0 + dT_s exp(-(xi + d)/ell) on THIS branch's own profile.
                //    This branch has no flux length -- its sample is the cell
                //    average at every depth -- so `ell` is its two-cell SLOPE
                //    along the inward normal, the same construction as the
                //    lateral criterion's own `ell` (Removal.H), one cell inward
                //    and therefore inside temp_mf's ghost layer (the read that
                //    criterion already makes). No inward decay (soaked or
                //    inverted) gives a zero bracket: the uniform-slab limit.
                //    ⚠ Declared: Removal.H's lateral criterion references with
                //    its flux length ell_use; the two lateral sites were never
                //    one model (this one is isothermal across the flake).
                auto Tref_at = [&] (Set::Scalar depth_n,
                                    int ii, int jj, int kk) -> Set::Scalar {
                    if (!cref_on_sp) return T0_at(ii, jj, kk);
                    if (nx_in == 0 && ny_in == 0)
                    {
                        const std::size_t c =
                            static_cast<std::size_t>(col_index_sp(i, j));
                        return CrefTAtDepth(c, (depth_n + cref_d_sp) * icth
                                               + Set::Scalar(cref_ktop[c] - ksurf) * dz,
                                            dz, sp_cref_vert, sp_cref_clamp);
                    }
                    ++sp_cref_lat;
                    const Set::Scalar T0c = T0_at(i, j, k);
                    const Set::Scalar dTc = T_a(i, j, k) - T0c;
                    const Set::Scalar dTb = T_a(i + nx_in, j + ny_in, k) - T0c;
                    if (!(dTc > 0.0) || !(dTb > 0.0) || dTc <= dTb)
                        return T_a(i, j, k);
                    const Set::Scalar ell_s = DX[0] / std::log(dTc / dTb);
                    if (!(ell_s > 0.0) || !std::isfinite(ell_s))
                        return T_a(i, j, k);
                    return T0c + dTc * std::exp(-(depth_n * icth + cref_d_sp) / ell_s);
                };
                auto sigma_at_depth = [&] (Set::Scalar depth) -> Set::Scalar {
                    const Set::Scalar dv = depth * icth;  // vertical sample depth
                    if (use_presc_s)
                    {
                        Set::Scalar z_phys = z_surf_face - dv;
                        return -presc_s_f(xc, yc, z_phys, time) + sp_conf_offset;
                    }
                    if (use_local_TE)
                    {
                        // P0: inside the top half-cell the profile is the
                        // surface lerp, not the flat cell average. rec_ok is
                        // false with the key off, so the pre-P0 expression
                        // below is reached unchanged and bit-for-bit.
                        if (rec_ok && dv <= 0.5 * dz)
                        {
                            // P2-2e: the drive's reference -- T0 at d = 0
                            // (bit-identical), T(depth + d) under d > 0.
                            const Set::Scalar Tref_h = Tref_at(depth, i, j, k);
                            if (cref_on_sp && T_recon_at(dv) < Tref_h) ++sp_cref_neg;
                            // P2-2b deliverable 4, law path.
                            if (te_g)
                                return te_relief
                                     * LawSigma(te_g, T_recon_at(dv),
                                                Tref_h)
                                     + sp_conf_offset;
                            // ⚠ pre-P2-2b expression, VERBATIM.
                            return te_relief * te_E * te_beta
                                   * (T_recon_at(dv) - Tref_h)
                                   / (1.0 - te_nu) + sp_conf_offset;
                        }
                        // 1-D-confined thermoelastic sigma_xx from T, same form
                        // and compression-positive convention as the h_spall
                        // scan (Removal.H). Sample temp_mf at the cell
                        // containing the requested depth along the INWARD FACE
                        // NORMAL, clamped into the box.
                        //
                        // D2q-2b: for the top face (nz_in = -1, nx_in = ny_in
                        // = 0) this is bit-for-bit the pre-D2q-2b expression.
                        // For a lateral face it steps sideways into the rock
                        // instead of downwards -- the same rule, evaluated
                        // where it currently is not. At a0 = 20 um against
                        // dx = 2 mm the step count is 0 for every sample, so
                        // the criterion is a function of THIS CELL alone.
                        const int nstep = static_cast<int>(std::floor(dv / dz));
                        int i_d = i + nx_in * nstep;
                        int j_d = j + ny_in * nstep;
                        int k_d = k + nz_in * nstep;
                        if (k_d > ksurf)    k_d = ksurf;
                        if (k_d < vbx_lo_k) k_d = vbx_lo_k;
                        if (i_d < dlo_pc[0]) i_d = dlo_pc[0];
                        if (i_d > dhi_pc[0]) i_d = dhi_pc[0];
                        if (j_d < dlo_pc[1]) j_d = dlo_pc[1];
                        if (j_d > dhi_pc[1]) j_d = dhi_pc[1];
                        // P2-2b deliverable 4, law path -- the SECOND branch
                        // of this site. ⚠ Both branches must convert together:
                        // leaving one on the scalar form would make the
                        // diagnostic disagree with itself depending only on
                        // whether the sample fell inside the top half-cell.
                        // P2-2e: the drive's reference -- T0 at the sampled
                        // cell at d = 0 (bit-identical); under d > 0 the
                        // gathered column (top face) or the profile model
                        // (lateral face). ⚠ The SAMPLE T_a(i_d, j_d, k_d) is
                        // unchanged: the stress formula is not, only its
                        // reference (packet guardrail).
                        const Set::Scalar Tref_d = Tref_at(depth, i_d, j_d, k_d);
                        if (cref_on_sp && T_a(i_d, j_d, k_d) < Tref_d) ++sp_cref_neg;
                        if (te_g)
                            return te_relief
                                 * LawSigma(te_g, T_a(i_d, j_d, k_d),
                                            Tref_d)
                                 + sp_conf_offset;
                        // ⚠ pre-P2-2b expression, VERBATIM.
                        return te_relief * te_E * te_beta
                               * (T_a(i_d, j_d, k_d) - Tref_d)
                               / (1.0 - te_nu) + sp_conf_offset;
                    }
                    int k_d = ksurf - static_cast<int>(std::floor(dv / dz));
                    if (k_d > ksurf)    k_d = ksurf;
                    if (k_d < vbx_lo_k) k_d = vbx_lo_k;
                    Set::Matrix sigma =
                        Numeric::Interpolate::NodeToCellAverage(sig_node, i, j, k_d, 0);
                    return -sigma(0, 0) + sp_conf_offset;
                };

                // P2-2 change 3, READ 1 of 2 (the VERTICAL flaw). The size term
                // is applied HERE, not in InitializeFlaws, because `ell_q` is a
                // per-column FLUX length and init runs once. Identity under the
                // default `weibull.size_volume = cell`.
                // ⚠ `ell_q` comes from the column's own reconstruction
                // (surf_ellq_col, built in BuildSurfaceTColumns); if it is not
                // available the factor stays 1 rather than guessing.
                Set::Scalar a_loc = use_weibull ? flaw_a(i, j, k) : a0_loc;
                if (sp_size_flux)
                {
                    const std::size_t c_sz =
                        static_cast<std::size_t>(col_index_sp(i, j));
                    const Set::Scalar lq = (c_sz < surf_ellq_col.size())
                        ? surf_ellq_col[c_sz] : Set::Scalar(0.0);
                    if (lq > 0.0) a_loc *= WeibullSizeFactor(lq * lq * lq);
                }

                Set::Scalar K_I_val =
                    Numeric::SpCriterion::K_I(sigma_at_depth, a_loc, n_seg);

                // P0 deliverable 3: K_Ic AT THE CRACK TIP, not at the cell.
                // The tip of the surface-breaking flaw sits at depth a_loc
                // (normal), i.e. vertical depth a_loc*icth, where the rock is
                // COLDER than the surface and therefore TOUGHER. With the key
                // off this is the cell average, unchanged.
                const Set::Scalar T_kic =
                    rec_ok ? T_recon_at(a_loc * icth) : Tsurf;

                Set::Scalar K_Ic_MPa;
                if (law == KIcLaw::GrimselTable)
                    // P2-2d d7, toughness site 1 of 4. ⚠ kic_scale_ multiplies
                    // uniformly here exactly as it does for Nasseri, so the
                    // R/scale degeneracy `(a2)` proves is unchanged by route.
                    K_Ic_MPa = kic_scale_ * kic_g->K_Ic_of(T_kic);
                else if (law == KIcLaw::NasseriTable)
                    K_Ic_MPa = Numeric::SpCriterion::KIcNasseri(T_kic, kic_scale_);
                else
                    K_Ic_MPa = Numeric::SpCriterion::KIcLinear(
                        T_kic, K_Ic0_, T_ref_kic_, T_melt_);

                Set::Scalar K_Ic_Pa = Numeric::SpCriterion::MPa_sqrtm_to_Pa_sqrtm(K_Ic_MPa);
                Sp_a(i, j, k) = (K_Ic_Pa > 0.0) ? K_I_val / K_Ic_Pa : 0.0;
            });
        }

        // P2-2e: this site's share of the constraint-reference witnesses.
        // Hand-reduced (extensive = false is not reduced by the Integrator), and
        // ADDED to the members BuildConstraintRefColumns zeroed this step --
        // Removal.H adds the depth scan's and the lateral criterion's share.
        if (cref_on_sp)
        {
            long cnt[4] = {sp_cref_vert, sp_cref_lat, sp_cref_clamp, sp_cref_neg};
            amrex::ParallelDescriptor::ReduceLongSum(cnt, 4);
            cref_n_vert  += static_cast<Set::Scalar>(cnt[0]);
            cref_n_lat   += static_cast<Set::Scalar>(cnt[1]);
            cref_n_clamp += static_cast<Set::Scalar>(cnt[2]);
            cref_n_neg   += static_cast<Set::Scalar>(cnt[3]);
        }

        // P0 witnesses. Hand-reduced (extensive = false is not reduced by the
        // Integrator) and written on every rank so thermo.dat is identical for
        // any rank count.
        if (recon_on)
        {
            long cnt[7] = {rec_n, rec_n_flat, rec_n_ellclamp, rec_n_Tclamp,
                           rec_n_fallback, rec_n_hfcap, rec_n_hfbnd};
            amrex::ParallelDescriptor::ReduceLongSum(cnt, 7);
            amrex::ParallelDescriptor::ReduceRealSum(rec_sum_dT);
            recon_cols        = static_cast<Set::Scalar>(cnt[0]);
            recon_flat        = static_cast<Set::Scalar>(cnt[1]);
            recon_ell_clamped = static_cast<Set::Scalar>(cnt[2]);
            recon_T_clamped   = static_cast<Set::Scalar>(cnt[3]);
            recon_fallback    = static_cast<Set::Scalar>(cnt[4]);
            recon_hf_capped   = static_cast<Set::Scalar>(cnt[5]);
            recon_hf_bounded  = static_cast<Set::Scalar>(cnt[6]);
            recon_dT_mean     = (cnt[0] > 0)
                                ? rec_sum_dT / static_cast<Set::Scalar>(cnt[0])
                                : 0.0;
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
        if (!SurfaceOwner(lev)) return;   // AMR: clusters live on the surface-owning level

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
        if (SurfaceOwner(lev) && cluster_count > 0 && sp_emit_clusters_csv &&
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
    // ⚠ `static` and taking the law by POINTER on purpose. The `Sp_field`
    // diagnostic site sits inside `amrex::ParallelFor` under a `[=, this]`
    // capture, where the per-column vectors and callables of the host path
    // cannot go (see Removal.H's note at the lateral criterion). A pointer
    // hoisted before the loop crosses that boundary; a 15.7 KB by-value copy
    // of the two property tables would be carried into every closure.
    //
    // ⚠⚠ THE SIGN. Returns COMPRESSION-POSITIVE for T_hot > T_cold, matching
    // `SpCriterion::K_I`'s convention and the double-negation the pre-P2-2b
    // expressions perform. `Removal.H:1833` returns `-sg(0,0)` for the same
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
    //
    // D2q-1: this is the ONE site that could NOT be made rem-authoritative
    // without changing behaviour, so it is left on phi and ESCALATED rather
    // than approximated. Replacing the phi window with "solid here, void
    // directly above" changed 24 plotfiles in `dev2d` -- and changed the
    // trajectory, not merely the written `surface` field: the differing files
    // migrate between the 40000cell / 48000cell / 56000cell output
    // directories. `input_2d_dev` sets no `surface.follow_mask`, so this mask
    // is a live decision path feeding SurfaceCellFlux, not a diagnostic.
    // See studies/d2q1_percellvoid/RESULTS.md section "Site 7".
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

    // D2q-2b items 6 + 8: LATERAL per-cell removal.
    //
    // Runs after UpdateSpAfterMechanics, which (with the key on) has already
    // written Sp for every exposed cell, scored along that cell's own face
    // normal. This pass only decides and flips.
    //
    // SCOPE -- and this is what makes criterion (f) exact rather than
    // approximate. Only cells that are NOT their column's top and that have an
    // exposed LATERAL face are considered here. The column top keeps its
    // existing path entirely (depth scan, flake coherence, h_col remainder).
    // A flat interior therefore has no candidate cells at all, so with the key
    // on its rate is not merely close to the key-off rate, it is produced by
    // the identical code. Anything (f) detects is a real difference, not drift.
    //
    // MPI CONSISTENCY -- no precedent in this codebase, so stated explicitly:
    //  * The criterion is LOCAL. At a0 = 20 um against dx = 2 mm every sample
    //    in the K_I integral lands in the scored cell (floor(dv/dz) = 0), so Sp
    //    is a function of that cell's own T and a_f. No column gather is needed
    //    -- unlike the top-cell path, where a column spans ranks and k_top must
    //    be reduced.
    //  * The only non-local read is the six-neighbour void test, taken from
    //    removed_mf's one-ghost layer refilled inside UpdateSpAfterMechanics.
    //  * Decide-then-flip. Sp was computed from the FROZEN pre-pass state, and
    //    this pass reads only Sp (an owned cell) when flipping. So no cell's
    //    fate depends on whether a neighbour was flipped first, and the result
    //    is independent of box decomposition and iteration order. Flipping in
    //    place while re-testing exposure would have made a rank's sweep order
    //    observable -- the same class of bug as the 0-ghost removed_mf seam
    //    defect found earlier in this project.
    //  * Ghosts are refilled after the flip, and the counters are reduced.
    void PerCellRemoval(int lev, Set::Scalar /*time*/, Set::Scalar /*dt*/)
    {
        pcr_lateral_cells  = 0.0;
        pcr_overhang_cells = 0.0;
        if (!spall_per_cell_removal || !SurfaceOwner(lev)) return;

        const amrex::Box& domain = geom[lev].Domain();
        const int dlo[3] = {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
        const int dhi[3] = {domain.bigEnd(0),   domain.bigEnd(1),   domain.bigEnd(2)};
        const amrex::Vector<int> k_top = ColumnTopSolid(lev);
        const int nx = domain.length(0), ilo = dlo[0], jlo = dlo[1];
        // Mirrors Removal.H:125 -- regime_field_mf exists only on this branch.
        const bool sp_weibull_spall =
            spall_enabled && spallation_model == SpallationModel::SpWeibull;

        Set::Scalar n_fired = 0.0;
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar>       const& rem = (*removed_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& Sp  = (*Sp_field_mf[lev]).const_array(mfi);
            amrex::Array4<Set::Scalar>       const& reg = (*regime_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& D   = (*D_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> regf;
            if (sp_weibull_spall) regf = (*regime_field_mf[lev]).array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (rem(i, j, k) > 0.5) continue;
                const int kt = k_top[(j - jlo) * nx + (i - ilo)];
                if (k == kt) continue;                  // the column path owns the top
                const bool lat =
                    ((i > dlo[0]) && rem(i - 1, j, k) > 0.5) ||
                    ((i < dhi[0]) && rem(i + 1, j, k) > 0.5) ||
                    ((j > dlo[1]) && rem(i, j - 1, k) > 0.5) ||
                    ((j < dhi[1]) && rem(i, j + 1, k) > 0.5);
                if (!lat) continue;
                if (!(Sp(i, j, k) >= 1.0)) continue;
                rem(i, j, k) = 1.0;
                reg(i, j, k) = 1.0;                     // spall regime
                if (sp_weibull_spall) regf(i, j, k) = 1.0;
                D(i, j, k) = 0.0;
                // H and T are deliberately NOT reset -- see Removal.H's note:
                // a void cell is insulating, its stored state is inert, and
                // resetting it both breaks the energy ledger's neutrality and
                // produces a spurious cold checkerboard.
                n_fired += 1.0;
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(n_fired);
        pcr_lateral_cells = n_fired;

        // The flip changed the void mask; refill before anything reads it.
        Util::RealFillBoundary(*removed_mf[lev], geom[lev]);

        // Item 8: overhang census -- void cells with SOLID directly above.
        // A one-cell fillet at the corner is expected. A DEEP undercut is a
        // finding to report, not a bug to suppress: nothing in this model
        // collapses an overhang, so whatever forms, stays.
        Set::Scalar n_over = 0.0;
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
                if (rem(i, j, k) > 0.5 && k < dhi[2] && rem(i, j, k + 1) <= 0.5)
                    n_over += 1.0;
        }
        amrex::ParallelDescriptor::ReduceRealSum(n_over);
        pcr_overhang_cells = n_over;
    }

    // D2q-2b item 0: flip the seeded index-space box to void, once, right
    // after InitPhi. Sets `removed` ONLY -- phi is deliberately left alone, so
    // the two genuinely disagree and every consumer is forced through the
    // authoritative mask that D2q-1 installed. That disagreement is the whole
    // point of the probe: with phi still saying "solid", any path that quietly
    // fell back to phi shows up immediately.
    //
    // No-op unless `spall.seed_void` was given.
    void ApplySeedVoid(int lev)
    {
        if (spall_seed_void.size() < 6 || spall_seed_void.size() % 6 != 0)
            return;
        long n = 0;
        for (std::size_t sb = 0; sb < spall_seed_void.size(); sb += 6)
        {
            const amrex::Box seed(
                amrex::IntVect{AMREX_D_DECL(spall_seed_void[sb + 0],
                                            spall_seed_void[sb + 1],
                                            spall_seed_void[sb + 2])},
                amrex::IntVect{AMREX_D_DECL(spall_seed_void[sb + 3],
                                            spall_seed_void[sb + 4],
                                            spall_seed_void[sb + 5])});
            for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & seed;
                if (!bx.ok()) continue;
                amrex::Array4<Set::Scalar> const& rem = (*removed_mf[lev]).array(mfi);
                const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                    if (rem(i, j, k) <= 0.5) { rem(i, j, k) = 1.0; ++n; }
            }
        }
        amrex::ParallelDescriptor::ReduceLongSum(n);
        Util::Message(INFO, "spall.seed_void: level ", lev, ", voided ", n, " cells");
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
        // Static-field ghost layers must be refilled on the next advance.
        if (lev >= static_cast<int>(micro_ghost_filled.size()))
            micro_ghost_filled.resize(lev + 1, 0);
        micro_ghost_filled[lev] = 0;
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
                // ⚠⚠ P2-2b SITE 4, resolution (ii): BUILD FROM THE LAW, never
                // from the scalar. `Table::Build` takes `const Material&` and is
                // polymorphic (`Table.H:30`), so a Grimsel object goes straight
                // in and the table then integrates rho(T)*c_p(T) instead of a
                // constant -- which is what routes rho and c_p into conduction
                // AND into the H->T inversion, from ONE curve.
                //
                // ⚠⚠ THE LATENT HEAT MUST BE COPIED ACROSS FIRST. Grimsel parses
                // only `T_ref`, so `T_m`/`fus_width`/`L_m`/`T_vap_lo`/`T_vap_hi`/
                // `L_v` sit at Material's base defaults and
                // `effective_rhocp`'s boxcars would contribute NOTHING -- the
                // latent heat vanishes silently. ⚠ Inert on both production
                // decks (neither sets a phase-change scalar), so the ledger
                // cannot see it; `unit/latent_heat` is the only witness.
                if (phases[p].grimsel_law)
                {
                    Numeric::Material::Grimsel g = phases[p].grimsel;
                    g.T_ref     = phases[p].T_ref;
                    g.T_m       = phases[p].T_m;
                    g.fus_width = phases[p].fus_width;
                    g.L_m       = phases[p].L_m;
                    g.T_vap_lo  = phases[p].T_vap_lo;
                    g.T_vap_hi  = phases[p].T_vap_hi;
                    g.L_v       = phases[p].L_v;
                    phase_tables[p].Build(g);
                }
                else
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
            // ⚠⚠ P2-2b SITE 5 -- THE ONE LEGITIMATE SCALAR -> FIELD WRITER, and
            // the ONLY place in the rock path that may read a per-phase scalar.
            // Every consumer reads the FIELD; this is where per-phase values
            // become per-cell ones.
            //
            // ⚠⚠ DO NOT DELETE THIS IN FAVOUR OF `UpdatePropertyFields` ALONE.
            // That function returns at `:5380` unless SOME phase sets
            // `law = grimsel`, and skips non-Grimsel phases at `:5412`. Under the
            // shipped default the fields would then be NEVER WRITTEN, and the
            // idle clock divides by `pin_rho*pin_Cp`. Together the two cover
            // every phase under every law: this at init, that per step.
            //
            // ⚠ LAW-AWARE so the fields hold the law's values at EVERY instant,
            // not only from the first Advance. Without this a grimsel deck
            // carries scalar values between Initialize and the first
            // UpdatePropertyFields -- a window nothing reads today, which is
            // exactly the kind of latent inconsistency that later becomes a bug.
            if (phases[p].grimsel_law)
            {
                rh[p]  = phases[p].grimsel.rho(phases[p].T_ref);
                Cpv[p] = phases[p].grimsel.cp (phases[p].T_ref);
            }
            else
            {
                rh[p]  = phases[p].rho;
                Cpv[p] = phases[p].Cp;
            }
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
        // P2-2 change 3: under `size_volume = flux` the size term is NOT baked
        // in here -- it is applied per face, per step, at the flaw READS, from
        // the flux length. Under the default it is the old factor, which is
        // exactly 1 on both decks anyway.
        const Set::Scalar vol_factor =
            sp_size_flux ? Set::Scalar(1.0)
                         : std::pow(V_cell / weibull_V0, inv_m);

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
    // AMR surface band (m): cells kept on the finest level below / above each
    // column's top solid cell. depth > 0 is required when max_level > 0.
    Set::Scalar surface_refine_depth  = 0.0;
    Set::Scalar surface_refine_height = 0.0;

    // S1: whole-domain energy ledger (EnergyLedger), thermo.dat vars.
    // E_H [J] is sum(H*V) over all level-0 cells (void included), E_src [J]
    // the cumulative applied energy, P_* [W] per step. The compensated-sum
    // state (ledger_E_H0_s/c) is internal.
    bool        ledger_enabled               = false;
    // C1: pinned-surface closure state (lev 0, ncols-sized host vectors,
    // identical on every rank). pin_T = top-cell T at the column's last spall
    // firing (0 = never fired), pin_t = that time, pin_rho / pin_Cp = the
    // per-cell rho and c_p (carried SEPARATELY -- P2-2b's rule) of the
    // top-cell phase then. pin_T_eff = the pin used by this step's H update
    // and ledger (0 = face form), built at the start of AdvanceMicrostructure.
    // NOT checkpointed: a restart falls back to face until each column fires.
    std::vector<Set::Scalar> pin_T, pin_t, pin_rho, pin_Cp, pin_T_eff;
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
    // D2c: march options and burner-body clearance. Host state (jet_fs_z0,
    // foot_clear_k) is NOT checkpointed, like the D2b state above: a restart
    // re-takes the original surface from the restart step's faces.
    Set::Scalar jet_D_e      = 0.0; // thermo: decay diameter used this step [m] (momentum)
    Set::Scalar jet_fs_cols  = 0.0; // thermo: marching free-surface columns
    Set::Scalar jet_m_ratio  = 0.0; // thermo: m_exit/m_stag
    Set::Scalar jet_T_mouth  = 0.0; // thermo: inlet T_gas of the first free-surface bin [K]
    Set::Scalar jet_P_neg    = 0.0; // thermo: sum of negative q_robin*dA [W] (count)
    std::vector<Set::Scalar> jet_fs_z0;   // per column face height at the first march [m]
    Set::Scalar foot_mech_vol  = 0.0;     // thermo: cumulative mechanically removed volume [m^3]
    // P2-2f: cumulative enthalpy that left the SOLID set with pad-cleared rock
    // [J]. Credited to ledger_err_solid; kept OUT of ledger_E_flake so the
    // heat-driven price (floor + wall) excludes mechanical removal on both
    // sides of the ratio. See Removal.H's pad-clearance block for why this is a
    // booking-only term and why the field must not be written.
    Set::Scalar foot_mech_E    = 0.0;     // thermo: cumulative pad-clearance enthalpy [J]
    Set::Scalar foot_mech_cols = 0.0;     // thermo: columns cleared this step
    std::vector<int> foot_clear_k;        // per column target top k (INT_MIN = none)
    // D2k: side-face flux on exposed vertical faces (surface_patch.side_face_flux).
    Set::Scalar side_P_face = 0.0;  // thermo: net power into side faces [W] (f applied)
    // P1e diagnostics, registered ONLY under surface_patch.side_face_staircase.
    Set::Scalar side_sc_mean = 1.0;  // area-weighted mean staircase factor
    Set::Scalar side_sc_min  = 1.0;  // smallest factor seen this step
    Set::Scalar side_sc_cols = 0.0;  // DISTINCT columns with a corrected face
    Set::Scalar side_P_in35  = 0.0;  // share of side_P_face landing at r < 35 mm
    // D2r-1 witness: side-face power delivered to columns ABOVE the
    // nozzle plane, and the face count. Registered ONLY when
    // side_face_factor_above is set, so unset runs keep their exact
    // thermo.dat -- which is what keeps the frozen A_f02 comparable.
    Set::Scalar side_P_face_above = 0.0;
    Set::Scalar side_w_mean = 0.0;  // D2x thermo: mean realised riser weight

    // D2z: per-column recession PHASE, maintained off the existing k_top_col
    // mask (side_face_flux = 1 forces surface.follow_mask = 1, so the mask is
    // always built). k_top_col is a domain-wide REDUCED vector identical on
    // every rank, so these are rank-consistent and regrid-safe by construction.
    //
    // ⚠ NOT the stored sub-cell remainder: under whole-cell removal phi_top is
    // single-valued (dz/2 - O(1e-16), drifting ~1.6e-16 over a run), so it
    // cannot express a lag at all. The phase must be TIME-based.
    amrex::Vector<int>         col_ktop_prev;   // last step's k_top per column
    amrex::Vector<Set::Scalar> col_t_last;      // time of that column's last void
    amrex::Vector<Set::Scalar> col_dt;          // its TRAILING flip interval
    amrex::Vector<int>         col_nvoid;       // voids so far
    std::size_t                step_split_nx = 0;  // domain nx, for the Step 3b assertion
    // Thermo, registered only under side_face_step_height (Goal 5).
    Set::Scalar delta_mean          = 0.0;  // TIME-WEIGHTED over the live faces
    Set::Scalar delta_fallback_cols = 0.0;
    Set::Scalar riser_area_ratio    = 0.0;  // assigned area / the other law's

    // phi_i = clamp((t - t_last_i)/dt_i, 0, 1), from the column's own TRAILING
    // interval. Fewer than two voids => no interval => phi = 1, so delta falls
    // back to the discrete cell count and g = 1, i.e. exactly today's
    // behaviour. The fallback is NAMED (delta_fallback_cols), never silent.
    AMREX_FORCE_INLINE
    static Set::Scalar ColPhi(const Set::Scalar* t_last, const Set::Scalar* dt_i,
                              const int* nvoid, int c, Set::Scalar time, int& fb)
    {
        if (nvoid[c] < 2 || !(dt_i[c] > 0.0)) { fb = 1; return 1.0; }
        const Set::Scalar p = (time - t_last[c]) / dt_i[c];
        return p < 0.0 ? Set::Scalar(0.0) : (p > 1.0 ? Set::Scalar(1.0) : p);
    }

    // delta = |n + phi_S - phi_T|, n = k_top[T] - k_top[S] >= 1 (S is the
    // SHORTER/deeper column, T the taller). BOTH phases -- reading only the
    // lower column's returns ~0 by construction, because a step is born the
    // instant the deeper column flips.
    AMREX_FORCE_INLINE
    static Set::Scalar PairDelta(const int* ktc, const Set::Scalar* t_last,
                                 const Set::Scalar* dt_i, const int* nvoid,
                                 int cT, int cS, Set::Scalar time, int& fb)
    {
        const int n = ktc[cT] - ktc[cS];
        const Set::Scalar pT = ColPhi(t_last, dt_i, nvoid, cT, time, fb);
        const Set::Scalar pS = ColPhi(t_last, dt_i, nvoid, cS, time, fb);
        const Set::Scalar d = Set::Scalar(n) + pS - pT;
        return d > 0.0 ? d : Set::Scalar(0.0);
    }

    // The heating area, in FOOTPRINTS (multiply by dx*dy).
    //   step   : delta                       -- a real vertical face
    //   wetted : sqrt(1 + delta^2) - 1       -- the EXACT ramp, which is
    //            RiserWeight's own formula with the continuous delta in place
    //            of the integer cell count. NOT delta^2/2 (small-step limit).
    // D2z Step 3b: how delta splits over the cells of an n-cell riser.
    // The taller column's exposed cells run from k_top(S)+1 to k_top(T), which
    // is exactly n cells. DERIVED, then asserted (AssertStepSplit):
    //
    //   n == 1 : the single cell carries the whole lag,  1 - phi_T + phi_S
    //   n >= 2 : top      (k == kT)        -> 1 - phi_T
    //            lowest   (k == kS + 1)    -> 1 + phi_S   <- phi_S folded in
    //            middle                    -> 1
    //   sum    = (1-phi_T) + (n-2) + (1+phi_S) = n - phi_T + phi_S = delta  ✓
    //
    // ⚠ The packet's split ((1-phi_T), 1 per middle, phi_S on the lowest) sums
    // to delta - 1: the phi_S piece geometrically belongs to a cell with NO
    // void neighbour, so it must be FOLDED INTO the lowest exposed cell rather
    // than occupy one of its own.
    AMREX_FORCE_INLINE
    static Set::Scalar StepCellWeight(int n, int k, int kT, int kS,
                                      Set::Scalar phiT, Set::Scalar phiS)
    {
        if (n <= 1) return 1.0 - phiT + phiS;
        if (k == kT)     return 1.0 - phiT;
        if (k == kS + 1) return 1.0 + phiS;
        return 1.0;
    }

    AMREX_FORCE_INLINE
    // ⚠⚠ P1e: THE STAIRCASE ORIENTATION CORRECTION.
    //
    // A staircase approximation of an inclined wall exposes MORE face area than
    // the wall has: on a Cartesian grid a round wall's diagonal column presents
    // two open faces where an axis column presents one, so the diagonals collect
    // ~sqrt(2) times the side heat for the same physical wall. The standard
    // correction is the ratio of the true face area to the staircase area, which
    // for a surface normal n in the horizontal plane is
    //
    //     factor = 1 / (|n_x| + |n_y|)
    //
    // and the normal's components are proportional to the summed RISER HEIGHTS
    // of the open faces on each axis. Hence hypot(hx, hy) / (hx + hy).
    //
    // ⚠ It is a DISCRETISATION FIX, NOT A THRESHOLD. Exactly 1 on a grid-aligned
    // wall, exactly 1 on a thin pillar (an opposing pair is not a sloped wall),
    // exactly 1/sqrt(2) on a 45 degree corner, and it converges with mesh.
    // Nothing in it is fitted, and by construction it lies in [1/sqrt(2), 1] for
    // any non-negative hx, hy -- so NO CLAMP IS NEEDED. If a clamp ever looks
    // necessary the inputs are wrong.
    //
    // ⚠ The arguments are RISER CELL COUNTS, not metres. The ratio is
    // scale-invariant -- dz cancels exactly between numerator and denominator --
    // so taking integers both removes a dependency on dz at the call sites and
    // avoids a rounding step. A face that is not open passes 0.
    static Set::Scalar StaircaseFactor(int n_xm, int n_xp, int n_ym, int n_yp)
    {
        const Set::Scalar hx = Set::Scalar(n_xm + n_xp);
        const Set::Scalar hy = Set::Scalar(n_ym + n_yp);
        const Set::Scalar sum = hx + hy;
        // No open face, or risers of zero height: there is no orientation to
        // correct, so the factor is exactly 1 and the area is untouched.
        if (!(sum > 0.0)) return 1.0;
        return std::hypot(hx, hy) / sum;
    }

    static Set::Scalar AreaOfDelta(Set::Scalar delta, bool wetted)
    {
        return wetted ? (std::sqrt(1.0 + delta * delta) - 1.0) : delta;
    }

    // g(delta, n, a) = (1 - e^{-delta a}) / (1 - e^{-n a}),  a = dz*v/alpha.
    //
    // D2z-2 deliverable 4. The factor corrects the GRID's lateral flux across a
    // step to what the TRUE step would pass. Under the quasi-steady recession
    // profile the grid's own drop across that face goes as 1 - e^{-n a},
    // because the grid has ALREADY RESOLVED the step as n cells; the true drop
    // goes as 1 - e^{-delta a}. The correction is their ratio.
    //
    // ⚠ D2z normalised by 1 - e^{-a} -- i.e. to n = 1 -- REGARDLESS of n, so for
    // n >= 2 it double-counted the depth the grid had already resolved and
    // AMPLIFIED the flux, up to 1/(1 - e^{-a}). That was an arithmetic error,
    // not a modelling choice, so there is no key and nothing to A/B:
    //
    //     a     n  delta   D2z form   derived
    //   1.29    2    2.3     1.31      1.026
    //   0.65    2    2.3     1.62      1.066
    //   0.65    3    3.2     1.83      1.020
    //   0.65    2    1.8     1.44      0.948   <- a CLAMP to 1 would be wrong
    //                                             in the other direction
    // Three properties, each verified before this was written:
    //   * n == 1 reproduces the D2z form to machine precision, so nothing
    //     measured inside the arbiter's r < 25 mm (0.0 % multi-cell) changes.
    //   * a -> 0 degenerates to delta/n, NOT min(delta, 1).
    //   * the never-fired fallback (phi = 1 on both columns => delta = n)
    //     returns EXACTLY 1 for any a, with no special case. The D2z form
    //     returned 1.67 at n = 3 and 2.90 at n = 7 (each at its own mesh's a).
    AMREX_FORCE_INLINE
    static Set::Scalar StepConductance(Set::Scalar delta, int n, Set::Scalar a)
    {
        const Set::Scalar nn = Set::Scalar(n > 1 ? n : 1);
        // a -> 0: (1 - e^{-delta a})/(1 - e^{-n a}) -> delta/n.
        if (!(a > 1.0e-12)) return delta < nn ? delta / nn : Set::Scalar(1.0);
        const Set::Scalar den = 1.0 - std::exp(-nn * a);
        if (!(den > 0.0)) return 1.0;
        Set::Scalar g = (1.0 - std::exp(-delta * a)) / den;
        return g < 0.0 ? Set::Scalar(0.0) : g;
    }
    Set::Scalar side_cols_above   = 0.0;
    Set::Scalar wall_step_q       = 0.0;  // mean k*dT/dx across the wall step [W/m^2]
    Set::Scalar wall_step_q_eff   = 0.0;  // D2z: the same mean under g(delta)
    Set::Scalar wall_step_n       = 0.0;  // faces in that mean
    // D2z-2: the same three over the step_diag_radius population, plus the
    // conductance the meter actually applied and the multi-cell fraction.
    Set::Scalar wall_step_n_in     = 0.0;
    Set::Scalar wall_step_q_in     = 0.0;
    Set::Scalar wall_step_q_eff_in = 0.0;
    Set::Scalar step_g_mean        = 0.0;
    Set::Scalar step_g_max         = 0.0;
    Set::Scalar wall_step_f_multi  = 0.0;
    Set::Scalar step_delta_faces   = 0.0;  // stepped FACES behind delta_mean
    Set::Scalar side_P_face_multi  = 0.0;  // side_P_face drawn by n >= 2 faces
    // D2aa: charged lateral area [m^2] and side power [W] by radial band.
    // D2ab deliverable 1: the wall gas temperature ACTUALLY applied this step.
    // Resolved once in BuildFlameColumns; the kernel, the march and the ledger
    // all read this member and never surface_patch.side_face_T. Host state,
    // NOT checkpointed, exactly like jet_T_rec.
    Set::Scalar side_T_eff  = -1.0; // thermo [K]
    Set::Scalar side_A_bore = 0.0, side_A_pads = 0.0, side_A_wall = 0.0;
    Set::Scalar side_P_bore = 0.0, side_P_pads = 0.0, side_P_wall = 0.0;
    Set::Scalar wall_step_f_multi_pads = 0.0;  // multi-cell fraction, pad ring
    Set::Scalar side_P_loss = 0.0;  // thermo: side-face surface losses [W] (f applied)
    Set::Scalar side_faces  = 0.0;  // thermo: number of exposed lateral faces heated

    // ---- P2-1: THE DERIVED WALL COEFFICIENT -------------------------------
    // Per-step state for the composite wall law. Computed ONCE at the top of
    // Advance (UpdateWallLawState) and read by all four roles, so the
    // criterion, the kernel, the march and the ledger cannot disagree about
    // which law ran -- which is the whole point of the step and what
    // unit/lateral_flux_lock's P2-1 block asserts.
    //
    // ⚠ Host state, NOT checkpointed, exactly like side_T_eff and jet_T_rec.
    // ⚠ The gas state is side_T_eff, which BuildFlameColumns writes AFTER the
    //   criterion runs, so at the criterion it is one step stale. That
    //   staleness is PRE-EXISTING (the criterion already heats with that same
    //   side_T_eff) and the lock declares T_gas equal across the four sites on
    //   that basis; the law inherits it and adds none of its own.
    Set::Scalar wall_uDe     = 0.0; // u*D_e, the jet's momentum product [m^2/s]
    Set::Scalar wall_rho     = 0.0; // gas density at the wall gas state [kg/m^3]
    Set::Scalar wall_nu      = 0.0; // kinematic viscosity [m^2/s]
    Set::Scalar wall_cp      = 0.0; // gas c_p [J/kgK] (surface_patch.jet_cp)
    Set::Scalar wall_h_duct  = 0.0; // the duct branch [W/m^2K]
    Set::Scalar wall_r_bore  = 0.0; // measured bore radius [m]
    // P2-3: the radius the height datum actually truncated its ray at, so the
    // datum's inner edge is reported rather than inferred from a key.
    Set::Scalar wall_datum_r = 0.0; // [m]
    Set::Scalar wall_sector  = 0.0; // simulated azimuthal fraction of the disc
    Set::Scalar wall_dz      = 0.0; // the cell height the datum counts in [m]
    bool        wall_law_ready = false;
    // k_corner per column, SAME indexing as ColumnTopSolid: (j-jlo)*nx+(i-ilo).
    amrex::Vector<int> wall_corner;
    // P2-3 deliverable 2: the wall march's per-bin gas temperature, indexed by
    // (k - wall_k0) where wall_k0 is the LOWEST above-plane face level this
    // step. ⚠ Host state, NOT checkpointed, exactly like side_T_eff.
    std::vector<Set::Scalar> wall_T_gas;
    int wall_k0 = 0;
    // P2-3 thermo, ALL PER-STEP STATE -- average them in a window, never sum
    // (standing rule 7). `jet_P_wall` is the enthalpy the stream lost to the
    // wall THIS step; `wall_T_mouth` the gas temperature after the last wall
    // bin; `wall_bins` the number of height bins; `wall_L_fit` the cooling
    // length fitted to wall_T_gas[k].
    Set::Scalar jet_P_wall   = 0.0;   // [W]
    Set::Scalar wall_T_mouth = 0.0;   // [K]
    Set::Scalar wall_bins    = 0.0;   // [-]
    Set::Scalar wall_L_fit   = 0.0;   // [m]
    // P2-4, all PER-STEP STATE except body_noLaw_steps (cumulative).
    Set::Scalar jet_P_body       = 0.0;  // [W] the cold body's debit this step
    Set::Scalar body_h_used      = 0.0;  // [W/m^2K] the duct scalar actually used
    Set::Scalar body_z_imm       = 0.0;  // [m] immersed length = wall_bins * dz
    Set::Scalar body_L           = 0.0;  // [m] mcp_wall / (h * P_q)
    Set::Scalar jet_P_retained   = 0.0;  // [W] entrained share, `fresh` only
    Set::Scalar body_noLaw_steps = 0.0;  // CUMULATIVE count -- difference it
    // P2-4: the regime-6 share of ledger_E_flake. CUMULATIVE, like its parent.
    Set::Scalar ledger_E_flake_wall = 0.0;  // [J]
    // thermo, registered ONLY under side_face_h_law = derived so the typed
    // default adds no column and the frozen tree is untouched.
    Set::Scalar side_h_used_min = 0.0, side_h_used_p05 = 0.0;
    Set::Scalar side_h_used_p50 = 0.0, side_h_used_p95 = 0.0;
    Set::Scalar side_h_used_max = 0.0, side_h_used_n   = 0.0;
    // P2-3 (a5): of the faces the DERIVED law reached, how many were ABOVE the
    // nozzle plane -- i.e. hole wall rather than floor cone flank. ⚠ This is
    // the witness the packet asked `side_h_used_n` for; that column is
    // band-wide and reads 8968 on a control leg whose above-plane power is
    // exactly 0, so it cannot discriminate.
    Set::Scalar side_h_used_n_above = 0.0;
    Set::Scalar wall_h_corner   = 0.0;  // the law at z = 0, r = r_bore
    Set::Scalar wall_z_fill     = 0.0;  // the wall-jet/duct switch height [m]
    Set::Scalar wall_law_on     = 0.0;  // 1 while the law is actually applied
    // P2-2 change 1: how many times the property fields were refreshed.
    Set::Scalar prop_refresh_n  = 0.0;
    // (b5)'s two-temperature probe reads these: the MEAN refreshed
    // property over the level, beside the mean temperature it was
    // evaluated at. A refresh that never happens leaves them at 0.
    Set::Scalar prop_kappa_mean = 0.0, prop_beta_mean = 0.0;
    Set::Scalar prop_E_mean     = 0.0, prop_T_mean    = 0.0;
    // =====================================================================
    // P2-2 change 4: THE STRESS-RELIEF FACTOR R.
    // =====================================================================
    // `spallation.sp.stress_relief`, default 1.0 (identity, so every frozen
    // study stays byte-identical). It multiplies the thermoelastic driving
    // stress at ALL SEVEN sites.
    //
    // ⚠⚠ WHY A RELIEF FACTOR AT ALL. With the measured properties the heated
    // layer carries 1.78x this rock's own compressive strength at Meier's onset
    // and 1.90x at Kant's, so an elastic criterion is outside its validity and
    // something relieves it. R stands in for that, bounded by the measurement:
    // R = UCS/sigma = 0.563 (Meier) / 0.527 (Kant), one material number to 7 %.
    //
    // ⚠⚠ RELIEVING ONE SITE RELIEVES THE WALL AND NOT THE FLOOR -- the exact
    // family `Removal.H:589-594` pre-registers, which has cost this programme
    // three packets. unit/lateral_flux_lock asserts that SEVEN sites take
    // exactly this member, so no site can be missed and an eighth stress
    // expression cannot appear un-relieved.
    Set::Scalar sp_relief = 1.0;
    // =====================================================================
    // P2-2e: THE CONSTRAINT DEPTH.
    // =====================================================================
    // `spallation.sp.constraint_depth` [m], default 0 = today, byte-identical.
    // Under d > 0 the thermal DRIVE at depth z below a face is referenced to the
    // temperature a fixed material length d BEHIND it, T(z + d), instead of to
    // the stored initial temperature. Limits: thin skin -> today's criterion
    // (the banked level stands); soaked face -> drive ~ q d / k, so a MINIMUM
    // FLUX q_min ~ k dT_c / d emerges; uniform slab -> 0; interior-peaked
    // (beam) profile -> negative drive and no firing.
    // ⚠⚠ WHAT CHANGES IS THE FAILURE DRIVE, NOT THE STRESS. In a half-space
    // with T = T(z) linear elasticity gives a surface stress set by T_s - T_far
    // whatever the gradient, so today's stress is elastically right. Under
    // d > 0 every quantity the three sites still call "sigma" is the MISMATCH
    // drive across a surface-parallel crack of lateral size ~d; the name is
    // kept and documented at each site (ACTIVE_STEP P2-2e, advisor log §17).
    // ⚠ A MODE SWITCH, not a knob: d = 0 is today's absolute criterion; any
    // d > 0 is the differential, which -> 0 as d -> 0. Small d is NOT "nearly
    // today" (offline: no Kant standoff fires at d <= 3 mm on Gotthard).
    // ⚠ A ROCK PROPERTY, set ONCE on Kant 2016's SOD 130 / 150 bracket and
    // never adjusted afterwards.
    Set::Scalar sp_cref_d = 0.0;
    // The per-step GATHERED reference column, built ONCE per step BEFORE either
    // firing site runs (BuildConstraintRefColumns, called just ahead of
    // UpdateSpAfterMechanics) so the Sp_field site and the depth scan read the
    // SAME values in the SAME step -- the two-firing-site invariant. Never
    // allocated at d = 0.
    amrex::Vector<Set::Scalar> cref_T_col; // [c*cref_colw + m] = T(i,j,cref_ktop[c]-m)
    amrex::Vector<int>         cref_ktop;  // the column top the buffer is indexed from
    amrex::Vector<int>         cref_nuse;  // usable cells: band, domain floor, solid run
    int                        cref_colw = 0;
    // Thermo witnesses, registered ONLY under d > 0 and appended LAST, so a d = 0
    // leg gains no column and every existing column keeps its position.
    // PER-STEP: reset in BuildConstraintRefColumns and written as-is (`false`),
    // each site hand-reducing what its own per-box loop accumulated.
    Set::Scalar cref_d_out = 0.0, cref_n_vert = 0.0, cref_n_lat = 0.0,
                cref_n_clamp = 0.0, cref_n_neg = 0.0;
    // =====================================================================
    // P2-2 change 3: THE STRESSED VOLUME.
    // =====================================================================
    // `weibull.size_volume = cell | flux`, default `cell` (identity).
    //
    // ⚠⚠ `weibull.m` IS INERT TODAY, ON BOTH DECKS. The size term is
    // pow(V_cell/V0, 1/m) and it is evaluated ONCE, in InitializeFlaws, where
    // V_cell/V0 == 1 EXACTLY: Meier runs 2 mm cells with V0 = 8e-9 (set by
    // argv), Kant 24^3 on 0.024 m with the 1e-9 default. The factor is 1 for
    // ANY m, so the modulus changes nothing -- "m = 20 gives exponent 0.13" is
    // not merely wrong here, it gives ZERO.
    //
    // ⚠⚠ AND THE ONSET IS STRUCTURALLY FLUX-BLIND. The stressed volume is one
    // MESH CELL where it must be the THERMAL LAYER, so the criterion cannot see
    // the flux at all and no modulus can give it a flux exponent.
    //
    // Under `flux` the size term becomes (ell_q^3/V0)^(1/m) and is applied
    // where the flaw is READ -- per face, per step -- because `ell_q` is a
    // per-column flux length and `InitializeFlaws` runs once.
    //
    // ⚠ `ell_q = k(T_s - T_0)/q` is the FLUX length, never kappa/v: at ONSET
    // v = 0, so kappa/v is both singular and circular (it needs the rate the
    // criterion is trying to predict). P1g-3 established ell_q as an identity.
    bool sp_size_flux = false;
    // ⚠⚠ P2-2d deliverable 8: `V_ref` IS A SEPARATE KEY FROM `weibull.V0`.
    // `V0` stays the MESH-INDEPENDENCE anchor read under `cell`, where
    // V_cell/V0 == 1 exactly on both decks. `flux` needs a different reference
    // -- the stressed volume is ell_q^3, not a cell -- and sharing V0 for both
    // jobs is what makes the cross-deck 8x difference in V0 (Meier 8e-9, Kant
    // 1e-9) leak into a quantity that has nothing to do with the mesh.
    // ⚠ -1 means UNSET, and `flux` with it unset ABORTS rather than silently
    // borrowing V0: an anchor taken from the wrong reference is exactly the
    // class of silent-default error this line of packets keeps hitting.
    Set::Scalar weibull_V_ref = -1.0;

    // The size factor for a stressed volume V. Returns 1 exactly under the
    // default, which is what keeps every frozen anchor byte-identical.
    Set::Scalar WeibullSizeFactor(Set::Scalar V) const
    {
        if (!sp_size_flux) return 1.0;
        // ⚠ P2-2d d8: `flux` reads V_ref, NOT V0. Parse aborts if it is unset,
        // so reaching here with a bad value is not a routine path.
        if (!(V > 0.0) || !(weibull_V_ref > 0.0)) return 1.0;
        return std::pow(V / weibull_V_ref, 1.0 / weibull_m);
    }
    // Per-column side power at the current flame vectors, for the enthalpy
    // march's budget. Empty unless side_face_flux and jet_closure = enthalpy.
    std::vector<Set::Scalar> side_P_col;
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
        // P0b deliverable 1: `surface` is a FOURTH form, not a branch that
        // overrides the other three. Revision 1 selected it whenever
        // spall.surface_T was applied, which would have made the `face` and
        // `cell` control legs identical to the scored leg by construction and
        // destroyed the comparison the gate depends on.
        enum class RobinForm { Cell, Face, Pinned, Surface };
        // Revision 7 deliverable 6, `surface_patch.side_face_closure`: which closure a LATERAL face takes under
        // robin_form = surface. `follow` (default) is the revision-5 rule --
        // the face Newton on the hole wall, the owning column's T_surf on a
        // floor riser. `cell` forces the pre-revision-5 cell average back onto
        // lateral faces ONLY, leaving the top faces on the surface form, so
        // the wall closure's effect can be measured on one binary with one key
        // instead of by differencing two packets. ⚠ Inert unless
        // robin_form = surface, so cell/face/pinned stay byte-identical.
        enum class SideTForm { Follow, CellAvg };

        bool enabled = false;
        Mode mode = Mode::PrescribedTemperature;
        RobinForm robin_form = RobinForm::Cell;
        SideTForm side_T_form = SideTForm::Follow;
        Set::Scalar pinned_idle_cycles = 2.0; // C1: robin_form = pinned only
        // D2k: flux on exposed vertical faces of a staircase wall (default-off;
        // side_face_flux = 0 is byte-identical to every pre-D2k result).
        bool side_face_flux = false;
        Set::Scalar side_face_factor = 1.0;
        // D2q-2b: cap on a side face's heat-transfer coefficient, applied only
        // when spall.per_cell_removal = 1. See the parse block for why 80.
        Set::Scalar side_face_h_max = 80.0;
        // D2r-1: the side factor for columns ABOVE the nozzle plane (s <= 0).
        // SENTINEL -1 = "use side_face_factor", so every run that does not set
        // it is byte-identical. The drain flows into the taller column's buried
        // interior; that column's top is above the plane, so its h is ~0 and
        // its side faces have received 0 x f = 0 in every run ever done. This
        // key is what lets a leg heat exactly those faces and nothing else.
        Set::Scalar side_face_factor_above = -1.0;
        // P2-3 deliverable 4, REVISED (user, 2026-10-13). The wall law's HEIGHT
        // DATUM truncates its inward ray at a radius. Shipped behaviour takes
        // that radius from `side_face_h_r_min` -- a TYPED number, 0.040 m, while
        // the measured bore radius on the P2-2f m13 leg is 0.0397679 m. The bore
        // wall therefore sits 0.23 mm INSIDE the datum's own edge, about a tenth
        // of a cell, and ColumnWallCorner skips it: every bore-wall face falls
        // back to the typed coefficient and the derived law never reaches the
        // surface this step exists to heat.
        //
        // ⚠⚠ AND THE TYPED RADIUS CANNOT SIMPLY BE LOWERED. It is read in THREE
        // coupled places -- the kernel's exhaust-vs-column choice (:2280), the
        // 7a/7b Newton-regime split (:2314), and this datum's inner truncation
        // (ColumnWallCorner) -- so moving it shifts z, and hence WallJetH, over
        // the WHOLE wall band, not just the new annulus. P2-1 recorded the datum
        // moving the answer 1.04-1.74x, which would confound exactly the
        // measurement P2-3 is making.
        //
        // `bore` takes the truncation from the per-step MEASURED bore radius
        // (wall_r_bore, already computed from the void mask and already logged)
        // less one cell of discretisation slack. The datum then follows the hole
        // as it widens, introduces no new tuning number, and is DECOUPLED from
        // the typed gate so the wall's effect stays attributable.
        // Default `typed` => byte-identical, and P2-1's accepted calibration is
        // untouched.
        bool side_face_h_datum_bore = false;
        // P2-3 deliverable 5: multiplies the DERIVED law's output at the one
        // resolver and nowhere else. 1.0 => byte-identical.
        Set::Scalar side_face_h_scale = 1.0;
        // P2-3 deliverable 2: march the gas UP THE BORE, depleting it against
        // the wall the way the floor march depletes it across the floor.
        // Default 0 => byte-identical.
        bool jet_wall_march = false;
        // P2-4: THE COLD BODY IN THE ANNULUS. P2-3 measured the rock wall as a
        // near-equilibrium surface -- it soaks to ~1118 K and the gas-to-wall
        // dT collapses, so the march sheds only 53 K where the mouth needs
        // ~410 K to stop firing. The rock cannot supply that at any
        // coefficient. The one surface in Meier's hole that CANNOT equilibrate
        // is the burner's water-cooled mantle (Ø 67 mm, 160 L/h), standing in
        // an 85-95 mm bore from the nozzle flange upward. Its radius is already
        // in this file as WJ_R_SHIELD, so the body adds NO new geometry number.
        // ⚠ It is a second one-sided debit per bin, not a face: it has no
        // area law, no criterion and no removal -- which is why SideFaceH's
        // role count and the cap-site count must stay where they are.
        bool jet_body = false;
        // the mantle's surface temperature [K]. ⚠ P6: 320 +/- 30 K moves
        // jet_P_body by +/-2.5 %, so this key carries no magnitude; the
        // immersion does, and that is MEASURED (body_z_imm), not assumed.
        Set::Scalar jet_body_T = 320.0;
        Set::Scalar jet_body_h_scale = 1.0;
        // which mass the WALL section marches. `jet` = the floor march's mass
        // including entrained; `fresh` = the nozzle's own mdot, on the reading
        // that the entrained pit gas stays below the plane and does not climb
        // the bore. Under `fresh` the retained share is booked so the budget
        // still closes -- see the closure below.
        enum class WallMdot { Jet, Fresh };
        WallMdot jet_wall_mdot = WallMdot::Jet;
        // D2x: heat a side face only where its column's radius from (x0, y0) is
        // >= this (0 = no restriction). The map's cone flank and the hole wall
        // need different treatment and one side term serves both.
        Set::Scalar side_face_r_min = 0.0;
        // D2z-2: radius gate for WallStepDiagnostic ONLY. 0.0 = whole domain =
        // D2z behaviour. Diagnostics-only: it selects which faces are averaged
        // and touches no flux, area or trajectory.
        Set::Scalar step_diag_radius = 0.0;
        // D2aa deliverable 2: an explicit heat-transfer pair for the BELOW-PLANE
        // lateral faces, replacing the top-face map they inherit today.
        // Sentinel -1 = inherit = exactly today's behaviour, so identity holds.
        // ⚠ BOTH-OR-NEITHER at parse: setting h without T_gas REFRIGERATES the
        // wall, because T defaults to the 293 K ambient. That has bitten this
        // line before and it is an abort, not a warning.
        Set::Scalar side_face_h = -1.0;
        Set::Scalar side_face_T = -1.0;
        // D2ab deliverable 1: surface_patch.side_face_T = live tracks the
        // march's OWN jet_T_exhaust instead of a frozen constant. D2aa fed the
        // wall 1463.5 K lifted from SW while WL and WL81's own exhaust ran at
        // 1601.3 and 1631.0 K -- gas that gives up less heat to the wall stays
        // hotter -- so the driving dT was 22 %/27 % too small and every D2aa
        // wall number was self-inconsistent. Default false = the constant form,
        // which is what reproduces D2aa.
        bool        side_face_T_live = false;
        // D2ab deliverable 2: the explicit (h, T) pair is the ANNULUS
        // coefficient, so it may only reach faces that are IN the annulus --
        // beyond the burner radius. Below r_min a face keeps the INHERITED
        // impinging map. Measured on D2aa's own region columns, 64.6 % of the
        // side power sat INSIDE foot_r_outer (bore staircase risers and pad-ring
        // step faces, in impinging flow, not in the annulus), so D2aa applied
        // the annulus coefficient to two thirds of the wrong faces.
        // 0.0 = today's behaviour (the pair reaches every lateral face), so an
        // unset run takes the identical branch and identity holds.
        Set::Scalar side_face_h_r_min = 0.0;
        // P2-1 deliverable 3: the wall coefficient's LAW.
        //   typed   (DEFAULT) -- side_face_h applied flat at every height, the
        //                        shipped behaviour; every frozen anchor from
        //                        d2v onward stays byte-identical.
        //   derived           -- the composite wall-jet/duct law, WallJetH().
        // ⚠ It REPLACES the pair's coefficient and is read only where that pair
        // is applied, so it is scoped by exactly the same radius gate.
        bool side_face_h_derived = false;
        // D2aa deliverable 1: open fraction of the pad ring [foot_r_inner,
        // foot_r_outer). 0 = fully shielded (no gas under the pads), 1 = today.
        // Sentinel -1 = inert. ⚠ A BOUND, not a fitted coefficient.
        Set::Scalar pad_shield_open = -1.0;
        // D2x: replace the side_face_factor fudge on BELOW-plane faces with the
        // derived per-edge riser area (see RiserWeight). f multiplies an AREA
        // inside h*A*(T_gas - T), where an area factor and an h factor are
        // algebraically indistinguishable, so no output can falsify it.
        bool side_face_area_exact = false;
        // D2z: the step-height closure. Pay each riser its TRUE sub-cell height
        // delta = |n + phi_S - phi_T| (BOTH columns' phases), and correct the
        // lateral face conductance across the step by
        //   g(delta) = (1 - exp(-delta*a)) / (1 - exp(-a)),  a = dz*v/alpha
        // with v = dz/dt_i from the column's own trailing flip interval and
        // alpha = kappa/(rho*cp) the thermal DIFFUSIVITY. g(1) = 1 exactly, so
        // delta = 1 is today's un-weighted face and the scheme degrades to
        // current behaviour rather than to zero. No tuning constant.
        bool side_face_step_height = false;
        // ⚠⚠ P1e: the staircase ORIENTATION correction, default OFF.
        // Multiplies the per-column side-face area sum by
        // hypot(hx, hy) / (hx + hy), the ratio of the true inclined area to the
        // staircase area, with hx/hy the summed riser heights of the open faces
        // on each axis. A DISCRETISATION FIX, not a threshold: exactly 1 on a
        // grid-aligned wall and on a thin pillar, 1/sqrt(2) on a 45 degree
        // corner, convergent with mesh, nothing fitted. See StaircaseFactor.
        // ⚠ NOT side_face_factor (a scalar user multiplier, default 1.0) -- a
        // different key with a different meaning. Do not conflate them.
        bool side_face_staircase = false;
        // Which law the HEAT sees. The conductance correction is identical on
        // both, so a leg pair isolates the heating law.
        //   step   : riser area = delta * dz * dy   (a real vertical face)
        //   wetted : the EXACT ramp (sqrt(1+delta^2) - 1) * dx * dy -- a phase
        //            lag on a smooth floor has no vertical face at all.
        // ⚠ NOT delta^2/2: that is the small-step limit and it OVERcharges
        // exactly where A's floor lives (2.0 vs a true 1.236 at delta = 2).
        enum class SideAreaLaw { Step, Wetted };
        SideAreaLaw side_face_area_law = SideAreaLaw::Step;
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
        Set::Scalar side_power_interval = 0.0;   // rev7 d2, 0 = off
        Set::Scalar jet_D = 0.0, jet_core_length = 5.0;
        bool jet_stag_decay = true, jet_entrained_mass = false;
        bool jet_T_ent_exhaust = false;     // D2b: jet_T_ent_mode = exhaust (or
                                            // D2e rec_fixed): recirculation on
        bool jet_rec_fixed = false;         // D2e: jet_T_ent_mode = rec_fixed
        Set::Scalar jet_T_rec_fixed = 0.0;  // D2e: its constant T_mix [K]
        // D2c: march options (see JetEnthalpyMarch); all default-off.
        bool jet_De_momentum = false;       // jet_decay_diameter = momentum
        Set::Scalar jet_De_ref = 0.0;       // full-jet D_e at jet_T_ent [m]
        bool jet_free_surface = false;
        Set::Scalar jet_fs_aspect = 1.0, jet_fs_exponent = 1.0;
        bool jet_bin_exponential = false;   // jet_bin_update = exponential
        bool jet_neg_count = false;         // jet_negative_flux = count
        // D2v: the confined h map ABOVE the nozzle plane. Under jet_closure =
        // enthalpy the jet branch hard-zeroes h for s <= 0 and leaves T_gas at
        // jet_T_ent. These two keys lift that, and BOTH are sentinels so a run
        // that sets neither is byte-identical.
        //   jet_above_h  = 1  -> evaluate the SAME h_expr for every s.
        //   jet_T_above >= 0  -> T_gas above the plane [K].
        // They are parsed as an inseparable pair: h without a hot T_gas
        // REFRIGERATES the wall (radiative + ambient losses are its only other
        // boundary conditions), which would read as a physics result.
        bool jet_above_h = false;
        Set::Scalar jet_T_above = -1.0;
        // D2b: burner-on-feet descent (nozzle_descent = feet).
        enum class FootRule { Pads, Max, Mean };
        bool feet = false;
        FootRule foot_rule = FootRule::Pads;
        Set::Scalar foot_r_inner = 0.0, foot_r_outer = 0.0, foot_standoff = 0.0;
        int foot_npads = 3;
        Set::Scalar foot_pad_quantile = 0.9;
        Set::Scalar nozzle_feed_max = std::numeric_limits<Set::Scalar>::infinity();
        bool foot_body_clearance = false;   // D2c

        // prescribed_T only: patch surface temperature at time t. Not valid
        // under convective_flame (T_f is left uncompiled in that mode).
        Set::Scalar T_at(Set::Scalar t) const { return T_f(t); }

        // D2aa deliverable 1. ⚠ THREE wrong homes before this became a helper,
        // each failing differently, which is why it is one:
        //   * beside pad_shield_open's range check it ran BEFORE the feet branch
        //     reads foot_r_inner/outer, so it aborted on exactly the legs that
        //     DID set the ring;
        //   * inside the feet branch it stopped firing at all, leaving the key
        //     SILENTLY INERT without feet -- worse than the abort;
        //   * at the end of Parse it was unreachable, because the expression
        //     path returns early at the `use_expr = true` exit and that is the
        //     path every study leg takes.
        // A guard must sit after everything it reads and on EVERY exit it has
        // to cover.
        static void CheckPadRing(const SurfacePatch &v)
        {
            if (v.pad_shield_open >= 0.0 && !(v.foot_r_outer > v.foot_r_inner))
                Util::Abort(INFO, "surface_patch.pad_shield_open needs a pad ring: "
                                  "foot_r_outer (", v.foot_r_outer, ") must exceed "
                                  "foot_r_inner (", v.foot_r_inner, "). The annulus "
                                  "is read under surface_patch.nozzle_descent = feet, "
                                  "so the shield needs a feet run.");
            // D2ab deliverable 1. Same both-exits home, for the same reason:
            // jet_closure is parsed BELOW the side_face block, so this cannot
            // live there, and the use_expr path returns early so it cannot live
            // at the tail of Parse either.
            //
            // ⚠ WITHOUT THE MARCH THERE IS NO jet_T_exhaust. `live` would then
            // sit on its first-step seed (jet_T_nozzle) for the whole run and
            // NOTHING in any output would say so -- the wall would simply be
            // fed a constant that nobody chose. That is the exact failure mode
            // this line of packets keeps paying for (the inert Meier lever, the
            // inert pad-ring guard), so it is an abort.
            if (v.side_face_T_live && !v.jet_enthalpy)
                Util::Abort(INFO, "surface_patch.side_face_T = live tracks the jet "
                                  "march's own jet_T_exhaust, which only exists under "
                                  "surface_patch.jet_closure = enthalpy. Without the "
                                  "march the key would SILENTLY hold its first-step "
                                  "seed for the entire run. Set jet_closure = enthalpy, "
                                  "or give side_face_T an explicit temperature in K.");
            // P2-1, and HERE for exactly the two reasons recorded above:
            // jet_closure parses BELOW the side_face block, and the use_expr
            // path returns early so the tail of Parse is unreachable.
            //
            // ⚠⚠ The derived law needs the burner's MOMENTUM and its exhaust
            // gas state. Both exist only under the enthalpy march; without it
            // the law would have no jet to derive from and would silently fall
            // back to the typed scalar it was meant to replace.
            if (v.side_face_h_derived && !v.jet_enthalpy)
                Util::Abort(INFO, "surface_patch.side_face_h_law = derived derives the "
                                  "wall coefficient from the burner's momentum and its "
                                  "exhaust gas state, both of which exist only under "
                                  "surface_patch.jet_closure = enthalpy. Without the "
                                  "march the law has no jet to derive from. Set "
                                  "jet_closure = enthalpy, or leave the law `typed`.");
            // P2-3, RE-HOMED HERE BY P2-4 -- it was unreachable in Parse because
            // BOTH values it reads are parsed later (see the note at the query).
            // Third member of the same family as the two above, same home, same
            // reason.
            if (v.jet_wall_march && !v.jet_enthalpy)
                Util::Abort(INFO, "surface_patch.jet_wall_march marches the ENTHALPY "
                    "jet's own stream up the bore and is read only inside "
                    "JetEnthalpyMarch, so without surface_patch.jet_closure = enthalpy "
                    "it would be SILENTLY INERT.");
            // P2-4: the body lives inside that march's bin loop.
            if (v.jet_body && !v.jet_enthalpy)
                Util::Abort(INFO, "surface_patch.jet_body is a sink inside the wall "
                    "march's bin loop, which runs only under "
                    "surface_patch.jet_closure = enthalpy.");
        }

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
            else if (robin_str == "surface") v.robin_form = RobinForm::Surface;
            else Util::Abort(INFO, "surface_patch.robin_form must be 'cell', 'face', "
                                   "'pinned' or 'surface', got: ", robin_str);
            // Revision 7 deliverable 6. Only meaningful under `surface`; read
            // unconditionally so ParmParse's strict unused-key check cannot
            // abort a leg that sets it while sweeping forms.
            {
                std::string side_str;
                pp.query_default("side_face_closure", side_str, std::string("follow"));
                if      (side_str == "follow") v.side_T_form = SideTForm::Follow;
                else if (side_str == "cell")   v.side_T_form = SideTForm::CellAvg;
                else Util::Abort(INFO, "surface_patch.side_face_closure must be "
                                       "'follow' or 'cell', got: ", side_str);
                if (v.side_T_form == SideTForm::CellAvg &&
                    v.robin_form != RobinForm::Surface)
                    Util::Message(INFO, "surface_patch.side_face_closure = cell is INERT "
                                        "unless robin_form = surface: it switches off "
                                        "the lateral surface closure, which only that "
                                        "form turns on.");
            }
            if (v.robin_form == RobinForm::Face && v.mode != Mode::ConvectiveFlame)
                Util::Abort(INFO, "surface_patch.robin_form = face requires "
                                  "surface_patch.mode = convective_flame");
            // C1: pinned-surface closure (see SurfaceCellFlux). A column that
            // has not fired for pinned_idle_cycles x (its own cell heating time)
            // falls back to the face form.
            if (v.robin_form == RobinForm::Pinned ||
                v.robin_form == RobinForm::Surface)
            {
                if (v.mode != Mode::ConvectiveFlame)
                    Util::Abort(INFO, "surface_patch.robin_form = ",
                                (v.robin_form == RobinForm::Pinned ? "pinned"
                                                                   : "surface"),
                                " requires surface_patch.mode = convective_flame");
                // ⚠ READ IT EVEN UNDER `surface`, WHERE IT IS INERT. ParmParse
                // strict-checks unused keys inside queryclass, so simply not
                // reading it aborts the run -- which is exactly what happened
                // to the first Meier rerun of this packet: every frozen leg
                // from d2v to d2ad sets pinned_idle_cycles = 1.0e9, so
                // switching one of them to robin_form = surface aborted at
                // parse before a single step ran. Deliverable 1.4 asks for the
                // key to be LOGGED as inert, not silently ignored; that is the
                // Message emitted with the other `surface` guards. Reading it
                // here is what makes the message reachable at all.
                pp.query_default("pinned_idle_cycles", v.pinned_idle_cycles, 2.0);
                if (!(std::isfinite(v.pinned_idle_cycles) && v.pinned_idle_cycles >= 0.0))
                    Util::Abort(INFO, "surface_patch.pinned_idle_cycles must be finite "
                                      "and >= 0, got: ", v.pinned_idle_cycles);
            }

            // D2k: side-face flux on exposed vertical faces of a staircase
            // wall. Default-off; with 0 the code path is byte-identical to
            // every pre-D2k result (no term, no new thermo column).
            {
                int sff = 0;
                pp.query_default("side_face_flux", sff, 0);
                if (sff != 0 && sff != 1)
                    Util::Abort(INFO, "surface_patch.side_face_flux must be 0 or 1, got: ", sff);
                v.side_face_flux = (sff == 1);
                if (v.side_face_flux)
                {
                    if (v.mode != Mode::ConvectiveFlame)
                        Util::Abort(INFO, "surface_patch.side_face_flux = 1 requires "
                                          "surface_patch.mode = convective_flame");
                    // D2x: query_default ADDS the key to the table (the D2c
                    // lesson), so "was it given by the user" must be captured
                    // BEFORE the first query_default that names it.
                    const bool f_given = pp.contains("side_face_factor");
                    pp.query_default("side_face_factor", v.side_face_factor, 1.0);
        // D2z-2: radius gate for the DRAIN DIAGNOSTIC only -- it selects which
        // faces WallStepDiagnostic averages and changes no flux, no area and no
        // trajectory. Default 0.0 = whole domain = D2z behaviour exactly, which
        // is what keeps the identity legs and the 2639-file tree gate valid.
        // ⚠ Deliberately NOT surface_patch.radius, which is 0.2 (the whole
        // domain) on the Meier input.
        pp.query_default("step_diag_radius", v.step_diag_radius, 0.0);
        // D2aa deliverable 2. ⚠ pp.contains() must be read BEFORE query_default,
        // which ADDS the key to the ParmParse table.
        {
            const bool h_given = pp.contains("side_face_h");
            const bool T_given = pp.contains("side_face_T");
            if (h_given != T_given)
                Util::Abort(INFO,
                    "surface_patch.side_face_h and surface_patch.side_face_T are "
                    "both-or-neither. Setting h without T_gas REFRIGERATES the wall: "
                    "side_face_T defaults to the 293 K ambient, so a lateral face "
                    "would be COOLED by a coefficient chosen to heat it. Give both, "
                    "or neither and inherit the top-face map.");
            pp.query_default("side_face_h", v.side_face_h, -1.0);
            // D2ab deliverable 1: side_face_T takes a temperature in K OR the
            // word `live`. Read as a STRING and converted here rather than
            // through the Scalar overload, because the Scalar path would reject
            // `live` with ALAMO's own parser error instead of this one.
            {
                std::string sfT;
                pp.query_default("side_face_T", sfT, std::string("-1"));
                if (sfT == "live") { v.side_face_T_live = true; v.side_face_T = -1.0; }
                else
                {
                    const char* b = sfT.c_str();
                    char* e = nullptr;
                    const double tv = std::strtod(b, &e);
                    while (e && (*e == ' ' || *e == '\t')) ++e;
                    if (!e || e == b || *e != '\0')
                        Util::Abort(INFO, "surface_patch.side_face_T must be a gas "
                            "temperature in K, or the word 'live' to track the "
                            "march's own jet_T_exhaust each step; got: ", sfT);
                    v.side_face_T = tv;
                }
            }
            // ⚠ contains() BEFORE query_default, which ADDS the key.
            const bool pad_given = pp.contains("pad_shield_open");
            pp.query_default("pad_shield_open", v.pad_shield_open, -1.0);
            if (pad_given
                && !(v.pad_shield_open >= 0.0 && v.pad_shield_open <= 1.0))
                Util::Abort(INFO, "surface_patch.pad_shield_open is the OPEN FRACTION "
                                  "of the pad ring and must lie in [0, 1]; got ",
                            v.pad_shield_open);
            if (h_given && !(v.side_face_h >= 0.0
                             && (v.side_face_T > 0.0 || v.side_face_T_live)))
                Util::Abort(INFO, "surface_patch.side_face_h must be >= 0 and "
                                  "surface_patch.side_face_T must be > 0 K (or "
                                  "'live'); got h = ",
                            v.side_face_h, ", T = ", v.side_face_T);
            // D2ab deliverable 2. ⚠ contains() BEFORE query_default, which ADDS
            // the key -- the trap that aborted every D2aa side_face_flux run.
            const bool rmin_given = pp.contains("side_face_h_r_min");
            pp.query_default("side_face_h_r_min", v.side_face_h_r_min, 0.0);
            // ⚠ Without the pair this key is SILENTLY INERT: it gates an
            // override that is not being applied. D2aa lost a leg to exactly
            // that failure mode (a guard in the wrong place stopped firing and
            // the key did nothing), so it is an abort, not a default.
            if (rmin_given && !h_given)
                Util::Abort(INFO, "surface_patch.side_face_h_r_min gates the explicit "
                    "(surface_patch.side_face_h, surface_patch.side_face_T) pair by "
                    "radius, so without that pair it is SILENTLY INERT -- it would "
                    "gate nothing and no output would show it. Set side_face_h and "
                    "side_face_T as well, or drop this key.");
            if (rmin_given && !(std::isfinite(v.side_face_h_r_min)
                                && v.side_face_h_r_min >= 0.0))
                Util::Abort(INFO, "surface_patch.side_face_h_r_min must be finite and "
                                  ">= 0 m, got: ", v.side_face_h_r_min);

            // P2-1 deliverable 3. Default `typed` == the shipped behaviour, so
            // an input that does not mention this key is byte-identical.
            std::string hlaw;
            pp.query_default("side_face_h_law", hlaw, std::string("typed"));
            if (hlaw == "derived") v.side_face_h_derived = true;
            else if (hlaw != "typed")
                Util::Abort(INFO, "surface_patch.side_face_h_law must be 'typed' or "
                                  "'derived', got: ", hlaw);
            // ⚠ Same reasoning as side_face_h_r_min above: the law replaces the
            // pair's coefficient at the pair's own call sites, so without the
            // pair it would be SILENTLY INERT -- the third inert key in this
            // line of packets. Abort, not a default.
            if (v.side_face_h_derived && !h_given)
                Util::Abort(INFO, "surface_patch.side_face_h_law = derived replaces the "
                    "coefficient of the explicit (surface_patch.side_face_h, "
                    "surface_patch.side_face_T) pair, and is read only where that pair "
                    "is applied, so without the pair it would be SILENTLY INERT. Set "
                    "side_face_h and side_face_T as well, or leave the law `typed`.");

            // P2-3 deliverable 4. Default `typed` == shipped.
            std::string hdat;
            pp.query_default("side_face_h_datum", hdat, std::string("typed"));
            if (hdat == "bore") v.side_face_h_datum_bore = true;
            else if (hdat != "typed")
                Util::Abort(INFO, "surface_patch.side_face_h_datum must be 'typed' or "
                                  "'bore', got: ", hdat);
            // ⚠⚠ P2-4 FIX, A DOUBLY-DEAD P2-3 ABORT. The query used to sit
            // BELOW this point, so `v.jet_wall_march` was still its `false`
            // default and the check never fired. Moving the query up exposed the
            // SECOND half of the same bug: `v.jet_enthalpy` is not parsed here
            // either -- `jet_closure` is read at :10511, far below the
            // side_face block, exactly as the D2ab note at :10021 records. So
            // the abort read TWO unparsed values and was dead twice over;
            // naming `jet_wall_march = 1` without `jet_closure = enthalpy` went
            // silently inert, the very failure it was written to prevent.
            // ⇒ the QUERY stays here (every P2-4 cross-key check below reads
            //   it, and those keys are all parsed in this block), and the
            //   jet_enthalpy-dependent abort MOVES to `CheckPadRing`, the
            //   post-parse both-exits home D2ab and P2-1 already use for
            //   precisely this reason.
            pp.query_default("jet_wall_march", v.jet_wall_march, false);
            if (v.side_face_h_datum_bore && !v.side_face_h_derived)
                Util::Abort(INFO, "surface_patch.side_face_h_datum = bore sets the "
                    "DERIVED wall law's height datum and is read nowhere else, so "
                    "without side_face_h_law = derived it would be SILENTLY INERT.");

            // P2-4 deliverable 1. Every key defaults to shipped behaviour and
            // aborts if named without the feature it extends -- the inert-key
            // discipline this line of packets has been bitten by four times.
            // ⚠⚠ `contains` MUST BE SAMPLED BEFORE `query_default`. ALAMO's
            // ParmParse REGISTERS a key when it is queried with a default, so
            // `pp.contains("jet_body_T")` AFTER the query is unconditionally
            // true and the guard below fired on every leg that never named the
            // key -- including the four P2-3 legs, which aborted rc 6. The
            // existing `h_given` flag for side_face_h is sampled the same way
            // for the same reason.
            const bool body_T_given = pp.contains("jet_body_T");
            const bool body_s_given = pp.contains("jet_body_h_scale");
            pp.query_default("jet_body", v.jet_body, false);
            pp.query_default("jet_body_T", v.jet_body_T, 320.0);
            pp.query_default("jet_body_h_scale", v.jet_body_h_scale, 1.0);
            {
                std::string wm;
                pp.query_default("jet_wall_mdot", wm, std::string("jet"));
                if (wm == "fresh") v.jet_wall_mdot = SurfacePatch::WallMdot::Fresh;
                else if (wm != "jet")
                    Util::Abort(INFO, "surface_patch.jet_wall_mdot must be 'jet' or "
                                      "'fresh', got: ", wm);
            }
            if (v.jet_body && !v.jet_wall_march)
                Util::Abort(INFO, "surface_patch.jet_body is a second sink inside the "
                    "WALL MARCH's bin loop and is read nowhere else, so without "
                    "surface_patch.jet_wall_march = 1 it would be SILENTLY INERT.");
            if (body_T_given && !v.jet_body)
                Util::Abort(INFO, "surface_patch.jet_body_T is the cold body's surface "
                    "temperature and is read only when surface_patch.jet_body = 1.");
            if (body_s_given && !v.jet_body)
                Util::Abort(INFO, "surface_patch.jet_body_h_scale scales the cold body's "
                    "coefficient and is read only when surface_patch.jet_body = 1.");
            if (!(std::isfinite(v.jet_body_h_scale) && v.jet_body_h_scale > 0.0))
                Util::Abort(INFO, "surface_patch.jet_body_h_scale must be finite and > 0 "
                    "(a silent 0 is how a sink becomes inert without saying so), got: ",
                    v.jet_body_h_scale);
            if (!(std::isfinite(v.jet_body_T) && v.jet_body_T > 0.0))
                Util::Abort(INFO, "surface_patch.jet_body_T must be finite and > 0 K, "
                    "got: ", v.jet_body_T);
            if (v.jet_wall_mdot == SurfacePatch::WallMdot::Fresh && !v.jet_wall_march)
                Util::Abort(INFO, "surface_patch.jet_wall_mdot selects which mass the "
                    "WALL MARCH carries and is read nowhere else, so without "
                    "surface_patch.jet_wall_march = 1 it would be SILENTLY INERT.");

            // P2-3 deliverable 5. ⚠ `jet_wall_march` is queried ABOVE now; see
            // the note there for why keeping it here made its abort dead.
            pp.query_default("side_face_h_scale", v.side_face_h_scale, 1.0);
            if (!(std::isfinite(v.side_face_h_scale) && v.side_face_h_scale > 0.0))
                Util::Abort(INFO, "surface_patch.side_face_h_scale must be finite and "
                                  "> 0, got: ", v.side_face_h_scale);
            if (v.side_face_h_scale != 1.0 && !v.side_face_h_derived)
                Util::Abort(INFO, "surface_patch.side_face_h_scale multiplies the "
                    "DERIVED law's output only, so without side_face_h_law = derived "
                    "it would be SILENTLY INERT.");
        }
                    if (!(std::isfinite(v.side_face_factor) && v.side_face_factor >= 0.0))
                        Util::Abort(INFO, "surface_patch.side_face_factor must be finite "
                                          "and >= 0, got: ", v.side_face_factor);
                    // D2q-2b item 7: cap on the coefficient a SIDE face may see.
                    // Applied only when spall.per_cell_removal = 1 (checked at
                    // the use site), so it is inert for every pre-D2q-2b run.
                    //
                    // The default 80 W/m2K sits BELOW h* -- the coefficient at
                    // which a wall reaches the 821 K firing point inside its
                    // 121 s dwell (116 at T_gas = 1400, 88.5 at the arbiter's
                    // 1600, 75 at 1750). So the wall provably cannot be heated
                    // into spalling on its own: any lateral removal has to come
                    // from the corner drain, which is the claim under test.
                    // DECLARED, never fitted -- a failure routes to the packet's
                    // decision rule, not to another value.
                    pp.query_default("side_face_factor_above",
                                     v.side_face_factor_above, -1.0);
                    if (v.side_face_factor_above >= 0.0 &&
                        !std::isfinite(v.side_face_factor_above))
                        Util::Abort(INFO, "surface_patch.side_face_factor_above must be finite "
                                          "and >= 0 (unset = use side_face_factor), got: ",
                                    v.side_face_factor_above);
                    pp.query_default("side_face_h_max", v.side_face_h_max, 80.0);
                    if (!(std::isfinite(v.side_face_h_max) && v.side_face_h_max > 0.0))
                        Util::Abort(INFO, "surface_patch.side_face_h_max must be finite "
                                          "and > 0, got: ", v.side_face_h_max);
                    // D2x: heat only faces at r >= side_face_r_min.
                    pp.query_default("side_face_r_min", v.side_face_r_min, 0.0);
                    if (!(std::isfinite(v.side_face_r_min) && v.side_face_r_min >= 0.0))
                        Util::Abort(INFO, "surface_patch.side_face_r_min must be finite "
                                          "and >= 0, got: ", v.side_face_r_min);
                    // D2x: the derived per-edge riser area.
                    int aex = 0;
                    pp.query_default("side_face_area_exact", aex, 0);
                    if (aex != 0 && aex != 1)
                        Util::Abort(INFO, "surface_patch.side_face_area_exact must be 0 "
                                          "or 1, got: ", aex);
                    v.side_face_area_exact = (aex == 1);
                    // Two area laws must never be set at once: f multiplies an
                    // AREA inside h*A*(T_gas - T), so a stray f would silently
                    // shadow the derived weight and nothing downstream could
                    // tell them apart.
                    if (v.side_face_area_exact && f_given)
                        Util::Abort(INFO, "surface_patch.side_face_area_exact = 1 replaces "
                                          "surface_patch.side_face_factor on below-plane "
                                          "faces; do not set both");
                    // D2z: the step-height closure. Same reasoning as above --
                    // it replaces the area factor outright, so a stray f or the
                    // D2x derived weight alongside it would shadow the derived
                    // delta and nothing downstream could tell them apart.
                    int sh = 0;
                    pp.query_default("side_face_step_height", sh, 0);
                    if (sh != 0 && sh != 1)
                        Util::Abort(INFO, "surface_patch.side_face_step_height must be 0 "
                                          "or 1, got: ", sh);
                    v.side_face_step_height = (sh == 1);
                    if (v.side_face_step_height && f_given)
                        Util::Abort(INFO, "surface_patch.side_face_step_height = 1 derives "
                                          "the riser area from the column phases and "
                                          "replaces surface_patch.side_face_factor; "
                                          "do not set both");
                    if (v.side_face_step_height && v.side_face_area_exact)
                        Util::Abort(INFO, "surface_patch.side_face_step_height = 1 and "
                                          "surface_patch.side_face_area_exact = 1 are two "
                                          "different area laws; set at most one");
                    // P1e. Default 0, so every frozen leg is byte-identical
                    // without it.
                    int sc = 0;
                    pp.query_default("side_face_staircase", sc, 0);
                    if (sc != 0 && sc != 1)
                        Util::Abort(INFO, "surface_patch.side_face_staircase must be "
                                          "0 or 1, got: ", sc);
                    v.side_face_staircase = (sc == 1);
                    // ⚠⚠ ALL FOUR AreaOfDelta SITES LIVE BEHIND
                    // side_face_step_height, so without it this key would be a
                    // SILENT NO-OP -- which is the single failure mode most
                    // likely to produce a clean-looking but meaningless leg.
                    // Refuse it instead.
                    if (v.side_face_staircase && !v.side_face_step_height)
                        Util::Abort(INFO, "surface_patch.side_face_staircase = 1 "
                            "requires surface_patch.side_face_step_height = 1: the "
                            "staircase correction multiplies the per-column area "
                            "sum at the four step-height AreaOfDelta sites, and "
                            "every one of them is inside that key's branch. "
                            "Without it the correction would silently do nothing.");
                    // The heating law. Parsed unconditionally so that naming it
                    // without the step key is refused rather than ignored --
                    // AMReX only WARNS about an unused ParmParse variable.
                    const bool law_given = pp.contains("side_face_area_law");
                    std::string law;
                    pp.query_default("side_face_area_law", law, std::string("step"));
                    if (law == "step")
                        v.side_face_area_law = SurfacePatch::SideAreaLaw::Step;
                    else if (law == "wetted")
                        v.side_face_area_law = SurfacePatch::SideAreaLaw::Wetted;
                    else
                        Util::Abort(INFO, "surface_patch.side_face_area_law must be 'step' "
                                          "or 'wetted', got: ", law);
                    if (law_given && !v.side_face_step_height)
                        Util::Abort(INFO, "surface_patch.side_face_area_law requires "
                                          "surface_patch.side_face_step_height = 1");
                }
                else if (pp.contains("side_face_h_max"))
                    Util::Abort(INFO, "surface_patch.side_face_h_max requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_factor_above"))
                    Util::Abort(INFO, "surface_patch.side_face_factor_above requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_factor"))
                    Util::Abort(INFO, "surface_patch.side_face_factor requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_r_min"))
                    Util::Abort(INFO, "surface_patch.side_face_r_min requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_area_exact"))
                    Util::Abort(INFO, "surface_patch.side_face_area_exact requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_step_height"))
                    Util::Abort(INFO, "surface_patch.side_face_step_height requires "
                                      "surface_patch.side_face_flux = 1");
                // ⚠ P1e: without this row the key would be SILENTLY IGNORED when
                // side_face_flux is off -- AMReX only WARNS about an unused
                // ParmParse variable -- which is the same silent-no-op failure the
                // step-height guard exists to prevent.
                else if (pp.contains("side_face_staircase"))
                    Util::Abort(INFO, "surface_patch.side_face_staircase requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_area_law"))
                    Util::Abort(INFO, "surface_patch.side_face_area_law requires "
                                      "surface_patch.side_face_flux = 1");
                // D2ab: the pair and its radius gate act only on side faces.
                else if (pp.contains("side_face_h"))
                    Util::Abort(INFO, "surface_patch.side_face_h requires "
                                      "surface_patch.side_face_flux = 1");
                else if (pp.contains("side_face_h_r_min"))
                    Util::Abort(INFO, "surface_patch.side_face_h_r_min requires "
                                      "surface_patch.side_face_flux = 1");
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
                // D2e: rec_fixed = the exhaust path with a constant T_mix =
                // jet_T_rec_fixed (a diagnostic that opens the recirculation
                // loop): same budget, P_recirc, thermo, T_ent floor and
                // dilution ambient; jet_T_rec_next is not used.
                if      (tem == "fixed")     v.jet_T_ent_exhaust = false;
                else if (tem == "exhaust")   v.jet_T_ent_exhaust = true;
                else if (tem == "rec_fixed") v.jet_T_ent_exhaust = v.jet_rec_fixed = true;
                else Util::Abort(INFO, "surface_patch.jet_T_ent_mode must be 'fixed', "
                                       "'exhaust' or 'rec_fixed', got: ", tem);
                if (v.jet_T_ent_exhaust && !v.jet_entrained_mass)
                    Util::Abort(INFO, "surface_patch.jet_T_ent_mode = ", tem, " requires "
                                      "surface_patch.jet_entrained_mass = 1");
                if (v.jet_rec_fixed)
                {
                    // (query_default adds the key, so check "given" first)
                    if (!pp.contains("jet_T_rec_fixed"))
                        Util::Abort(INFO, "surface_patch.jet_T_ent_mode = rec_fixed requires "
                                          "surface_patch.jet_T_rec_fixed [K]");
                    pp.query_default("jet_T_rec_fixed", v.jet_T_rec_fixed, 0.0);
                    if (!std::isfinite(v.jet_T_rec_fixed) || v.jet_T_rec_fixed < v.jet_T_ent)
                        Util::Abort(INFO, "surface_patch.jet_T_rec_fixed must be finite and "
                                          ">= jet_T_ent (", v.jet_T_ent, " K), got: ",
                                    v.jet_T_rec_fixed);
                }
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
                // D2v: the above-plane pair. Each aborts without the other,
                // and jet_T_above must exceed jet_T_ent, because h against a
                // cold gas is a heat SINK on the wall, not a weaker source.
                int jah = 0;
                pp.query_default("jet_above_h", jah, 0);
                if (jah != 0 && jah != 1)
                    Util::Abort(INFO, "surface_patch.jet_above_h must be 0 or 1, got: ", jah);
                v.jet_above_h = (jah == 1);
                // (query_default adds the key, so read "given" first)
                const bool ta_given = pp.contains("jet_T_above");
                pp.query_default("jet_T_above", v.jet_T_above, -1.0);
                if (v.jet_above_h && !ta_given)
                    Util::Abort(INFO, "surface_patch.jet_above_h = 1 requires "
                                      "surface_patch.jet_T_above [K]: h above the plane "
                                      "without a hot T_gas refrigerates the wall");
                if (ta_given && !v.jet_above_h)
                    Util::Abort(INFO, "surface_patch.jet_T_above requires "
                                      "surface_patch.jet_above_h = 1");
                if (v.jet_above_h && !(std::isfinite(v.jet_T_above) &&
                                       v.jet_T_above > v.jet_T_ent))
                    Util::Abort(INFO, "surface_patch.jet_T_above must be finite and > "
                                      "jet_T_ent (", v.jet_T_ent, " K), got: ", v.jet_T_above);
                pp.query_default("jet_dr",        v.jet_dr,    0.0); // 0 -> max(dx, dy)
                // Optional per-bin profile log <plot_file>_jet_profile.csv every
                // jet_profile_interval seconds of simulated time (0 = off).
                pp.query_default("jet_profile_interval", v.jet_profile_interval, 0.0);
                // rev7 d2: per-column side-power log, seconds of simulated
                // time between rows. 0 = off, so no existing leg gains a file.
                pp.query_default("side_power_interval", v.side_power_interval, 0.0);
                if (!(std::isfinite(v.side_power_interval) &&
                      v.side_power_interval >= 0.0))
                    Util::Abort(INFO, "surface_patch.side_power_interval must be "
                                      "finite and >= 0, got ", v.side_power_interval);
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
                                      "jet_negative_flux", "jet_above_h", "jet_T_above"})
                    if (pp.contains(k))
                        Util::Abort(INFO, "surface_patch.", k, " requires "
                                          "surface_patch.jet_closure = enthalpy");
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
                    // D2c: burner-body clearance (rock above a pad is removed
                    // mechanically at the next removal pass; see FeetDescent).
                    int cl = 0;
                    pp.query_default("foot_body_clearance", cl, 0);
                    if (cl != 0 && cl != 1)
                        Util::Abort(INFO, "surface_patch.foot_body_clearance must be 0 or 1, got: ", cl);
                    v.foot_body_clearance = (cl == 1);
                    if (v.foot_body_clearance && v.foot_rule != FootRule::Pads)
                        Util::Abort(INFO, "surface_patch.foot_body_clearance = 1 requires "
                                          "surface_patch.foot_rule = pads");
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
                // D2aa ⚠ THE EARLY RETURN. Every expression-based run -- which
                // is every study leg and the unit test -- leaves Parse here, so
                // a guard placed after it is not a guard at all (the D2z-0
                // lesson, repeated). Checked on BOTH exits.
                CheckPadRing(v);
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

            // D2aa: the non-expression exit takes the same check.
            CheckPadRing(v);
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
        // P2-2 change 2: the per-phase property LAW.
        //   constant (DEFAULT) -- the scalar kappa/beta/E/mu/rho/Cp above, the
        //                         shipped behaviour; every frozen study is
        //                         byte-identical.
        //   grimsel            -- the MEASURED Grimsel / Central Aare laws,
        //                         Numeric::Material::Grimsel, including the
        //                         tabulated thermal strain (change 5).
        // ⚠ The scalars stay parsed and are what `constant` uses, so turning
        // the law on is the only thing that changes behaviour.
        bool grimsel_law = false;
        Numeric::Material::Grimsel grimsel;
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
    // heterogeneous thermoelastic update. Cell-centred, 1-ghost; sampled by
    // CellToNodeAverage in UpdateModel.
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

    // D2z-1: whole-cell removal is UNCONDITIONAL -- there is no key and no
    // half-cell mode. Under computed-T the depth scan's surface region is flat
    // across the whole top half-cell (the T_face = T_top fallback), so at a
    // marginal firing the K_I crossing landed just past dz/2 and a column fired
    // TWICE per cell. Now a firing voids the top cell (h = dz), so an event in
    // removal_events.csv IS a removal.
    //   * The depth scan is RETAINED under prescribed_T. That is CONSERVATISM,
    //     not a justified split: D2z-0 takeaway 7 found the same dz/2 landing
    //     in a prescribed_T control with a real gradient, which points at the
    //     GATE LOCATION (Sp_top >= 1 tested at the cell CENTRE while depth is
    //     measured from the FACE) rather than the lerp. Keeping it untouched
    //     keeps Hu and Kant unchanged; fixing it is a Step 16d question.
    //   * Applies to the SPALL regime only. Vapor sets its own h_v and the bit
    //     is drawworks-set; neither had the half-cell artefact.
    // Registered ONLY when spall.surface_normal = 1, so it appears in exactly
    // the four inputs it describes and changes no other thermo.dat header.
    // Mean 1/cos(theta) over each step's spall-firing columns = the factor by
    // which recession WOULD have been scaled before whole-cell discarded the
    // scan's result. 1.0 means nothing was lost.
    Set::Scalar sn_recession_superseded = 0.0;

    // D2q-2b item 5: per-cell removal (default off). The two preconditions
    // are cross-checked in Parse() once every block has been read, since they
    // live in other parse blocks.
    bool spall_per_cell_removal = false;
    // Item 8: thermo counters, registered only when the key is on.
    // pcr_overhang_cells = void cells with solid directly above (the undercut
    // census -- a one-cell fillet is expected, a deep undercut is a FINDING,
    // not a bug, because nothing in the model collapses one).
    // pcr_lateral_cells  = cells removed this step through a lateral face.
    Set::Scalar pcr_overhang_cells = 0.0;
    Set::Scalar pcr_lateral_cells  = 0.0;

    // D2q-2b item 0: debug-only seeded interior void, `ilo jlo klo ihi jhi khi`
    // in index space. Empty = off (the default), and ApplySeedVoid returns
    // immediately. See the parse block for why it carries no diagnostic.
    std::vector<int> spall_seed_void;

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
    // ⚠⚠ P2-2b: THE MECHANICAL HALF OF THE GRIMSEL SET, KEYED SEPARATELY FROM
    // THE THERMAL HALF. `microstructure.phaseN.law = grimsel` routes lambda,
    // rho and c_p into conduction; this routes E and the thermal strain into
    // the CRITERION. They must be independently switchable or the Kant ladder
    // collapses -- `K1` (thermal only) and `K2` (+ mechanical) would be the
    // same run, and `(c1)` could not say which property moved the onset.
    // Default false ⇒ the criterion keeps its pre-P2-2b scalar arithmetic bit
    // for bit even when the thermal law is on.
    bool sp_stress_law_grimsel = false;
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
    // P0: per-cell witnesses for spall.surface_T = reconstructed, written at
    // each column's scored surface cell in UpdateSpAfterMechanics. Registered
    // only when the key is on, so the 2639-file identity tree -- which never
    // sets it -- gains no plotfile variable.
    Set::Field<Set::Scalar> T_surf_field_mf;
    Set::Field<Set::Scalar> ell_eff_field_mf;
    // P0 thermo.dat witnesses (registered only under the key, so no pre-P0
    // thermo.dat gains or loses a column). ⚠ Registered `extensive = false`,
    // which the Integrator does NOT reduce across ranks -- so these are summed
    // by hand onto every rank before assignment (D2z-0 lesson 2).
    Set::Scalar recon_cols        = 0.0;  // columns reconstructed this step
    Set::Scalar recon_flat        = 0.0;  // ... with no usable gradient
    Set::Scalar recon_ell_clamped = 0.0;  // ... where ell_eff >= dz/4 bound
    Set::Scalar recon_T_clamped   = 0.0;  // ... where the gas ceiling bound
    Set::Scalar recon_fallback    = 0.0;  // ... skipped: no solid second cell
    Set::Scalar recon_hf_capped   = 0.0;  // ... where h_f hit h_spall_max
    Set::Scalar recon_hf_bounded  = 0.0;  // ... where h_f hit spall.flake_h_max
    Set::Scalar recon_dT_mean     = 0.0;  // mean (T_surf - T_top) [K]
    // P0b deliverable 3.
    Set::Scalar sp_flake_c_f = 0.4;
    // Revision 7 deliverable 1: the EMPIRICAL flake-thickness bound, in metres.
    // 0 = inert. ⚠ A STOPGAP for a gradient-blind criterion; see ReconFromX.
    Set::Scalar sp_flake_h_max = 0.0;
    // P1. Default `absolute` + `flake` keeps every existing leg byte-identical,
    // P1b: unsupported hot rock comes off. Default OFF, so every frozen leg is
    // byte-identical without it. The rule keys on CURRENT roughness, which is
    // what makes it self-extinguishing: a smooth wall reads 1 exposed side and
    // a smooth floor reads 0, so as the floor levels the rule runs out of work.
    // ⚠ A PHYSICAL LENGTH [m], never a cell count -- a cell count would make a
    // mesh quantity set a physical one (the `ell >= dz/4` mistake). At 1 cell
    // the ordinary staircase phase lag (5-9 %, a converged-floor feature) puts
    // one-cell steps everywhere and the rule would fire a third of a HEALTHY
    // floor; measured offline, 6 mm fires 0.0 % of rev 7's floor and 35-41 % of
    // P1's stalled one. See tests/.../p1b_support_rule/exposure.py.
    // ⚠⚠ P1f: the LATERAL mode. Both default OFF.
    // The vertical criterion reads a per-COLUMN 1-D profile, and after P1e's area
    // correction a diagonal wall column still carries 1.4142 face-equivalents of
    // an axis column -- the same arc with sqrt(2) fewer columns. So no per-column
    // criterion can be azimuthally isotropic on a round wall, however the area is
    // corrected. These keys put the criterion on the FACE instead.
    bool sp_lateral_drive   = false;  // spall.lateral_drive
    bool sp_lateral_removal = false;  // spall.lateral_removal (regime 6)
    bool sp_lateral_subcell = false;  // spall.lateral_subcell (P1g-1 banking)
    bool sp_lateral_thermal = false;  // spall.lateral_thermal (P1g-2 the sink)
    // ⚠⚠ P1g-3: the flake thickness, derived instead of picked.
    bool sp_flake_len_flux     = false;  // spall.flake_length  = flux
    bool sp_flake_frac_weibull = false;  // spall.flake_fraction = weibull
    // ⚠⚠ W_CONE'S PROVENANCE, since it is the one new constant. It is the MEDIAN
    // of the cone-weighted failure density
    //      F(w) = 1 - e^{-w} (1 + w + w^2/2)
    // which is Rauenzahn Eq. 4.96-4.100 / Walsh 2014 Eq. 3: the failure
    // probability integrated over the volume of a cone whose cross-section grows
    // with depth through the aspect ratio, dV ~ x^2 dx. F(w) = 1/2 gives
    // w = 2.67 bare and w = 2.35 under Rauenzahn's normalisation.
    // ⚠ A MEDIAN, not a percentile anyone chose, and NOT fitted -- choosing w to
    // make a comparison land is forbidden by this step's guardrails. The 2.35
    // normalised value is used because the modulus it pairs with (`weibull.m`)
    // is Rauenzahn/Wilkinson's. See validation/spall_thickness/REFERENCES.md.
    static constexpr Set::Scalar W_CONE_MEDIAN = 2.35;
    // deliverable 7's two counters: LEN fell back to `twocell` because the face
    // carried no flame coefficient at all, or because the net flux was <= 0.
    // ⚠ The fallback is to `twocell`, NEVER to a clamp -- a clamp would hide the
    // absence behind a plausible number.
    Set::Scalar sp_len_fb_h0 = 0.0;
    Set::Scalar sp_len_fb_q0 = 0.0;
    // the LATERAL site's own pair -- reported separately, because the wall's
    // flux comes from a different source than the floor's and the two sites
    // fall back for different reasons.
    Set::Scalar lat_len_fb_h0_n = 0.0;
    Set::Scalar lat_len_fb_q0_n = 0.0;
    // P1g-2 diagnostics, registered only under spall.lateral_thermal.
    Set::Scalar lat_shift_n       = 0.0;  // faces whose row was shifted
    Set::Scalar lat_shift_refuse  = 0.0;  // refused: row shorter than the minimum
    Set::Scalar lat_shift_E       = 0.0;  // enthalpy the shift took out [J]
    Set::Scalar lat_shift_E_est   = 0.0;  // the INDEPENDENT analytic estimate [J]
    Set::Scalar lat_shift_ratio   = 0.0;  // displaced row volume / banked volume
    Set::Scalar lat_overlap_n     = 0.0;  // cells written by BOTH writers in a step
    Set::Scalar lat_orphan_n      = 0.0;  // banks left referring to no solid cell
    // P1g-2c deliverable 1: THE TRUTH DIAGNOSTIC. lat_shift_E is what the shift
    // INTENDED to move; lat_shift_real is what the FIELD actually lost, summed
    // over the touched cells before and after the write-back. ⚠⚠ The two
    // existing ratios (lat_shift_ratio, lat_shift_E_est/lat_shift_E) are both
    // built from the intended shift, so both read ~1 while the field loses 15 %
    // -- which is exactly why nothing caught the overwrite.
    Set::Scalar lat_shift_real    = 0.0;  // enthalpy the FIELD lost [J]
    Set::Scalar lat_shift_realfrac= 0.0;  // lat_shift_real / lat_shift_E
    Set::Scalar lat_touch_multi   = 0.0;  // cells written by MORE THAN ONE row
    // P1g-2c deliverable 2: without this a clean ratio cannot be told from a run
    // that simply had no corners -- which is how the old 45 deg test passed.
    Set::Scalar lat_twoaxis_n     = 0.0;  // shifting faces with BOTH axes active
    // P1g-2c deliverable 6: lat_shift_refuse conflated three outcomes under one
    // name that reads as a failure. Split: a genuinely bad decay length, the
    // VOID-at-gather case (the ~25,000, where the bank voided the face cell so
    // the row gathers length 0), and a real but too-short row.
    Set::Scalar lat_ref_ell       = 0.0;  // ell <= 0
    Set::Scalar lat_ref_void      = 0.0;  // runlen == 0: face cell already void
    Set::Scalar lat_ref_short     = 0.0;  // runlen in 1..ROWMIN-1
    // ⚠⚠ P1g-1: THE LATERAL RECESSION ACCUMULATOR, as PERSISTENT PER-COLUMN
    // state -- the same idiom as lat_score / lat_kfire / k_top_global, which are
    // per-column Vectors reduced across ranks rather than MultiFabs.
    //
    // Why not a MultiFab. The accumulator must CARRY ITS REMAINDER PAST A VOID:
    // once a cell voids, the lateral criterion skips it forever (it tests
    // rm > 0.5 first), so a per-cell accumulator can never void more than once
    // and the banking law would read 0.5 cells per firing instead of 0.75. The
    // remainder therefore has to live with the recession TRACK, not the cell.
    // Keying it to the column also avoids the face normal entirely, which is the
    // blocker that rules out mirroring the vertical remap: a two-open-face corner
    // cell has two normals and no principled split between them.
    // ⚠ The honest cost of keying to the column: recession at different k in one
    // column shares one accumulator. On the riser geometry the column recedes as
    // a single track so this is exact; on a re-entrant wall it would conflate two.
    // ⚠ NOT a MultiFab also keeps every plotfile's component list unchanged, so
    // (a1)'s 2639 files cannot move for a reporting reason.
    amrex::Vector<Set::Scalar> lat_accum_col;   // banked depth per column [m]
    amrex::Vector<Set::Scalar> lat_nfire_col;   // that column's lateral firings
    // P1g-1 diagnostics, registered only under spall.lateral_subcell so a
    // keys-off run keeps its exact thermo.dat column set and owes no re-pin.
    Set::Scalar lat_bank_n      = 0.0;  // lateral firings that BANKED this step
    Set::Scalar lat_bank_hf     = 0.0;  // sum of h_f_lat banked this step [m]
    Set::Scalar lat_accum_sum   = 0.0;  // sum of the accumulator over solid cells [m]
    Set::Scalar lat_void_colmax = 0.0;  // max regime-6 voids in ONE column this step
    // ⚠ lat_void_colmax is the cascade witness named in the packet's inertness
    // list: voiding a cell and pushing the un-flaked remainder's enthalpy into
    // the neighbour could fire it next step and run away. It is a MAX, not a
    // mean -- a mean would average a cascade away.
    // P1f diagnostics, registered only under spall.lateral_drive.
    Set::Scalar lat_cand_n     = 0.0;  // (cell, face) pairs that reached K_I/K_Ic
    Set::Scalar lat_fires      = 0.0;  // columns whose lateral criterion fired
    Set::Scalar lat_voided_n   = 0.0;  // cells actually voided as regime 6
    Set::Scalar lat_no_faceT_n = 0.0;  // faces skipped for want of a face T
    Set::Scalar lat_ell_med    = 0.0;  // mean ell_lat over firing columns [m]
    // ⚠⚠ P1f deliverable 12: THE SUB-THRESHOLD SCORE. lat_q / lat_score are
    // written only inside `if (score >= 1.0 ...)` (Removal.H:675), so a face
    // that did NOT fire has no logged score and lat_q reads 0.0 by
    // initialisation -- which is exactly where a bound needs checking. These
    // two are taken BEFORE that branch, over every candidate, so (w)'s upper
    // bound and the (v) ladder have a population to read. lat_Tface_max is the
    // face temperature OF THE CANDIDATE ATTAINING lat_qmax, kept paired across
    // ranks (a separate max on each would decouple them).
    Set::Scalar lat_qmax       = 0.0;  // max K_I/K_Ic over ALL candidates [-]
    Set::Scalar lat_Tface_max  = 0.0;  // T_face of the candidate attaining it [K]
    // P1b diagnostics, registered only under spall.support_rule.
    // P0b witnesses, registered only under robin_form = surface.
    // surf_side_share is deliverable 1.2's requirement: side faces keep the
    // cell average, so the Meier rerun is only interpretable if the share of
    // P_robin they carry is on the record.
    Set::Scalar surf_cols_robin = 0.0;
    Set::Scalar surf_q_mean     = 0.0;
    Set::Scalar surf_side_share = 0.0;
    // P0b rev5 deliverable 9 bin 2: the side Robin power split by whether the
    // receiving cell is its own column's top cell (scored next step) or buried
    // beneath it. Registered with surf_side_share, so no run that does not
    // reconstruct gains a thermo column.
    Set::Scalar surf_side_P_top    = 0.0;
    Set::Scalar surf_side_P_buried = 0.0;
    // Revision 7 deliverable 2: the same power split at side_face_h_r_min --
    // the LIVE gate -- rather than at foot_r_inner/outer, which are 0 under
    // `prescribed` descent and sent everything into the `wall` bin.
    Set::Scalar surf_side_P_rin    = 0.0;   // r <  side_face_h_r_min
    Set::Scalar surf_side_P_rout   = 0.0;   // r >= side_face_h_r_min
    bool side_power_header_written = false;
    long side_power_last = -1;
    // P0b (c''): cumulative enthalpy carried out of the domain by flakes, the
    // sink the remap creates. ⚠ The existing ledger_E_H sums H over ALL cells
    // INCLUDING void ones (whose H is deliberately frozen rather than reset),
    // so removal has never been a sink in its eyes. This is what closes it.
    Set::Scalar ledger_E_flake = 0.0;
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

    // P0b deliverable 1.1: ONE surface temperature per column per step, built
    // once at the top of Advance from the same T_old the conduction kernel
    // will use, and read by BOTH criterion sites AND the Robin BC. That is the
    // whole point of the packet -- P0 left the criterion on a reconstructed
    // surface while the BC stayed on the cell average, and the two errors that
    // used to cancel stopped cancelling.
    std::vector<Set::Scalar> surf_T_col;     // T_surf [K] (0 = not built)
    std::vector<Set::Scalar> surf_ell_col;   // ell_eff in use [m]
    // ⚠⚠ P1g-3: WITHOUT THIS ARRAY THE ell_q COLUMN LOGS A CONSTANT ZERO.
    // SharedSurfaceRecon rebuilds a SurfaceRecon from these per-column arrays,
    // so a field that has no array here silently keeps its default -- the
    // physics used the flux length (through surf_hf_col) while the log reported
    // 0 on every one of 13,064 events. A derived quantity must be carried in
    // the same per-column state as the quantity it derives.
    std::vector<Set::Scalar> surf_ellq_col;  // the flux-based length used [m]
    std::vector<Set::Scalar> surf_ell_raw;   // this step's raw estimate
    std::vector<Set::Scalar> surf_hf_col;    // c_f * ell_eff [m]
    std::vector<Set::Scalar> surf_tbar_col;  // flake average [K]
    std::vector<Set::Scalar> surf_flag_col;  // 1 flat | 2 ell clamp | 4 T cap
    // The ONLY persistent per-column estimator state left after rev5: the last
    // x = dz/ell_eff that had content, for the flat-hold (1.3). ncols is the
    // DOMAIN column count, which no regrid changes, so this needs sizing but
    // no index repair. ⚠ rev5 deleted surf_delta (declared and zeroed, never
    // read -- deliverable 11) and surf_fired_prev with the skip-after-firing.
    std::vector<Set::Scalar> surf_x_hold;    // last valid x = dz/ell_eff

    // Previous step's ell_eff per column, for the step-to-step variation the
    // packet asks to be measured rather than assumed small. Diagnostic only;
    // sized lazily in Removal.H (ncols is the domain column count, which no
    // regrid changes).
    amrex::Vector<Set::Scalar> rc_ell_prev;

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

    // ------------------------------------------------------------------
    // P0: a real surface temperature for the firing criterion.
    //
    // ⚠ THE DEFECT THIS REPLACES. Under `spall.stress_source =
    // local_thermoelastic` (the Meier and Kant configuration) BOTH criterion
    // sites read the TOP CELL AVERAGE as the surface temperature:
    //   * UpdateSpAfterMechanics -- sigma_xx is flat at T_a(i,j,k) over the
    //     whole cell (a0 = 4-20 um << dz, so every K_I sample lands in it) and
    //     K_Ic is evaluated at the same value;
    //   * MMWSpalling/Removal.H:846-858 -- the depth scan's surface lerp has
    //     T_face = T_top, i.e. the top HALF-CELL is isothermal at the average.
    // For a steady ablation front T(z) = T0 + (T_s - T0) e^{-z/l}, l = kappa/v,
    // the cell average and the true surface are related by
    //     T_s = T0 + (Tbar - T0) * (dz/l) / (1 - e^{-dz/l}),
    // so the criterion's "surface" is mesh- AND rate-dependent: at Tbar = 823 K
    // the true surface is 1215 K at dz = 2 mm and 1003 K at 1 mm. D2d's
    // measured +43 % ROP from 2 mm to 1 mm is that ratio.
    //
    // THE FIX (spall.surface_T = reconstructed). For an exponential profile
    // the top two cell AVERAGES have the exact ratio e^{dz/l}, so l follows
    // from the profile alone -- no recession rate, no flame coupling, no
    // one-step lag, no persistent state. Everything below is a function of
    //     x = dz / l_eff = log((T_top - T0) / (T_top_minus_1 - T0)).
    //
    // ⚠ T0 IS phases[0].T_ref. That is the temperature at which the
    // thermoelastic driving stress is zero, and in every input this programme
    // uses it equals the undisturbed rock temperature (293.15 K in both
    // validation/kant/* and validation/meier/*). It is NOT a new constant and
    // it is NOT tunable here. If a future input sets T_ref away from the
    // initial rock temperature, the ratio loses its meaning -- the
    // `nsolid_fb` / `flat` witnesses in the event log are what would show it.
    // ⚠ THREE modes, not two, and the third is what makes the second
    // FALSIFIABLE. A threshold criterion fires when ITS OWN temperature
    // variable reaches the threshold, so "the reconstruction fires at 462 K on
    // every mesh" is nearly a tautology on its own -- the cell-average
    // criterion fires at 462 K of CELL AVERAGE on every mesh just as flatly.
    // What actually differs is the PHYSICAL STATE at firing: the cell-average
    // criterion fires the rock at a true surface temperature that moves with
    // dz. `diagnostic` computes and logs the reconstruction while leaving the
    // criterion on the cell average, so that number can be MEASURED on the
    // control leg instead of inferred.
    enum class SurfaceTMode { CellAverage, Diagnostic, Reconstructed };
    SurfaceTMode sp_surface_T_mode = SurfaceTMode::CellAverage;
    // Compute + log. False => every path below is skipped and the run is
    // byte-identical to pre-P0.
    bool sp_surface_T_recon = false;
    // Feed it to the criterion. False under `diagnostic`, so the criterion
    // stays bit-for-bit the pre-P0 criterion while the witness is written.
    bool sp_surface_T_apply = false;
    // P0b: `spall.removal_energy = flake` recedes phi by the PHYSICAL flake
    // thickness h_f = min(l_eff, dz) instead of forcing h_c = dz, and voids
    // the cell only when phi crosses zero (the existing mechanism). Requires
    // surface_T = reconstructed (l_eff has no other source).
    bool sp_removal_energy_flake = false;


    // The ONE place the reconstruction is written. Both criterion sites call
    // this, so they cannot drift apart. Pure function of its arguments --
    // no state, no fields, no MPI -- which is also what makes it unit-testable
    // and decomposition-independent.
    //
    //   x       = log((T_top - T0)/(T_1 - T0))          = dz / ell_eff
    //   T_surf  = T0 + (T_top - T0) * x/(1 - e^{-x})    -> T_top as x -> 0
    //   h_f     = min(ell_eff, dz) = dz / max(x, 1)
    //   T_bar   = T0 + (T_surf - T0) * (1 - e^{-u})/u,  u = h_f/ell_eff
    //                                                    = min(x, 1)
    //
    // Guards, in order: a non-decaying or non-finite pair degrades to the cell
    // average (flat = true, factor exactly 1.0 -- no correction exactly where
    // none is available); x is capped at 4 (ell_eff >= dz/4) so a near-tie in
    // the two averages cannot produce an unbounded surface; and T_surf is
    // capped at the configuration's gas temperature when one exists.
    // The raw two-cell estimate x = dz/ell_eff. Returns 0 when the pair carries
    // no usable decaying gradient ("flat"), which the CALLER resolves -- under
    // P0b it holds the column's last valid smoothed x (deliverable 1.3) rather
    // than dropping to the cell average, because dropping gives each such
    // column a x1.3 flux kick and 12.8 % of Meier firings are flat.
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

    // Everything from a given x. Split out from ReconstructSurfaceT so the
    // per-column builder can substitute a smoothed or held x (deliverables
    // 1.3, 5.1, 5.2) without duplicating a line of the algebra.
    // ⚠ P1g-3: `ell_flux` and `frac_over` are appended with defaults that
    // reproduce the shipped expression exactly, so every existing caller and
    // every run with both new keys absent is byte-identical. (b1) witnesses it.
    //   ell_flux <= 0  => use the two-cell ell_eff (the shipped length)
    //   frac_over <= 0 => use c_f            (the shipped fraction)
    static SurfaceRecon ReconFromX(Set::Scalar T_top, Set::Scalar T0,
                                   Set::Scalar dz, Set::Scalar x,
                                   Set::Scalar T_gas_cap, Set::Scalar c_f,
                                   Set::Scalar h_f_max,
                                   Set::Scalar h_f_bound = 0.0,
                                   Set::Scalar ell_flux = 0.0,
                                   Set::Scalar frac_over = 0.0)
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
        // is a MESH-INDEPENDENT multiple of the surface rise (P6). With P0's
        // min(ell_eff, dz) the multiple was mesh-dependent and collapsed to the
        // cell average wherever ell_eff exceeded dz.
        // ⚠⚠ P1g-3 deliverables 5 + 6, APPLIED HERE SO BOTH SITES SHARE IT.
        // length: the flux-based heated depth ell_q = k (T_s - T0) / q_net when
        // the caller could compute one, else the shipped two-cell ell_eff.
        // fraction: w/m from the cone-weighted median when FRAC is on, else c_f.
        // ⚠ `ell_q` is reported through r.ell_q so the event log can carry the
        // length that was actually used, not the one that was available.
        const Set::Scalar ell_used = (ell_flux > 0.0) ? ell_flux : r.ell_eff;
        const Set::Scalar frac_use = (frac_over > 0.0) ? frac_over : c_f;
        r.ell_q = (ell_flux > 0.0) ? ell_flux : 0.0;
        r.h_f   = (x > 0.0) ? frac_use * ell_used : dz;
        // ⚠⚠ ell_eff IS CLAMPED FROM BELOW (>= dz/4) BUT NOT FROM ABOVE, so a
        // nearly-flat column -- x small but not zero -- gives an UNBOUNDED
        // h_f. This is not hypothetical: the Meier rerun aborted at t = 129 s
        // with a column trying to void 287 cells (574 mm, more than the whole
        // block) off a two-cell ratio of 1.0014. The assertion caught it.
        // The cap is spall.h_spall_max, the key that already defines how deep
        // this model looks for a flake at all. ⚠ THAT IS NOT THE dz CLAMP THE
        // PACKET FORBIDS: dz is a MESH quantity and clamping there is what
        // made T_bar mesh-dependent in P0; h_spall_max is PHYSICAL and
        // mesh-independent, so the cap cannot reintroduce mesh dependence.
        // It binds on nearly-flat columns only -- h_f is 0.4 mm against a
        // 15 mm cap on every 1-D column leg -- and its bind fraction is
        // reported (flag bit 8).
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
        // in this model says a preheated slab should spall less. The real fix
        // is a gradient-aware criterion (buckling / energy release with a
        // flake thickness), which is a P1 DESIGN ITEM. This key buys the
        // Meier comparison; it does not earn it.
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

    static SurfaceRecon ReconstructSurfaceT(Set::Scalar T_top, Set::Scalar T_1,
                                            Set::Scalar T0, Set::Scalar dz,
                                            Set::Scalar T_gas_cap,
                                            Set::Scalar c_f,
                                            Set::Scalar h_f_max)
    {
        return ReconFromX(T_top, T0, dz,
                          RawSurfaceX(T_top, T_1, T0), T_gas_cap, c_f, h_f_max);
    }

    // Upper bound on the gas temperature this configuration can present to the
    // rock, used only as a ceiling on the reconstruction. Time-independent and
    // read from Parse'd values, so it introduces NO dependence on
    // BuildFlameColumns -- which runs AFTER the criterion (see the call-order
    // note in Removal.H). Returns 0 when no gas temperature is defined, which
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

        // P2-2e: THE CONSTRAINT DEPTH. ⚠⚠ Queried BEFORE the sp_weibull early
        // return on purpose: ALAMO does not abort on an unused top-level key
        // (standing rule 1), so a key read only on the sp_weibull branch would
        // be SILENTLY INERT on a damage_law deck. ⚠ Every input of the abort
        // below is already parsed -- the P2-4 lesson that an abort reading an
        // unparsed member is dead code: spall.enabled and spall.stress_source in
        // ParseSpallSettings (Parse :164), surface_patch at :139, and
        // spallation.model just above, all before this function (:206).
        pp.query_default("spallation.sp.constraint_depth", value.sp_cref_d, 0.0);
        if (!(std::isfinite(value.sp_cref_d) && value.sp_cref_d >= 0.0))
            Util::Abort(INFO, "spallation.sp.constraint_depth must be finite and "
                              ">= 0 m, got: ", value.sp_cref_d);
        if (value.sp_cref_d > 0.0)
        {
            // The SAME conjunction Removal.H forms as `sp_weibull_use_T_face`
            // (sp_weibull_spall && (prescribed_T || local_thermoelastic)).
            const bool presc_T_face = value.surface_patch.enabled &&
                value.surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
            const bool use_T_face = value.spall_enabled &&
                value.spallation_model == SpallationModel::SpWeibull &&
                (presc_T_face || value.sp_local_thermoelastic);
            if (!use_T_face)
                Util::Abort(INFO, "spallation.sp.constraint_depth > 0 references the "
                    "criterion's drive to the temperature d BEHIND the face, and only "
                    "the T-derived criterion path reads it (spallation.model = "
                    "sp_weibull, spall.enabled = 1, and spall.stress_source = "
                    "local_thermoelastic or surface_patch.mode = prescribed_T). "
                    "Without that path the key would be SILENTLY INERT. Got d = ",
                    value.sp_cref_d, " m");
        }

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
        // runs at :206 and `value.phases` is not resized until :293, so a
        // `phases.empty()` test at this point is ALWAYS true and aborts every
        // legal leg. The check sits with the other cross-key aborts after the
        // microstructure block instead.

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
        // ⚠⚠ P1g-3 deliverable 4: `weibull.m` IS A MATERIAL CONSTANT, NOT A
        // WEIBULL-MODE KEY, and it is parsed OUTSIDE the mode guard for that
        // reason. Inside it, a run with `weibull.enabled = 0` never executed the
        // parse and silently took the MEMBER default 15.0 (:9608) -- which is
        // exactly what `unit/lateral_drive/input:120` does, while Meier's
        // `input_drilling:124` sets 20.0. The flake fraction below reads `m`, so
        // the fixture and the validation leg would have derived thicknesses
        // 33 % apart for no stated reason.
        // ⚠ IDENTITY-SAFE BY CONSTRUCTION: the query_default's fallback (15.0)
        // equals the member default at :9608, so any run that does not set the
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
            // P2-2 change 3. Default `cell` == the shipped behaviour.
            {
                std::string svol;
                pp.query_default("weibull.size_volume", svol,
                                 std::string("cell"));
                if (svol == "flux") value.sp_size_flux = true;
                else if (svol != "cell")
                    Util::Abort(INFO, "weibull.size_volume must be 'cell' or "
                                      "'flux', got: ", svol);
                // ⚠⚠ P2-2d d8: the flux reference, separate from V0.
                pp.query_default("weibull.V_ref", value.weibull_V_ref, -1.0);
                if (value.sp_size_flux
                    && !(std::isfinite(value.weibull_V_ref)
                         && value.weibull_V_ref > 0.0))
                    Util::Abort(INFO,
                        "weibull.size_volume = flux requires weibull.V_ref > 0."
                        " It is NOT weibull.V0: V0 is the mesh-independence "
                        "anchor (V_cell/V0 == 1 on both decks), while the flux "
                        "mode's stressed volume is ell_q^3. Borrowing V0 would "
                        "import the 8x cross-deck difference in V0 into a "
                        "quantity that has nothing to do with the mesh. Got: ",
                        value.weibull_V_ref);
            }

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
