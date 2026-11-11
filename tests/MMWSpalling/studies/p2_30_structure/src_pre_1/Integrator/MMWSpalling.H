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
#include "Numeric/SpCriterion.H"
#include "Numeric/Material/Material.H"
#include "Numeric/Material/Zhang.H"
#include "Numeric/Material/Grimsel.H"
#include "Numeric/Material/Constant.H"
#include "Numeric/Material/Table.H"
#include "BC/Constant.H"



#include <algorithm>
#include <iomanip>
#include <cmath>
// The face list's ONE WRITE collects per-cell deltas in a small ordered map before
// touching the field (FaceRemoval.H, included INSIDE the class body, so this
// cannot live there).
#include <map>
#include <array>
#include <queue>       // P2-11b: the surface path's Dijkstra
#include <functional>  // std::greater
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

        // P2-29: `solver` (explicit) and `alpha` (0 0) are retired at their defaults (the table above): one heat
        // path, forward Euler; the thermal expansion is each phase's beta.

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

            // P2-27 stage 3 (census C1b): D, kappa_eff and k_eff -- the damage
            // field, kappa_phase under another name, and a uniform bulk average --
            // are no longer registered: every later plotfile component moves up
            // three places (a declared re-pin of both gates, component-audited).
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
            // P2-28: surface.follow_mask is retired at 1 (the table above): the mask column top is the
            // face list's `top` class, and removal requires the list, so the switch IS the list's.
            value.surface_follow_mask = value.surface_patch.face_list;
            value.ledger_enabled      = (le != 0);
            // P0b deliverable 1: robin_form = surface. Same prerequisites as
            // pinned -- it is the same closure with the surface temperature
            // computed rather than remembered -- PLUS the reconstruction
            // itself, without which the form has no surface to sit on.
            if (value.surface_patch.enabled &&
                value.surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
            {
                // P2-28: surf_cols_robin / surf_q_mean (the column reconstruction's witnesses) are deleted
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

            // D2z-1 Goal 4 (revised): NO firings_per_void column. Under
            // whole-cell it is 1.000 by construction, so it carries no
            // information, and it appended after the eight jet_* columns and
            // broke jet_enthalpy criterion (a). The invariant lives in the
            // n_voided == 1 assertion in the column path (deleted, P2-28) instead.
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

protected:
    void Initialize(int lev) override
    {
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
        if (list_state()) InitPhi(lev);
        if (list_state()) ApplySeedVoid(lev);
        InitSurfaceMask(lev);
        // InitializeMicrostructure must run before ConvertTtoH on the
        // microstructure path because the H(T) inversion is per-phase.
        if (microstructure_enabled) InitializeMicrostructure(lev, /*announce=*/true);
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
        if (microstructure_enabled) CheckTimeStepLimit();
        HeatConduction::TimeStepBegin(a_time, a_step);
        // Step 7c/7d repair for the mechanics AMR path. With mechanics off the
        // thermal/removal AMR path keeps its fields across steps (regrid
        // repair happens in Regrid), so skip the per-step re-initialisation.
        if (microstructure_enabled && max_level > 0 &&
            m_type != Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);
    }

    // ⚠ P2-27: THE EXPLICIT TIME-STEP LIMIT, PER LEVEL. AdvanceMicrostructure is
    // forward Euler (H = H_old + dt*(div(kappa grad T) + Q - q_loss/dz)), stable
    // for dt <= dx^2/(2 d alpha) in d dimensions; harmonic faces, voids, the GB
    // resistance and the latent boxcars only LOWER the effective alpha, so the
    // bound uses each phase's plain diffusivity. Checked once per level, at the
    // first TimeStepBegin the level exists, against THAT level's own dt (ALAMO
    // subcycles: dt[l] = dt[l-1]/nsubsteps[l]) and its own dx -- never the level-0
    // dt against a finer dx, which is 2x too strict per level. alpha_max is the
    // largest phase diffusivity at the COLDEST temperature the run starts from,
    // T_cold = min(phase T_ref, losses.T_amb when set): Grimsel's lambda/(rho c_p)
    // DECREASES with T, so the cold end binds (297.5 K: 1.6746e-6 m^2/s, 0.2430 s
    // at 1.5625 mm); the property table's 200 K floor would reject production.
    // Over the limit, a face-list deck (the production path) ABORTS and any other
    // deck WARNS: the user, 2026-11-09 -- seven unit tests run a uniform-temperature
    // field at dt = 1 s with stand-in properties, where the scheme cannot go
    // unstable (studies/p2_27_cleanup/Z0.txt).
    // P2-28: a level whose top carries a PRESCRIBED patch has a top cell with 2d - 1 neighbour faces plus the
    // half-cell Dirichlet face (conductance 2 kc / dz): its diagonal is (2d + 1) alpha / dx^2, so the bound there is
    // dx^2/((2d + 1) alpha) (Hu at 3.125 mm and 280.15 K: 0.769 s against 0.897 s).
    void CheckTimeStepLimit()
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (lev < static_cast<int>(dt_limit_checked.size()) && dt_limit_checked[lev]) continue;
            if (lev >= static_cast<int>(dt_limit_checked.size())) dt_limit_checked.resize(lev + 1, 0);
            dt_limit_checked[lev] = 1;
            Set::Scalar a_max = -1.0, T_arg = 0.0;
            int p_arg = 0;
            for (int p = 0; p < static_cast<int>(phases.size()); ++p)
            {
                const MineralPhase& ph = phases[p];
                const Set::Scalar T_cold = (T_amb > 0.0) ? std::min(ph.T_ref, T_amb) : ph.T_ref;
                // ⚠ P2-29: `kappa` means two things -- Grimsel::kappa is the DIFFUSIVITY, Zhang::kappa
                // (Material's) the CONDUCTIVITY [W/(m K)], so the zhang diffusivity is formed here
                const Set::Scalar a = ph.grimsel_law ? ph.grimsel.kappa(T_cold)
                                    : ph.zhang_law   ? ph.zhang.kappa(T_cold)
                                                       / (ph.zhang.rho(T_cold) * ph.zhang.cp(T_cold))
                                                     : ph.kappa / (ph.rho * ph.Cp);
                if (a > a_max) { a_max = a; p_arg = p; T_arg = T_cold; }
            }
            const Set::Scalar* DX = geom[lev].CellSize();
            Set::Scalar dx = DX[0];
            for (int d = 1; d < AMREX_SPACEDIM; ++d) dx = std::min(dx, DX[d]);
            const bool presc_top = surface_patch.enabled && SurfaceOwner(lev) &&
                                   surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
            const int n_den = 2 * AMREX_SPACEDIM + (presc_top ? 1 : 0);
            const Set::Scalar dt_lim = dx * dx / (Set::Scalar(n_den) * a_max);
            Util::Message(INFO, "explicit time-step limit, level ", lev, ": dt = ", dt[lev],
                          " s, dt_lim = dx^2/(", n_den, " alpha) = ", dt_lim,
                          " s (dx ", dx, " m; phase ", p_arg, ", alpha ", a_max, " m^2/s at T_cold ",
                          T_arg, " K)");
            if (!(dt[lev] > dt_lim)) continue;
            if (surface_patch.face_list)
                Util::Abort(INFO, "level ", lev, ": timestep dt = ", dt[lev], " s exceeds the explicit "
                            "limit dt_lim = dx^2/(", n_den, " alpha) = ", dt_lim, " s (dx ", dx,
                            " m; phase ", p_arg, ", alpha ", a_max, " m^2/s at T_cold ", T_arg,
                            " K): the forward-Euler conduction step would amplify, not damp, the shortest "
                            "wavelength. Lower `timestep` (P2-27)");
            Util::Warning(INFO, "level ", lev, ": timestep dt = ", dt[lev], " s exceeds the explicit limit "
                          "dt_lim = ", dt_lim, " s (phase ", p_arg, ", alpha ", a_max, " m^2/s at T_cold ",
                          T_arg, " K); harmless only while the temperature field stays uniform. Not on the "
                          "face list, so the run continues (P2-27)");
        }
    }

    void TimeStepComplete(Set::Scalar a_time, int a_step) override
    {
        HeatConduction::TimeStepComplete(a_time, a_step);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepComplete(a_time, a_step);
        if (microstructure_enabled && max_level > 0 &&
            m_type != Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::Regrid(lev, time);
        InitSurfaceMask(lev);
        if (microstructure_enabled)
        {
            InitializeMicrostructure(lev, /*announce=*/false);
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
        // on the cell-centred heat conduction (AdvanceMicrostructure), but the
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
        // P2-29: the homogeneous branch (the global alpha[] x eta eigenstrain) is deleted: every deck states
        // its rock as phases (microstructure.nphases is required), so the branch above is the only one.
    }

    // P0b: declared here, ahead of its first use as a parameter type in
    // the column reconstruction (deleted, P2-28). An in-class member function BODY is a
    // complete-class context, but its DECLARATION is not.
    struct SurfaceRecon
    {
        Set::Scalar T_surf  = 0.0;   // reconstructed surface temperature [K]
        Set::Scalar ell_eff = 0.0;   // recovered decay length [m]
        // the face list's old tangent depth k_c (T_s - T0) / q_net, logged
        // only (FaceRemoval.H); 0 everywhere else
        Set::Scalar ell_q   = 0.0;
        Set::Scalar h_f     = 0.0;   // flake thickness min(ell_eff, dz) [m]
        Set::Scalar T_bar   = 0.0;   // flake-average temperature [K]
        bool ell_clamped    = false; // the ell_eff >= dz/4 floor bound
        bool h_f_capped     = false; // h_f hit the spall.h_spall_max SEARCH depth
        bool h_f_bounded    = false; // h_f hit spall.flake_h_max (deliverable 1)
        bool T_clamped      = false; // the gas-temperature ceiling bound
        bool flat           = false; // no usable gradient; T_surf == T_top
    };

    // P2-21: the column optical depth (ColumnOpticalDepth), built once per step
    // per path when the loss table varies and the beam is on. Indexed
    // [col*nz + (k - klo)] as ColumnSolidAboveCount; tau_top is the optical depth
    // from the column's top solid face to the cell's top face, tau_bot to its
    // bottom face (tau_top + a_k dz, the SAME a_k the gather used). k_last is the
    // column's lowest solid cell (the transmitted power's exit), INT_MIN if none.
    // Declared here, ahead of its first use as a parameter type (as SurfaceRecon).
    struct BeamTau
    {
        amrex::Vector<Set::Scalar> top, bot;
        amrex::Vector<int>         k_last;
        // P2-24: the column's reflection R_col = R + (R_melt - R) f_L(T_top), one per column (ncols); used at
        // P_inc, P_trans and every cell's (1 - R) prefactor in that column
        amrex::Vector<Set::Scalar> Rcol;
        bool built = false;
    };
    // P2-21: a cell's (tau_top, tau_bot) for the kernels, captured by value in a
    // ParallelFor. Unbuilt => -1, -1 => BeamSource's closed form.
    struct TauView
    {
        const Set::Scalar* top = nullptr;
        const Set::Scalar* bot = nullptr;
        int ilo = 0, jlo = 0, klo = 0, nx = 0;
        std::size_t nz = 0;
        AMREX_GPU_HOST_DEVICE std::size_t at(int i, int j, int k) const
        {
            return (static_cast<std::size_t>(j - jlo) * nx + (i - ilo)) * nz + (k - klo);
        }
        AMREX_GPU_HOST_DEVICE Set::Scalar t(int i, int j, int k) const { return top ? top[at(i, j, k)] : -1.0; }
        AMREX_GPU_HOST_DEVICE Set::Scalar b(int i, int j, int k) const { return bot ? bot[at(i, j, k)] : -1.0; }
        // P2-24: the column's R when built, else the deck's scalar R (byte for byte)
        const Set::Scalar* rc = nullptr;
        AMREX_GPU_HOST_DEVICE Set::Scalar R(int i, int j, Set::Scalar R_scalar) const
        {
            return rc ? rc[static_cast<std::size_t>(j - jlo) * nx + (i - ilo)] : R_scalar;
        }
    };
    TauView ViewOf(const BeamTau& bt, int lev) const
    {
        TauView v;
        if (!bt.built) return v;
        const amrex::Box& d = geom[lev].Domain();
        v.top = bt.top.data();
        v.bot = bt.bot.data();
        v.rc  = bt.Rcol.empty() ? nullptr : bt.Rcol.data();
        v.ilo = d.smallEnd(0); v.jlo = d.smallEnd(1); v.klo = d.smallEnd(2);
        v.nx  = d.length(0);   v.nz  = static_cast<std::size_t>(d.length(2));
        return v;
    }

    // Override HeatConduction::Advance with the enthalpy method: ONE heat path, AdvanceMicrostructure
    // (P2-29: the const-alpha, material and implicit kernels are deleted; every deck states its rock as phases).
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        // P2-9 rev 2 item 7: the step's wall time (a thermo column only under
        // the face list; measured always, read by nothing else)
        const double fr_step_t0 = amrex::second();
        // P2-2 change 1: the property fields, refreshed at the cell
        // temperature BEFORE the criterion and the conduction dispatch read
        // them. No-op unless a phase sets law = grimsel or (P2-29) zhang.
        UpdatePropertyFields(lev);
        Mechanics<Model::Solid::Affine::Isotropic>::Advance(lev, time, dt);
        // Step 10/11, P2-28: THE removal step -- the face list's (FaceRemoval.H), the only one: the
        // criterion on every listed line and the recession, BEFORE the H/Temp swaps below, so what it
        // writes is the "old" state the thermal advance consumes.
        RemovalStep(lev, time, dt);
        std::swap(*H_mf[lev], *H_old_mf[lev]);
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        AdvanceMicrostructure(lev, time, dt);
        // Step 6b: every path evolves H_mf as the conserved state, so
        // dispatch a single InvertHtoT here. The microstructure variant
        // selects the per-phase Material::Table by cell phase ID.
        InvertHtoT(lev);
        {
            double ts = amrex::second() - fr_step_t0;
            amrex::ParallelDescriptor::ReduceRealMax(ts);
            fr_t_step = static_cast<Set::Scalar>(ts);
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
        // Perf (2026-09-18): phase, grain_id and kappa_phase are written only
        // by InitializeMicrostructure (Initialize / Regrid), so their ghost
        // layers are filled once after each (re)initialisation instead of
        // every step. (D, refilled every step until P2-27, is deleted with the damage lineage.)
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
        // P2-28: and under the face list's evaluate-only mode (nothing is removed;
        // `removed` reads 0, the list and its column tops are built as with removal)
        const bool void_aware = list_state();
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
        // D2q-1: the beam path counts SOLID cells, not geometric distance.
        // Only needed when the beam is actually on (BeamSource returns 0 at
        // P = 0), which keeps the large beam-free runs free of the gather.
        const amrex::Vector<int> solid_above_col =
            (void_aware && P_now > 0.0 && !k_top_col.empty())
            ? ColumnSolidAboveCount(lev, k_top_col) : amrex::Vector<int>();
        const int* sab    = solid_above_col.data();
        const int  nz_dom = dhi[2] - dlo[2] + 1;
        const int  klo_d  = dlo[2];
        // P2-21: the column optical depth, through solid cells only from each
        // column's top solid cell when void-aware (the solid-count convention),
        // else from the domain top (unbuilt unless the loss table varies and the
        // beam is on: then the kernel keeps the closed form)
        const BeamTau btau = BuildBeamTau(lev, P_now, void_aware ? k_top_col : amrex::Vector<int>());
        const TauView tv = ViewOf(btau, lev);

        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar bw0 = beam.omega0, blam = beam.lambda;
        const Set::Scalar bx0 = beam.x0, by0 = beam.y0, bz0 = beam.z0;
        const Set::Scalar bR  = beam.R;
        const Set::Scalar bn  = beam.n_idx;   // P2-22: in-rock spread
        const int         bprof = beam.profile;
        const Set::Scalar brad  = beam.radius;
        const Set::Scalar b_alpha = beam.alpha;   // P2-20: computed from beam.eps_r, beam.tan_delta, beam.lambda
        const Set::Scalar eps_lo = eps_rad, eps_hi = eps_rad_hi, T_eps_b = T_eps_break;
        const Set::Scalar sig_SB = sigma_SB, h_c = h_conv, T_a = T_amb;
        const Set::Scalar inv_dz = 1.0 / DX[2];
        const Set::Scalar inv_dx2 = 1.0 / (DX[0] * DX[0]);
        const Set::Scalar inv_dy2 = 1.0 / (DX[1] * DX[1]);
        const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);
        const Set::Scalar h_gb0   = h_gb_conductance;
        const Set::Scalar dx_ax[3] = {DX[0], DX[1], DX[2]};
        const Set::Scalar inv_dx  = 1.0 / DX[0];

        // Step 11b: the Hu onset validation uses a homogeneous one-phase
        // microstructure (until P2-27 so its DP damage fields existed). Apply the same
        // surface patch here as the deleted material kernel did (P2-29): prescribed_T injects the
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
        // P2-8 F1: the face list, after the column tops and BEFORE the march
        // (BuildFlameColumns), which bills every listed face.
        // P2-23: the list runs WITHOUT a flame too. `owner` is kept (the finest
        // level owns the surface).
        // P2-28: and with a PRESCRIBED patch: its top faces carry the half-cell
        // Dirichlet term (FaceListHeat sets T_presc), billed in the face loop below.
        const bool fl_on = surface_patch.face_list && owner && void_aware;
        if (fl_on) BuildFaceList(lev, k_top_col);
        if (sp_on) BuildFlameColumns(lev, time, dt, k_top_col);
        // P2-23: with the patch OFF BuildFlameColumns -- which bills the list's
        // heat, gives every listed cell its T_s and every face its gas
        // temperature -- does not run, so its LIST half runs here, in its order,
        // without the flame: h = 0 on every face, T_s, T_gas = 0.
        else if (fl_on)
        {
            FaceListHeat(lev, 0.0, time, k_top_col);
            FrAssignTs(lev);
            FaceListTg(lev);
        }
        const Set::Scalar* fh  = sp_on ? flame_h.data()  : nullptr;
        const Set::Scalar* fTg = sp_on ? flame_Tg.data() : nullptr;
        // P2-8 F5: the list every kernel face reads.
        const FaceEntry* flp = fl_on ? fl_faces.data() : nullptr;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face;

        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> T_old = temp_old_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> H_old = H_old_mf   .Patch(lev, mfi);
            Set::Patch<Set::Scalar>       H     = H_mf       .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> gid   = grain_id_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> kp    = kappa_mf   .Patch(lev, mfi);
            Set::Patch<const Set::Scalar> surf  = surface_mf.Patch(lev, mfi);
            // When void_aware is false, rem is bound to surface_mf as a
            // harmless placeholder that the kernel never reads.
            Set::Patch<const Set::Scalar> rem   =
                void_aware ? removed_mf.Patch(lev, mfi) : surface_mf.Patch(lev, mfi);
            amrex::Array4<const int> fli;
            if (fl_on) fli = fl_index->const_array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tc = T_old(i, j, k);
                Set::Scalar kc = kp(i, j, k);
                const bool self_void = void_aware && rem(i, j, k) > 0.5;
                if (self_void) kc = 0.0;
                int gc = static_cast<int>(gid(i, j, k));

                auto k_face = [&](int di, int dj, int dk, int axis) {
                    int ii = i + di, jj = j + dj, kk = k + dk;
                    // Domain-boundary face: kappa_mf has no BC
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
                    Set::Scalar k_n = kp(ii, jj, kk);
                    int g_n = static_cast<int>(gid(ii, jj, kk));   // Step 6c: the GB hook keys on grain_id
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
                        // P2-27: h_gb0 (1 - D_face) with D = 0 on every deck since
                        // the damage lineage was retired
                        inv += 1.0 / (h_gb0 * dx_ax[axis]);
                    }
                    Set::Scalar kf = 1.0 / inv;
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
                        // ⚠⚠ P2-28: THE HALF-CELL DIRICHLET TERM. The held temperature sits on the
                        // cell's top FACE, dz/2 from its centre: conductance 2 kc / dz, so the
                        // divergence term is 2 kc (T_patch - T_c) / dz^2. Until P2-28 this read
                        // kc (...) / dz^2 -- the held value a full cell away, half the boundary
                        // conductance (-12 to -20 % in the top cell's rise at 3.125 mm; first
                        // order, so refinement hid it). Under the face list the term is the
                        // LIST's (T_presc on the top face, the loop below), so it is not added
                        // here twice.
                        if (sp_presc)
                        {
                            if (!fl_on) div_kgrad += 2.0 * kc * (sp_T - Tc) * inv_dz2;
                        }
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
                                    q_in, q_out, q_rob, q_los);
                // P2-8 F4: the column-top +z face keeps the surface form, x w
                // (gain AND losses: both scale with the face's area).
                int fl_f = 0, fl_n = 0;
                if (fl_on && !self_void) { fl_f = fli(i, j, k, 0); fl_n = fli(i, j, k, 1); }
                // P2-9: a LISTED cell's every face -- the top included -- is
                // heated below at the cell's one T_s, x w (FrFaceFlux).
                if (fl_n > 0) { q_in = 0.0; q_out = 0.0; }
                if (in_patch) div_kgrad += q_in * inv_dz;


                // P2-9 R1: EVERY listed face -- top, up, side, down -- at its
                // OWN h (F2), the CELL's T_s, x its weight (F3). Risers are
                // part of the surface, not surfaces of their own.
                for (int q = fl_f; q < fl_f + fl_n; ++q)
                {
                    // P2-28: a PRESCRIBED face -- the column term's own expression (2 kc
                    // (T_presc - T_c) / dz^2, the same kc, T_c and evaluation time), x w.
                    // ⚠ w FIRST, in the column term's product order: the compiler fuses the
                    // last multiply with the add (FMA), so w (P2-9's area weight; 1.0 exactly
                    // on a flat top) must not round the product before it -- a parenthesised
                    // w * (...) differed from the column term by one ulp in H at 8 of 76 cells.
                    if (flp[q].T_presc > 0.0)
                    {
                        div_kgrad += flp[q].w * 2.0 * kc * (flp[q].T_presc - Tc) * inv_dz2;
                        continue;
                    }
                    Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                    FrFaceFlux(flp[q], fTg ? flp[q].Tg : Set::Scalar(0.0), kc, Tc, DX[0],
                               a_in, a_out, a_rob, a_los);
                    div_kgrad += flp[q].w * (a_in - a_out) * inv_dx;
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
                                                       tv.R(i, j, bR), bprof, brad, b_alpha, x, y, z,
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
                                                         : Set::Scalar(-1.0),
                                                       tv.t(i, j, k), tv.b(i, j, k), bn);

                const Set::Scalar q_loss = q_out;

                // Conservative H update; T is recovered in InvertHtoT
                // using the per-phase Numeric::Material::Table. Detached
                // (void) cells exchange no heat and get no source, so their
                // H stays frozen at its value when voided (removal never
                // writes H).
                H(i, j, k) = H_old(i, j, k) + dt * (div_kgrad + Q - q_loss * inv_dz);
            });
        }

        // E1: standing beam-deposition closure counter (thermo.dat).
        if (P_now > 0.0) BeamEnergyBalance(lev, time, dt, k_top_col, btau);
        // S1: whole-domain energy ledger (thermo.dat), after the H update.
        // P2-8: the per-class robin power, the kernel's terms summed by class.
        if (fl_on) FaceListDiagnostics(lev);
        if (ledger_enabled) EnergyLedger(lev, time, dt, k_top_col);
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
            // P2-28: a PRESCRIBED patch on the face list (the only constant-flame list: a
            // convective list needs the h_expr map, Parse) -- the list's half, in the
            // flameless branch's order: every face's h and T_presc, every cell's T_s,
            // every face's gas temperature
            if (surface_patch.face_list)
            {
                FaceListHeat(lev, 0.0, time, k_top_col);
                FrAssignTs(lev);
                FaceListTg(lev);
            }
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
                // set by JetEnthalpyMarch (T_ent here, a harmless finite value;
                // above the plane the march never visits the column, so it is
                // final there). P2-26 (census B3): the D2v above-plane pair and
                // the D2aa pad shield that reached this site are deleted.
                const Set::Scalar hj = (s > 0.0)
                                     ? surface_patch.h_f(x, y, r, s, time)
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
        // P2-8 F2: every listed face's h, at its centre, BEFORE the march.
        if (surface_patch.face_list)
        {
            FaceListHeat(lev, z_n, time, k_top_col);
            // P2-9 R1: every listed cell's T_s, before the march bills them
            FrAssignTs(lev);
        }
        if (surface_patch.jet_enthalpy) JetEnthalpyMarch(lev, time, k_top_col);
        if (surface_patch.face_list) FaceListTg(lev);
    }

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
        // P2-29: pads is the feet rule (max / mean, closed with no live deck, are deleted)
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

    // The solid fraction at first use: 1 on solid, 0 on void (after any seed).
    void FrInitFraction(int lev)
    {
        if (fr_f_init) return;
        for (amrex::MFIter mfi(*fr_f_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<Set::Scalar> const& ff = (*fr_f_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& rm = (*removed_mf[lev]).const_array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                ff(i, j, k) = (rm(i, j, k) > 0.5) ? 0.0 : 1.0; });
        }
        fr_f_init = true;
    }

    // The fraction's ghost rule: box to box and periodic by FillBoundary; out
    // of the domain, ABOVE the top is gas (0); a SYMMETRY PLANE (the patch
    // centre lies on that side) is mirrored; every other side takes the
    // nearest in-domain cell.
    void FillFractionGhosts(int lev)
    {
        amrex::MultiFab& Fm = *fr_f_mf[lev];
        Fm.FillBoundary(geom[lev].periodicity());
        const amrex::Box& dom = geom[lev].Domain();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar* PHI = geom[lev].ProbHi();
        const Set::Scalar tol = 1.0e-9 * (PHI[0] - PLO[0]);
        const bool mir_lo[3] = {std::abs(surface_patch.x0 - PLO[0]) <= tol,
                                std::abs(surface_patch.y0 - PLO[1]) <= tol, false};
        const bool mir_hi[3] = {std::abs(surface_patch.x0 - PHI[0]) <= tol,
                                std::abs(surface_patch.y0 - PHI[1]) <= tol, false};
        const bool per[3] = {geom[lev].isPeriodic(0), geom[lev].isPeriodic(1),
                             geom[lev].isPeriodic(2)};
        for (amrex::MFIter mfi(Fm, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box gb = mfi.fabbox();
            amrex::Array4<Set::Scalar> const& a = Fm.array(mfi);
            const amrex::Dim3 lo = amrex::lbound(gb), hi = amrex::ubound(gb);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int p[3] = {i, j, k};
                bool out = false, gas = false;
                for (int d = 0; d < 3; ++d)
                {
                    if (per[d]) continue;
                    const int dl = dom.smallEnd(d), dh = dom.bigEnd(d);
                    if (p[d] < dl) { out = true; p[d] = mir_lo[d] ? 2 * dl - 1 - p[d] : dl; }
                    else if (p[d] > dh)
                    {
                        out = true;
                        if (d == 2) gas = true;
                        p[d] = mir_hi[d] ? 2 * dh + 1 - p[d] : dh;
                    }
                }
                if (!out) continue;
                if (gas) { a(i, j, k) = 0.0; continue; }
                if (!gb.contains(amrex::IntVect(AMREX_D_DECL(p[0], p[1], p[2]))))
                    Util::Abort(INFO, "FillFractionGhosts: a ghost's source cell (", p[0], ", ",
                                p[1], ", ", p[2], ") is outside its own fab -- a box narrower "
                                "than 3 cells at a domain side");
                a(i, j, k) = a(p[0], p[1], p[2]);
            }
        }
    }

    // F3 -- THE AREA RULE. The one weight function.
    static Set::Scalar FaceWeight(const Set::Vector& g, Set::Scalar& nx,
                                  Set::Scalar& ny, Set::Scalar& nz, int& fb)
    {
        const Set::Scalar gn = std::sqrt(g(0) * g(0) + g(1) * g(1) + g(2) * g(2));
        if (!(gn > 0.0)) { nx = ny = nz = 0.0; fb = 1; return 1.0; }
        nx = -g(0) / gn; ny = -g(1) / gn; nz = -g(2) / gn;
        fb = 0;
        return 1.0 / (std::abs(nx) + std::abs(ny) + std::abs(nz));
    }

    // F2 -- THE HEAT RULE. The one heat resolver: h at a face centre.
    Set::Scalar FaceH(Set::Scalar x, Set::Scalar y, Set::Scalar s, Set::Scalar time)
    {
        const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
        const Set::Scalar r2 = dxp * dxp + dyp * dyp;
        const Set::Scalar r = std::sqrt(r2);
        return FaceHr(x, y, r, s, time);
    }
    // P2-11b: the same with r given -- face_heat = path passes r := s_path,
    // s := H.
    Set::Scalar FaceHr(Set::Scalar x, Set::Scalar y, Set::Scalar r, Set::Scalar s,
                       Set::Scalar time)
    {
        Set::Scalar h = 0.0;
        if (surface_patch.jet_enthalpy)
            h = (s > 0.0) ? surface_patch.h_f(x, y, r, s, time) : Set::Scalar(0.0);
        else
            h = surface_patch.h_f(x, y, r, s, time);
        if (!(std::isfinite(h) && h >= 0.0))
            Util::Abort(INFO, "surface_patch.face_list: h_expr must be finite and >= 0 "
                        "at every face centre; got h = ", h, " at x = ", x, ", y = ", y,
                        ", s = ", s);
        return h;
    }

    // F1 + F3 -- THE LIST BUILDER. Once per step, before BuildFlameColumns.
    void BuildFaceList(int lev, const amrex::Vector<int>& k_top_col_in)
    {
        fl_ready = false;
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1), klo = domain.smallEnd(2);
        const int ihi = domain.bigEnd(0),   jhi = domain.bigEnd(1),   khi = domain.bigEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        // cubic cells: checked once at Parse (CheckInterlocks, P2-27)
        const amrex::Vector<int> k_top_col =
            k_top_col_in.empty() ? ColumnTopSolid(lev) : k_top_col_in;
        const Set::Scalar r2max = surface_patch.radius * surface_patch.radius;

        // P2-9 rev 2 item 4: F3's field is the SOLID FRACTION (3 ghosts, the
        // fraction's own ghost rule), and its normal is the smoothed one.
        const amrex::BoxArray& ba = removed_mf[lev]->boxArray();
        const amrex::DistributionMapping& dm = removed_mf[lev]->DistributionMap();
        FrInitFraction(lev);
        FillFractionGhosts(lev);

        // the local faces, packed FL_NF doubles each
        std::vector<Set::Scalar> loc;
        for (amrex::MFIter mfi(*removed_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& rem = (*removed_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& Sa  = (*fr_f_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& T_old = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& kp  = (*kappa_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (rem(i, j, k) > 0.5) continue;
                const Set::Scalar x = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (dxp * dxp + dyp * dyp > r2max) continue;
                bool open[6];
                open[0] = (i > ilo) && rem(i - 1, j, k) > 0.5;
                open[1] = (i < ihi) && rem(i + 1, j, k) > 0.5;
                open[2] = (j > jlo) && rem(i, j - 1, k) > 0.5;
                open[3] = (j < jhi) && rem(i, j + 1, k) > 0.5;
                open[4] = (k > klo) && rem(i, j, k - 1) > 0.5;
                open[5] = (k == khi) || rem(i, j, k + 1) > 0.5;
                if (!(open[0] || open[1] || open[2] || open[3] || open[4] || open[5]))
                    continue;
                const int col = (j - jlo) * nx + (i - ilo);
                // P2-9 rev 2: THE normal (FrSmoothedGradient), the one function
                // every consumer's n and w come from
                const Set::Vector g = FrSmoothedGradient(Sa, i, j, k, DX[0]);
                Set::Scalar fnx = 0.0, fny = 0.0, fnz = 0.0;
                int fb = 0;
                const Set::Scalar w = FaceWeight(g, fnx, fny, fnz, fb);
                // the conductivity exactly as the kernel forms it (kappa_phase;
                // the damage factor exp(-a_D D) is retired, P2-27)
                const Set::Scalar K = kp(i, j, k);
                for (int d = 0; d < 6; ++d)
                {
                    if (!open[d]) continue;
                    const int cls = (d == 5) ? ((k == k_top_col[static_cast<std::size_t>(col)])
                                                ? FL_TOP : FL_UP)
                                  : (d == 4) ? FL_DOWN : FL_SIDE;
                    const Set::Scalar pk[FL_NF] = {
                        Set::Scalar(i), Set::Scalar(j), Set::Scalar(k), Set::Scalar(d),
                        Set::Scalar(cls), Set::Scalar(col), Set::Scalar(fb), w,
                        fnx, fny, fnz, T_old(i, j, k), K};
                    loc.insert(loc.end(), pk, pk + FL_NF);
                }
            }
        }
        // gather: disjoint rank segments, summed -- identical on every rank
        const int np = amrex::ParallelDescriptor::NProcs();
        const int me = amrex::ParallelDescriptor::MyProc();
        std::vector<int> cnt(static_cast<std::size_t>(np), 0);
        cnt[static_cast<std::size_t>(me)] = static_cast<int>(loc.size() / FL_NF);
        amrex::ParallelDescriptor::ReduceIntSum(cnt.data(), np);
        std::size_t off = 0, tot = 0;
        for (int p = 0; p < np; ++p)
        {
            if (p == me) off = tot;
            tot += static_cast<std::size_t>(cnt[static_cast<std::size_t>(p)]);
        }
        std::vector<Set::Scalar> all(tot * FL_NF, 0.0);
        std::copy(loc.begin(), loc.end(), all.begin() + static_cast<std::ptrdiff_t>(off * FL_NF));
        if (tot > 0)
            amrex::ParallelDescriptor::ReduceRealSum(all.data(), static_cast<int>(tot * FL_NF));
        fl_faces.assign(tot, FaceEntry());
        for (std::size_t q = 0; q < tot; ++q)
        {
            const Set::Scalar* p = all.data() + q * FL_NF;
            FaceEntry& e = fl_faces[q];
            e.i = static_cast<int>(p[0]); e.j = static_cast<int>(p[1]);
            e.k = static_cast<int>(p[2]); e.dir = static_cast<int>(p[3]);
            e.cls = static_cast<int>(p[4]); e.col = static_cast<int>(p[5]);
            e.fb = static_cast<int>(p[6]); e.w = p[7];
            e.nx = p[8]; e.ny = p[9]; e.nz = p[10]; e.T = p[11]; e.K = p[12];
        }
        std::sort(fl_faces.begin(), fl_faces.end(),
                  [](const FaceEntry& a, const FaceEntry& b) {
                      if (a.col != b.col) return a.col < b.col;
                      if (a.k != b.k) return a.k < b.k;
                      return a.dir < b.dir; });
        // per column: the entry range
        fl_col_first.assign(ncols, 0);
        fl_col_count.assign(ncols, 0);
        Set::Scalar n_cls[4] = {0.0, 0.0, 0.0, 0.0}, n_fb = 0.0;
        for (std::size_t q = 0; q < fl_faces.size(); ++q)
        {
            const FaceEntry& e = fl_faces[q];
            const std::size_t c = static_cast<std::size_t>(e.col);
            if (fl_col_count[c] == 0) fl_col_first[c] = static_cast<int>(q);
            ++fl_col_count[c];
            n_cls[e.cls] += 1.0;
            if (e.fb) n_fb += 1.0;
        }
        fl_n_top = n_cls[FL_TOP]; fl_n_up = n_cls[FL_UP];
        fl_n_side = n_cls[FL_SIDE]; fl_n_down = n_cls[FL_DOWN];
        fl_w_fallback = n_fb;
        // per cell: the entry range (a cell's faces are contiguous by the sort)
        if (!fl_index || fl_index->boxArray() != ba || fl_index->DistributionMap() != dm)
            fl_index = std::make_unique<amrex::iMultiFab>(ba, dm, 2, 0);
        fl_index->setVal(0);
        if (fl_w_mf[lev]) fl_w_mf[lev]->setVal(0.0);
        for (amrex::MFIter mfi(*fl_index, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<int> const& fi = fl_index->array(mfi);
            amrex::Array4<Set::Scalar> fw;
            if (fl_w_mf[lev]) fw = (*fl_w_mf[lev]).array(mfi);
            for (std::size_t q = 0; q < fl_faces.size(); ++q)
            {
                const FaceEntry& e = fl_faces[q];
                if (!bx.contains(amrex::IntVect(AMREX_D_DECL(e.i, e.j, e.k)))) continue;
                if (fi(e.i, e.j, e.k, 1) == 0) fi(e.i, e.j, e.k, 0) = static_cast<int>(q);
                fi(e.i, e.j, e.k, 1) += 1;
                if (fw) fw(e.i, e.j, e.k) += e.w;
            }
        }
        fl_ready = true;
    }

    // F2 over the list. Called by BuildFlameColumns once z_n is known and
    // BEFORE the march, which bills every face at this h.
    void FaceListHeat(int lev, Set::Scalar z_n, Set::Scalar time,
                      const amrex::Vector<int>& k_top_col)
    {
        if (!fl_ready) return;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        long mism = 0;
        // P2-23: no flame (the patch off): every listed face is billed at h = 0
        // -- its radiative and convective losses only. No FaceH, and no
        // flame_h to compare against (it is empty).
        if (!surface_patch.enabled)
        {
            for (auto& e : fl_faces) e.h = 0.0;
        }
        // P2-28: a PRESCRIBED patch -- no flame: every face h = 0, and a TOP face (the
        // column top's +z face) whose column centre lies inside the patch radius holds
        // T_expr(t), the column term's own test (x, y at the cell centre, r^2 <= R^2).
        // Side and roof faces are never prescribed: the patch is a flat top boundary.
        else if (surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature)
        {
            const Set::Scalar Tp = surface_patch.T_at(time);
            const Set::Scalar r2 = surface_patch.radius * surface_patch.radius;
            for (auto& e : fl_faces)
            {
                e.h = 0.0;
                e.T_presc = -1.0;
                if (e.cls != FL_TOP) continue;
                const Set::Scalar x = PLO[0] + (Set::Scalar(e.i) + 0.5) * DX[0];
                const Set::Scalar y = PLO[1] + (Set::Scalar(e.j) + 0.5) * DX[1];
                const Set::Scalar dxp = x - surface_patch.x0, dyp = y - surface_patch.y0;
                if (dxp * dxp + dyp * dyp <= r2) e.T_presc = Tp;
            }
            for (const auto& e : fl_faces)
                if (e.cls == FL_TOP) flame_h[static_cast<std::size_t>(e.col)] = e.h;
        }
        // P2-11b: face_heat = path -- every face from the surface path; the
        // column's flame_h (read only by column diagnostics under the list)
        // takes its top face's value, so fl_h_top_mismatch stays an identity.
        else if (surface_patch.face_heat_path)
        {
            FaceListPath(lev, time, k_top_col);
            for (const auto& e : fl_faces)
                if (e.cls == FL_TOP) flame_h[static_cast<std::size_t>(e.col)] = e.h;
        }
        else
        for (auto& e : fl_faces)
        {
            Set::Scalar ox = 0.5, oy = 0.5;
            if (e.dir == 0) ox = 0.0; else if (e.dir == 1) ox = 1.0;
            if (e.dir == 2) oy = 0.0; else if (e.dir == 3) oy = 1.0;
            // ⚠ the column-centre and top-face expressions are BuildFlameColumns'
            // own, token for token, so the top face reproduces flame_h exactly.
            const Set::Scalar x = (ox == 0.5) ? PLO[0] + (Set::Scalar(e.i) + 0.5) * DX[0]
                                              : PLO[0] + (Set::Scalar(e.i) + ox) * DX[0];
            const Set::Scalar y = (oy == 0.5) ? PLO[1] + (Set::Scalar(e.j) + 0.5) * DX[1]
                                              : PLO[1] + (Set::Scalar(e.j) + oy) * DX[1];
            Set::Scalar z_f;
            if (e.dir == 5)      z_f = PLO[2] + Set::Scalar(e.k + 1) * DX[2];
            else if (e.dir == 4) z_f = PLO[2] + Set::Scalar(e.k) * DX[2];
            else                 z_f = PLO[2] + (Set::Scalar(e.k) + 0.5) * DX[2];
            e.h = FaceH(x, y, z_n - z_f, time);
            if (e.cls == FL_TOP
                && !(e.h == flame_h[static_cast<std::size_t>(e.col)])) ++mism;
        }
        fl_h_top_mismatch = static_cast<Set::Scalar>(mism);
        if (fl_h_mf[lev])
        {
            fl_h_mf[lev]->setVal(0.0);
            for (amrex::MFIter mfi(*fl_index, false); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const int> const& fi = fl_index->const_array(mfi);
                amrex::Array4<Set::Scalar> const& fh = (*fl_h_mf[lev]).array(mfi);
                const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                    for (int q = fi(i, j, k, 0); q < fi(i, j, k, 0) + fi(i, j, k, 1); ++q)
                        fh(i, j, k) = std::max(fh(i, j, k), fl_faces[static_cast<std::size_t>(q)].h);
            }
        }
    }

    // The per-class robin power, from the replicated list with the kernel's
    // own inputs (T_old, kc, each face's h and gas temperature) -- the kernel's
    // deposit term by term, summed by class. Every rank computes the same.
    void FaceListDiagnostics(int lev)
    {
        fl_P_top = fl_P_up = fl_P_side = fl_P_down = 0.0;
        if (!fl_ready) return;
        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar A = DX[0] * DX[1];
        Set::Scalar P[4] = {0.0, 0.0, 0.0, 0.0};
        for (const auto& e : fl_faces)
        {
            Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
            // P2-28: a PRESCRIBED face -- the kernel's half-cell Dirichlet flux 2 K (T_presc - T) / dz
            if (e.T_presc > 0.0)
            {
                P[e.cls] += e.w * (2.0 * e.K * (e.T_presc - e.T) / DX[2]) * A;
                continue;
            }
            // P2-9: the kernel's own term, FrFaceFlux, for every class
            // (P2-11b: at the face's own gas temperature)
            FrFaceFlux(e, e.Tg, e.K, e.T, DX[0], q_in, q_out, q_rob, q_los);
            P[e.cls] += e.w * q_rob * A;
        }
        fl_P_top = P[FL_TOP]; fl_P_up = P[FL_UP]; fl_P_side = P[FL_SIDE]; fl_P_down = P[FL_DOWN];
        const Set::Scalar Ps = fl_P_top + fl_P_up + fl_P_side + fl_P_down;
        // P2-23: the march's identity exists only with a march (the patch on);
        // flameless, jet_P_face is never set and the identity is vacuous.
        // P2-28: a prescribed patch has no flame and no march either
        const bool presc = surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        fl_id_march = (surface_patch.enabled && !presc && std::abs(Ps) > 0.0)
                    ? std::abs(jet_P_face - Ps) / std::abs(Ps) : 0.0;
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
    static void SurfaceCellFlux(bool in_patch, bool face_form, Set::Scalar surf_c,
                                Set::Scalar h, Set::Scalar T_gas,
                                Set::Scalar kc, Set::Scalar dz, Set::Scalar Tc,
                                Set::Scalar eps_lo, Set::Scalar eps_hi, Set::Scalar T_eps_b,
                                Set::Scalar sig_SB, Set::Scalar h_c, Set::Scalar T_a,
                                Set::Scalar& q_in, Set::Scalar& q_out,
                                Set::Scalar& q_robin, Set::Scalar& q_loss,
                                Set::Scalar T_surf_in = 0.0,
                                Set::Scalar* T_face_out = nullptr)
    {
        // T_face_out (P1f deliverable 0): the temperature of the face the Robin
        // flux was evaluated at, read by the face list (FaceRemoval.H). `T_surf_in`
        // is an INPUT and the face form's `Tf` is otherwise a local, so this is
        // the only way out for it.
        //
        // ⚠⚠ EACH BRANCH WRITES THE TEMPERATURE ITS OWN q_robin USED, and the CELL
        // FORM WRITES NOTHING. There q_robin = h*(T_gas - Tc) is built on the CELL
        // AVERAGE, which is not a face temperature; writing Tc would hand the
        // reader a plausible number that is the wrong quantity.
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
        // ⚠ A caller that passes nothing gets T_surf_in = 0, so `cell` and
        // `face` stay byte-identical (the face list passes its cell's T_s).
        if (in_patch && T_surf_in > 0.0)
        {
            const Set::Scalar Ts    = T_surf_in;
            const Set::Scalar eps_s = (Ts <= T_eps_b) ? eps_lo : eps_hi;
            q_robin = h * (T_gas - Ts);
            q_loss  = eps_s * sig_SB * (Ts * Ts * Ts * Ts - T_a * T_a * T_a * T_a)
                    + h_c * (Ts - T_a);
            q_in    = q_robin - q_loss;
            q_out   = 0.0;
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
            q_in    = G * (Tf - Tc);
            q_out   = 0.0;
            q_robin = h * (T_gas - Tf);
            q_loss  = q_robin - q_in;
            // P1f: the Newton's own root -- the temperature this q_robin used.
            // ⚠ Inferred on a LINEAR half-cell (G = 2k/dz) while the criterion
            // reads an exponential; the two disagree once ell_lat < dz/2, which
            // criterion (k2) scores rather than assumes away.
            if (T_face_out) *T_face_out = Tf;
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
    //   ledger_err_raw = (E_H - E_H0 - E_src)/max(|E_src|, 1e-30)
    // P_beam reuses the E1 BeamEnergyBalance sum (same BeamSource calls);
    // P_robin/P_loss re-evaluate SurfaceCellFlux with the H update's inputs
    // (prescribed_T patch flux, if any, is counted in P_robin).
    // EXACT only with Neumann-0 outer temperature BCs on every face: a
    // Dirichlet/Robin domain BC puts a flux through the boundary faces that
    // the ledger does not count. The ledger state (ledger_started, E_H0,
    // E_src) is NOT checkpointed: a restarted run re-baselines at its first
    // step. E_H is a compensated (Neumaier) sum so the
    // difference E_H - E_H0 keeps the small per-step increments.
    // surface_missing_cols counts columns that have a mask-top solid cell but
    // no surface_mf cell (the phi-window dropout at phi_top ~ dz ties).
    void EnergyLedger(int lev, Set::Scalar time, Set::Scalar dt,
                      const amrex::Vector<int>& k_top_col)
    {
        if (lev != 0) return;
        // P2-28: the kernel's own switch (list_state: the evaluate-only list included)
        const bool void_aware  = list_state();
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

        const bool        sp_on    = surface_patch.enabled;
        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
        const bool        sp_face  = surface_patch.robin_form ==
                                     SurfacePatch::RobinForm::Face;
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
        // Per-column side power, for the wall/riser split that no plotfile can
        // give (38 fields, none of them a side flux). Sized only when the key
        // is on.
        const bool sidecsv = surface_patch.side_power_interval > 0.0;
        std::vector<Set::Scalar> sidecol;
        // P2-8: under the face list the ledger mirrors the LIST, face for face.
        // P2-23: the kernel's switch (the ledger is single-level, so `owner` is
        // true). P2-28: a prescribed patch is on the list too (its top faces' T_presc)
        const bool        fl_on_l = surface_patch.face_list && void_aware && fl_ready;
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
            amrex::Array4<const Set::Scalar> const& surf  = (*surface_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& rem   =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
            amrex::Array4<const int> fli_l;
            if (fl_on_l) fli_l = fl_index->const_array(mfi);
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
                const Set::Scalar kc = kp(i, j, k);
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
                        // P2-28: the kernel's half-cell Dirichlet term, mirrored; under the
                        // face list it is the list's (T_presc, below), not added twice
                        if (sp_presc) { if (!fl_on_l) P_robin += 2.0 * kc * (sp_T - Tc) * inv_dz2 * V; }
                        else          in_patch = true;
                    }
                }
                Set::Scalar q_in = 0.0, q_out = 0.0, q_rob = 0.0, q_los = 0.0;
                SurfaceCellFlux(in_patch, sp_face, surf_c,
                                sp_on ? flame_h[col]  : Set::Scalar(0.0),
                                sp_on ? flame_Tg[col] : Set::Scalar(0.0), kc, dz, Tc,
                                eps_rad, eps_rad_hi, T_eps_break, sigma_SB, h_conv, T_amb,
                                q_in, q_out, q_rob, q_los);
                // P2-8 F4: the top face x w, exactly as the kernel weights it.
                int fl_f = 0, fl_n = 0;
                if (fl_on_l) { fl_f = fli_l(i, j, k, 0); fl_n = fli_l(i, j, k, 1); }
                // P2-9: a listed cell's faces all come through FrFaceFlux below.
                if (fl_n > 0) { q_rob = 0.0; q_los = 0.0; }
                P_robin += q_rob * A;
                P_loss  += q_los * A;
                // P2-9 R1: every listed face, mirrored from the kernel.
                for (int q = fl_f; q < fl_f + fl_n; ++q)
                {
                    const FaceEntry& e = fl_faces[static_cast<std::size_t>(q)];
                    // P2-28: a PRESCRIBED face -- the kernel's term x V (a top face: no side share),
                    // in the column mirror's product order (w first)
                    if (e.T_presc > 0.0)
                    {
                        P_robin += e.w * 2.0 * kc * (e.T_presc - Tc) * inv_dz2 * V;
                        continue;
                    }
                    Set::Scalar a_in = 0.0, a_out = 0.0, a_rob = 0.0, a_los = 0.0;
                    FrFaceFlux(e, e.Tg, kc, Tc, DX[0], a_in, a_out, a_rob, a_los);
                    const Set::Scalar Af = e.w * DX[0] * dz;
                    P_robin      += a_rob * Af;
                    P_loss       += a_los * Af;
                    if (e.cls == FL_TOP) continue;
                    P_robin_side += a_rob * Af;
                    if (sidecsv && col >= 0 && col < static_cast<int>(sidecol.size()))
                        sidecol[col] += a_rob * Af;
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(EH_s);
        amrex::ParallelDescriptor::ReduceRealSum(EH_c);
        amrex::ParallelDescriptor::ReduceRealSum(EHsol_s);
        amrex::ParallelDescriptor::ReduceRealSum(EHsol_c);
        amrex::ParallelDescriptor::ReduceRealSum(P_robin);
        amrex::ParallelDescriptor::ReduceRealSum(P_loss);
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
        if (fl_on_l)
        {
            const Set::Scalar Ps = fl_P_top + fl_P_up + fl_P_side + fl_P_down;
            fl_id_ledger = (std::abs(Ps) > 0.0) ? std::abs(P_robin - Ps) / std::abs(Ps) : 0.0;
        }
        // P0b witnesses for robin_form = surface -- a FLAME patch's witnesses
        // (registered only with the patch on; they read flame_h, which is empty
        // without it: P2-23)
        if (sp_on && surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
        {
            amrex::ParallelDescriptor::ReduceRealSum(P_robin_side);
            surf_side_share = (std::abs(P_robin) > 0.0)
                              ? P_robin_side / P_robin : 0.0;
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
        }
        ledger_P_loss  = P_loss;
        ledger_E_src  += dt * (P_beam + P_robin - P_loss);
        ledger_E_H     = EH_s + EH_c;
        const Set::Scalar dE = (EH_s - ledger_E_H0_s) + (EH_c - ledger_E_H0_c);
        ledger_err_raw = (dE - ledger_E_src) / std::max(std::abs(ledger_E_src), Set::Scalar(1.0e-30));
        // P0b (c''): with the flake sink, the SOLID enthalpy must close against
        // the source minus what left with the rock.
        ledger_E_H_solid = EHsol_s + EHsol_c;
        const Set::Scalar dE_sol = (EHsol_s - ledger_E_H0_sol_s)
                                 + (EHsol_c - ledger_E_H0_sol_c);
        // P2-2f: foot_mech_E is the pad-clearance sink -- enthalpy that left the
        // SOLID set mechanically. Without it this row cannot close on any leg
        // with feet + foot_body_clearance, which is one cause and not two: the
        // gap between ledger_err_raw and ledger_err_solid is the flake credit the
        // TOTAL row omits by design, and the accepted CTRL5 shows a larger gap
        // while closing here to 7.78e-04.
        // P2-24: and the vapour's enthalpy (the face list's regime 2: H V of each vaporised cell)
        ledger_err_solid = (dE_sol - (ledger_E_src - ledger_E_flake - foot_mech_E - ledger_E_vapor))
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
    // Called FIRST in Advance (where P2-1's wall-law state also ran), so the
    // criterion and the kernel/march/ledger all see ONE set by
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
    // No-op, and no traversal, unless some phase sets `law = grimsel` or (P2-29) `zhang`.
    // P2-27 stage 4: also called at the end of InitializeMicrostructure with
    // `step = false` -- a law phase's scalars are refused, so the t = 0 fields
    // come from THIS writer (the law at the start temperature) instead of the
    // scalars; the per-step counters (prop_refresh_n, prop_*_mean) are not
    // touched then, so no thermo column moves.
    void UpdatePropertyFields(int lev, bool step = true)
    {
        if (!microstructure_enabled) return;
        bool any = false;
        for (const auto& ph : phases) any = any || ph.grimsel_law || ph.zhang_law;
        if (!any) return;
        if (lev >= static_cast<int>(phase_mf.size()) || !phase_mf[lev]) return;

        const int n_phases = static_cast<int>(phases.size());
        std::vector<char> use(n_phases, 0);
        // P2-29: 2 = a zhang phase, refreshed in its THERMAL fields only (the mechanics are its scalars)
        for (int p = 0; p < n_phases; ++p) use[p] = phases[p].grimsel_law ? 1 : (phases[p].zhang_law ? 2 : 0);
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
            // the host. It is one pass over the level per step.
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                if (!use[p]) continue;
                if (use[p] == 2)
                {
                    // P2-29: Zhang -- kappa_phase is a CONDUCTIVITY [W/(m K)] (Material::kappa)
                    const Numeric::Material::Zhang& z = phases[p].zhang;
                    const Set::Scalar Tz = T(i, j, k);
                    if (!(Tz > 0.0) || !std::isfinite(Tz)) continue;
                    kapa (i, j, k) = z.kappa(Tz);
                    rho_a(i, j, k) = z.rho(Tz);
                    Cp_a (i, j, k) = z.cp(Tz);
                    continue;
                }
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
        if (!step) return;
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

    // E1: global topmost solid (removed <= 0.5) cell index per (i, j)
    // column at level lev, INT_MIN for fully removed columns. Index is
    // (j - dlo_y)*nx + (i - dlo_x). MPI-reduced, so identical on every rank
    // and independent of the z box decomposition. D2q-1: every column-top
    // scan in the code now uses this rule -- `removed` alone -- so the column path (deleted, P2-28)
    // PASS 1 no longer additionally requires phi >= 0.
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

    // P2-21: THE COLUMN OPTICAL DEPTH for a temperature-dependent loss
    // (beam.tan_delta_T). Mirrors ColumnSolidAboveCount: every solid cell's
    // a_k*dz, a_k = beam.alpha_at(T_old), into an ncols*nz vector (-1 marks a
    // void cell; each cell is owned by exactly one box, so a sum over ranks
    // gathers the column), ONE ReduceRealSum, then per column from the top solid
    // cell down the running sums tau_top (0 at the top cell) and tau_bot =
    // tau_top + a_k*dz -- the SAME a_k*dz, never recomputed from another T.
    // The CENTRAL column's temperatures (the column containing (x0, y0); it can
    // span several z-boxes, so it is never assumed to sit on one rank) ride in
    // the same vector, so the deciders read the T_old the gather used.
    // k_top_col empty => every column's top is the domain top and every cell is
    // solid (the material paths, and the microstructure path without removal);
    // the beam then enters at the domain top, so beam.z0 must be that top.
    // Also sets the thermo deciders beam_tau_max ... bc_ell_pk.
    // Cost: one (ncols*nz + nz)-double allreduce per step, ONLY on a deck whose
    // table varies and only while the beam is on (BuildBeamTau).
    BeamTau ColumnOpticalDepth(int lev, const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1);
        const int klo = domain.smallEnd(2), khi = domain.bigEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t NZ = static_cast<std::size_t>(domain.length(2));
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const Set::Scalar dz = DX[2];
        const bool void_aware = !k_top_col.empty();
        const int kmin = std::numeric_limits<int>::min();
        if (!void_aware && std::abs(beam.z0 - geom[lev].ProbHi()[2]) > 1.0e-9 * dz)
            Util::Abort(INFO, "beam.tan_delta_T: the column optical depth starts at the domain top (",
                        geom[lev].ProbHi()[2], " m); beam.z0 = ", beam.z0, " must equal it");
        auto top_of = [&](std::size_t c) { return void_aware ? k_top_col[c] : khi; };

        // a column is LIT when the beam reaches its top cell (Pi > 0 for any P > 0;
        // a Gaussian's waist only grows with depth, a uniform disk's radius is fixed)
        std::vector<char> lit(ncols, 0);
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = top_of(c);
            if (kt == kmin) continue;
            const int i = ilo + static_cast<int>(c % nx), j = jlo + static_cast<int>(c / nx);
            lit[c] = BeamIncidentIntensity(1.0, beam.omega0, beam.lambda, beam.x0, beam.y0, beam.z0,
                                           beam.profile, beam.radius,
                                           PLO[0] + (Set::Scalar(i) + 0.5) * DX[0],
                                           PLO[1] + (Set::Scalar(j) + 0.5) * DX[1],
                                           PLO[2] + (Set::Scalar(kt) + 0.5) * dz) > 0.0;
        }
        const int ic = std::min(nx - 1, std::max(0, static_cast<int>(std::floor((beam.x0 - PLO[0]) / DX[0]))));
        const int jc = std::min(ny - 1, std::max(0, static_cast<int>(std::floor((beam.y0 - PLO[1]) / DX[1]))));
        const std::size_t c0 = static_cast<std::size_t>(jc) * nx + ic;

        // [ncols*NZ: a_k dz | NZ: the central column's T | ncols: each column's TOP solid cell T (P2-24)]
        std::vector<Set::Scalar> g(ncols * NZ + NZ + ncols, 0.0);
        const std::size_t o_top = ncols * NZ + NZ;
        const Set::Scalar mf = beam.melt_factor;
        Set::Scalar amax = 0.0;
        for (amrex::MFIter mfi(*temp_old_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T = (*temp_old_mf[lev]).const_array(mfi);
            amrex::Array4<const Set::Scalar> const& rem =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*temp_old_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(j - jlo) * nx + (i - ilo);
                const std::size_t p = c * NZ + (k - klo);
                if (c == c0) g[ncols * NZ + (k - klo)] = T(i, j, k);
                if (void_aware && rem(i, j, k) > 0.5) { g[p] = -1.0; continue; }
                if (k == top_of(c)) g[o_top + c] = T(i, j, k);
                // P2-24: the melt's absorption on the liquid fraction (1 + (m - 1) f_L); m = 1 skips it
                const Set::Scalar a = (mf != 1.0)
                    ? beam.alpha_at(T(i, j, k)) * (1.0 + (mf - 1.0) * MeltFraction(T(i, j, k)))
                    : beam.alpha_at(T(i, j, k));
                g[p] = a * dz;
                if (lit[c] && a > amax) amax = a;
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(g.data(), static_cast<int>(g.size()));
        amrex::ParallelDescriptor::ReduceRealMax(amax);

        BeamTau out;
        out.top.assign(ncols * NZ, -1.0);
        out.bot.assign(ncols * NZ, -1.0);
        out.k_last.assign(ncols, kmin);
        Set::Scalar tmax = 0.0, emax = 0.0;
        for (std::size_t c = 0; c < ncols; ++c)
        {
            const int kt = top_of(c);
            if (kt == kmin) continue;
            Set::Scalar acc = 0.0, e_top = 1.0, s = 0.0;
            for (int k = kt; k >= klo; --k)
            {
                const std::size_t p = c * NZ + (k - klo);
                if (g[p] < 0.0) continue;              // void: no material, no attenuation
                out.top[p] = acc;
                out.bot[p] = acc + g[p];
                acc = out.bot[p];
                out.k_last[c] = k;
                // the kernel's own pair: e^{-tau_top} of this cell IS e^{-tau_bot}
                // of the solid cell above (the same double), so the sum telescopes
                const Set::Scalar e_bot = std::exp(-out.bot[p]);
                s += e_top - e_bot;
                e_top = e_bot;
            }
            if (!lit[c] || out.k_last[c] == kmin) continue;
            tmax = std::max(tmax, acc);
            emax = std::max(emax, std::abs(s + e_top - 1.0));
        }
        beam_tau_max     = tmax;
        beam_alpha_max   = amax;
        beam_col_err_max = emax;

        // P2-24: the column's reflection from its top solid cell's liquid fraction (only when the melt keys
        // differ from the no-op: unbuilt => every site reads the scalar R, byte for byte), and the surface
        // record of the melt: the synthetic radiometer over the DOE launcher's view, the melt spot's area
        const Set::Scalar R_s = beam.R, R_l = beam.R_melt;
        if (beam.melt_active())
        {
            out.Rcol.assign(ncols, R_s);
            for (std::size_t c = 0; c < ncols; ++c)
                if (top_of(c) != kmin) out.Rcol[c] = R_s + (R_l - R_s) * MeltFraction(g[o_top + c]);
        }
        {
            // the radiometer's view: the 32.5 mm launcher aperture (DOE DE-EE0005504 L992-1002), r <= 16.25 mm
            // from the beam axis. eps_137 = (1 - R_melt) f_L + (1 - R)(1 - f_L): the DOE's own melt absorptance
            // (0.7) and the solid's from its permittivity, both from the one dielectric description. A LOG.
            const Set::Scalar r_view = 0.01625;
            Set::Scalar n_v = 0.0, sTr = 0.0, sT = 0.0, A_m = 0.0;
            for (std::size_t c = 0; c < ncols; ++c)
            {
                if (top_of(c) == kmin) continue;
                const Set::Scalar Tt = g[o_top + c];
                const Set::Scalar fL = MeltFraction(Tt);
                if (fL > 0.0) A_m += DX[0] * DX[1];
                const int i = ilo + static_cast<int>(c % nx), j = jlo + static_cast<int>(c / nx);
                const Set::Scalar dxc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0] - beam.x0;
                const Set::Scalar dyc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1] - beam.y0;
                if (dxc * dxc + dyc * dyc > r_view * r_view) continue;
                n_v += 1.0;
                sT  += Tt;
                sTr += ((1.0 - R_l) * fL + (1.0 - R_s) * (1.0 - fL)) * Tt;
            }
            bc_T_rad  = (n_v > 0.0) ? sTr / n_v : std::numeric_limits<Set::Scalar>::quiet_NaN();
            bc_T_view = (n_v > 0.0) ? sT / n_v : std::numeric_limits<Set::Scalar>::quiet_NaN();
            bc_A_melt = A_m;
        }

        // the central column: its contiguous solid run from the top, T from the
        // gathered T_old; depths are cell centres below the top solid face
        const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
        bc_T_peak = bc_z_peak = bc_z_melt = bc_tau = bc_ell_p = bc_ell_pk = nan;
        bc_f_L_peak = bc_z_liq = bc_R_col = nan;
        const int kt0 = top_of(c0);
        if (kt0 != kmin)
        {
            // P2-24: the central column's reflection, from its top solid cell (the form every site reads)
            bc_R_col = R_s + (R_l - R_s) * MeltFraction(g[o_top + c0]);
            // P2-24: the ROCK PHASE's solidus on the microstructure path (phase.T_m left the beam deck)
            const Set::Scalar Tm = !phases.empty() ? phases[0].T_m : Set::Scalar(1.0e30);
            std::vector<Set::Scalar> Tcol;
            for (int k = kt0; k >= klo && g[c0 * NZ + (k - klo)] >= 0.0; --k)
                Tcol.push_back(g[ncols * NZ + (k - klo)]);
            const int nmax = static_cast<int>(Tcol.size()) - 1;
            int d_pk = 0;
            for (int d = 0; d <= nmax; ++d)
            {
                if (Tcol[d] > Tcol[d_pk]) d_pk = d;
                if (!std::isfinite(bc_z_melt) && Tcol[d] >= Tm) bc_z_melt = (Set::Scalar(d) + 0.5) * dz;
            }
            bc_T_peak = Tcol[d_pk];
            bc_z_peak = (Set::Scalar(d_pk) + 0.5) * dz;
            // P2-24: the liquid fraction at the peak, and the deepest fully liquid cell of the column
            bc_f_L_peak = MeltFraction(bc_T_peak);
            for (int d = 0; d <= nmax; ++d)
                if (MeltFraction(Tcol[d]) >= 1.0) bc_z_liq = (Set::Scalar(d) + 0.5) * dz;
            if (out.k_last[c0] != kmin) bc_tau = out.bot[c0 * NZ + (out.k_last[c0] - klo)];
            // P2-19's helpers on the rise above the ambient (the Zhang decks' IC);
            // the face sample is the top cell's rise (no surface reconstruction here)
            auto rise_at_d = [&](int d) { return Tcol[d] - T_amb; };
            Set::Scalar z_pk = nan;
            bc_ell_p = ProfileEFoldDepth(rise_at_d, Tcol[0] - T_amb, dz, nmax);
            ProfilePeak(rise_at_d, Tcol[0] - T_amb, dz, nmax, z_pk, bc_ell_pk);
        }
        out.built = true;
        return out;
    }

    // P2-21: the column optical depth when it applies (the table varies AND the
    // beam is on this step), else an unbuilt BeamTau (every kernel then keeps the
    // closed form) with the deciders cleared on a table deck.
    BeamTau BuildBeamTau(int lev, Set::Scalar P_now, const amrex::Vector<int>& k_top_col)
    {
        beam_P_beam = P_now;   // P2-24: the schedule's power now (the DOE's denominator)
        if (beam.alpha_varies() && P_now > 0.0) return ColumnOpticalDepth(lev, k_top_col);
        if (beam.alpha_varies())
        {
            const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
            beam_tau_max = beam_alpha_max = beam_col_err_max = 0.0;
            bc_T_peak = bc_z_peak = bc_z_melt = bc_tau = bc_ell_p = bc_ell_pk = nan;
            // P2-24: with the beam OFF the surface still cools; on a BEAM deck the radiometer and the melt spot
            // are read anyway (a top-cell gather, no optical depth). A flame deck (no beam power ever) skips it.
            if (beam.can_deposit()) SurfaceMeltRecord(lev, k_top_col);
        }
        return BeamTau();
    }

    // P2-24: THE ROCK'S LIQUID FRACTION f_L(T) -- ONE phase (SR 13): phase 0's fusion band on the microstructure
    // path (P2-29: the only path); 0 when no band is declared (T_m >= 1e29).
    // f_L = clamp((T - (T_m - w/2)) / w, 0, 1); a sharp band (w <= 0) is the step at T_m. The beam's melt factor
    // and reflection, the face list's molten / lid exclusions and the melt records all read THIS.
    Set::Scalar MeltFraction(Set::Scalar T) const
    {
        Set::Scalar Tm = 1.0e30, w = 0.0;
        if (!phases.empty()) { Tm = phases[0].T_m; w = phases[0].fus_width; }
        if (!(Tm < 1.0e29)) return 0.0;
        if (!(w > 0.0)) return (T >= Tm) ? 1.0 : 0.0;
        const Set::Scalar f = (T - (Tm - 0.5 * w)) / w;
        return (f <= 0.0) ? 0.0 : ((f >= 1.0) ? 1.0 : f);
    }

    // P2-24: the beam-OFF surface record (after the beam is switched off the melt still cools): bc_T_rad,
    // bc_T_view, bc_A_melt from each column's top solid cell, as ColumnOpticalDepth forms them
    void SurfaceMeltRecord(int lev, const amrex::Vector<int>& k_top_col)
    {
        const amrex::Box& domain = geom[lev].Domain();
        const int ilo = domain.smallEnd(0), jlo = domain.smallEnd(1), khi = domain.bigEnd(2);
        const int nx = domain.length(0), ny = domain.length(1);
        const std::size_t ncols = static_cast<std::size_t>(nx) * ny;
        const Set::Scalar* DX  = geom[lev].CellSize();
        const Set::Scalar* PLO = geom[lev].ProbLo();
        const bool void_aware = !k_top_col.empty();
        const int kmin = std::numeric_limits<int>::min();
        std::vector<Set::Scalar> Ttop(ncols, 0.0);
        for (amrex::MFIter mfi(*temp_old_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& T = (*temp_old_mf[lev]).const_array(mfi);
            const amrex::Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(j - jlo) * nx + (i - ilo);
                const int kt = void_aware ? k_top_col[c] : khi;
                if (kt != kmin && kt >= lo.z && kt <= hi.z) Ttop[c] = T(i, j, kt);
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(Ttop.data(), static_cast<int>(Ttop.size()));
        const Set::Scalar R_s = beam.R, R_l = beam.R_melt, r_view = 0.01625;   // as ColumnOpticalDepth
        Set::Scalar n_v = 0.0, sTr = 0.0, sT = 0.0, A_m = 0.0;
        for (std::size_t c = 0; c < ncols; ++c)
        {
            if ((void_aware ? k_top_col[c] : khi) == kmin) continue;
            const Set::Scalar fL = MeltFraction(Ttop[c]);
            if (fL > 0.0) A_m += DX[0] * DX[1];
            const int i = ilo + static_cast<int>(c % nx), j = jlo + static_cast<int>(c / nx);
            const Set::Scalar dxc = PLO[0] + (Set::Scalar(i) + 0.5) * DX[0] - beam.x0;
            const Set::Scalar dyc = PLO[1] + (Set::Scalar(j) + 0.5) * DX[1] - beam.y0;
            if (dxc * dxc + dyc * dyc > r_view * r_view) continue;
            n_v += 1.0;
            sT  += Ttop[c];
            sTr += ((1.0 - R_l) * fL + (1.0 - R_s) * (1.0 - fL)) * Ttop[c];
        }
        const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
        bc_T_rad  = (n_v > 0.0) ? sTr / n_v : nan;
        bc_T_view = (n_v > 0.0) ? sT / n_v : nan;
        bc_A_melt = A_m;
        bc_f_L_peak = bc_z_liq = bc_R_col = nan;
    }

    // E1: per-step closure of the microstructure-path beam deposition.
    //
    //   P_inc   = sum over lit columns of (1-R)*Pi*dx*dy at the top solid cell
    //   P_dep   = sum over solid cells of Q*dx*dy*dz  (same BeamSource call)
    //   P_trans = (1-R)*Pi*exp(-a*s_bot)*dx*dy at k == domain bottom
    //   closure : P_dep + P_trans == P_inc
    //
    // Exact (to round-off) only for CONSTANT alpha and the UNIFORM profile.
    // A temperature-dependent alpha deposits through the cumulative optical
    // depth tau = sum_k a_k*dz (P2-21: ColumnOpticalDepth), never a per-cell
    // a_k applied to the whole path above (which over-deposits). A Gaussian profile adds an
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
    // P2-21: with a built BeamTau (a varying loss table) every cell's source is
    // the optical-depth pair the kernel used and P_trans leaves at each column's
    // LOWEST solid cell with e^{-tau_bot} there; the per-column pairs telescope
    // (beam_col_err_max), but the closure keeps a residual of the size
    // (d/z_R)^2 x the centre-to-edge contrast in alpha: Pi is the cell's own, and
    // a Gaussian's Pi changes with depth -- the lateral power a column sum cannot
    // follow (studies/p2_21_loss_table/step0). Unbuilt: the closed form, as before.
    void BeamEnergyBalance(int lev, Set::Scalar time, Set::Scalar dt,
                           const amrex::Vector<int>& k_top_col,
                           const BeamTau& btau)
    {
        if (!SurfaceOwner(lev)) return;
        // P2-28: the kernel's own switch (list_state)
        const bool void_aware = list_state();
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
        const Set::Scalar b_alpha = beam.alpha;   // P2-20: computed from beam.eps_r, beam.tan_delta, beam.lambda
        // D2q-1: solid-cell count above k, as in AdvanceMicrostructure.
        const amrex::Vector<int> solid_above =
            (void_aware && P_now > 0.0 && !k_top_col.empty())
            ? ColumnSolidAboveCount(lev, k_top_col) : amrex::Vector<int>();
        const int* sab    = solid_above.data();
        const int  nz_dom = domain.length(2);
        const TauView tv = ViewOf(btau, lev);

        Set::Scalar P_inc = 0.0, P_dep = 0.0, P_trans = 0.0, n_viol = 0.0;
        for (amrex::MFIter mfi(*H_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            // Placeholders as in AdvanceMicrostructure when not void-aware.
            // (P2-20: alpha no longer depends on T, so no temperature is read
            // here -- the call is valid on every advance path.)
            amrex::Array4<const Set::Scalar> const& rem =
                void_aware ? (*removed_mf[lev]).const_array(mfi) : (*surface_mf[lev]).const_array(mfi);
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
                // P2-22: `path` is -1 off the void-aware path (BeamSource then takes the
                // geometric z0 - z as the solid path); the incident intensity must use the
                // same solid path as the deposit or P_inc and P_dep disagree (smoke: 1.43).
                const Set::Scalar s_rock = (path >= 0.0) ? path : (beam.z0 - z);
                const Set::Scalar Pi = BeamIncidentIntensity(P_now, beam.omega0, beam.lambda,
                                                             beam.x0, beam.y0, beam.z0,
                                                             beam.profile, beam.radius, x, y, z,
                                                             beam.n_idx, s_rock);
                if (Pi == 0.0) continue;

                // P2-24: the column's R (the melt's reflection over a molten top), the deck's scalar R unbuilt
                const Set::Scalar Rc = tv.R(i, j, beam.R);
                const Set::Scalar Q = BeamSource(P_now, beam.omega0, beam.lambda,
                                                 beam.x0, beam.y0, beam.z0, Rc,
                                                 beam.profile, beam.radius, b_alpha,
                                                 x, y, z, dz, path, tv.t(i, j, k), tv.b(i, j, k),
                                                 beam.n_idx);
                P_dep += Q * V;
                if (at_col_top) P_inc += (1.0 - Rc) * Pi * A;
                if (btau.built)
                {
                    // P2-21: the gathered value at the column's exit cell
                    if (k == btau.k_last[static_cast<std::size_t>(col_b)])
                        P_trans += (1.0 - Rc) * Pi * std::exp(-tv.b(i, j, k)) * A;
                }
                else if (k == klo_dom)
                {
                    const Set::Scalar a = b_alpha;
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
                                  Set::Scalar a,
                                  Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                  Set::Scalar dz_cell = 0.0,
                                  Set::Scalar dz_abs = -1.0,
                                  Set::Scalar tau_top = -1.0,
                                  Set::Scalar tau_bot = -1.0,
                                  Set::Scalar n_rock = 1.0)
    {
        if (P == 0.0) return 0.0; // perf: beam off (Pi would be 0 anyway)
        Set::Scalar dz = z0 - z;
        // P2-22: the solid path to this cell (dz when the beam origin is the
        // rock surface) sets the in-rock spread; n_rock 1 keeps the geometric dz.
        Set::Scalar path = (dz_abs >= 0.0) ? dz_abs : dz;
        Set::Scalar Pi = BeamIncidentIntensity(P, w0, lam, x0, y0, z0,
                                               profile, radius, x, y, z, n_rock, path);
        if (Pi == 0.0) return 0.0;
        // P2-21: a VARYING coefficient. The column's optical depth to this cell's
        // top and bottom faces (ColumnOpticalDepth) replaces a*s: the Beer-Lambert
        // integral of a(T) along the solid path, so the column's pairs telescope to
        // 1 - e^{-tau_N}. Pi stays the cell's own (the closed form's, unchanged).
        // tau_top < 0 (the default) is the closed form below, byte for byte.
        if (tau_top >= 0.0)
            return (1.0 - R) * Pi * (std::exp(-tau_top) - std::exp(-tau_bot)) / dz_cell;
        // P2-20: a = beam.alpha, the POWER coefficient computed from the
        // dielectric description (no temperature dependence: user, 2026-11-05)
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
    // P2-22: n_rock = sqrt(eps_r) and s_rock = the solid path from the surface to
    // this point. A Gaussian beam crossing a plane interface into a dielectric
    // keeps its waist and propagates with the Rayleigh range n*zR, i.e. the waist
    // law is evaluated at the REDUCED distance d_air + s_rock / n. Written as
    // dz - s_rock*(1 - 1/n) so that n = 1 reproduces the geometric dz BIT FOR BIT
    // (s_rock * 0.0 = 0.0; every eps_r 1 fixture is byte-identical).
    static Set::Scalar BeamIncidentIntensity(Set::Scalar P, Set::Scalar w0, Set::Scalar lam,
                                             Set::Scalar x0, Set::Scalar y0, Set::Scalar z0,
                                             int profile, Set::Scalar radius,
                                             Set::Scalar x, Set::Scalar y, Set::Scalar z,
                                             Set::Scalar n_rock = 1.0, Set::Scalar s_rock = 0.0)
    {
        Set::Scalar dz = z0 - z;
        if (dz < 0.0 || P == 0.0) return 0.0;
        Set::Scalar dz_red = dz - s_rock * (1.0 - 1.0 / n_rock);
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
            Set::Scalar w  = w0 * std::sqrt(1.0 + (dz_red / zR) * (dz_red / zR));
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
        // P2-28: the AMR surface band (surface.refine_depth / refine_height) is retired with the column path
    }

private:
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

        // ------------------------------------------------------------------
        // P2-28: spall.surface_T is retired at `reconstructed` (the table in Parse): the face list's
        // R1 IS the reconstruction; site A's cell-average and diagnostic forms are deleted.

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

        // P2-28: spall.removal_energy is retired at `flake` (the table in Parse): the face list's R3
        // IS the flake remap.

        if (!value.spall_enabled) return;

        // P2-27 (census C1a): the damage-law removal's keys (spall.damage_threshold,
        // xi_s, C_h, alpha, t_spall, prescribed_ratio.*, sample, shallow_depth) are
        // retired with it -- see the abort list and the retired-at-default table.

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

    // The removal code lives in two class-body fragments, textually included
    // here: Shared.H (the helpers and the step prologue) and FaceRemoval.H (the
    // face list, the ONE removal path since P2-28 deleted the column path).
    #include "Integrator/MMWSpalling/Shared.H"
    #include "Integrator/MMWSpalling/FaceRemoval.H"

    // Build the global-material H<->T table once. Microstructure tables
    // are built per phase in InitializeMicrostructure(); see phase_tables.

    void ConvertTtoH(int lev)
    {
        // P2-29: one path -- the per-cell phase table (this->phase_tables[pid]); the global material's table and
        // the const-alpha closed form are deleted
        const Numeric::Material::Table *ptabs = phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& H    = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Hold = (*H_old_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph = (*phase_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tv = T(i, j, k);
                Set::Scalar Hv = 0.0;
                Set::Scalar ls = 1.0, ll = 0.0, lv = 0.0;
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                Hv = ptabs[p].H_of_T(Tv);
                Set::Scalar dummyT;
                ptabs[p].Invert(Hv, dummyT, ls, ll, lv);
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
                H(i,j,k)    = Hv;
                Hold(i,j,k) = Hv;
            });
        }
    }

    void InvertHtoT(int lev)
    {
        // P2-29: one path -- the per-cell phase table
        const Numeric::Material::Table *ptabs = phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& H = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& T  = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph = (*phase_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Hv = H(i, j, k);
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                Set::Scalar Tv, ls, ll, lv;
                ptabs[p].Invert(Hv, Tv, ls, ll, lv);
                T(i,j,k) = Tv;
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
            });
        }
    }

    void InitSurfaceMask(int lev)
    {
        if (list_state())
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
            InitializeMicrostructure(lev, announce);
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
    void InitializeMicrostructure(int lev, bool announce = true)
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
                // P2-29: a zhang phase -- the same rule, one rule for both laws (Zhang's own latent
                // scalars are refused at the parse; the phase's are the ones that count)
                else if (phases[p].zhang_law)
                {
                    Numeric::Material::Zhang z = phases[p].zhang;
                    z.T_ref     = phases[p].T_ref;
                    z.T_m       = phases[p].T_m;
                    z.fus_width = phases[p].fus_width;
                    z.L_m       = phases[p].L_m;
                    z.T_vap_lo  = phases[p].T_vap_lo;
                    z.T_vap_hi  = phases[p].T_vap_hi;
                    z.L_v       = phases[p].L_v;
                    phase_tables[p].Build(z);
                }
                else
                    phase_tables[p].Build(phase_materials[p]);
                // P2-24: THE TABLE THE RUN INVERTS, printed once per build when a latent band is declared: its
                // band temperatures and enthalpies [J/m^3, H(T_ref) = 0] (b1 reads the fusion jump from here)
                if ((phases[p].T_m < 1.0e29 || phases[p].T_vap_lo < 1.0e29)
                    && amrex::ParallelDescriptor::IOProcessor())
                {
                    const Numeric::Material::Table& tb = phase_tables[p];
                    const Set::Scalar fl = phases[p].T_m - 0.5 * phases[p].fus_width;
                    const Set::Scalar fh = phases[p].T_m + 0.5 * phases[p].fus_width;
                    std::cout << std::setprecision(17)
                              << "P2-24 phase " << p << " enthalpy table [J/m^3]: H(fus_lo " << fl << " K) = "
                              << tb.H_fus_lo() << ", H(fus_hi " << fh << " K) = " << tb.H_fus_hi()
                              << ", H(vap_lo " << phases[p].T_vap_lo << " K) = " << tb.H_vap_lo()
                              << ", H(vap_hi " << phases[p].T_vap_hi << " K) = " << tb.H_vap_hi()
                              << std::setprecision(6) << std::endl;
                }
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
            // That function returns early unless SOME phase sets
            // `law = grimsel`, and skips non-Grimsel phases. Under the
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
            else if (phases[p].zhang_law)   // P2-29
            {
                rh[p]  = phases[p].zhang.rho(phases[p].T_ref);
                Cpv[p] = phases[p].zhang.cp (phases[p].T_ref);
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
        // P2-27 stage 4: a LAW phase's scalars are refused (Parse), so its property
        // fields come from the one per-cell writer -- the law at the start
        // temperature, what the first Advance writes -- not from the scalars
        // above (no-op unless some phase sets law = grimsel; no counter moves).
        UpdatePropertyFields(lev, /*step=*/false);
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

    Set::Scalar eps_rad      = 0.0;     // emissivity below T_eps_break  (-)
    Set::Scalar eps_rad_hi   = 0.0;     // emissivity above T_eps_break  (-)
    Set::Scalar T_eps_break  = 1.0e30;  // breakpoint for ε(T)           (K)
    Set::Scalar sigma_SB = 0.0; // Stefan-Boltzmann constant       (W/m²K⁴)
    Set::Scalar h_conv   = 0.0; // convective heat-transfer coeff. (W/m²K)
    Set::Scalar T_amb    = 0.0; // ambient temperature             (K)

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

    // P2-21: the deciders, thermo columns registered ONLY when
    // beam.alpha_varies() && beam.can_deposit() (every flat and every beam-OFF
    // deck keeps its header). Per step, non-extensive; set by ColumnOpticalDepth.
    Set::Scalar beam_tau_max     = 0.0;  // max over lit columns of tau_bot at the last solid cell
    Set::Scalar beam_alpha_max   = 0.0;  // max a_k over lit solid cells [1/m]
    Set::Scalar beam_col_err_max = 0.0;  // max over lit columns |sum_k (e^-tt - e^-tb) + e^-tN - 1|
    Set::Scalar bc_T_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();  // central column [K]
    Set::Scalar bc_z_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();  // [m below the top face]
    Set::Scalar bc_z_melt = std::numeric_limits<Set::Scalar>::quiet_NaN();  // shallowest T >= T_m [m]
    Set::Scalar bc_tau    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // central tau_bot at k_last
    Set::Scalar bc_ell_p  = std::numeric_limits<Set::Scalar>::quiet_NaN();  // P2-19 face e-fold [m]
    Set::Scalar bc_ell_pk = std::numeric_limits<Set::Scalar>::quiet_NaN();  // P2-19 peak e-fold [m]
    // P2-24: the melt and vapour record (same registration condition as the bc_* above)
    Set::Scalar beam_P_beam = 0.0;   // the schedule's power now [W] (the DOE's denominator; P_inc is after R)
    Set::Scalar bc_T_rad    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // synthetic radiometer eps T [K]
    Set::Scalar bc_T_view   = std::numeric_limits<Set::Scalar>::quiet_NaN();  // view mean of T_top [K]
    Set::Scalar bc_A_melt   = 0.0;   // top solid cells with f_L > 0 [m^2]
    Set::Scalar bc_f_L_peak = std::numeric_limits<Set::Scalar>::quiet_NaN();  // f_L at the central peak [-]
    Set::Scalar bc_z_liq    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // deepest central f_L >= 1 [m]
    Set::Scalar bc_R_col    = std::numeric_limits<Set::Scalar>::quiet_NaN();  // the central column's R_col [-]
    Set::Scalar bc_n_vapor     = 0.0;   // cumulative surface cells vaporised (the face list's vapour rule)
    Set::Scalar bc_n_vapor_int = 0.0;   // cumulative cell-steps an INTERIOR cell sat past H_vap_hi

    // S1: surface.follow_mask applies the convective_flame patch and the
    // surface losses at the `removed`-mask column top (ColumnTopSolid)
    // instead of the fixed domain top / the phi-window surface_mf.
    bool surface_follow_mask = false;

    // S1: whole-domain energy ledger (EnergyLedger), thermo.dat vars.
    // E_H [J] is sum(H*V) over all level-0 cells (void included), E_src [J]
    // the cumulative applied energy, P_* [W] per step. The compensated-sum
    // state (ledger_E_H0_s/c) is internal.
    bool        ledger_enabled               = false;
    // D2a: per-column flame (h, T_gas) and nozzle stand-off s for the current
    // AdvanceMicrostructure call (BuildFlameColumns). Constants path: h_conv /
    // T_flame everywhere. Read by the kernel, EnergyLedger and the face list.
    std::vector<Set::Scalar> flame_h, flame_Tg, flame_s;
    Set::Scalar patch_min_standoff = 0.0; // thermo: min s over in-patch live columns [m]
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
    // sides of the ratio. See the column path (deleted, P2-28)'s pad-clearance block for why this is a
    // booking-only term and why the field must not be written.
    Set::Scalar foot_mech_E    = 0.0;     // thermo: cumulative pad-clearance enthalpy [J]
    Set::Scalar foot_mech_cols = 0.0;     // thermo: columns cleared this step
    std::vector<int> foot_clear_k;        // per column target top k (INT_MIN = none)
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
    // ⚠⚠ RELIEVING ONE SITE RELIEVES THE WALL AND NOT THE FLOOR -- a family
    // that has cost this programme three packets. unit/lateral_flux_lock
    // asserts that every surviving stress expression takes exactly this
    // member, so no site can be missed and a new stress expression cannot
    // appear un-relieved.
    Set::Scalar sp_relief = 1.0;
    // P2-15: the buckling thickness condition (spall.buckle_span [m], 0 = off; face list only).
    // A fired line is withdrawn unless C_L x <= L_span (FaceRemoval.H, BuckleRatio). Witnesses
    // registered only under the key; per step, reset and added once in RecedeAlongLines.
    Set::Scalar sp_buckle_span = 0.0;
    // P2-23: spall.buckle_log_csv -- one row per line B evaluated, on
    // thermo-output steps, into <plot_file>_buckle_log.csv (FaceRemoval.H)
    bool sp_buckle_log_csv = false;
    bool sp_buckle_log_header_written = false;
    Set::Scalar bk_span_out = 0.0, bk_n_eval = 0.0, bk_n_block = 0.0, bk_n_inf = 0.0,
                bk_n_block_floor = 0.0, bk_n_block_wall = 0.0, bk_n_block_roof = 0.0,
                bk_max_pass = 0.0;
    long        jet_profile_last = -1; // last jet_profile_interval index written
    bool        jet_profile_header_written = false;

    bool        ledger_started               = false;
    Set::Scalar ledger_E_H0_s                = 0.0;
    Set::Scalar ledger_E_H0_c                = 0.0;
    Set::Scalar ledger_E_H                   = 0.0;
    Set::Scalar ledger_E_src                 = 0.0;
    Set::Scalar ledger_err_raw               = 0.0;
    Set::Scalar ledger_P_beam                = 0.0;
    Set::Scalar ledger_P_robin               = 0.0;
    Set::Scalar ledger_P_loss                = 0.0;
    Set::Scalar ledger_surface_missing_cols  = 0.0;

    // Step 6d / Step 15a: surface patch driver. A circular patch on the +z
    // face applies one of two heat-load modes inside the patch footprint:
    //   mode = prescribed_T (default, Step 6d) — holds T = T_expr(t) on the top
    //     FACE and injects the half-cell Dirichlet flux
    //     q = 2*kappa*(T_patch - T_cell)/dz as a divergence contribution
    //     2*kappa*(T_patch - T_cell)/dz^2 on the in-patch top row (P2-28
    //     corrected the factor 2; P2-29 deleted the material kernel that still
    //     had the half term). Under the face list the top faces carry it.
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
        // P2-28: `pinned` (the column pin rule) is retired -- it aborts by value.
        enum class RobinForm { Cell, Face, Surface };
        // P2-8 group A: `side_face_closure` (SideTForm) DELETED -- `cell` was
        // set only by an abort test; the lateral closure is always `follow`.

        bool enabled = false;
        Mode mode = Mode::PrescribedTemperature;
        RobinForm robin_form = RobinForm::Cell;
        // P2-8: THE SURFACE FACE LIST (`surface_patch.face_list`, default 0).
        // With it on, every solid face touching void (all six directions) is
        // listed once per step and heated, weighted and judged from that one
        // list (BuildFaceList / FaceListHeat / FaceWeight). 0 => every site
        // runs its pre-P2-8 code.
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
            // patch ON they stay required exactly as before; with it OFF the mask
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

            // S1: Robin discretisation for convective_flame (microstructure
            // path only; see SurfaceCellFlux). `cell` (default) applies
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
                // P2-28: pinned_idle_cycles governed the pin rule, retired with robin_form = pinned.
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

            // T_ext / h_conv / T_flame. Under convective_flame, h_conv is the
            // flame heat-transfer coefficient h_fl and T_flame the flame
            // temperature: q = h_fl*(T_flame - T_surface). Under prescribed_T
            // they are reserved/unused.

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
            // simulated time between rows. 0 = off, so no existing leg gains a
            // file. P2-26 (census B1): parsed under EVERY jet_closure -- it sat
            // inside the enthalpy branch, so a jet_closure = none deck (CAL)
            // could not set it although its consumer is closure-independent.
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
                    // P2-29: surface_patch.foot_rule is retired at `pads` (the table in Parse)
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
    bool microstructure_enabled = false;
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

    // Step 6c: the grain-boundary contact conductance (MICRO; 0 = off).
    Set::Scalar h_gb_conductance   = 0.0;

    // Step 10: surface advancement and spall detachment. Default-off.
    // `phi = z_surface(x,y) - z_cell` is a level-set state; cells with
    // phi >= 0 are material, phi < 0 are void. `removed` is the binary
    // void mask. `spall_event`, `spall_thickness`, and `RoP_spall` are
    // per-step diagnostics reset at the start of each spall update.
    bool spall_enabled = false;

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

    // D2q-2b item 0: debug-only seeded interior void, `ilo jlo klo ihi jhi khi`
    // in index space. Empty = off (the default), and ApplySeedVoid returns
    // immediately. See the parse block for why it carries no diagnostic.
    std::vector<int> spall_seed_void;

    // Step 16: Sp-vs-depth h_spall scan parameters (sp_weibull branch only).
    // `spall.h_spall_max` caps the downward column scan in metres (safety
    // stop; in practice the K_I < K_Ic crossing always fires first). The
    // segment count `spall.h_spall_n_segments` falls back to `sp_n_segments`
    // when set to 0 (default), so the depth scan uses the same Simpson
    // resolution as `site A (deleted, P2-28)`'s surface Sp.
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
    // Distinct from the scalar `regime_mf` (registered earlier; kept).
    Set::Field<Set::Scalar> regime_field_mf;

    bool removal_enabled() const { return spall_enabled || vapor_enabled; }
    // ⚠ P2-28: THE FACE LIST'S EVALUATE-ONLY MODE -- the list on with spall and vapour off (the onset decks: Kant, the
    // pressure sweep, Hu). The state the list reads (phi, removed, regime, RoP) exists and the criterion runs on every
    // listed line; nothing is removed. Until P2-28 the list refused removal off (Parse), so list_state() equals
    // removal_enabled() on every deck that ran before it.
    bool list_state() const { return removal_enabled() || surface_patch.face_list; }
    bool evaluate_only() const { return surface_patch.face_list && !removal_enabled(); }

    // Step 12: spallation-model toggle. `none` (default) loads no criterion;
    // `sp_weibull` activates the LEFM Sp = K_I/K_Ic(T) onset criterion.
    // (`damage_law` was retired at P2-27 and aborts by value.)
    enum class SpallationModel { None, SpWeibull };
    SpallationModel spallation_model = SpallationModel::None;

    // Sp criterion parameters (sp_weibull branch only). Deterministic
    // surface flaw length a0 (m); Weibull randomization arrives in Step 15.
    Set::Scalar sp_a0 = 0.0;
    int         sp_n_segments = 64;     // composite Simpson segments

    // Lithostatic confining pressure (Pa, Kant sign convention: positive =
    // compressive). Applied as +p·nu/(1-nu) offset in both sigma_xx lambdas
    // (site A (deleted, P2-28) and the k_i_at_zcrack depth scan), NOT as an
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
    // under robin_form = surface). P2-28: surf_cols_robin / surf_q_mean, the
    // column reconstruction's witnesses, are deleted.
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
    // (spall.removal_events_csv; sp_weibull + spall.enabled, level 0). One
    // row per (step, column) whose surface moved: time, col_i, col_j, regime,
    // k_top, T_top, a_f, Sp_top, h_scan, h_applied, n_voided. Unlike the
    // column path's retired edge-triggered CSVs, it is a complete removal record.
    bool sp_emit_removal_events_csv = false;
    bool sp_removal_events_csv_header_written = false;

    // ------------------------------------------------------------------
    // P0: a real surface temperature for the firing criterion.
    //
    // ⚠ THE DEFECT THIS REPLACES. Under `spall.stress_source =
    // local_thermoelastic` (the Meier and Kant configuration) BOTH criterion
    // sites read the TOP CELL AVERAGE as the surface temperature:
    //   * site A (deleted, P2-28) -- sigma_xx is flat at T_a(i,j,k) over the
    //     whole cell (a0 = 4-20 um << dz, so every K_I sample lands in it) and
    //     K_Ic is evaluated at the same value;
    //   * the column path (deleted, P2-28)'s depth scan (T_at_zv) -- its surface lerp has
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
    // P2-28: the three P0 flags (sp_surface_T_recon / _apply, sp_removal_energy_flake) are deleted:
    // the face list always reconstructs (R1) and always removes a flake (R3).


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

    // Everything from a given x. Split out of the reconstruction so the
    // per-column builder can substitute a smoothed or held x (deliverables
    // 1.3, 5.1, 5.2) without duplicating a line of the algebra.
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
        // is a MESH-INDEPENDENT multiple of the surface rise (P6). With P0's
        // min(ell_eff, dz) the multiple was mesh-dependent and collapsed to the
        // cell average wherever ell_eff exceeded dz.
        r.h_f   = (x > 0.0) ? c_f * r.ell_eff : dz;
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

    // Upper bound on the gas temperature this configuration can present to the
    // rock, used only as a ceiling on the reconstruction. Time-independent and
    // read from Parse'd values, so it introduces NO dependence on
    // BuildFlameColumns -- which runs AFTER the criterion (see the call-order
    // note in the column path (deleted, P2-28)). Returns 0 when no gas temperature is defined, which
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
        // P2-27: the default `none` is what the 48 decks that never named the key
        // always had -- no criterion loaded (no Sp_field, sp.a0 not required).
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
        // runs; the offset is only applied in the sp_weibull sigma_xx
        // lambdas (site A (deleted, P2-28) + Step 16d's k_i_at_zcrack).
        pp.query_default("confining.p", value.confining_p, 0.0);
        if (!std::isfinite(value.confining_p))
            Util::Abort(INFO, "confining.p must be finite, got ",
                        value.confining_p);

        // P2-15: the buckling span, read on EVERY deck, BEFORE the sp_weibull early return:
        // ALAMO does not abort on an unused top-level key (standing rule 1), so a key read
        // only on one branch would be silently inert on the other. Its cross-key requirements
        // are checked in Parse after the face list and spall.surface_T are both parsed.
        pp.query_default("spall.buckle_span", value.sp_buckle_span, 0.0);
        if (!(std::isfinite(value.sp_buckle_span) && value.sp_buckle_span >= 0.0))
            Util::Abort(INFO, "spall.buckle_span must be finite and >= 0 m, got: ",
                        value.sp_buckle_span);
        // P2-23: the blocked-line record (default off; a LOG, read by nothing).
        // Checked against the span in Parse, beside the span's own face-list check.
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

        // Step 15b: Weibull per-cell flaw-length distribution. Default-off.
        // Only parsed under sp_weibull (we already returned for none);
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
        // exactly what P1f's lateral-drive fixture did, while Meier's
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

            pp.query_default("weibull.seed", value.weibull_seed, 0);

            // P2-28: weibull.detach_mode is retired at per_face, weibull.A_crit aborts (Parse)
        }
    }
};
} // namespace Integrator
#endif
