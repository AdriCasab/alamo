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

    #include "Integrator/MMWSpalling/Parse.H"

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

    #include "Integrator/MMWSpalling/Model.H"

    #include "Integrator/MMWSpalling/Heat.H"

    #include "Integrator/MMWSpalling/Flame.H"

    #include "Integrator/MMWSpalling/FaceList.H"

    #include "Integrator/MMWSpalling/Ledger.H"

    #include "Integrator/MMWSpalling/BeamProps.H"

    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow) override
    {
        // Reuses HeatConduction's tagging (∇T) and Mechanics' tagging.
        HeatConduction::TagCellsForRefinement(lev, a_tags, time, ngrow);
        Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev, a_tags, time, ngrow);
        // P2-28: the AMR surface band (surface.refine_depth / refine_height) is retired with the column path
    }

private:
    #include "Integrator/MMWSpalling/ParseRemoval.H"

    // The removal code lives in two class-body fragments, textually included
    // here: Shared.H (the helpers and the step prologue) and FaceRemoval.H (the
    // face list, the ONE removal path since P2-28 deleted the column path).
    #include "Integrator/MMWSpalling/Shared.H"
    #include "Integrator/MMWSpalling/FaceRemoval.H"

    #include "Integrator/MMWSpalling/Fields.H"

    #include "Integrator/MMWSpalling/SurfacePatch.H"

    #include "Integrator/MMWSpalling/State.H"
};
} // namespace Integrator
#endif
