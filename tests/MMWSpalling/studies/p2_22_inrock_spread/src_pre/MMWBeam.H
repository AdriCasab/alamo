#ifndef NUMERIC_MMWBEAM_H
#define NUMERIC_MMWBEAM_H

// Beer-Lambert + Gaussian millimetre-wave beam source (main.tex eqs. 1-3).
//
//   Q̇(x) = (1 - R) P(x) α exp(-α (z_0 - z))                   (eq. 1)
//   P_GB(x,y,z) = (2 P_0)/(π ω_0²) (ω_0/ω(z))² exp(-2 r²/ω(z)²)   (eq. 2)
//   ω(z) = ω_0 sqrt(1 + (λ (z_0 - z)/(π ω_0²))²)             (eq. 3)
//
// P2-20: ONE DIELECTRIC DESCRIPTION. The rock's relative permittivity ε′
// (`beam.eps_r`) and loss tangent tan δ (`beam.tan_delta`) are declared once
// (the rock file); the POWER absorption coefficient and the normal-incidence
// reflection are COMPUTED from them and the deck's wavelength:
//   α = 2π √ε′ tan δ / λ   [1/m]   (main.tex eq. alpha_def; its "low-loss"
//                                    π tan δ √ε′/λ is the FIELD coefficient)
//   R = ((√ε′ − 1)/(√ε′ + 1))²     (Fresnel, normal incidence)
// P2-21: an optional `beam.tan_delta_T` ratio table makes the loss depend on
// temperature (alpha_at); with it the deposition is the column optical-depth
// sum (MMWSpalling::ColumnOpticalDepth). Absent or one pair: flat, as P2-20.
// The typed `beam.alpha_expr` and
// `beam.R` are retired (MMWSpalling::Parse aborts naming them). Surface
// position z_0 is constant in step 3; surface advancement is step 11.
//
// See in-main-tex-you-will-quizzical-treasure.md, step 3.

#include <cmath>
#include <string>
#include <vector>

#include "Set/Set.H"
#include "Util/Util.H"
#include "Unit/Unit.H"
#include "IO/ParmParse.H"

namespace Numeric
{
class MMWBeam
{
public:
    static constexpr const char* name = "mmwbeam";

    MMWBeam() = default;

    static void Parse(MMWBeam& v, IO::ParmParse& pp)
    {
        // Beam power [W]; default zero so the integrator reduces to step 2 if
        // no beam block is supplied.
        pp.query_default("P0",     v.P0,     std::string("0.0"),   Unit::Power());
        // Beam waist at waveguide exit [m].
        pp.query_default("omega0", v.omega0, std::string("1.0"),   Unit::Length());
        // MMW wavelength [m].
        pp.query_default("lambda", v.lambda, std::string("1.0e-3"),Unit::Length());
        // Beam-axis location (x0, y0) and surface depth z0 [m].
        pp.query_default("x0",     v.x0,     std::string("0.0"),   Unit::Length());
        pp.query_default("y0",     v.y0,     std::string("0.0"),   Unit::Length());
        pp.query_default("z0",     v.z0,     std::string("0.0"),   Unit::Length());
        // P2-20: the dielectric description, read AFTER `lambda` (every input
        // of the two derived members below is parsed before they are formed).
        // Dimensionless, through the same strtod path as `lambda`. Defaults:
        // vacuum (eps_r 1, tan_delta 0) => alpha 0, R 0.
        pp.query_default("eps_r",     v.eps_r,     std::string("1.0"), Unit::Less());
        pp.query_default("tan_delta", v.tan_delta, std::string("0.0"), Unit::Less());
        if (!(v.eps_r >= 1.0))
            Util::Abort(INFO, "beam.eps_r (relative permittivity) must be >= 1, got ", v.eps_r);
        if (!(v.tan_delta >= 0.0))
            Util::Abort(INFO, "beam.tan_delta (loss tangent) must be >= 0, got ", v.tan_delta);
        // The POWER coefficient (exp(-alpha s) on the power) and the Fresnel
        // reflection. The volumetric source is scaled by (1 − R), so the power
        // absorbed on a column of length L integrates to (1 − R) P_0 (1 − e^{-αL}).
        // ⚠ Evaluated as 2 pi sqrt(eps_r) * (tan_delta / lambda), the order kept by the
        // parentheses (no reassociation without fast-math): algebraically the formula
        // above, and the one order in which every beam-ON fixture's tan_delta string
        // (eps_r 1) reproduces its old typed alpha BIT-FOR-BIT (studies/
        // p2_20_beam_deposition/Z0.txt), so the frozen gate decks do not move.
        v.alpha = 2.0 * M_PI * std::sqrt(v.eps_r) * (v.tan_delta / v.lambda);
        {
            const Set::Scalar n = std::sqrt(v.eps_r);
            v.R = ((n - 1.0) / (n + 1.0)) * ((n - 1.0) / (n + 1.0));
        }
        // P2-21: the TEMPERATURE DEPENDENCE of the loss, as a dimensionless RATIO
        // on the declared cold value: pairs `T r` (T in K, strictly increasing;
        // r = tan_delta(T) / tan_delta > 0), linear between pairs, clamped at
        // both ends. alpha_at(T) = alpha * r(T); alpha itself and its parse
        // order are untouched. Absent or ONE pair => alpha_varies() is false and
        // every consumer keeps the closed form (byte identity on flat decks).
        // The deposition with a varying coefficient is the column optical-depth
        // sum (MMWSpalling::ColumnOpticalDepth), never alpha_at in the closed form.
        {
            std::vector<Set::Scalar> tab;
            pp.queryarr("tan_delta_T", tab);
            if (tab.size() % 2 != 0)
                Util::Abort(INFO, "beam.tan_delta_T must be PAIRS `T r` (K, ratio); got ", tab.size(), " numbers");
            v.tdT_T.clear();
            v.tdT_r.clear();
            for (std::size_t i = 0; i + 1 < tab.size(); i += 2)
            {
                if (!(tab[i + 1] > 0.0))
                    Util::Abort(INFO, "beam.tan_delta_T: every ratio must be > 0, got ", tab[i + 1], " at T = ", tab[i]);
                if (!v.tdT_T.empty() && !(tab[i] > v.tdT_T.back()))
                    Util::Abort(INFO, "beam.tan_delta_T: T must be strictly increasing, got ", tab[i],
                                " after ", v.tdT_T.back());
                v.tdT_T.push_back(tab[i]);
                v.tdT_r.push_back(tab[i + 1]);
            }
        }

        // Radial intensity profile. `gaussian` (default) is the focused-beam
        // P_GB(x,y,z) above (eqs. 2-3). `uniform` is a flat-top disk of radius
        // `beam.radius`: I = P_0/(π radius²) for r ≤ radius, 0 outside, with no
        // Gaussian waist evolution. The uniform profile models an impinging
        // flame jet (flat stagnation-zone heat flux) far better than a laser
        // Gaussian -- it avoids the central hot-spot that overheats the surface.
        std::string prof = "gaussian";
        pp.query_default("profile", prof, "gaussian");
        if      (prof == "gaussian")                    v.profile = Gaussian;
        else if (prof == "uniform" || prof == "tophat") v.profile = Uniform;
        else Util::Abort(INFO, "beam.profile must be 'gaussian' or 'uniform', got ", prof);
        // Disk radius [m]; required and used only for the uniform profile.
        pp.query_default("radius", v.radius, std::string("0.0"), Unit::Length());
        if (v.profile == Uniform && v.radius <= 0.0)
            Util::Abort(INFO, "beam.profile = uniform requires beam.radius > 0");

        // Optional piecewise-constant power schedule. If present, P0 above is
        // ignored and the active power is `schedule.P[i]` whenever
        // `t < schedule.t_end[i]`. Both arrays must have equal length and
        // `t_end` must be strictly increasing.
        pp.queryarr("schedule.t_end", v.sched_t);
        pp.queryarr("schedule.P",     v.sched_P);
        if (v.sched_t.size() != v.sched_P.size())
            Util::Abort(INFO, "MMWBeam.schedule.t_end and .P must have the same length");

        v.parsed = true;
    }

    // Active beam power at simulation time `t`. With no schedule supplied
    // this returns the constant `P0`.
    Set::Scalar power_at(Set::Scalar t) const
    {
        if (sched_t.empty()) return P0;
        for (std::size_t i = 0; i < sched_t.size(); ++i)
            if (t < sched_t[i]) return sched_P[i];
        return sched_P.back();
    }

    bool active() const { return parsed && P0 != 0.0; }

    // P2-21: the loss table varies (>= 2 pairs). Only then is the column
    // optical-depth sum built, and only then are its thermo columns registered.
    bool alpha_varies() const { return tdT_T.size() >= 2; }

    // P2-21: can this deck deposit at all? Known at parse time: with a schedule
    // the schedule's powers decide (P0 is ignored by power_at), else P0. A
    // beam-OFF deck that includes the rock file (and so the loss table) keeps
    // its thermo header.
    bool can_deposit() const
    {
        if (sched_t.empty()) return P0 > 0.0;
        for (Set::Scalar p : sched_P) if (p > 0.0) return true;
        return false;
    }

    // P2-21: the power absorption coefficient at temperature T [1/m]. With a
    // flat or absent table this is `alpha` itself.
    Set::Scalar alpha_at(Set::Scalar T) const
    {
        if (!alpha_varies()) return alpha;
        if (!(T > tdT_T.front())) return alpha * tdT_r.front();
        if (!(T < tdT_T.back()))  return alpha * tdT_r.back();
        std::size_t i = 1;
        while (tdT_T[i] < T) ++i;            // tdT_T[i-1] < T <= tdT_T[i]
        const Set::Scalar f = (T - tdT_T[i - 1]) / (tdT_T[i] - tdT_T[i - 1]);
        return alpha * (tdT_r[i - 1] + (tdT_r[i] - tdT_r[i - 1]) * f);
    }

    Set::Scalar P0     = 0.0;
    Set::Scalar omega0 = 1.0;
    Set::Scalar lambda = 1.0e-3;
    Set::Scalar x0     = 0.0;
    Set::Scalar y0     = 0.0;
    Set::Scalar z0     = 0.0;
    Set::Scalar eps_r     = 1.0;   // relative permittivity epsilon' (P2-20)
    Set::Scalar tan_delta = 0.0;   // loss tangent (P2-20)
    Set::Scalar alpha     = 0.0;   // POWER absorption coefficient [1/m], computed
    Set::Scalar R         = 0.0;   // normal-incidence reflection, computed
    std::vector<Set::Scalar> tdT_T;  // P2-21: loss-table temperatures [K]
    std::vector<Set::Scalar> tdT_r;  // P2-21: tan_delta(T) / tan_delta [-]

    enum ProfileType { Gaussian = 0, Uniform = 1 };
    int         profile = Gaussian; // radial intensity profile
    Set::Scalar radius  = 0.0;      // flat-top disk radius [m] (uniform only)

    std::vector<Set::Scalar> sched_t; // interval end times [s]
    std::vector<Set::Scalar> sched_P; // power per interval  [W]

private:
    bool parsed = false;
};
} // namespace Numeric

#endif
