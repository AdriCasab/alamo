#ifndef NUMERIC_MMWBEAM_H
#define NUMERIC_MMWBEAM_H

// Beer-Lambert + Gaussian millimetre-wave beam source (main.tex eqs. 1-3).
//
//   Q̇(x,T) = P(x) α(T) exp(-α(T)(z_0 - z))                   (eq. 1)
//   P_GB(x,y,z) = (2 P_0)/(π ω_0²) (ω_0/ω(z))² exp(-2 r²/ω(z)²)   (eq. 2)
//   ω(z) = ω_0 sqrt(1 + (λ (z_0 - z)/(π ω_0²))²)             (eq. 3)
//
// α(T) is parsed from an AMReX expression in the single variable T (so
// constants like "10" or "10*(1 + 0.001*(T-300))" both work). Surface
// position z_0 is constant in step 3; surface advancement is step 11.
//
// See in-main-tex-you-will-quizzical-treasure.md, step 3.

#include <cmath>
#include <string>
#include <vector>

#include "AMReX_Parser.H"
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
        // Surface reflectivity R ∈ [0, 1). The volumetric Beer-Lambert source
        // is scaled by (1 − R) so the absorbed power on a column of length L
        // integrates to (1 − R) · P_0 · (1 − e^{-αL}). Default 0 keeps the
        // earlier "incident == absorbed" convention used in steps 1–4.
        pp.query_default("R",      v.R,      0.0);

        // Attenuation α(T) in m^-1 — parsed as an expression of T so a
        // temperature-dependent loss tangent can be plugged in without code
        // changes. Default 10 m^-1.
        std::string alpha_expr = "10";
        pp.query_default("alpha_expr", alpha_expr, "10");
        v.alpha_parser = amrex::Parser(alpha_expr);
        v.alpha_parser.registerVariables({"T"});
        v.alpha_f = v.alpha_parser.compile<1>();

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

    // Beer-Lambert volumetric source [W/m³]. Returns 0 above the surface.
    AMREX_GPU_HOST_DEVICE
    Set::Scalar Q(Set::Scalar x, Set::Scalar y, Set::Scalar z, Set::Scalar T) const
    {
        Set::Scalar dz = z0 - z;
        if (dz < 0.0) return 0.0;
        Set::Scalar zR = M_PI * omega0 * omega0 / lambda;
        Set::Scalar w  = omega0 * std::sqrt(1.0 + (dz / zR) * (dz / zR));
        Set::Scalar r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        Set::Scalar P_int = (2.0 * P0 / (M_PI * omega0 * omega0))
                          * (omega0 / w) * (omega0 / w)
                          * std::exp(-2.0 * r2 / (w * w));
        Set::Scalar a = alpha_f(T);
        return (1.0 - R) * P_int * a * std::exp(-a * dz);
    }

    bool active() const { return parsed && P0 != 0.0; }

    Set::Scalar P0     = 0.0;
    Set::Scalar omega0 = 1.0;
    Set::Scalar lambda = 1.0e-3;
    Set::Scalar x0     = 0.0;
    Set::Scalar y0     = 0.0;
    Set::Scalar z0     = 0.0;
    Set::Scalar R      = 0.0;   // surface reflectivity (0 = full absorption)

    amrex::Parser alpha_parser;
    amrex::ParserExecutor<1> alpha_f;

    std::vector<Set::Scalar> sched_t; // interval end times [s]
    std::vector<Set::Scalar> sched_P; // power per interval  [W]

private:
    bool parsed = false;
};
} // namespace Numeric

#endif
