#ifndef NUMERIC_SPCRITERION_H
#define NUMERIC_SPCRITERION_H

#include "Set/Set.H"
#include "Util/Util.H"

#include <AMReX.H>

#include <cmath>

// Step 12: LEFM Sp = K_I / K_Ic(T) onset criterion. Edge-crack-in-half-space
// stress intensity via the Tada weight function (Tada-Paris-Irwin 1985 /
// Kant 2017 eq. 6):
//
//     f_K(xi) = (1.3 - 0.3 xi^(5/4)) / sqrt(1 - xi^2)
//     K_I    = 2 sqrt(a/pi) * integral_0^1 sigma_xx(xi*a) f_K(xi) dxi
//
// The 1/sqrt(1-xi^2) singularity at xi -> 1 is integrable. We remove it
// exactly by the substitution xi = sin(theta), dxi = cos(theta) dtheta,
// 1/sqrt(1-xi^2) dxi = dtheta:
//
//     integral_0^1 g(xi)/sqrt(1-xi^2) dxi = integral_0^{pi/2} g(sin(theta)) dtheta
//
// where g(xi) = 1.3 - 0.3 xi^(5/4). The result is a smooth integral in
// theta in [0, pi/2] handled by composite Simpson's rule.
//
// `K_Ic(T)` ships with two laws:
//   nasseri_table — piecewise linear on the Nasseri 2007 Westerly granite
//                   five-point table (clamped at endpoints) with an optional
//                   uniform multiplicative `scale` for cross-rock use (e.g.
//                   1.05 for Central Aare granite where K_Ic(RT) = 1.5 vs
//                   1.43 MPa*m^0.5).
//   linear        — K_Ic(T) = K_Ic0 * (1 - (T - T_ref)/T_melt), main.tex
//                   eq. 30 fallback for sensitivity sweeps.
namespace Numeric
{
namespace SpCriterion
{

// Tada weight g(xi) without the integrable singularity (the trig
// substitution above absorbs the 1/sqrt(1-xi^2) factor into the dtheta
// integration measure).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar TadaG(Set::Scalar xi)
{
    if (xi <= 0.0) return 1.3;
    Set::Scalar xi54 = std::pow(xi, 1.25);
    return 1.3 - 0.3 * xi54;
}

// Composite Simpson on theta in [0, pi/2] using xi = sin(theta). Requires
// an even number of segments; the caller passes n_segments.
//
// Template parameter SigmaXX is any callable
//   Set::Scalar sigma_xx(Set::Scalar depth)
// returning the in-plane tensile stress at the requested depth below the
// surface in metres.
template <typename SigmaXX>
AMREX_FORCE_INLINE
Set::Scalar K_I(SigmaXX sigma_at_depth, Set::Scalar a, int n_segments = 64)
{
    if (!(a > 0.0)) return 0.0;
    if (n_segments < 2) n_segments = 2;
    if (n_segments % 2 != 0) ++n_segments;

    constexpr Set::Scalar pi = 3.141592653589793238462643383279502884;
    const Set::Scalar h = (0.5 * pi) / Set::Scalar(n_segments);

    Set::Scalar acc = 0.0;
    for (int k = 0; k <= n_segments; ++k)
    {
        Set::Scalar theta = h * Set::Scalar(k);
        Set::Scalar xi    = std::sin(theta);
        Set::Scalar gxi   = TadaG(xi);
        Set::Scalar sxx   = sigma_at_depth(xi * a);
        Set::Scalar w;
        if      (k == 0 || k == n_segments) w = 1.0;
        else if (k % 2 == 1)                w = 4.0;
        else                                w = 2.0;
        acc += w * sxx * gxi;
    }
    acc *= h / 3.0;
    return 2.0 * std::sqrt(a / pi) * acc;
}

// Nasseri 2007 Table 4 (Westerly granite) K_Ic(T) tabulation. Temperatures
// are stored internally in Kelvin. Endpoint behaviour: clamped (constant
// extrapolation) below RT and above 850 C. The `scale` argument multiplies
// every table value uniformly for cross-rock use.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar KIcNasseri(Set::Scalar T, Set::Scalar scale = 1.0)
{
    // (T_K, K_Ic_MPa_sqrtm) anchor points (Nasseri 2007 + revised-model-approach.md s11).
    constexpr int N = 5;
    constexpr Set::Scalar T_K[N]  = {
        298.15,   // ~25 C (RT)
        523.15,   // 250 C
        723.15,   // 450 C
        923.15,   // 650 C
        1123.15   // 850 C
    };
    constexpr Set::Scalar Kic[N] = { 1.43, 1.35, 0.98, 0.43, 0.22 };

    // Clamp below the first point and above the last.
    if (T <= T_K[0])     return scale * Kic[0];
    if (T >= T_K[N - 1]) return scale * Kic[N - 1];
    for (int i = 0; i < N - 1; ++i)
    {
        if (T <= T_K[i + 1])
        {
            Set::Scalar f = (T - T_K[i]) / (T_K[i + 1] - T_K[i]);
            return scale * ((1.0 - f) * Kic[i] + f * Kic[i + 1]);
        }
    }
    return scale * Kic[N - 1];
}

// Linear fallback (main.tex eq. 30). Returns K_Ic0 below T_ref and clamps
// at 0 above T_ref + T_melt to avoid negative toughness.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar KIcLinear(Set::Scalar T,
                      Set::Scalar K_Ic0,
                      Set::Scalar T_ref,
                      Set::Scalar T_melt)
{
    if (!(T_melt > 0.0)) return K_Ic0;
    Set::Scalar v = K_Ic0 * (1.0 - (T - T_ref) / T_melt);
    if (v < 0.0) v = 0.0;
    return v;
}

// SI MPa*m^0.5 -> Pa*m^0.5 helper (table is in MPa*m^0.5; the integral is
// in Pa). Apply once when comparing K_I (Pa*m^0.5) to K_Ic (table units).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar MPa_sqrtm_to_Pa_sqrtm(Set::Scalar v) { return v * 1.0e6; }

} // namespace SpCriterion
} // namespace Numeric

#endif
