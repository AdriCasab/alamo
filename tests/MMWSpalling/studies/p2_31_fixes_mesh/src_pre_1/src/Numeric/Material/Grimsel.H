#ifndef NUMERIC_MATERIAL_GRIMSEL_H
#define NUMERIC_MATERIAL_GRIMSEL_H

// Grimsel / Central Aare granite — the MEASURED property set for the rock the
// Meier and Kant experiments were actually run on. P2-2 changes 2 and 5.
//
// EVERY law and table here is sourced, by line, from
//   tests/MMWSpalling/validation/grimsel_properties/
//       grimsel_granite_properties_reference.md
// and the offline twin that verifies them against that file is
//   tests/MMWSpalling/studies/p2_2_props/grimsel.py
// The bracketed numbers below are that file's line numbers. Standing rule 8: a
// number is backed only if both the file it was read from and the code that
// reduced it are in the repository.
//
//   theta = T - 273.15   (T in Kelvin; the reference's fits are in degC)
//
//   c_p(theta)    = 1370 - 178000/(theta + 271)          J/(kg K)      [43]
//   H(theta)      = 1370(theta-20) - 178000 ln((theta+271)/291)        [48]
//   lambda(theta) = 2000/(theta + 563)                   W/(m K)       [89]
//     `lambda.cap_theta` (default -1 = off) holds lambda at its value at that
//     theta and above: the P2-24 melt heat-transport SENSITIVITY (Zhang's
//     kappa cap at 1273 K is the recorded alternative to the extrapolated fit).
//   rho           = 2650                                 kg/m^3       [141]
//
// ⚠ kappa is NEVER set independently of lambda/(rho c_p) [112]: this class
//   exposes lambda (the CONDUCTIVITY, which the integrator writes to its
//   conductivity field) and `kappa()` returns the DIFFUSIVITY formed from it,
//   exactly as the reference demands. Its three fits agree to 5.0 % against a 6.5 %
//   uncertainty propagated from c_p 768+/-48 and rho 2650+/-50.
//   ⚠ UNIT TRAP: Zhang::kappa and Constant::kappa are CONDUCTIVITIES [W/(m K)];
//   only this class's kappa() is a diffusivity [m^2/s].
//
// ⚠ H(theta) here (`H_of`, J/kg, zero at 20 degC) is ANALYTIC and is the exact
//   integral of c_p. The integrator does not call it: a run's H<->T inversion
//   is the phase's Table, a trapezoid of rho*c_p built from this class's
//   rho() and cp() (Table::Build).
//
// =========================================================================
// THE THERMAL STRAIN — P2-2 change 5, and the reason this class exists.
// =========================================================================
// ⚠⚠ THE STRESS IS A FORMULA, NOT A NUMBER. Three independent attempts at the
// one thermoelastic stress in this rock produced three different wrong answers,
// because the formula was never written down beside the result:
//
//       sigma(T_s) = E(T_s) * INTEGRAL[T_0..T_s] alpha(theta) dtheta / (1 - nu)
//
//   * `alpha(T_s) * dT` OVERSTATES it: alpha TRIPLES between 300 and 500 degC,
//     so its value at T_s is far above its mean over the range.
//   * the accumulated strain with COLD stiffness UNDERSTATES it: E halves by
//     600 degC.
//   * `INTEGRAL(E*alpha dtheta)` is wrong as well: ELASTICITY HAS NO PATH
//     MEMORY. Stress is the CURRENT stiffness times the CURRENT elastic strain.
//
// ⇒ accumulate the STRAIN, then multiply by the stiffness AT T_s. The strain is
// tabulated as a cumulative integral the same way the enthalpy table is, and
// `unit/lateral_flux_lock` asserts that `sigma_thermal` below has exactly ONE
// caller in the integrator (`LawSigma`). If anyone writes `alpha(T)*dT` at one site, the
// relief factor R silently absorbs a factor of two -- invisibly, because R is
// fitted. That is gate (b6).
//
// alpha, E/E_0 and K_Ic/K_Ic,0 come from the reference's Table 3 [229-238].
// ⚠⚠ alpha is the reference's ONE [ESTIMATE] (+/-30 % above 300 degC), and part
// of the differential expansion opens grain-boundary microcracks rather than
// loading the flake [172] -- so `alpha_cap_theta` holds it at a stated value,
// which is the (c7) sensitivity.

#include <cmath>
#include <string>
#include "IO/ParmParse.H"
#include "Numeric/Material/Material.H"

namespace Numeric { namespace Material {

class Grimsel : public Material
{
public:
    static constexpr const char *name = "grimsel";

    // Table 3 [229-238]: theta [degC], alpha [1e-6/K], E/E_0, K_Ic/K_Ic,0.
    // ⚠ TEN rows. The first offline pass TYPED this table and silently omitted
    // 550 and 650 degC; its verifier still "passed", because a check over the
    // rows you typed cannot find a row you left out. Keep all ten.
    static constexpr int  NT3 = 10;
    static constexpr Set::Scalar T3_th[NT3] =
        {  20.0, 100.0, 200.0, 300.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0 };
    static constexpr Set::Scalar T3_alpha[NT3] =
        {   8.4,   8.7,   9.1,   9.5,  16.4,  19.9,  23.4,  26.9,  28.3,  28.4 };
    static constexpr Set::Scalar T3_EoE0[NT3] =
        {  1.00,  1.00,  1.00,  0.93,  0.85,  0.75,  0.65,  0.58,  0.50,  0.50 };
    static constexpr Set::Scalar T3_KIc[NT3] =
        {  1.00,  1.00,  0.90,  0.82,  0.78,  0.77,  0.75,  0.65,  0.60,  0.60 };

    Grimsel() { BuildTables(); }

    Grimsel(IO::ParmParse &pp, std::string a_prefix)
    {
        pp.queryclass<Grimsel>(a_prefix, *this);
    }

    static void Parse(Grimsel &v, IO::ParmParse &pp)
    {
        // c_p [43]. c_p(25) = 768 +/- 48 J/(kg K).
        pp.query_default("cp.a",   v.cp_a,   1370.0);
        pp.query_default("cp.b",   v.cp_b,   178000.0);
        pp.query_default("cp.th0", v.cp_th0, 271.0);
        // lambda [89], first heating, intact rock.
        pp.query_default("lambda.a",   v.lam_a,   2000.0);
        pp.query_default("lambda.th0", v.lam_th0, 563.0);
        // hold lambda at its value at this theta and above. -1 = off (the fit
        // extrapolated through the melt; P2-24's production setting).
        pp.query_default("lambda.cap_theta", v.lam_cap_theta, -1.0);
        // rho [141].
        pp.query_default("rho",   v.rho_0,  2650.0);
        // E_0 and nu [180, 195]. ⚠ mu moves with E: mu = E/(2(1+nu)).
        pp.query_default("E0",    v.E_0,    53.3e9);
        pp.query_default("nu",    v.nu,     0.30);
        // K_Ic,0 [209], and the AT-TEMPERATURE ratio of Table 3 -- NOT the
        // post-cooling shape, which reads 32-51 % low at the spall temperature.
        pp.query_default("K_Ic0", v.K_Ic0,  1.5);
        // ⚠⚠ (c7): hold alpha at its value at this theta and above. -1 = off
        // (face value). The reference's [172] says to keep the bulk alpha at
        // the low side near T_s because part of the expansion opens
        // microcracks instead of loading the flake.
        pp.query_default("alpha_cap_theta", v.alpha_cap_theta, -1.0);
        v.BuildTables();
    }

    // ---- the Material interface -----------------------------------------
    Set::Scalar rho(Set::Scalar /*T*/) const override { return rho_0; }

    Set::Scalar cp(Set::Scalar T) const override
    {
        const Set::Scalar th = T - 273.15;
        return cp_a - cp_b / (th + cp_th0);
    }

    // ⚠ the DIFFUSIVITY, formed as lambda/(rho c_p) [112], never set apart.
    Set::Scalar kappa(Set::Scalar T) const override
    {
        return lambda(T) / (rho(T) * cp(T));
    }

    // ---- the laws the criterion and the stress need ----------------------
    Set::Scalar lambda(Set::Scalar T) const       // W/(m K) [89]
    {
        Set::Scalar th = T - 273.15;
        if (lam_cap_theta > 0.0 && th > lam_cap_theta) th = lam_cap_theta;
        return lam_a / (th + lam_th0);
    }

    // The analytic enthalpy [48], exactly the integral of cp above.
    Set::Scalar H_of(Set::Scalar T) const
    {
        const Set::Scalar th = T - 273.15;
        return cp_a * (th - 20.0)
             - cp_b * std::log((th + cp_th0) / (20.0 + cp_th0));
    }

    Set::Scalar alpha_of(Set::Scalar T) const      // 1/K, Table 3
    {
        Set::Scalar th = T - 273.15;
        if (alpha_cap_theta > 0.0 && th > alpha_cap_theta) th = alpha_cap_theta;
        return 1.0e-6 * Interp(th, T3_alpha);
    }

    // Pa, Table 3, read from the DENSE UNIFORM table (P2-2b deliverable 4).
    //
    // ⚠⚠ WHY THIS ONE NEEDS NO ANALYTIC REMAINDER AND `strain_of` DOES.
    // `strain_of` interpolates a CUMULATIVE INTEGRAL of a piecewise-linear
    // alpha, which is piecewise QUADRATIC, so linear interpolation between
    // nodes is wrong at second order (it read 4.5e-9 off a dense quadrature).
    // `E/E_0` is itself piecewise LINEAR in theta, and every one of Table 3's
    // kinks (20/100/200/300/400/450/500/550/600/650 degC) is an integer that
    // lands exactly on this 1 K grid -- so no cell straddles a kink and linear
    // interpolation reproduces `Interp` to round-off at EVERY theta, not just
    // at the nodes. ⇒ the dense table is a speed change, not an accuracy one,
    // which is what lets `(b2)` hold at 1e-12.
    Set::Scalar E_of(Set::Scalar T) const
    {
        return E_0 * DenseAt(T - 273.15, EoE0_tab);
    }

    Set::Scalar mu_of(Set::Scalar T) const         // Pa
    {
        return E_of(T) / (2.0 * (1.0 + nu));
    }

    Set::Scalar K_Ic_of(Set::Scalar T) const       // MPa m^0.5, AT TEMPERATURE
    {
        return K_Ic0 * Interp(T - 273.15, T3_KIc);
    }

    // =====================================================================
    // THE STRAIN TABLE — the single source for every stress site.
    // =====================================================================
    // Cumulative INTEGRAL[20 degC .. theta] alpha dtheta, trapezoidal on a
    // uniform 1 K grid from 20 to 1000 degC; between nodes `strain_of` adds
    // the analytic remainder. ⚠ 1 K spacing makes the table agree with a dense offline
    // quadrature to better than 1e-12 -- gate (b6) asks for 1e-10.
    // ⚠ Called from the constructor AND from the end of `Parse`:
    // `alpha_cap_theta` (which reshapes the strain table) is parsed
    // AFTER default construction, so a table built only once would silently
    // ignore it (`E0` needs no rebuild: the stiffness table holds the RATIO
    // E/E_0). `(c7)` is exactly the `alpha_cap_theta` sensitivity, so a
    // stale strain table would make that row measure nothing.
    void BuildTables()
    {
        // the stiffness ratio on the same uniform grid as the strain
        for (int i = 0; i < NSTR; ++i)
            EoE0_tab[i] = Interp(STR_TH0 + Set::Scalar(i) * STR_DTH, T3_EoE0);

        strain_tab[0] = 0.0;
        for (int i = 1; i < NSTR; ++i)
        {
            const Set::Scalar th0 = STR_TH0 + Set::Scalar(i - 1) * STR_DTH;
            const Set::Scalar th1 = STR_TH0 + Set::Scalar(i) * STR_DTH;
            const Set::Scalar a0 = AlphaAtTheta(th0), a1 = AlphaAtTheta(th1);
            strain_tab[i] = strain_tab[i - 1] + 0.5 * (a0 + a1) * STR_DTH;
        }
    }

    // INTEGRAL[20 degC .. T] alpha dtheta, dimensionless.
    //
    // ⚠⚠ THE REMAINDER IS ANALYTIC, NOT INTERPOLATED. alpha is piecewise
    // LINEAR in theta (Table 3), so the cumulative integral is piecewise
    // QUADRATIC -- interpolating the cumulative table linearly between nodes
    // is wrong at second order and read 4.5e-9 off a dense quadrature, which
    // fails gate (b6)'s 1e-10 bar. Table 3's kinks all land on integer theta
    // and the grid is 1 K, so no cell straddles a kink and the exact
    // contribution from the last node is one trapezoid:
    //        strain(theta) = tab[i] + (alpha(theta_i) + alpha(theta))/2 * d
    // which is exact to round-off at every theta, not just at the nodes.
    Set::Scalar strain_of(Set::Scalar T) const
    {
        const Set::Scalar th = T - 273.15;
        if (th <= STR_TH0) return 0.0;
        const Set::Scalar x = (th - STR_TH0) / STR_DTH;
        int i = static_cast<int>(x);
        if (i >= NSTR - 1) return strain_tab[NSTR - 1];
        if (i < 0) i = 0;
        const Set::Scalar th_i = STR_TH0 + Set::Scalar(i) * STR_DTH;
        const Set::Scalar d = th - th_i;
        return strain_tab[i]
             + 0.5 * (AlphaAtTheta(th_i) + AlphaAtTheta(th)) * d;
    }

    // ⚠⚠ THE STRESS. Every site that needs a thermoelastic stress in this
    // material calls THIS, so there is one formula and one place to be wrong.
    // `T_s` is the hot face, `T_0` the unheated reference.
    Set::Scalar sigma_thermal(Set::Scalar T_s, Set::Scalar T_0) const
    {
        const Set::Scalar eps = strain_of(T_s) - strain_of(T_0);
        return E_of(T_s) * eps / (1.0 - nu);
    }

    // Parsed constants, all defaulted to the reference's recommendations.
    Set::Scalar cp_a = 1370.0, cp_b = 178000.0, cp_th0 = 271.0;
    Set::Scalar lam_a = 2000.0, lam_th0 = 563.0, lam_cap_theta = -1.0;
    Set::Scalar rho_0 = 2650.0;
    Set::Scalar E_0 = 53.3e9, nu = 0.30;
    Set::Scalar K_Ic0 = 1.5;
    Set::Scalar alpha_cap_theta = -1.0;

private:
    static constexpr int NSTR = 981;          // 20 .. 1000 degC at 1 K
    static constexpr Set::Scalar STR_TH0 = 20.0;
    static constexpr Set::Scalar STR_DTH = 1.0;
    Set::Scalar strain_tab[NSTR] = {0.0};
    Set::Scalar EoE0_tab[NSTR]   = {0.0};

    // Linear read of a dense uniform table on the strain grid, clamped at both
    // ends exactly as `Interp` clamps (below 20 degC -> the first node; above
    // 1000 -> the last, and Table 3 is already flat above 650).
    static Set::Scalar DenseAt(Set::Scalar th, const Set::Scalar *tab)
    {
        if (th <= STR_TH0) return tab[0];
        const Set::Scalar x = (th - STR_TH0) / STR_DTH;
        int i = static_cast<int>(x);
        if (i >= NSTR - 1) return tab[NSTR - 1];
        if (i < 0) i = 0;
        return tab[i] + (x - Set::Scalar(i)) * (tab[i + 1] - tab[i]);
    }

    Set::Scalar AlphaAtTheta(Set::Scalar th) const
    {
        Set::Scalar t = th;
        if (alpha_cap_theta > 0.0 && t > alpha_cap_theta) t = alpha_cap_theta;
        return 1.0e-6 * Interp(t, T3_alpha);
    }

    // Table 3 lookup, clamped at both ends.
    static Set::Scalar Interp(Set::Scalar th, const Set::Scalar *y)
    {
        if (th <= T3_th[0])       return y[0];
        if (th >= T3_th[NT3 - 1]) return y[NT3 - 1];
        for (int i = 0; i < NT3 - 1; ++i)
            if (th <= T3_th[i + 1])
            {
                const Set::Scalar f = (th - T3_th[i]) / (T3_th[i + 1] - T3_th[i]);
                return (1.0 - f) * y[i] + f * y[i + 1];
            }
        return y[NT3 - 1];
    }
};

}} // namespace Numeric::Material

#endif
