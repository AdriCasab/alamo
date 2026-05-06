#ifndef NUMERIC_MATERIAL_ZHANG_H
#define NUMERIC_MATERIAL_ZHANG_H

// Zhang/Oglesby granite material model — the polynomial fits used in
// Zhang et al. (2023) Table 1 / Eq. 10. All coefficients are parsed from
// input so the user can substitute different rock types without touching
// code; sensible defaults reproduce Zhang's granite.
//
// theta = T - 273.15   (T in Kelvin)
//
// rho(T)   = a_rho * theta + b_rho                       (T <= T_rho_cap)
//          = rho_cap                                     (T >  T_rho_cap)
// kappa(T) = 1 / (a_k + b_k * theta)                     (T <= T_kappa_cap)
//          = kappa(T_kappa_cap)                          (T >  T_kappa_cap)
// cp(T)    = a_cp_lo*theta^2 + b_cp_lo*theta + c_cp_lo   (T <= T_cp_break)
//          = a_cp_hi*theta + b_cp_hi                     (T >  T_cp_break)

#include <string>
#include "IO/ParmParse.H"
#include "Numeric/Material/Material.H"

namespace Numeric { namespace Material {

class Zhang : public Material
{
public:
    static constexpr const char *name = "zhang";

    Zhang() = default;

    Zhang(IO::ParmParse &pp, std::string a_prefix)
    {
        pp.queryclass<Zhang>(a_prefix, *this);
    }

    static void Parse(Zhang &v, IO::ParmParse &pp)
    {
        // Density (Zhang Table 1, Eq. 10): rho = a_rho*theta + b_rho, capped.
        pp.query_default("rho.a",    v.a_rho,    -0.2877);
        pp.query_default("rho.b",    v.b_rho,    2884.0);
        pp.query_default("rho.cap",  v.rho_cap,  2300.0);
        pp.query_default("rho.T_cap",v.T_rho_cap,1273.0);

        // Conductivity (Zhang Table 1): kappa = 1/(a_k + b_k*theta), capped.
        pp.query_default("kappa.a",     v.a_kappa,    0.3154);
        pp.query_default("kappa.b",     v.b_kappa,    3.8e-4);
        pp.query_default("kappa.T_cap", v.T_kappa_cap,1273.0);

        // Specific heat (Zhang Table 1): piecewise about T_cp_break.
        pp.query_default("cp.a_lo",  v.a_cp_lo,  -2.73e-4);
        pp.query_default("cp.b_lo",  v.b_cp_lo,  0.894);
        pp.query_default("cp.c_lo",  v.c_cp_lo,  788.0);
        pp.query_default("cp.a_hi",  v.a_cp_hi,  0.184);
        pp.query_default("cp.b_hi",  v.b_cp_hi,  1249.0);
        pp.query_default("cp.T_break", v.T_cp_break, 1510.5);

        // Phase-change parameters (stored on the Material base).
        pp.query_default("T_ref",     v.T_ref,     273.15);
        pp.query_default("T_m",       v.T_m,       1510.5);
        pp.query_default("fus_width", v.fus_width, 200.0);
        pp.query_default("L_m",       v.L_m,       3.4e5);
        pp.query_default("T_vap_lo",  v.T_vap_lo,  3233.0);
        pp.query_default("T_vap_hi",  v.T_vap_hi,  3503.0);
        pp.query_default("L_v",       v.L_v,       4.8e6);
    }

    Set::Scalar rho(Set::Scalar T) const override
    {
        if (T >= T_rho_cap) return rho_cap;
        return a_rho * (T - 273.15) + b_rho;
    }

    Set::Scalar kappa(Set::Scalar T) const override
    {
        Set::Scalar Tc = (T < T_kappa_cap) ? T : T_kappa_cap;
        return 1.0 / (a_kappa + b_kappa * (Tc - 273.15));
    }

    Set::Scalar cp(Set::Scalar T) const override
    {
        Set::Scalar theta = T - 273.15;
        if (T <= T_cp_break) return a_cp_lo * theta * theta + b_cp_lo * theta + c_cp_lo;
        return a_cp_hi * theta + b_cp_hi;
    }

    // Density polynomial coefficients.
    Set::Scalar a_rho     = -0.2877;
    Set::Scalar b_rho     = 2884.0;
    Set::Scalar rho_cap   = 2300.0;
    Set::Scalar T_rho_cap = 1273.0;

    // Conductivity polynomial coefficients.
    Set::Scalar a_kappa     = 0.3154;
    Set::Scalar b_kappa     = 3.8e-4;
    Set::Scalar T_kappa_cap = 1273.0;

    // Specific-heat polynomial coefficients.
    Set::Scalar a_cp_lo   = -2.73e-4;
    Set::Scalar b_cp_lo   = 0.894;
    Set::Scalar c_cp_lo   = 788.0;
    Set::Scalar a_cp_hi   = 0.184;
    Set::Scalar b_cp_hi   = 1249.0;
    Set::Scalar T_cp_break = 1510.5;
};

}} // namespace Numeric::Material

#endif
