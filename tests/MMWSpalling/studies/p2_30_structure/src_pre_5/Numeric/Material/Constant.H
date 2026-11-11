#ifndef NUMERIC_MATERIAL_CONSTANT_H
#define NUMERIC_MATERIAL_CONSTANT_H

// Constant-property material: rho, c_p, kappa do not depend on T. Used by
// the microstructure path to give every mineral phase the same H<->T
// machinery as the global Zhang material (Step 6b). Phase-change scalars
// (T_m, fus_width, L_m, T_vap_lo, T_vap_hi, L_v) are inherited from the
// Material base; defaults disable phase change so it behaves like a pure
// linear H = rho*Cp*(T - T_ref) law.

#include "IO/ParmParse.H"
#include "Numeric/Material/Material.H"

namespace Numeric { namespace Material {

class Constant : public Material
{
public:
    static constexpr const char *name = "constant";

    Constant() = default;

    Constant(Set::Scalar rho_, Set::Scalar cp_, Set::Scalar kappa_)
        : m_rho(rho_), m_cp(cp_), m_kappa(kappa_) {}

    Constant(IO::ParmParse &pp, std::string a_prefix)
    {
        pp.queryclass<Constant>(a_prefix, *this);
    }

    static void Parse(Constant &v, IO::ParmParse &pp)
    {
        pp.query_required("rho",   v.m_rho);
        pp.query_required("cp",    v.m_cp);
        pp.query_required("kappa", v.m_kappa);
        // Phase-change scalars on the Material base; defaults disable phase
        // change so a constant-property rock acts as pure linear
        // H = rho*cp*(T - T_ref).
        pp.query_default("T_ref",     v.T_ref,     273.15);
        pp.query_default("T_m",       v.T_m,       1.0e30);
        pp.query_default("fus_width", v.fus_width, 0.0);
        pp.query_default("L_m",       v.L_m,       0.0);
        pp.query_default("T_vap_lo",  v.T_vap_lo,  1.0e30);
        pp.query_default("T_vap_hi",  v.T_vap_hi,  1.0e30);
        pp.query_default("L_v",       v.L_v,       0.0);
    }

    Set::Scalar rho  (Set::Scalar) const override { return m_rho;   }
    Set::Scalar cp   (Set::Scalar) const override { return m_cp;    }
    Set::Scalar kappa(Set::Scalar) const override { return m_kappa; }

    Set::Scalar m_rho   = 0.0;
    Set::Scalar m_cp    = 0.0;
    Set::Scalar m_kappa = 0.0;
};

}} // namespace Numeric::Material

#endif
