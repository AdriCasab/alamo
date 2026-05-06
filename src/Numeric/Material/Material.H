#ifndef NUMERIC_MATERIAL_MATERIAL_H
#define NUMERIC_MATERIAL_MATERIAL_H

// Abstract base for temperature-dependent thermophysical materials. The
// ALAMO MMW spalling integrator dispatches κ(T), ρ(T), c_p(T) through this
// interface; phase-change scalars (T_m, T_v, latent heats, smoothing widths)
// are stored on the base so they are uniform across material implementations.
//
// Step 4b only needs Zhang/Oglesby granite. New rock types or per-mineral
// models slot in by deriving a new subclass and adding it to MMWSpalling's
// `select_default<...>("material", ...)` factory call.

#include "Set/Set.H"

namespace Numeric { namespace Material {

class Material
{
public:
    virtual ~Material() = default;

    virtual Set::Scalar rho  (Set::Scalar T) const = 0;   // [kg/m^3]
    virtual Set::Scalar cp   (Set::Scalar T) const = 0;   // [J/(kg K)]
    virtual Set::Scalar kappa(Set::Scalar T) const = 0;   // [W/(m K)]

    // ρ · c_p with the smoothed latent-heat boxcars baked in. This is what
    // implicit time-steppers want as their thermal-mass coefficient, since
    // it is the slope dH/dT used by the H↔T table. Equal to ρ·c_p outside
    // the fusion / vaporisation windows.
    Set::Scalar effective_rhocp(Set::Scalar T) const
    {
        Set::Scalar cpe = cp(T);
        if (fus_width > 0.0) {
            Set::Scalar fus_lo = T_m - 0.5 * fus_width;
            Set::Scalar fus_hi = T_m + 0.5 * fus_width;
            if (T >= fus_lo && T <= fus_hi) cpe += L_m / fus_width;
        }
        if (T_vap_hi > T_vap_lo) {
            if (T >= T_vap_lo && T <= T_vap_hi)
                cpe += L_v / (T_vap_hi - T_vap_lo);
        }
        return rho(T) * cpe;
    }

    // Phase-change scalars (defaults disable phase change).
    Set::Scalar T_ref     = 273.15;   // [K] reference T at which H = 0
    Set::Scalar T_m       = 1.0e30;   // [K] fusion centre temperature
    Set::Scalar fus_width = 0.0;      // [K] full width of fusion smoothing
    Set::Scalar L_m       = 0.0;      // [J/kg] latent heat of fusion
    Set::Scalar T_vap_lo  = 1.0e30;   // [K] start of vaporisation range
    Set::Scalar T_vap_hi  = 1.0e30;   // [K] end of vaporisation range
    Set::Scalar L_v       = 0.0;      // [J/kg] latent heat of vaporisation
};

}} // namespace Numeric::Material

#endif
