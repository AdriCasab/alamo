#ifndef NUMERIC_COHESIVEZONE_H
#define NUMERIC_COHESIVEZONE_H

#include "Set/Set.H"

#include <AMReX.H>

#include <cmath>

// Bilinear normal-opening cohesive-zone law (Step 9). The envelope is
// parameterised by the elastic stiffness `K_n`, the opening at peak `delta_c`,
// and the opening at full failure `delta_max`. The peak traction is
// `f_t = K_n * delta_c`. The fracture energy under the envelope is
// `G_c = 0.5 * f_t * delta_max`. Calibration helpers convert between (K_n,
// G_c, f_t) and (delta_c, delta_max).
//
// Irreversibility uses the secant from the origin to the envelope point at
// the historical maximum opening: t = (1 - D(delta_hist)) * K_n * delta_now,
// with D = 1 - t_env(delta_hist) / (K_n * delta_hist).
namespace Numeric
{
namespace CohesiveZone
{

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar EnvelopeTraction(Set::Scalar delta,
                             Set::Scalar K_n,
                             Set::Scalar delta_c,
                             Set::Scalar delta_max)
{
    if (!(delta > 0.0)) return 0.0;
    if (delta < delta_c) return K_n * delta;
    if (delta >= delta_max) return 0.0;
    Set::Scalar f_t = K_n * delta_c;
    Set::Scalar denom = delta_max - delta_c;
    if (!(denom > 0.0)) return 0.0;
    return f_t * (delta_max - delta) / denom;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar EnvelopeDamage(Set::Scalar delta_hist,
                           Set::Scalar K_n,
                           Set::Scalar delta_c,
                           Set::Scalar delta_max)
{
    if (!(delta_hist > delta_c)) return 0.0;
    if (delta_hist >= delta_max) return 1.0;
    Set::Scalar t_env = EnvelopeTraction(delta_hist, K_n, delta_c, delta_max);
    Set::Scalar t_lin = K_n * delta_hist;
    if (!(t_lin > 0.0)) return 0.0;
    Set::Scalar D = 1.0 - t_env / t_lin;
    if (D < 0.0) D = 0.0;
    if (D > 1.0) D = 1.0;
    return D;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar SecantTraction(Set::Scalar delta_now,
                           Set::Scalar delta_hist,
                           Set::Scalar K_n,
                           Set::Scalar delta_c,
                           Set::Scalar delta_max)
{
    if (!(delta_now > 0.0)) return 0.0;
    Set::Scalar D = EnvelopeDamage(delta_hist, K_n, delta_c, delta_max);
    return (1.0 - D) * K_n * delta_now;
}

}
}

#endif
