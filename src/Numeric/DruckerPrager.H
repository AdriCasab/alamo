#ifndef NUMERIC_DRUCKERPRAGER_H
#define NUMERIC_DRUCKERPRAGER_H

#include "Set/Set.H"

#include <AMReX.H>

#include <cmath>

namespace Numeric
{
namespace DruckerPrager
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar I1(const Set::Matrix& sigma)
{
    return sigma.trace();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Matrix Deviatoric(const Set::Matrix& sigma)
{
    return sigma - (I1(sigma) / Set::Scalar(AMREX_SPACEDIM))
                 * Set::Matrix::Identity();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar J2(const Set::Matrix& sigma)
{
    Set::Matrix s = Deviatoric(sigma);
    return 0.5 * (s.transpose() * s).trace();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar SqrtJ2(const Set::Matrix& sigma)
{
    Set::Scalar j2 = J2(sigma);
    return (j2 > 0.0) ? std::sqrt(j2) : 0.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar AlphaFromPhiRadians(Set::Scalar phi_rad)
{
    Set::Scalar s = std::sin(phi_rad);
    return (2.0 * s) / (std::sqrt(3.0) * (3.0 - s));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar KFromCohesionPhiRadians(Set::Scalar cohesion, Set::Scalar phi_rad)
{
    Set::Scalar s = std::sin(phi_rad);
    Set::Scalar c = std::cos(phi_rad);
    return (6.0 * cohesion * c) / (std::sqrt(3.0) * (3.0 - s));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar YieldFunction(const Set::Matrix& sigma,
                          Set::Scalar alpha_dp,
                          Set::Scalar k_dp)
{
    return SqrtJ2(sigma) + alpha_dp * I1(sigma) - k_dp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar Psi(const Set::Matrix& sigma,
                Set::Scalar alpha_dp,
                Set::Scalar k_dp)
{
    return YieldFunction(sigma, alpha_dp, k_dp) / k_dp;
}
}
}

#endif
