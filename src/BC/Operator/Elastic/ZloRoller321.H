//
// Hu-style zlo roller for 3D static thermoelastic validation runs.
//
// The bottom support is u_z = 0 everywhere on zlo, with homogeneous
// tangential tractions evaluated against the physical zlo normal even on
// bottom edges and corners. Two small bottom-corner anchors remove the
// remaining rigid-body x/y translations and rotation about z.
//

#ifndef BC_OPERATOR_ELASTIC_ZLOROLLER321_H
#define BC_OPERATOR_ELASTIC_ZLOROLLER321_H

#include "BC/Operator/Elastic/Elastic.H"

namespace BC
{
namespace Operator
{
namespace Elastic
{
class ZloRoller321 : public Elastic
{
public:
    static constexpr const char* name = "zlo_roller_321";

    ZloRoller321() = default;

    ZloRoller321(IO::ParmParse &pp, std::string name): ZloRoller321()
    { pp_queryclass(name, *this); }

    ~ZloRoller321() override = default;

    using Elastic::Init;

    void
    Init(amrex::FabArray<amrex::BaseFab<Set::Vector>> *a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
        (void)a_homogeneous;
        amrex::Box domain(a_geom.Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
        for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(2);
            bx = bx & domain;
            amrex::Array4<Set::Vector> const& rhs = a_rhs->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (isBoundary(i, j, k, lo, hi))
                    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                        rhs(i, j, k)(dir) = 0.0;
            });
        }
    }

    void
    Init(amrex::MultiFab *a_rhs,
        const amrex::Geometry &a_geom,
        bool a_homogeneous = false) const override
    {
        (void)a_homogeneous;
        amrex::Box domain(a_geom.Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
        for (amrex::MFIter mfi(*a_rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(2);
            bx = bx & domain;
            amrex::Array4<amrex::Real> const& rhs = a_rhs->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (isBoundary(i, j, k, lo, hi))
                    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                        rhs(i, j, k, dir) = 0.0;
            });
        }
    }

    std::array<Type, AMREX_SPACEDIM>
    getType(const int &i, const int &j, [[maybe_unused]] const int &k,
        const amrex::Box &domain) override
    {
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
        if (!isBoundary(i, j, k, lo, hi))
            return {AMREX_D_DECL(Type::None, Type::None, Type::None)};

#if AMREX_SPACEDIM == 3
        if (k == lo.z)
        {
            if (i == lo.x && j == lo.y)
                return {Type::Displacement, Type::Displacement, Type::Displacement};
            if (i == hi.x && j == hi.y)
                return {Type::Traction, Type::Displacement, Type::Displacement};
            return {Type::Traction, Type::Traction, Type::Displacement};
        }
        return {Type::Traction, Type::Traction, Type::Traction};
#else
        return {AMREX_D_DECL(Type::Traction, Type::Traction, Type::Traction)};
#endif
    }

    AMREX_FORCE_INLINE
    Set::Vector
    operator()(const Set::Vector &u,
        const Set::Matrix &gradu,
        const Set::Matrix &sigma,
        const int &i, const int &j, const int &k,
        const amrex::Box &domain) override
    {
        (void)gradu;
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
        if (!isBoundary(i, j, k, lo, hi)) return Set::Vector::Zero();

        std::array<Type, AMREX_SPACEDIM> bc_type = getType(i, j, k, domain);
        Set::Vector n = outwardNormal(i, j, k, lo, hi);
#if AMREX_SPACEDIM == 3
        if (k == lo.z) n = Set::Vector(0.0, 0.0, -1.0);
#endif
        return apply(bc_type, u, sigma, n);
    }

    static void
    Parse(ZloRoller321 &value, IO::ParmParse &pp)
    {
        (void)value;
        (void)pp;
#if AMREX_SPACEDIM != 3
        Util::Abort(INFO, "zlo_roller_321 is only defined for 3D elastic solves");
#endif
    }

private:
    AMREX_FORCE_INLINE
    static bool
    isBoundary(const int &i, const int &j, [[maybe_unused]] const int &k,
        const amrex::Dim3 &lo, const amrex::Dim3 &hi)
    {
        return AMREX_D_TERM(i == lo.x || i == hi.x,
            || j == lo.y || j == hi.y,
            || k == lo.z || k == hi.z);
    }

    AMREX_FORCE_INLINE
    static Set::Vector
    outwardNormal(const int &i, const int &j, [[maybe_unused]] const int &k,
        const amrex::Dim3 &lo, const amrex::Dim3 &hi)
    {
        Set::Vector n = Set::Vector::Zero();
        int count = 0;

        if (i == lo.x)      { n(0) = -1.0; count++; }
        else if (i == hi.x) { n(0) =  1.0; count++; }
#if AMREX_SPACEDIM >= 2
        if (j == lo.y)      { n(1) = -1.0; count++; }
        else if (j == hi.y) { n(1) =  1.0; count++; }
#endif
#if AMREX_SPACEDIM == 3
        if (k == lo.z)      { n(2) = -1.0; count++; }
        else if (k == hi.z) { n(2) =  1.0; count++; }
#endif

        if (count == 2)      n *= SQRT2INV;
        else if (count == 3) n *= SQRT3INV;
        return n;
    }

    AMREX_FORCE_INLINE
    static Set::Vector
    apply(const std::array<Type, AMREX_SPACEDIM> &bc_type,
        const Set::Vector &u,
        const Set::Matrix &sigma,
        const Set::Vector &n)
    {
        Set::Vector f = Set::Vector::Zero();
        Set::Vector traction = sigma * n;
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            if (bc_type[dir] == Type::Displacement)
                f(dir) = u(dir);
            else if (bc_type[dir] == Type::Traction)
                f(dir) = traction(dir);
        }
        return f;
    }
};
}
}
}

#endif
