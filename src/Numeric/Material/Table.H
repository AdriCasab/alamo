#ifndef NUMERIC_MATERIAL_TABLE_H
#define NUMERIC_MATERIAL_TABLE_H

// Tabulated H<->T inverter for a single Material. Owns the integration of
//   H(T) = integral_{T_ref}^T rho(T') * cp_eff(T') dT'
// where cp_eff bakes in the smoothed latent-heat boxcars from the Material
// base (fusion / vaporisation widths and specific latent heats). The table
// is built once per material instance and used by:
//
//   - the global `material` path (single Zhang/whatever instance), and
//   - the microstructure path (one Table per phase, indexed by phase ID).
//
// Phase fractions Lambda_S/L/V are derived from H, not T, so a sharp
// fusion plateau (fus_width = 0) where T jumps discontinuously across the
// latent interval still resolves the partial-melt fraction correctly.

#include <vector>
#include "Set/Set.H"
#include "Numeric/Material/Material.H"

namespace Numeric { namespace Material {

class Table
{
public:
    Table() = default;

    // Build the H(T) table for `mat` over [T_min, T_max] with step dT. The
    // table is shifted so H(mat.T_ref) == 0.
    void Build(const Material &mat,
               Set::Scalar T_min = 200.0,
               Set::Scalar T_max = 5000.0,
               Set::Scalar dT    = 0.5)
    {
        m_T_min = T_min;
        m_dT    = dT;
        m_N     = static_cast<int>((T_max - T_min) / dT) + 1;
        m_T.assign(m_N, 0.0);
        m_H.assign(m_N, 0.0);

        const Set::Scalar fus_lo = mat.T_m - 0.5 * mat.fus_width;
        const Set::Scalar fus_hi = mat.T_m + 0.5 * mat.fus_width;

        m_T[0] = T_min;
        m_H[0] = 0.0;
        for (int i = 1; i < m_N; ++i)
        {
            Set::Scalar T_lo = T_min + (i - 1) * dT;
            Set::Scalar T_hi = T_min + i * dT;
            m_T[i] = T_hi;
            // Trapezoidal integration of dH/dT = rho * cp_eff (cp_eff
            // includes the latent-heat boxcars via Material::effective_rhocp).
            Set::Scalar avg = 0.5 * (mat.effective_rhocp(T_lo)
                                    + mat.effective_rhocp(T_hi));
            m_H[i] = m_H[i - 1] + avg * dT;
        }
        Set::Scalar H_ref = H_of_T(mat.T_ref);
        for (auto &h : m_H) h -= H_ref;

        // Cache plateau breakpoints in H so phase fractions are inferable
        // from H alone (works for sharp plateaux too).
        m_fus_lo_T = fus_lo;
        m_fus_hi_T = fus_hi;
        m_vap_lo_T = mat.T_vap_lo;
        m_vap_hi_T = mat.T_vap_hi;
        m_T_m      = mat.T_m;
        m_T_v      = mat.T_vap_lo; // sharp-vap reference T (unused if width>0)

        m_H_fus_lo = H_of_T(fus_lo);
        m_H_fus_hi = H_of_T(fus_hi);
        // For sharp fusion (fus_width == 0) the integrator does not put a
        // boxcar pulse anywhere, so add the latent-heat jump explicitly.
        if (mat.fus_width <= 0.0 && mat.L_m > 0.0)
            m_H_fus_hi = m_H_fus_lo + mat.rho(mat.T_m) * mat.L_m;
        m_sharp_fus = (mat.fus_width <= 0.0) && (mat.L_m > 0.0);

        m_H_vap_lo = H_of_T(mat.T_vap_lo);
        m_H_vap_hi = H_of_T(mat.T_vap_hi);
        if (mat.T_vap_hi <= mat.T_vap_lo && mat.L_v > 0.0)
            m_H_vap_hi = m_H_vap_lo + mat.rho(mat.T_vap_lo) * mat.L_v;
        m_sharp_vap = (mat.T_vap_hi <= mat.T_vap_lo) && (mat.L_v > 0.0);
    }

    bool built() const { return m_N > 0; }

    Set::Scalar H_of_T(Set::Scalar T) const
    {
        if (T <= m_T_min)            return m_H.front();
        if (T >= m_T[m_N - 1])       return m_H.back();
        Set::Scalar f = (T - m_T_min) / m_dT;
        int i = static_cast<int>(f);
        Set::Scalar w = f - i;
        return m_H[i] * (1.0 - w) + m_H[i + 1] * w;
    }

    Set::Scalar T_of_H(Set::Scalar H) const
    {
        // Pin to plateau temperatures for sharp jumps so the inverse is
        // well-defined inside the latent interval.
        if (m_sharp_fus && H >= m_H_fus_lo && H <= m_H_fus_hi) return m_T_m;
        if (m_sharp_vap && H >= m_H_vap_lo && H <= m_H_vap_hi) return m_vap_lo_T;
        if (H <= m_H.front()) return m_T.front();
        if (H >= m_H.back())  return m_T.back();
        int lo = 0, hi = m_N - 1;
        while (hi - lo > 1)
        {
            int mid = (lo + hi) / 2;
            if (m_H[mid] <= H) lo = mid; else hi = mid;
        }
        Set::Scalar denom = m_H[hi] - m_H[lo];
        Set::Scalar w = (denom > 0.0) ? (H - m_H[lo]) / denom : 0.0;
        return m_T[lo] * (1.0 - w) + m_T[hi] * w;
    }

    // Combined inversion: returns T plus the phase fractions inferred from
    // H. Use this everywhere in lieu of T->Lambda inference, because phase
    // fractions inside a latent interval are not recoverable from T alone
    // when the plateau is sharp.
    void Invert(Set::Scalar H,
                Set::Scalar &T,
                Set::Scalar &LS, Set::Scalar &LL, Set::Scalar &LV) const
    {
        T = T_of_H(H);
        if (m_H_fus_hi > m_H_fus_lo && H >= m_H_fus_lo && H <= m_H_fus_hi)
        {
            LL = (H - m_H_fus_lo) / (m_H_fus_hi - m_H_fus_lo);
            LS = 1.0 - LL; LV = 0.0; return;
        }
        if (m_H_vap_hi > m_H_vap_lo && H >= m_H_vap_lo && H <= m_H_vap_hi)
        {
            LV = (H - m_H_vap_lo) / (m_H_vap_hi - m_H_vap_lo);
            LL = 1.0 - LV; LS = 0.0; return;
        }
        if (H < m_H_fus_lo)       { LS = 1.0; LL = 0.0; LV = 0.0; }
        else if (H < m_H_vap_lo)  { LS = 0.0; LL = 1.0; LV = 0.0; }
        else                      { LS = 0.0; LL = 0.0; LV = 1.0; }
    }

private:
    std::vector<Set::Scalar> m_T, m_H;
    Set::Scalar m_T_min = 0.0, m_dT = 0.0;
    int         m_N     = 0;

    Set::Scalar m_fus_lo_T = 0.0, m_fus_hi_T = 0.0;
    Set::Scalar m_vap_lo_T = 0.0, m_vap_hi_T = 0.0;
    Set::Scalar m_T_m = 0.0, m_T_v = 0.0;
    Set::Scalar m_H_fus_lo = 0.0, m_H_fus_hi = 0.0;
    Set::Scalar m_H_vap_lo = 0.0, m_H_vap_hi = 0.0;
    bool        m_sharp_fus = false, m_sharp_vap = false;
};

}} // namespace Numeric::Material

#endif
