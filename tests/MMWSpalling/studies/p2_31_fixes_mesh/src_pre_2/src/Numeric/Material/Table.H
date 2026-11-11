#ifndef NUMERIC_MATERIAL_TABLE_H
#define NUMERIC_MATERIAL_TABLE_H

// Tabulated H<->T inverter for a single Material. Owns the integration of
//   H(T) = integral_{T_ref}^T rho(T') * cp_eff(T') dT'
// where cp_eff bakes in the smoothed latent-heat boxcars from the Material
// base (fusion / vaporisation widths and specific latent heats). The
// integrator builds one Table per mineral phase, once, from that phase's law
// (Constant, Grimsel or Zhang), indexed by phase ID. A SHARP latent band
// (zero width, L > 0) is a jump of rho*L in H at one temperature, carried by
// H_of_T / T_of_H rather than by the trapezoid (see Build).
//
// Phase fractions Lambda_S/L/V are derived from H, not T, so a sharp
// fusion plateau (fus_width = 0), where T stays pinned at T_m while H crosses
// the latent interval, still resolves the partial-melt fraction correctly.

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

        // ⚠⚠ P2-29: A SHARP latent band (zero width) is a JUMP in H at one temperature, which the trapezoid below
        // cannot carry: m_H stays the sensible enthalpy and the jump is added by H_of_T above the band and removed
        // by T_of_H above its top. Carried by the cached band top alone, a liquid cell would read L/c_p too
        // hot (100 K on the Stefan fixture) and ConvertTtoH would drop rho*L.
        // The jumps are set before ANY H_of_T call (H_ref, the band ends) so every cached value carries them.
        m_T_m      = mat.T_m;
        m_vap_lo_T = mat.T_vap_lo;
        m_sharp_fus = (mat.fus_width <= 0.0) && (mat.L_m > 0.0);
        m_sharp_vap = (mat.T_vap_hi <= mat.T_vap_lo) && (mat.L_v > 0.0);
        m_fus_jump  = m_sharp_fus ? mat.rho(mat.T_m) * mat.L_m : 0.0;
        m_vap_jump  = m_sharp_vap ? mat.rho(mat.T_vap_lo) * mat.L_v : 0.0;

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
        m_vap_hi_T = mat.T_vap_hi;
        m_T_v      = mat.T_vap_lo; // sharp-vap reference T (cached only: T_of_H pins to m_vap_lo_T, nothing reads this)

        m_H_fus_lo = H_of_T(fus_lo);
        m_H_fus_hi = H_of_T(fus_hi);
        // For sharp fusion (fus_width == 0) the trapezoid above does not put a
        // boxcar pulse anywhere, so the band's top is its bottom plus the jump
        // (H_of_T(T_m) itself excludes it: the jump applies strictly above T_m).
        if (m_sharp_fus)
            m_H_fus_hi = m_H_fus_lo + m_fus_jump;

        m_H_vap_lo = H_of_T(mat.T_vap_lo);
        m_H_vap_hi = H_of_T(mat.T_vap_hi);
        if (m_sharp_vap)
            m_H_vap_hi = m_H_vap_lo + m_vap_jump;

        // Perf (2026-09-18): uniform-in-H seed table for T_of_H. For a
        // query H the seed lo = m_inv_lo[j] satisfies m_H[lo] <= H, and a
        // short forward walk reaches the unique bracket m_H[lo] <= H <
        // m_H[lo+1] that the bisection in T_of_H's fallback branch returns,
        // so the interpolated T is bit-identical to it. Saves ~14
        // data-dependent bisection steps per cell per step.
        m_inv_N   = 4 * m_N;
        m_inv_H0  = m_H.front();
        m_inv_dH  = (m_H.back() - m_H.front()) / m_inv_N;
        m_inv_lo.assign(m_inv_N + 1, 0);
        for (int j = 0, lo = 0; j <= m_inv_N; ++j)
        {
            const Set::Scalar Hj = m_inv_H0 + j * m_inv_dH;
            while (lo + 1 < m_N - 1 && m_H[lo + 1] <= Hj) ++lo;
            m_inv_lo[j] = lo;
        }
    }

    bool built() const { return m_N > 0; }
    // P2-24: the cached latent-band enthalpies [J/m^3, per volume, H(T_ref) = 0] -- read by the face list's
    // vapour rule (a cell at H >= H_vap_hi() has received its whole vaporisation enthalpy) and printed once
    // per run as the table the run inverts (its enthalpy jump across the fusion band).
    Set::Scalar H_fus_lo() const { return m_H_fus_lo; }
    Set::Scalar H_fus_hi() const { return m_H_fus_hi; }
    Set::Scalar H_vap_lo() const { return m_H_vap_lo; }
    Set::Scalar H_vap_hi() const { return m_H_vap_hi; }

    Set::Scalar H_of_T(Set::Scalar T) const
    {
        Set::Scalar h;
        if (T <= m_T_min)            h = m_H.front();
        else if (T >= m_T[m_N - 1])  h = m_H.back();
        else
        {
            Set::Scalar f = (T - m_T_min) / m_dT;
            int i = static_cast<int>(f);
            Set::Scalar w = f - i;
            h = m_H[i] * (1.0 - w) + m_H[i + 1] * w;
        }
        // P2-29: a sharp band's jump, strictly above its temperature (no-op unless sharp)
        if (m_sharp_fus && T > m_T_m)      h += m_fus_jump;
        if (m_sharp_vap && T > m_vap_lo_T) h += m_vap_jump;
        return h;
    }

    Set::Scalar T_of_H(Set::Scalar H) const
    {
        // Pin to plateau temperatures for sharp jumps so the inverse is
        // well-defined inside the latent interval.
        if (m_sharp_fus && H >= m_H_fus_lo && H <= m_H_fus_hi) return m_T_m;
        if (m_sharp_vap && H >= m_H_vap_lo && H <= m_H_vap_hi) return m_vap_lo_T;
        // P2-29: above a sharp band's top the sensible table is read without the band's jump (no-op unless sharp)
        if (m_sharp_vap && H > m_H_vap_hi) H -= m_vap_jump;
        if (m_sharp_fus && H > m_H_fus_hi) H -= m_fus_jump;
        if (H <= m_H.front()) return m_T.front();
        if (H >= m_H.back())  return m_T.back();
        int lo = 0;
        if (m_inv_N > 0)
        {
            int j = static_cast<int>((H - m_inv_H0) / m_inv_dH);
            if (j < 0) j = 0; else if (j > m_inv_N) j = m_inv_N;
            lo = m_inv_lo[j];
            // Hj <= H, so m_H[lo] <= Hj <= H; walk to the last such index.
            while (lo + 1 < m_N - 1 && m_H[lo + 1] <= H) ++lo;
        }
        else
        {
            int hi = m_N - 1;
            while (hi - lo > 1)
            {
                int mid = (lo + hi) / 2;
                if (m_H[mid] <= H) lo = mid; else hi = mid;
            }
        }
        const int hi = lo + 1;
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
    // Uniform-in-H seed table for T_of_H (see Build()).
    std::vector<int> m_inv_lo;
    int              m_inv_N  = 0;
    Set::Scalar      m_inv_H0 = 0.0, m_inv_dH = 0.0;

    Set::Scalar m_fus_lo_T = 0.0, m_fus_hi_T = 0.0;
    Set::Scalar m_vap_lo_T = 0.0, m_vap_hi_T = 0.0;
    Set::Scalar m_T_m = 0.0, m_T_v = 0.0;
    Set::Scalar m_H_fus_lo = 0.0, m_H_fus_hi = 0.0;
    Set::Scalar m_H_vap_lo = 0.0, m_H_vap_hi = 0.0;
    bool        m_sharp_fus = false, m_sharp_vap = false;
    Set::Scalar m_fus_jump = 0.0, m_vap_jump = 0.0;   // P2-29: rho*L of a SHARP band [J/m^3], else 0
};

}} // namespace Numeric::Material

#endif
