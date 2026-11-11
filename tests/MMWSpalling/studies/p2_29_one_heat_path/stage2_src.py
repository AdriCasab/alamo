#!/usr/bin/env python3
"""P2-29 stage 2 -- THE SEPARATE HEAT KERNELS DELETED (Goals 6-9), as one transform of the stage-2 pre-image
(S/src_pre_2/src/Integrator/MMWSpalling.H -> src/Integrator/MMWSpalling.H). Every anchor must occur exactly once; each
deleted span's first and last lines are printed for review.
  6  AdvanceConstAlpha / AdvanceMaterial / AdvanceMaterialImplicit, the dispatch's three branches, SolverKind and its
     members, material / material_table / BuildEnthalpyTable / SetPhaseFractionsConstAlpha, the global and closed-form
     branches of ConvertTtoH / InvertHtoT, UpdateModel's homogeneous branch, therm_alpha + its parse, the material
     fallbacks, the material + microstructure abort, the struct comment about the half term
  7  a deck without microstructure.nphases aborts naming P2-29
  8  retired keys: by presence material.* / implicit.* / phase.{T_ref,T_m,T_v,lh_m,lh_v}; at a default solver =
     explicit, surface_patch.foot_rule = pads, alpha = 0 0; the multi-token compare becomes token-wise numeric
  9  foot_rule max / mean deleted (the feed cap stays)
Usage: stage2_src.py   (writes src/Integrator/MMWSpalling.H; prints the line count before/after)"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PRE = os.path.join(HERE, "src_pre_2", "src", "Integrator", "MMWSpalling.H")
OUT = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
s = open(PRE).read()
n0 = s.count("\n")


def sub(a, b):
    global s
    assert s.count(a) == 1, "anchor not unique (%d): %r" % (s.count(a), a[:90])
    s = s.replace(a, b)


def span_end(start):
    """index just past the closing brace (and its newline) of the block whose first '{' follows `start`"""
    i = s.index("{", start)
    depth = 0
    while True:
        c = s[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return s.index("\n", i) + 1
        i += 1


def cut(first, last_sig, label):
    """delete from the line holding `first` to the end of the function whose signature is `last_sig`"""
    global s
    assert s.count(first) == 1 and s.count(last_sig) == 1, label
    a = s.rfind("\n", 0, s.index(first)) + 1
    b = span_end(s.index(last_sig))
    seg = s[a:b]
    L = seg.rstrip("\n").split("\n")
    print("  cut %-46s %4d lines: [%s] .. [%s]" % (label, len(L), L[0].strip()[:60], L[-1].strip()[:40]))
    s = s[:a] + s[b:]


# ---- 8: the retired keys (presence) after P2-28's, and the nphases abort (7)
sub('''        for (const char* family : {"bit", "spallation.sp.prescribed_stress",
                                   "spallation.sp.prescribed_T_surface"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                Util::Abort(INFO, k, " was retired in P2-28 (the column removal path is deleted: "
                            "the face list is the one removal path); remove it from the input");
''', '''        for (const char* family : {"bit", "spallation.sp.prescribed_stress",
                                   "spallation.sp.prescribed_T_surface"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                Util::Abort(INFO, k, " was retired in P2-28 (the column removal path is deleted: "
                            "the face list is the one removal path); remove it from the input");
        // P2-29 (census D): THE SEPARATE HEAT KERNELS are retired -- one heat path, AdvanceMicrostructure. The
        // global material (`material.*`), its implicit solver (`implicit.*`) and the const-alpha closed form's
        // phase scalars abort by PRESENCE; a deck states its rock as one phase.
        for (const char* k : {"phase.T_ref", "phase.T_m", "phase.T_v", "phase.lh_m", "phase.lh_v"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-29 (the const-alpha path is deleted: one heat path); state "
                            "the rock as one phase -- microstructure.nphases = 1, microstructure.phase0.{T_m, "
                            "fus_width, L_m, ...}; remove it from the input");
        for (const char* family : {"material", "implicit"})
            for (const std::string& k : amrex::ParmParse::getEntries(family))
                Util::Abort(INFO, k, " was retired in P2-29 (the separate heat kernels and the implicit solver are "
                            "deleted: one heat path); state the rock as one phase -- microstructure.nphases = 1, "
                            "microstructure.phase0.law = constant | grimsel | zhang; remove it from the input");
        // P2-29 (census D): a deck without microstructure.nphases ran a separate kernel; there is none now
        if (!pp.contains("microstructure.nphases"))
            Util::Abort(INFO, "microstructure.nphases is required since P2-29 (one heat path: the const-alpha and "
                        "material kernels are deleted): state the rock as one phase -- microstructure.mode = "
                        "expression, microstructure.phase_expr = 0, microstructure.nphases = 1, "
                        "microstructure.phase0.law = constant | grimsel | zhang");
''')
# ---- 8: the table rows at a default, and the token-wise numeric compare
sub('''            {"surface.follow_mask",  "1",
             "P2-28: the face list's `top` class IS the mask column top, and removal requires the list"},
        };''', '''            {"surface.follow_mask",  "1",
             "P2-28: the face list's `top` class IS the mask column top, and removal requires the list"},
            // P2-29 (census D): the separate kernels' selectors, each retired at the value the one heat path runs
            {"solver",                  "explicit",
             "P2-29: the implicit solver is deleted; the one heat path is forward Euler"},
            {"alpha",                   "0 0",
             "P2-29: the homogeneous eigenstrain is deleted; a phase's beta is the thermal expansion"},
            {"surface_patch.foot_rule", "pads",
             "P2-29: the max / mean feet rules were closed (no live deck); pads is the feet rule"},
        };''')
sub('''        // P2-27: a ONE-token value that parses fully as a double is compared as
        // a double, so `0.0` equals the table's "0" (P2-26 DR-3: the string
        // compare refused 622 record argvs that wrote a default as 0.0); any
        // other value (a word, several tokens) is compared as a string.
        auto one_double = [] (const std::string& s, double& d) {
            if (s.empty() || s.find(' ') != std::string::npos) return false;
            char* end = nullptr;
            d = std::strtod(s.c_str(), &end);
            return end == s.c_str() + s.size();
        };''', '''        // P2-27: a ONE-token value that parses fully as a double is compared as
        // a double, so `0.0` equals the table's "0" (P2-26 DR-3: the string
        // compare refused 622 record argvs that wrote a default as 0.0).
        // P2-29: TOKEN-WISE -- several tokens compare as doubles too, the same
        // count and each equal, so `0.0 0.0` equals the table's "0 0"; a word,
        // or a different count, is compared as a string (and differs).
        auto doubles = [] (const std::string& s, std::vector<double>& d) {
            std::istringstream in(s);
            std::string t;
            d.clear();
            while (in >> t)
            {
                char* end = nullptr;
                const double x = std::strtod(t.c_str(), &end);
                if (end != t.c_str() + t.size()) return false;
                d.push_back(x);
            }
            return !d.empty();
        };''')
sub('''            double d_got = 0.0, d_old = 0.0;
            if (one_double(got, d_got) && one_double(r.value, d_old) && d_got == d_old) continue;''',
    '''            std::vector<double> d_got, d_old;
            if (doubles(got, d_got) && doubles(r.value, d_old) && d_got == d_old) continue;''')
# ---- 6: the parses
sub('''        // Cache the diffusivity for use in Advance (HeatConduction's alpha is
        // private; cheaper than another queryclass).
        pp.query_required("hc.heat.alpha", value.therm_alpha, Unit::ThermalDiffusivity());

        // Phase-change parameters (single material in step 2; per-mineral in
        // step 5). Defaults disable phase change.
        pp.query_default("phase.T_ref", value.T_ref, 0.0);
        pp.query_default("phase.T_m",   value.T_m,   1.0e30);
        pp.query_default("phase.T_v",   value.T_v,   1.0e30);
        pp.query_default("phase.lh_m",  value.lh_m,  0.0);  // L_m / c_p [K]
        pp.query_default("phase.lh_v",  value.lh_v,  0.0);  // L_v / c_p [K]

''', '')
sub('''        // Optional material model. When omitted, the constant-α path of
        // step 2 is used (driven by hc.heat.alpha). When present, Advance
        // takes the variable-property flux-divergence path and the
        // enthalpy table is built from the material's ρ(T)·c_p(T).
        if (pp.contains("material.type"))
            pp.select_default<Numeric::Material::Zhang,
                              Numeric::Material::Constant>("material", value.material);

''', '')
sub('''        // Solver choice for the diffusion term: "explicit" forward Euler
        // (default; matches steps 2-4) or "implicit" backward Euler (built
        // on amrex::MLABecLaplacian + MLMG). The source Q and surface loss
        // q_loss stay explicit in both — they are not stiff once diffusion
        // is. Currently only the material path supports "implicit"; the
        // constant-α path is always explicit.
        std::string solver_str;
        pp.query_default("solver", solver_str, std::string("explicit"));
        if      (solver_str == "explicit") value.solver_kind = SolverKind::Explicit;
        else if (solver_str == "implicit") value.solver_kind = SolverKind::Implicit;
        else Util::Abort(INFO, "solver must be 'explicit' or 'implicit', got: ", solver_str);
        pp.query_default("implicit.tol_rel",  value.implicit_tol_rel,  1.0e-10);
        pp.query_default("implicit.tol_abs",  value.implicit_tol_abs,  0.0);
        pp.query_default("implicit.max_iter", value.implicit_max_iter, 100);
        pp.query_default("implicit.verbose",  value.implicit_verbose,  0);

        // Thermal-expansion eigenstrain coefficient per model component.
        pp_queryarr("alpha", value.alpha);

''', '''        // P2-29: `solver` (explicit) and `alpha` (0 0) are retired at their defaults (the table above): one heat
        // path, forward Euler; the thermal expansion is each phase's beta.

''')
# ---- 6: Initialize
sub('''        // Step 6: simultaneous `material` + `microstructure` is not yet
        // composed (the heterogeneous path uses constant per-phase ρ, c_p, κ;
        // the material path uses T-dependent polynomial ρ(T), c_p(T), κ(T)).
        // Reject the combination until a composition rule is defined.
        if (material && microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: simultaneous `material` and `microstructure` blocks "
                              "are not supported (Step 6 has no composition rule). Drop one.");

        if (material && !material_table.built()) BuildEnthalpyTable();   // build once
        HeatConduction::Initialize(lev);''', '''        HeatConduction::Initialize(lev);''')
# ---- 6: UpdateModel's homogeneous branch (after the microstructure branch's return)
a = s.index("        // Homogeneous path: existing alpha[]/eta[] thermal eigenstrain.")
b = s.index("    // P0b: declared here, ahead of its first use as a parameter type in")
seg = s[a:b]
assert seg.rstrip().endswith("}") and "F0 += (eta(i, j, k, n) * alpha[n])" in seg
print("  cut %-46s %4d lines" % ("UpdateModel's homogeneous branch", seg.count("\n")))
s = s[:a] + "        // P2-29: the homogeneous branch (the global alpha[] x eta eigenstrain) is deleted: every deck states\n" \
            "        // its rock as phases (microstructure.nphases is required), so the branch above is the only one.\n" \
            "    }\n\n" + s[b:]
sub('''            return;
        }

        // P2-29: the homogeneous branch''', '''            return;
        }
        // P2-29: the homogeneous branch''')
sub('''        // Step 8b: prescribed-temperature surface_patch acts as a flux source
        // on the cell-centred heat conduction (AdvanceMaterial), but the''',
    '''        // Step 8b: prescribed-temperature surface_patch acts as a flux source
        // on the cell-centred heat conduction (AdvanceMicrostructure), but the''')
# ---- 6: the dispatch and the three kernels
sub('''    // Override HeatConduction::Advance with the enthalpy method.
    // Branch: microstructure      → AdvanceMicrostructure (the production path);
    //         material + implicit → AdvanceMaterialImplicit (BE on diffusion);
    //         material + explicit → AdvanceMaterial (forward Euler);
    //         no material         → AdvanceConstAlpha (step 4 path).''',
    '''    // Override HeatConduction::Advance with the enthalpy method: ONE heat path, AdvanceMicrostructure
    // (P2-29: the const-alpha, material and implicit kernels are deleted; every deck states its rock as phases).''')
sub('''        if      (microstructure_enabled)                          AdvanceMicrostructure(lev, time, dt);
        else if (material && solver_kind == SolverKind::Implicit) AdvanceMaterialImplicit(lev, time, dt);
        else if (material)                                        AdvanceMaterial(lev, time, dt);
        else                                                      AdvanceConstAlpha(lev, time, dt);''',
    '''        AdvanceMicrostructure(lev, time, dt);''')
cut("    void AdvanceConstAlpha(int lev, Set::Scalar time, Set::Scalar dt)",
    "    void AdvanceMaterialImplicit(int lev, Set::Scalar time, Set::Scalar dt)", "the three kernels")
sub('''        // surface patch here as in AdvanceMaterial(): prescribed_T injects the''',
    '''        // surface patch here as the deleted material kernel did (P2-29): prescribed_T injects the''')
# ---- 6: the material fallbacks
sub('''            const Set::Scalar Tm = !phases.empty() ? phases[0].T_m : (material ? material->T_m : T_m);''',
    '''            const Set::Scalar Tm = !phases.empty() ? phases[0].T_m : Set::Scalar(1.0e30);''')
sub('''        if (!phases.empty()) { Tm = phases[0].T_m; w = phases[0].fus_width; }
        else if (material)   { Tm = material->T_m; w = material->fus_width; }''',
    '''        if (!phases.empty()) { Tm = phases[0].T_m; w = phases[0].fus_width; }''')
sub('''    // path, the global material's on a material path; 0 when no band is declared (T_m >= 1e29).''',
    '''    // path (P2-29: the only path); 0 when no band is declared (T_m >= 1e29).''')
# ---- 6: BuildEnthalpyTable + SetPhaseFractionsConstAlpha, and the converters' other branches
cut("    void BuildEnthalpyTable()", "    void SetPhaseFractionsConstAlpha(Set::Scalar T,", "BuildEnthalpyTable + SetPhaseFractionsConstAlpha")
a = s.index("    void ConvertTtoH(int lev)")
b = span_end(s.index("    void InvertHtoT(int lev)"))
old = s[a:b]
assert "SetPhaseFractionsConstAlpha" in old and "gtab" in old
new = '''    void ConvertTtoH(int lev)
    {
        // P2-29: one path -- the per-cell phase table (this->phase_tables[pid]); the global material's table and
        // the const-alpha closed form are deleted
        const Numeric::Material::Table *ptabs = phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        for (amrex::MFIter mfi(*H_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& H    = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& Hold = (*H_old_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& T = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph = (*phase_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Tv = T(i, j, k);
                Set::Scalar Hv = 0.0;
                Set::Scalar ls = 1.0, ll = 0.0, lv = 0.0;
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                Hv = ptabs[p].H_of_T(Tv);
                Set::Scalar dummyT;
                ptabs[p].Invert(Hv, dummyT, ls, ll, lv);
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
                H(i,j,k)    = Hv;
                Hold(i,j,k) = Hv;
            });
        }
    }

    void InvertHtoT(int lev)
    {
        // P2-29: one path -- the per-cell phase table
        const Numeric::Material::Table *ptabs = phase_tables.data();
        const int n_phases = static_cast<int>(phase_tables.size());
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& H = (*H_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& T  = (*temp_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LS = (*Lambda_S_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LL = (*Lambda_L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& LV = (*Lambda_V_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> ph = (*phase_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) {
                Set::Scalar Hv = H(i, j, k);
                int p = static_cast<int>(ph(i, j, k));
                if (p < 0 || p >= n_phases) p = 0;
                Set::Scalar Tv, ls, ll, lv;
                ptabs[p].Invert(Hv, Tv, ls, ll, lv);
                T(i,j,k) = Tv;
                LS(i,j,k) = ls; LL(i,j,k) = ll; LV(i,j,k) = lv;
            });
        }
    }
'''
print("  rewrite ConvertTtoH + InvertHtoT: %d -> %d lines (the phase-table branch kept, the rest deleted)"
      % (old.count("\n"), new.count("\n")))
s = s[:a] + new + s[b:]
# ---- 6: the members
sub('''    Set::Scalar therm_alpha = NAN; // κ/(ρ c_p), cached from hc.heat.alpha
    Set::Scalar T_ref = 0.0;
    Set::Scalar T_m   = 1.0e30;
    Set::Scalar T_v   = 1.0e30;
    Set::Scalar lh_m  = 0.0;
    Set::Scalar lh_v  = 0.0;

''', '')
sub('''    Numeric::Material::Material *material = nullptr; // optional T-dependent props

    // H↔T lookup table for the global `material` path. Step 6b moved the
    // implementation into a reusable Numeric::Material::Table; the
    // microstructure path uses one Table per phase (see phase_tables).
    Numeric::Material::Table material_table;

    enum class SolverKind { Explicit, Implicit };
    SolverKind  solver_kind        = SolverKind::Explicit;
    Set::Scalar implicit_tol_rel   = 1.0e-10;
    Set::Scalar implicit_tol_abs   = 0.0;
    int         implicit_max_iter  = 100;
    int         implicit_verbose   = 0;

    std::vector<Set::Scalar> alpha; // thermoelastic eigenstrain coefficient per component

''', '')
# ---- 6: the struct comment about the half term
sub('''    //     2*kappa*(T_patch - T_cell)/dz^2 on the in-patch top row (the
    //     microstructure kernel; P2-28 corrected the factor 2 there -- the
    //     material kernel's AdvanceMaterial term is Zhang's path, D, and still
    //     reads kappa*(...)/dz^2). Under the face list the top faces carry it.''',
    '''    //     2*kappa*(T_patch - T_cell)/dz^2 on the in-patch top row (P2-28
    //     corrected the factor 2; P2-29 deleted the material kernel that still
    //     had the half term). Under the face list the top faces carry it.''')
# ---- 9: foot_rule max / mean
sub('''        const bool pads = sp.foot_rule == SurfacePatch::FootRule::Pads;
        const int npads = pads ? sp.foot_npads : 1;
        std::vector<int> pad(zf.size(), 0);
        if (pads && npads > 1)''', '''        // P2-29: pads is the feet rule (max / mean, closed with no live deck, are deleted)
        const int npads = sp.foot_npads;
        std::vector<int> pad(zf.size(), 0);
        if (npads > 1)''')
sub('''            Set::Scalar hp = 0.0;
            if (sp.foot_rule == SurfacePatch::FootRule::Max)
                hp = v.back();
            else if (sp.foot_rule == SurfacePatch::FootRule::Mean)
            {
                for (Set::Scalar z : v) hp += z;
                hp /= Set::Scalar(v.size());
            }
            else
            {
                long idx = static_cast<long>(std::ceil(sp.foot_pad_quantile * Set::Scalar(v.size()))) - 1;
                idx = std::max(0L, std::min(idx, static_cast<long>(v.size()) - 1));
                hp = v[static_cast<std::size_t>(idx)];
            }''', '''            Set::Scalar hp = 0.0;
            {
                long idx = static_cast<long>(std::ceil(sp.foot_pad_quantile * Set::Scalar(v.size()))) - 1;
                idx = std::max(0L, std::min(idx, static_cast<long>(v.size()) - 1));
                hp = v[static_cast<std::size_t>(idx)];
            }''')
sub('''        enum class FootRule { Pads, Max, Mean };
        bool feet = false;
        FootRule foot_rule = FootRule::Pads;''', '''        bool feet = false;''')
sub('''                    std::string rule;
                    pp.query_default("foot_rule", rule, std::string("pads"));
                    if      (rule == "pads") v.foot_rule = FootRule::Pads;
                    else if (rule == "max")  v.foot_rule = FootRule::Max;
                    else if (rule == "mean") v.foot_rule = FootRule::Mean;
                    else Util::Abort(INFO, "surface_patch.foot_rule must be 'pads', 'max' or "
                                           "'mean', got: ", rule);
''', '''                    // P2-29: surface_patch.foot_rule is retired at `pads` (the table in Parse)
''')
sub('''                    v.foot_body_clearance = (cl == 1);
                    if (v.foot_body_clearance && v.foot_rule != FootRule::Pads)
                        Util::Abort(INFO, "surface_patch.foot_body_clearance = 1 requires "
                                          "surface_patch.foot_rule = pads");''', '''                    v.foot_body_clearance = (cl == 1);''')
# ---- the implicit solver's includes
sub('''#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
''', '')

open(OUT, "w").write(s)
print("MMWSpalling.H: %d -> %d lines" % (n0, s.count("\n")))
left = re.findall(r"\b(AdvanceConstAlpha|AdvanceMaterialImplicit|AdvanceMaterial|SolverKind|solver_kind|material_table|"
                  r"therm_alpha|BuildEnthalpyTable|SetPhaseFractionsConstAlpha|FootRule|implicit_tol_\w+)\b", s)
print("names left (code or comments): %s" % (sorted(set(left)) or "none"))
