#!/usr/bin/env python3
"""P2-31 stage 2a (goal 9): every item of the dead list z0 proved dead is removed (Z0.txt L8); the vapour keys retire at
their defaults (Parse's table); microstructure_enabled, its always-true tests and the seven unreachable aborts go.
Comments that named a removed identifier are reworded (R13 reads comments too). Every replacement must match exactly
once or NOTHING is written. Usage: stage2a_src.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
I = "src/Integrator/MMWSpalling/"

EDITS = {
    "src/Integrator/MMWSpalling.H": [
        ("        if (microstructure_enabled) InitializeMicrostructure(lev, /*announce=*/true);\n",
         "        InitializeMicrostructure(lev, /*announce=*/true);\n"),
        ("        if (microstructure_enabled) CheckTimeStepLimit();\n", "        CheckTimeStepLimit();\n"),
        ("""        if (microstructure_enabled && max_level > 0 &&
            m_type != Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);""",
         """        if (max_level > 0 &&
            m_type != Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable)
            RefreshMicrostructureDerivedFields(/*announce=*/false);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);"""),
        ("""        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepComplete(a_time, a_step);
        if (microstructure_enabled && max_level > 0 &&""",
         """        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepComplete(a_time, a_step);
        if (max_level > 0 &&"""),
        ("""        InitSurfaceMask(lev);
        if (microstructure_enabled)
        {
            InitializeMicrostructure(lev, /*announce=*/false);
            // H and T were interpolated independently; re-derive T and the
            // phase fractions from the conserved H on the remade level.
            if (lev > 0) InvertHtoT(lev);
        }""",
         """        InitSurfaceMask(lev);
        InitializeMicrostructure(lev, /*announce=*/false);
        // H and T were interpolated independently; re-derive T and the
        // phase fractions from the conserved H on the remade level.
        if (lev > 0) InvertHtoT(lev);"""),
    ],
    I + "Model.H": [
        ("""        if (microstructure_enabled)
        {
            // Step 7/alpha-beta: heterogeneous mechanics.""",
         """        {
            // Step 7/alpha-beta: heterogeneous mechanics."""),
        ("""        // the tangent depth k_c (T_s - T0) / q_net [m], set on the face list's
        // fallback lines (FaceRemoval.H) and read by nothing; 0 everywhere else
        Set::Scalar ell_q   = 0.0;
""", ""),
    ],
    I + "FaceRemoval.H": [
        ("profile -- NOT the tangent depth ell_q = k_c (T_s - T0) / q_net, which is",
         "profile -- NOT the tangent depth k_c (T_s - T0) / q_net, which is"),
        ("""the window has no depth: no recession (fr_nodepth). q_net is logged;
                    // the tangent depth ell_q is stored on the record and read by nothing.
""", """the window has no depth: no recession (fr_nodepth). q_net is logged.
"""),
        ("""                    const Set::Scalar q_net = a3 - a4;
                    const Set::Scalar k_c = Gc[2 * l + 1];
                    o.q_net = q_net;
                    if (q_net > 0.0 && Tf > T0l && k_c > 0.0)
                        o.rec.ell_q = k_c * (Tf - T0l) / q_net;   // the tangent depth (stored; read by nothing)
                    depth_from_Ts(Tf);""",
         """                    o.q_net = a3 - a4;
                    depth_from_Ts(Tf);"""),
    ],
    I + "Shared.H": [
        ("""    // nu is outside (-1, 1). `ok` = phase 0 exists and either the rock law is
    // selected or E > 0 and mu > 0 (nu is 0 otherwise). ⚠ No caller reads `ok`:
    // the one caller (RecedeAlongLines, FaceRemoval.H) takes the members as
    // they are.
    struct ElasticConsts { Set::Scalar E = 0.0, beta = 0.0, nu = 0.0, conf = 0.0; bool ok = false; };""",
         """    // nu is outside (-1, 1). nu is 0 unless phase 0 exists and either the rock
    // law is selected or E > 0 and mu > 0.
    struct ElasticConsts { Set::Scalar E = 0.0, beta = 0.0, nu = 0.0, conf = 0.0; };"""),
        ("""        // phase's thermal law is grimsel too (Parse), so they read 0 here and
        // `ok` does not ask for them: every consumer takes the law then""",
         """        // phase's thermal law is grimsel too (Parse), so they read 0 here and
        // the test below does not ask for them: every consumer takes the law then"""),
        ("""        if (!law && !(phases[0].E > 0.0 && phases[0].mu > 0.0)) return c;
        c.ok = true;
""", """        if (!law && !(phases[0].E > 0.0 && phases[0].mu > 0.0)) return c;
"""),
        ("""    // The guards: removal on (P2-28: or the face list's evaluate-only mode, whose
    // criterion runs here and removes nothing), the microstructure path, and the
    // surface-owning level only (AMR: coarser levels receive phi/removed by
    // average-down). Returns false when the step has nothing to do on this level.
    bool RemovalStepBegins(int lev)
    {
        if (!list_state()) return false;
        Util::Assert(INFO, TEST(microstructure_enabled));
        if (!SurfaceOwner(lev)) return false;""",
         """    // The guards: removal on (P2-28: or the face list's evaluate-only mode, whose
    // criterion runs here and removes nothing) and the surface-owning level only
    // (AMR: coarser levels receive phi/removed by average-down). Returns false
    // when the step has nothing to do on this level.
    bool RemovalStepBegins(int lev)
    {
        if (!list_state()) return false;
        if (!SurfaceOwner(lev)) return false;"""),
    ],
    I + "ParseRemoval.H": [
        ("""        value.vapor_enabled = (on != 0);
        if (!value.vapor_enabled) return;

        if (pp.contains("vapor.A_beam"))
        {
            value.vapor_A_beam_override = true;
            pp.query_required("vapor.A_beam", value.vapor_A_beam);
            if (!(std::isfinite(value.vapor_A_beam) && value.vapor_A_beam > 0.0))
                Util::Abort(INFO, "vapor.A_beam must be finite and > 0, got ",
                            value.vapor_A_beam);
        }
        if (pp.contains("vapor.T0"))
        {
            value.vapor_T0_override = true;
            pp.query_required("vapor.T0", value.vapor_T0);
            if (!std::isfinite(value.vapor_T0))
                Util::Abort(INFO, "vapor.T0 must be finite, got ", value.vapor_T0);
        }

        int presc = 0;
        pp.query_default("vapor.prescribed_rop.enabled", presc, 0);
        value.vapor_prescribed_rop_enabled = (presc != 0);
        if (value.vapor_prescribed_rop_enabled)
        {
            std::string expr;
            pp.query_required("vapor.prescribed_rop.expr", expr);
            value.vapor_rop_parser = amrex::Parser(expr);
            value.vapor_rop_parser.registerVariables({"x", "y", "z", "t"});
            value.vapor_rop_f = value.vapor_rop_parser.compile<4>();
        }
    }""",
         """        value.vapor_enabled = (on != 0);
        // P2-31: the beam area, the reference temperature and the prescribed recession rate the old vapour removal
        // read are retired at their defaults (Parse's table): the vapour rule is the face list's (P2-24, FaceRemoval.H)
        // and reads none of them
    }"""),
    ],
    I + "Parse.H": [
        ("""            {"surface_patch.foot_rule", "pads",
             "P2-29: the max / mean feet rules were closed (no live deck); pads is the feet rule"},
""", """            {"surface_patch.foot_rule", "pads",
             "P2-29: the max / mean feet rules were closed (no live deck); pads is the feet rule"},
            // P2-31: the keys of the vapour removal the face list's vapour rule (P2-24) replaced -- parsed, read by
            // nothing
            {"vapor.A_beam",                 "*",
             "P2-31: parsed into a member nothing read (the vapour rule is the face list's, P2-24)"},
            {"vapor.T0",                     "*",
             "P2-31: parsed into a member nothing read (the vapour rule is the face list's, P2-24)"},
            {"vapor.prescribed_rop.enabled", "0",
             "P2-31: the prescribed recession rate had no reader once the face list's vapour rule replaced it (P2-24)"},
            {"vapor.prescribed_rop.expr",    "*",
             "P2-31: read only under the retired vapor.prescribed_rop.enabled = 1"},
"""),
        ("""            for (auto &p : value.phases) p.fraction /= fsum;

            value.microstructure_enabled = true;
""", """            for (auto &p : value.phases) p.fraction /= fsum;
"""),
        ("""        if (value.vapor_enabled && !value.vapor_A_beam_override)
        {
            // Default beam-area convention: Gaussian effective area
            //   A_beam = 0.5 * pi * omega0^2
            // (matches the (1-R)*P0/A_eff intensity scaling used by the
            // Gaussian source, where peak intensity is 2P/(pi*omega0^2)).
            constexpr Set::Scalar pi = 3.141592653589793238462643383279502884;
            value.vapor_A_beam = 0.5 * pi * value.beam.omega0 * value.beam.omega0;
            if (!(std::isfinite(value.vapor_A_beam) && value.vapor_A_beam > 0.0))
                Util::Abort(INFO, "vapor: default A_beam is non-positive (",
                            value.vapor_A_beam, "); set vapor.A_beam explicitly");
        }
        if (value.vapor_enabled && !value.vapor_T0_override)
            value.vapor_T0 = value.T_amb;

        // P2-27: LAST""", """        // P2-27: LAST"""),
        ("""        if (value.spall_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spall.enabled requires microstructure.nphases");

""", ""),
        ("""        if (value.vapor_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: vapor.enabled requires microstructure.nphases");

        if (value.spallation_model == SpallationModel::SpWeibull &&
            !value.microstructure_enabled)
            Util::Abort(INFO, "MMWSpalling: spallation.model = sp_weibull requires "
                              "microstructure.nphases (Step 15 Weibull samples will key "
                              "on grain-boundary faces)");

        if (value.ledger_enabled && !value.microstructure_enabled)
            Util::Abort(INFO, "energy_ledger.enabled requires "
                              "the microstructure path (microstructure.nphases)");
""", ""),
        ("""        if (value.surface_patch.enabled &&
            value.surface_patch.robin_form == SurfacePatch::RobinForm::Face &&
            !value.microstructure_enabled)
            Util::Abort(INFO, "surface_patch.robin_form = face is implemented on "
                              "the microstructure path only");
""", ""),
        ("""        // P0b deliverable 1: robin_form = surface needs the microstructure
        // path, the face list (its reconstruction) and the criterion.
        if (value.surface_patch.enabled &&
            value.surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
        {
            if (!value.microstructure_enabled)
                Util::Abort(INFO, "surface_patch.robin_form = surface is implemented "
                                  "on the microstructure path only");
""", """        // P0b deliverable 1: robin_form = surface needs the face list (its
        // reconstruction) and the criterion.
        if (value.surface_patch.enabled &&
            value.surface_patch.robin_form == SurfacePatch::RobinForm::Surface)
        {
"""),
        ("""            if (!value.microstructure_enabled)
                Util::Abort(INFO, "surface_patch.face_list = 1 requires the "
                            "microstructure path (its consumers live in "
                            "AdvanceMicrostructure)");
""", ""),
    ],
    I + "State.H": [
        ("""    bool microstructure_enabled = false;
""", ""),
        ("""    // ⚠ vapor.A_beam, vapor.T0 and vapor.prescribed_rop.* are still parsed
    // into the members below, and nothing outside the parse reads them.
    bool vapor_enabled = false;
    bool vapor_A_beam_override = false;
    Set::Scalar vapor_A_beam = 0.0;
    bool vapor_T0_override = false;
    Set::Scalar vapor_T0 = 0.0;
    bool vapor_prescribed_rop_enabled = false;
    amrex::Parser vapor_rop_parser;
    amrex::ParserExecutor<4> vapor_rop_f;
""", """    bool vapor_enabled = false;
"""),
        ("""    // only). Default-off; `weibull.enabled = 1` requires sp_weibull +
    // microstructure (the latter is already enforced by the sp_weibull +
    // !microstructure_enabled abort). `flaw_a_mf` is cell-centred: one""",
         """    // only). Default-off; `weibull.enabled = 1` requires sp_weibull (every
    // deck has the microstructure: microstructure.nphases is required).
    // `flaw_a_mf` is cell-centred: one"""),
        ("""        // the microstructure requirement is enforced by the existing
        // sp_weibull + !microstructure_enabled abort. With""",
         """        // every deck has the microstructure (microstructure.nphases is
        // required). With"""),
    ],
    I + "BeamProps.H": [
        ("""    {
        if (!microstructure_enabled) return;
        bool any = false;""", """    {
        bool any = false;"""),
    ],
    I + "Fields.H": [
        ("        if (microstructure_enabled && max_level > 0)\n", "        if (max_level > 0)\n"),
    ],
    "src/Numeric/Material/Grimsel.H": [
        ("""// ⚠ H(theta) here (`H_of`, J/kg, zero at 20 degC) is ANALYTIC and is the exact
//   integral of c_p. The integrator does not call it: a run's H<->T inversion
//   is the phase's Table, a trapezoid of rho*c_p built from this class's
//   rho() and cp() (Table::Build).""",
         """// ⚠ H(theta) [48] is the exact integral of c_p, and the class does not form it:
//   a run's H<->T inversion is the phase's Table, a trapezoid of rho*c_p built
//   from this class's rho() and cp() (Table::Build)."""),
        ("""    // The analytic enthalpy [48], exactly the integral of cp above.
    Set::Scalar H_of(Set::Scalar T) const
    {
        const Set::Scalar th = T - 273.15;
        return cp_a * (th - 20.0)
             - cp_b * std::log((th + cp_th0) / (20.0 + cp_th0));
    }

""", ""),
    ],
    "src/Numeric/MMWBeam.H": [
        ("""    bool active() const { return parsed && P0 != 0.0; }

""", ""),
    ],
    "src/Numeric/Material/Table.H": [
        ("""        m_fus_lo_T = fus_lo;
        m_fus_hi_T = fus_hi;
        m_vap_hi_T = mat.T_vap_hi;
        m_T_v      = mat.T_vap_lo; // sharp-vap reference T (cached only: T_of_H pins to m_vap_lo_T, nothing reads this)

""", ""),
        ("""    Set::Scalar m_fus_lo_T = 0.0, m_fus_hi_T = 0.0;
    Set::Scalar m_vap_lo_T = 0.0, m_vap_hi_T = 0.0;
    Set::Scalar m_T_m = 0.0, m_T_v = 0.0;""",
         """    Set::Scalar m_vap_lo_T = 0.0;
    Set::Scalar m_T_m = 0.0;"""),
    ],
}


def main():
    new = {}
    for rel, reps in EDITS.items():
        p = os.path.join(ROOT, rel)
        t = open(p).read()
        for old, rep in reps:
            n = t.count(old)
            if n != 1:
                sys.exit("%s: expected 1 match, found %d of:\n%s" % (rel, n, old[:300]))
            t = t.replace(old, rep, 1)
        new[p] = t
    for p, t in new.items():
        open(p, "w").write(t)
    print("2a applied to %d files: %s" % (len(new), ", ".join(sorted(os.path.relpath(p, ROOT) for p in new))))


if __name__ == "__main__":
    main()
