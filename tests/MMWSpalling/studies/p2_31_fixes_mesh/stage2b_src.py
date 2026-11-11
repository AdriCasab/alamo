#!/usr/bin/env python3
"""P2-31 stage 2b (goal 10): `regime_field` duplicates `regime` at every write site (1 spall, 2 vapour, 4 pad clearance,
FaceRemoval.H) and has NO live test reader (Z0.txt L15), so it is removed: the registration, the member, the three
writes, the per-step reset. The one reader outside L15's scope (extensions/hu_spall_onset_mmwbeam_lefm/test, an ungated
informational print that tolerates the field's absence) is pointed at `regime`, which holds the same values.
Every replacement must match exactly once or NOTHING is written."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
I = "src/Integrator/MMWSpalling/"

EDITS = {
    I + "State.H": [
        ("""    // Step 17: cell-resolved per-firing-cell regime diagnostic under
    // sp_weibull. {0 none, 1 spall, 2 vap, 4 pad clearance}. Registered only when
    // spallation_model == SpWeibull && (spall_enabled || vapor_enabled).
    // Distinct from the scalar `regime_mf` (registered earlier; kept).
    Set::Field<Set::Scalar> regime_field_mf;
""", ""),
    ],
    I + "Shared.H": [
        ("""        // Step 17: per-step reset for the per-cell regime field under
        // sp_weibull. The field is registered only when
        // spallation_model == SpWeibull && (spall_enabled || vapor_enabled).
        if (spallation_model == SpallationModel::SpWeibull &&
            (spall_enabled || vapor_enabled))
        {
            regime_field_mf[lev]->setVal(0.0);
        }
""", ""),
    ],
    I + "Parse.H": [
        ("""            // Step 17: cell-resolved regime field under sp_weibull,
            // registered after the scalar `regime_mf` (kept, in its order).
            // The face list writes the SAME code to both at every removal
            // (1 spall, 2 vapour, 4 pad clearance; FaceRemoval.H).
            if (value.spallation_model == SpallationModel::SpWeibull &&
                (value.spall_enabled || value.vapor_enabled))
            {
                value.RegisterNewFab(value.regime_field_mf, nullptr, 1, 0,
                                     "regime_field", true);
            }
""", ""),
    ],
    I + "FaceRemoval.H": [
        ("    // state: fr_f 0, removed 1, regime / regime_field 4, RoP n dz / dt; the\n",
         "    // state: fr_f 0, removed 1, regime 4, RoP n dz / dt; the\n"),
        ("            amrex::Array4<Set::Scalar> const& regf = (*regime_field_mf[lev]).array(mfi);\n", ""),
        ("                        regf(i, j, k) = 4.0;\n", ""),
        ("                    amrex::Array4<Set::Scalar> const& regf= (*regime_field_mf[lev]).array(mfi);\n", ""),
        ("                            regf(p[0], p[1], p[2]) = 1.0;\n", ""),
        ("                            regf(c.i, c.j, c.k) = 2.0;\n", ""),
    ],
    "tests/MMWSpalling/extensions/hu_spall_onset_mmwbeam_lefm/test": [
        ("""        if "regime_field" in fields and out["n_removed"] > 0:
            reg = to_numpy(ad["regime_field"])""",
         """        if "regime" in fields and out["n_removed"] > 0:   # P2-31: regime_field (its duplicate) is removed
            reg = to_numpy(ad["regime"])"""),
        ("""        if "regime_field" in {fld[1] for fld in ds_local.field_list}:
            reg = to_numpy(ad_local["regime_field"])""",
         """        if "regime" in {fld[1] for fld in ds_local.field_list}:
            reg = to_numpy(ad_local["regime"])"""),
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
    print("2b applied to %d files: %s" % (len(new), ", ".join(sorted(os.path.relpath(p, ROOT) for p in new))))


if __name__ == "__main__":
    main()
