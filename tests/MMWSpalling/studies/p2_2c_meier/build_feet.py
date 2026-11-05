#!/usr/bin/env python3
"""P2-2c deliverable 3 + a3 — build the feet configuration and RESOLVE EVERY KEY.

THE BUILD IS BY SUBTRACTION, NOT ADDITION (guardrail): take `CTRL5.argv.sh`,
REMOVE `surface_patch.nozzle_z0` and `surface_patch.nozzle_feed`, SWITCH
`surface_patch.nozzle_descent` from `prescribed` to `feet`, and ADD the seven
remaining feet keys from `input_feet_d2c:166-173`.

⚠⚠ a3 IS A GATE: every key must resolve to a read site AND its guard. A leg
that aborts or runs inert after an overnight launch is what this prevents.

Usage: build_feet.py            # prints the ledger and writes feet_argv.txt
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
CTRL5 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_1_wallh",
                     "CTRL5.argv.sh")
FEET_DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                         "sp_meier_pilot", "input_feet_d2c")
SRC = [os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H"),
       os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "Removal.H")]

DROP = ("surface_patch.nozzle_z0", "surface_patch.nozzle_feed")
GRIMSEL = ["microstructure.phase0.law=grimsel",
           "spall.stress_law=grimsel",
           "spallation.sp.kic.law=grimsel_table",
           "spallation.sp.kic.scale=1.0",
           "spallation.sp.stress_relief=0.531294"]


def feet_keys():
    """The eight feet keys, READ FROM THE DECK -- never transcribed."""
    out = []
    for ln in open(FEET_DECK).read().splitlines()[165:173]:
        ln = ln.split("#")[0].strip()
        if "=" in ln:
            k, v = [x.strip() for x in ln.split("=", 1)]
            out.append("%s=%s" % (k, v))
    return out


def build():
    # ⚠⚠ `CTRL5.argv.sh` IS A SHELL COMMAND LINE, NOT A KEY LIST, AND NAIVE
    # SPLITTING BREAKS IT TWO WAYS.
    #   (1) QUOTES: `surface_patch.h_expr` is single-quoted (its value is a
    #       parser expression), so a whitespace split leaves a leading `'` on
    #       the token -- it fails every startswith test and would be handed to
    #       ALAMO verbatim. `shlex` removes them correctly.
    #   (2) MULTI-VALUE KEYS: `amr.n_cell=60 60 300`,
    #       `geometry.prob_lo=0.0 0.0 -0.1` and `geometry.prob_hi=...` carry
    #       THREE values each. A split on whitespace keeps only the first and
    #       orphans the rest as bare tokens with no `=`, which the key filter
    #       then DROPS. The 60 s leg aborted with
    #       "Invalid domain specified: 0 < y < 0.12" because prob_lo had lost
    #       its y and z components. ⇒ re-attach every bare token to the key
    #       that precedes it.
    import shlex
    toks = shlex.split(open(CTRL5).read())
    keys = []
    def is_number(x):
        try:
            float(x); return True
        except ValueError:
            return False
    for t in toks:
        # ⚠⚠ A NEGATIVE VALUE LOOKS LIKE A FLAG. `geometry.prob_lo=0.0 0.0 -0.1`
        # carries z_lo = -0.1 -- 100 mm of rock BELOW the origin -- and a filter
        # that drops anything starting with `-` eats it, leaving prob_lo with
        # two components and silently truncating the domain. Treat a token as a
        # flag only if it is not a number.
        if (t.startswith("-") and not is_number(t)) or \
           t.endswith(("mpirun", "caffeinate")):
            continue
        if "=" in t:
            keys.append(t)
        elif keys and not os.path.sep in t:
            keys[-1] = keys[-1] + " " + t     # a further value of the last key
    keys = [k for k in keys if "=" in k]
    # strip the harness's own leg plumbing; the runner supplies these
    keys = [k for k in keys if not k.startswith(("plot_file=", "stop_time="))]
    dropped = [k for k in keys if k.split("=")[0] in DROP]
    keys = [k for k in keys if k.split("=")[0] not in DROP]
    fk = feet_keys()
    # the SWITCH, in place -- not an append
    switched = None
    for i, k in enumerate(keys):
        if k.startswith("surface_patch.nozzle_descent="):
            switched = (k, "surface_patch.nozzle_descent=feet")
            keys[i] = switched[1]
    added = [k for k in fk
             if not k.startswith("surface_patch.nozzle_descent")]
    return keys + added + GRIMSEL, dropped, switched, added


def resolve(key):
    """(read-site, guard) for a key, derived against live source.

    ⚠ TWO SPELLINGS MUST BOTH MATCH. A key is parsed either by its FULL dotted
    name in one string literal (`"spallation.sp.kic.law"`) or by a PREFIX
    concatenation (`pfx + ".law"`), and a resolver that only tries the last
    segment misses the first kind while a resolver that only tries the full
    name misses the second. The first version here tried only the bare
    segment and reported five live keys as NOT FOUND."""
    bare = key.split(".")[-1]
    pats = [r'"%s"' % re.escape(key),            # full dotted name
            r'"%s"' % re.escape(bare),           # bare segment
            r'"\.%s"' % re.escape(bare)]         # pfx + ".bare"
    site, guard = None, None
    for p in SRC:
        txt = open(p).read().splitlines()
        for i, ln in enumerate(txt):
            if any(re.search(q, ln) for q in pats) and \
               re.search(r"query|contains|select", ln):
                site = "%s:%d" % (os.path.basename(p), i + 1)
                # the nearest enclosing `if (...)` above the parse
                for j in range(i, max(0, i - 40), -1):
                    m = re.match(r"\s*if \((.+?)\)\s*$", txt[j])
                    if m and "query" not in m.group(1):
                        guard = m.group(1)[:58]
                        break
                return site, guard
    return site, guard


def main():
    keys, dropped, switched, added = build()
    print("P2-2c deliverable 3 / a3 -- THE FEET CONFIGURATION AND ITS LEDGER")
    print("=" * 74)
    print("\nbuilt from %s by SUBTRACTION:" % os.path.relpath(CTRL5, ROOT))
    for k in dropped:
        print("  REMOVED  %s" % k)
    if switched:
        print("  SWITCHED %s  ->  %s" % switched)
    for k in added:
        print("  ADDED    %s   (input_feet_d2c:166-173)" % k)
    for k in GRIMSEL:
        print("  ADDED    %s   (the Grimsel laws)" % k)
    print("\n  total keys: %d" % len(keys))

    print("\n[a3] THE REACHABILITY LEDGER -- read site AND guard")
    # ⚠⚠ THE `guard` COLUMN IS A HEURISTIC AND MUST BE READ AS ONE. It reports
    # the nearest enclosing `if (...)` above the parse, which for a key parsed
    # inside a long block picks up an UNRELATED condition -- e.g. it attributes
    # `pp.contains("microstructure.grain_expr")` to `phase0.law`, which is
    # simply the previous `if` in the same block. A wrong guard is worse than
    # no guard, so the AUTHORITATIVE part of this row is the explicit guard
    # CONJUNCTION block below, which tests the three conditions that can
    # actually abort a feet leg. Treat the column as a pointer, not a fact.
    print("  %-44s %-26s %s" % ("key", "read site", "guard (heuristic)"))
    unresolved = []
    for k in sorted(set(added + GRIMSEL)):
        name = k.split("=")[0]
        site, guard = resolve(name)
        print("  %-44s %-26s %s" % (name, site or "** NOT FOUND **",
                                    (guard or "(unguarded / top level)")[:40]))
        if site is None:
            unresolved.append(name)
    print()
    if unresolved:
        print("  FAIL  %d key(s) do not resolve: %s"
              % (len(unresolved), ", ".join(unresolved)))
    else:
        print("  ok    every added key resolves to a read site")

    # ⚠⚠ THE CONJUNCTION THE GUARDS REQUIRE, checked explicitly.
    print("\n  THE GUARD CONJUNCTION (what feet mode demands):")
    mmw = open(SRC[0]).read()
    checks = [
        ("nozzle_z0/nozzle_feed MUST BE ABSENT under feet",
         'pp.contains("nozzle_z0") || pp.contains("nozzle_feed")' in mmw,
         not any(k.split("=")[0] in DROP for k in keys)),
        ("feet REQUIRES surface_patch.h_expr",
         "nozzle_descent = feet requires the flame" in mmw,
         any(k.startswith("surface_patch.h_expr") for k in keys)),
        ("foot_body_clearance REQUIRES feet",
         "foot_body_clearance requires" in mmw,
         any(k.startswith("surface_patch.nozzle_descent=feet") for k in keys)),
    ]
    ok_all = True
    for name, guard_exists, satisfied in checks:
        print("    %s guard present=%s  argv satisfies=%s  %s"
              % ("ok  " if (guard_exists and satisfied) else "FAIL",
                 guard_exists, satisfied, name))
        ok_all = ok_all and guard_exists and satisfied

    out = os.path.join(HERE, "feet_argv.txt")
    with open(out, "w") as f:
        f.write("\n".join(keys) + "\n")
    print("\n  written: %s" % os.path.relpath(out, ROOT))
    return 0 if (not unresolved and ok_all) else 1


if __name__ == "__main__":
    sys.exit(main())
