#!/usr/bin/env python3
"""check_shared.py DECK [ARGV-FILE] [--allow KEY ...] [--run DIR]

Does the run use the shared Grimsel files as they are written? Four checks (P2-13 a2):
  (i)   no key defined in rock_grimsel.inputs or model_form_flame.inputs is set in the deck's own
        input or on its argv, except keys passed with --allow (a control leg's DECLARED override);
  (ii)  the run's `metadata` resolves every shared key (last value wins) to the file's value
        (an --allow key: to its override value). A key the run never queried (e.g. weibull.* under a declared
        weibull.enabled=0) is absent from `metadata`; it passes only if the run's log (run.log beside the run
        directory) lists it among AMReX's "Unused ParmParse Variables" WITH the file's value -- proof the file
        was read and the key was simply not consulted on that branch;
  (iii) microstructure.phase0.T_ref == hc.ic.constant.value (the surface reconstruction's T0 is T_ref);
  (iv)  under weibull.enabled = 1, weibull.V0 equals the deck's cell volume to 1e-12 relative.
The run directory is the resolved plot_file (argv over deck), or --run DIR. (ii)-(iv) read its
`metadata`, so they need a run (one step is enough). Prints one line per check; exit 1 on a failure.
"""
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
SHARED = [os.path.join(HERE, "rock_grimsel.inputs"), os.path.join(HERE, "model_form_flame.inputs")]


def parse(path):
    """ParmParse-style: [(key, value-string)] in file order, comments stripped; FILE lines kept as ('FILE', path)"""
    out = []
    for ln in open(path):
        ln = ln.split("#", 1)[0].strip()
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        out.append((k.strip(), " ".join(v.split())))
    return out


def parse_argv(path):
    toks = shlex.split(open(path).read())
    return [tuple(t.split("=", 1)) for t in toks if "=" in t and not t.startswith("-")]


def parse_meta(path):
    d = {}
    for ln in open(path):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            d[k.strip()] = " ".join(v.split())        # last value wins
    return d


def same(a, b):
    ta, tb = a.split(), b.split()
    if len(ta) != len(tb):
        return False
    for x, y in zip(ta, tb):
        try:
            fx, fy = float(x), float(y)
            if not (fx == fy or abs(fx - fy) <= 1e-12 * max(abs(fx), abs(fy))):
                return False
        except ValueError:
            if x != y:
                return False
    return True


def main(argv):
    allow, run, pos = set(), None, []
    i = 0
    while i < len(argv):
        if argv[i] == "--allow":
            i += 1
            while i < len(argv) and not argv[i].startswith("--"):
                allow.add(argv[i]); i += 1
            continue
        if argv[i] == "--run":
            run = argv[i + 1]; i += 2
            continue
        pos.append(argv[i]); i += 1
    if not pos:
        sys.exit(__doc__)
    deck = pos[0]
    deck_kv = parse(deck)
    arg_kv = parse_argv(pos[1]) if len(pos) > 1 else []
    shared = {}
    for f in SHARED:
        for k, v in parse(f):
            shared[k] = v
    bad = []

    def line(tag, ok, msg):
        print("  %-5s %s  %s" % (tag, "PASS" if ok else "FAIL", msg))
        if not ok:
            bad.append(tag)

    print("check_shared %s%s%s" % (os.path.relpath(deck, ROOT), "  argv " + os.path.relpath(pos[1], ROOT)
                                   if len(pos) > 1 else "", "  allow " + " ".join(sorted(allow)) if allow else ""))
    inc = [os.path.normpath(os.path.join(ROOT, v)) for k, v in deck_kv if k == "FILE"]
    line("incl", all(os.path.normpath(f) in inc for f in SHARED),
         "the deck includes both shared files (%d FILE lines)" % len(inc))

    # (i)
    own = [(k, "deck") for k, v in deck_kv if k in shared and k != "FILE"]
    own += [(k, "argv") for k, v in arg_kv if k in shared]
    undeclared = [(k, w) for k, w in own if k not in allow]
    unused = sorted(a for a in allow if a not in [k for k, _ in own])
    line("(i)", not undeclared and not unused,
         "no undeclared override of a shared key (%d declared%s%s)"
         % (len([1 for k, _ in own if k in allow]),
            "; UNDECLARED: " + ", ".join("%s (%s)" % kw for kw in undeclared) if undeclared else "",
            "; --allow names a key nobody sets: " + ", ".join(unused) if unused else ""))

    # the run
    if run is None:
        pf = dict(deck_kv).get("plot_file")
        for k, v in arg_kv:
            if k == "plot_file":
                pf = v
        run = pf if pf is None or os.path.isabs(pf) else os.path.join(ROOT, pf)
    meta_p = os.path.join(run or "", "metadata")
    if not os.path.exists(meta_p):
        line("(ii)", False, "no metadata at %s (run the deck first)" % meta_p)
        print("check_shared FAIL")
        return 1
    meta = parse_meta(meta_p)

    # (ii)
    want = dict(shared)
    for k, v in arg_kv + [kv for kv in deck_kv if kv[0] in allow]:
        if k in allow:
            want[k] = v
    unused = {}
    logp = os.path.join(os.path.dirname(os.path.normpath(run)), "run.log")
    if os.path.exists(logp):
        tail = open(logp, errors="replace").read().split("Unused ParmParse Variables:")
        if len(tail) > 1:
            for ln in tail[-1].splitlines():
                if "::" in ln and "[" in ln:
                    k = ln.split("::", 1)[1].split("(", 1)[0].strip()
                    unused[k] = " ".join(ln.rsplit("[", 1)[1].rstrip("]").split())
    absent = [k for k in want if k not in meta]
    unread = [k for k in absent if k in unused and same(unused[k], want[k])]
    miss = [k for k in absent if k not in unread]
    wrong = ["%s: %s != %s" % (k, meta[k], want[k]) for k in want if k in meta and not same(meta[k], want[k])]
    line("(ii)", not miss and not wrong,
         "metadata resolves all %d shared keys to the file's value%s%s%s"
         % (len(want) - len(unread), " (+%d not queried on this branch, read from the file: %s)"
            % (len(unread), ", ".join(unread)) if unread else "", "; MISSING: " + ", ".join(miss) if miss else "", "; WRONG: " + "; ".join(wrong) if wrong else ""))

    # (iii)
    tr, ic = meta.get("microstructure.phase0.T_ref"), meta.get("hc.ic.constant.value")
    line("(iii)", tr is not None and ic is not None and same(tr, ic),
         "phase0.T_ref %s == hc.ic.constant.value %s" % (tr, ic))

    # (iv)
    if meta.get("weibull.enabled", "0").split()[0] in ("1", "true"):
        lo = [float(x) for x in meta["geometry.prob_lo"].split()]
        hi = [float(x) for x in meta["geometry.prob_hi"].split()]
        n = [int(x) for x in meta["amr.n_cell"].split()]
        vc = 1.0
        for a, b, m in zip(lo, hi, n):
            vc *= (b - a) / m
        v0 = float(meta["weibull.V0"])
        line("(iv)", abs(v0 - vc) <= 1e-12 * vc, "weibull.V0 %.6g == V_cell %.6g (rel %.1e)" % (v0, vc, abs(v0 - vc) / vc))
    else:
        line("(iv)", True, "weibull.enabled = 0: V0 not read (n/a)")
    print("check_shared %s" % ("PASS" if not bad else "FAIL: " + " ".join(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
