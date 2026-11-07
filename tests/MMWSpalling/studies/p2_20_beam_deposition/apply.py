#!/usr/bin/env python3
"""P2-20 Goal 5 -- the decks carry the dielectric description, never the typed alpha / R.

  usage: apply.py --dry      list every line that would go or come; change nothing
         apply.py --apply    do it; writes DECK_EDITS.txt beside this script

DECKS = unit/, validation/, extensions/ (input*, *.inputs) + studies/s1_surface_resolution/input (a gate deck), the
same scope as P2-18's a2_apply.py. NOT touched: output dirs, retired_units, */pre/*, src_pre, */legs/*, argv.sh
records, study run.py scripts (records) -- except s1's run_sweep.py, which the tree gate runs (c).
The Zhang/Oglesby deck is edited by Goal 6 (by hand, recorded there), not here.

  (a) a beam-OFF deck (P0 = 0 and no schedule): its beam.alpha_expr / beam.R lines go (inert at P = 0)
  (b) a beam-ON FIXTURE: beam.alpha_expr -> beam.eps_r = 1.0 + beam.tan_delta = <z0's fixture string> with the
      comment "# fixture: alpha = <old> 1/m by construction (eps_r 1 => R 0)"; its beam.R line (0) goes
  (c) s1's run_sweep.py: beam.alpha_expr=<a> on argv -> beam.eps_r=1.0 beam.tan_delta=<fixture string> (lambda is
      the s1 deck's default 1.0e-3)
  (d) the 13 Hu beam extensions: the Barre pair (eps_r 5.02, tan_delta 0.01065 => alpha 68.5 1/m, R 0.147 at their
      2.19 mm), commented as typed-then-replaced; unscored and not run here
  (e) validation/shared/rock_grimsel.inputs: the dielectric block (class LIT, Barre, flat in T)
"""
import importlib.machinery
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
ld = importlib.machinery.SourceFileLoader("z0_p220", os.path.join(HERE, "z0.py"))
Z0 = importlib.util.module_from_spec(importlib.util.spec_from_loader("z0_p220", ld))
ld.exec_module(Z0)

SKIP = ("/output", "/legs/", "/run/", ".old", "retired_units", "/pre/", "src_pre", "inputs_pre", "/restore/",
        "/zhang_pre", "zhang_oglesby")
RE_A = re.compile(r"^\s*beam\.alpha_expr\s*=\s*(.*?)\s*(#.*)?$")
RE_R = re.compile(r"^\s*beam\.R\s*=\s*(.*?)\s*(#.*)?$")
S1_SWEEP = os.path.join(T, "studies", "s1_surface_resolution", "run_sweep.py")
ROCK = os.path.join(T, "validation", "shared", "rock_grimsel.inputs")
BARRE = ("beam.eps_r     = 5.02      # P2-20: Barre granite (DOE DE-EE0005504 p. 52: n 2.24); typed alpha 100 /"
         " R 0.2 replaced -- unscored, not run\n",
         "beam.tan_delta = 0.01065   # P2-20: Barre, 30 % single-pass absorption in 25 mm at 28 GHz (DOE p. 44)"
         " => alpha 68.5 1/m, R 0.147 at lambda 2.19 mm\n")
ROCK_BLOCK = '''
# --- the dielectric description (P2-20): the beam's absorption and reflection ---------
beam.eps_r     = 5.02       # LIT (another rock: Barre granite): eps' = n^2, n = 2.24 (Oglesby, Woskov et al.,
                            #   DOE final report DE-EE0005504, Dec 2014, OSTI 1169951, p. 52)
beam.tan_delta = 0.01065    # LIT (Barre): 30 % single-pass absorption through 25 mm of solid granite at 28 GHz
                            #   (DOE p. 44/52) => POWER alpha = -ln(0.7)/0.025 ~ 14 1/m; tan_delta = alpha lambda/(2 pi n)
                            # The code computes alpha = 2 pi sqrt(eps_r) tan_delta / lambda (the POWER coefficient)
                            #   and R = ((sqrt eps_r - 1)/(sqrt eps_r + 1))^2 = 0.147 (the DOE's Fresnel 0.14) with
                            #   each deck's own beam.lambda (MMWBeam.H). FLAT with temperature by the user's decision
                            #   (2026-11-05). To be replaced by a Grimsel measurement.
'''


def decks():
    out = []
    # P2-21 (P2-20 review DR-1): studies/ too -- two live study decks kept the retired keys and aborted
    for sub in ("unit", "validation", "extensions", "studies"):
        for dp, _ds, fs in os.walk(os.path.join(T, sub)):
            if any(s in dp + "/" for s in SKIP):
                continue
            for f in fs:
                if f.startswith("input") or f.endswith(".inputs"):
                    if not f.endswith((".py", ".md", ".txt", ".csv")):
                        out.append(os.path.join(dp, f))
    out.append(os.path.join(T, "studies", "s1_surface_resolution", "input"))
    return sorted(set(out))


def kv(text):
    d = {}
    for ln in text.splitlines():
        ln = ln.split("#")[0]
        if "=" in ln:
            k, v = ln.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def beam_on(d):
    try:
        p0 = float(d.get("beam.P0", "0").split()[0])
    except ValueError:
        p0 = 1.0
    return p0 != 0.0 or "beam.schedule.P" in d


def plan():
    edits, files = [], {}
    for p in decks():
        text = open(p, errors="replace").read()
        if not (re.search(r"^\s*beam\.alpha_expr\s*=", text, re.M) or re.search(r"^\s*beam\.R\s*=", text, re.M)):
            continue
        rel = os.path.relpath(p, ROOT)
        d = kv(text)
        on = beam_on(d)
        hu = "/extensions/hu_spall_onset_mmwbeam" in p
        lam_s = d.get("beam.lambda", "1.0e-3").split()[0]
        new = []
        for ln in text.splitlines(keepends=True):
            ma, mr = RE_A.match(ln), RE_R.match(ln)
            if not (ma or mr):
                new.append(ln)
                continue
            if mr:
                if on and not hu:
                    assert float(mr.group(1)) == 0.0, "%s: a beam-ON fixture with R != 0" % rel
                edits.append((rel, "-", ln.rstrip("\n")))
                continue
            # alpha_expr
            edits.append((rel, "-", ln.rstrip("\n")))
            if not on:
                continue
            if hu:
                for b in BARRE:
                    new.append(b)
                    edits.append((rel, "+", b.rstrip("\n")))
                continue
            a_old = float(ma.group(1))                   # a pure number on every fixture (asserted by float())
            s, exact, k, got = Z0.fixture_string(a_old, lam_s)
            assert exact, "%s: no bit-exact tan_delta for alpha %r at lambda %s" % (rel, a_old, lam_s)
            for b in ("beam.eps_r     = 1.0\n",
                      "beam.tan_delta = %s   # fixture: alpha = %s 1/m by construction (eps_r 1 => R 0); lambda %s\n"
                      % (s, ma.group(1), lam_s)):
                new.append(b)
                edits.append((rel, "+", b.rstrip("\n")))
        files[p] = "".join(new)
    # (c) s1's sweep script
    sw = open(S1_SWEEP).read()
    old = '''               f"beam.alpha_expr={src['alpha']!r}",\n'''
    assert sw.count(old) == 1
    alphas = sorted({2000.0, 32000.0, 100.0})
    tds = {a: Z0.fixture_string(a, "1.0e-3") for a in alphas}
    assert all(v[1] for v in tds.values())
    sw = sw.replace(old, '''               "beam.eps_r=1.0", f"beam.tan_delta={TD_FIXTURE[src['alpha']]}",\n''')
    anchor = "CASES = {"
    assert sw.count(anchor) == 1
    tdl = "".join("    %r: %r,\n" % (a, tds[a][0]) for a in alphas)
    sw = sw.replace(anchor, "# P2-20: beam.alpha_expr is retired. Each beam case's alpha is reproduced BIT-FOR-BIT by eps_r = 1\n"
                    "# (=> R 0) and this tan_delta at the s1 deck's default lambda 1.0e-3 (MMWBeam.H: alpha =\n"
                    "# 2 pi sqrt(eps_r) * (tan_delta / lambda); studies/p2_20_beam_deposition/Z0.txt).\n"
                    "TD_FIXTURE = {\n" + tdl + "}\n\n" + anchor)
    edits.append((os.path.relpath(S1_SWEEP, ROOT), "-", old.strip()))
    edits.append((os.path.relpath(S1_SWEEP, ROOT), "+", '"beam.eps_r=1.0", f"beam.tan_delta={TD_FIXTURE[src[\'alpha\']]}"'
                  " + TD_FIXTURE " + repr({a: tds[a][0] for a in alphas})))
    files[S1_SWEEP] = sw
    # (e) the rock file
    rk = open(ROCK).read()
    assert "beam.eps_r" not in rk
    anchor = "# --- the surface's emissivity (radiation to the surroundings) --------------\n"
    assert rk.count(anchor) == 1
    files[ROCK] = rk.replace(anchor, ROCK_BLOCK.lstrip("\n") + "\n" + anchor)
    for ln in ROCK_BLOCK.strip("\n").splitlines():
        edits.append((os.path.relpath(ROCK, ROOT), "+", ln))
    return edits, files


def main(a):
    edits, files = plan()
    by = {}
    for rel, sg, ln in edits:
        by.setdefault(rel, []).append("  %s %s" % (sg, ln))
    L = ["P2-20 deck edits (%d files, %d lines)" % (len(by), len(edits)), "=" * 100]
    for rel in sorted(by):
        L.append(rel)
        L += by[rel]
    txt = "\n".join(L)
    if "--apply" in a:
        for p, t in files.items():
            open(p, "w").write(t)
        open(os.path.join(HERE, "DECK_EDITS.txt"), "w").write(txt + "\n")
        print("applied: %d files (DECK_EDITS.txt)" % len(files))
    else:
        print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
