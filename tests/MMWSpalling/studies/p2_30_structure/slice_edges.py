#!/usr/bin/env python3
"""P2-30: what the token identity cannot see at a slice edge (a stranded comment, macro or access specifier).
In the parent, the non-blank line above each NEW slice include is not a // line, an AMREX_GPU_* line or template<; each
slice ends on } / }; / ; and begins on a comment, an AMREX_GPU_* line or a member; neither access specifier is inside a
slice. The two older fragments (Shared.H, FaceRemoval.H) keep the explanatory comment above their includes (M:4737-4739),
which is not a stranded member comment; they are named and skipped.
Exit 1 on any violation."""
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))
DIR = os.path.join(ROOT, "src", "Integrator", "MMWSpalling")
OLD = ("Shared.H", "FaceRemoval.H")


def main():
    P = open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")).read().split("\n")
    bad = []
    names = []
    for i, ln in enumerate(P):
        m = re.match(r'\s*#include "Integrator/MMWSpalling/([^"]+)"', ln)
        if not m:
            continue
        nm = m.group(1)
        if nm in OLD:
            print("  skipped (older fragment, its explanatory comment stands): %s" % nm)
            continue
        names.append(nm)
        j = i - 1
        while j >= 0 and P[j].strip() == "":
            j -= 1
        above = P[j].strip()
        if above.startswith("//") or above.startswith("AMREX_GPU") or above.startswith("template<"):
            bad.append("%s: the line above its include in the parent is %r" % (nm, above[:70]))
        L = open(os.path.join(DIR, nm)).read().split("\n")
        while L and L[-1] == "":
            L.pop()
        if not (L[-1].strip() in ("}", "};") or L[-1].rstrip().endswith(";")):
            bad.append("%s ends on %r" % (nm, L[-1].strip()[:70]))
        first = L[0].strip()
        if not (first.startswith("//") or first.startswith("AMREX_GPU") or
                re.match(r"(static|void|struct|template|enum|std::|Set::|int|bool|const|BeamTau|amrex::|inline)\b", first)):
            bad.append("%s begins on %r" % (nm, first[:70]))
        for k, s in enumerate(L, 1):
            if re.match(r"\s*(public|protected|private)\s*:\s*(//.*)?$", s) and not s.startswith("        "):
                # an access label at member depth (a struct's own label sits deeper and belongs to that struct)
                bad.append("%s:%d carries an access specifier: %r" % (nm, k, s.strip()))
    specs = [i + 1 for i, s in enumerate(P) if re.match(r"^(public|protected|private)\s*:", s)]
    print("  %d new slices: %s" % (len(names), ", ".join(names)))
    print("  access specifiers in the parent at lines %s (%s)" % (specs, ", ".join(P[i - 1].strip() for i in specs)))
    if len(specs) != 3 or [P[i - 1].strip() for i in specs] != ["public:", "protected:", "private:"]:
        bad.append("the parent should hold public / protected / private: found %s" % specs)
    for b in bad:
        print("  VIOLATION  " + b)
    print("VERDICT", "no edge violation" if not bad else "%d violation(s)" % len(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
