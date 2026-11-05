"""Map every line of the integrator sources to its enclosing member function (heuristic).
A function starts at a 4- or 8-space-indented signature line whose body opens with '{'
at the same indentation. Output: work/funcmap.tsv  file  start  end  name"""
import re, os, json
ROOT = "/Users/tzetze20/amr_tools/alamo/"
FILES = ["src/Integrator/MMWSpalling.H", "src/Integrator/MMWSpalling/Removal.H",
         "src/Integrator/MMWSpalling/FaceRemoval.H"]
SIG = re.compile(r'^( {4}| {8})(?!return|if\b|for\b|while\b|switch\b|else\b|//|\*|case\b|using\b|#)'
                 r'[A-Za-z_][\w:<>,\s\*&\[\]]*?\b(~?\w+)\s*\(')
KEYW = {"if", "for", "while", "switch", "return", "sizeof", "static_assert", "Abort",
        "Message", "Warning", "Assert", "AMREX_ASSERT", "define", "defined"}
def scan(f):
    L = open(ROOT + f).read().split("\n")
    starts = []
    for i, l in enumerate(L):
        m = SIG.match(l)
        if not m or l.rstrip().endswith(";") or m.group(2) in KEYW:
            continue
        ind = m.group(1)
        # body must open with '{' at the same indentation within the next 12 lines
        for j in range(i, min(i + 14, len(L))):
            s = L[j]
            if s.rstrip().endswith(";") and j > i and not s.strip().startswith("//") and "(" not in s and ")" not in s:
                break
            if s.startswith(ind + "{") or (j == i and s.rstrip().endswith("{")) or \
               (s.rstrip().endswith("{") and j > i and s.startswith(ind) and not s.startswith(ind + " ")):
                starts.append((i + 1, m.group(2), len(ind)))
                break
            if s.rstrip().endswith(";") and j > i:
                break
    # end of a function = line where a '}' at the same indentation closes it
    out = []
    for (s, name, ind) in starts:
        e = None
        for j in range(s, len(L)):
            if L[j].startswith(" " * ind + "}") and not L[j].startswith(" " * (ind + 1)):
                e = j + 1
                break
        out.append((f, s, e or len(L), name, ind))
    return out
allf = []
for f in FILES:
    allf += scan(f)
with open(os.path.join(os.path.dirname(__file__), "funcmap.tsv"), "w") as fh:
    for r in allf:
        fh.write("%s\t%d\t%d\t%s\t%d\n" % r)
print(len(allf))
