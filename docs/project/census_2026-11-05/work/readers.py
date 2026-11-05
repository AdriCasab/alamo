"""For every key in keys.py, find the reads of the member(s) it sets, outside the parse
functions. Comments and string literals are stripped first. Output: readers.json"""
import re, json, os, collections
import keys
ROOT = "/Users/tzetze20/amr_tools/alamo/"
H = os.path.dirname(os.path.abspath(__file__))
FILES = {"MMWSpalling.H": "src/Integrator/MMWSpalling.H",
         "Removal.H": "src/Integrator/MMWSpalling/Removal.H",
         "FaceRemoval.H": "src/Integrator/MMWSpalling/FaceRemoval.H"}
# parse ranges (MMWSpalling.H): Parse, ParseDamage, ParseCohesive, ParseSpall, ParseVapor,
# ParseBit, CheckPadRing, SurfacePatch::Parse, ParseSpallationModel
PARSE = [(88, 1126), (5013, 5110), (5273, 5324), (5400, 5811), (5813, 5847), (5854, 5877),
         (7887, 7895), (7897, 8356), (9137, 9372)]
SP_STRUCT = (7786, 8357)
funcs = collections.defaultdict(list)
for l in open(os.path.join(H, "funcmap.tsv")):
    f, s, e, n, ind = l.rstrip("\n").split("\t")
    funcs[os.path.basename(f)].append((int(s), int(e), n))
def func_of(fb, ln):
    best = None
    for s, e, n in funcs[fb]:
        if s <= ln <= e and (best is None or s > best[0]):
            best = (s, n)
    return best[1] if best else None
def strip(line):
    line = re.sub(r'"(\\.|[^"\\])*"', '""', line)
    i = line.find("//")
    return line[:i] if i >= 0 else line
SRC = {fb: [strip(x) for x in open(ROOT + p).read().split("\n")] for fb, p in FILES.items()}
def in_parse(fb, ln):
    return fb == "MMWSpalling.H" and any(a <= ln <= b for a, b in PARSE)
def pats(kind, mem):
    if kind == "m":
        if mem == "alpha":
            return [re.compile(r'(?<![\w\.>])alpha\s*\[')]
        return [re.compile(r'(?<![\w\.>])' + mem + r'\b(?!\s*\()')]
    if kind == "sp":
        return [re.compile(r'\b(?:surface_patch|sp)\.' + mem + r'\b')]
    if kind == "ph":
        return [re.compile(r'phases\[[^\]]*\]\.' + mem + r'\b(?!\s*\()'),
                re.compile(r'(?<![\w\.])(?:p|ph|phase|mp)\.' + mem + r'\b(?!\s*\()')]
    if kind == "beam":
        return [re.compile(r'\bbeam\.' + mem + r'\b')]
    return []
out = {}
for r in keys.K:
    hits = []
    for mem in r["members"]:
        for p in pats(r["kind"], mem):
            for fb, L in SRC.items():
                for i, line in enumerate(L):
                    ln = i + 1
                    if p.search(line):
                        fn = func_of(fb, ln)
                        tag = "parse" if in_parse(fb, ln) else ("decl" if fn is None else "read")
                        hits.append(dict(mem=mem, file=fb, line=ln, func=fn, tag=tag, text=line.strip()[:160]))
            # SurfacePatch members used bare inside the struct's own methods
            if r["kind"] == "sp":
                q = re.compile(r'(?<![\w\.>])' + mem + r'\b')
                L = SRC["MMWSpalling.H"]
                for ln in range(SP_STRUCT[0], SP_STRUCT[1] + 1):
                    if q.search(L[ln - 1]) and not any(a <= ln <= b for a, b in PARSE):
                        hits.append(dict(mem=mem, file="MMWSpalling.H", line=ln, func="SurfacePatch(struct)",
                                         tag="decl", text=L[ln - 1].strip()[:160]))
    out[r["key"]] = hits
json.dump(out, open(os.path.join(H, "readers.json"), "w"), indent=0)
# summary
for r in keys.K:
    h = [x for x in out[r["key"]] if x["tag"] == "read"]
    c = collections.Counter(x["func"] for x in h)
    print(r["key"], len(h), dict(c))
