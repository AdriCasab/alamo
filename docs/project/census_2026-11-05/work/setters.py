"""Who sets each key: production legs (effective values), kept unit tests, validation
decks, tree-gate inputs. Output: setters.json"""
import os, re, glob, json, shlex
import keys
ROOT = "/Users/tzetze20/amr_tools/alamo/"
T = ROOT + "tests/MMWSpalling/"
H = os.path.dirname(os.path.abspath(__file__))

def deck_pairs(path, seen=None):
    """(key, value) in order, FILE includes expanded in place (AMReX: later wins)."""
    out = []
    for ln in open(path, errors="replace"):
        s = ln.split("#", 1)[0].strip()
        if not s or "=" not in s:
            continue
        k, v = (x.strip() for x in s.split("=", 1))
        if k == "FILE":
            out += deck_pairs(ROOT + v)
            continue
        out.append((k, v))
    return out
def argv_pairs(path):
    toks = shlex.split(open(path).read())
    i = next(n for n, t in enumerate(toks) if "/tests/MMWSpalling/" in t and "=" not in t and not t.endswith("g++"))
    base = deck_pairs(toks[i])
    return base + [tuple(t.split("=", 1)) for t in toks[i + 1:] if "=" in t]
PROD_LEGS = {
    "HOS_R": argv_pairs(T + "studies/p2_16_span_rule/legs/HOS_R/argv.sh"),
    "CAL": argv_pairs(T + "studies/p2_16_span_rule/legs/CAL/argv.sh"),
}
PROD_FILES = ["validation/shared/rock_grimsel.inputs", "validation/shared/model_form_flame.inputs",
              "validation/kant/grimsel_calibration/input", "validation/hoeser/input_lab"]

def textfiles_unit():
    fs = []
    for d in sorted(glob.glob(T + "unit/*/")):
        for f in sorted(os.listdir(d)):
            p = d + f
            if os.path.isdir(p) or f.endswith((".log", ".pyc", ".txt", ".md")):
                continue
            if f.startswith("input") or f == "test" or f.endswith(".py"):
                fs.append(p)
    return fs
def files_validation():
    fs = []
    for sub in ["hu", "zhang", "meier", "rossi", "kant", "hoeser", "p0_steady_column",
                "spall_thickness", "grimsel_properties"]:
        for p in glob.glob(T + "validation/" + sub + "/**/input*", recursive=True):
            if "/output" in p or os.path.isdir(p) or p.endswith((".log", ".md", ".txt", ".py", ".csv")):
                continue
            fs.append(p)
    return sorted(fs)
GATE_TESTS = ["robin_face", "robin_pinned", "robin_jet", "jet_enthalpy", "robin_feet",
              "beam_void_closure", "scalar_flaw", "spall_event", "jet_d2c"]
def files_gate():
    fs = []
    for t in GATE_TESTS:
        for f in sorted(os.listdir(T + "unit/" + t)):
            if f.startswith("input") or f == "test":
                fs.append(T + "unit/" + t + "/" + f)
    fs += [T + "validation/meier/sp_meier_pilot/input_2d_dev",
           T + "studies/s1_surface_resolution/input", T + "studies/s1_surface_resolution/run_sweep.py"]
    return fs

def key_regex(key):
    if "phase<i>" in key:
        k = re.escape(key).replace(re.escape("phase<i>"), r"phase\d+")
    else:
        k = re.escape(key)
    return [re.compile(r'(?<![\w.])(' + k + r')\s*=\s*("[^"]*"|\'[^\']*\'|[^\s"\',\]\)#]+)'),
            re.compile(r'["\'](' + k + r')["\']\s*:\s*("[^"]*"|\'[^\']*\'|[^\s,}]+)')]
def scan(files, key, py_only_strings=False):
    hits = []
    pats = key_regex(key)
    for p in files:
        try:
            L = open(p, errors="replace").read().split("\n")
        except Exception:
            continue
        isdeck = os.path.basename(p).startswith("input") or p.endswith(".inputs")
        for i, ln in enumerate(L):
            s = ln.split("#", 1)[0] if isdeck else ln
            if isdeck and s.strip().startswith("#"):
                continue
            for pat in pats:
                for m in pat.finditer(s):
                    hits.append((p.replace(T, ""), i + 1, m.group(1), m.group(2).strip("\"'")))
    return hits

def nondefault(val, default):
    d = default.strip()
    if d.upper().startswith("REQUIRED") or d.startswith("ABORT") or d.startswith("absent"):
        return True
    try:
        return abs(float(val) - float(d)) > 1e-12 * max(1.0, abs(float(d)))
    except ValueError:
        pass
    if d in ("inf",) and val in ("inf",):
        return False
    if d.startswith("=") or " " in d or "(" in d or "/" in d:
        return None   # not comparable
    return val != d

out = {}
UF, VF, GF = textfiles_unit(), files_validation(), files_gate()
for r in keys.K:
    key = r["key"]
    pk = re.compile(key_regex(key)[0].pattern.split(r'\s*=')[0].replace("(?<![\\w.])(", "^(") + "$")
    prod = {}
    for leg, pairs in PROD_LEGS.items():
        vals = [(k, v) for k, v in pairs if pk.match(k)]
        prod[leg] = vals[-1][1] if vals else None
        if len(vals) > 1:
            prod[leg + "_all"] = [v for _, v in vals]
    pfile = scan([T + f for f in PROD_FILES], key)
    u = scan(UF, key)
    v = scan(VF, key)
    g = scan(GF, key)
    def ann(h):
        return [dict(file=f, line=l, key=k, val=val, nd=nondefault(val, r["default"])) for f, l, k, val in h]
    out[key] = dict(prod=prod, prod_files=ann(pfile), unit=ann(u), val=ann(v), gate=ann(g))
json.dump(out, open(os.path.join(H, "setters.json"), "w"), indent=0)
print(len(UF), len(VF), len(GF))
for t in ["UF", "VF", "GF"]:
    pass
open(os.path.join(H, "files_scanned.txt"), "w").write(
    "UNIT\n" + "\n".join(x.replace(T, "") for x in UF) + "\nVALIDATION\n" + "\n".join(x.replace(T, "") for x in VF)
    + "\nGATE\n" + "\n".join(x.replace(T, "") for x in GF) + "\n")
