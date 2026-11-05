import json, collections, keys
S = json.load(open("setters.json"))
R = json.load(open("readers.json"))
def short(p):
    return (p.replace("validation/", "").replace("unit/", "u/").replace("studies/s1_surface_resolution/", "s1/")
             .replace("sp_meier_pilot/", "").replace("sp_kant_onset", "kant_onset"))
def uniq(vals, n=4):
    v = list(dict.fromkeys(vals))
    return ",".join(v[:n]) + (",..." if len(v) > n else "")
def who(key):
    s = S[key]
    p = s["prod"]
    prod = "; ".join("%s=%s" % (k, p[k]) for k in ("HOS_R", "CAL") if p.get(k) is not None) or "—"
    if p.get("CAL_all"): prod += " (CAL: %s, last wins)" % "→".join(p["CAL_all"])
    if p.get("HOS_R_all"): prod += " (HOS_R: %s, last wins)" % "→".join(p["HOS_R_all"])
    def grp(lst, by_dir):
        d = collections.OrderedDict()
        for x in lst:
            parts = x["file"].split("/")
            g = parts[1] if by_dir else short(x["file"])
            d.setdefault(g, []).append(x["val"])
        return d
    u = grp(s["unit"], True)
    uu = "%d: %s" % (len(u), " ".join("%s(%s)" % (g, uniq(v)) for g, v in u.items())) if u else "—"
    v = grp(s["val"], False)
    vv = "%d: %s" % (len(v), " ".join("%s(%s)" % (g, uniq(x)) for g, x in v.items())) if v else "—"
    gnd = [x for x in s["gate"] if x["nd"] is not False]
    g = grp(gnd, False)
    gg = " ".join("%s(%s)" % (k, uniq(x)) for k, x in g.items()) if g else "—"
    return prod, uu, vv, gg
def readers(key):
    h = [x for x in R[key] if x["tag"] == "read"]
    c = collections.OrderedDict()
    for x in h:
        f = {"MMWSpalling.H": "M", "Removal.H": "R", "FaceRemoval.H": "FR"}[x["file"]]
        c.setdefault(x["func"], []).append("%s:%d" % (f, x["line"]))
    return h, c
if __name__ == "__main__":
    for r in keys.K:
        prod, u, v, g = who(r["key"])
        h, c = readers(r["key"])
        print("##", r["key"])
        print("  prod:", prod); print("  unit:", u); print("  val :", v); print("  gate:", g)
