#!/usr/bin/env python3
"""Census part 4 helper (read-only): per study directory under tests/MMWSpalling/studies and extensions, which
key families its recorded argv (legs/**/argv.sh, any argv*.sh / *argv*.txt) and its drivers (run*.py, *.sh) name.
Static scan only: a retired key (P2-10 list, P2-16 constraint_depth) in a recorded argv means that leg ABORTS
today; a strict surface_patch.* key the parser no longer reads also aborts (queryclass strict check).
Writes STUDIES_SCAN.tsv next to this script."""
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
SP_KEYS = set(open(os.path.join(HERE, "sp_keys.txt")).read().split())
P210 = re.compile(r"\b(spall\.lateral_(?:drive|removal|subcell|thermal)|spall\.per_cell_removal|"
                  r"surface_patch\.(?:side_face_(?:flux|h|T|h_r_min|h_law|h_datum|factor_above|area_law|"
                  r"step_height|staircase)|step_diag_radius|jet_wall_march|jet_wall_mdot|jet_body(?:_T|_h_scale)?))=")
P216 = re.compile(r"\bspallation\.sp\.constraint_depth=")
FAM = {
    "DAMAGE": re.compile(r"\b(damage\.(?:enabled|ic\.type|[a-z_]+)|cohesive\.[a-z_]+|spallation\.model=damage_law|"
                         r"spall\.sample=(?:column_max|shallow_max)|spall\.(?:prescribed_ratio|alpha|t_spall|xi_s|C_h|"
                         r"damage_threshold))="),
    "ZHANG": re.compile(r"\b(material\.type|solver=implicit)"),
    "FACE": re.compile(r"surface_patch\.face_list=1"),
    "surface_normal": re.compile(r"spall\.surface_normal=1"),
    "melt_aware": re.compile(r"spall\.melt_aware=1"),
    "flake_coh": re.compile(r"spall\.flake_coherence_length=(?!0(?:\.0*)?\b)"),
    "flake_len_flux": re.compile(r"spall\.flake_length=flux|spall\.flake_fraction=weibull"),
    "pad_shield": re.compile(r"surface_patch\.pad_shield_open="),
    "jet_above": re.compile(r"surface_patch\.jet_above_h=1|surface_patch\.jet_T_above="),
    "rec_fixed": re.compile(r"jet_T_ent_mode=rec_fixed"),
    "robin_cell_face": re.compile(r"surface_patch\.robin_form=(?:cell|face)\b"),
    "driving_stress": re.compile(r"spall\.driving_(?:stress|ref|ref_depth)=|spall\.support_"),
}


def tokens_of(path):
    try:
        txt = open(path, errors="replace").read()
    except OSError:
        return ""
    return txt


rows = []
dirs = sorted(glob.glob(os.path.join(T, "studies", "*"))) + sorted(glob.glob(os.path.join(T, "extensions", "*")))
for d in dirs:
    if not os.path.isdir(d):
        continue
    argvs = [p for p in glob.glob(os.path.join(d, "**", "argv*.sh"), recursive=True)
             + glob.glob(os.path.join(d, "**", "*argv*.txt"), recursive=True)
             if "/output" not in p.replace(d, "")[:0] + p]
    drivers = [p for p in glob.glob(os.path.join(d, "**", "*.py"), recursive=True)
               + glob.glob(os.path.join(d, "**", "*.sh"), recursive=True)
               + glob.glob(os.path.join(d, "**", "input*"), recursive=True)
               if os.path.isfile(p) and "__pycache__" not in p and "/run/" not in p and not p.endswith(".log")]
    drivers = [p for p in drivers if p not in argvs]
    a_txt = "\n".join(tokens_of(p) for p in argvs)
    d_txt = "\n".join(tokens_of(p) for p in drivers)
    # normalise "key = value" in input files and quoted "key=value" strings in drivers to key=value
    norm = lambda s: re.sub(r"([A-Za-z_][A-Za-z_0-9.]*)\s*=\s*", r"\1=", s)
    a_n, d_n = norm(a_txt), norm(d_txt)
    n_abort_210 = sum(1 for p in argvs if P210.search(norm(tokens_of(p))))
    n_abort_216 = sum(1 for p in argvs if P216.search(norm(tokens_of(p))))
    strict = set()
    for m in re.finditer(r"surface_patch\.([A-Za-z_0-9]+)=", a_n):
        k = m.group(1)
        if k not in SP_KEYS and not P210.search("surface_patch.%s=" % k):
            strict.add(k)
    n_face = sum(1 for p in argvs if FAM["FACE"].search(norm(tokens_of(p))))
    fams = []
    for name, rx in FAM.items():
        if name == "FACE":
            continue
        a, b = bool(rx.search(a_n)), bool(rx.search(d_n))
        if a or b:
            fams.append(name + ("" if a else "(driver)"))
    rows.append((os.path.relpath(d, T), len(argvs), n_abort_210, n_abort_216, ",".join(sorted(strict)),
                 n_face, " ".join(fams), bool(P210.search(d_n)), bool(P216.search(d_n))))

with open(os.path.join(HERE, "STUDIES_SCAN.tsv"), "w") as f:
    f.write("dir\tn_argv\targv_abort_P2-10\targv_abort_P2-16\tstrict_sp_unknown\targv_face_list=1\tfamilies\t"
            "driver_names_P2-10\tdriver_names_P2-16\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
for r in rows:
    print("\t".join(str(x) for x in r))
