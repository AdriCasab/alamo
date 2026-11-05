#!/usr/bin/env python3
"""Model card, pass 1: every resolved input key of each representative deck, read from the run's own AMReX
`metadata` file (# PARAMETERS block). A key listed twice (input file, then the command line) resolves to the LAST
value, as ParmParse does; the overridden values are kept in the `overridden` column. Writes RAW_KEYS.tsv
(key x deck) and DECKS.tsv (which run each deck label resolves to, its binary date and git hash).

Run: python3 docs/project/model_card/raw_keys.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
S = "tests/MMWSpalling/studies"
V = "tests/MMWSpalling/validation"
# label, metadata path, what it represents
DECKS = [
    ("KANT_onset", f"{V}/kant/sp_kant_onset_meiercfg/output/R500_1.0/run/metadata",
     "Kant 2017 onset anchor (flat flame, R500 1 mm leg) -- the rock-law calibration deck"),
    ("KANT_dref", f"{S}/p2_2e_cref/legs/G2.0_130/run/metadata",
     "Kant 2016 minimum-flux ladder, Gotthard properties, d_ref 2.0 mm, SOD 130 (P2-2e, not reviewed)"),
    ("HOSER_47", f"{S}/p2_11_heat/legs/G1/run/metadata",
     "Hoeser 47 MPa series deck, G1 (P2-11, accepted) -- today's Hoeser deck"),
    ("HOSER_0", f"{S}/p2_7_hoeser_heat/legs/HG/run/metadata",
     "Hoeser 0 MPa at her measured heat, Goldstein average profile (P2-7 HG, accepted)"),
    ("HOSER_0_S", f"{S}/whatif_uptake/legs/s100/run/metadata",
     "Hoeser 0 MPa what-if, HG with h x 0.10 (whatif_uptake, no /plan)"),
    ("MEIER", f"{S}/p2_2c_meier/ctrl8/m13/run/metadata",
     "Meier pilot, feet, Grimsel properties, m = 13 (P2-2f rev 2, accepted) -- last quotable Meier deck"),
    ("MEIER_dref", f"{S}/p2_2e_cref/legs/SD/run/metadata",
     "Meier with d_ref 2.0 mm (P2-2e SD, stalled at the pads, not reviewed)"),
    ("MMW_ZHANG", f"{V}/zhang/zhang_oglesby/output/metadata",
     "Zhang/Oglesby MMW thermal validation (beam on, no removal)"),
    ("HU_onset", f"{V}/hu/hu_spall_onset/output/granite2/metadata",
     "Hu granite-2 prescribed-temperature onset"),
]


def parse(path):
    vals, comp, inpar = {}, {}, False
    for ln in open(os.path.join(ROOT, path), errors="replace"):
        s = ln.rstrip("\n")
        if s.startswith("# PARAMETERS"):
            inpar = True
            continue
        if s.startswith("# RUN DETAILS"):
            inpar = False
            continue
        if s.startswith("#") or not s.strip():
            continue
        m = re.match(r"^(\S+?)(\(nvals = \d+\))?\s*(=|::)\s*(.*)$", s)
        if not m:
            continue
        k, v = m.group(1), m.group(4).strip()
        if inpar:
            vals.setdefault(k, []).append(v)
        else:
            comp[k] = v
    return vals, comp


def main():
    table, info = {}, []
    for lab, p, what in DECKS:
        if not os.path.exists(os.path.join(ROOT, p)):
            info.append((lab, "MISSING " + p, "", "", what))
            continue
        vals, comp = parse(p)
        table[lab] = vals
        info.append((lab, p, comp.get("Compilation_Date", "?"), comp.get("Git_commit_hash", "?"), what))
    keys = sorted(set(k for v in table.values() for k in v))
    labs = [d[0] for d in DECKS if d[0] in table]
    with open(os.path.join(HERE, "RAW_KEYS.tsv"), "w") as f:
        f.write("key\tdiffers\t" + "\t".join(labs) + "\toverridden\n")
        for k in keys:
            row = [table[l][k][-1] if k in table[l] else "-" for l in labs]
            present = [r for r in row if r != "-"]
            differs = "DIFF" if len(set(present)) > 1 else ""
            ov = "; ".join("%s:%s" % (l, "|".join(table[l][k][:-1])) for l in labs
                           if k in table[l] and len(table[l][k]) > 1 and len(set(table[l][k])) > 1)
            f.write("\t".join([k, differs] + row + [ov]) + "\n")
    with open(os.path.join(HERE, "DECKS.tsv"), "w") as f:
        f.write("label\tmetadata\tcompiled\tgit\trepresents\n")
        for r in info:
            f.write("\t".join(r) + "\n")
    print("decks %d, keys %d -> RAW_KEYS.tsv, DECKS.tsv" % (len(labs), len(keys)))


if __name__ == "__main__":
    main()
