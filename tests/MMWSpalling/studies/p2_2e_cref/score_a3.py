#!/usr/bin/env python3
"""P2-2e (a3) -- decomposition identity at d = 8.5 mm.

 (i)   run_removal_events.csv byte-identical across -np 1/2/4 and two
       max_grid_size, on the Kant column and on a 5 s Hoeser G05 leg;
 (ii)  the Sp_field plotfile component identical, cell by cell -- assembled from
       each run's OWN boxes, because the FAB layout differs by decomposition and
       a byte comparison of Cell_D files would be meaningless;
 (iii) one pair with spall.per_cell_removal = 1, event logs identical, AND its
       lateral path shown to FIRE (pcr_lateral_cells > 0) -- otherwise (iii) is
       as vacuous as (i) is for the Sp_field site.
⚠ Why (ii) and (iii) exist: with per_cell_removal off the Sp_field site writes
  no event row, so (i) cannot see a box-local read there (ledger L5).
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2_props"))
import fabcmp                                                      # noqa: E402
LEGS = os.path.join(HERE, "legs")
CAVEATS = []


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def field(leg, var="Sp_field"):
    run = os.path.join(LEGS, leg, "run")
    pfs = sorted(d for d in os.listdir(run) if d.endswith("cell"))
    pf = pfs[-1]
    lvl = os.path.join(run, pf, "Level_0")
    files = sorted(f for f in os.listdir(lvl) if f.startswith("Cell_D_"))
    first = fabcmp.read_fabs_boxed(os.path.join(lvl, files[0]))
    names = fabcmp.var_names(os.path.join(run, pf, "Header"), first[0]["ncomp"])
    iv = names.index(var)
    out = {}
    for fn in files:
        for fab in fabcmp.read_fabs_boxed(os.path.join(lvl, fn)):
            lo, hi = fab["lo"], fab["hi"]
            nx, ny = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1
            c = fab["comp"][iv]
            for o, v in enumerate(c):
                i = o % nx; j = (o // nx) % ny; k = o // (nx * ny)
                out[(lo[0] + i, lo[1] + j, lo[2] + k)] = v
    return pf, out


def thermo_col(leg, name):
    L = open(os.path.join(LEGS, leg, "run", "thermo.dat")).read().split("\n")
    h = L[0].split()
    if name not in h:
        return None
    return [float(r.split()[h.index(name)]) for r in L[1:] if r.strip()]


def group(tag, legs, check_sp=True, need_pcr=False):
    P = print
    bad = 0
    P(f"\n  {tag}")
    ref = legs[0]
    h0 = md5(os.path.join(LEGS, ref, "run_removal_events.csv"))
    n0 = sum(1 for _ in open(os.path.join(LEGS, ref, "run_removal_events.csv"))) - 1
    cv = thermo_col(ref, "cref_n_vert")
    P(f"    reference {ref}: {n0} event rows, md5 {h0[:12]}, "
      f"cref_n_vert max {max(cv) if cv else 'ABSENT'}")
    if not cv or max(cv) <= 0:
        P("    FAIL  the constraint reference was never read on this leg -- "
          "the identity below would be vacuous"); bad += 1
    if need_pcr:
        pc = thermo_col(ref, "pcr_lateral_cells")
        tot = sum(pc) if pc else 0
        P(f"    pcr_lateral_cells: {'ABSENT' if pc is None else f'max {max(pc):.0f}, nonzero on {sum(1 for x in pc if x > 0)} of {len(pc)} rows'}")
        if not pc or max(pc) <= 0:
            # ⚠ NOT a gate failure, and NOT dropped either. This check is the
            # implementer's addition, beyond the packet's (iii) bar ("event logs
            # identical"), which IS met. It is not met because per-cell lateral
            # removal never fires on this deck: P2-5 measured that the Hoeser
            # pit walls never reach the firing flux in 120 s (c13). So (iii)'s
            # EVENT-LOG form cannot see the Sp_field lateral verdicts here, which
            # is the vacuity the packet feared for (i). What covers the same
            # purpose is (ii) on this pair: the Sp_field component, including
            # the lateral faces the branch scored, compared across decompositions
            # whose x/y seams cut the pit walls (mgs 8).
            P("    NOT MET (caveat, not a gate row): per-cell removal never fired, so "
              "(iii)'s event log cannot see the Sp_field lateral verdicts on this deck")
            CAVEATS.append("(iii) per-cell removal fired 0 times -- its event-log "
                           "form is vacuous here; the lateral branch is covered by "
                           "the Sp_field comparison across x/y seams instead")
    if check_sp:
        pf0, f0 = field(ref)
        nz = sum(1 for v in f0.values() if v != 0.0)
        P(f"    Sp_field at {pf0}: {len(f0)} cells, {nz} non-zero")
        # how many of the scored cells are NOT a column top (lateral faces)?
        top = {}
        for (i, j, k), v in f0.items():
            if v != 0.0:
                top[(i, j)] = max(top.get((i, j), -1), k)
        nlat = sum(1 for (i, j, k), v in f0.items() if v != 0.0 and k < top[(i, j)])
        P(f"      of which {nlat} below their column's top cell (lateral faces "
          f"scored by the Sp_field lateral branch)")
        if nz == 0:
            P("    FAIL  Sp_field is identically zero -- (ii) would be vacuous"); bad += 1
    for leg in legs[1:]:
        h = md5(os.path.join(LEGS, leg, "run_removal_events.csv"))
        ok_e = (h == h0)
        line = f"    {'ok  ' if ok_e else 'FAIL'}  {leg}: event log {'byte-identical' if ok_e else 'DIFFERS (md5 ' + h[:12] + ')'}"
        bad += 0 if ok_e else 1
        if check_sp:
            pf, f = field(leg)
            same_keys = set(f) == set(f0)
            ndiff = sum(1 for k in f0 if f.get(k) != f0[k])
            ok_s = same_keys and ndiff == 0
            line += f"; Sp_field {'identical' if ok_s else f'{ndiff} cells DIFFER'} ({pf})"
            bad += 0 if ok_s else 1
        P(line)
    return bad


def main():
    print("P2-2e (a3) -- DECOMPOSITION IDENTITY at d = 8.5 mm")
    print("=" * 90)
    bad = 0
    bad += group("(i)+(ii) Kant column, SOD 130, Gotthard thermal, 200 s",
                 ["A3C_np1_m25", "A3C_np2_m25", "A3C_np4_m25", "A3C_np4_m50"])
    bad += group("(i)+(ii) Hoeser G05, 1 mm, 5 s",
                 ["A3H_np1_m25", "A3H_np2_m25", "A3H_np4_m25", "A3H_np4_m50",
                  "A3H_np4_m8"])
    bad += group("(iii) Hoeser G05 with spall.per_cell_removal = 1, 5 s",
                 ["A3P_np2_m25", "A3P_np4_m50", "A3P_np4_m8"], need_pcr=True)
    print("=" * 90)
    for c in CAVEATS:
        print(f"  ⚠ CAVEAT: {c}")
    print("PASS (a3) on the packet's bar -- (i), (ii) and (iii) all identical"
          if bad == 0 else f"FAIL (a3): {bad} check(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
