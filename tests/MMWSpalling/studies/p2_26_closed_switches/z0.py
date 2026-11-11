#!/usr/bin/env python3
"""P2-26 z0 (offline gate) -- everything step zero must establish before the first edit.

  1. the census items (CODE_CENSUS_2026-11-05 §2 B1, B2, B3, C2, C3) re-resolved to TODAY's lines: item, symbol,
     file:line, action (every item found, or struck with a reason)
  2. the deck-edit list (decks.py): every setter of the twelve keys, DEFAULT/LIVE, EDIT/RECORD/ARGV; every LIVE
     EDIT line is a C2/C3 key; zero B3 deck setters
  3. the gate-header grep: `len_fb_h0`, `sn_recession_superseded`, `ell_q` in the gate outputs on disk (the p2_25
     sweep: hash.sh over them equals witness/p2_25.md5)
  4. the all-cell melt read: Temp and Lambda_L per-box maxima from EVERY plotfile header (Level_0/Cell_H) of the six
     melt-aware gate tests; T_top over every firing in their event logs; the IEEE facts MeltSkeleton(0) == 1.0 and
     1.0 * x == x
  5. the flake-cap scan: every deck's effective spall.flake_h_max against its effective spall.h_spall_max (FILE =
     includes resolved); every spall-off deck that sets h_spall_max (B1 moves that parse above the spall_enabled
     return, so such a deck would START reading it); the L1 probe's discrimination (the check reads the member
     initialiser 0.010 today)
  6. the pre-registered movers: the 233 paths of beam_void_closure + dev2d in witness/p2_18.md5 (MOVERS_PRE.txt)
  7. units of every deleted key (L9); src_pre_1/ complete and equal to src/ (copied before any edit)
Usage: z0.py > Z0.txt   (exit 1 on a failed bar)
"""
import glob
import hashlib
import importlib.util
import os
import re
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
SRC = {"M": "src/Integrator/MMWSpalling.H", "R": "src/Integrator/MMWSpalling/Removal.H",
       "F": "src/Integrator/MMWSpalling/FaceRemoval.H", "S": "src/Integrator/MMWSpalling/Shared.H"}
_spec = importlib.util.spec_from_file_location("decks", os.path.join(HERE, "decks.py"))
decks = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(decks)
BAD = []


def need(ok, msg):
    print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
    if not ok:
        BAD.append(msg)


def text(k):
    return open(os.path.join(ROOT, SRC[k])).read().split("\n")


# ------------------------------------------------------------------ 1. items ----
# (census item, file key, regex, action). Regexes name a SYMBOL or a distinctive phrase, never a line number.
ITEMS = [
    ("B1 flake_h_max check (reads h_spall_max)", "M", r'value\.sp_flake_h_max >= value\.spall_h_spall_max',
     "keep; the h_spall_max parse moves above it"),
    ("B1 h_spall_max parse", "M", r'pp\.query_default\("spall\.h_spall_max"', "move above the flake_h_max check"),
    ("B1 the spall-off return", "M", r'if \(!value\.spall_enabled\) return;', "unchanged (the parse moves above it)"),
    ("B1 side_power_interval parse (enthalpy branch)", "M", r'pp\.query_default\("side_power_interval"',
     "move out of the jet_closure = enthalpy branch"),
    ("B1 fl_wtop (written, never read)", "M", r'\bfl_wtop\b', "delete"),
    ("B1 sp_surface_T_mode + SurfaceTMode", "M", r'\bsp_surface_T_mode\b|\bSurfaceTMode\b',
     "delete (keep sp_surface_T_recon/_apply and the surfT parse)"),
    ("B2 'P1b: the acceptance test ...' (unfinished)", "M", r'P1b: the acceptance test for the support rule', "delete"),
    ("B2 Advance header above SurfaceRecon", "M", r'Override HeatConduction::Advance with the enthalpy method',
     "move above Advance (and name the microstructure branch)"),
    ("B2 C1 pin header above BuildSurfaceTColumns", "M", r'C1: per-column pin used by this step',
     "move above BuildPinEffective"),
    ("B2 SurfaceRobinT header above SharedSurfaceRecon", "M", r'P0b deliverable 1\.1: the per-column surface temperature',
     "move above SurfaceRobinT"),
    ("B2 orphan 'P1 -- THE REFERENCE TEMPERATURE'", "M", r'P1 -- THE REFERENCE TEMPERATURE, READ FROM THE EXPONENTIAL',
     "delete"),
    ("B2 SurfaceCellFlux header (stranded)", "M", r'S1: per-cell surface heat flux \[W/m\^2\]',
     "move above SurfaceCellFlux"),
    ("B2 RiserWeight tombstone", "M", r'RiserWeight \(D2x, side_face_area_exact\) DELETED', "delete"),
    ("B2 'P1f deliverable 0' lateral-criterion rationale", "M", r'P1f deliverable 0 -- THE FACE TEMPERATURE, OUT',
     "rewrite: T_face_out is read by FrFaceFlux; the lateral criterion is gone"),
    ("B2 ColumnTopSolid header above UpdatePropertyFields", "M", r'E1: global topmost solid \(removed <= 0\.5\)',
     "move above ColumnTopSolid"),
    ("B2 'P1 -- THE GRADIENT-AWARE DRIVING STRESS' (no parse)", "M", r'P1 -- THE GRADIENT-AWARE DRIVING STRESS',
     "delete"),
    ("B2 GatherColumnSurface header 'phi >= 0'", "M", r'Fill k_top_col\[ncols\] with the topmost live',
     "fix: removed-only since D2q-1"),
    ("B2 UpdateSpAfterMechanics header above ColumnSurfaceCos", "M", r'Sp = K_I/K_Ic\(T\) surface diagnostic for',
     "move above UpdateSpAfterMechanics"),
    ("B2 LawSigma 'ParallelFor under [=, this]'", "M", r'amrex::ParallelFor` under a', "fix: LoopOnCpu, [&]"),
    ("B2 duplicated 'Step 5: Voronoi' line", "M", r'Step 5: Voronoi mineral microstructure', "delete one of two"),
    ("B2 P1/P1b support-rule member text", "M", r'P1\. Default `absolute` \+ `flake` keeps|P1b diagnostics, registered only',
     "delete"),
    ("B2 dead local pc", "M", r'int pc = static_cast<int>\(ph\(i, j, k\)\)', "delete with its (void)"),
    ("B2 dead local p_n", "M", r'\bp_n\b', "delete"),
    ("B2 (void) q_rob/q_los", "M", r'\(void\)\s*q_rob;\s*\(void\)\s*q_los;', "delete (and the unused bindings)"),
    ("B2 (void) dt of BuildSurfaceTColumns", "M", r'\(void\) dt;', "delete the parameter"),
    ("B2 (void) z_n of FaceListPath", "M", r'\(void\) z_n;', "delete the parameter"),
    ("B2 redundant (void) pair in TimeStepComplete", "M", r'^\s*\(void\)\s*a_time;$', "delete (both used on the lines above)"),
    ("B3 pad_shield_open", "M", r'\bpad_shield_open\b', "delete parse, member, uses, ph_abort"),
    ("B3 CheckPadRing", "M", r'\bCheckPadRing\b', "delete"),
    ("B3 jet_above_h / jet_T_above", "M", r'\bjet_above_h\b|\bjet_T_above\b', "delete parse, members, uses, ph_abort"),
    ("B3 rec_fixed + jet_T_rec_fixed", "M", r'\bjet_rec_fixed\b|\bjet_T_rec_fixed\b|"rec_fixed"',
     "delete; jet_T_ent_mode keeps fixed / exhaust"),
    ("B3 flake_length=flux / flake_fraction=weibull", "M",
     r'\bsp_flake_len_flux\b|\bsp_flake_frac_weibull\b|"spall\.flake_length"|"spall\.flake_fraction"',
     "delete parse, members, ell_flux branch, fr_abort clause"),
    ("B3 len_fb_h0/q0 thermo columns", "M", r'\bsp_len_fb_[hq]0\b', "delete (registered only under flux)"),
    ("B3 ReconFromX ell_flux/frac_over params", "M", r'\bell_flux\b|\bfrac_over\b', "delete (default args only)"),
    ("B3 the ell_q event column (flux only) + rc_ellq", "R", r'\brc_ellq\b|sp_flake_len_flux|sp_flake_frac_weibull',
     "delete with the switch"),
    ("B3 WeibullSizeFactor + size_volume + V_ref", "M", r'\bWeibullSizeFactor\b|\bsp_size_flux\b|\bweibull_V_ref\b',
     "delete; every caller returns to the bare flaw"),
    ("B3 WeibullSizeFactor caller (thickness scan)", "R", r'\bWeibullSizeFactor\b|\bsp_size_flux\b', "delete"),
    ("C2 surface_normal parse/member/registration", "M",
     r'"spall\.surface_normal"|\bspall_surface_normal\b|\bsn_recession_superseded\b', "delete"),
    ("C2 ColumnSurfaceCos + SP_SN_COS_MIN", "M", r'\bColumnSurfaceCos\b|\bSP_SN_COS_MIN\b', "delete"),
    ("C2 site-A /cos sampling (sn_cos_p, icth)", "M", r'\bsn_cos_p\b|\bicth\b|\bsn_cos_col\b', "delete (icth = 1)"),
    ("C2 column-path uses", "R", r'\bspall_surface_normal\b|\bsn_recession_superseded\b|\bsn_cos_p\b|\bsn_cos_col\b|'
     r'\bsn_active\b|\bColumnSurfaceCos\b', "delete"),
    ("C2 DepthScanH sn_cth/sn_icth", "S", r'\bsn_cth\b|\bsn_icth\b', "delete the parameters (1.0 everywhere)"),
    ("C2 DepthScanH call sites", "R", r'\bsn_cth\b|\bsn_icth\b', "drop the arguments"),
    ("C2 the face list's literal 1.0, 1.0", "F", r'DepthScanH\(nprof, 0, -\(L\.Lg \+ 1\), dx, 1\.0, 1\.0',
     "drop the two literals"),
    ("C2 flake_coherence_length", "M", r'"spall\.flake_coherence_length"|\bspall_flake_coh_len\b', "delete"),
    ("C2 the coherence cap", "R", r'\bspall_flake_coh_len\b|\bcoh_cells\b', "delete"),
    ("C3 melt_aware + lambda_crit + MeltSkeleton", "M",
     r'"spall\.melt_aware"|"melt\.lambda_crit"|\bmelt_aware\b|\bmelt_lambda_crit\b|\bMeltSkeleton\b', "delete"),
    ("C3 the Lambda_L gather, LL_at_didx, f_skel", "R", r'\bLL_col_prof\b|\bLL_col_p\b|\bLLcol_i\b|\bLL_at_didx\b|'
     r'\bf_skel\b|\bmelt_aware\b|\bMeltSkeleton\b', "delete"),
    ("C3 ScanSigmaT's f_skel parameter", "S", r'\bf_skel\b', "delete the parameter"),
    ("C3 the face list's literal f_skel 1.0", "F", r'ScanSigmaT\(scan_g, [^;]*?, 1\.0,', "drop the literal"),
]


def items():
    print("1. CENSUS ITEMS, re-resolved to today's lines (M MMWSpalling.H, R Removal.H, F FaceRemoval.H, S Shared.H)")
    for item, fk, rx, action in ITEMS:
        L = text(fk)
        hits = [i for i, ln in enumerate(L, 1) if re.search(rx, ln)]
        where = "%s:%s" % (fk, ",".join(str(h) for h in hits[:14]) + (" ... (%d)" % len(hits) if len(hits) > 14 else ""))
        print("  %-4s %-54s %-44s -> %s" % ("ok" if hits else "GONE", item, where if hits else "-", action))
        need(bool(hits), "%s resolves to a present symbol" % item) if not hits else None
    print("  (every item found)" if not BAD else "")


# ------------------------------------------------------------------ 2. decks ----
def deck_list():
    print("\n2. DECK-EDIT LIST (decks.py; EDIT = kept decks, RECORD = study/retired/pre decks, ARGV = live test scripts)")
    rows = decks.setters()
    for r in rows:
        print("  %-6s %-66s %4d  %-30s %-16s %-7s %s" % r)
    edit = [r for r in rows if r[0] == "EDIT"]
    files = sorted({r[1] for r in rows if r[0] in ("EDIT", "MOVED")})
    print("  distinct EDIT/MOVED files: %d" % len(files))
    need(all(r[6] in ("C2", "C3") for r in edit if r[5] == "LIVE"),
         "every LIVE line on a kept deck is a C2/C3 key (L10)")
    need(not [r for r in edit if r[6] == "B3"], "zero B3 setters on any kept deck (P1: stage 1 moves no deck)")
    argv_b3 = sorted({(r[1], r[3]) for r in rows if r[0] == "ARGV" and r[6] == "B3"})
    print("  B3 setters in live test scripts (the legs this packet drops or restates): %s" % argv_b3)
    gate = ("unit/beam_void_closure/input", "validation/meier/sp_meier_pilot/input_2d_dev",
            "unit/robin_pinned/input_pinned", "unit/robin_jet/input_jet", "unit/jet_enthalpy/input_quarter",
            "unit/jet_enthalpy/input_strip", "unit/robin_feet/input_feet", "unit/robin_feet/input_column",
            "unit/scalar_flaw/input", "studies/s1_surface_resolution/input")
    live_gate = sorted({(r[1], r[3]) for r in edit if r[5] == "LIVE" and r[1] in gate})
    print("  LIVE lines on GATE decks: %s" % live_gate)
    need({k for _, k in live_gate if k in ("spall.surface_normal", "spall.flake_coherence_length")} and
         {f for f, k in live_gate if k in ("spall.surface_normal", "spall.flake_coherence_length")}
         == {"unit/beam_void_closure/input", "validation/meier/sp_meier_pilot/input_2d_dev"},
         "the C2 LIVE gate decks are exactly beam_void_closure and dev2d (the two pre-registered movers)")
    return rows


# --------------------------------------------------------- 3. gate headers ----
GATE_DIRS = ["unit/robin_face/output", "unit/robin_pinned/output", "unit/robin_jet/output", "unit/jet_enthalpy/output",
             "unit/beam_void_closure/output", "unit/scalar_flaw/output", "unit/spall_event/output",
             "unit/spall_event/output_no_event", "unit/robin_feet/output", "validation/meier/sp_meier_pilot/output/dev2d",
             "studies/s1_surface_resolution/output_smoke"]


def headers():
    print("\n3. GATE HEADERS (the outputs on disk; hash.sh over them == witness/p2_25.md5, checked below)")
    tmp = os.path.join(HERE, ".z0_now.md5")
    subprocess.run(["bash", os.path.join(W, "hash.sh"), tmp], check=True)
    same = open(tmp).read() == open(os.path.join(W, "p2_25.md5")).read()
    os.remove(tmp)
    need(same, "the gate outputs on disk hash to witness/p2_25.md5 (binary c052c348's sweep)")
    hit = {"len_fb_h0": [], "sn_recession_superseded": [], "ell_q": []}
    nth = nev = 0
    for d in GATE_DIRS:
        for p in glob.glob(os.path.join(T, d, "**", "thermo.dat"), recursive=True):
            if ".old" in p:
                continue
            nth += 1
            h = open(p).readline().split()
            for k in ("len_fb_h0", "sn_recession_superseded"):
                if k in h:
                    hit[k].append(os.path.relpath(p, T))
        for p in glob.glob(os.path.join(T, d, "**", "*_events.csv"), recursive=True):
            nev += 1
            if "ell_q" in open(p).readline().strip().split(","):
                hit["ell_q"].append(os.path.relpath(p, T))
    print("  %d thermo.dat headers, %d event-log headers read" % (nth, nev))
    for k, v in hit.items():
        print("  %-24s in %d: %s" % (k, len(v), ", ".join(v)))
    need(not hit["len_fb_h0"], "no gate header carries len_fb_h0 (L4)")
    need(not hit["ell_q"], "no gate event log carries the flux-only ell_q column")
    tops = sorted({p.split("/output")[0] for p in hit["sn_recession_superseded"]})
    # ⚠ dev2d writes NO thermo.dat and no event log (input_2d_dev sets no amr.thermo.*; its 64 gate files are 16
    # plotfiles x 4 Cell_D fabs), so the column can leave only beam_void_closure's header. The packet's b1 clause
    # "both legs' thermo.dat" is restated to "beam_void_closure's thermo.dat; dev2d moves through Cell_D only".
    d2 = os.path.join(T, "validation/meier/sp_meier_pilot/output/dev2d")
    print("  dev2d writes thermo.dat: %s; event logs: %d" % (os.path.exists(os.path.join(d2, "thermo.dat")),
          len(glob.glob(d2 + "*_events.csv"))))
    need(tops == ["unit/beam_void_closure"] and not os.path.exists(os.path.join(d2, "thermo.dat")),
         "sn_recession_superseded is in exactly one gate header, beam_void_closure's; dev2d (the other surface_normal "
         "gate deck) writes no thermo.dat (L5, restated)")


# -------------------------------------------------------------- 4. melt read ----
MELT_TESTS = ["unit/robin_pinned/output", "unit/robin_jet/output", "unit/jet_enthalpy/output", "unit/robin_feet/output",
              "unit/scalar_flaw/output", "studies/s1_surface_resolution/output_smoke"]
SOLIDUS = 1473.0 - 25.0     # phase0.T_m 1473.0, fus_width 50.0 => the mushy band starts at 1448 K


def cell_h_max(pf, comps):
    """per-component max over every box, from Level_0/Cell_H's min/max block"""
    hdr = open(os.path.join(pf, "Header")).read().split("\n")
    n = int(hdr[1])
    names = hdr[2:2 + n]
    idx = {c: names.index(c) for c in comps if c in names}
    L = open(os.path.join(pf, "Level_0", "Cell_H")).read().split("\n")
    blocks = [i for i, ln in enumerate(L) if re.fullmatch(r"\d+,\d+", ln.strip())]
    if len(blocks) < 2:
        return None
    nb = int(L[blocks[1]].split(",")[0])
    mx = [list(map(float, L[blocks[1] + 1 + b].rstrip(",").split(","))) for b in range(nb)]
    return {c: max(row[i] for row in mx) for c, i in idx.items()}


def melt():
    print("\n4. THE ALL-CELL MELT READ (six melt-aware gate tests; every plotfile's per-box maxima)")
    print("   solidus = phase0.T_m 1473.0 - fus_width/2 (50.0) = %.0f K" % SOLIDUS)
    allT = allL = -1.0
    tot = 0
    for d in MELT_TESTS:
        pfs = [p for p in glob.glob(os.path.join(T, d, "**", "*cell"), recursive=True)
               if os.path.isdir(p) and ".old" not in p and os.path.exists(os.path.join(p, "Level_0", "Cell_H"))]
        mT = mL = -1.0
        for pf in pfs:
            r = cell_h_max(pf, ("Temp", "Lambda_L"))
            mT, mL = max(mT, r["Temp"]), max(mL, r["Lambda_L"])
        ttop = -1.0
        for p in glob.glob(os.path.join(T, d, "**", "*_removal_events.csv"), recursive=True):
            L = open(p).read().split("\n")
            h = L[0].split(",")
            c = next((h.index(x) for x in ("T_top", "Ttop", "T_surf") if x in h), None)
            if c is not None:
                for ln in L[1:]:
                    f = ln.split(",")
                    if len(f) > c:
                        try:
                            ttop = max(ttop, float(f[c]))
                        except ValueError:
                            pass
        tot += len(pfs)
        allT, allL = max(allT, mT), max(allL, mL)
        print("  %-44s %4d plotfiles  max Temp %8.2f K  max Lambda_L %.3g  max event T_top %s"
              % (d, len(pfs), mT, mL, ("%.2f K" % ttop) if ttop > 0 else "-"))
    need(allT < SOLIDUS and allL == 0.0,
         "hottest cell %.2f K < solidus %.0f K and Lambda_L == 0 in every cell of every plotfile (%d) -- L7"
         % (allT, SOLIDUS, tot))
    one = 1.0 - 0.0 / 0.4
    probe = [655.4 * 1.23e6, -3.7e-9, 1.0e300, 5e-324]
    need(one == 1.0 and struct.pack("<d", one) == struct.pack("<d", 1.0) and all(one * x == x for x in probe),
         "IEEE: MeltSkeleton(0) = 1 - 0/0.4 == 1.0 bit-exactly, and 1.0 * x == x (the stress is unchanged)")


# ------------------------------------------------------------ 5. flake cap ----
def resolved(path, seen=None):
    """{key: value} with FILE = includes resolved (an include's keys first, the deck's own after -- the order
    ParmParse reads them; a later duplicate wins in a query_default)"""
    seen = seen or set()
    out = {}
    for ln in open(path, errors="replace").read().split("\n"):
        m = decks.a2.LINE.match(ln)
        if not m:
            continue
        k, v = m.group(1), m.group(2).strip()
        if k == "FILE":
            inc = os.path.join(ROOT, v)
            if os.path.exists(inc) and inc not in seen:
                seen.add(inc)
                out.update(resolved(inc, seen))
        else:
            out[k] = v
    return out


def flake_cap():
    print("\n5. THE FLAKE-CAP SCAN (B1/L1): effective spall.flake_h_max vs spall.h_spall_max, includes resolved")
    files = sorted(set(decks.a2.deck_files()) | {os.path.join(T, "studies", "relief_diag", "U", "input_cell")})
    bad_cap, spalloff = [], []
    for p in files:
        kv = resolved(p)
        rel = os.path.relpath(p, T)
        fh = float(kv.get("spall.flake_h_max", "0"))
        hm = float(kv.get("spall.h_spall_max", "0.010"))
        on = kv.get("spall.enabled", "0").strip() in ("1",)
        if fh > 0:
            print("  %-62s flake_h_max %-8g h_spall_max %-6g spall.enabled %d  %s"
                  % (rel, fh, hm, on, "ok" if fh < hm else "VIOLATES"))
            if fh >= hm:
                bad_cap.append(rel)
        if not on and "spall.h_spall_max" in kv and float(kv["spall.h_spall_max"]) != 0.010:
            spalloff.append((rel, kv["spall.h_spall_max"]))
    print("  test argv overrides: face_heat BASE flake_h_max=0.0015 on validation/hoeser/input_lab (h_spall_max %s via "
          "its includes); surface_T_recon's refusal leg flake_h_max=0.015 (h_spall_max 0.015: aborts before and after)"
          % resolved(os.path.join(T, "validation", "hoeser", "input_lab")).get("spall.h_spall_max", "0.010"))
    need(not bad_cap, "every live flake cap is below its deck's scan depth (the fix changes no live deck)")
    print("  spall-off decks that set a non-default h_spall_max (would START reading it): %s" % (spalloff or "none"))
    need(not spalloff, "no spall-off deck sets h_spall_max away from 0.010 (moving the parse above the return is inert)")
    # L1 probe discrimination: today the check reads the member initialiser 0.010
    init = 0.010
    for fh, hm, what in ((0.02, 0.015, "the packet's probe"), (0.012, 0.015, "discriminating (old aborts, new runs)"),
                         (0.008, 0.005, "discriminating (old runs, new aborts)")):
        old, new = fh >= init, fh >= hm
        print("  L1 probe flake_h_max=%g h_spall_max=%g: before %s, after %s  -- %s%s"
              % (fh, hm, "ABORTS" if old else "runs", "ABORTS" if new else "runs", what,
                 "" if old != new else "  ** NOT DISCRIMINATING"))


# ---------------------------------------------------------------- 6. movers ----
def movers():
    print("\n6. PRE-REGISTERED MOVERS (stage 2): every gate path of beam_void_closure and dev2d")
    L = [ln.split(None, 1)[1] for ln in open(os.path.join(W, "p2_18.md5")).read().split("\n") if ln.strip()]
    mv = [p for p in L if p.startswith("unit/beam_void_closure/output") or "/output/dev2d" in p]
    t0 = [p for p in mv if re.search(r"/0+(?:cell|node)/Level_0/", p)]
    bvc = sum(1 for p in mv if "beam_void_closure" in p)
    open(os.path.join(HERE, "MOVERS_PRE.txt"), "w").write(
        "P2-26 pre-registered stage-2 movers (written by z0.py BEFORE any edit): the 233 gate paths of "
        "beam_void_closure + dev2d. A differing path outside this list fails b1. t = 0 plotfiles (cannot move "
        "through C2, which acts at a firing) are marked T0.\n" +
        "\n".join(("T0 " if p in t0 else "   ") + p for p in mv) + "\n")
    print("  beam_void_closure %d + dev2d %d = %d paths (t = 0 plotfiles among them: %d); MOVERS_PRE.txt"
          % (bvc, len(mv) - bvc, len(mv), len(t0)))
    need(len(mv) == 233 and bvc == 169, "the mover list is exactly 169 + 64 = 233 paths")
    print("  C3 movers pre-registered: NONE (section 4)")


# ------------------------------------------------------------- 7. units, pre ----
def units_pre():
    print("\n7. UNITS OF THE DELETED KEYS (L9) and the pre-images")
    for k, (it, dflt, u) in decks.KEYS.items():
        print("  %-32s %-3s default %-8s [%s]" % (k, it, dflt, u))
    print("  spall.flake_h_max [m], spall.h_spall_max [m] (B1, moved); surface_patch.side_power_interval [s] (B1, moved)")
    pre = os.path.join(HERE, "src_pre_1")
    ok = os.path.isdir(pre)
    rows = []
    for k, rel in sorted(SRC.items()) + [("B", "src/Numeric/MMWBeam.H")]:
        a = os.path.join(ROOT, rel)
        b = os.path.join(pre, rel)
        ha = hashlib.md5(open(a, "rb").read()).hexdigest()
        hb = hashlib.md5(open(b, "rb").read()).hexdigest() if os.path.exists(b) else "MISSING"
        rows.append(ha == hb)
        print("  src_pre_1/%-44s %s" % (rel, "== src/" if ha == hb else ("DIFFERS from src/ (%s)" % hb)))
    need(ok and all(rows), "src_pre_1/ holds every header the stages touch, equal to src/ (copied before any edit)")
    bp = os.path.join(HERE, "BINARY_pre.txt")
    print("  BINARY_pre.txt: %s" % (open(bp).read().strip() if os.path.exists(bp) else "MISSING"))
    need(os.path.exists(bp), "BINARY_pre.txt written")


def main():
    print("P2-26 z0 -- step zero, offline\n" + "=" * 100)
    items()
    deck_list()
    headers()
    melt()
    flake_cap()
    movers()
    units_pre()
    print("=" * 100)
    print("z0 %s%s" % ("PASS" if not BAD else "FAIL", "" if not BAD else " (%d): %s" % (len(BAD), "; ".join(BAD))))
    return 0 if not BAD else 1


if __name__ == "__main__":
    sys.exit(main())
