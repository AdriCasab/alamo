#!/usr/bin/env python3
"""P1g-2b ord 1 (a5): emit `WALL_FINDINGS.md` — only what is SETTLED about the
hole-wall coefficient, and explicitly what is not.

⚠⚠ THIS NOTE EXISTS TO RECORD A DEFERRAL AND TWO OPEN ITEMS, not to supply a
coefficient. P2-0 proposed replacing the typed 169.3 with a duct correlation; its
own ord 1 refuted the premise. What is settled is recorded; what is
convention-dependent is recorded AS a band with unknown sign; and the reason the
wall-jet law is deferred is the BRAKE, not the lever.

⚠ STANDING RULE 8. None of the wall-jet literature is in `validation/`, so every
external number here is cited as **a reviewer's estimate** from
`Claude_markdowns/2026-10-06a.md`, not as a backed value. The duct figures ARE
backed — they come from `studies/p2_0_annulus/annulus_h.py`, which is imported
here rather than retyped, so one implementation produces both notes.

⚠⚠ WHAT REV 1 GOT WRONG, AND WHY IT MATTERS. Rev 1 declared the friction
convention settled by comparing the wall-jet corner value (268) against
251 W/m²K read out of Leg B's own `h_expr`. But that 251 is **Martin multiplied by
a hand-set inverse-standoff decay** (`0.09/s`, s clamped at 0.09 inside Martin) —
i.e. a fitted model closure, not a correlation. Martin *proper* at r = 51 mm is
~529, which lands on the OTHER convention's ~536. **A physics estimate is
validated against the correlation, never against a fitted map.**
"""
import csv
import math
import os
import statistics as st
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
sys.path.insert(0, os.path.join(ROOT, "tests/MMWSpalling/studies/p2_0_annulus"))
CTRL2 = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g3_thickness/MEIER_CTRL2")

# The duct implementation already scored in P2-0's (a4). Imported, not retyped.
import annulus_h as duct                                    # noqa: E402

H_TYPED = 169.3          # surface_patch.side_face_h, the constant under question
R_GATE = 0.040           # side_face_h_r_min -- which faces take wall treatment
DZ, KTOP0 = 0.002, 299
T_CUT = 100.0            # the packet's void census starts after 100 s

# ---------------------------------------------------------------------------
# REVIEWER'S ESTIMATE (Claude_markdowns/2026-10-06a.md §7 item 1, the
# density-corrected central table). ⚠ NOT in validation/ -- rule 8.
# h [W/m²K] on the Bradshaw-Gee / Launder-Rodi c_f with the 1/2 in the
# definition. The OTHER convention (c_f on rho*U^2, no 1/2) DOUBLES every row.
# ---------------------------------------------------------------------------
WJ = ((0.0, 268.0),      # the corner, floor side (radial wall jet at r = R)
      (25.0, 190.0), (50.0, 164.0), (100.0, 132.0),
      (150.0, 113.0), (200.0, 100.0))
CONV = {"Launder-Rodi as used (with the 1/2)": 1.0,
        "c_f on rho*U^2 (no 1/2)": 2.0}

# Carried from the planner, NOT reproduced here: converting a coefficient ratio
# into a bore change needs a dh -> dOe transfer function this step does not have.
LEVER_MM = (-3.9, +12.8)
LEVER_BORE = (98.3, 115.0)
BORE_NOW = 102.2


def wj_h(z_mm, scale=1.0):
    """linear interpolation of the reviewer's profile, in mm above the floor."""
    if z_mm <= WJ[0][0]:
        return WJ[0][1] * scale
    if z_mm >= WJ[-1][0]:
        return WJ[-1][1] * scale
    for (a, va), (b, vb) in zip(WJ, WJ[1:]):
        if a <= z_mm <= b:
            return (va + (vb - va) * (z_mm - a) / (b - a)) * scale
    return WJ[-1][1] * scale


def void_heights():
    """height above the advancing floor for every lateral VOID, plus the
    final-state bore profile near the front.

    ⚠⚠ "VOID" means `n_voided >= 1`, NOT every regime-6 firing. Under
    `lateral_subcell=1` a firing banks its thickness and voids a cell only when
    the accumulator crosses dx, so firings outnumber voids. Counting firings here
    would overstate the census by ~35 % and would also mis-describe the mechanism
    the banking row exists to measure. Both counts are returned and both are
    reported.

    The floor at time t is the deepest removal reached up to t, so events are
    swept in time order and the running maximum is the floor.

    ⚠ The near-front bore is taken from the FINAL state -- the hole's width as a
    function of height above the LAST floor -- not as a running maximum over the
    whole run. A running max asks "was the hole ever this wide this close to the
    floor", which a wall that widens long after the floor has passed answers YES
    to; the scope question is about the FRESHLY cut hole at the front, which is
    narrow. The two measures differ by ~3.5x and only the final-state one speaks
    to scope.
    """
    rows = []
    with open(os.path.join(CTRL2, "FH_removal_events.csv")) as fh:
        for r in csv.DictReader(fh):
            try:
                t = float(r["time"])
                kt = int(float(r["k_top"]))
                rd = float(r.get("riser_depth") or 0.0)
                ci, cj = int(r["col_i"]), int(r["col_j"])
                reg = int(float(r["regime"]))
                nv = float(r.get("n_voided") or 0.0)
            except (ValueError, KeyError, TypeError):
                continue
            d = (KTOP0 - kt) * DZ * 1e3 + rd * 1e3
            rr = DZ * math.hypot(ci + 0.5, cj + 0.5) * 1e3
            rows.append((t, d, rr, reg, nv))
    rows.sort(key=lambda x: x[0])
    floor = 0.0
    hs, hs_fire = [], []
    depth_rmax = {}
    for t, d, rr, reg, nv in rows:
        floor = max(floor, d)
        if reg == 6 and t > T_CUT:
            hs_fire.append(floor - d)
            if nv >= 1.0:
                hs.append(floor - d)
        k = int(d // 5) * 5          # bore vs ABSOLUTE depth, for the final state
        if rr > depth_rmax.get(k, 0.0):
            depth_rmax[k] = rr
    # final-state profile: height above the LAST floor
    bore_by_h = {}
    for k, rr in depth_rmax.items():
        h = floor - k
        if h < 0:
            continue
        b = int(h // 5) * 5
        bore_by_h[b] = max(bore_by_h.get(b, 0.0), 2 * rr)
    return hs, hs_fire, bore_by_h, floor


def main():
    P = print
    hs, hs_fire, bore_by_h, floor = void_heights()
    hs_sorted = sorted(hs)
    q = lambda p: hs_sorted[min(len(hs_sorted) - 1, int(p * len(hs_sorted)))]
    near = sorted((b, v) for b, v in bore_by_h.items() if b <= 25)

    P("# P1g-2b `(a5)` — the hole-wall coefficient: what is settled, what is not")
    P("")
    P("Generated by `wall_findings.py`. Offline; no source edited, nothing run.")
    P("")
    P("## ⚠⚠ The deferral, and its reason")
    P("")
    P("**Wall heating is not the width-setting feedback, so the wall-jet law is")
    P("deferred.** The reason is the **brake**, not the lever:")
    P("")
    P("- P2-0's physics case was a *feedback*: a duct coefficient rises steeply as")
    P("  the hole narrows, so the bore becomes an equilibrium instead of a runaway.")
    P("- In a **radial wall jet** that feedback nearly vanishes. Momentum is")
    P("  conserved over an expanding annulus, so `v_m ∝ 1/R` while `Re` stays")
    P("  roughly constant, hence `c_f` barely moves and `h ∝ 1/R` — a **linear**")
    P("  brake, against the duct's steep one.")
    P("- ⚠ The brake is **convention-independent**: changing the `c_f` definition is")
    P("  a constant multiplier, so every *ratio* survives it. This is the one")
    P("  wall-jet quantity that can be recorded as settled.")
    P("")
    P("Measured on one diameter pair at a time, since a ratio needs a stated pair:")
    P("")
    P("| diameter pair | duct brake (backed) | wall-jet brake `1/R` | ratio |")
    P("|---|---|---|---|")
    mdot = 4.0 * 0.0037826388888888888 / 0.3364     # from annulus_h's (a1)
    for d0, d1 in ((84.0, 121.8), (102.0, 121.8), (84.0, 102.0)):
        h0, _, _, _ = duct.annulus(d0 * 1e-3 / 2, mdot, 1514.0, "fit")
        h1, _, _, _ = duct.annulus(d1 * 1e-3 / 2, mdot, 1514.0, "fit")
        P("| Ø%.0f → Ø%.1f | **%.2f×** | **%.2f×** | %.1f |"
          % (d0, d1, h0 / h1, d1 / d0, (h0 / h1) / (d1 / d0)))
    P("")
    P("⚠ **The packet's \"1.2× against the duct's 4.04×\" compares two different")
    P("diameter pairs** — 4.04 is Ø84→Ø121.8 and 1.2 is Ø102→Ø121.8. On a common")
    P("pair the gap is 4.04 vs 1.45, or 1.75 vs 1.19. The conclusion is unchanged")
    P("either way and is **qualitative**: the duct brake is 2.8–3× the wall jet's.")
    P("")
    P("## ✅ Settled: the wall is a jet, not a duct")
    P("")
    P("The faces that actually take wall treatment sit inside an impinging/radial")
    P("wall jet, not in developed annular duct flow. The duct value is the")
    P("**far-field limit misapplied**:")
    P("")
    hd, Re, Dh, _ = duct.annulus(BORE_NOW * 1e-3 / 2, mdot, 1514.0, "fit")
    P("- duct Gnielinski at Ø%.1f mm reproduces **%.1f W/m²K** (backed, from"
      % (BORE_NOW, hd))
    P("  `annulus_h.py`), and the reviewer's note recovers it as the `z > z_fill`")
    P("  limit rather than contradicting it — it is simply not where the faces are.")
    P("- hydraulic diameter there is **%.1f mm**, so the live wall band sits a few"
      % (1e3 * Dh))
    P("  `D_h` below the nozzle, inside the developing jet.")
    P("")
    P("## ⚠⚠ NOT settled: the friction convention, and therefore every absolute h")
    P("")
    P("The open item is **Launder & Rodi 1981's definition of `c_f`**. Taken as")
    P("`c_f = τ_w/(½ρU_m²) = 0.0315 Re^−0.182` the corner is ~268 W/m²K; if the")
    P("original defines `c_f` on `ρU_m²` with no ½, **every value doubles**. That")
    P("single definition is the largest lever in the whole estimate.")
    P("")
    P("| z above floor [mm] | %s | %s |"
      % tuple(CONV.keys()))
    P("|---|---|---|")
    for z, _ in WJ:
        P("| %s | %.0f | %.0f |"
          % ("corner (r = R)" if z == 0 else "%.0f" % z,
             wj_h(z, 1.0), wj_h(z, 2.0)))
    P("")
    P("**The discriminating cross-check — and rev 1 got it backwards:**")
    P("")
    P("| comparand | value | what it is |")
    P("|---|---|---|")
    P("| Martin **proper**, r = 51 mm | **~529** | the correlation itself |")
    P("| Leg B `h_expr` at r = 51 mm | **~251** | ⚠ Martin **×** a hand-set `0.09/s` decay — a FITTED MAP |")
    P("| wall jet, with the ½ | **268** | matches the fitted map |")
    P("| wall jet, no ½ | **~536** | matches Martin proper to **1.3 %** |")
    P("")
    P("⇒ the corner is **268 or ~536**, and it stays open until Launder & Rodi")
    P("1981 is read. ⚠ Rev 1 declared it settled by matching 268 against the 251,")
    P("i.e. **against a model closure containing a hand-set decay**. A physics")
    P("estimate is validated against the correlation, never against a fitted map.")
    P("")
    P("## ⚠⚠ NOT settled: the lever's SIGN")
    P("")
    P("Where the model actually cuts the wall, weighted by its own lateral voids:")
    P("")
    P("- **%d** lateral **voids** after %.0f s (from **%d** regime-6 firings — ⚠ the"
      % (len(hs), T_CUT, len(hs_fire)))
    P("  two differ because banking voids a cell only when the accumulator crosses")
    P("  `dx`, so **%.0f %%** of firings bank without voiding);"
      % (100.0 * (1 - len(hs) / len(hs_fire))))
    P("  height above the advancing floor p05 **%.0f** / p50 **%.0f** / p95 **%.0f** mm;"
      % (q(0.05), q(0.50), q(0.95)))
    P("  **%.1f %%** in 50–100 mm."
      % (100.0 * sum(1 for x in hs if 50 <= x < 100) / len(hs)))
    P("- **%d (%.1f %%)** within 25 mm of the floor."
      % (sum(1 for x in hs if x < 25), 100.0 * sum(1 for x in hs if x < 25) / len(hs)))
    P("")
    P("⚠⚠ **That zero is SCOPE, not physics.** In the final state (floor at")
    P("**%.0f mm**) the freshly cut hole just above the front is:" % floor)
    P("")
    P("| height above the final floor [mm] | bore Ø [mm] | radius [mm] | inside the %.0f mm gate? |"
      % (1e3 * R_GATE))
    P("|---|---|---|---|")
    for b, v in near:
        P("| %d–%d | **%.1f** | %.1f | %s |"
          % (b, b + 5, v, v / 2, "**yes**" if v / 2 < 1e3 * R_GATE else "no"))
    P("")
    P("⇒ no face within ~25 mm of the front is in lateral scope at all, so the")
    P("absence of firings there cannot be a statement about heat transfer.")
    P("⚠ Measured on the **final state**, not as a running maximum over the run: a")
    P("running max answers \"was the hole ever this wide this close to the floor\",")
    P("which a wall that widens long after the floor has passed answers yes to")
    P("(it gives ~90 mm at 0–10 mm). Only the final-state profile speaks to scope.")
    P("")
    P("The %.0f mm median is mostly the depth the hole needs to pass" % q(0.50))
    P("r = 40 mm, plus ~10 mm of thermal lag. ⇒ a wall-jet law changes nothing at")
    P("the corner **unless the 40 mm gate is also lifted**, and lifting it puts the")
    P("law onto the pad-ring risers that feed the floor rate (D2w/D2x/D2o-0) —")
    P("a rate-containment decision, not a side effect.")
    P("")
    P("Void-weighted coefficient against the typed %.1f, both conventions:" % H_TYPED)
    P("")
    P("| convention | void-weighted mean h | ratio to %.1f |" % H_TYPED)
    P("|---|---|---|")
    for lbl, sc in CONV.items():
        m = st.mean(wj_h(x, sc) for x in hs)
        P("| %s | **%.0f** | **%.2f×** |" % (lbl, m, m / H_TYPED))
    P("")
    P("⇒ the law is **below** the constant on one convention and **above** it on")
    P("the other. The planner's bore translation of that band is **%+.1f to %+.1f mm**"
      % LEVER_MM)
    P("(bore **%.1f or %.1f mm** against today's %.1f) — ⚠ **carried from the**"
      % (LEVER_BORE[0], LEVER_BORE[1], BORE_NOW))
    P("**planner, not reproduced here**, because converting a coefficient ratio into")
    P("a bore change needs a `dh → dØ` transfer function this step does not have.")
    P("")
    P("## What this note backs — and the four things not to quote")
    P("")
    P("✅ **Backs:** the wall is a jet not a duct; the duct **%.1f** reproduces as"
      % hd)
    P("the far-field limit; the **brake** is linear (`1/R`) against the duct's steep")
    P("one, by a factor of **2.8–3× on a common diameter pair**, and that ratio is")
    P("convention-independent; ⇒ **wall heating is not the width lever**, which is")
    P("why the implementation is deferred.")
    P("")
    P("⚠ **Do NOT quote:**")
    P("")
    P("1. **\"169.3 is right in the mean\"** — true on one convention only.")
    P("2. **\"−3.9 mm\"** as the lever — the band spans a sign change.")
    P("3. any absolute wall-jet `h` as settled — all of them double under the")
    P("   other `c_f` definition.")
    P("4. **\"the coefficient is why nothing fires at the corner\"** — that is the")
    P("   40 mm **gate**, measured above.")
    P("")
    P("**Open item, named:** Launder & Rodi 1981's `c_f` definition. ⚠ Every")
    P("external number above is **a reviewer's estimate**")
    P("(`Claude_markdowns/2026-10-06a.md` §2–§3, §7) — none is in `validation/`,")
    P("standing rule 8. The duct figures are backed by `annulus_h.py`, imported here.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
