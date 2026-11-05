#!/usr/bin/env python3
"""P2-2f — the price, the frozen floor, and (a3)'s cross-check. ONE module, so
`run_ctrl8.py` and the standalone (a3) gate cannot drift apart.

⚠⚠ WHAT REV 1 GOT WRONG, AND WHY EACH FIX IS WHAT IT IS.

(1) THE QUARTER-DOMAIN FACTOR WAS APPLIED TO ONE SIDE OF A RATIO.
    `run_ctrl8.py:184` multiplied the volume by 4 while `:196` divided the
    still-quarter-domain `ledger_E_flake` by it, so the price came out 4x low.
    ⇒ A PRICE IS A RATIO: SCALE BOTH SIDES OR NEITHER. Here: NEITHER. Both the
    event log and `ledger_E_flake` are one quadrant, so their ratio is already
    the full-domain specific energy and the x4 cancels exactly. Standing rule 3
    ("x4 for any Meier VOLUME") still governs a volume REPORTED on its own,
    which is why `c3`'s volume carries the factor and `c1`'s does not.

(2) THE FLOOR WAS SUMMED ACROSS EVERY REGIME. Regime 4 (pad clearing) is
    MECHANICAL and its 3,272 cells inflated the heat-driven floor by 4.3 %.
    ⚠⚠ Counting it gives 1.216 GJ/m^3, INSIDE the 1.21 +/- 0.05 band, while
    counting correctly gives 1.270, outside it: THE WRONG METHOD YIELDS THE
    BETTER-LOOKING NUMBER. That is the whole reason this is scored from
    regime labels and not from a total.

(3) THE WALL CHARGE WAS READ OFF A PER-STEP COLUMN. `lat_bank_hf` at its final
    value is 0 (per-step, max 0.1365), so the wall charge was identically zero
    and the reduction was floor-only. ⚠ `lat_bank_hf x 6.0226` is NOT the fix
    either: ARCHIVE_DONE L18670 states 6.0226 is "a sampling factor derived
    from the event log" OF P1g-2c'S RUN -- it is not a constant and does not
    transfer. ⇒ the wall path is counted from `n_voided` on regime 6, like the
    floor, which needs no factor at all.

⇒ EVERYTHING IS COUNTED FROM `n_voided`. One cell is dx*dy*dz; `n_voided` is
how many cells that event voided. (a3) is the cross-check that this agrees with
the independent `h_applied` route on the floor, where both are available.
⚠ The wall path has NO second route: regime 6 BANKS thickness, so its
`h_applied` is 0 by construction (standing rule: an assertion over a quantity
that is zero by construction proves nothing). (a3) is therefore a floor-only
cross-check, and that is stated rather than hidden.

⚠ NUMERATOR CAVEAT, CARRIED INTO THE OUTPUT. `ledger_E_flake` is cumulative
with FOUR `+=` terms (Removal.H:3102 flake+boundary+support, :3535 shift,
:3543 lateral) -- ARCHIVE_DONE L18675. It is NOT "the floor flake energy", so a
miss of a per-cent or two is a numerator question before it is physics.
⚠ P2-2f adds `foot_mech_E` as a FIFTH booked sink, but deliberately NOT into
`ledger_E_flake`: pad-clearance energy is mechanical, and folding it in would
put pad energy in the numerator while the denominator excludes pad volume --
defect (1) again, in the other variable.
"""

CELL_DX = 0.002      # amr.n_cell = 60 60 300 over [0,0.12]^2 x [-0.1,0.5]
CELL_DZ = 0.002
V_CELL = CELL_DX * CELL_DX * CELL_DZ          # 8e-9 m^3, one quadrant

R_FLOOR_FALLBACK = 0.040     # used only if no column fires in the final window

# regime labels, from Removal.H
R_FLOOR, R_PAD, R_WALL = 1, 4, 6


def _f(r, k):
    v = r.get(k)
    if v in (None, ""):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _i(r, k):
    v = _f(r, k)
    return None if v is None else int(v)


def col_radius(i, j):
    """Cell-centre radius of column (i, j). Quarter domain, axis at (0, 0)."""
    xc = (i + 0.5) * CELL_DX
    yc = (j + 0.5) * CELL_DX
    return (xc * xc + yc * yc) ** 0.5


def charged_volume(ev):
    """Heat-driven charged volume from `n_voided`, by path.

    Returns a dict. ⚠ Volumes are ONE QUADRANT and unscaled -- see the module
    docstring on why the price must not scale them.
    """
    out = {}
    for name, reg in (("floor", R_FLOOR), ("wall", R_WALL), ("pad", R_PAD)):
        rows = [r for r in ev if _i(r, "regime") == reg]
        nv = sum(_i(r, "n_voided") or 0 for r in rows)
        ha = sum(_f(r, "h_applied") or 0.0 for r in rows)
        out[name] = {
            "events": len(rows),
            "cells": nv,
            "vol_nvoided": nv * V_CELL,
            # the independent route: h_applied is a LENGTH per event, so a
            # column's charged volume is h * dx * dy.
            "vol_happlied": ha * CELL_DX * CELL_DX,
            "h_sum": ha,
        }
    # the charge: heat-driven paths only. Pad clearing is mechanical and is
    # excluded from BOTH sides of the price.
    out["charged"] = out["floor"]["vol_nvoided"] + out["wall"]["vol_nvoided"]
    return out


def a3_cross_check(ev, bar=1.0e-3):
    """(ok, rel, detail) -- floor volume from `n_voided` vs from `h_applied`.

    bar is a FRACTION (1e-3 = 0.1 %), the packet's (a3) bar.
    """
    v = charged_volume(ev)["floor"]
    a, b = v["vol_nvoided"], v["vol_happlied"]
    if a <= 0.0 and b <= 0.0:
        return False, None, ("floor path is EMPTY on both routes -- the check "
                             "would compare two zeros and prove nothing")
    rel = abs(a - b) / max(a, b)
    detail = ("floor: %d events, %d cells => %.6f L from n_voided vs %.6f L "
              "from h_applied (sum h = %.4f m); rel %.3e (bar %.1e). "
              "⚠ FLOOR ONLY: regime 6 banks thickness so its h_applied is 0 "
              "by construction and has no second route."
              % (v["events"], v["cells"], 1e3 * a, 1e3 * b, v["h_sum"],
                 rel, bar))
    return rel < bar, rel, detail


def frozen_floor(ev, t_end=None, silent_s=30.0, active_s=50.0):
    """Deliverable 4's FROZEN floor fraction.

    A floor column is FROZEN if it fired at least once and was then silent for
    the final `silent_s`. The POPULATION is set by radius: columns whose last
    firing falls in the final `active_s` define the active radius, and every
    column that ever fired inside that radius is scored.

    ⚠⚠ REV 1 SCORED `foot_carry_cols`, WHICH IS NOT THIS. MMWSpalling.H:2884 /
    :8879 define it as "columns at or above their pad height" -- LOAD-BEARING
    CONTACT, whose maximum is 162 = 100 %. A "< 10 %" bar on it measures how
    far the feet have sunk, not whether the floor stopped advancing. Rev 1's
    13.6 % is withdrawn.

    ⚠ WHAT THE RADIUS CUT EXCLUDES, AND WHY.
      - Columns OUTSIDE the active radius. These are hole wall and the rock
        beyond the bore. Many fired once early, when the floor passed them, and
        are silent thereafter BY DESIGN -- the floor has moved on and left them
        behind as wall. Scoring them would count the entire drilled sidewall as
        "frozen" and the fraction would approach 1 on a perfectly healthy leg.
      - Columns that NEVER fired. They were never part of the advancing floor,
        so "stopped advancing" is not defined for them. Including them would
        make the fraction a statement about domain width.
      Both exclusions therefore remove columns for which the bar's question is
      meaningless, not columns that would fail it.

    Returns a dict, or None if nothing fired.
    """
    last = {}
    for r in ev:
        if _i(r, "regime") != R_FLOOR:
            continue
        i, j, t = _i(r, "col_i"), _i(r, "col_j"), _f(r, "time")
        if i is None or j is None or t is None:
            continue
        k = (i, j)
        if t > last.get(k, -1.0):
            last[k] = t
    if not last:
        return None
    t_last = max(last.values())
    if t_end is None:
        t_end = t_last

    active = [k for k, t in last.items() if t >= t_end - active_s]
    if active:
        r_active = max(col_radius(*k) for k in active)
        r_from = "columns firing in the final %.0f s" % active_s
    else:
        r_active = R_FLOOR_FALLBACK
        r_from = "FALLBACK %.3f m -- no column fired in the final %.0f s" % (
            R_FLOOR_FALLBACK, active_s)

    pop = [k for k in last if col_radius(*k) <= r_active]
    frozen = [k for k in pop if last[k] < t_end - silent_s]
    return {
        "t_end": t_end,
        "r_active": r_active,
        "r_from": r_from,
        "n_pop": len(pop),
        "n_frozen": len(frozen),
        "n_fired_total": len(last),
        "n_active": len(active),
        "frac": (len(frozen) / len(pop)) if pop else None,
    }


def firing_T(ev):
    """The c2 population: reconstructed surface temperature of the firings that
    HAVE one. Returns (list, note).

    ⚠⚠ REV 1'S 851.64 K WAS COMPUTED OVER ZERO-VALUED ROWS. Measured on the
    m = 13 log: ALL 3,272 regime-4 rows and 1,858 of 1,865 regime-6 rows carry
    `T_surf` = 0 exactly -- `Removal.H`'s pad-clearance comment says so in as
    many words ("other fields 0"), and regime 6 reports its face temperature in
    `T_face_lat`, not here. Averaging those 5,130 zeros into the sample pulled
    the median from 852.61 K to 851.64 K. Both clear the 845 +/- 10 bar, so the
    verdict does not move, but 851.64 is withdrawn: it is a statistic over a
    population that includes rows where the column does not exist.

    ⚠ This is the ZERO-COUNT TRAP in a new guise -- not a count that came back
    zero, but zeros inside the sample. The population is therefore stated
    explicitly (regime 1, `T_surf` > 0) rather than taken as "every row".
    """
    ts = [_f(r, "T_surf") for r in ev if _i(r, "regime") == R_FLOOR]
    ts = [t for t in ts if t is not None and t > 0.0]
    n_all = sum(1 for r in ev if _f(r, "T_surf") is not None)
    n_zero = sum(1 for r in ev
                 if (_f(r, "T_surf") or 0.0) <= 0.0)
    note = ("population: regime %d with T_surf > 0, n = %d. ⚠ the log has %d "
            "rows carrying a T_surf field of which %d are EXACTLY 0 (all pad "
            "rows + the banking wall rows); averaging those in is what gave "
            "rev 1 its withdrawn median on the m = 13 leg (851.64 vs 852.61)."
            % (R_FLOOR, len(ts), n_all, n_zero))
    return ts, note


def clamp_fraction(ev, T_clamp=873.15):
    """(n_above, n_total, fraction) of firings whose reconstructed T_surf sits
    ABOVE the temperature where both Grimsel tables go flat.

    ⚠⚠ c2 REPORTS THIS BESIDE THE MEDIAN OR IT OVERCLAIMS. Both tables repeat
    their last value at 600 and 650 C, so the last real interval is 823-873 K
    and nothing above 873.15 K is scored by either property law. A median that
    transfers while a quarter of the population is outside the tables' validity
    is a weaker result than the median alone suggests.
    """
    ts, _ = firing_T(ev)        # ⚠ the SAME population as the median, or the
    if not ts:                  #   two numbers c2 prints describe different
        return 0, 0, None       #   samples and the pair is not interpretable.
    n = sum(1 for t in ts if t > T_clamp)
    return n, len(ts), n / len(ts)


def pctl(xs, q):
    if not xs:
        return None
    s = sorted(xs)
    k = min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))
    return s[k]
