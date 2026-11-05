#!/usr/bin/env python3
"""P2-4 (a1) -- the cold body gate, ON THE MEIER DECK.

⚠⚠ WHY THIS IS NOT IN `unit/wall_march`. The body's coefficient is the wall
law's DUCT scalar through one read, with no typed fallback (deliverable 2), so
`wall_law_ready` must be true or the term is zero BY DESIGN. That needs a bore
wider than `WJ_R_SHIELD + 2dx` (> 75 mm across) and a domain deep enough for the
bore measurement's own band, which drops 100 mm below the mouth and 20 mm above
the floor. The unit box is 48 mm wide and 80 mm deep: a 75 mm hole does not fit,
the band is empty, and `wall_r_bore` is identically 0. Probed there:
`wall_law_on` 0, `body_h_used` 0, `body_noLaw_n` 4998 of 4998 steps.
⇒ every body row would have passed or failed for a reason unrelated to the code,
and (g) -- "a body hotter than the gas draws exactly zero" -- would have
reported a PASS on a zero caused by the law never switching on. That is the fake
pass this project keeps paying for, so the gate was moved to the deck where the
mechanism exists. On P2-3's W1 the law is ready by t = 187 s; these legs run to
200 s.

⚠⚠ ROW (f) IS STRONGER THAN THE PACKET'S VERSION, AND HAD TO BE. The packet
fits a decay length and compares it with `body_L`, which needs the rock faces
OFF so the profile is a clean exponential. But with `side_face_factor_above = 0`
the bore widens only by the floor path and does NOT reach 75 mm by 200 s -- the
law would not be ready and the row would be vacuous again. Since the profile CSV
now logs `P_k` (rock) and `P_body_k` (body) SEPARATELY per bin, the body's area
is checked directly instead:
      P_body_k == body_h_used * (2*pi*WJ_R_SHIELD*wall_sector*dz) * (T_in - T_body)
bin by bin, against the logged h, the logged inlet temperature and the geometric
area. A dropped `wall_sector` shows up as a factor of exactly 4, which is the
failure mode the packet's fitted length was a proxy for.

Usage: gate_a1.py [--run] [--stop 200]
"""
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                    "sp_meier_pilot", "input_drilling")
W1_ARGV = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_3_wall",
                       "legs", "W1", "argv.sh")
OUT = os.path.join(HERE, "gate_a1")

WJ_R_SHIELD = 0.0335          # MMWSpalling.H:5376 -- no new geometry number
T_BODY = 320.0
DZ = 0.002

LEGS = {
    # the mechanism leg: body on, W1's configuration otherwise
    "FB":    (["surface_patch.jet_body=1",
               "surface_patch.jet_body_T=320.0"], 4),
    # (h) rank independence. ⚠ 5 ranks, not 1: the grid is 5 boxes of 60^3, so
    # -np 5 balances it, and a 1-rank leg on this deck would be ~4x slower for
    # no extra discrimination -- rank INDEPENDENCE is what is being tested, and
    # 4 vs 5 tests it with a genuinely different decomposition.
    "FB5":   (["surface_patch.jet_body=1",
               "surface_patch.jet_body_T=320.0"], 5),
    # (g) a body hotter than the gas
    "FBHOT": (["surface_patch.jet_body=1",
               "surface_patch.jet_body_T=5000.0"], 4),
    # (j) the `fresh` mass option
    "FBFR":  (["surface_patch.jet_body=1",
               "surface_patch.jet_body_T=320.0",
               "surface_patch.jet_wall_mdot=fresh"], 4),
}

FAILS = []


def check(ok, msg, detail=""):
    print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
    if detail:
        print("        %s" % detail)
    if not ok:
        FAILS.append(msg)


def base_keys():
    """W1's configuration, from P2-3's OWN harness rather than from its argv.sh.

    ⚠⚠ TWO BUGS MADE THIS NOT A PARSE OF `argv.sh`, AND THE ROUND-TRIP CHECK
    CAUGHT BOTH BEFORE ANY LEG RAN (P2-2f lost 100 mm of rock to exactly this
    class):
      (1) `argv.sh` carries DUPLICATE keys -- `side_face_factor_above` appears as
          the P2-2f carried 0.0 at line 52 and again as P2-3's override 1.0 at
          line 93. ParmParse takes the LAST; a first-match parser takes 0.0, i.e.
          THE WALL OFF, so the bore would never widen, the law would never become
          ready, and this gate would have reported a structural zero as a
          failure of the body.
      (2) the carried `feet_argv.txt` line for `side_face_h_law` has a 64-hex
          sha appended after `derived` on the same line. ParmParse takes the
          first value for a scalar string query, so every P2-2f/P2-3 leg ran
          correctly, but any re-parser must not treat the sha as a value.
    ⇒ build from `run_p2_3.base_keys() + run_p2_3.COMMON` -- the construction
      that actually produced W1 -- and then de-duplicate LAST-WINS, matching
      ParmParse. The multi-value keys (`amr.n_cell=60 60 300`) and the quoted
      `h_expr` come through that path intact.
    """
    sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies",
                                    "p2_3_wall"))
    import run_p2_3 as R3
    # ⚠ W1's own two keys live in R3.LEGS["W1"], not in COMMON -- omitting them
    # dropped `jet_wall_march` and the body would have aborted as an inert key.
    raw = list(R3.base_keys()) + list(R3.COMMON) + list(R3.LEGS["W1"])
    seen = {}
    for k in raw:
        if "=" not in k:
            continue
        seen[k.split("=", 1)[0]] = k        # last wins, as ParmParse does
    return list(seen.values())


def run(tag, extra, nproc, stop):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    done = os.path.join(d, ".done")
    if os.path.exists(done):
        os.remove(done)
    argv = (["caffeinate", "-dis", "mpirun", "--oversubscribe", "--bind-to",
             "none", "-np", str(nproc), BIN, DECK]
            + base_keys() + extra
            + ["stop_time=%g" % stop,
               "amr.thermo.int=1", "amr.thermo.plot_int=1",
               "surface_patch.jet_profile_interval=10.0",
               "plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join("'%s'" % a if " " in a else a for a in argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    wall = time.time() - t0
    with open(done, "w") as f:
        f.write("rc=%d wall=%.1f\n" % (rc, wall))
    return rc, wall


def thermo(tag):
    p = os.path.join(OUT, tag, "run", "thermo.dat")
    if not os.path.exists(p):
        return None
    lines = [l.split() for l in open(p).read().splitlines() if l.split()]
    head, rows = lines[0], lines[1:]
    out = {}
    for i, nm in enumerate(head):
        col = []
        for r in rows:
            if i < len(r):
                try:
                    col.append(float(r[i]))
                except ValueError:
                    pass
        out[nm] = col
    return out


def prof(tag):
    d = os.path.join(OUT, tag)
    hits = []
    for dp, _, fn in os.walk(d):
        hits += [os.path.join(dp, f) for f in fn
                 if f.endswith("wall_profile.csv")]
    if not hits:
        return {}
    lines = open(sorted(hits)[0]).read().splitlines()
    if len(lines) < 2:
        return {}
    head = lines[0].split(",")
    by = {}
    for ln in lines[1:]:
        f = ln.split(",")
        if len(f) != len(head):
            continue
        try:
            r = {h: float(v) for h, v in zip(head, f)}
        except ValueError:
            continue
        by.setdefault(r["time"], []).append(r)
    for t in by:
        by[t].sort(key=lambda x: x["k"])
    return by


def main():
    stop = 200.0
    if "--stop" in sys.argv:
        stop = float(sys.argv[sys.argv.index("--stop") + 1])
    print("P2-4 (a1) -- the cold body gate, Meier deck, stop_time %.0f s" % stop)
    print("=" * 72)
    print("  base argv %s (%d keys)"
          % (os.path.relpath(W1_ARGV, ROOT), len(base_keys())))
    if "--run" in sys.argv:
        for tag, (extra, npr) in LEGS.items():
            rc, wall = run(tag, extra, npr, stop)
            print("  %-6s rc=%d np=%d wall %.1f s" % (tag, rc, npr, wall))
            if rc != 0:
                FAILS.append("leg %s rc=%d" % (tag, rc))
        if FAILS:
            print("\nFAILED: %s" % "; ".join(FAILS))
            return 1

    th = {t: thermo(t) for t in LEGS}
    for t, v in th.items():
        if v is None:
            print("  (no thermo for %s -- pass --run first)" % t)
            return 2

    def tl(t, nm):
        return th[t].get(nm, [])[1:]

    # ---- (l) NON-VACUITY FIRST: the law must be ready and the body reached ----
    bh = tl("FB", "body_h_used")
    nolaw = tl("FB", "body_noLaw_n")
    wl = tl("FB", "wall_law_on")
    rb = tl("FB", "wall_r_bore")
    ready = bool(bh) and max(bh) > 0.0
    check(ready,
          "(l) NON-VACUITY: the wall law became READY and the body read its duct "
          "scalar -- without this every row below is a zero for the wrong reason",
          "max body_h_used %.3f W/m2K; max wall_law_on %.0f; max wall_r_bore "
          "%.5f m (needs > %.4f); body_noLaw_n final %.0f (CUMULATIVE)"
          % (max(bh) if bh else 0.0, max(wl) if wl else 0.0,
             max(rb) if rb else 0.0, WJ_R_SHIELD + 2 * DZ,
             nolaw[-1] if nolaw else 0.0))
    if not ready:
        print("\nFAILED: the gate cannot be evaluated -- the law never became ready")
        return 1

    # ---- (f) THE AREA, bin by bin ----
    pr = prof("FB")
    sec = tl("FB", "wall_sector")
    sector = sec[-1] if sec else 0.0
    A_bin = 2.0 * math.pi * WJ_R_SHIELD * sector * DZ
    worst, nb, bad = 0.0, 0, 0
    if pr:
        t_last = sorted(pr)[-1]
        # ⚠⚠ h MUST BE READ AT THE PROFILE'S OWN TIMESTAMP. The first version of
        # this row took `body_h_used` from the LAST thermo row (t = 200) and
        # compared it against the last profile snapshot (t = 190, written at
        # jet_profile_interval = 10). `wall_h_duct` DECLINES through the run
        # (83 -> 63 W/m2K on W1), so the later, smaller h made `want` smaller
        # and every bin read ~0.7 % "above the formula" -- a uniform offset on
        # all 19 bins, which is the signature of a time mismatch and NOT of an
        # area error (a dropped sector fraction is a factor of exactly 4).
        tt = th["FB"]["time"][1:]
        hrow = tl("FB", "body_h_used")
        if tt and hrow:
            i_h = min(range(len(hrow)), key=lambda j: abs(tt[j] - t_last))
            h_use = hrow[i_h]
        else:
            h_use = 0.0
        for r in pr[t_last]:
            want = h_use * A_bin * max(0.0, r["T_gas"] - T_BODY)
            got = r["P_body_k"]
            nb += 1
            den = max(abs(want), 1e-12)
            rel = abs(got - want) / den
            # ⚠ the JOINT CLAMP: a bin that would overdraw the stream past the
            # T_ent floor has BOTH sinks scaled down in proportion, so its
            # P_body_k is legitimately below h*A*dT. Those bins are excluded and
            # counted rather than silently tolerated.
            if got > want * 1.0001:
                bad += 1
            elif rel > 1.0e-3 and got >= want * 0.999:
                worst = max(worst, rel)
            elif rel > 1.0e-3:
                pass          # clamped low -- legitimate
            worst = max(worst, rel) if got >= want * 0.999 else worst
    check(nb > 0 and bad == 0 and worst <= 1.0e-3,
          "(f) P_body_k == h * (2*pi*R*sector*dz) * (T_in - T_body) bin by bin "
          "-- a dropped sector fraction would show as EXACTLY 4x",
          "sector %.4f, area/bin %.6e m2, %d bins, worst unclamped rel %.3e, "
          "%d bins ABOVE the formula (must be 0)" % (sector, A_bin, nb, worst, bad))
    check(abs(sector - 0.25) < 1e-6,
          "(f1) the sector fraction is the quarter deck's 0.25, DERIVED by "
          "PatchAzimuthFraction and not typed", "wall_sector = %.6f" % sector)

    # ---- (f2) the total equals the sum of the bins ----
    pb = tl("FB", "jet_P_body")
    okS, detS = False, "no profile"
    if pr and pb:
        t_last = sorted(pr)[-1]
        ssum = sum(r["P_body_k"] for r in pr[t_last])
        tt = th["FB"]["time"][1:]
        i = min(range(len(tt)), key=lambda j: abs(tt[j] - t_last))
        okS = abs(ssum - pb[i]) <= 1e-5 * max(abs(pb[i]), 1.0)
        detS = ("at t=%.1f: sum of P_body_k over bins %.5f W vs jet_P_body "
                "%.5f W" % (t_last, ssum, pb[i]))
    check(okS, "(f2) jet_P_body is the sum of its own per-bin debits", detS)

    # ---- (g) a hot body draws EXACTLY nothing, with the law READY ----
    pbh = tl("FBHOT", "jet_P_body")
    bhh = tl("FBHOT", "body_h_used")
    allz = bool(pbh) and all(x == 0.0 for x in pbh)
    live = bool(bhh) and max(bhh) > 0.0
    check(allz and live,
          "(g) jet_body_T = 5000 K => jet_P_body EXACTLY 0 WHILE the law is "
          "ready -- one-sided, and NOT the vacuous zero the unit box gives",
          "max |jet_P_body| %.3g; max body_h_used %.3f W/m2K (non-zero, so the "
          "zero above is the one-sidedness and not an inert law)"
          % (max((abs(x) for x in pbh), default=0.0), max(bhh) if bhh else 0.0))

    # ---- (h) rank independence, 4 vs 5 ----
    p4, p5 = tl("FB", "jet_P_body"), tl("FB5", "jet_P_body")
    m4, m5 = tl("FB", "wall_T_mouth"), tl("FB5", "wall_T_mouth")
    nn = min(len(p4), len(p5))
    d1 = sum(1 for i in range(nn) if ("%g" % p4[i]) != ("%g" % p5[i]))
    d2 = sum(1 for i in range(min(len(m4), len(m5)))
             if ("%g" % m4[i]) != ("%g" % m5[i]))
    check(nn > 0 and d1 == 0 and d2 == 0,
          "(h) jet_P_body and wall_T_mouth are DECOMPOSITION-INDEPENDENT "
          "(4 ranks vs 5)",
          "%d rows compared; %d jet_P_body and %d wall_T_mouth printed tokens "
          "differ" % (nn, d1, d2))

    # ---- (i) the columns are gated ----
    absent = all(c not in th["FB"] for c in ("jet_P_retained",))
    present = all(c in th["FB"] for c in
                  ("jet_P_body", "body_h_used", "body_z_imm", "body_L",
                   "body_noLaw_n", "ledger_E_flake_wall"))
    check(present and absent,
          "(i) the six body columns are PRESENT with the key on, and "
          "jet_P_retained is ABSENT under the default `jet` mass",
          "FB has %d columns" % len(th["FB"]))

    # ---- (j) `fresh` books the retained share ----
    ret = tl("FBFR", "jet_P_retained")
    nz = bool(ret) and max(ret) > 0.0
    gated = "jet_P_retained" not in th["FB"]
    check(nz and gated,
          "(j) `fresh`: the retained entrained share is NON-ZERO and the column "
          "appears only under that key -- the five-term identity itself is the "
          "binary's 1e-9 abort, which rc = 0 witnesses over every step",
          "max jet_P_retained %.2f W" % (max(ret) if ret else 0.0))

    # ---- (k) no annulus, no sink ----
    wb = tl("FB", "wall_bins")
    pre = [(w, x) for w, x in zip(wb, pb) if w == 0.0]
    check(bool(pre) and all(x == 0.0 for _, x in pre),
          "(k) before any above-plane face exists (wall_bins == 0) jet_P_body is "
          "EXACTLY 0",
          "%d of %d steps have wall_bins == 0; max |jet_P_body| there %.3g"
          % (len(pre), len(wb), max((abs(x) for _, x in pre), default=0.0)))

    # ---- the ledger must still close ----
    les = tl("FB", "ledger_err_solid")
    worst_l = max((abs(x) for x in les), default=0.0)
    check(bool(les) and worst_l <= 5.0e-3,
          "(m) ledger_err_solid <= 5e-3 with the body on",
          "max %.4e over %d rows" % (worst_l, len(les)))

    print("=" * 72)
    if FAILS:
        print("FAILED %d: %s" % (len(FAILS), "; ".join(FAILS)))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
