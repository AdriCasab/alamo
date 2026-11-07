#!/usr/bin/env python3
"""P2-22 -- the in-rock beam spread (Rayleigh range n*zR inside the dielectric) and test #19 at its TRUE spot.

DOE DE-EE0005504 Appendix Table 2: #19 standoff 0.31" = 7.9 mm (not the "1/2-2 aperture diameters" P2-21 scanned).
HE11 from the 32.5 mm aperture ~ Gaussian 1/e^2 diameter 25.4 mm (report L978) => w0 12.7 mm, zR = pi w0^2/lambda =
47.3 mm => w at the rock = 12.7 sqrt(1 + (7.9/47.3)^2) = 12.88 mm (COMPUTED; see RESULTS.txt).
Legs (the deck, the block and Fig. 34's measured power ramp exactly as P2-21's T-legs; only omega0 differs):
  A   the 17-pair Hartlieb table, omega0 0.01288, implicit dt 2 s, 400 s       (P2-21 T1a's twin at the true spot)
  F   A + beam.tan_delta_T="295.15 1.0 1268.15 1.0" (two EQUAL pairs = flat alpha 14 on the sum path; the control)
  X   A + solver=explicit timestep=0.03 amr.plot_int=333 stop_time=max(60, ceil(1.25 x A's onset))  (the fine-dt check)
  C1  the Zhang deck + one pair (flat), amr.node.any=1 amr.plot_int=5: P2-21's C1 on the new binary (its byte
      reference against P2-20's Z1 is EXPECTED to move: the deck has eps_r 5.02)
Usage: run.py LEG [LEG ...] | --argv LEG ...
"""
import importlib.machinery, importlib.util, math, os, shlex, sys
HERE = os.path.dirname(os.path.abspath(__file__))
P21 = os.path.join(HERE, "..", "p2_21_loss_table", "run.py")
ld = importlib.machinery.SourceFileLoader("p21", P21)
p21 = importlib.util.module_from_spec(importlib.util.spec_from_loader("p21", ld)); ld.exec_module(p21)
p21.LEGS_DIR = os.path.join(HERE, "legs")
extra21 = p21.extra   # the P2-21 leg table, kept for C1

W0_AP, STANDOFF, LAM = 12.7e-3, 7.9e-3, 1.070687e-2      # m; report L978, Table 2 (#19), 28 GHz
ZR = math.pi * W0_AP ** 2 / LAM
W_ROCK = W0_AP * math.sqrt(1.0 + (STANDOFF / ZR) ** 2)
OMEGA = "%.5f" % W_ROCK                                  # 0.01288

def extra(name):
    if name == "A":
        return p21.T19 + ["beam.omega0=" + OMEGA] + p21.BEAM
    if name == "F":
        return extra("A") + ["beam.tan_delta_T=295.15 1.0 1268.15 1.0"]
    if name == "X":
        t_on = p21.onset_time("A")
        if t_on is None:
            raise SystemExit("X needs A's onset, and A never melted")
        stop = max(60.0, math.ceil(1.25 * t_on))
        return extra("A") + ["solver=explicit", "timestep=0.03", "amr.plot_int=333", "stop_time=%g" % stop]
    if name == "C1":
        return extra21("C1")
    if name.startswith("SMOKE_"):
        return extra(name[6:]) + ["stop_time=20.0", "amr.n_cell=32 32 16", "amr.plot_int=100000"]
    raise KeyError(name)

p21.extra = extra
if __name__ == "__main__":
    if sys.argv[1:2] == ["--spot"]:
        print("w0 %.4f mm, zR %.2f mm, standoff %.1f mm => w(rock) %.3f mm => beam.omega0=%s" %
              (W0_AP * 1e3, ZR * 1e3, STANDOFF * 1e3, W_ROCK * 1e3, OMEGA)); sys.exit(0)
    sys.exit(p21.main(sys.argv[1:]))
