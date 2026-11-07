#!/usr/bin/env python3
"""P2-21 z0: DOE DE-EE0005504 Fig. 34 (p. 54, test #19, virgin Barre) re-read from the PDF in this folder.

The LOWER panel carries both traces on one frame: the incident power (black, left axis 0-4 kW) and the 28 GHz
test-chamber pick-off signal (red, right axis 0-0.15 V). The page is rendered at 400 dpi (pdftoppm, poppler), the
frame found as the bounding box of the long dark axis lines, and the axes taken as the frame: x 0 at the left edge
to 2200 s at the right (the tick labels 0 / 500 / ... / 2000 sit on that scale), y 0 at the bottom edge to 4 kW /
0.15 V at the top (the "2" and "3" kW labels and the "0.10" V label check it). Per 10 s window the MEDIAN pixel row
of each colour is the reading (the traces are noisy bands a few pixels thick). Usage: fig34_digitise.py > FIG34.txt
"""
import os
import subprocess
import tempfile

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "DE-EE0005504_final_report.pdf")
DPI = 400


def page():
    d = tempfile.mkdtemp()
    subprocess.check_call(["pdftoppm", "-f", "54", "-l", "54", "-r", str(DPI), "-png", PDF, os.path.join(d, "p")])
    f = [x for x in os.listdir(d) if x.endswith(".png")][0]
    return np.asarray(Image.open(os.path.join(d, f)).convert("RGB")).astype(int)


def frame(img, y0, y1):
    """the lower panel's frame inside rows [y0, y1): columns / rows that are dark over most of the panel"""
    sub = img[y0:y1]
    dark = (sub.sum(axis=2) < 200)
    cols = np.where(dark.sum(axis=0) > 0.6 * (y1 - y0))[0]
    rows = np.where(dark.sum(axis=1) > 0.6 * (cols.max() - cols.min()))[0]
    return cols.min(), cols.max(), y0 + rows.min(), y0 + rows.max()


def main():
    img = page()
    H = img.shape[0]
    # the lower panel lies in the page's 36-60 % band (the upper panel ends above 33 %)
    xl, xr, yt, yb = frame(img, int(0.36 * H), int(0.60 * H))
    r, g, b = img[..., 0], img[..., 1], img[..., 2]
    red = (r > 170) & (g < 110) & (b < 110)
    blk = (r < 70) & (g < 70) & (b < 70)
    t_of = lambda x: 2200.0 * (x - xl) / (xr - xl)
    P_of = lambda y: 4.0 * (yb - y) / (yb - yt)
    V_of = lambda y: 0.15 * (yb - y) / (yb - yt)
    print(__doc__.split("\n\n")[0])
    print("frame (px at %d dpi): x %d-%d, y %d-%d; scale %.3f px/s, %.2f px/kW" % (
        DPI, xl, xr, yt, yb, (xr - xl) / 2200.0, (yb - yt) / 4.0))
    print("\n  t [s]  P_inc [kW]  pick-off [V]  V/kW (the transmitted FRACTION, arb.)")
    rows = []
    for t0 in range(150, 700, 10):
        x0 = int(xl + (t0 + 1) * (xr - xl) / 2200.0)
        x1 = int(xl + (t0 + 9) * (xr - xl) / 2200.0)
        ys_r = np.where(red[yt + 3:yb - 3, x0:x1])[0] + yt + 3
        ys_b = np.where(blk[yt + 3:yb - 3, x0:x1])[0] + yt + 3
        P = P_of(np.median(ys_b)) if ys_b.size > 20 else float("nan")
        V = V_of(np.median(ys_r)) if ys_r.size > 20 else float("nan")
        rows.append((t0 + 5, P, V))
        print("  %5d  %10.2f  %12.4f  %8.4f" % (t0 + 5, P, V, V / P if P > 0.3 else float("nan")))
    on = [x for x in rows if x[1] > 0.5 and x[2] > 0]
    t_on = on[0][0]
    # the collapse: the first window whose pick-off falls below 0.6 x the window before it (that window already
    # holds the vertical drop, so "just before the collapse" is the window before it)
    k_c = next(i for i in range(1, len(on)) if on[i][2] < 0.6 * on[i - 1][2])
    pre = on[:k_c]
    f0 = np.mean([x[2] / x[1] for x in pre[:2]])
    fr = [(t, (V / P) / f0) for t, P, V in pre]
    tmin, fmin = min(fr, key=lambda z: z[1])
    tmax, fmax = max(fr, key=lambda z: z[1])
    print("\nbeam on: the first window with both traces is %d s (incident %.2f kW); the collapse starts in the "
          "%d s window" % (t_on, on[0][1], on[k_c][0]))
    print("just before the collapse (%d s): incident %.2f kW; mean incident from beam-on to there %.2f kW" % (
        pre[-1][0], pre[-1][1], np.mean([x[1] for x in pre])))
    print("transmitted FRACTION (V/kW) relative to its first two beam-on windows: minimum %.2f at %d s, maximum %.2f "
          "at %d s, %.2f just before the collapse" % (fmin, tmin, fmax, tmax, fr[-1][1]))
    print("collapse: pick-off %.4f V at %d s -> %.4f V at %d s" % (pre[-1][2], pre[-1][0], on[k_c + 6][2],
                                                                   on[k_c + 6][0]))


if __name__ == "__main__":
    main()
