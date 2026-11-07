#!/usr/bin/env python3
"""DOE DE-EE0005504 Fig. 34 (p. 54, test #19) UPPER panel: the 137 GHz radiometer temperature (blue, left axis
0-2500 C, "uncorrected for emissivity", averaged over the radiometer's view) and the 28 GHz pick-off (red, right
axis, same trace as the lower panel). Same method as fig34_digitise.py (400 dpi render, frame = the long dark axis
lines, x 0-2200 s across the frame, y 0-2500 C bottom to top; per 10 s window the MEDIAN blue pixel row).
Usage: fig34_top_digitise.py > FIG34_TOP.txt"""
import os, subprocess, tempfile
import numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "DE-EE0005504_final_report.pdf"); DPI = 400

def page():
    d = tempfile.mkdtemp()
    subprocess.check_call(["pdftoppm", "-f", "54", "-l", "54", "-r", str(DPI), "-png", PDF, os.path.join(d, "p")])
    f = [x for x in os.listdir(d) if x.endswith(".png")][0]
    return np.asarray(Image.open(os.path.join(d, f)).convert("RGB")).astype(int)

def frame(img, y0, y1):
    sub = img[y0:y1]; dark = (sub.sum(axis=2) < 200)
    cols = np.where(dark.sum(axis=0) > 0.6 * (y1 - y0))[0]
    rows = np.where(dark.sum(axis=1) > 0.6 * (cols.max() - cols.min()))[0]
    return cols.min(), cols.max(), y0 + rows.min(), y0 + rows.max()

img = page(); H = img.shape[0]
xl, xr, yt, yb = frame(img, int(0.05 * H), int(0.34 * H))
r, g, b = img[..., 0], img[..., 1], img[..., 2]
blue = (b > 150) & (r < 110) & (g < 130)
red = (r > 170) & (g < 110) & (b < 110)
t_of = lambda x: 2200.0 * (x - xl) / (xr - xl)
T_of = lambda y: 2500.0 * (yb - y) / (yb - yt)
print(__doc__.split("\n\n")[0])
print("frame (px): x %d-%d, y %d-%d; %.3f px/s, %.3f px/C" % (xl, xr, yt, yb, (xr - xl) / 2200.0, (yb - yt) / 2500.0))
print("\n  t [s]  T_137GHz [C, uncorrected]   (red pick-off, arb. row)")
rows = []
for t0 in range(0, 2200, 10):
    x0 = int(xl + (t0 + 1) * (xr - xl) / 2200.0); x1 = int(xl + (t0 + 9) * (xr - xl) / 2200.0)
    ys = np.where(blue[yt + 3:yb - 3, x0:x1])[0] + yt + 3
    T = T_of(np.median(ys)) if ys.size > 10 else float("nan")
    rows.append((t0 + 5, T))
for t, T in rows:
    if t % 50 == 5 or (200 <= t <= 600):
        print("  %5d  %8.0f" % (t, T))
def mean(a, b): return np.nanmean([T for t, T in rows if a <= t <= b])
print("\npre-melt (beam-on 215 s to 495 s): mean %.0f C, max %.0f C" % (mean(215, 495), np.nanmax([T for t, T in rows if 215 <= t <= 495])))
jump = [(t, T) for t, T in rows if 495 <= t <= 600]
print("the jump: " + ", ".join("%d s %.0f C" % (t, T) for t, T in jump[::2]))
for a, b in ((600, 800), (800, 1200), (1200, 1500), (1500, 1800)):
    print("plateau %4d-%4d s: mean %.0f C" % (a, b, mean(a, b)))
print("after beam-off (1800 s): " + ", ".join("%d s %.0f C" % (t, T) for t, T in rows if t in (1805, 1855, 1905, 2005, 2105)))
