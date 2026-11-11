#!/usr/bin/env python3
"""P2-9 unit/face_removal -- generate ./input (provenance; run once, output committed).

The line rule needs the whole Hoeser lab configuration (robin_form = surface,
surface.follow_mask, reconstructed surface T, flake removal, microstructure), so
the input is the lab deck `validation/hoeser/input_lab` with the P2-7/what-if base
leg's command-line keys (`studies/whatif_uptake/legs/s100/argv.sh`) applied on top,
key for key -- as unit/face_list does -- and then made a FIXTURE (the packet's
`## Legs`):
    weibull.enabled = 0, spallation.sp.a0 = 8.0e-6        (one flaw, every cell)
    surface_patch.jet_closure = none, T_flame_expr = 1866.0 (a uniform gas)
    surface_patch.h_expr = 700.0                          (a uniform coefficient)
    24 x 24 x 64 cells at 1 mm, insulated sides (the deck's Neumann-0), 20 s,
    the nozzle plane 63.2 mm above the block, the patch covering it.
The jet keys (`surface_patch.jet_*` but jet_closure) and the side-face keys
(`surface_patch.side_face_*`, group B and side_face_flux, and the side-power
log interval it needs) and the four spall.lateral_* keys (which require the
side-face staircase) are DROPPED: under
jet_closure = none they abort, there is no jet, and the line rule reads no group-B or
lateral key (a flat key-0 floor has no side face either). The mesh, the key, the seeded
staircase and the step are given per leg by ./test.
Keys are unique in the output (the base's value replaced in place).
"""
import os
import shlex

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser", "input_lab")
BASE = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "whatif_uptake", "legs",
                    "s100", "argv.sh")

OVERRIDE = {
    "amr.n_cell": "24 24 64",
    "amr.max_grid_size": "16",
    "amr.blocking_factor": "2",
    "geometry.prob_lo": "0.0 0.0 0.0",
    "geometry.prob_hi": "0.024 0.024 0.064",
    "surface_patch.radius": "1.0",
    "surface_patch.nozzle_z0": "0.1272",           # block top 0.064 + 63.2 mm
    "surface_patch.jet_closure": "none",
    "surface_patch.T_flame_expr": "1866.0",          # h_expr and T_flame_expr go together
    "surface_patch.h_expr": "700.0",
    "weibull.enabled": "0",
    "spallation.sp.a0": "8.0e-6",
    "timestep": "0.008",
    "stop_time": "20.0",
    "amr.plot_int": "125",                         # every 1 s at 8 ms
    "amr.thermo.int": "1",
    "amr.thermo.plot_int": "1",
    "plot_file": "tests/MMWSpalling/unit/face_removal/output/run",
}


# ---------------------------------------------------------------------------
# P2-9 rev 2: THE CURVED FIXTURES (key 1, full domain, no symmetry plane: the
# patch centre is the domain centre). The base input above with, per leg, the
# keys of curved_keys() and the seeds of curved_seeds():
#   SPH   a closed spherical cavity R 24 mm centred in a 64 mm cube;
#   BORE  a cylinder R 12 mm, 24 mm deep, flat floor, open to the top, in a
#         48 mm cube.
# The ORIGINAL TOP SURFACE is lidded -- h = 0.001 on faces within 0.5 mm of the
# rock-top plane (s < 63.7 mm; the nozzle 63.2 mm above it), as the Hoeser
# flange is -- so that only the cavity recedes; every other face h = 700.
# Seeds are by cell centre: inside the sphere / inside the bore.
# ---------------------------------------------------------------------------
CURVED = {
    "SPH":  dict(L=0.064, R=0.024, kind="sphere"),
    "BORE": dict(L=0.048, R=0.012, depth=0.024, kind="bore"),
}


def curved_keys(name, dx):
    c = CURVED[name]
    n = int(round(c["L"] / dx))
    return {"amr.n_cell": "%d %d %d" % (n, n, n),
            "geometry.prob_hi": "%r %r %r" % (c["L"], c["L"], c["L"]),
            "amr.max_grid_size": str(16 if dx < 0.0015 else 8),
            "surface_patch.x0": repr(c["L"] / 2), "surface_patch.y0": repr(c["L"] / 2),
            "surface_patch.nozzle_z0": repr(c["L"] + 0.0632),
            "surface_patch.h_expr": "if(s<0.0637,0.001,700.0)",
            "timestep": "0.008", "amr.plot_int": "125"}


def curved_seeds(name, dx):
    """6-int boxes, one per column (i, j) that the cavity crosses."""
    import math
    c = CURVED[name]
    n = int(round(c["L"] / dx))
    xc = c["L"] / 2
    boxes = []
    for i in range(n):
        for j in range(n):
            x, y = (i + 0.5) * dx, (j + 0.5) * dx
            r2 = (x - xc) ** 2 + (y - xc) ** 2
            if c["kind"] == "sphere":
                if r2 >= c["R"] ** 2:
                    continue
                hz = math.sqrt(c["R"] ** 2 - r2)
                ks = [k for k in range(n) if abs((k + 0.5) * dx - xc) < hz]
            else:
                if r2 >= c["R"] ** 2:
                    continue
                ks = [k for k in range(n) if (k + 0.5) * dx > c["L"] - c["depth"]]
            if ks:
                boxes.append((i, j, min(ks), i, j, max(ks)))
    return boxes


def main():
    keys, order = {}, []
    for ln in open(DECK):
        s = ln.split("#", 1)[0].strip()
        if not s or "=" not in s:
            continue
        k, v = (x.strip() for x in s.split("=", 1))
        if k not in keys:
            order.append(k)
        keys[k] = v
    toks = shlex.split(open(BASE).read())
    start = next(n for n, t in enumerate(toks) if t.endswith("input_lab")) + 1
    for t in toks[start:]:
        k, v = t.split("=", 1)
        if k not in keys:
            order.append(k)
        keys[k] = v
    for k in [k for k in order if (k.startswith("surface_patch.jet_")
              and k != "surface_patch.jet_closure")
              or k.startswith("surface_patch.side_face_")
              or k == "surface_patch.side_power_interval"
              or k.startswith("spall.lateral_")]:
        order.remove(k)
        del keys[k]
    for k, v in OVERRIDE.items():
        if k not in keys:
            order.append(k)
        keys[k] = v
    with open(os.path.join(HERE, "input"), "w") as f:
        f.write("# P2-9 unit/face_removal -- GENERATED by ./make_input.py (do not edit by hand):\n")
        f.write("# validation/hoeser/input_lab + studies/whatif_uptake/legs/s100/argv.sh,\n")
        f.write("# made a fixture: uniform gas (T_flame 1866 K, h 700), one flaw (8 um),\n")
        f.write("# a 24 x 24 x 64 mm block at 1 mm. Mesh, key and seeds per leg in ./test.\n")
        for k in order:
            f.write(f"{k} = {keys[k]}\n")
    print(f"wrote {len(order)} keys")


if __name__ == "__main__":
    main()
