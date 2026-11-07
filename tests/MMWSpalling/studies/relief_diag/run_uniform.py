#!/usr/bin/env python3
"""THE RESOLUTION LADDER for the operator's surface stress. Under UNIFORM heating of the whole top face of a
block whose sides are free, the exact answer at the centre (many skin depths from any edge) is FULL constraint:
sigma_op / sigma_1D = 1. Whatever the operator reads below 1 at the centre is its discretisation error for a skin
of ell/dx cells. Fixture: unit/face_removal/input (uniform gas 1866 K, h 700, 24 x 24 x 64 mm block, one 8 um flaw),
legs at 2 / 1 / 0.5 mm, each as Lc (local law on the cell T) and M (operator), same keys as run.py's CELL legs.
Read with reduce_uniform.py: the centre-column ratio vs time and vs ell/dx."""
import hashlib, os, shlex, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(HERE, "U", "input_cell")   # the fixture minus its pinned_idle_cycles line (unread under the cell form, strict scope)
MESH = {"2": ("12 12 32", "0.008", "8", "2"), "1": ("24 24 64", "0.008", "16", "3"), "05": ("48 48 128", "0.004", "32", "4")}
COMMON = ["surface_patch.face_list=0", "spall.buckle_span=0", "spall.surface_T=cell_average",
          "spall.removal_energy=cell", "surface_patch.robin_form=cell", "stop_time=1.0", "amr.thermo.plot_int=1"]
LOCAL = ["spall.stress_source=local_thermoelastic", "el.type=disable"]
MECH = ["spall.stress_source=mechanics", "el.type=static", "el.interval=10", "el.time_evolving=1",
        "el.bc.type=zlo_roller_321", "el.print_model=0", "el.solver.verbose=2", "el.solver.nriters=1",
        "el.solver.max_iter=200", "el.solver.bottom_max_iter=200", "el.solver.tol_rel=1.0e-5",
        "el.solver.tol_abs=1.0e-16", "el.solver.bottom_solver=cg", "el.solver.normalize_ddw=1",
        "el.zero_out_displacement=1"]

def run(name, extra):
    d = os.path.join(HERE, "U", name); os.makedirs(d, exist_ok=True)
    argv = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, DECK] + COMMON + extra + ["plot_file=" + os.path.join(d, "run")]
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(a) for a in argv) + "\n")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    line = "U/%s rc=%d wall %.0f s" % (name, rc, time.time() - t0)
    print(line, flush=True); open(os.path.join(HERE, "RUN.log"), "a").write(line + "\n")

if __name__ == "__main__":
    for dxk in (sys.argv[1:] or ["2", "1", "05"]):
        n, dt, mgs, mcl = MESH[dxk]
        mesh = ["amr.n_cell=" + n, "timestep=" + dt, "amr.max_grid_size=" + mgs, "amr.plot_int=" + str(int(round(0.2 / float(dt))))]
        run("Lc" + dxk, mesh + LOCAL)
        run("M" + dxk, mesh + MECH + ["el.max_coarsening_level=" + mcl])
