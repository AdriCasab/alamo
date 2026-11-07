#!/usr/bin/env python3
"""THE RELIEF DIAGNOSTIC (2026-11-05, advisor-run at the user's request). Is R = 0.531 geometry or material?
Two legs on ONE temperature history: the Hoeser lab flame (P2-16 HOS_R argv) on a FULL 100 x 100 mm block, spot
centred, free sides, bottom roller, 2 mm cells, removal ON through the column path (the production surface form needs the mask, the mask needs removal; face_list=0 because the face list refuses the operator), to 1.0 s (HOS_R fires 0.84 s).
  L  spall.stress_source = local_thermoelastic, el.type = disable, R = 0.531294 (the production form)
  Lc spall.stress_source = local_thermoelastic on the CELL temperature (spall.surface_T = cell_average)
  M  spall.stress_source = mechanics (the elastic operator; that site applies no relief), el.type = static,
     spall.surface_T = cell_average (the reconstruction is refused with the operator)
Read with reduce.py: sigma_op / sigma_1D = R * Sp_M / Sp_L column by column (same flaw, same K_Ic, same T)."""
import hashlib, os, shlex, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser", "input_lab")
HX = ("if((r>0.02)and(s<0.0633),0.001,1189.763953*533*(533+0.3575*44*pow(r/0.0079,1.285))"
      "/pow(533+44*pow(r/0.0079,1.285),2)*pow(0.0632/max(0.0632,s),1.3185305455))")
COMMON = ["beam.P0=0.0", "surface_patch.mode=convective_flame", "surface_patch.pinned_idle_cycles=1.0e9",
          "surface_patch.side_power_interval=10.0", "surface_patch.enabled=1", "surface_patch.x0=0.05",
          "surface_patch.y0=0.05", "surface_patch.radius=0.08", "surface_patch.nozzle_descent=prescribed",
          "surface_patch.nozzle_feed=0.0", "surface_patch.nozzle_collision_radius=0.004",
          "surface_patch.jet_closure=enthalpy", "surface_patch.jet_stagnation=decay",
          "surface_patch.jet_entrained_mass=1", "surface_patch.jet_T_ent_mode=fixed", "surface_patch.jet_T_ent=293.15",
          "surface_patch.jet_mdot=3.4e-3", "surface_patch.jet_cp=1250.0", "surface_patch.jet_D=0.0079",
          "surface_patch.jet_T_nozzle=1866.0", "surface_patch.jet_decay_diameter=momentum",
          "surface_patch.jet_De_ref=0.003745", "surface_patch.jet_bin_update=exponential",
          "surface_patch.jet_negative_flux=count", "surface_patch.jet_profile_interval=1.0", "confining.p=0.0",
          "amr.thermo.plot_int=6", "amr.thermo.int=1", "amr.max_grid_size=24", "timestep=0.008", "weibull.V0=1.0e-9",
          "surface_patch.jet_core_length=30.0", "surface_patch.nozzle_z0=0.2132", "surface_patch.h_expr=" + HX,
          # 48^3 with 24-cell boxes: coarsenable by 2 four times (48/24/12/6/3), which multigrid needs; the first
          # attempt at 50^3 / 25-cell boxes could coarsen only once and MLMG failed on the first heated solve
          "amr.n_cell=48 48 48", "geometry.prob_lo=0.0 0.0 0.05", "geometry.prob_hi=0.1 0.1 0.15",
          "surface_patch.face_list=0", "spall.buckle_span=0", "stop_time=1.0", "amr.plot_int=25"]
LOCAL = ["spall.stress_source=local_thermoelastic", "el.type=disable"]
# el.interval = 10: one elastic solve per 80 ms (the stress follows the temperature quasi-statically; at
# el.interval = 1 a single solve took > 1 min on 125k cells). verbose = 1 logs the MLMG iterations.
MECH = ["spall.stress_source=mechanics", "el.type=static", "el.interval=10", "el.time_evolving=1",
        "el.bc.type=zlo_roller_321", "el.print_model=0", "el.solver.verbose=2", "el.solver.nriters=1",
        "el.solver.max_iter=200", "el.solver.bottom_max_iter=200", "el.solver.tol_rel=1.0e-5",
        "el.solver.tol_abs=1.0e-16", "el.solver.bottom_solver=cg", "el.solver.normalize_ddw=1",
        "el.zero_out_displacement=1", "el.max_coarsening_level=4"]

def run(name, extra):
    d = os.path.join(HERE, name); os.makedirs(d, exist_ok=True)
    common = list(COMMON)
    if "surface_patch.robin_form=cell" in extra:   # read under the surface form only; strict scope otherwise
        common = [t for t in common if not t.startswith("surface_patch.pinned_idle_cycles=")]
    argv = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, DECK] + common + extra + ["plot_file=" + os.path.join(d, "run")]
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(a) for a in argv) + "\n")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    line = "%s rc=%d wall %.0f s" % (name, rc, time.time() - t0)
    print(line, flush=True); open(os.path.join(HERE, "RUN.log"), "a").write(line + "\n")
    return rc

if __name__ == "__main__":
    open(os.path.join(HERE, "BINARY.txt"), "w").write(hashlib.sha256(open(BIN, "rb").read()).hexdigest() + "\n")
    legs = sys.argv[1:] or ["L", "M"]
    # Lc: the local law on the CELL temperature (the operator's stress is the cell's, so the ratio is
    # read Lc vs M, like for like); L keeps the production reconstruction for the onset reading.
    # The interlock chain (MMWSpalling::Parse): robin_form = surface needs surface_T = reconstructed, the flake
    # removal energy needs it too, and the operator REFUSES it. So both comparison legs run the cell Robin form on
    # the cell temperature with the cell energy form: the same thermal BC in both, and the stress ratio at the same
    # cell temperature is the geometric relief. L alone keeps the production form, for the onset reading.
    CELL = ["spall.surface_T=cell_average", "spall.removal_energy=cell", "surface_patch.robin_form=cell"]
    FINE = ["amr.n_cell=96 96 96", "amr.max_grid_size=24"]      # 1.04 mm; still coarsenable 96/48/24/12/6/3
    for n in legs:
        base = {"L": LOCAL, "Lc": LOCAL + CELL, "M": MECH + CELL}[n.rstrip("1")]
        run(n, base + (FINE if n.endswith("1") else []))
