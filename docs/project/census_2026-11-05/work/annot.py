# Hand-written classification, guard and production status per key.
# (class, guard for readers to execute, production status, notes/evidence)
G_MS = "microstructure.nphases present"
G_COL = "removal on & surface_patch.face_list = 0 (Removal.H after the :543 face-list dispatch)"
G_FL = "surface_patch.face_list = 1 (FaceRemoval.H / face-list builders)"
G_SPW = "spallation.model = sp_weibull"
G_JET = "surface_patch.enabled & h_expr & jet_closure = enthalpy (JetEnthalpyMarch, called BuildFlameColumns:2234)"
G_FEET = "surface_patch.nozzle_descent = feet (FeetDescent, called BuildFlameColumns:2123)"
G_MECH = "el.type != disable (UpdateModel runs from Base::Mechanics::TimeStepBegin, Base/Mechanics.H:167,176)"
G_CONST = "no material.type & no microstructure (AdvanceConstAlpha / const-alpha branch of ConvertTtoH, InvertHtoT)"
G_MAT = "material.type given (AdvanceMaterial / AdvanceMaterialImplicit; Advance:1458-1459)"
G_BEAM = "beam power > 0 at time t (beam.P0 or schedule); readers copy it every step, BeamSource is 0 at P=0"
A = {}
def a(k, cls, guard, prod, note=""):
    A[k] = (cls, guard, prod, note)

a("hc.heat.alpha", "PROD", "AdvanceConstAlpha: " + G_CONST + "; Removal.H:297: spall.enabled & !spall.alpha (runs on EVERY removal step, before the :543 dispatch)",
  "set 6.9022e-7 both legs; READ every step at Removal.H:297 only as the h_spall positivity abort (:299-301); value otherwise inert on the face list",
  "REQUIRED twice: HeatConduction.H:79 and MMWSpalling.H:153. The rock file calls it inert; it is not removable while Removal.H:299 aborts on alpha_eff <= 0")
for k in ["phase.T_ref", "phase.T_m", "phase.T_v", "phase.lh_m", "phase.lh_v"]:
    a(k, "ZHANG", G_CONST + ("; phase.T_ref also UpdateModel:1362 homogeneous eigenstrain (" + G_MECH + ", no microstructure)" if k == "phase.T_ref" else ""),
      "not set; not read (microstructure path)", "const-alpha (Step 2-4) thermal path; only unit/stefan + Hu thermoelastic decks set any")
for k in ["losses.epsilon", "losses.sigma_SB", "losses.h_conv", "losses.T_amb", "losses.epsilon_high", "losses.T_eps_break"]:
    a(k, "PROD", "none (every conduction kernel, the jet march, the reconstruction, the face list, the ledger)",
      "set on both legs (rock file / model-form / decks); read every step", "")
A["losses.T_eps_break"] = ("PROD", "none", "set 1.0e30 (= default) on both legs via rock file; read every step",
                           "the set value is the default: the two-emissivity switch is unused on production; key could be retired with losses.epsilon_high")
a("material.type", "ZHANG", G_MAT, "not set", "select_default<Zhang,Constant>: also parses material.zhang.* / material.constant.* (see grouped rows)")
a("solver", "ZHANG", G_MAT + " (Advance:1458 tests solver_kind only with material)", "set 'explicit' (= default) on both legs; reader never reached (no material)", "")
for k in ["implicit.tol_rel", "implicit.tol_abs", "implicit.max_iter", "implicit.verbose"]:
    a(k, "ZHANG", G_MAT + " & solver = implicit", "not set", "only zhang/zhang_oglesby sets any")
a("alpha", "MECH", G_MECH + " & no microstructure (UpdateModel homogeneous branch :1361-1398)",
  "set '0.0 0.0' both legs; never read (el.type = disable, microstructure on)", "top-level eigenstrain array, ThermoElastic convention")
a("microstructure.nphases", "PROD", "none (presence selects the microstructure path everywhere)", "set 1 both legs", "presence-keyed: contains() at :292 gates ~45 other keys (report ii)")
a("microstructure.mode", "PROD", G_MS, "set expression both legs; read at InitializeMicrostructure", "")
a("microstructure.number_of_grains", "MICRO", G_MS + " & mode = voronoi", "not set", "Voronoi")
a("microstructure.seed", "MICRO", G_MS + " & mode = voronoi", "not set", "Voronoi")
a("microstructure.phase_expr", "PROD", G_MS + " & mode = expression", "set '0' both legs", "")
a("microstructure.grain_expr", "MICRO", G_MS + " & mode = expression", "not set", "grain topology (Step 6c)")
a("microstructure.phase<i>.name", "DEAD", "-", "set 'grimsel' both legs; never read", "parsed :338 into MineralPhase::name; no reader anywhere (grep '\\.name\\b' finds only the parse). Strict? NO: top-level, so deleting the parse makes decks that set it silently inert, not aborting")
a("microstructure.phase<i>.fraction", "PROD", G_MS, "set 1.0 both legs; read InitializeMicrostructure:7099 (phase sampling CDF)", "REQUIRED; single phase makes the value irrelevant (normalised at :418-421)")
a("microstructure.phase<i>.kappa", "PROD", G_MS, "set 1.5 both legs; read at init (field write, overwritten each step by UpdatePropertyFields under law=grimsel) and BuildSurfaceTColumns:3009 fallback only if the field is unavailable",
  "REQUIRED; value inert under law=grimsel (rock file comment). Removal.H:702 read is column path only")
a("microstructure.phase<i>.beta", "PROD", G_MS, "set 8.0e-6; read ElasticConstants (Removal.H:36) -> pB, unused by ScanSigmaT/BuckleRatio when RockLaw() is non-null; init field write",
  "REQUIRED; value inert under spall.stress_law = grimsel")
a("microstructure.phase<i>.E", "PROD", G_MS, "set 30.0e9; read ElasticConstants Removal.H:35,37,43 -- the E>0 & mu>0 test gates ec.ok, and UpdateSpAfterMechanics:6023 ABORTS without it",
  "REQUIRED; value inert under law (nu from grimsel.nu since P2-14a) but positivity is load-bearing")
a("microstructure.phase<i>.mu", "PROD", G_MS, "set 11.538e9; same as E (Removal.H:37,43; abort :6023)", "REQUIRED; value inert under law")
a("microstructure.phase<i>.rho", "PROD", G_MS, "set 2750.0; read InitializeMicrostructure:7033 (phase_materials, unused under law) and :7233 (non-law branch only)", "REQUIRED; value inert under law=grimsel")
a("microstructure.phase<i>.Cp", "PROD", G_MS, "set 790.0; as rho (:7034, :7234 non-law only)", "REQUIRED; value inert under law=grimsel")
a("microstructure.phase<i>.alpha_attenuation", "DEAD", "-", "set 2000.0 both legs; never read",
  "parsed :347-348 into MineralPhase::alpha_atten; zero readers in the 3 files (the beam reads beam.alpha_expr). REQUIRED by the parser, so deleting it means deleting the query_required and the deck lines; MODEL_CARD §12 item 11 already says 'never read'")
a("microstructure.phase<i>.A_D", "DAMAGE", G_MS + " & damage.enabled (UpdateDamageAfterMechanics:5244)", "set 0.0 both legs (REQUIRED by parser); never read on production", "REQUIRED on every microstructure deck although only the damage law reads it")
a("microstructure.phase<i>.law", "PROD", G_MS, "set grimsel both legs; read UpdatePropertyFields:4516,4522, InitializeMicrostructure:7057,7226", "")
a("microstructure.phase<i>.T_ref", "PROD", G_MS, "set 293.15 (= IC); read BuildSurfaceTColumns:2935, UpdateSpAfterMechanics:6029, FrAssignTs FR:228, RecedeAlongLines FR:468, init", "the stress-free reference T0 of the reconstruction")
for k in ["microstructure.phase<i>.T_m", "microstructure.phase<i>.fus_width", "microstructure.phase<i>.L_m",
          "microstructure.phase<i>.T_vap_lo", "microstructure.phase<i>.T_vap_hi", "microstructure.phase<i>.L_v"]:
    a(k, "BEAM", G_MS + " (table build, InitializeMicrostructure:7037-7066; T_vap_lo/L_m/L_v also Removal.H:1667-1676 vapour regime, " + G_COL + " & vapor.enabled)",
      "not set; READ on production at default (latent heat copied into the Grimsel table, :7061-7066) -- default = no phase change", "melt/vapour latent heat; unit/latent_heat is the only witness under the law (comment :7052-7055)")
for k in ["microstructure.phase<i>.alpha_beta.enabled", "microstructure.phase<i>.alpha_beta.T0", "microstructure.phase<i>.alpha_beta.width",
          "microstructure.phase<i>.alpha_beta.dV", "microstructure.phase<i>.alpha_beta.linear_strain"]:
    a(k, "MECH", G_MS + "; physics effect only through UpdateModel (" + G_MECH + "); InitializeMicrostructure also writes the plotted thermal_free_strain field",
      "not set; read at init at default (enabled=0 => zero strain)", "quartz alpha-beta eigenstrain; only unit/alpha_beta_transition sets it" + ("; linear_strain is set by NO input" if k.endswith("linear_strain") else ""))
GR = {"cp.a": "Grimsel::cp:125, H_of:145", "cp.b": "Grimsel::cp:125, H_of:146", "cp.th0": "Grimsel::cp:125, H_of:146",
      "lambda.a": "Grimsel::lambda:138 (kappa:129)", "lambda.th0": "Grimsel::lambda:138",
      "rho": "Grimsel::rho:120", "E0": "Grimsel::E_of:170 (mu_of, sigma_thermal)", "nu": "Grimsel::mu_of:175, sigma_thermal:242; ElasticConstants Removal.H:43 (law->nu)",
      "K_Ic0": "Grimsel::K_Ic_of:180", "alpha_cap_theta": "Grimsel::alpha_of:152, AlphaAtTheta:276, BuildTables (strain table)"}
CALL = {"cp.a": "UpdatePropertyFields:4558, table build :7067, init :7229", "lambda.a": "UpdatePropertyFields:4554",
        "rho": "UpdatePropertyFields:4557, table build :7067, init :7228", "E0": "UpdatePropertyFields:4555-4556, LawSigma (via RockLaw)",
        "nu": "UpdatePropertyFields:4556, LawSigma, ElasticConstants", "K_Ic0": "ScanKIcPa Removal.H:114, UpdateSpAfterMechanics:6049 (kic.law = grimsel_table)",
        "alpha_cap_theta": "UpdatePropertyFields:4561-4563 (strain_of/alpha_of), LawSigma, BuckleRatio"}
for nm, rd in GR.items():
    caller = CALL.get(nm, CALL.get(nm.split(".")[0] + ".a", ""))
    a("microstructure.phase<i>.grimsel." + nm, "PROD",
      G_MS + "; thermal half needs phase<i>.law = grimsel, criterion half needs spall.stress_law = grimsel (RockLaw:6614)" + ("; kic.law = grimsel_table" if nm == "K_Ic0" else ""),
      "set (rock file) both legs; read every step", "reader in Grimsel.H: " + rd + "; called from " + caller + ". Strict queryclass scope (MMWSpalling.H:372), parsed for every phase whatever the law")
a("microstructure.phase<i>.grimsel.T_ref", "DEAD", "-", "not set; never read",
  "Grimsel.H:115 parses Material::T_ref, but Grimsel's methods never read it and the only consumer (Table::Build, Table.H:57) receives a COPY whose T_ref is overwritten by phase<i>.T_ref at MMWSpalling.H:7060 before :7067. Set by no input anywhere")
a("conductivity.kappa_damage_alpha", "DAMAGE", G_MS + " (AdvanceMicrostructure:1857, JetEnthalpyMarch:2500, PatchDiagnostics, BuildFaceList:3711, EnergyLedger, RefreshKappaEff, RecedeAlongLines FR:453)",
  "set 0.0 both legs; READ every step (exp(-0*D) = 1: inert)", "damage-modified conductivity (Step 6/8)")
a("conductivity.h_gb0", "MICRO", G_MS + " (AdvanceMicrostructure:1861 grain-boundary conductance)", "set 0.0 both legs; read every step, inert at 0", "")
a("damage.ic.type", "DAMAGE", G_MS + " (InitializeMicrostructure reads ic_damage)", "not set", "select<IC::Constant, IC::Expression>: also damage.ic.constant.value / damage.ic.expression.* (grouped rows)")
a("surface.follow_mask", "PROD", G_MS, "set 1 both legs; read AdvanceMicrostructure:1827, EnergyLedger:4069", "required by robin_form = surface and face_list = 1")
a("energy_ledger.enabled", "DIAG", G_MS + " (AdvanceMicrostructure:2081 calls EnergyLedger; :1829 also triggers the k_top gather, already on under follow_mask)",
  "set 1 both legs; runs every step", "thermo witnesses only; also hosts the side_power CSV (EnergyLedger:4108)")
for k in ["surface.refine_depth", "surface.refine_height"]:
    a(k, "CLOSED", "amr.max_level > 0 & removal on (TagCellsForRefinement:4991)", "not set; reader returns at :4991 on max_level 0",
      "WEAKEST VERDICT IN THIS FILE, unconfirmed by C.7.3/ARCHIVE headings. AMR surface band. Evidence: studies/amr_band/RESULTS.md (G1/G2 PASS, wall-time ratio 1.00 - no speed gain, 'speed is the next item'); the face list refuses amr.max_level > 0 (MMWSpalling.H:1042); set by NO input or test (only amr_band/run.py argv). Not in C.7.3 - the weakest CLOSED verdict here")
for nm in ["P0", "omega0", "lambda", "x0", "y0", "z0", "R", "profile", "radius", "alpha_expr", "schedule.t_end", "schedule.P"]:
    a("beam." + nm, "BEAM", G_BEAM + (" (beam.power_at, MMWBeam.H:93-98)" if nm.startswith(("P0", "schedule")) else ""),
      "set on both legs (P0 = 0, so inert)" if not nm.startswith("schedule") else "not set",
      ("also the default of bit.%s (MMWSpalling.H:5873-5875)" % nm if nm in ("x0", "y0", "z0") else "")
      + ("also the default vapor.A_beam (MMWSpalling.H:1119)" if nm == "omega0" else "")
      + ("schedule only in zhang/zhang_oglesby" if nm.startswith("schedule") else ""))
a("damage.enabled", "DAMAGE", "none (read Regrid, UpdateModel, InitializeMicrostructure, UpdateDamageAfterMechanics:5164)", "set 0 (= default) both legs", "sp_weibull aborts with it on (:9172-9175)")
for k in ["damage.phi_deg", "damage.cohesion", "damage.psi0", "damage.n", "damage.m", "damage.gb_multiplier"]:
    a(k, "DAMAGE", "damage.enabled (UpdateDamageAfterMechanics:5164)", "not set", "")
a("damage.stiffness_floor", "DAMAGE", "damage.enabled & " + G_MECH + " (UpdateModel:1292)", "not set", "")
a("damage.prescribed_stress.enabled", "DAMAGE", "damage.enabled", "not set", "")
for c in ["xx", "yy", "zz", "xy", "yz", "zx"]:
    a("damage.prescribed_stress.%s_expr" % c, "DAMAGE", "damage.enabled & damage.prescribed_stress.enabled", "not set", "")
for k in ["cohesive.enabled", "cohesive.K_n", "cohesive.delta_c", "cohesive.delta_max", "cohesive.G_c", "cohesive.f_t",
          "cohesive.prescribed_delta.enabled", "cohesive.prescribed_delta.expr"]:
    a(k, "DAMAGE", "cohesive.enabled (UpdateCohesiveAfterMechanics:5328)", "not set",
      ("G_c/f_t are read only when a delta is missing (:5293); set by NO input (gb_cohesive gives both deltas)" if k in ("cohesive.G_c", "cohesive.f_t") else "only unit/gb_cohesive sets it"))
a("spall.enabled", "PROD", "none", "set 1 both legs", "")
a("spall.stress_source", "PROD", G_SPW + " (UpdateSpAfterMechanics:6016); Removal.H:732 column path", "set local_thermoelastic both legs", "face list requires local_thermoelastic (:1039-1041)")
a("spall.stress_law", "PROD", "none (RockLaw:6614, every criterion site)", "set grimsel both legs", "")
a("spall.surface_normal", "CLOSED", G_SPW + " (UpdateSpAfterMechanics:6086); Removal.H:554 column path", "set 0 (= default) both legs; read at default",
  "Evidence: CODE_CENSUS C.7.3 row 'spall.surface_normal' (Step 20; recession half superseded by whole-cell removal, D2z-1 ARCHIVE L13103) and its own parse Message :5494-5501; the face list refuses it (:1038). ON THE TREE GATE: unit/beam_void_closure and meier/input_2d_dev set 1")
a("spall.melt_aware", "BEAM", G_COL + " & " + G_SPW + " (Removal.H:1027,1373,1473)", "set 0 (= default) both legs; never read on the face list", "P1 melt skeleton (ARCHIVE L2771); face list refuses it (:1037, 'the MMW step'). Tree gate: robin_pinned/jet/feet/jet_enthalpy inputs set 1")
a("melt.lambda_crit", "BEAM", "spall.melt_aware & " + G_COL + " (MeltSkeleton:8614 via Removal.H:1477)", "not set", "")
a("spall.clusters_csv", "DIAG", G_SPW + " & weibull.enabled (UpdateSpClusters:6535)", "not set", "")
a("spall.h_col_events_csv", "DIAG", G_COL + " & " + G_SPW + " (Removal.H:1851)", "not set", "")
a("spall.removal_events_csv", "DIAG", "removal on; face list: RecedeAlongLines FR:864; column: Removal.H:772,2447", "set 1 both legs; written every removal step", "parse requires sp_weibull + spall.enabled (:780-786)")
a("spall.surface_T", "PROD", G_SPW + " (BuildSurfaceTColumns:2868, UpdateSpAfterMechanics:6096,6161); Removal.H:807-2038 column path", "set reconstructed both legs", "face list requires reconstructed (:1047-1049)")
a("spall.flake_c_f", "PROD", "spall.surface_T != cell_average (BuildSurfaceTColumns) / " + G_FL, "set 0.4 both legs", "")
a("spall.flake_h_max", "PROD", "spall.surface_T != cell_average / " + G_FL + " (FR:368-369 gather window, FR:505,536)", "set 0.0015 both legs",
  "⚠ parse check :5610-5612 compares it with spall_h_spall_max BEFORE spall.h_spall_max is parsed (:5752) -> it sees the member default 0.010, not the deck's value")
a("spall.flake_length", "CLOSED", "BuildSurfaceTColumns:2975 (recon on); Removal.H:2488,2540 column path", "set twocell (= default) both legs; read at default",
  "Evidence: ARCHIVE L18497 'P1g-3 ... REJECTED ON PROVENANCE' and C.7.3 row; face list refuses flux (:1052-1054). Set by no unit/validation input")
a("spall.flake_fraction", "CLOSED", "BuildSurfaceTColumns:3016; Removal.H:2488,2540 column path", "set cf (= default) both legs; read at default",
  "Evidence: ARCHIVE L18497 (P1g-3 REJECTED) and C.7.3 row; face list refuses weibull (:1052-1054). Set by no unit/validation input")
a("spall.removal_energy", "PROD", G_COL + " (Removal.H:985,2038,2343,2382); on the face list NO runtime reader", "set flake both legs; enforced only by the parse interlock :1050-1051 (face list requires flake)",
  "on the production path this key is a parse-time requirement only")
for k in ["spall.damage_threshold", "spall.xi_s"]:
    a(k, "DAMAGE", "spall.enabled (copied at Removal.H:291-292 every removal step); used only by the damage-law column path", "not set; copied, unused on the face list", "")
for k in ["spall.C_h", "spall.alpha", "spall.t_spall"]:
    a(k, "DAMAGE", "spall.enabled (Removal.H:294-305: h_spall computed and abort-checked EVERY removal step, before the face-list dispatch)",
      "not set; h_spall computed every step on production and discarded", "damage-law h_spall = C_h sqrt(alpha t)")
a("spall.prescribed_ratio.enabled", "DAMAGE", "spall.enabled (Removal.H:309) -> damage-law column path", "not set", "")
a("spall.prescribed_ratio.expr", "DAMAGE", "spall.enabled & prescribed_ratio.enabled", "not set", "")
a("spall.sample", "DAMAGE", G_COL + " (Removal.H:566)", "set top_cell (= default) both legs; never read on the face list", "model-form comment says 'column path only; inert on the face list'")
a("spall.shallow_depth", "DAMAGE", G_COL + " & spall.sample = shallow_max", "not set", "Hu end-to-end only")
a("spall.h_spall_max", "PROD", G_FL + " (FR:369,378,505,535) & BuildSurfaceTColumns", "set 0.015 both legs", "")
a("spall.h_spall_n_segments", "PROD", G_FL + " (FR:473) / column (Removal.H:666)", "not set; READ on production at default 0 = defer to spallation.sp.n_segments",
  "set by NO input or test anywhere: deletion = hard-code the fallback")
a("spall.flake_coherence_length", "CLOSED", G_COL + " (Removal.H:1808)", "set 0.0 (= default) both legs; never read on the face list",
  "Evidence: CODE_CENSUS C.7.3 row 'superseded for the feet line (set to 0 there)'; face list refuses > 0 (:1055-1056). ON THE TREE GATE: beam_void_closure and meier/input_2d_dev set 0.004")
a("spall.seed_void", "DIAG", "removal on (Initialize:1151 -> ApplySeedVoid:6816)", "not set", "debug-only seeded interior void (D2q-2b, comment :5768-5779); unit/seed_void, unit/face_list")
for k in ["vapor.enabled", "vapor.A_beam", "vapor.T0", "vapor.prescribed_rop.enabled", "vapor.prescribed_rop.expr"]:
    a(k, "BEAM", ("none" if k == "vapor.enabled" else "vapor.enabled") + " (Removal.H:324-340; InitPhi:6876)",
      "set 0 (= default) both legs" if k == "vapor.enabled" else "not set",
      "face list refuses vapor (:1033)" + ("; prescribed_rop set by NO input or test" if "prescribed_rop" in k else ""))
for k in ["bit.enabled", "bit.feed_rate", "bit.radius", "bit.x0", "bit.y0", "bit.z0", "bit.start_time"]:
    a(k, "BEAM", "bit.enabled & " + G_COL + " (Removal.H:1903-1971)", "set 0 (= default) both legs" if k == "bit.enabled" else "not set",
      "Step 19 drawworks bit (regime 3); face list refuses it (:1034); used by meier/input_drilling(_25mm)")
a("spallation.model", "PROD", "none", "set sp_weibull both legs", "")
a("confining.p", "PROD", "ElasticConstants (Removal.H:44-45) -> every criterion site", "set 0.0 both legs (Kant pressure legs vary it)", "parsed before the sp_weibull early return (:9151)")
a("spall.buckle_span", "PROD", G_FL + " (RecedeAlongLines FR:346,488,623); parse requires face_list & surface_T = reconstructed (:947-956)", "HOS_R 0.012; CAL 0.0 (rock file 0.012 overridden later in the Kant deck)", "P2-15/P2-17")
a("spallation.sp.a0", "PROD", G_SPW + "; read only when weibull.enabled = 0 (FR:452, Removal.H:1125,1167, UpdateSpAfterMechanics:6042)",
  "set 8.0e-6 both legs; REQUIRED (:9177) but inert under weibull.enabled = 1", "rock file: 'inert under Weibull'")
a("spallation.sp.n_segments", "PROD", G_SPW, "set 128 both legs", "")
a("spallation.sp.stress_relief", "PROD", G_SPW + " (UpdateSpAfterMechanics:6019, ScanSigmaT:139,143)", "set 0.531294 both legs", "CAL-ONCE")
a("spallation.sp.kic.law", "PROD", G_SPW, "set grimsel_table both legs", "")
a("spallation.sp.kic.scale", "PROD", G_SPW + " & kic.law in {nasseri_table, grimsel_table}", "set 1.0 (= default) both legs", "")
for k in ["spallation.sp.kic.K_Ic0", "spallation.sp.kic.T_ref", "spallation.sp.kic.T_melt"]:
    a(k, "MECH", G_SPW + " & kic.law = linear (ScanKIcPa:120-124, UpdateSpAfterMechanics:6050-6052)", "not set (silently ignored under grimsel_table: parsed only under linear)",
      "the Step 12 linear K_Ic law: kant/sp_kant_onset/input_verification(_hfl), sp_kant_pressure_sweep only")
for k in ["spallation.sp.prescribed_stress.enabled", "spallation.sp.prescribed_stress.xx_expr",
          "spallation.sp.prescribed_T_surface.enabled", "spallation.sp.prescribed_T_surface.expr"]:
    a(k, "MECH", G_SPW + " (UpdateSpAfterMechanics:6006-6009)" + ("" if k.endswith("enabled") else " & its .enabled"), "not set",
      "criterion-site-A test inputs (unit/sp_lefm)")
a("weibull.enabled", "PROD", G_SPW, "set 1 both legs", "")
a("weibull.m", "PROD", G_SPW + " (parsed outside the weibull.enabled guard since P1g-3, :9288)", "set 13 both legs", "")
a("weibull.a0_gb", "MICRO", G_SPW + " & weibull.enabled (InitializeFlaws:7418, grain-boundary cells only)", "set 40.0e-6 both legs; REQUIRED (:9294) but inert: rock is grain-free", "")
a("weibull.a0_ig", "PROD", G_SPW + " & weibull.enabled", "set 8.0e-6 both legs (CAL-ONCE)", "")
a("weibull.V0", "PROD", G_SPW + " & weibull.enabled (InitializeFlaws:7451)", "set 1.0e-9 both legs (= V_cell, factor 1)", "")
a("weibull.size_volume", "CLOSED", G_SPW + " & weibull.enabled", "set cell (= default) both legs; read at default",
  "Evidence: ARCHIVE L19170 P2-2d ('size_volume = flux stays default-off ... unvalidated, per-ladder anchor, no target', L19252); face list refuses flux (:1052-1054); the only non-default setter is the lateral_flux_lock abort leg")
a("weibull.V_ref", "CLOSED", G_SPW + " & weibull.enabled & size_volume = flux (WeibullSizeFactor:7757)", "not set; reader returns 1 before reading it",
  "Evidence: as weibull.size_volume (P2-2d L19250). Set by NO input")
a("weibull.seed", "PROD", G_SPW + " & weibull.enabled", "set 12345 both legs", "")
a("weibull.detach_mode", "PROD", G_SPW + " & weibull.enabled (UpdateSpClusters:6465, InitializeFlaws:7458)", "set per_face (= default) both legs", "")
a("weibull.A_crit", "MECH", G_SPW + " & weibull.enabled & detach_mode = connected_cluster (UpdateSpClusters:6466-6473; column detach Removal.H:1162)",
  "not set", "Step 15c cluster-area detach rule; kant/sp_kant_onset/input_weibull_cluster only; face list refuses it (:1057-1061)")
a("surface_patch.enabled", "PROD", "none", "HOS_R 1 (input_lab 0, argv 1); CAL 1", "strict scope: enabled = 0 makes every other surface_patch.* key an abort")
a("surface_patch.face_list", "PROD", "surface_patch.enabled", "set 1 both legs", "")
a("surface_patch.face_heat", "FLAME-CAP", "face_list = 1 (FaceListTg:2246, JetEnthalpyMarch:2639, FaceListHeat:3858)", "not set (position)",
  "X1 surface-path heat rule: P2-11b accepted (ARCHIVE L24468); deleted and RESTORED by the user in P2-12 (L24861) - NOT closed. Only unit/face_heat sets it")
a("surface_patch.barrier_expr", "FLAME-CAP", "face_heat = path (FaceListPath:3411,3421)", "not set", "as face_heat (P2-11b/P2-12)")
a("surface_patch.type", "DEAD", "-", "not set", "consumed :7940, compared with 'prescribed_temperature' under prescribed_T only (:7952-7955); no reader. STRICT scope: 17 Hu decks set it, so removing the parse turns those decks into aborts - delete the deck lines in the same step (CODE_CENSUS C.7.4)")
a("surface_patch.mode", "PROD", "surface_patch.enabled", "set convective_flame both legs", "")
a("surface_patch.robin_form", "PROD", "surface_patch.enabled", "set surface both legs", "")
a("surface_patch.pinned_idle_cycles", "FLAME-CAP", "robin_form = pinned (PinRule:3111 via BuildPinEffective:3092 and JetEnthalpyMarch:2605)",
  "set 1.0e9 both legs; INERT under robin_form = surface (Parse Message :770-774 says so)",
  "C1 pinned closure (ARCHIVE L4422), the Meier feet decks' form (input_feet*, pinned); production moved to robin_form = surface at P0b (ARCHIVE L14061; C.7.1 'surface on every study since P0b'). The brief lists 'pinned idle' under FLAME-CAP. Parsed under pinned AND surface (strict scope :7976-7993), so the prod decks' 1.0e9 must go with the parse. TREE GATE exercises the pinned form (robin_pinned, robin_jet, robin_feet, jet_enthalpy)")
a("surface_patch.pad_shield_open", "CLOSED", "surface_patch.enabled & h_expr (BuildFlameColumns:2116, FaceHr:3686); effective only with feet geometry (CheckPadRing:7889)",
  "not set", "Evidence: ARCHIVE L13601 'D2aa ... Wall lever measured and exhausted'; C.7.3 row 'D2aa bounded; P1 confirmed, P2 refuted'. Only setter: unit/face_heat R_pad abort leg")
a("surface_patch.face", "DEAD", "-", "set zhi (= default, the only legal value) both legs", "consumed :8014 and checked == 'zhi'; no reader. STRICT scope: 2 prod decks + 30 validation + 12 unit decks set it -> delete deck lines with the parse")
for k in ["surface_patch.x0", "surface_patch.y0", "surface_patch.radius"]:
    a(k, "PROD", "surface_patch.enabled", "set both legs", "REQUIRED")
a("surface_patch.T_expr", "ZHANG", "mode = prescribed_T (UpdateModel:1240, AdvanceMaterial:1554, AdvanceMicrostructure:1878, EnergyLedger:4090, Removal.H:685 column)",
  "not set", "Hu prescribed-T patch (Step 6d), the global-material path. ⚠ consumed and IGNORED under convective_flame (:8034): silently inert")
a("surface_patch.T_ext", "DEAD", "-", "not set", "parsed :8041 into SurfacePatch::T_ext; zero readers. STRICT scope: 27 validation (Hu) + 4 unit decks set it")
a("surface_patch.h_expr", "PROD", "surface_patch.enabled & mode = convective_flame", "set both legs (HOS_R map; CAL 500.0)", "")
a("surface_patch.T_flame_expr", "PROD", "h_expr given & jet_closure = none (BuildFlameColumns:2196)", "CAL 1673.15; HOS_R absent (forbidden under enthalpy, :8066-8068)", "")
a("surface_patch.jet_closure", "PROD", "h_expr given", "HOS_R enthalpy; CAL none", "")
for k, v in [("jet_mdot", "3.4e-3"), ("jet_T_nozzle", "1866.0"), ("jet_cp", "1250.0"), ("jet_D", "0.0079"), ("jet_T_ent", "293.15 (= default)"),
             ("jet_core_length", "30.0"), ("jet_stagnation", "decay (= default)"), ("jet_entrained_mass", "1"), ("jet_T_ent_mode", "fixed (= default)"),
             ("jet_decay_diameter", "momentum"), ("jet_De_ref", "0.003745"), ("jet_bin_update", "exponential"), ("jet_negative_flux", "count")]:
    a("surface_patch." + k, "PROD", G_JET, "HOS_R " + v + "; CAL n/a (jet_closure = none)", "")
A["surface_patch.jet_T_ent_mode"] = ("PROD", G_JET, "HOS_R fixed (= default); CAL n/a",
    "values: fixed (prod) | exhaust (feet line: robin_feet, meier/input_feet*) | rec_fixed (D2e diagnostic, see jet_T_rec_fixed)")
a("surface_patch.jet_T_rec_fixed", "CLOSED", G_JET + " & jet_T_ent_mode = rec_fixed (JetEnthalpyMarch:2571)", "not set",
  "Evidence: ARCHIVE L8974 D2e and C.7.3 row 'jet_T_ent_mode = rec_fixed: D2e diagnostic'. Only setter: unit/jet_d2c (KT/K15 legs - ON THE TREE GATE)")
for k in ["jet_free_surface", "jet_fs_aspect", "jet_fs_exponent"]:
    a("surface_patch." + k, "FLAME-CAP", G_JET + " & nozzle_descent = feet" + ("" if k == "jet_free_surface" else " & jet_free_surface = 1"), "not set",
      "D2c steady closures (ARCHIVE L7335); meier/input_feet_d2c; face_heat = path refuses it (:987). TREE GATE: jet_d2c")
a("surface_patch.jet_above_h", "CLOSED", "h_expr given (BuildFlameColumns:2156, FaceHr:3678)", "not set",
  "Evidence: ARCHIVE L12600 D2v (the funnel was the map) and C.7.3 row; ARCHIVE L18980 'above-plane h is hard-zeroed unless jet_above_h ... set on M2 only'. Only setter: unit/face_heat R_above abort leg")
a("surface_patch.jet_T_above", "CLOSED", "jet_above_h = 1 (BuildFlameColumns:2190-2191)", "not set", "as jet_above_h (D2v); pair is both-or-neither (:8189-8196)")
a("surface_patch.jet_dr", "FLAME-CAP", G_JET + " (JetEnthalpyMarch:2509); face_heat = path (FaceListPath:3364)", "not set; READ on HOS_R at default 0 = max(dx, dy)", "march bin width; set only by jet_enthalpy / jet_d2c tests")
a("surface_patch.jet_profile_interval", "DIAG", G_JET + " (JetEnthalpyMarch:2677)", "HOS_R 1.0 (writes _jet_profile.csv)", "")
a("surface_patch.side_power_interval", "DIAG", G_JET + " (parse branch) & energy_ledger.enabled & robin_form = surface (EnergyLedger:4108,4250)", "HOS_R 10.0 (writes _side_power.csv)",
  "rev7 d2 side-power CSV (the brief lists side power under FLAME-CAP; it is a CSV-interval switch, so DIAG); survives the P2-10 side-face deletion. ⚠ parsed inside the enthalpy branch, so a jet_closure = none deck cannot set it (strict abort)")
a("surface_patch.nozzle_descent", "PROD", "h_expr given (BuildFlameColumns:2122); Removal.H:431 (feet & foot_body_clearance)", "set prescribed (= default) both legs", "selects prescribed vs feet")
a("surface_patch.foot_body_clearance", "FLAME-CAP", "feet (FeetDescent:2316; Removal.H:431-432 runs BEFORE the face-list dispatch)", "not set",
  "D2c/P2-2f (ARCHIVE L7335, L19275); face list refuses it (regime 4, :1035-1036). TREE GATE: jet_d2c")
for k in ["h_conv", "T_flame"]:
    a("surface_patch." + k, "FLAME-CAP", "surface_patch.enabled & no h_expr (BuildFlameColumns:2108-2109 constant flame; AdvanceMaterial:1555-1556)" + ("; T_flame also SurfaceGasTempCap:9121" if k == "T_flame" else ""),
      "not set (forbidden with h_expr, :8259)", "constant-coefficient flame (Step 15a); TREE GATE: robin_face, robin_pinned, s1 smoke; Hu decks")
a("surface_patch.nozzle_z0", "PROD", "h_expr & nozzle_descent = prescribed (BuildFlameColumns:2124)", "HOS_R 0.2132; CAL 0.1272", "REQUIRED under prescribed")
a("surface_patch.nozzle_feed", "PROD", "h_expr & nozzle_descent = prescribed (BuildFlameColumns:2124)", "set 0.0 (= default) both legs", "")
for k in ["foot_r_inner", "foot_r_outer", "foot_standoff", "foot_rule", "foot_npads", "foot_pad_quantile", "nozzle_feed_max"]:
    a("surface_patch." + k, "FLAME-CAP", G_FEET + ("; also BuildFlameColumns:2182,2220 / FaceHr:3687 under pad_shield_open >= 0; foot_r_outer also JetEnthalpyMarch:2544 under jet_free_surface" if k in ("foot_r_inner", "foot_r_outer") else ""),
      "not set", "D2b feet descent (ARCHIVE L6390); TREE GATE: robin_feet")
a("surface_patch.nozzle_collision_radius", "PROD", "h_expr given (BuildFlameColumns:2120-2121)", "set 0.004 both legs", "standing rule 6")
