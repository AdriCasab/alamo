# The key table, written by hand from the parse sites (MMWSpalling.H unless noted).
# kind: m = class member (bare name), sp = SurfacePatch member, ph = MineralPhase member,
#       beam = MMWBeam member, grim = Grimsel member, none = no stored member (consumed only)
# (key, parse lines, default, members, kind, parse_guard)
M = "MMWSpalling.H"
K = []
def k(key, lines, default, members, kind="m", guard=""):
    K.append(dict(key=key, lines=lines, default=default, members=members, kind=kind, guard=guard))

# ---- retired (abort) -------------------------------------------------------
P210 = ["spall.lateral_drive","spall.lateral_removal","spall.lateral_subcell","spall.lateral_thermal",
 "spall.per_cell_removal","surface_patch.side_face_flux","surface_patch.side_face_h",
 "surface_patch.side_face_T","surface_patch.side_face_h_r_min","surface_patch.side_face_h_law",
 "surface_patch.side_face_h_datum","surface_patch.side_face_factor_above",
 "surface_patch.side_face_area_law","surface_patch.side_face_step_height",
 "surface_patch.side_face_staircase","surface_patch.step_diag_radius","surface_patch.jet_wall_march",
 "surface_patch.jet_wall_mdot","surface_patch.jet_body","surface_patch.jet_body_T",
 "surface_patch.jet_body_h_scale"]
for i, key in enumerate(P210):
    k(key, ["%s:%d (list :97-118), :119 contains -> abort :120-122" % (M, 98 + i)], "ABORT", [], "retired")
k("spallation.sp.constraint_depth", [M + ":127 (list), :128 contains -> abort :129-132"], "ABORT", [], "retired")

# ---- MMWSpalling::Parse body ------------------------------------------------
k("hc.heat.alpha", [M+":153 (also HeatConduction via queryclass hc :136)"], "REQUIRED", ["therm_alpha"])
k("phase.T_ref", [M+":157"], "0.0", ["T_ref"])
k("phase.T_m", [M+":158"], "1.0e30", ["T_m"])
k("phase.T_v", [M+":159"], "1.0e30", ["T_v"])
k("phase.lh_m", [M+":160"], "0.0", ["lh_m"])
k("phase.lh_v", [M+":161"], "0.0", ["lh_v"])
k("losses.epsilon", [M+":165"], "0.0", ["eps_rad"])
k("losses.sigma_SB", [M+":166"], "0.0", ["sigma_SB"])
k("losses.h_conv", [M+":167"], "0.0", ["h_conv"])
k("losses.T_amb", [M+":168"], "0.0", ["T_amb"])
k("losses.epsilon_high", [M+":172"], "=losses.epsilon", ["eps_rad_hi"])
k("losses.T_eps_break", [M+":173"], "1.0e30", ["T_eps_break"])
k("material.type", [M+":179 contains, :180 select_default<Zhang,Constant>"], "absent (zhang if block present)", ["material"], guard="")
k("solver", [M+":260"], "explicit", ["solver_kind"])
k("implicit.tol_rel", [M+":264"], "1.0e-10", ["implicit_tol_rel"])
k("implicit.tol_abs", [M+":265"], "0.0", ["implicit_tol_abs"])
k("implicit.max_iter", [M+":266"], "100", ["implicit_max_iter"])
k("implicit.verbose", [M+":267"], "0", ["implicit_verbose"])
k("alpha", [M+":270 pp_queryarr"], "absent (empty)", ["alpha"])
G_MS = "microstructure.nphases present (:292)"
k("microstructure.nphases", [M+":292 contains", M+":332"], "absent (=> microstructure off)", ["phases", "microstructure_enabled"])
k("microstructure.mode", [M+":299"], "voronoi", ["micro_mode"], guard=G_MS)
k("microstructure.number_of_grains", [M+":306"], "REQUIRED if mode=voronoi", ["num_grains"], guard=G_MS+" & mode=voronoi")
k("microstructure.seed", [M+":307"], "1", ["micro_seed"], guard=G_MS+" & mode=voronoi")
k("microstructure.phase_expr", [M+":312"], "REQUIRED if mode=expression", ["phase_parser","phase_f"], guard=G_MS+" & mode=expression")
k("microstructure.grain_expr", [M+":321 contains", M+":323"], "absent", ["grain_parser","grain_f","has_grain_expr"], guard=G_MS+" & mode=expression")
PH = "microstructure.phase<i>"
k(PH+".name", [M+":338"], "phase<i>", ["name"], "ph", G_MS)
for j,(nm, mem) in enumerate([("fraction","fraction"),("kappa","kappa"),("beta","beta"),("E","E"),("mu","mu"),
                              ("rho","rho"),("Cp","Cp")]):
    k(PH+"."+nm, [M+":%d" % (340+j)], "REQUIRED", [mem], "ph", G_MS)
k(PH+".alpha_attenuation", [M+":347"], "REQUIRED", ["alpha_atten"], "ph", G_MS)
k(PH+".A_D", [M+":349"], "REQUIRED", ["A_D"], "ph", G_MS)
k(PH+".law", [M+":356"], "constant", ["grimsel_law"], "ph", G_MS)
for j,(nm,mem,d) in enumerate([("T_ref","T_ref","273.15"),("T_m","T_m","1.0e30"),("fus_width","fus_width","0.0"),
                                ("L_m","L_m","0.0"),("T_vap_lo","T_vap_lo","1.0e30"),("T_vap_hi","T_vap_hi","1.0e30"),
                                ("L_v","L_v","0.0")]):
    k(PH+"."+nm, [M+":%d" % (387+j)], d, [mem], "ph", G_MS)
k(PH+".alpha_beta.enabled", [M+":400"], "0", ["alpha_beta_enabled","alpha_beta_linear_strain"], "ph", G_MS)
k(PH+".alpha_beta.T0", [M+":401"], "846.15", ["alpha_beta_T0"], "ph", G_MS)
k(PH+".alpha_beta.width", [M+":403"], "10.0", ["alpha_beta_width"], "ph", G_MS)
k(PH+".alpha_beta.dV", [M+":405"], "0.0063", ["alpha_beta_dV", "alpha_beta_linear_strain"], "ph", G_MS)
k(PH+".alpha_beta.linear_strain", [M+":411"], "cbrt(1+dV)-1", ["alpha_beta_linear_strain"], "ph", G_MS)
GR = "Numeric/Material/Grimsel.H"
for (nm, ln, d, mem) in [("cp.a",96,"1370.0","cp_a"),("cp.b",97,"178000.0","cp_b"),("cp.th0",98,"271.0","cp_th0"),
                         ("lambda.a",100,"2000.0","lam_a"),("lambda.th0",101,"563.0","lam_th0"),
                         ("rho",103,"2650.0","rho_0"),("E0",105,"53.3e9","E_0"),("nu",106,"0.30","nu"),
                         ("K_Ic0",109,"1.5","K_Ic0"),("alpha_cap_theta",114,"-1.0","alpha_cap_theta"),
                         ("T_ref",115,"293.15","T_ref")]:
    k(PH+".grimsel."+nm, ["%s:%d (queryclass at %s:372, every phase, unconditional)" % (GR, ln, M)], d, [mem], "grim", G_MS)
k("conductivity.kappa_damage_alpha", [M+":489"], "0.0", ["kappa_damage_alpha"], guard=G_MS)
k("conductivity.h_gb0", [M+":491"], "0.0", ["h_gb_conductance"], guard=G_MS)
k("damage.ic.type", [M+":495 contains, :496 select<IC::Constant,IC::Expression>"], "absent (D=0)", ["ic_damage"], guard=G_MS)
k("surface.follow_mask", [M+":659"], "0", ["surface_follow_mask"])
k("energy_ledger.enabled", [M+":660"], "0", ["ledger_enabled"])
k("surface.refine_depth", [M+":687"], "0.0", ["surface_refine_depth"])
k("surface.refine_height", [M+":688"], "0.0", ["surface_refine_height"])
# ---- beam (queryclass "beam" :138, Numeric/MMWBeam.H) ----------------------
BM = "Numeric/MMWBeam.H"
for (nm, ln, d, mem) in [("P0",39,"0.0","P0"),("omega0",41,"1.0","omega0"),("lambda",43,"1.0e-3","lambda"),
                         ("x0",45,"0.0","x0"),("y0",46,"0.0","y0"),("z0",47,"0.0","z0"),("R",52,"0.0","R"),
                         ("profile",61,"gaussian","profile"),("radius",66,"0.0","radius"),
                         ("alpha_expr",74,"10","alpha_f"),("schedule.t_end",83,"empty","sched_t"),
                         ("schedule.P",84,"empty","sched_P")]:
    k("beam."+nm, ["%s:%d (queryclass %s:138)" % (BM, ln, M)], d, [mem], "beam")
# ---- ParseDamageSettings ------------------------------------------------------
k("damage.enabled", [M+":5016"], "0", ["damage_enabled"])
k("damage.phi_deg", [M+":5021 (required if enabled)", M+":5026 (default otherwise)"], "0.0 / REQUIRED if damage.enabled", ["damage_phi_deg","damage_phi_rad","damage_alpha_dp","damage_k_dp"])
k("damage.cohesion", [M+":5022", M+":5027"], "1.0 / REQUIRED if damage.enabled", ["damage_cohesion","damage_k_dp"])
k("damage.psi0", [M+":5029"], "0.0", ["damage_psi0"])
k("damage.n", [M+":5030"], "1.0", ["damage_n"])
k("damage.m", [M+":5031"], "1.0", ["damage_m"])
k("damage.gb_multiplier", [M+":5032"], "1.0", ["damage_gb_multiplier"])
k("damage.stiffness_floor", [M+":5033"], "1.0e-6", ["damage_stiffness_floor"])
k("damage.prescribed_stress.enabled", [M+":5036"], "0", ["damage_prescribed_stress_enabled"])
for j, c in enumerate(["xx","yy","zz","xy","yz","zx"]):
    k("damage.prescribed_stress.%s_expr" % c, [M+":5045 (lambda), call :%d" % (5051 + 2*j)], "0.0",
      ["damage_s%s_f" % c, "damage_s%s_parser" % c], guard="damage.prescribed_stress.enabled = 1")
# ---- ParseCohesiveSettings ----------------------------------------------------
GC = "cohesive.enabled = 1"
k("cohesive.enabled", [M+":5276"], "0", ["cohesive_enabled"])
k("cohesive.K_n", [M+":5280"], "REQUIRED if enabled", ["cohesive_K_n"], guard=GC)
k("cohesive.delta_c", [M+":5289 contains", M+":5291"], "f_t/K_n", ["cohesive_delta_c"], guard=GC)
k("cohesive.delta_max", [M+":5290 contains", M+":5292"], "2 G_c/f_t", ["cohesive_delta_max"], guard=GC)
k("cohesive.G_c", [M+":5296"], "REQUIRED unless both deltas", ["cohesive_delta_max"], guard=GC+" & a delta missing")
k("cohesive.f_t", [M+":5297"], "REQUIRED unless both deltas", ["cohesive_delta_c","cohesive_delta_max"], guard=GC+" & a delta missing")
k("cohesive.prescribed_delta.enabled", [M+":5314"], "0", ["cohesive_prescribed_enabled"], guard=GC)
k("cohesive.prescribed_delta.expr", [M+":5319"], "REQUIRED if enabled", ["cohesive_delta_f","cohesive_delta_parser"], guard=GC+" & prescribed_delta.enabled")
# ---- ParseSpallSettings -------------------------------------------------------
k("spall.enabled", [M+":5403"], "0", ["spall_enabled"])
k("spall.stress_source", [M+":5415"], "mechanics", ["sp_local_thermoelastic"])
k("spall.stress_law", [M+":5438"], "scalar", ["sp_stress_law_grimsel"])
k("spall.surface_normal", [M+":5474"], "0", ["spall_surface_normal"])
k("spall.melt_aware", [M+":5519"], "0", ["melt_aware"])
k("melt.lambda_crit", [M+":5521"], "0.4", ["melt_lambda_crit"])
k("spall.clusters_csv", [M+":5532"], "0", ["sp_emit_clusters_csv"])
k("spall.h_col_events_csv", [M+":5535"], "0", ["sp_emit_h_col_events_csv"])
k("spall.removal_events_csv", [M+":5540"], "0", ["sp_emit_removal_events_csv"])
k("spall.surface_T", [M+":5552"], "cell_average", ["sp_surface_T_mode","sp_surface_T_recon","sp_surface_T_apply"])
k("spall.flake_c_f", [M+":5599"], "0.4", ["sp_flake_c_f"])
k("spall.flake_h_max", [M+":5606"], "0.0", ["sp_flake_h_max"])
k("spall.flake_length", [M+":5626"], "twocell", ["sp_flake_len_flux"])
k("spall.flake_fraction", [M+":5632"], "cf", ["sp_flake_frac_weibull"])
k("spall.removal_energy", [M+":5640"], "cell", ["sp_removal_energy_flake"])
GS = "spall.enabled = 1 (early return :5681)"
k("spall.damage_threshold", [M+":5683"], "0.99", ["spall_damage_threshold"], guard=GS)
k("spall.xi_s", [M+":5684"], "8.0", ["spall_xi_s"], guard=GS)
k("spall.C_h", [M+":5685"], "1.0", ["spall_C_h"], guard=GS)
k("spall.alpha", [M+":5687 contains", M+":5690"], "absent (derived)", ["spall_alpha","spall_alpha_override"], guard=GS)
k("spall.t_spall", [M+":5692 contains", M+":5695"], "absent (derived)", ["spall_t_spall","spall_t_spall_override"], guard=GS)
k("spall.prescribed_ratio.enabled", [M+":5715"], "0", ["spall_prescribed_ratio_enabled"], guard=GS)
k("spall.prescribed_ratio.expr", [M+":5720"], "REQUIRED if enabled", ["spall_ratio_f","spall_ratio_parser"], guard=GS+" & prescribed_ratio.enabled")
k("spall.sample", [M+":5732"], "top_cell", ["spall_sample_mode"], guard=GS)
k("spall.shallow_depth", [M+":5743"], "0.006", ["spall_shallow_depth"], guard=GS)
k("spall.h_spall_max", [M+":5752"], "0.010", ["spall_h_spall_max"], guard=GS)
k("spall.h_spall_n_segments", [M+":5757"], "0", ["spall_h_spall_n_segments"], guard=GS)
k("spall.flake_coherence_length", [M+":5763"], "0.0", ["spall_flake_coh_len"], guard=GS)
k("spall.seed_void", [M+":5789 queryarr"], "absent", ["spall_seed_void"], guard=GS)
# ---- ParseVaporSettings --------------------------------------------------------
GV = "vapor.enabled = 1 (early return :5818)"
k("vapor.enabled", [M+":5816"], "0", ["vapor_enabled"])
k("vapor.A_beam", [M+":5820 contains", M+":5823"], "0.5 pi omega0^2", ["vapor_A_beam","vapor_A_beam_override"], guard=GV)
k("vapor.T0", [M+":5828 contains", M+":5831"], "losses.T_amb", ["vapor_T0","vapor_T0_override"], guard=GV)
k("vapor.prescribed_rop.enabled", [M+":5837"], "0", ["vapor_prescribed_rop_enabled"], guard=GV)
k("vapor.prescribed_rop.expr", [M+":5842"], "REQUIRED if enabled", ["vapor_rop_f","vapor_rop_parser"], guard=GV+" & prescribed_rop.enabled")
# ---- ParseBitSettings ------------------------------------------------------------
GB = "bit.enabled = 1 (early return :5859)"
k("bit.enabled", [M+":5857"], "0", ["bit_enabled"])
k("bit.feed_rate", [M+":5865"], "REQUIRED if enabled", ["bit_feed_rate"], guard=GB)
k("bit.radius", [M+":5869"], "REQUIRED if enabled", ["bit_radius"], guard=GB)
k("bit.x0", [M+":5873"], "beam.x0", ["bit_x0"], guard=GB)
k("bit.y0", [M+":5874"], "beam.y0", ["bit_y0"], guard=GB)
k("bit.z0", [M+":5875"], "beam.z0", ["bit_z0"], guard=GB)
k("bit.start_time", [M+":5876"], "0.0", ["bit_start_time"], guard=GB)
# ---- ParseSpallationModel ----------------------------------------------------------
k("spallation.model", [M+":9140"], "damage_law", ["spallation_model"])
k("confining.p", [M+":9151"], "0.0", ["confining_p"])
k("spall.buckle_span", [M+":9160"], "0.0", ["sp_buckle_span"])
GW = "spallation.model = sp_weibull (early return :9164)"
k("spallation.sp.a0", [M+":9177"], "REQUIRED under sp_weibull", ["sp_a0"], guard=GW)
k("spallation.sp.n_segments", [M+":9181"], "64", ["sp_n_segments"], guard=GW)
k("spallation.sp.stress_relief", [M+":9184"], "1.0", ["sp_relief"], guard=GW)
k("spallation.sp.kic.law", [M+":9197"], "nasseri_table", ["kic_law"], guard=GW)
k("spallation.sp.kic.scale", [M+":9224"], "1.0", ["kic_scale"], guard=GW+" & kic.law in {nasseri_table, grimsel_table}")
k("spallation.sp.kic.K_Ic0", [M+":9231"], "REQUIRED if kic.law=linear", ["kic_K_Ic0"], guard=GW+" & kic.law = linear")
k("spallation.sp.kic.T_ref", [M+":9232"], "REQUIRED if kic.law=linear", ["kic_T_ref"], guard=GW+" & kic.law = linear")
k("spallation.sp.kic.T_melt", [M+":9233"], "REQUIRED if kic.law=linear", ["kic_T_melt"], guard=GW+" & kic.law = linear")
k("spallation.sp.prescribed_stress.enabled", [M+":9243"], "0", ["sp_prescribed_sigma_enabled"], guard=GW)
k("spallation.sp.prescribed_stress.xx_expr", [M+":9248"], "REQUIRED if enabled", ["sp_prescribed_sigma_f","sp_prescribed_sigma_parser"], guard=GW+" & prescribed_stress.enabled")
k("spallation.sp.prescribed_T_surface.enabled", [M+":9255"], "0", ["sp_prescribed_T_surface_enabled"], guard=GW)
k("spallation.sp.prescribed_T_surface.expr", [M+":9260"], "REQUIRED if enabled", ["sp_prescribed_T_surface_f","sp_prescribed_T_surface_parser"], guard=GW+" & prescribed_T_surface.enabled")
k("weibull.enabled", [M+":9275"], "0", ["weibull_enabled"], guard=GW)
k("weibull.m", [M+":9288"], "15.0", ["weibull_m"], guard=GW)
GWE = GW + " & weibull.enabled = 1"
k("weibull.a0_gb", [M+":9294"], "REQUIRED if weibull.enabled", ["weibull_a0_gb"], guard=GWE)
k("weibull.a0_ig", [M+":9302"], "a0_gb/5", ["weibull_a0_ig"], guard=GWE)
k("weibull.V0", [M+":9315"], "1.0e-9", ["weibull_V0"], guard=GWE)
k("weibull.size_volume", [M+":9322"], "cell", ["sp_size_flux"], guard=GWE)
k("weibull.V_ref", [M+":9329"], "-1.0", ["weibull_V_ref"], guard=GWE)
k("weibull.seed", [M+":9343"], "0", ["weibull_seed"], guard=GWE)
k("weibull.detach_mode", [M+":9351"], "per_face", ["weibull_detach_mode"], guard=GWE)
k("weibull.A_crit", [M+":9363"], "REQUIRED if connected_cluster", ["weibull_A_crit"], guard=GWE+" & detach_mode = connected_cluster")
# ---- SurfacePatch::Parse (queryclass "surface_patch", strict scope) ---------------
SPE = "surface_patch.enabled = 1 (early return :7902)"
def s(nm, lines, d, mem, guard=SPE):
    k("surface_patch."+nm, [M+":"+x for x in lines], d, mem, "sp", guard)
s("enabled", ["7900"], "0", ["enabled"], guard="")
s("face_list", ["7906"], "0", ["face_list"])
s("face_heat", ["7912 contains","7918"], "absent (position)", ["face_heat_path"], SPE+" & face_list = 1")
s("barrier_expr", ["7924 contains","7930"], "absent", ["barrier_f","barrier_parser","has_barrier"], SPE+" & face_heat = path")
s("type", ["7940"], "prescribed_temperature", [])
s("mode", ["7944"], "prescribed_T", ["mode"])
s("robin_form", ["7962"], "cell", ["robin_form"])
s("pinned_idle_cycles", ["7993"], "2.0", ["pinned_idle_cycles"], SPE+" & robin_form in {pinned, surface}")
s("pad_shield_open", ["8004 contains","8005"], "-1.0", ["pad_shield_open"])
s("face", ["8014"], "zhi", [])
s("x0", ["8018"], "REQUIRED", ["x0"])
s("y0", ["8019"], "REQUIRED", ["y0"])
s("radius", ["8020"], "REQUIRED", ["radius"])
s("T_expr", ["8027 (required, prescribed_T)","8034 (default, otherwise)"], "REQUIRED if mode=prescribed_T, else ''", ["T_f","T_parser"])
s("T_ext", ["8041"], "0.0", ["T_ext"])
s("h_expr", ["8049 contains","8263"], "absent", ["h_f","h_parser","use_expr"])
s("T_flame_expr", ["8050 contains","8269"], "absent", ["Tg_f","Tg_parser"], SPE+" & h_expr given & jet_closure = none")
s("jet_closure", ["8056"], "none", ["jet_enthalpy"])
JE = SPE + " & jet_closure = enthalpy"
s("jet_mdot", ["8069"], "REQUIRED under enthalpy", ["jet_mdot"], JE)
s("jet_T_nozzle", ["8070"], "REQUIRED under enthalpy", ["jet_T_nozzle"], JE)
s("jet_cp", ["8071"], "REQUIRED under enthalpy", ["jet_cp"], JE)
s("jet_D", ["8072"], "REQUIRED under enthalpy", ["jet_D"], JE)
s("jet_T_ent", ["8073"], "293.15", ["jet_T_ent"], JE)
s("jet_core_length", ["8077 contains","8078"], "5.0", ["jet_core_length"], JE)
s("jet_stagnation", ["8080"], "decay", ["jet_stag_decay"], JE)
s("jet_entrained_mass", ["8086"], "0", ["jet_entrained_mass"], JE)
s("jet_T_ent_mode", ["8094"], "fixed", ["jet_T_ent_exhaust","jet_rec_fixed"], JE)
s("jet_T_rec_fixed", ["8110 contains","8113"], "REQUIRED if rec_fixed", ["jet_T_rec_fixed"], JE+" & jet_T_ent_mode = rec_fixed")
s("jet_decay_diameter", ["8126","8229 (else-branch abort)"], "nozzle", ["jet_De_momentum"], JE)
s("jet_De_ref", ["8139","8144 contains (abort)","8229"], "REQUIRED if momentum", ["jet_De_ref"], JE+" & jet_decay_diameter = momentum")
s("jet_free_surface", ["8150","8229"], "0", ["jet_free_surface"], JE)
s("jet_fs_aspect", ["8156","8163 contains (abort)","8229"], "1.0", ["jet_fs_aspect"], JE+" & jet_free_surface = 1")
s("jet_fs_exponent", ["8157","8163 contains (abort)","8229"], "1.0", ["jet_fs_exponent"], JE+" & jet_free_surface = 1")
s("jet_bin_update", ["8167","8229"], "explicit", ["jet_bin_exponential"], JE)
s("jet_negative_flux", ["8173","8229"], "drop", ["jet_neg_count"], JE)
s("jet_above_h", ["8182","8229"], "0", ["jet_above_h"], JE)
s("jet_T_above", ["8187 contains","8188","8229"], "-1.0", ["jet_T_above"], JE)
s("jet_dr", ["8200"], "0.0", ["jet_dr"], JE)
s("jet_profile_interval", ["8203"], "0.0", ["jet_profile_interval"], JE)
s("side_power_interval", ["8206"], "0.0", ["side_power_interval"], JE)
s("nozzle_descent", ["8237"], "prescribed", ["feet"])
EX = SPE + " & h_expr/T_flame_expr given"
s("foot_body_clearance", ["8248 contains (abort if !feet)","8308"], "0", ["foot_body_clearance"], EX+" & nozzle_descent = feet")
s("h_conv", ["8259 contains (abort with expr)","8342"], "0.0", ["h_conv"], SPE+" & no h_expr (constant flame)")
s("T_flame", ["8259 contains (abort with expr)","8343"], "0.0", ["T_flame"], SPE+" & no h_expr (constant flame)")
s("nozzle_z0", ["8276 contains (abort under feet)","8318"], "REQUIRED (expr & prescribed)", ["nozzle_z0"], EX+" & nozzle_descent = prescribed")
s("nozzle_feed", ["8276 contains (abort under feet)","8319"], "0.0", ["nozzle_feed"], EX+" & nozzle_descent = prescribed")
s("foot_r_inner", ["8279"], "REQUIRED under feet", ["foot_r_inner"], EX+" & feet")
s("foot_r_outer", ["8280"], "REQUIRED under feet", ["foot_r_outer"], EX+" & feet")
s("foot_standoff", ["8281"], "REQUIRED under feet", ["foot_standoff"], EX+" & feet")
s("foot_rule", ["8289"], "pads", ["foot_rule"], EX+" & feet")
s("foot_npads", ["8295"], "3", ["foot_npads"], EX+" & feet")
s("foot_pad_quantile", ["8298"], "0.9", ["foot_pad_quantile"], EX+" & feet")
s("nozzle_feed_max", ["8301"], "inf", ["nozzle_feed_max"], EX+" & feet")
s("nozzle_collision_radius", ["8329"], "radius (prescribed) / foot_r_inner (feet)", ["nozzle_collision_radius"], EX)
