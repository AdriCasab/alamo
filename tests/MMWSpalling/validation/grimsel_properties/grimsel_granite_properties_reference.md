Grimsel / Central Aare granite — temperature-dependent properties and spallation temperature
=========================================================================================

Reference data for the ALAMO MMWSpalling Meier (flame-jet) and Kant validation cases.
Compiled 2026-10-06 by the reviewer from the sources in §10. Long-form discussion and
the derivations are in `Claude_markdowns/2026-10-06b.md`; this file carries the laws, the
numbers, and the provenance of each one. No code was changed; the model inputs quoted
below were read from `tests/MMWSpalling/validation/meier/sp_meier_pilot/input_drilling` and
the resolved `metadata` of the CTRL4_FIXED leg (studies/p1g2c_shiftfix).

Conventions: θ is temperature in °C, T in kelvin (T = θ + 273.15). "Measured" means the
number comes from an experiment on the named rock in the cited paper. "[estimate]" means
it was assembled here from mineral data or from other granites. "[recalled]" means a
number quoted from memory of a paper that could not be re-opened today; treat it as
provisional until the paper is on disk.

The rock
--------
Central Aare granite, sampled at Chüenzentennlen ~1 km north of the Grimsel Test Site
(Kant et al. 2017a): 46 % quartz, 31 % K-feldspar, 20 % plagioclase, 3 % biotite, minor
chlorite; grain sizes from micrometres to ~4 mm (K-feldspar up to 4 mm, plagioclase and
biotite ~1 mm, quartz ~0.9 mm). Meier's drilling tests and the Grimsel Test Site datasets
(Krietsch et al. 2018) concern the same massif: Central Aar granite and the adjacent
Grimsel granodiorite. Kant & von Rohr 2016's threshold tests used Gotthard granite
(Gurtnellen) and Bethel granite (Vermont).

Model inputs today (for comparison)
-----------------------------------
    rho = 2750 kg/m^3            constant
    cp  = 790 J/(kg K)           constant   (microstructure.phase0.Cp)
    k   = 1.5 W/(m K)            constant   (microstructure.phase0.kappa)  => alpha_diff 6.9e-7 m^2/s
    E   = 30 GPa (per brief §9)  constant
    Weibull: a0_ig = 8 um, m = 20
    firing temperature in the Meier legs: 757-767 K (threshold-pinned, flux-independent)


1. Specific heat capacity  c_p(θ)   — MEASURED on this rock
-----------------------------------------------------------
Law (Kant et al. 2017a, GRL, Table 2; laser flash ASTM E1461, 25-500 °C, 25 K steps,
three samples, ±5 % device accuracy, all data inside a 10 % band; identical on heating,
cooling and later cycles, i.e. no hysteresis):

    c_p(θ) = 1370 - 178 000 / (θ + 271)          J/(kg K)        valid 25-500 °C
    c_p(25 °C) measured = 768 ± 48 J/(kg K)

Enthalpy from a 20 °C datum (closed form, use this for the H <-> T inversion):

    H(θ) = 1370 (θ - 20) - 178 000 ln((θ + 271) / 291)      J/kg

Cross-checks (independent rocks or assemblages):
  - Whittington, Hofmeister & Nabelek 2009 (Nature 458:319; erratum: 10^-6 -> 10^6), crustal
    assemblage 30 % quartz, 60 % albite, 5 % phlogopite, 5 % annite, molar mass 221.78 g/mol:
        C_P = 199.50 + 0.0857 T - 5.0e6 / T^2      J/(mol K)   T <= 846 K
        C_P = 229.32 + 0.0323 T - 47.9e6 / T^2     J/(mol K)   T >= 846 K
    Divide by 0.22178 kg/mol. Agrees with the Kant law to < 2 % between 300 and 550 °C.
  - Zhang/Oglesby granite law already in validation/zhang: c_p = -2.73e-4 θ^2 + 0.894 θ + 788
    (agrees to 3 %).
  - Eibenstock granite (Acta Geotechnica 2020): 706 -> 912 J/(kg K) from 25 to 600 °C.
  - Vosteen & Schellschmidt 2003, Fig. 4: mean magmatic-rock c_p ~750 at 0 °C, ~1000 at 300 °C.

Recommendation: use the Kant law to 773 K (measured range) and extrapolate to 873 K (3 %
from Whittington at 823 K). Above the α-β quartz transition (846 K) β-quartz has a lower
c_p and Whittington's second branch drops to ~880 J/(kg K); the λ-peak at the transition
(a few kJ/kg) is in neither law. Irrelevant unless T_s is pushed above ~850 K.

Consequence: mean c_p 20 -> 487 °C = 1005 J/(kg K) (model 790); enthalpy of removed rock
at the present firing temperature 1.24 GJ/m^3 with rho = 2650 (model 1.01), 1.43 GJ/m^3 at
823 K. Energy-limited removal rate scales inversely.
**Table 1 — specific heat capacity c_p [J/(kg·K)] and enthalpy from 20 °C**

|T [K]|θ [°C]|Kant 2017 fit|Whittington 2009|Zhang (repo ref)|H(θ) Kant [kJ/kg]|mean c_p 20→θ|
|---|---|---|---|---|---|---|
|293|20|758|750|806|-0|nan|
|373|100|890|882|875|66|829|
|473|200|992|982|956|161|894|
|573|300|1058|1052|1032|263|941|
|673|400|1105|1110|1102|372|979|
|723|450|1123|1136|1135|427|994|
|773|500|1139|1161|1167|484|1009|
|823|550|1153|1184|1197|541|1022|
|873|600|1166|878|1226|599|1034|
|923|650|1177|915|1254|658|1044|


2. Thermal conductivity λ(θ) and diffusivity κ(θ)   — MEASURED on this rock
---------------------------------------------------------------------------
Laws (Kant et al. 2017a, Table 2):

    λ_1(θ)  = 2000 / (θ + 563)        W/(m K)    first heating (intact rock)
    λ_2+(θ) = 3400 / (θ + 1385)       W/(m K)    after one excursion to 500 °C (cracked)
    κ_1(θ)  = 487 / (θ + 264)         mm^2/s     first heating
    κ_2+(θ) = 526 / (θ + 400)         mm^2/s     after cracking
    measured at 25 °C: λ = 3.247 ± 0.066 W/(m K), κ = 1.656 ± 0.09 mm^2/s, rho = 2550 kg/m^3

The irreversible 25 % drop at room temperature after the first cycle is thermal
microcracking along grain boundaries (crack density 1.6 -> 5.8 mm/mm^2 after 500 °C,
mean crack length 130 um; Kant 2017a §5). No further change on later cycles below 500 °C
(thermal Kaiser effect). c_p shows no such change.

General law for any crystalline rock given its room-temperature value (Vosteen &
Schellschmidt 2003, Phys. Chem. Earth 28:499, Eq. 4-5; T in °C, 0-500 °C):

    λ(T) = λ(0) / (0.99 + T (a - b / λ(0)))      a = 0.0030 ± 0.0015,  b = 0.0042 ± 0.0006
    λ(0) = 0.53 λ(25) + 0.5 sqrt(1.13 λ(25)^2 - 0.42 λ(25))

Seeded with 3.25 W/(m K) it reproduces λ_1 within 5 % (Table 2). The Zhang law
κ = 1/(0.3154 + 3.8e-4 θ) has the same shape within 10 %. Eibenstock granite: 3.65 -> 1.8
W/(m K) from 25 to 600 °C (Acta Geotechnica 2020).

Which curve: λ_1 for everything below the flake layer (rock heated once); λ_2+ for the
top millimetre that has already cracked. Recommendation: λ_1 as the single law; carry
λ_2+ as a sensitivity. κ follows from λ/(rho c_p) and must not be set independently.

Model comparison: 1.5 W/(m K) is the intact value at ~1060 °C. Cold rock is 3.2-3.4 W/(m K)
and 1.66 mm^2/s (2.4x the model's diffusivity); at the firing temperature (~490 °C) the
measured 1.9 W/(m K) and 0.64 mm^2/s are close to the model. The steady ablation rate does
not depend on λ; the thermal-layer thickness κ/v, surface gradient, flake-thickness
length ℓ, start-up transient and wall heating do.
**Table 2 — thermal conductivity λ [W/(m·K)] and diffusivity κ [mm²/s]**

|T [K]|θ [°C]|λ Kant intact|λ Kant cracked|λ V&S from 3.25|λ Zhang|κ Kant intact|κ Kant cracked|
|---|---|---|---|---|---|---|---|
|293|20|3.43|2.42|3.26|3.10|1.72|1.25|
|373|100|3.02|2.29|2.87|2.83|1.34|1.05|
|473|200|2.62|2.15|2.50|2.56|1.05|0.88|
|573|300|2.32|2.02|2.21|2.33|0.86|0.75|
|673|400|2.08|1.90|1.98|2.14|0.73|0.66|
|723|450|1.97|1.85|1.88|2.06|0.68|0.62|
|773|500|1.88|1.80|1.80|1.98|0.64|0.58|
|823|550|1.80|1.76|1.72|1.91|0.60|0.55|
|873|600|1.72|1.71|1.64|1.84|0.56|0.53|
|923|650|1.65|1.67|1.58|1.78|0.53|0.50|


3. Density   — MEASURED
-----------------------
    Central Aar granite        2660 ± 24 kg/m^3    Krietsch et al. 2018 (Keusen et al. 1989 data)
    Grimsel granodiorite       2706 ± 14           Krietsch et al. 2018
    Kant 2017a thermal samples 2550                Kant et al. 2017a (five bulk measurements)
    Gotthard / Bethel granite  2594 / 2566         Kant & von Rohr 2016 (validation/kant)
Recommendation: rho = 2650 ± 50 kg/m^3. Thermal volume change to 500 °C ~1 % (Kant 2017a);
neglect, or derive from α below. The model's 2750 is at the top of the range.


4. Linear thermal expansion α(θ)   — [estimate], mineral mixing at 46 % quartz
-------------------------------------------------------------------------------
No measured α(θ) for Central Aare granite was obtained (Kant et al. 2017b use one; the
paper is behind the Wiley wall). Assembled from the minerals; the rock's α is dominated by
quartz, whose data are excellent.

Quartz (α-quartz):
  - 20 °C: 13.55e-6 /K along a, 7.43e-6 /K along c (recommended values; arXiv 1612.07049 §III
    after Bechmann & Kosinski and Sosman; IEEE-UFFC "Constants of alpha quartz").
  - Polynomials referenced to 0 °C, valid to ~300 °C (Bechmann & Kosinski):
        a_11 = 13.71 + 6.5e-3 θ - 1.9e-6 θ^2     (x1e-6 /°C)   perpendicular to c
        a_33 =  7.48 + 2.9e-3 θ - 1.5e-6 θ^2     (x1e-6 /°C)   parallel to c
  - Toward the α-β transition (573.1 °C, Rosenholtz & Smith 1941) the expansion diverges;
    linear strains of ~1.5 % (a) and ~0.9 % (c) between 25 °C and the transition
    [recalled: Kihara 1990 Eur. J. Mineral. 2:63; Carpenter et al. 1998 Am. Mineral. 83:2],
    i.e. mean isotropic α ≈ 24e-6 over 25-570 °C and ≈ 35e-6 over 300-570 °C; a further
    0.63 % volume jump at the transition (Carpenter et al. 1998, via secondary quotation).
  - Used here: polynomials to 300 °C, then a linear ramp of the isotropic α to 53e-6 at 570 °C
    (reproduces the integrated strain).
Feldspars: 5-8e-6 /K, weakly T-dependent (Fei 1995, AGU Reference Shelf 2); used 5.5 + 3e-3 θ (x1e-6).

Rock-level checks: granites 5.5-8e-6 at room temperature (Britannica 8 ± 3; Eibenstock
granite 6-8e-6); Nanan granite (12 % quartz) +62 % at 400 °C, +95 % at 500 °C, +175 % at
600 °C relative to 200 °C (Acta Geod. Geophys. 2022); several granites average 35-43e-6 over
25-1000 °C (J. Therm. Anal. Calorim. 2015). All attribute the rise above ~300 °C to quartz.

Result (Table 3): α_rock = 8.4e-6 at 20 °C (ETH group uses ~8e-6), 9.5 at 300, 16 at 400,
23 at 500, 27 at 550 °C. Uncertainty ±10 % below 300 °C, ±30 % above. Caution: part of the
high-T differential expansion is taken up by grain-boundary microcracks (Kant 2017a §5)
and does not load the flake; keep the bulk α at the low side near T_s.


5. Elastic modulus and Poisson's ratio
--------------------------------------
Room temperature (MEASURED; Krietsch et al. 2018, from Keusen et al. 1989):
    Central Aar granite    E = 53.3 ± 11.0 GPa   ν = 0.37 ± 0.12   UCS 169 ± 37 MPa
    Grimsel granodiorite   E = 47.3 ± 15.4 GPa   ν = 0.33 ± 0.15   UCS 117 ± 48 MPa
    dynamic (30 MPa confinement)  E = 65-81 GPa
    Kant & von Rohr 2016: Gotthard 54 GPa, Bethel 58 GPa (validation/kant)
The model's E = 30 GPa is 1.8x low; σ_thermal ∝ E, so any flaw-size fit absorbs this.

Temperature dependence [band]:
  - Measured at temperature: Mount Yan granite E = 29.98 - 0.0145 θ GPa to 400 °C (-19 %),
    ν = 0.1185 + 4e-4 θ (ScientificWorldJournal 2014, PMC3980840); Geofluids 2022 real-time
    series 25-700 °C: E and UCS rise slightly to 200-300 °C then fall, clearly above 400 °C.
  - Measured after cooling (overstates loss; cracks stay open): Eibenstock granite 62 / 60 /
    45 / 20 GPa at 25 / 400 / 600 / 800 °C; coarse granites lose ~60 % by 600 °C.
  - Heuzé 1983 (Int. J. Rock Mech. 20:3) review: little change to ~300 °C, steepening decline
    through the quartz transition.
Recommendation: E(θ)/E_0 = 1.0 to 200 °C, 0.85 at 400, 0.65 at 500, 0.5 at 600 °C (±0.1);
ν = 0.25-0.30 rising slowly. The product E·α (Table 3) is flat to 300 °C and 1.7-1.8x at
400-550 °C: for a given ΔT the thermal stress near the firing temperature is ~1.8x the
constant-property value, which lowers T_s unless a0 is refitted.


6. Fracture toughness and tensile strength
------------------------------------------
Room temperature (MEASURED, other granites):
    Westerly granite    1.79 ± 0.02 MPa m^1/2   Meredith & Atkinson 1985 (via Chandler 2017)
    Lanhelin granite    1.52                    Chandler et al. 2017
    Fujian granite      1.31                    Geothermal Energy 2022
    heat-treatment baselines 1.2-1.7            RSOS 2019; Sci. Rep. 2025
    Central Aare: 1-1.5 (Stanford 2018 summary, unverified); Grimsel granodiorite 1.2-2.3
    anisotropic [recalled: Dutler et al. 2018]
Recommendation: K_Ic,0 = 1.5 ± 0.3 MPa m^1/2.
Tensile strength (MEASURED): Central Aar 9.06 ± 1.48 MPa; Grimsel granodiorite 6-11 MPa
(Krietsch et al. 2018).

Temperature dependence measured AT temperature (the relevant kind):
  - Chandler et al. 2017, Lanhelin granite pre-treated to 150 °C (temperature effect only):
        22 °C 1.52 | 54 °C 1.65 | 86 °C 1.58 | 118 °C 1.49 | 150 °C 1.39 | 175 °C 1.37 | 200 °C 1.31
  - Geothermal Energy 2022, Fujian granite, 0.001 mm/min:
        25 °C 1.31 | 100 1.10 | 200 1.03 | 300 1.01 | 400 0.92 | 500 0.99
        fit: K_Ic = 1.051 + 0.489 exp(-θ / 86.4)
  - Meredith & Atkinson 1985, Westerly 20-400 °C: initial increase (crack blunting) then
    gradual decrease (crack coalescence); values not retrieved.
After cooling (damage curves, not applicable ahead of the front): 1.21 -> 0.81 (400 °C) ->
0.32 (600 °C) (RSOS 2019); 1.74 -> 0.59 (400 °C) -> 0.34 (600 °C) (Sci. Rep. 2025).
Recommendation: K_Ic(θ)/K_Ic,0 = 1.0 to 100 °C, 0.90 at 200, 0.82 at 300, 0.78 at 400,
0.75 at 500 °C, ~0.6 at 573 °C (±0.1).
**Table 3 — mechanical laws (estimates, §4–§6): α_rock [10⁻⁶/K] at 46 % quartz, E/E₀, K_Ic/K_Ic,0, and E·α relative to 20 °C**

|T [K]|θ [°C]|α_rock|E/E₀|K_Ic/K_Ic,0|E·α rel.|
|---|---|---|---|---|---|
|293|20|8.4|1.00|1.00|1.00|
|373|100|8.7|1.00|1.00|1.04|
|473|200|9.1|1.00|0.90|1.08|
|573|300|9.5|0.93|0.82|1.04|
|673|400|16.4|0.85|0.78|1.66|
|723|450|19.9|0.75|0.77|1.78|
|773|500|23.4|0.65|0.75|1.81|
|823|550|26.9|0.58|0.65|1.84|
|873|600|28.3|0.50|0.60|1.69|
|923|650|28.4|0.50|0.60|1.69|


7. Spallation surface temperature T_s(q)
----------------------------------------
Measurements:
  - Wilkinson & Tester 1993 (Rock Mech. Rock Eng. 26:29-62), pyrometry under flame jets
    (values read from the thesis Figs 1.8/1.9, brief §9):
        Barre granite:    480 / 560 / 680 / 830 / 930 °C at 0.3 / 0.6 / 1.2 / 2.0 / 2.9 MW/m^2
        Westerly granite: 450 / 600 / 650 / 740 °C     at 0.3 / 1.2 / 1.5 / 2.9 MW/m^2
    Conclusion of the paper: mechanically determined Weibull parameters do not describe
    spallation; the T_s(q) slope fits m ≈ 8 (Barre), 13 (Westerly).
  - Kant & von Rohr 2016 (validation/kant): Gotthard and Bethel granite, methane-air burner,
    sensor-plate flux 0.45-0.87 MW/m^2, h_fl 630-1010 W/(m^2 K), flame 1100-1450 °C: mean
    spalling surface ≈ 600 °C, failure plane ≈ 550 °C, minimum 500 ± 50 °C, ≈100 °C drop per
    detachment, minimum h ≈ 500 W/(m^2 K); Gotthard stops spalling between SOD 130 and 140 mm.
  - Kant et al. 2017b (JGR 122:1805): 550-600 °C for granitic rocks; Δx = 200 um (Höser 2016).
  - Rossi 2018 Stanford summary (unverified): onset at ~400-600 °C surface temperature.

Law: Weibull flaw statistics give ΔT_s ∝ q^(3/(m+3)) (Rauenzahn 1986 Eq. 4.91; brief §9).
Power-law fits to Wilkinson's points: exponent 0.306 for Barre (m = 6.8), 0.229 for
Westerly (m = 10.1). Kant's Gotthard points (ΔT ≈ 440-560 K at 0.45-0.87 MW/m^2) lie on or
slightly above the Westerly line. The repo's m = 20 gives exponent 0.13 (nearly flat).

Recommendation (Grimsel / Central Aare):
    T_s(q) = 293 + 530 (q / 0.6 MW m^-2)^0.25   K,   ±50 K
    => 800 K at 0.5 MW/m^2 (Meier), 823 K at 0.6, 875 K at 0.87, 895 K at 1.0, 1010 K at 2.0
Model: 757-767 K, flux-independent — 40-60 K low at Meier's flux, 130 K low at 1 MW/m^2.
The exponent matters more than the level: it is what makes the removal rate respond less
than linearly to flux. In the criterion this is the Weibull modulus (8-12 for the T_s(q)
slope, not the mechanical 20) together with the flaw size; fit both jointly on the Kant
legs (0.45-0.87 MW/m^2) and check at Meier's 0.5, with E, α, K_Ic from §4-§6 in place
first, otherwise a0 absorbs the property errors.
**Table 4 — spallation surface temperature T_s vs heat flux q**

|q [MW/m²]|Barre fit (Wilkinson)|Westerly fit (Wilkinson)|recommended Grimsel law|
|---|---|---|---|
|0.3|741 K|722 K|739 K|
|0.45|800 K|764 K|786 K|
|0.5|816 K|776 K|800 K|
|0.6|846 K|796 K|823 K|
|0.87|913 K|841 K|875 K|
|1.0|939 K|858 K|895 K|
|1.5|1025 K|913 K|960 K|
|2.0|1092 K|955 K|1009 K|
|2.9|1188 K|1014 K|1079 K|



8. Spall thickness (for completeness)
-------------------------------------
Höser 2016 (ETH thesis): mean spall ≈ 200 um, Central Aare granite (used by Kant 2017b).
Rossi et al. 2020 (J. Pet. Sci. Eng. 190:107005, field test): populations at ~250 um and ~1 mm.
Meier thesis: cuttings < 0.1 mm average, p90 < 0.3 mm.  Wilkinson 1993: 0.4 mm at 6e-4 m/s.
Rauenzahn 1986: Weibull median thickness 2.67 ℓ/m.  Model: floor flakes 0.2 mm (right);
wall flakes capped at 1.5 mm (not).


9. Suggested order of adoption
------------------------------
  1. c_p(θ) + rho = 2650 (no criterion coupling; re-run Zhang/Hu/Kant thermal regressions).
  2. λ(θ) = λ_1 (κ derived).
  3. E = 53 GPa, ν ≈ 0.3 static values.
  4. E(θ), α(θ), K_Ic(θ) together with a joint refit of (a0, m) on the 21 Kant legs.
  5. Check T_s(q) at Meier's flux against §7; then the wall flake cap against §8.


10. Sources and how each was consulted
--------------------------------------
Read in full (PDF text on disk during compilation):
  [K17a] Kant M.A., Ammann J., Rossi E., Madonna C., Höser D., Rudolf von Rohr Ph. (2017)
         Thermal properties of Central Aare granite for temperatures up to 500 °C: irreversible
         changes due to thermal crack formation. Geophys. Res. Lett. 44:771-776.
         doi:10.1002/2016GL070990. Open access: ETH research collection.
         -> §1, §2, §3 laws; composition; crack densities.
  [K17c] Kant M.A. et al. (2017) Thermal spallation drilling, an alternative drilling technology
         for deep heat mining — performance analysis... Stanford Geothermal Workshop SGW 2017.
         https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2017/Kant2.pdf
         -> drilling-velocity model; parameters Δx = 200 um, λ = 2, κ = 0.7 mm^2/s, Θ_SP = 500 °C.
  [VS03] Vosteen H.-D., Schellschmidt R. (2003) Influence of temperature on thermal conductivity,
         thermal capacity and thermal diffusivity for different types of rock. Phys. Chem. Earth
         28:499-509.  -> §2 general law (Eq. 4-5), §1 Fig. 4 check.
  [Ch17] Chandler M.R., Meredith P.G., Brantut N., Crawford B.R. (2017) Effect of temperature on
         the fracture toughness of anisotropic shale and other rocks. Geol. Soc. Lond. Spec. Publ.
         454. http://www.normalesup.org/~brantut/files/Chandler_et_al_GSSP_2017.pdf
         -> §6 Lanhelin K_Ic(θ) table; Westerly 1.79 ± 0.02 and the Meredith & Atkinson trend.
  [Gu19] Gupta A., Choudhary H., Dixit S.A. (2019/20) Mean velocity scaling in plane turbulent
         wall jets. arXiv 1905.00248 (J. Fluid Mech.). -> not rock; used in the wall-jet note.
  [aQ]   arXiv 1612.07049, Material properties of α-quartz relevant for X-ray optics, §III.
         -> §4 quartz α at 20 °C and TPRC conductivity.
  [KvR16] Kant M.A., Rudolf von Rohr Ph. (2016) Minimal required boundary conditions for the
         thermal spallation process of granitic rocks. Int. J. Rock Mech. Min. Sci. 84:177-186.
         -> via tests/MMWSpalling/validation/kant (Table 2 properties, Table 4 flame BCs, §5 thresholds).
  [Zh]   Zhang/Oglesby granite material model -> via tests/MMWSpalling/validation/zhang.

Read as article page / abstract (numbers as displayed by the publisher page):
  [Kr18] Krietsch H., Doetsch J., Dutler N., Jalali M., Gischig V., Loew S., Amann F. (2018)
         Comprehensive geological dataset describing a crystalline rock mass for hydraulic
         stimulation experiments. Sci. Data 5:180269. PMC6259022. -> §3, §5, §6 statics.
  [GE22] Effect of real-time high temperature and loading rate on mode I fracture toughness of
         granite. Geothermal Energy 10 (2022), doi:10.1186/s40517-022-00225-3. -> §6 real-time K_Ic.
  [Eib20] Laboratory testing and numerical simulation of properties and thermal-induced cracking
         of Eibenstock granite at elevated temperatures. Acta Geotech. 15 (2020),
         doi:10.1007/s11440-020-00926-8. -> §1, §2, §4, §5 cross-checks.
  [Nan22] Effects of high temperature on the linear thermal expansion coefficient of Nanan granite.
         Acta Geod. Geophys. 57 (2022), doi:10.1007/s40328-022-00375-7. -> §4 relative rises.
  [MY14] Thermophysical and mechanical properties of granite and its effects on borehole stability
         in high temperature and three-dimensional stress. ScientificWorldJournal 2014, PMC3980840.
         -> §5 real-time E(θ), ν(θ).
  [IEEE] IEEE-UFFC, Constants of alpha quartz (Bechmann & Kosinski polynomials). -> §4.
  [Ro18] Rossi E. et al. (2018) Stanford Geothermal Workshop (CTMD rock-bit interaction).
         https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2018/Rossi.pdf — the property
         summary attributed to it above could not be confirmed in the PDF text; treat as unverified.

Via search-engine summaries only (numbers quoted, paper not opened):
  [WHN09] Whittington A.G., Hofmeister A.M., Nabelek P.I. (2009) Nature 458:319-321 (+ erratum
         Nature 458, "10^-6 should be 10^6"). -> §1 C_P equations; cross-checked against [K17a].
  [RSOS19] Fracture failure of granite after varied durations of thermal treatment: an experimental
         study. R. Soc. Open Sci. 6:190144 (2019). -> §6 post-treatment K_Ic series.
  [SR25] Investigation of fracture properties of mode I fracture in heat-treated granite. Sci. Rep.
         (2025). -> §6 post-treatment K_Ic.
  [JTAC15] Thermal expansion behaviour of granites. J. Therm. Anal. Calorim. (2015),
         doi:10.1007/s10973-015-4996-z. -> §4 (35-43e-6 averaged 25-1000 °C).
  [Geo22] Effects of real-time high temperature and loading rate on deformation and strength
         behavior of granite. Geofluids 2022:9426378. -> §5 trend.
  [RS41] Rosenholtz J.L., Smith D.T. (1941) Linear thermal expansion and inversions of quartz,
         var. rock crystal. Am. Mineral. 26:103-109. -> §4 transition temperature 573.1 °C.
  [Ca98] Carpenter M.A. et al. (1998) Am. Mineral. 83:2-22 — 0.63 % volume jump at the transition.
  [He83] Heuzé F.E. (1983) High-temperature mechanical, physical and thermal properties of
         granitic rocks — a review. Int. J. Rock Mech. Min. Sci. 20:3-10. -> §5 trend only.

From the programme's own notes (thesis OCR and earlier reading, brief 2026-09-30a §9):
  [WT93] Wilkinson M.A., Tester J.W. (1993) Experimental measurement of surface temperatures
         during flame-jet induced thermal spallation. Rock Mech. Rock Eng. 26:29-62. -> §7 data.
  [Ra86] Rauenzahn R.M. (1986) MIT thesis; Rauenzahn & Tester (1989) Int. J. Rock Mech. Min. Sci.
         26:381. -> Weibull T_s(q) and thickness laws.
  [Hö16] Höser D. (2016) ETH thesis. -> 200 um mean spall.
  [Ro20] Rossi E., Jamali S., Saar M.O., Rudolf von Rohr Ph. (2020) Field test of a Combined
         Thermo-Mechanical Drilling technology, Mode I. J. Pet. Sci. Eng. 190:107005. -> spall sizes.
  [Me]   Meier thesis, Ch. 8 (validation/meier). -> cuttings sizes, rates, geometry.

Recalled, not re-opened (provisional):
  [Ki90] Kihara K. (1990) An X-ray study of the temperature dependence of the quartz structure.
         Eur. J. Mineral. 2:63-77. -> §4 a/c strains to the transition.
  [Du18] Dutler N. et al. (2018) — Grimsel granodiorite anisotropic K_Ic. -> §6.
  [K17b] Kant M.A., Rossi E., Madonna C., Höser D., Rudolf von Rohr Ph. (2017) A theory on thermal
         spalling of rocks with a focus on thermal spallation drilling. J. Geophys. Res. Solid Earth
         122:1805-1815, doi:10.1002/2016JB013800. Abstract only (550-600 °C; h > 500 W/(m^2 K)).

Not obtained: Balme et al. 2004 (JVGR 132:159); Nasseri et al. 2007 (IJRMMS 44:601);
Cooper & Simmons 1977; Fei 1995 tables; a Grimsel-specific α(θ) and K_Ic. These would
narrow §4 and §6 from ±30 % to ±10 %; none would change a recommendation.
