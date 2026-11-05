# Physical origin of the thermal-spallation onset temperature in granites, and the high-temperature material data needed to model it

Research note, 2026-11-02. Scope: three pre-registered questions about the onset temperature at which
granites (Grimsel / Central Aare in particular) begin to spall under flame-jet, laser or hydrothermal
heating, and the temperature-dependent property data a thermoelastic/LEFM onset criterion needs.

Conventions used throughout:
- **[MEASURED]** = an experiment on the named rock in the cited paper. **[MODELLED]** = output of a model.
  **[ESTIMATE]** = assembled from other data (e.g. mineral mixing). **[fig]** = read off a published figure
  by the note-taker (±5 % of axis span). **[via X]** = the number was taken from a secondary source X that
  quotes the primary; the primary was not opened here.
- "Hot" = tested at temperature. "Cooled" = heated, cooled to room temperature, then tested.
- Local files: `franck_material/` and `tests/MMWSpalling/validation/hoeser/*.md` under
  `/Users/tzetze20/amr_tools/alamo/`. The project's own property compilation is
  `tests/MMWSpalling/validation/grimsel_properties/grimsel_granite_properties_reference.md` ("the project
  reference" below).
- Context from the brief: the model fires at 843.76 K = 570.6 °C; its α(θ) table (ESTIMATE, ±30 % above
  300 °C) triples between 300 and 500 °C and plateaus by 600 °C; its E/E₀ halves by 600 °C; both tables
  stop at 650 °C.

---

## Q1. Is the spallation onset of granites a material temperature pinned by the quartz α–β transition (573 °C)?

### Takeaway
No. Every direct onset measurement on a granite that resolves individual spalls (Rauenzahn & Tester 1989;
Wilkinson & Tester 1993; Kant & von Rohr 2016; Höser 2016) lies in a band of roughly **400–630 °C**, flat with
heat flux over 0.2–7.4 MW/m², and the one rock that is almost pure quartz (Sioux quartzite, 93 %) spalls
*lowest* of all, at **~355 °C**. The only dataset whose range crosses 573 °C (Kant et al. 2017, Central Aare,
553–694 °C) measures the failure-plane temperature *after* the spall leaves and is above its authors' own
model. The quartz-transition explanation was a Russian hypothesis of the 1950s; the Canadian Mines Branch
found no petrographic trace of it in 1964 and Rauenzahn & Tester formally rejected it in 1989. Quartz
content does govern *spallability*, but through expansion mismatch and stiffness, not the inversion.

### Cited Findings

**Compiled onset / spall temperatures (rock, quartz content, method, hot measurement unless noted):**

| Rock | Quartz | Onset or spall T | Flux range | Method | Source |
|---|---|---|---|---|---|
| Sioux quartzite | ~93 % (fine, 0.2–0.5 mm) | **331–364 °C**, mean ≈ 355 °C, no flux dependence; onset = steady-state | 0.7–2.3 MW/m² | spinning CO₂ laser, IR scanner [MEASURED] | Wilkinson & Tester 1993, Table 5 and p. 54 ([on disk](/Users/tzetze20/amr_tools/alamo/franck_material/SurfaceTemperature.pdf)) |
| Barre granite (VT) | 27.2 % (Chayes 1950, via Beentjes 2018) | onset ΔT_s **465 / 435 / 445 / 472 K** above ambient ⇒ T_s ≈ 450–500 °C; "quite small, if any" flux dependence; **all below 573 °C even with maximum uncertainty** | 0.68–1.19 MW/m² | propane–O₂ torch, 2000 fps video + calorimetry [MEASURED] | Rauenzahn & Tester 1989, IJRMMS 26:381–399, Fig. 9 and p. 397 ([on disk](/Users/tzetze20/amr_tools/alamo/franck_material/RockFailureMechanism.pdf)) |
| Westerly granite (RI) | ~30 % (Violay et al. 2017) | onset ΔT_s **473 / 446 / 423 / 455 K** (same series) | same | same | same |
| Barre + Westerly | as above | **global onset 400–600 °C, "relatively insensitive to the applied heat flux level"**; steady-state surface T rises 500→930 °C (Barre), 480→740 °C (Westerly) | 0.2–2.9 MW/m² | torch + laser, IR scanner [MEASURED] | Wilkinson & Tester 1993, Rock Mech. Rock Eng. 26:29–62, Fig. 15 vs Figs 11–14 ([notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/wilkinson_tester_1993_notes.md)) |
| Gotthard granite (Gurtnellen) | n/r | θ_SP ≈ **600 °C** before, θ_FP ≈ **550 °C** after each spall; individual spalls 550–700 °C; 620→570 °C over SOD 50→130 mm (no clear trend); minimum θ_SP ≈ 500 °C | 0.45–0.87 MW/m² | 10 kW multi-orifice burner, 1.5 mm pyrometer at 1 kHz [MEASURED] | Kant & Rudolf von Rohr 2016, IJRMMS 84:177–186, Figs 10–13 ([on disk](/Users/tzetze20/amr_tools/alamo/franck_material/MinimalBoundaryConditions.pdf)) |
| Bethel White granite (VT) | n/r | same ≈ 600 / 550 °C; 630→540 °C over SOD 50→150 mm | same | same | same |
| Grimsel granite/granodiorite (Chuenzentennlen quarry) | ~30–46 % (see Q2) | **430–600 °C** one IR frame before the first spall: ≈ 500 (23 kW), 545 (28), 460 (33), 555 °C (38 kW), ±25–45 °C; "any trend was not recorded" with firing rate | 4.3–7.4 MW/m² | 7.9 mm CH₄–air jet, Optris PI450 IR camera, ε = 0.82 [MEASURED, fig] | Höser 2016, Diss. ETH 23896, Ch. 4.2, Fig. 4.11 ([notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/hoeser_2016_thesis_notes.md)) |
| Central Aare granite (same quarry) | 46 % (Kant et al. 2017a) | **θ_FP 553–694 °C, mean ≈ 630 ± 40 °C at 0 MPa**; ≈ 575 ± 30 at 27 MPa; ≈ 560 ± 35 at 48 MPa (failure-plane T *after* the spall leaves). Their LEFM model with a = 20 µm predicts 475 (390–560) °C [MODELLED] | ~1,400 °C flame, Höser flange | pyrometer 1.5 mm / 1 kHz [MEASURED, fig] | Kant, Rossi, Madonna, Höser, Rudolf von Rohr 2017, JGR Solid Earth 122:1805–1815, Fig. 7, Table 2 ([on disk](/Users/tzetze20/amr_tools/alamo/franck_material/JGR%20Solid%20Earth%20-%202017%20-%20Kant%20-%20A%20theory%20on%20thermal%20spalling%20of%20rocks%20with%20a%20focus%20on%20thermal%20spallation%20drilling.pdf)) |
| Webatuck dolomite | 0 % (≈ 80 % MgCO₃) | surface 457–645 °C at 1.24–3.92 MW/m² but "unclear whether true spallation was occurring"; surface disintegrated, likely calcination | laser [MEASURED] | Wilkinson & Tester 1993, p. 51 and Table 5 |
| Barre granite, hydrothermal | 27.2 % | comminution at supercritical-water jet **535–560 °C**, 22.5–27 MPa, below the atmospheric onset minima; preferential removal of quartz grains; "since the impinging jet temperature was below 573 °C, the α-β quartz phase transition should not play a role" | ~0.5 MW/m² class | SCW jet, autoclave [MEASURED] | Beentjes, Bender, Hillson, Tester 2018, SGW-43, [SGP-TR-213](https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2018/Beentjes.pdf) |
| Gneiss (Meier/Brkic hydrothermal flames) | n/r | only a quoted "surface temperature (~500 °C) as stated by Rauenzahn"; no onset measurement | — | — | Meier & Brkic WGC 2015 ([on disk](/Users/tzetze20/amr_tools/alamo/franck_material/Meier_Brkic_WGC2015_HydrothermalFlames_Gneiss.pdf)) |
| Gabbro, Rorschach sandstone, Paradiso granite | n/r | qualitative screening only: granites spalled most, **sandstone melted** | 23–54 kW | flame [MEASURED, qualitative] | Höser 2016, p. 41 |

- Rauenzahn & Tester 1989 report that Gray's compressive-strength criterion ΔT_s = (1−ν)σ_c/(β_r E),
  applied to the 13 Mines-Branch samples, gives an "essentially constant" ΔT_s, and they tabulate ΔT_s by rock
  *category* (their Table 2, refs 14–17) while warning that broad categories hide wide property scatter —
  i.e. the onset is set by thermoelastic properties, not by a phase temperature — [Rauenzahn & Tester 1989, pp. 385–386 (on disk)](/Users/tzetze20/amr_tools/alamo/franck_material/RockFailureMechanism.pdf).
- **Who proposed the quartz-transition mechanism and who refuted it:** the α–β inversion was "once thought to be
  important in causing spallation" in the Soviet thermal-drilling literature — Rauenzahn & Tester cite
  Vasiliev A. P., *The thermal method of drilling in hard rocks*, in *Destruction of Coals and Rocks*
  (Terpigor'yev & Protodyakonov, eds), Ugletekhizdat, Moscow (1958) as ref. [28]. Soles & Geller (Mines
  Branch Canada, Tech. Bull. TB 53, 1964) "examined the effect of rock physical properties and composition on
  susceptibility to jet piercing and concluded ... that no single rock property can predict or determine
  spallability ... because **no petrographic evidence for the α–β transition in quartz at 573 °C or the
  decomposition of dolomite was ever found, those transformations ... are certainly of little or no
  consequence**." Rauenzahn & Tester add: "A common argument for the relatively excellent spallability of
  quartzitic species has been the effect of the rapid increase in thermal expansion coefficient of quartz
  near the α–β phase transition temperature of 573 °C [28]. As more evidence against such a claim, the results
  of this spallation temperature study add to previous data of spallation temperatures that are all below
  this equilibrium condition. In fact, none of the temperatures calculated from the propane experiments were
  above 573 °C, even with maximum levels of uncertainty included." —
  [Rauenzahn & Tester 1989, pp. 385 and 397 (on disk)](/Users/tzetze20/amr_tools/alamo/franck_material/RockFailureMechanism.pdf).
- Thirumalai (US Bureau of Mines, 1969): "the process of failure by thermal spalling in predominantly
  quartz-containing rock (Sioux Quartzite) takes place at temperatures considerably below the high-low
  inversion of quartz"; "spalling behavior in a rock appears to depend on the thermal expansion and
  shear-strain characteristics of the rock" — [Thirumalai, *Process of thermal spalling behavior in rocks — an exploratory study*, 11th US Symp. Rock Mech. (ARMA-69-0705), OSTI 5804838](https://www.osti.gov/biblio/5804838) (abstract via search summary; page itself did not load).
- Soles & Geller found "an increased propensity to spall with the volume fraction of quartz" —
  [Walsh et al., LLNL-TR-632239 (2013), OSTI 1088469](https://www.osti.gov/servlets/purl/1088469); Rauenzahn & Tester paraphrase: "Quartz and nepheline promote spallation, while mafic, micaceous and other soft, pliable minerals have an adverse effect on spallability [29 = Freeman, Sawdye & Mumpton 1963, Colo. Sch. Mines Q. 58:225–252]". Wilkinson & Tester: "spallation drilling works best in rocks with high quartz (SiO₂) content, such as granite, quartzite, and some sandstones" (p. 30, on disk).
- Preston (1926) and Preston & White (1934, *J. Am. Ceram. Soc.* 17:137–144) proposed the **buckling-of-a-compressed-skin** mechanism for brick spalling, which the LLNL group (Walsh, Lomov, Roberts) adopts as Figure 1 of each report; none of the four LLNL documents on disk invokes the quartz transition (grep for "573", "alpha-beta", "transition" returns nothing) — [Walsh et al. SGW 2012, LLNL-CONF-523558 (on disk)](/Users/tzetze20/amr_tools/alamo/franck_material/Walsh_Lomov_Roberts_SGW2012_GrainScaleFailure_LLNL-CONF-523558.pdf); [LLNL-TR-632239 (on disk)](/Users/tzetze20/amr_tools/alamo/franck_material/Walsh_2013_ModelingThermallyInducedFailure_LLNL-TR-632239.pdf).
- The LLNL grain-scale model (quartz/plagioclase/K-spar grains, Drucker–Prager + Johnson–Cook) "fail[s] to
  produce significant spalling at fluid temperatures less than ~300 °C above the ambient rock temperature",
  and damage grows continuously above that with no special temperature [MODELLED] — Walsh, Lomov & Roberts,
  SGW 2012, LLNL-CONF-523558, "Results and Discussion" (on disk, link above). Their Table 1 mineral
  properties: bulk modulus 37.0 / 50.8 / 53.7 GPa, Poisson 0.08 / 0.26 / 0.28, conductivity 7.7 / 1.5 / 1.5 W/m·K for quartz / plagioclase / K-spar.
- Rauenzahn & Tester's Weibull theory gives a flux exponent 3/(m+3) for the *steady* spalling temperature
  (0.13 at m = 20) and predicts ΔT_s 400–550 °C for 1–10 MW/m² with E 45 GPa, β 1e-5, σ₀ 70 MPa, m 20
  [MODELLED]; Wilkinson & Tester separate the flat *onset* (weakest flaw) from the rising *steady-state*
  surface mean (overheating of partly buckled flakes), and note Sioux quartzite shows no such rise because
  it has no soft phase — [notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/rauenzahn_tester_1989_notes.md), [notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/wilkinson_tester_1993_notes.md).
- Kant & von Rohr 2016 Fig. 14 regress their Gotthard/Bethel points together with Wilkinson's Barre/Westerly
  and Rauenzahn's data and state "the spallation temperature increases linearly with heat flux" (≈ 530 →
  950 °C over 0.35–2.4 MW/m²) — but the high-flux points are Wilkinson's *steady-state* means on other rocks
  by another method, so this regression mixes the two quantities above — [Kant & von Rohr 2016 (on disk)](/Users/tzetze20/amr_tools/alamo/franck_material/MinimalBoundaryConditions.pdf); contradicted on the onset by Wilkinson & Tester 1993 Fig. 15 and Höser 2016 Fig. 4.11.
- Kant et al. 2017 attribute their measured θ_FP sitting *above* the model to the first spall rarely detaching
  inside the 1.5 mm pyrometer spot; pressure lowers the measured onset by ≈ 70 °C over 48 MPa (model −45 °C)
  — [Kant 2017 notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/kant_2017_jgr_notes.md).
- Potter Drilling "identified the onset of hydrothermal spallation at a water jet temperature of 700 °C"
  (fluid temperature, not rock surface) — quoted in [Beentjes et al. 2018](https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2018/Beentjes.pdf).
- Mechanistic reason quartz matters without the inversion: volumetric CTE of quartz 4.98e-5 /°C vs feldspar
  ≈ 2e-5 /°C (Robertson 1988, via Beentjes 2018); Rauenzahn & Tester's inclusion model gives ~10 MPa tension
  around a quartz grain for a 40 °C rise near a free surface — [Rauenzahn & Tester 1989 Eq. 6–8, Fig. 7 (on disk)](/Users/tzetze20/amr_tools/alamo/franck_material/RockFailureMechanism.pdf).
- Acoustic emission on a quartz-bearing rock shows "a strong peak of microcracking at the phase transition
  temperature for quartz (∼573 °C) superimposed upon a background of microcracking due to thermal expansion",
  and a second peak at ~800 °C from hornblende/chlorite dehydroxylation; permeability/porosity jump in the
  500–600 °C range [MEASURED, slow furnace heating] — [Glover, Baud, Darot, Meredith, Boon, Le Ravalec 1995, Geophys. J. Int. 120:775–782](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-246X.1995.tb01852.x). This is *bulk* damage under uniform slow heating, not a surface onset.

### Inferences
- Across six granites and two methods, the first-spall temperature is a **material property band of ≈ 430–630 °C
  that does not scale with flux**; the band's lower half is below 573 °C and the quartz-richest rock (Sioux
  quartzite) spalls ~200 °C below it. A transition-pinned onset would put quartzite *at or above* 573 °C and
  would move with quartz content toward it; the data move the other way (more quartz ⇒ lower onset, because
  of higher stiffness, higher expansion mismatch, and no soft phase to relieve stress).
- The project's 570.6 °C firing temperature is inside Höser's Grimsel band (430–600 °C) and Kant 2017's
  Central Aare band (553–694 °C), but its coincidence with 573 °C is numerical, not mechanistic: the
  reconstructed onset is where σ = E·∫α dθ/(1−ν) reaches the flaw criterion, and the *estimated* α-tripling
  between 300 and 500 °C is what lands it there (see Q2).
- The high end of Kant 2017 (to 694 °C) is a failure-plane temperature biased upward by spot sampling; the
  pre-spall surface temperatures θ_SP in Kant 2016 (≈ 600 °C mean) and Höser (≈ 500–555 °C) are the better
  analogue of a model's firing temperature.
- Glover's 573 °C AE peak is real but describes bulk microcrack generation in a uniformly heated sample over
  minutes; at spallation heating rates (hundreds of K/s, 0.1–1 mm hot layer) the surface has already failed by
  buckling before the inversion is reached, which is exactly what Rauenzahn, Thirumalai and Beentjes observe.

### Gaps
- No onset measurement for basalt, gabbro or gneiss resolved to a surface temperature was found; Höser's
  screening is qualitative. The brief's "Kant 2016 Gotthard 570–620 °C, Bethel 540–630 °C" corresponds to the
  SOD-dependent θ_SP ranges in Kant & von Rohr 2016 Fig. 11 as read in the project notes.
- Soles & Geller 1964 (TB 53), Gray 1965 and Vasiliev 1958 were not opened; their content is known here only
  through Rauenzahn & Tester 1989 and LLNL-TR-632239.
- The OSTI landing page for Thirumalai 1969 did not load (socket error twice); the quoted sentence comes from
  a search-engine rendering of the OSTI abstract.
- Quartz contents of Gotthard and Bethel granites were not found.

---

## Q2. Does a measured thermal-expansion curve exist for Grimsel / Central Aare granite through 20–800 °C? What measured curves exist for comparable granites, and what is the modal quartz of Grimsel granite?

### Takeaway
No dilatometric α(θ) for Grimsel or Central Aare granite was found: the ETH thermal-property paper on Central
Aare granite (Kant et al. 2017, GRL) measured diffusivity, conductivity and heat capacity only, to 500 °C, and
the project's α table is explicitly an **estimate** from mineral mixing. Measured granite curves (Westerly,
Stripa, Climax quartz monzonite, Eibenstock, eight Egyptian granites, Czech granites, Nanan granite) agree on
the shape: ≈ 6–12 ×10⁻⁶/K near room temperature, rising to ≈ 10–15 ×10⁻⁶/K by 325 °C and ≈ 22–25 ×10⁻⁶/K by
475 °C, with a sharp dilatometer peak at ≈ 570 °C and a second step at ≈ 870 °C. Central Aare granite is
reported with 46 % quartz by the ETH drilling group; Grimsel granodiorite from the Grimsel Test Site with
≈ 30–33 % quartz.

### Cited Findings

**Grimsel / Central Aare — what exists**
- Kant, Ammann, Rossi, Madonna, Höser, Rudolf von Rohr 2017, *Geophys. Res. Lett.* 44:771–776, "Thermal
  properties of Central Aare granite for temperatures up to 500 °C: Irreversible changes due to thermal crack
  formation": thermal diffusivity, heat capacity and conductivity 25–500 °C over three heating/cooling cycles;
  after the first cycle diffusivity and conductivity dropped to ~75 % of the initial room-temperature value,
  heat capacity showed no irreversible decay; correlations are given. **No thermal expansion measurement is in
  the abstract.** [MEASURED, no α] — [DOI 10.1002/2016GL070990](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016gl070990) (abstract via OpenAlex; full text paywalled).
- The project reference states: "No measured α(θ) for Central Aare granite was obtained ... Assembled from the
  minerals; the rock's α is dominated by quartz" and tabulates α_rock = 8.4, 9.5, 16.4, 23.4, 26.9, 28.3,
  28.4 ×10⁻⁶/K at 20, 300, 400, 500, 550, 600, 650 °C with "±10 % below 300 °C, ±30 % above" [ESTIMATE] —
  [grimsel_granite_properties_reference.md §4, Table 3](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/grimsel_properties/grimsel_granite_properties_reference.md).
- Room-temperature α used by the ETH group for Grimsel granodiorite: 8.5e-6 /K (Höser 2016 Table 3.1; Höser et
  al. 2018 IJRMMS 106:311–318); Kant et al. 2017 Table 1: 7.8–8.2e-6 /K (literature values, Nagra reports) —
  [Höser notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/hoeser_2016_thesis_notes.md); [Kant 2017 notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/kant_2017_jgr_notes.md).

**Measured α(θ) for comparable granites (dilatometry or triaxial displacement)**
- **Westerly and Stripa granite, hot, 19–350 °C, 6–55 MPa confinement** [MEASURED]: "α typically increased from
  6 to 12 ×10⁻⁶ °C⁻¹ at 40 °C to 10–15 ×10⁻⁶ °C⁻¹ at 325 °C"; increasing P from 6 to 55 MPa lowered α by
  1–5 ×10⁻⁶; both rocks isotropic in α; microcracks control the response — Heard & Page 1982, *J. Geophys.
  Res.* 87(B11):9340–9348, [DOI 10.1029/JB087iB11p09340](https://doi.org/10.1029/JB087iB11p09340).
- **Climax quartz monzonite, hot, 19–500 °C, ≤ 55 MPa** [MEASURED]: α "increased from 8 to 11 ×10⁻⁶ °C⁻¹ at
  40 °C to 22 to 25 ×10⁻⁶ °C⁻¹ at 475 °C", most values above the crack-free aggregate; Sudbury gabbro 6–11
  ×10⁻⁶ with "very little net increase" 19–300 °C — Page & Heard 1981, LLNL report (UNT Digital Library record,
  via OpenAlex); companion JGR paper above.
- **Westerly granite and other igneous rocks, 25–550 °C** [MEASURED]: mean volume CTE 17.8 ×10⁻⁶/°C (Fairfax
  diabase) to 33.6 ×10⁻⁶/°C (Westerly granite) — Richter & Simmons 1974, *Int. J. Rock Mech. Min. Sci.*
  11:403–411, [DOI 10.1016/0148-9062(74)91111-5](https://doi.org/10.1016/0148-9062(74)91111-5) (numbers via a search-engine summary of the abstract; not independently opened). Cooper & Simmons 1977, *Earth Planet. Sci. Lett.* 36:404–412, "The effect of cracks on the thermal expansion of rocks", [DOI 10.1016/0012-821X(77)90065-6](https://doi.org/10.1016/0012-821X(77)90065-6) — no abstract retrievable.
- **Eight commercial granitic rocks (Egypt), dilatometer NETZSCH DIL 402 PC, 5 °C/min, 30–1000 °C** [MEASURED]:
  modal quartz 24.5–37.7 %; linear α 10.97–14.7 ×10⁻⁶/K up to 300 °C; "the first sharp increase ... was
  detected around 570 °C, which represented the quartz mineral (α–β transition phase)", confirmed by DSC —
  Alzahrani, Lasheen & Rashwan 2022, *Materials* 15(6):2041, [DOI 10.3390/ma15062041](https://pmc.ncbi.nlm.nih.gov/articles/PMC8950072/).
- **Czech granites, 25–1000 °C** [MEASURED]: mean linear α 35–43 ×10⁻⁶/K over the whole range — Plevová,
  Vaculíková, Kožušníková, Ritz, Simha Martynková 2015, *J. Therm. Anal. Calorim.* 123:1555–1561,
  [DOI 10.1007/s10973-015-4996-z](https://link.springer.com/article/10.1007/s10973-015-4996-z) (value via search summary; article paywalled).
- **Eibenstock granite (Saxony), to 1000 °C** [MEASURED]: "linear thermal expansion coefficient ... shows a
  steep increase around 573 °C as well as at 870 °C"; UCS, E and tensile strength measured before and after
  treatment — Wang, Konietzky, Frühwirt & Dai 2020, *Acta Geotechnica*, [DOI 10.1007/s11440-020-00926-8](https://link.springer.com/article/10.1007/s11440-020-00926-8) (abstract via OpenAlex).
- **Nanan granite (12 % quartz)** [MEASURED]: linear α +62 % at 400 °C, +95 % at 500 °C, +175 % at 600 °C
  relative to 200 °C — Zhu, Yang, Wang, Tian, Jiang & Dou 2022, *Acta Geod. Geophys.*, [DOI 10.1007/s40328-022-00375-7](https://doi.org/10.1007/s40328-022-00375-7) [via project reference §4; no abstract on OpenAlex].
- Westerly granite AE/velocity during heating to 450 °C [MEASURED]: P-wave velocity −50 % at 450 °C on the
  first cycle and −40 % after cooling, i.e. extensive microcracking already *below* the inversion — Griffiths,
  Lengliné, Heap, Baud & Schmittbuhl 2018, *JGR Solid Earth*, [DOI 10.1002/2017JB015191](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB015191).
- Quartz single-crystal data used by the project reference [MEASURED, crystal]: 13.55e-6 /K ⊥c and 7.43e-6 /K
  ∥c at 20 °C (Bechmann–Kosinski polynomials); ~1.5 % (a) and ~0.9 % (c) linear strain from 25 °C to the
  573.1 °C inversion (Rosenholtz & Smith 1941; Kihara 1990 *Eur. J. Mineral.* 2:63; Carpenter et al. 1998
  *Am. Mineral.* 83) ⇒ mean isotropic α ≈ 24e-6 over 25–570 °C and ≈ 35e-6 over 300–570 °C; feldspars
  5–8e-6 /K weakly T-dependent (Fei 1995) — [project reference §4](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/grimsel_properties/grimsel_granite_properties_reference.md) (marked "recalled" there; primaries not opened here).

**Modal quartz of Grimsel / Central Aare granite**
- Central Aare granite (ETH drilling group, Chüenzentennlen quarry at the Grimsel Pass): **46 % quartz, 31 %
  K-feldspar, 20 % plagioclase, 3 % biotite**, minor retrograde chlorite; K-feldspar to ~4 mm, quartz ~0.9 mm
  [MEASURED] — Kant et al. 2017, GRL (as quoted in the [project reference](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/grimsel_properties/grimsel_granite_properties_reference.md) and in a [search summary of the GRL paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016gl070990)).
- Grimsel granodiorite, Grimsel Test Site cores, XRD (ETH Zürich): **quartz ~32 %, plagioclase ~38 %, K-feldspar
  ~20 %, biotite ~8 %** [MEASURED] — Wenning, Madonna, de Haller & Burg 2018, *Solid Earth* 9:683–698,
  [DOI 10.5194/se-9-683-2018](https://doi.org/10.5194/se-9-683-2018), as tabulated in [Selvadurai et al. 2019, Solid Earth 10:2001](https://se.copernicus.org/articles/10/2001/2019/); other reports give ~30 % quartz / 30 % plagioclase / 25 % K-feldspar / 15 % biotite, or 37 / 33 / 17 / 6 % (search summary citing Grimsel granodiorite studies).
- Primary petrography of the site: Keusen, Ganguin, Schuler & Buletti 1989, *Grimsel Test Site: Geology*,
  Nagra Technical Report NTB 87-14E, Geotest Zollikofen — [OSTI record](https://www.osti.gov/etdeweb/biblio/6997735) (not opened). Krietsch et al. 2018, *Sci. Data* 5:180269 compile Keusen's rock-mechanical values: Central Aar granite E 53.3 ± 11.0 GPa, ν 0.37 ± 0.12, UCS 169 ± 37 MPa, ρ 2660 ± 24 kg/m³, tensile 9.06 ± 1.48 MPa; Grimsel granodiorite E 47.3 ± 15.4 GPa, ν 0.33 ± 0.15, UCS 117 ± 48 MPa, tensile 6–11 MPa [MEASURED, room T] — [DOI 10.1038/sdata.2018.269](https://doi.org/10.1038/sdata.2018.269) (values via project reference §5–6).
- Rock identity: Höser's "Grimsel granite/granodiorite" and Kant's "Central Aare granite" both come from the
  Chuenzentennlen quarry; the project's `Grimsel.H` is this rock — [Kant 2017 notes](/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/validation/hoeser/kant_2017_jgr_notes.md).

### Inferences
- The project's estimated α table matches the *shape* of the measured granite curves (flat to ~300 °C, then
  rising) but its magnitude in 400–550 °C (16–27 ×10⁻⁶/K) should be compared with the hot-measured
  Climax quartz monzonite (22–25 ×10⁻⁶/K at 475 °C) and the mixing estimate's own ±30 %: it is plausible,
  not confirmed. Its plateau at 600–650 °C is an artefact of stopping the quartz ramp at 570 °C; real
  dilatometer curves show a *spike* at ≈ 570 °C followed by a drop to a lower β-quartz value and a second step
  at ≈ 870 °C (Eibenstock, Egyptian granites).
- Part of the measured rock α above ~300 °C is crack opening (Heard & Page; Griffiths), which relieves rather
  than loads a confined surface skin; a stress model should use the lower, pressure-confined values (Heard &
  Page's 55 MPa curve) or cap α, as the project reference recommends.
- With 46 % (Central Aare) vs 32 % (Grimsel granodiorite) quartz from the same massif, a mineral-mixing α for
  this rock carries at least a ±20 % compositional uncertainty on top of the mixing-rule uncertainty.

### Gaps
- No dilatometric α(θ) for Grimsel/Central Aare granite was found in any source; a search for "Grimsel granite
  thermal expansion dilatometer Nagra" returned none. The Nagra heater-test report (OSTI 5456594) may contain
  an in-situ expansion coefficient but was not opened.
- Hartlieb, Toifl, Kuchar, Meisels & Antretter 2016, *Minerals Eng.* 91:34–41, [DOI 10.1016/j.mineng.2015.11.008](https://doi.org/10.1016/j.mineng.2015.11.008): abstract only ("phase transitions (e.g. α-β quartz) influence the texture and stability" of granite, sandstone, basalt, 25–1000 °C); the specific α values, specific heat and diffusivity curves and rock origins could not be retrieved (ScienceDirect, Leoben repository and Hindawi secondary all refused).
- Full α(θ) tables of Richter & Simmons 1974, Cooper & Simmons 1977, Plevová 2015, Wang 2020 and Zhu 2022 were
  not accessible; only the abstract-level numbers above.
- Which rock Glover et al. 1995 used is not stated in the abstract.

---

## Q3. Which measured datasets give granite E, tensile strength, K_Ic and compressive strength vs temperature to 1000 °C, hot vs cooled, and where is the softening / brittle–ductile transition?

### Takeaway
Hot-tested granite data are sparse above 500 °C and show gradual loss: E and K_Ic fall ≈ 20–25 % by 400 °C
(Mount Yan, Fujian, Lanhelin granites), UCS and E rise slightly to 200–300 °C then fall, and under confinement
Westerly granite stays brittle to 700 °C (Wong 1982) and to 900 °C (Violay et al. 2017), becoming ductile only
at 1000 °C. After-cooling data (the bulk of the literature) collapse between 400 and 600 °C — K_Ic to ~25 % and
tensile strength to ~27 % by 600 °C, E to ~30 % by 800 °C — and the often-quoted "brittle–ductile transition at
500–600 °C" comes from those cooled uniaxial tests. The project's E/E₀ = 0.5 at 600 °C and the K_Ic ratio 0.6
are after-cooling-like numbers; no hot-tested granite dataset above 500 °C with tabulated values was reachable.

### Cited Findings

**Hot (tested at temperature)**
- **Mount Yan granite (North China, 1000 m), uniaxial, 20–400 °C, held 2 h** [MEASURED]: E = 29.98 − 0.0145·T
  GPa (28.5, 26.9, 25.2, 23.6 GPa at 100, 200, 300, 400 °C ⇒ −19 % at 400 °C), ν = 0.1185 + 0.0004·T,
  UCS = 64.05 − 0.0284·T − 0.0001·T² MPa (62.0, 59.5, 55.6, 50.4 MPa) — Wang Yu, Liu Bao-lin, Zhu Hai-yan et
  al. 2014, *Sci. World J.* 2014:650683, [DOI 10.1155/2014/650683](https://pmc.ncbi.nlm.nih.gov/articles/PMC3980840/).
- **Granite, real-time 25–700 °C, uniaxial, three loading rates** [MEASURED]: UCS and E "increase first and then
  decrease with the increase of temperature"; UCS "clearly increases at 200 °C and decreases gradually when the
  temperature exceeds 300 °C"; at 600 °C the loading-rate effect "decreases significantly"; real-time results
  lack the "thermal hardening" seen in quenched granite — Zhang Y., Zhang F., Yang K., Cai Z. 2022, *Geofluids*
  2022:9426378, [DOI 10.1155/2022/9426378](https://onlinelibrary.wiley.com/doi/10.1155/2022/9426378) (abstract; the value table was not retrievable).
- **Fujian granite (Shuitou; albite 39.6 %, biotite 39.2 %, quartz 16.8 %, tremolite 4.5 %), NSCB K_Ic at real-time
  25–500 °C** [MEASURED]: K_Ic ≈ 2.46 → ≈ 1.23 MPa·m⁰·⁵ (25 → 500 °C) at 0.1 mm/min (≈ 50 % remaining), ≈ 2.37
  → ≈ 1.30 at 0.01 mm/min, ≈ 2.31 → ≈ 1.75 at 0.001 mm/min (≈ 76 %); monotonic decrease; intergranular cracks
  25–300 °C, transgranular (feldspar) above 400 °C; biotite ductility moderates the loss; less damage than
  water-quenched specimens above 200 °C — Yang K., Zhang F., Meng F., Hu D., Tan X. 2022, *Geothermal Energy*
  10:13, [DOI 10.1186/s40517-022-00225-3](https://link.springer.com/article/10.1186/s40517-022-00225-3). ⚠ The project reference records this same paper as 1.31 → 0.99 MPa·m⁰·⁵ (25 → 500 °C); the two readings conflict and the ratio (~0.5–0.76) is the safer quantity.
- **Lanhelin granite, K_Ic at 22–200 °C after 150 °C pre-treatment** [MEASURED]: 1.52 (22 °C), 1.65 (54), 1.58
  (86), 1.49 (118), 1.39 (150), 1.37 (175), 1.31 MPa·m⁰·⁵ (200 °C) — Chandler et al. 2017 [via project reference §6].
- **Westerly granite and Black gabbro, double-torsion K_Ic and subcritical crack growth, hot 20–400 °C**
  [MEASURED]: fracture toughness decreases with temperature "due to the gradual opening of grain-grain
  boundaries in response to thermal stresses" after an initial increase (crack blunting); Westerly room-T K_Ic
  1.79 ± 0.02 MPa·m⁰·⁵ — Meredith & Atkinson 1985, *Phys. Earth Planet. Inter.* 39:33–51,
  [DOI 10.1016/0031-9201(85)90113-X](https://doi.org/10.1016/0031-9201(85)90113-X) (trend via search summary and project reference; values not retrieved).
- **Westerly and Stripa granite, hot E and K, 19–350 °C, 6–55 MPa** [MEASURED]: "From 19° to 350 °C, E decreased
  by as much as a factor of 2 and K decreased by 2 to 3 times, depending on P"; from 6 to 55 MPa E increased 3–6×
  — i.e. the loss is crack-controlled and largely removed by confinement — Heard & Page 1982, *JGR*
  87:9340–9348, [DOI 10.1029/JB087iB11p09340](https://doi.org/10.1029/JB087iB11p09340). Climax quartz monzonite: E and K "decreased by up to a factor of approx. 80 % from 19 to 500 °C" isobarically — Page & Heard 1981 (OpenAlex record).
- **Granite tensile strength, real-time miniature Brazilian** [MEASURED]: "increasing below 300 °C before
  decreasing, with sharp declines at 400–500 °C and 600–700 °C"; damage from dehydration, phase transformation
  and differential expansion — Li X., Tian Y., Li Z., Heng S., Zhang X., Liu B. 2025, *Int. J. Min. Sci.
  Technol.*, [DOI 10.1016/j.ijmst.2025.07.003](https://doi.org/10.1016/j.ijmst.2025.07.003) (abstract only).
- **Westerly granite, triaxial, hot, to 700 °C and 400 MPa, 95 experiments** [MEASURED]: "For dry samples, the
  effects of temperature and strain rate on the failure stress are not significant relative to that of
  pressure"; "temperature was shown to be a stabilizing factor"; class I and II post-failure behaviour
  observed; shear fracture energy ~10⁴ J/m² up to 250 MPa and 700 °C — Wong 1982, *Mech. Mater.* 1:3–17,
  [DOI 10.1016/0167-6636(82)90020-5](https://www.sciencedirect.com/science/article/abs/pii/0167663682900205).
- **Westerly granite (30 % quartz, 30 % oligoclase, 30 % microcline, 10 % biotite), triaxial, hot 600–1000 °C,
  130 MPa confining / 100 MPa effective, 10⁻⁵ s⁻¹** [MEASURED]: peak differential stress **815 MPa at 600 °C →
  719 MPa at 900 °C**, "considerably lower" at 1000 °C; "brittle up to 900 °C (shear fracture formation) but
  ductile at 1000 °C"; dilatant only to 800 °C; thermal shock did not change strength; extrapolated to
  geological strain rates the BDT is 400 ± 100 °C — Violay, Heap, Acosta & Madonna 2017, *Sci. Rep.* 7:7705,
  [DOI 10.1038/s41598-017-08108-5](https://www.nature.com/articles/s41598-017-08108-5).

**Cooled (heated, cooled, then tested)**
- **Eibenstock granite** [MEASURED]: E ≈ 62 / 60 / 45 / 20 GPa at 25 / 400 / 600 / 800 °C (≈ 97 / 73 / 32 %);
  tensile strength "nearly vanishes after a temperature treatment of 1000 °C" — Wang, Konietzky, Frühwirt & Dai
  2020, *Acta Geotechnica*, [DOI 10.1007/s11440-020-00926-8](https://link.springer.com/article/10.1007/s11440-020-00926-8) (E values via project reference §5; tensile statement via search summary).
- **Pingyi granite (Shandong), K_Ic, natural cooling, 20–900 °C** [MEASURED]: K_Ic 1.21 / 1.23 / 0.81 / 0.32 /
  0.20 MPa·m⁰·⁵ at 25 / 200 / 400 / 600 / 800 °C (≈ 100 / 102 / 67 / 26 / 17 %), −89 % at 900 °C; "brittle-ductile
  transition temperature range is identified as 500 °C to 600 °C"; attributed to thermal crack density, not the
  quartz transition — Niu, Wu, Chen & Zhao 2025, *Sci. Rep.* 15:13679, [DOI 10.1038/s41598-025-98787-2](https://www.nature.com/articles/s41598-025-98787-2).
- **Granite, K_Ic after treatment (RSOS)** [MEASURED]: 1.74 → 0.59 (400 °C) → 0.34 MPa·m⁰·⁵ (600 °C) — *R. Soc.
  Open Sci.* 6:190144 (2019), [DOI 10.1098/rsos.190144](https://royalsocietypublishing.org/rsos/article/6/7/190144/95502/Fracture-failure-of-granite-after-varied-durations) [via project reference §6; page refused].
- **Westerly granite, CCNBD K_Ic after treatment to 850 °C** [MEASURED]: toughness decreases with treatment
  temperature while fracture roughness increases (negative correlation); P-wave velocity decrease tracks K_Ic —
  Nasseri, Schubnel & Young 2007, *IJRMMS* 44:601–616, [DOI 10.1016/j.ijrmms.2006.09.008](https://doi.org/10.1016/j.ijrmms.2006.09.008); Nasseri, Tatone, Grasselli & Young 2009, *Pure Appl. Geophys.* 166:801–822, [DOI 10.1007/s00024-009-0476-3](https://link.springer.com/article/10.1007/s00024-009-0476-3) (values not retrievable).
- **Granite tensile strength after treatment** [MEASURED, granite unnamed in summary]: 10.59 (25 °C), 10.11
  (200), 7.83 (400), 2.84 MPa (600 °C) (≈ 27 % at 600 °C) — search summary of [Sci. Rep. 2024 s41598-024-59470-0](https://www.nature.com/articles/s41598-024-59470-0) and related Brazilian-test papers; "Brazilian tensile strength and P-wave velocity exhibit a clear decline beyond a temperature threshold of 450–600 °C".
- **Granite UCS after treatment** [MEASURED, unnamed in summary]: 35.5 % retained at 800 °C, 6.8 % at 1000 °C;
  "600 °C is the brittle–ductile transition critical temperature of granite stress–strain curves" — search
  summaries of [Processes 13:1082 (2025)](https://doi.org/10.3390/pr13041082) and related cooled-specimen studies.
- **Thermally damaged granite, characteristic temperatures** [MEASURED]: physical properties change
  significantly at **400, 800 and 1000 °C**; albite and K-feldspar decompose toward anorthoclase at high T;
  quartz is the only mineral whose XRD intensity does not fall — Peng, Li, Peng, Gan & Wang 2024, *Sci. Rep.*,
  [DOI 10.1038/s41598-024-65752-4](https://www.nature.com/articles/s41598-024-65752-4).
- **Strathbogie granite (Australia)**: Shao, Ranjith, Wasantha & Chen 2015, *Geothermics* 54:96–108,
  [DOI 10.1016/j.geothermics.2014.11.005](https://doi.org/10.1016/j.geothermics.2014.11.005) (tested at temperature, uniaxial, to 800 °C per secondary summaries); Kumari, Ranjith, Perera, Shao, Chen & Lashin 2016, *Geothermics* 65:44–59, [DOI 10.1016/j.geothermics.2016.07.002](https://doi.org/10.1016/j.geothermics.2016.07.002) (in-situ stress + temperature); Kumari et al. 2017, *Eng. Geol.* 229:31–44, [DOI 10.1016/j.enggeo.2017.09.012](https://doi.org/10.1016/j.enggeo.2017.09.012) (cooling treatments) — **numerical values could not be retrieved** (all abstracts empty on OpenAlex; ScienceDirect and MDPI pages refused). Yang, Tian & Ranjith 2017, *Energies* 10:756, reproduces the Strathbogie UCS vs T in a PFC model: intergranular cracking rises "rapidly and linearly" above an inflection temperature [MODELLED] — [DOI 10.3390/en10060756](https://doi.org/10.3390/en10060756).
- Yang, Ranjith, Jing, Tian & Ju 2017, *Geothermics* 65:180–197, "An experimental investigation on thermal damage
  and failure mechanical behavior of granite after exposure to different high temperature treatments",
  [DOI 10.1016/j.geothermics.2016.09.008](https://doi.org/10.1016/j.geothermics.2016.09.008) — cooled specimens; values not retrieved.
- Heuzé 1983, *Int. J. Rock Mech. Min. Sci.* 20:3–10, "High-temperature mechanical, physical and thermal
  properties of granitic rocks — a review" (modulus, Poisson's ratio, tensile and compressive strength,
  viscosity, expansion, density, permeability, melting, specific heat, conductivity, diffusivity),
  [DOI 10.1016/0148-9062(83)91609-1](https://doi.org/10.1016/0148-9062(83)91609-1) — content not retrievable here; the project reference summarises it as "little change to ~300 °C, steepening decline through the quartz transition".
- Homand-Etienne & Houpert 1989, *IJRMMS* 26:125–134, "Thermally induced microcracking in granites:
  characterization and analysis", [DOI 10.1016/0148-9062(89)90001-6](https://doi.org/10.1016/0148-9062(89)90001-6) — crack-density vs treatment temperature; no abstract retrievable.

### Inferences
- **Shape of the hot curves:** E and K_Ic decline roughly linearly from room temperature at ≈ 5–6 % per 100 °C
  (Mount Yan E; Fujian K_Ic at the fast rate), with no knee at 573 °C in any hot dataset reached; UCS has a
  small maximum at 200–300 °C (Zhang 2022; Wang 2014 shows only decline). Extrapolating the hot linear trends
  gives ≈ 70–75 % of E and ≈ 50–65 % of K_Ic at 570 °C — higher than the project's 0.55 and 0.6 — but this is an
  extrapolation past the last hot datum (400–500 °C), not a measurement.
- **Shape of the cooled curves:** flat or slightly rising to 200–300 °C, knee at ≈ 400 °C, collapse to
  25–35 % by 600 °C (K_Ic, tensile) and to ≈ 30 % (E) / 35 % (UCS) by 800 °C, near-zero tensile strength and
  < 10 % UCS at 1000 °C. These curves measure *accumulated crack damage after cooling*, which the material
  ahead of a spallation front has not experienced; they overstate the loss relevant to onset.
- **Softening / brittle–ductile transition:** under confinement granite is brittle to 700 °C (Wong) and 900 °C
  (Violay) at laboratory strain rates; the "500–600 °C BDT" is an unconfined, cooled-specimen stress–strain
  shape change produced by crack damage. Thermal softening relevant to a spall skin (hot, confined laterally,
  loaded in seconds) is therefore gradual, with Violay's 12 % strength drop from 600 to 900 °C the best hot
  number above 500 °C.
- The project's E/E₀ = 0.5 and K_Ic/K_Ic,0 = 0.6 at 600 °C sit between the hot and cooled literature; the
  plateau above 600 °C has no measured basis either way.

### Gaps
- No hot-tested granite dataset with tabulated E, UCS or K_Ic above 500 °C was reachable: Zhang 2022 (to 700 °C)
  and Shao 2015 (to 800 °C) exist but their tables were behind refused pages; Pang et al. 2025/2026, *Deep
  Underground Sci. Eng.*, [DOI 10.1002/dug2.70008](https://doi.org/10.1002/dug2.70008), compares real-time hot vs subsequently cooled granite directly but the page refused.
- Hot tensile-strength values (Li et al. 2025) and Meredith & Atkinson 1985's K_Ic(T) values were not retrieved.
- Nasseri 2007/2009 Westerly K_Ic values at 250/450/650/850 °C: not retrieved.
- Heuzé 1983 content: citation only.
- A permission denial blocked bulk PDF downloads (arXiv microwave review quoting Hartlieb, Zhang 2022 PDF, MDPI
  Energies PDF, Acta Geotechnica PDF); those were not retried.
- No Grimsel/Central Aare property-vs-temperature data beyond Kant 2017 GRL's thermal properties to 500 °C were
  found; Dutler/Nejati 2018 (ETH report, [DOI 10.3929/ethz-b-000289969](https://doi.org/10.3929/ethz-b-000289969)) give room-temperature anisotropic K_Ic and tensile strength only.

---

## One-line verdicts

- **Q1:** The evidence decides **against** a quartz-transition-pinned onset: granites first spall at ≈ 430–630 °C
  independent of flux, 93 %-quartz quartzite spalls at ≈ 355 °C, Barre/Westerly onsets are all below 573 °C
  within uncertainty, and the hypothesis (Vasiliev 1958) was rejected on petrographic grounds by Soles & Geller
  1964 and Rauenzahn & Tester 1989; quartz raises spallability through expansion mismatch and stiffness, not the inversion.
- **Q2:** **No measured α(θ) exists for Grimsel/Central Aare granite** (Kant 2017 GRL measured only diffusivity,
  conductivity and c_p to 500 °C); the model's α table is an estimate, and measured analogues (Westerly/Stripa
  6–12 → 10–15 ×10⁻⁶/K to 325 °C hot; Climax 8–11 → 22–25 ×10⁻⁶/K at 475 °C hot; Egyptian granites 11–15 ×10⁻⁶/K
  to 300 °C with a spike at ≈ 570 °C) support its shape to ±30 % but not its plateau; quartz is 46 % (Central
  Aare, ETH) or ≈ 32 % (Grimsel granodiorite, GTS cores).
- **Q3:** Hot-tested granite loses E and K_Ic only ≈ 20–25 % by 400 °C and stays brittle under confinement to
  900 °C; the steep 400–600 °C collapse and the "600 °C brittle–ductile transition" belong to cooled specimens,
  so the model's E/E₀ = 0.5 and K_Ic/K_Ic,0 = 0.6 at 600 °C are unsupported by any hot measurement and no hot data
  above 500 °C with values could be reached.
