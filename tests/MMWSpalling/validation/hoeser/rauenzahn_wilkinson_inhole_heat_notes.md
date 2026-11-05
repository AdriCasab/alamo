<!-- Copied 2026-10-28 from a research agent's scratchpad report. Sources: Rauenzahn 1986 MIT PhD (DSpace 1721.1/14884) and Wilkinson 1989 MIT PhD (DSpace 1721.1/14031), image scans OCR'd with tesseract; page numbers = printed pages. Items marked [inference] are the agent's, not the theses'. -->

# In-hole heat transfer: Rauenzahn (1986) and Wilkinson (1989)

Sources: `litheat/rauenzahn.pdf` (R), `litheat/wilkinson.pdf` (W). PDF page number = printed page number in both
(checked: R p.379 prints "379", W p.225 prints "225"). Fresh OCR: `litheat/ocr2/r/p-NNN.txt` (R pp.17-460),
`litheat/ocr2/w/p-NNN.txt` (W pp.11-262). Figure renders that were read by eye are in `litheat/ocr2/fig/`
(`r379rot.png` = Fig 6.16, `rot386.png` = Fig 6.20, `rot373.png`/`rot377.png` = Figs 6.12/6.15, `r380/r381` = Figs 6.17/6.18,
`r429/r430` = Figs 7.14/7.15, `rot408` = Fig 7.5, `w228` = W Fig 6.1, `w239` = W Fig 6.6, `t404`/`t426` = Tables 7.1/7.3).

Read-off uncertainty: unless stated, a value read from a figure is good to about ±3 % of full scale on y and ±0.1-0.2
of the x tick spacing. "[inference]" marks my own derivation, not something the thesis says.

---------------------------------------------------------------------------------------------------------------------

## 0. Headline answers

* **The in-hole k-ε flux profile is peaked, then has a long, slowly falling tail.** In R Fig 6.16 (p.379) the flux peaks at
  about 0.25 R_noz, falls to about 1/2 of the peak by about 2.2 R_noz and to about 1/3 by about 3.5 R_noz. It then flattens:
  0.33 at 4-5.5 R_noz, about 0.18 at 7-10 R_noz, and about 0.12 at 14.6 R_noz, the end of the curve. Both meshes agree to
  about 10 % (p.378).
* **The cavity tail falls 2-4.5 times more slowly than the flat-plate tail from the same code and the same jet**
  (R Fig 6.20 against Fig 6.12, both with p_jet = 0.44 MPa and Z/R = 10) [inference: the ratio is mine, read off two figures].
* **The x axis of the in-hole plots is very probably distance along the surface, not radius** [inference]. Figs 6.16 and 6.20
  are captioned "spanwise distribution … along curved wellbore surface" and their axes run to 15 R_noz. The same code's hole
  radius at those conditions is only about 8-11 R_noz (Fig 6.17).
* **The only measured in-hole flux** is copper calorimetry in an ellipsoidal mock hole (R §7.3, Table 7.3). It was fitted with
  **Q(s) = Q_max exp(−s/s₀), where s is the distance along the curved face from the centreline** (Eq 7.5, p.416). This is
  exactly an "exposure rule". The fit gives s₀ = 3.0-5.2 cm and Q_max = 5.3-7.4 MW/m². The fits are poor and partly
  time-dependent.
* **Hole radius in both models comes from the self-consistency condition u_n = V_dr cos θ = Q/((ρC)_r ΔT_s), not from a
  flux threshold.** The radius is wherever the wall turns near-vertical (the mesh outlet is past θ = 85°). No minimum
  heat flux appears in either thesis: **NOT FOUND**.
* **R_h ≈ 0.45 × standoff** is not stated as a sentence in either thesis (**NOT FOUND verbatim**). It is the slope of the plotted
  simulation and experiment points over Z/R_noz = 20-70: W Fig 6.6 (p.239) and R Fig 7.14 (p.429).

---------------------------------------------------------------------------------------------------------------------

## 1. Rauenzahn's k-ε CFD of the jet in the drilled hole

### 1.1 Model, geometry and boundary conditions

* The model is axisymmetric, compressible and explicit (Lagrangian, then donor-cell advection, then diffusion), with standard
  k-ε (Launder-Spalding) and wall functions. The wall functions are extended for compressibility (Van Driest) and roughness
  (Dipprey-Sabersky): §6.3-6.6, pp.277-371.
* Domain: "the inlet to the flow field … is physically represented by the condition at the sonic throat of the drill nozzle"
  (p.275). "The outlet is forced to extend past the point on the hole where θᵢ = 85°" (p.369).
* Mesh: "The standoff distance … is usually divided into 10 to 40 cells … and the mesh is comprised of 50 to 80 cells
  axially" (p.369). There are "usually two cells across the fluid inlet from the sonic nozzle" (p.298). The cells are nominally
  orthogonal with "tents" on the curved wall (Fig 6.3b, Fig 6.5, pp.293-299).
* Rock-wall heat flux (Eq 6.116, p.354; restated p.362):
  Q_w = (T_p − T_w + ½u_p²/C_p) ρ C_p C_μ^{1/4} k_p^{1/2} / [ (Pr_t/κ) ln(F y_p/h_c) + 5.19 Pr^{0.44} Re_c^{0.2} ] × (1/A)[sin⁻¹((2A²−B)/(B²+4A²)^{1/2}) + sin⁻¹(B/(B²+4A²)^{1/2})].
  The roughness is h_c ≈ 0.5 mm, so "Re_c is about 50 to 100 for this application" (p.348), with E = 32.6/Re_c.
* The wall temperature is the Weibull spallation temperature, Eq 4.91 (p.354):
  T_w = T_ro + [(Q/(ρC_p)_r)³ ((1−ν)σ₀/(β_r E))^m (2(.693)/(πC_L²)) (m/α_r)³]^{1/(m+3)},
  "i.e. (T_w − T_ro) ∝ Q_w^{3/(m+3)}" (p.355). This is iterated cell by cell with Eq 6.116.
* Inlet turbulence is set to 20 % of the mean kinetic energy with l = D/3, and it "after 2 diameters … has all but vanished"
  (p.365).
* "Nearly all of the fluid forms a wall jet, and the boundary conditions on the drill housing are relatively unimportant"
  (p.300).
* Dimensionless inputs: p_jet/p₀, T_jet/T_ro and Z_dr/R_dr (p.364-367). Results are reported as R_h/R_dr and
  <St> = <h>/(ρC_p v)_jet (Figs 6.17/6.18).

### 1.2 Flux distribution along the cavity surface: Fig 6.16 (p.379)

The figure is titled "Comparison of spanwise distribution of heat transfer along curved wellbore surface for two different
nominal mesh cell sizes (1 cm and 0.5 cm)".

* y axis: Q_w/[(ρC_p v)_jet (T_jet − T₀)]. x axis: "distance from jet centerline / nozzle radius".
* Legend values: ρ_jet = 0.8 kg/m³, I_jet = 2.1×10⁶ J/kg, C_p = 1390 J/kg·K, v_jet = 840 m/s, T₀ = 300 K, γ = 1.264.
* Implied T_jet/T_ro = I/(C_v T₀) = 2.1e6/(1100·300) ≈ 6.4, about 1940 K [inference]. These values match Browning's
  full-scale Barre torch: Table 7.5 (p.441) gives T_jet = 1940 K, R_dr = 3.81 cm, p_jet = 0.46 MPa, and W Table 6.1 (p.241)
  gives ρ = 0.8 kg/m³, U = 837 m/s, R_noz = 3.81 cm [inference that this is the case plotted].
* **The standoff and the hole radius for Fig 6.16 are NOT stated.**

| s/R_noz (x axis) | 0.5 cm mesh: Q/(ρC_pv ΔT) | relative to peak | 1 cm mesh | relative to peak |
|---|---|---|---|---|
| 0 | 0.0073 | 0.92 | 0.0064 | 1.00 |
| 0.25 | **0.0079 (peak)** | 1.00 | 0.0064 | 1.00 |
| 0.5 | 0.0076 | 0.96 | 0.0064 | 1.00 |
| 0.9 | 0.0067 | 0.85 | 0.0059 | 0.92 |
| 1.2 | 0.0061 | 0.77 | 0.0054 | 0.84 |
| 1.5 | 0.0054 | 0.68 | 0.0049 | 0.77 |
| 2.0 | 0.0043 | 0.54 | 0.0041 | 0.64 |
| 2.4 | 0.0036 | 0.46 | 0.0040 | 0.63 |
| 3.0 | 0.0029 | 0.37 | 0.0034 | 0.53 |
| 3.4 | 0.0026 | 0.33 | 0.0030 | 0.47 |
| 4.0-4.6 | 0.0025 | 0.32 | 0.0027 | 0.42 |
| 5.2 | 0.0022 | 0.28 | 0.0024 | 0.38 |
| 5.8 | 0.0019 | 0.24 | 0.0022 | 0.34 |
| 6.3 | 0.0015 (step) | 0.19 | 0.0021 | 0.33 |
| 7-8 | 0.0014 | 0.18 | 0.0019-0.0020 | 0.30 |
| 9-10 | 0.0013 | 0.16 | 0.0016 | 0.25 |
| 12 | 0.0012 | 0.15 | 0.0012 | 0.19 |
| 14.5 | 0.00095 | 0.12 | 0.0010 | 0.16 |

* Width: the half-maximum is at about 2.2 R_noz on the 0.5 cm mesh and about 2.7 R_noz on the 1 cm mesh. One third of the
  peak is reached at about 3.4 R_noz on the 0.5 cm mesh, which confirms the coordinator's note of "~1/3 by 3 nozzle radii".
  It is reached at about 6 R_noz on the 1 cm mesh.
* Text on the mesh test: "the consequence of refining the mesh cell size from 1 cm to 0.5 cm (nominal) … is seen to be less
  than 10 percent" (p.378). This refers to integrated quantities; the local values differ by up to 30 % at r ≈ 7.
* At the end of the curve (14.6 R_noz, close to the outlet) the flux is about 12-16 % of the peak.

### 1.3 Same code, low against high pressure in the hole: Fig 6.20 (p.386)

* Caption: "spanwise distribution of heat transfer rates for the two cases shown qualitatively in Fig. 6.19 (high pressure
  and low pressure jets at same relative standoff distance)".
* Z_dr/R_dr = 10. The y axis is Q_w/Q_stag and the x axis is "distance from centerline (jet axis)/nozzle radius".
* **T_jet and R_noz are NOT stated.**

| x/R_noz | p_jet = 0.44 MPa (Q_stag = 8.9 MW/m²) | p_jet = 1.77 MPa (Q_stag = 15.5 MW/m²) |
|---|---|---|
| 0 | 1.00 | 1.00 |
| 0.6 | 1.00 | 1.03 |
| 1.0 | 0.90 | 1.06 |
| 1.7 | 0.75 | 1.15 |
| 2.0 | 0.62 | 1.18 |
| 2.4 | 0.62 | **1.60 (peak)** |
| 3.0 | 0.53 | 1.17 |
| 3.4 | 0.51 | 1.12 |
| 4.0 | 0.45 | **1.47** |
| 5.0 | 0.36 | 1.17 |
| 6.0 | 0.34 | 1.00 |
| 7.0 | 0.30 | 0.80 |
| 8.0 | 0.28 | 0.72 |
| 9.0 | 0.25 | 0.58 |
| 10 | 0.24 | 0.52 |
| 12 | 0.20 | 0.45 |
| 15 | 0.14 | 0.43 |

* Low pressure: "the decay in Q_w proceeds monotonically, as in most subsonic jets on flat plates" (p.385).
* High pressure (p_jet/p₀ > 9): "flow reversal suddenly becomes favored, and a bubble of slowly recirculating, rather hot
  gas is trapped in the hole bottom … for reasons similar to those responsible for slot jet-concave trough bubble flows
  (Dyban and Mazur (1970))" (p.385). The flux maximum then moves off-axis, to 2.4-4 R_noz.
* "If the standoff distance is increased to near 20 … the highest heat flux again returns to near the jet centerline"
  (p.385-387). Streamline sketches are in Fig 6.19 (p.382).

### 1.4 Flat plate, same code and same jet: Figs 6.12 and 6.15 (pp.373, 377)

Fig 6.12 conditions: ρ_jet = 0.71 kg/m³, p_jet = 4.4×10⁵ Pa, **R_jet = 0.02 m**, γ = 1.26, Z_dr = 20 cm (Z/R = 10).
The plate temperature is 294 K, which "enhanced heat transfer rates … the effective temperature driving force is augmented
by about 50 percent" (p.372). The y axis is St_loc = Q_w/[(ρC_pv)_jet(T_jet − T_w)] against r/R_noz.

| r/R_noz | Fig 6.12 St (Z/R = 10) | relative to stagnation | Fig 6.15 St (Z = 10 cm, Z/R = 5) | relative to stagnation |
|---|---|---|---|---|
| 0 | 0.0082 | 1.00 | 0.0094 | 1.00 |
| 0.75 | 0.0078 | 0.95 | 0.0098 | 1.04 |
| 1.25 | 0.0060 | 0.73 | 0.0075 | 0.80 |
| 1.75 | 0.0044 | 0.54 | 0.0053 | 0.56 |
| 2.25 | 0.0032 | 0.39 | 0.0037 | 0.39 |
| 2.75 | 0.0024 | 0.29 | 0.0027 | 0.29 |
| 3.2 | 0.0019 | 0.23 | 0.0020 | 0.21 |
| 4.2 | 0.0013 | 0.16 | 0.0014 | 0.15 |
| 5.2 | 0.0010 | 0.12 | 0.0011 | 0.12 |
| 6.2 | 0.0008 | 0.10 | 0.0009 | 0.10 |
| 7.2 | 0.0006 | 0.07 | 0.0007 | 0.07 |
| 8.7 | 0.00045 | 0.055 | — | — |

* Piesik's method is also plotted. It gives a constant core at 0.0035-0.0043 and a spike at about 1.05-1.3 R_noz reaching
  0.018 and 0.0135. The text says "Within Piesik's 'constant core' region … the heat transfer rate is less than half of the
  present model prediction … the asymptotic values at longer distances along the plate agree very well" (p.372).
* **Width scaling in the flat-plate runs** [inference]: the relative profile is almost identical at Z/R = 5 and Z/R = 10.
  The half-maximum is at about 1.8 R_noz in both. So for Z/R ≤ 10 the width scales with the **nozzle radius**, not the
  standoff. The thesis's own explanation is the supersonic core: "the equivalent 'potential' core length for most supersonic
  jets often extends to at least 20 nozzle diameters" (p.427).

### 1.5 Cavity against flat plate at the same jet [inference: ratio of Fig 6.20 low-p to Fig 6.12]

Both runs use p_jet = 0.44 MPa and Z/R = 10. Values are normalised to stagnation.

| x/R_noz | cavity (Fig 6.20) | flat plate (Fig 6.12) | cavity/flat |
|---|---|---|---|
| 2 | 0.62 | 0.47 | 1.3 |
| 3 | 0.53 | 0.26 | 2.0 |
| 4 | 0.45 | 0.17 | 2.6 |
| 6 | 0.34 | 0.10 | 3.4 |
| 8.7 | 0.27 | 0.055 | 4.9 |

* The cavity profile is slightly wider near stagnation, with half-maximum at about 2.7 against 1.8 R_noz.
* The cavity tail is much fatter: an e-fold of about 10 R_noz between 6 and 15 R_noz, against roughly a power-law decay
  r^-1.5 to r^-2 on the plate.
* Caveats:
  - The two normalisations differ in the plate temperature (294 K against the spall temperature).
  - The plate is smooth and the cavity rough.
  - T_jet for Fig 6.20 is not given.
  - If the Fig 6.20 x axis is surface distance (see 1.6), comparing it with the plate's radius is exactly the "exposure"
    comparison: at the same distance from stagnation, the cavity keeps 2-5 times more relative flux.

### 1.6 Is the in-hole abscissa radius or surface distance? [inference]

Captions say "spanwise distribution … along curved wellbore surface" (Fig 6.16) and "spanwise distribution" (Fig 6.20).
The axis labels say "distance from jet centerline/nozzle radius". For Z/R = 10 and p_jet/p₀ = 4.4, the same code's hole
radius is R_h/R ≈ 7-8 (Fig 6.17), yet the Fig 6.20 low-p curve continues smoothly to 15. Radius cannot exceed R_h, so the
abscissa must be distance along the surface (the curve continuing up the near-vertical wall toward the outlet), or at
least include the wall. **The thesis does not say this explicitly.**

### 1.7 Peak, mouth value and fall-off

| Case | Peak location | Value at the end of the curve (near the outlet, θ ≥ 85°) |
|---|---|---|
| Low-pressure jets (Fig 6.16, Fig 6.20 at 0.44 MPa) | on axis or within 0.25 R_noz | 12-16 % of peak at about 14.6 R_noz (Fig 6.16); 14 % at 15 R_noz (Fig 6.20) |
| High-pressure jet (Fig 6.20 at 1.77 MPa) | off-axis, 2.4 and 4 R_noz | 43 % of stagnation at 15 R_noz |

* "<Q_w>, the heat flux that computationally controls the drilling velocity, is not particularly sensitive to the exact
  value of hole radius since the majority of heat transfer occurs near the hole bottom" (p.367).
* "Computationally, over 90 percent of the actively spalling surface is encompassed within the farthest axial cells of the
  grid" (p.42).
* **No explicit hole-mouth (collar) value is given anywhere: NOT FOUND.** The domain ends just above θ = 85°.

### 1.8 Wilkinson's CFD distribution: W Fig 6.1 (p.228)

* W used a second-order finite-volume scheme with RK4 and the same k-ε plus wall functions (W §5.2.5, pp.167-174;
  Table 5.1 p.150).
* Fig 6.1 is the "Effect of grid resolution on Stanton number distribution": SOD = 60 R_noz, "a prescribed ellipsoidal
  boundary shape" (p.227). The y axis is St×10³ and the x axis is "Radius", presumably r/R_noz. St is
  Q/((ρC_pU)_jet (T_jet − T_s)) (W Eq 2.9, p.68).

| r/R_noz | 0 | 3 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 25 | 26.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| St×10³ | 0.78 | 0.75 | 0.70 | 0.63 | 0.56 | 0.47 | 0.38 | 0.31 | 0.25 | 0.19 | 0.14 | 0.08 | 0.045 | 0.02 |
| relative to stagnation | 1 | 0.96 | 0.90 | 0.81 | 0.72 | 0.60 | 0.49 | 0.40 | 0.32 | 0.24 | 0.18 | 0.10 | 0.06 | 0.03 |

* The two grids (2701 and 5225 nodes) differ by less than 2 % (p.227).
* At SOD = 60 the profile is broad: half-maximum at about 13.8 R_noz, about 0.23 SOD [inference]. Here the width clearly
  scales with **standoff** through jet spreading, not with nozzle size. Compare about 1.8-2.7 R_noz at Z/R = 5-10 in R.
* The fall toward 0 near r ≈ 27 is where the prescribed ellipse turns vertical. Plotted against radius, the whole wall is
  compressed onto r → R_h. **No surface-distance plot exists in W: NOT FOUND.**
* W gives no flux profile for a converged (self-consistent) hole: NOT FOUND. Only <St> against R_h is shown
  (Figs 6.2-6.5, 6.7).
* W does report a recirculation artefact from a modified k-ε model: a recirculation zone "starting at the stagnation point
  and having a radius of approximately 6 nozzle radii … predicted penetration rates were two orders of magnitude less than
  those routinely observed" (p.214). The standard k-ε model was kept.

### 1.9 Does the profile width scale with nozzle size or with standoff? (coordinator's priority)

* At Z/R ≤ 10, R's flat-plate runs give an identical relative profile at Z/R = 5 and 10. The width is about 1.8 R_noz to
  half-maximum, so it is set by the nozzle radius and the supersonic core.
* In R's cavity runs at Z/R = 10 the half-width is about 2.2-2.7 R_noz.
* At Z/R = 60 (W Fig 6.1) the half-width is about 14 R_noz, about 0.23 Z, so it is set by jet spreading with standoff.
* Neither thesis states a width law: **NOT FOUND as a statement.** R notes that the potential core is at least 20 D
  (p.427), and that for Z/D above about 20 "the heat transfer coefficient should decline rapidly since subsonic jet behavior
  should then dominate" (p.427).
* [inference] The crossover between nozzle-scaled and standoff-scaled width lies somewhere between Z/R ≈ 10 and 60. The
  measured holes (Z/R = 22-70) are in or past that transition.

---------------------------------------------------------------------------------------------------------------------

## 2. Measured heat flux in a mock hole (R §7.3, pp.406-427)

### 2.1 Apparatus

* "A portion of the copper block … was machined to replicate the hole geometry. An ellipsoidal shape with major and minor
  axes of 15 cm and 7.5 cm, respectively (Fig. 7.5, Table 7.2), was formed" (p.407).
* From Fig 7.5 (p.408): mouth radius 1.5 in (3.8 cm), block 7.5 in × 5 in (19.1 × 12.7 cm), cavity bottom 2 in (5 cm)
  above the base. The cavity is therefore about 7.7 cm deep with a 3.8 cm mouth radius, i.e. semi-axes 7.5 cm (depth) by
  3.75 cm (radius).
* "the geometry of this mock 'hole' was selected prior to any field drilling experiments, and an educated guess concerning
  hole shape was necessarily abstracted from computer simulations" (p.407).
* There were 10 type-K thermocouples (Table 7.2, p.409). TC1 is on the centreline 4.45 cm below the bottom face. Others are
  at 1.27-10.16 cm off-axis and 2.5-11.4 cm deep.
* The torch was the small propane-O₂ unit: 0.64 cm throat, so R_noz = 0.32 cm (p.392). O₂ flow was 0.3 m³/min STP and
  propane 0.06 m³/min (Table 7.3).
* The heat flux was inverted with the FE conduction code AYER (Fig 7.9, p.415).

### 2.2 Assumed spatial form, which is the key point for an exposure rule

* "Q(s) = Q_max exp(−s/s₀) (7.5), where, in this case, s denotes the distance along the curved face as tracked from the
  centerline locus, not the radial distance from the axis of symmetry (Fig. 7.9). This assumption is consistent with the
  compiled experimental evidence of Hrycak (1981), which supports the remarkable similarity among boundary layer wall-jet
  flows resulting from impingement on flat, convex, and concave objects" (p.416).
* The fits were poor: "the postulated form of the heat flux distribution on the ellipsoidal surface consistently
  overestimates the rate of temperature increase at late times. All thermocouple positions … exhibited this anomaly"
  (p.419).
* Alternative forms were tried and "did not prove any more useful" (p.419): a flat top for s < s* with an exponential
  beyond (7.6), and a Gaussian beyond s* (7.7).
* The final fit added a linear decay in time: Q(s) = Q_max(1 − d_t t/t_tot) exp(−s/s₀) (Eq 7.8, p.419). Figs 7.12-7.13
  (pp.421-422) show the fits.

### 2.3 Results: Table 7.3 (p.426), with derived columns added

| Standoff Z (cm) | Z/R_noz | Q_max (MW/m²) | s₀ (cm) | d_t | Q_max range (MW/m²) | s₀/R_noz [inf.] | s₀/Z [inf.] | Q(mouth)/Q_max [inf.] |
|---|---|---|---|---|---|---|---|---|
| 7.6 | 24 | 5.6 | 5.1 | .05 | 5.3-5.6 | 16 | 0.67 | 0.17 |
| 7.6 | 24 | 7.1 | 4.0 | .32 | 4.8-7.1 | 12.6 | 0.53 | 0.10 |
| 10.2 | 32 | 5.3 | 5.2 | .21 | 4.2-5.3 | 16.4 | 0.51 | 0.17 |
| 10.2 | 32 | 7.4 | 3.0 | .38 | 4.6-7.4 | 9.4 | 0.29 | 0.05 |
| 10.2 | 32 | 5.4 | 4.8 | .15 | 4.6-5.4 | 15.1 | 0.47 | 0.15 |
| 12.7 | 40 | 5.3 | 3.0 | .48 | 2.8-5.3 | 9.4 | 0.24 | 0.05 |
| 15.2 | 48 | 5.4 | 3.7 | .43 | 3.1-5.4 | 11.7 | 0.24 | 0.09 |

* Mouth column: s_mouth ≈ 9.1 cm, the quarter-perimeter of a 7.5 × 3.75 cm ellipse (Ramanujan formula) [inference]. The
  mouth fraction is exp(−9.1/s₀).
* "The final heat transfer data … for maximum heating rate (Q_max) are relatively invariant with standoff distance"
  (p.425).
* "the only firm conclusion to be drawn is that the heat transfer rate is somewhere between 3 to 6 MW/m²" (p.427).
* The time decay of 20-40 % was attributed most plausibly to falling propane flow on cool days (pp.424-425).
* [inference] s₀ shows no clear trend with standoff, given the scatter. The most decay-free run (d_t = .05) has s₀ = 5.1 cm
  at Z = 7.6 cm.
* Versus copper, the CFD "computed heat transfer rates usually matched calibrated values to within 20 percent" (Abstract,
  p.2), and "the Stanton numbers match extremely well, except for the experimental independence of <St>_eff,2 on standoff
  distance" (p.431, Fig 7.15).
* The rock-derived <St>_eff,1 (heat balance) is 3-4× lower (p.431). W attributes about half of that gap to the spall
  temperature (W conclusion 4, p.253).

### 2.4 Flat-plate copper calibration (R Ch 5)

Not re-read in detail here. R used Anderson and Stresino's radial exponential on a flat copper block (p.412, p.416) for
the subsonic propane welding torch. In R Table 6.3 (p.314) it is listed as Q_loc = Q_stag exp[−9.2 (r/D_F) Re_jet^{−.34}],
with D_F the flame diameter at impingement.

---------------------------------------------------------------------------------------------------------------------

## 3. Hole shape and R_h ≈ 0.45 · standoff

### 3.1 Experimental data: R Table 7.1 (p.404)

Small propane-O₂ torch, R_dr = 0.32 cm, p_jet = 0.45 MPa, T_jet = 2830 K (p.405).

| Granite | V_dr (m/h) | Z (cm) | R_h (cm) | Z/R_dr | R_h/R_dr | R_h/Z [inf.] |
|---|---|---|---|---|---|---|
| Barre | 2.74 | 21.6 | 8.9 | 70 | 28 | 0.41 |
| Barre | 3.25 | 12.7 | 7.6 | 42 | 24 | 0.60 |
| Barre | 4.11 | 10.2 | 6.4-7.0 | 34 | 20-22 | 0.66 |
| Barre | 4.57 | 10.2 | 5.7-6.4 | 34 | 18-20 | 0.59 |
| Barre | 5.49 | 7.6 | 3.8-4.4 | 26 | 12-14 | 0.54 |
| Barre | 6.31 | 7.6 | 3.8 | 26 | 12 | 0.50 |
| Barre | 6.77 | 6.4 | 3.2-3.5 | 22 | 10-11 | 0.52 |
| Barre | 7.31 | 7.0 | 3.5 | 24 | 11-12 | 0.50 |
| Westerly | 3.39 | 17.1 | 7.6-8.3 | 56 | 24-26 | 0.46 |
| Westerly | 4.57 | 15.2 | 6.4-7.0 | 50 | 20-22 | 0.44 |
| Westerly | 5.62 | 10.2 | 5.7-6.4 | 34 | 18-20 | 0.59 |
| Westerly | 6.31 | 8.9 | 4.4-5.1 | 30 | 14-16 | 0.53 |

* Measured R_h/Z = 0.41-0.66, mean about 0.53 [inference].
* Steady state was judged by an invariant radius over the last third to half of about 45 cm of drilling (p.399).
* The standoff was measured with a tape after shutdown (p.399).

### 3.2 Simulated R_h against standoff

**R Fig 7.14 (p.429)**: p_jet/p₀ ≈ 4.4, T_jet/T_ro ≈ 9.3. The simulation triangles are (Z/R, R_h/R) ≈ (20, 10.3),
(30, 14.5), (40, 19.4), (51, 23.1). The drawn line runs from about (14, 8) to (59, 27.5), a slope of about 0.43. R_h/Z is
0.52, 0.48, 0.49, 0.45 [inference].
* Text: "Agreement in trends and actual numerical values between numerical and experimental hole radii (Fig. 7.14) is
  extremely close. The predicted diameter, especially as standoff increases, tends to slightly underestimate the actual
  drilled hole size, but the scatter in the experimental data encompasses code-generated results" (p.428).

**W Fig 6.6 (p.239)**: the model with mass injection and C_p(T). The triangles are (20, 9.8), (30, 12.4), (40, 17.0),
(50, 21.9), (60, 26.6), (70, 32.4). The line runs from about (20, 8.4) to (70, 31.6), a slope of about 0.46 with an
intercept of about −0.8 R_noz. R_h/Z is 0.49, 0.41, 0.43, 0.44, 0.44, 0.46 [inference]. Rauenzahn's points (circles) are
about 1-2 R_noz higher.
* Text: "good agreement is obtained between predicted and experimental non-dimensional hole radii as a function of
  non-dimensional stand-off distances from 20 to 70" (p.238).
* "R_h ≈ 0.45·SOD" as a phrase: **NOT FOUND in either thesis.** It is a fair summary of the slope of W Fig 6.6 and R Fig 7.14
  over SOD 20-70 at p_jet/p₀ = 4.4.

**R Fig 6.17 (p.380)**: the parametric map. Values at T_jet/T_ro ≈ 6.4:

| p_jet/p₀ | Z/R = 10 → R_h/R | Z/R = 20 → R_h/R |
|---|---|---|
| 4.4 | 8.0 | 11.1 |
| 8.8 | about 12.6 | about 13 |
| 13.3 | 19.9 | 14.3 |
| 17.7 | about 20.5* | 17.3 |

* *Marked "numerical difficulties" (p.385).
* R_h falls slightly with T_jet. At Z/R = 10, R_h/Z = 0.8 for p 4.4 and about 2 for p 13.3.
* "large holes will be favored at large nozzle-to-rock separations for low pressure jets, a trend that reverses as higher
  pressure operation is adopted" (p.384).
* So the 0.45 slope is specific to low-pressure (4.4 atm) jets at Z/R ≥ 20.

### 3.3 What sets the hole radius in the model

**Rauenzahn, Eq 6.134 (p.368, checked on the page image):**
"forcing each boundary cell to conform to a self-consistent hole shape through adjustment of the inclination angles of
the boundary cells indirectly determines overall hole size: u_ri ≡ V_dr cos θ_i = f(Q_wi, T_wi); R_h = Σ_cell l_i/cos θ_i"
(the formula is quoted as printed).
* Local recession (Eq 4.87, p.370; Eq 1.18, p.45): u_n = Q/((ρC_p)_r ΔT_s). Here "θ_i is the angle of inclination of the
  normal to the surface from vertical".
* V_dr is "computed by arbitrarily averaging local Q's over all cells with boundaries inclined to the horizontal by less
  than 45 degrees" (p.370). Each update relaxes cos θ toward u_n/V_dr (Eq 6.136, p.371).
* "The outlet is forced to extend past the point on the hole where θ_i = 85°" (p.369).

**Wilkinson, Eqs 5.83-5.87 (pp.215-216):** the same consistency rule, U_sr = V_dr cos θ. Here "V_dr is calculated as an
average of predicted penetration rates for all points from a radius of one-half nozzle diameters to the point where
θ = 45°" (p.215). The boundary nodes move along the normal by ω R_vr (ω = 2), and the standoff is held fixed.

**Mechanism [inference from these equations]:**
* There is no threshold. The wall angle follows cos θ(s) = [Q(s)/ΔT_s(Q(s))] / [<Q>/ΔT_s(<Q>)]_bottom.
* With ΔT_s ∝ Q^{3/(m+3)} and m = 20, u_n ∝ Q^{0.87}, so cos θ ≈ (Q/<Q>_bottom)^{0.87}.
* The wall goes near-vertical (θ = 85°, cos θ = 0.087) where Q(s) falls to about 6 % of the bottom-average flux.
  The bottom average is taken over θ < 45°, so it is already below the stagnation peak.
* The hole radius is therefore the integral of sin θ along s. It is set by the whole decay shape of Q(s), especially its
  long tail, and not by any minimum flux.
* A flux tail that never reaches zero means the wall never becomes exactly vertical. R_h is then defined where the mesh
  ends, just past θ = 85°.
* In R Fig 6.16 the tail is still 12-16 % of peak at 14.6 R_noz. This is consistent with the wall still receding (θ < 90°)
  where the domain ends.

**W: what controls R_h.**
* A uniform change of the flux level does not change R_h: "the predicted hole radii are unchanged when the effect of a
  variable heat capacity is included. This is because the use of a variable heat capacity lowers the heat fluxes everywhere
  along the rock boundary by approximately the same amount. Therefore, the variation of heat flux with radius, and the
  resulting predicted hole shape, should not be substantially different" (W p.231).
* The shape of the distribution does change R_h: "The mass injection model changes the predicted hole radii and Stanton
  numbers because most of the effect of mass injection occurs near the stagnation point where the penetration rate … is
  highest. This changes the heat transfer distribution along the rock surface and results in a change in the predicted
  hole shape" (W p.235).
* [inference] Only the normalised shape Q(s)/Q_bottom sets R_h. That makes R_h ≈ 0.45 Z a statement about the width of
  the profile relative to standoff.

### 3.4 Hole profile shape (depth against radius)

* **No converged hole profile z(r) is plotted in either thesis: NOT FOUND.**
* Shapes are only schematic: R Fig 6.1 (p.276), Figs 6.5/1.4 (pp.299/43), Fig 6.19 (p.382); W Fig 2.7 (p.67), Fig 1.12 grid
  (p.44).
* Starting shapes were ellipsoidal (W p.216, p.226). The mock copper hole was an ellipse with 2:1 depth-to-radius, as
  "abstracted from computer simulations" (R p.407).
* Field observation: "Since the hole diameters produced ranged from 7.5 to 15 cm and few inconsistencies were observed,
  steady-state conditions presumably prevailed" (R p.399).

---------------------------------------------------------------------------------------------------------------------

## 4. Heat transfer formula used for the shape calculation

* **Flux:** there is no correlation; the flux comes from CFD plus wall functions. R Eq 6.116 is quoted in 1.1.
  W Eq 5.37 (p.171) is the same form:
  Q_w = (C_p(T_p − T_w) + ½u_p²) ρ_p C_μ^{1/4} √k_p · β / [Pr_t/κ ln(32.6 y⁺/h⁺ + 1) + St_k⁻¹],
  with St_k⁻¹ = 5.19 h⁺^{0.2} Pr^{0.44} (Eq 5.38), Pr_t = 0.86 and κ = 0.41. The driving force uses the recovery enthalpy
  (R p.345).
* **Gas temperature:**
  - The local cell gas temperature T_p from the CFD is the driving temperature, not a correlation T_gas.
  - Inlet conditions are fixed at the throat. For the small-scale runs, T_jet = 2830 K, p_jet = 4.4 p₀, Ma = 1.01 and
    Re = 20,000 (W Fig 5.4, p.175; R Table 7.1 notes p.405). Full-scale Browning runs use T_jet = 1940 K (W Table 6.1).
  - The combustion chamber was taken at the adiabatic flame temperature of about 3200 K (R p.401).
  - The wall gas runs at about 900 K: "near the rock wall, where the temperature is typically about 900 K" (W p.229).
* **Wall temperature:** T_w is the Weibull spallation temperature, R Eq 4.91 / W Eq 3.11. W refits the Weibull parameters
  (W conclusion 1, p.252).
* **Reported dimensionless forms:**
  - R: <St> = <h>/(ρC_p v)_jet, the heat transferred within the 45° limit over the jet enthalpy flux (p.383). Results are
    in Figs 6.18/1.6: <St> ≈ 0.001-0.0145, rising with T_jet and falling with pressure.
  - W: St = Q/((ρC_pU)_jet (T_jet − T_s)) (Eq 2.6-2.9, p.68). The surface average is <St> = (1/i₄₅)Σ St_i/cos θ_i over
    nodes up to θ = 45° (p.226).
* **St against SOD:**
  - R Fig 7.15 (p.430): simulated <St> at about Z/R = 22, 32, 41, 50 is about 7.2e-3, 4.1e-3, 2.3e-3, 1.6e-3.
  - The copper-calibrated <St>_eff,2 runs from about 3.1e-3 at Z/R = 24 to 2.0e-3 at 50, falling gently. The heat-balance
    values <St>_eff,1 fall from 1.6e-3 at 22 to 4e-4 at 56.
  - W Fig 2.8 (p.70) and Fig 4.16 (p.144) give the same St against SOD. The figure pages were not re-read here.
* **No fitted St(SOD) or St(position) correlation exists in either thesis: NOT FOUND.**
  - W recommends one as future work: "suitable relationships must be devised for estimating the variation of heat transfer
    with distance along the rock boundary from the stagnation point, since most correlations proposed for subsonic
    impinging jets predict only stagnation point Nusselt number values" (W p.262).
  - W Fig 6.8 (p.245) gives St against Re_noz at SOD 20: roughly linear, varying only 20 % over Re 80k-700k (W p.244).
* **The only position law used anywhere** is the exponential in surface distance for the copper calibration, R Eq 7.5
  (Section 2).

---------------------------------------------------------------------------------------------------------------------

## 5. Minimum heat flux / spallation stopping on walls

* **No minimum heat flux, flux threshold, or wall-spallation cut-off exists in either model: NOT FOUND.** The searches
  covered "minimum heat/flux", "threshold" and "cease" over the full OCR of R pp.17-460 and W pp.11-262. The wall recedes
  for any Q > 0 (Section 3.3).
* Related statements:
  - Upper-temperature cut-off, not a minimum flux: "Above a rock-specific temperature, spallation ceases to occur. Calaman
    and Rolseth attributed this to the rock becoming plastic at high temperatures" (W p.82).
  - "a high heating rate would cause premature fusion and promptly halt the process" in semi-spallable rock (R p.130).
  - Linde field data: "a higher piercing speed will produce a narrower hole for a given burner geometry. An almost linear
    increase in hole diameter is seen as drill diameter increases (Figs. 3.10 and 3.11)" (R p.127). This is consistent with
    the cos θ rule: slower advance means more wall recession.
  - Evidence that walls keep spalling far from the axis: LANL "a cavity with a maximum diameter of 1.6 m was formed by
    leaving the burner head stationary in a pre-drilled 0.6 m diameter hole (Figure 3.5)" (W p.86).
  - Onset is statistical: the global onset temperature "represents the low temperature limit of the statistical
    distribution" (W p.141). Steady-spall temperatures are more than 2× onset at fluxes above 1 MW/m² (W p.252-254).
* Wall heat loss: "heat conducted and lost to the surroundings can be almost 25 percent of the heat returned by rock spalls"
  over 10-15 min tests (R p.433, Table 7.4 p.434, Fig 7.16: R_h = 0.15 m, 1 m heated length, ΔT 500 K).

---------------------------------------------------------------------------------------------------------------------

## 6. Concave / cavity enhancement relative to flat plates

* **Wall-jet similarity (the basis for using surface distance s):**
  - "After transition of the boundary layer to full turbulence, which typically occurs within one or two radii from the jet
    axis, the velocity profile in the wall jet is fairly well documented (Fig. 6.7, after Hrycak (1981)), even for concave
    and convex impingement surfaces" (R p.306).
  - Inner-layer exponent n = 10-15 (Eq 6.37); v_max decay in Eq 6.38. Table 6.1 (p.309) gives the constants: plate
    y_{½}/b ≈ .22, δ/b ≈ .01(?); concave hemisphere .223 and .056; convex .23 and .01 (OCR, partly illegible).
* **Concave-hemisphere correlation, Livingood and Gauntner (1973a), Eq 6.43 (pp.315-316)**, valid for Z_n/D ≤ 7 and
  13000 < Re < 88000:
  - Nu_stag = 1.04 Re_D^{.554} (D/D_H)^{.3} (Z_n/D)^{.004}
  - Nu_z = 2.04 Re_D^{.586} (D/D_H)^{1.1} (Z_n/D)^{−.023} ((Z_n − z)/D_H)^{?}. The final exponent is not legible on the
    page image.
  - Nu̅ = 2.98 Re_D^{.585} (D/D_H)^{1.1} (Z_n/D)^{−.007}
  - "No local maxima in heat transfer were noticed, unlike flat plate data or that for plane jets entering troughs"
    (p.315). "Their results compare quite well with those for a semi-cylindrical trough when relative wetted areas are
    taken into account" (p.316).
* **Trapped recirculation in concave targets:**
  - "if a slot jet discharges gas toward a target with highly angled sides at somewhat larger standoff distances (Z_n/D > 8),
    the stagnation zone appears to be trapped and the flow forced to recirculate" (p.316; Fig 6.8, p.317: wedges with
    θ < 45° and Z/B > 6, troughs with Z/B = 11).
  - "Only Hrycak (1981) has noticed such behavior with a series of round jets but at Z_n/D = 8 and D/D_plate = 0.1"
    (p.318).
* **R's overall judgement:**
  - "general heat transfer behavior of a single jet exhausting into an arbitrarily shaped concave body has not been
    adequately examined … the use of flat plate heat transfer predictions for subsonic flow is a gross oversimplification
    of the complexities of flow and heat transfer in the rock cavity" (p.318).
  - "If the trough diameter was many times larger than the nozzle diameter, then the stagnation zone … should resemble that
    for a flat plate, except for a slight variation in inner boundary layer thickness" (p.316).
  - "comparable experimental data for concave targets have not, as yet, been reported" (p.412).
* **Code-predicted cavity effects:**
  - The trapped-bubble flow reversal at p_jet/p₀ > 9 and Z/R < about 20 gives off-axis peaks up to 1.6× stagnation
    (Fig 6.20, p.385). The jet is "deflected to the side of the hole by the strong shock standing at the hole bottom"
    (p.48).
  - "an optimal standoff distance is predicted at high pressures and may well exist for all supersonic flame-jets directed
    into enclosed, concave cavities" (p.48).
  - At low pressure the cavity tail holds 2-5× the flat-plate relative flux at 3-9 R_noz from stagnation (Section 1.5,
    [inference]).
  - The code "reliably predicts the correct range of heat transfer for a supersonic jet impinging on a curved, concave
    body, as affirmed by calibration" (p.446).
* **No quantitative cavity/flat enhancement factor is stated by either author: NOT FOUND.** The factors in 1.5 are mine.

---------------------------------------------------------------------------------------------------------------------

## 7. Caveats that bear on using these numbers

1. R's scheme is first-order upwind. "Numerical dispersion accounts for some of the lack of definition" (p.372). W's
   second-order scheme raises fluxes by about 30 % (W p.258). The peak width in R may be numerically broadened.
2. Rock-derived Stanton numbers are 3-4× below the CFD (R p.431). W closes the gap with refitted Weibull T_s (about 2×),
   C_p(T) (−30 %) and blowing (−35 %) (W pp.252-259). W's hole radii stay within ±5 % (W p.259).
3. The copper calibration fits are poor and time-dependent (R pp.419-427). s₀ is a fit parameter of an assumed form, not a
   measured profile.
4. The x axis of R Figs 6.16/6.20 is taken as surface distance by inference only (Section 1.6).
