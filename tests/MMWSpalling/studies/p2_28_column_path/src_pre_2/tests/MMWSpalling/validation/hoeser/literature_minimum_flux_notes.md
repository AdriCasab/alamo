# Literature on the minimum-flux / flux-dependence question (advisor, 2026-10-14)

Question: has anyone else met the problem that a surface-temperature spall criterion carries no heat-flux
dependence (no minimum flux, no steady-state surface-T rise), and what length did they use? Found: the
flux dependence in the literature comes from TWO mechanisms of OPPOSITE sign, and every analytical model
references the stress to the far-field/initial temperature exactly as ours does.

## A. Models that reference stress to the far field and get flux dependence from the STRESSED VOLUME
- **Dey & Kranz 1985/1987 → Rauenzahn & Tester 1989 → Wilkinson & Tester 1993 → Walsh & Lomov 2014**
  (LLNL-TR-632239 §4, Eqs 14–19, on disk `franck_material/Walsh_2013_…pdf`): Weibull "median condition"
  G(σ) = 0.5 over the cone volume that would remove a surface point (disk spalls, aspect ratio C_L);
  quasi-steady "subliming" profile **T = ΔT_s·exp(−v x/α)** (ℓ = α/v); **σ = βEΔT/(1−ν)** with ΔT to the
  far field; closed form Q(v) and **T_s = T_r + Q/(ρ c_p v)**. ⇒ steady-state T_s RISES with Q (the stressed
  volume shrinks as ω = vm/α grows) — the Wilkinson rise — but **no minimum flux** (at low Q the volume
  grows and spalling gets EASIER). Walsh & Lomov 2014 (IJRMMS 70:375): spall aspect ratio is SIZE-DEPENDENT
  (C_L 1 below 100 µm, 5 above, field + grain-scale sims); putting that FIXED LENGTH into the Weibull model
  "changes the predicted surface temperature significantly … as the scaled thermal penetration distance
  approaches the spall thickness associated with the aspect-ratio transition". ⇒ a fixed physical length is
  what carries flux dependence; the same structural fact as our constraint depth.
- **Kant et al. 2017 JGR** (on disk): LEFM Sp with a = 20 µm, full lateral constraint, stress to the initial
  T — the group that MEASURED the minimum flux (Kant 2016) did not put it in their criterion either.
- **Seagraves, ARMA-2019-1799** "A mechanistic model for thermal spallation of geological materials"
  (paywalled, onepetro): continuous thermally-activated surface recession + fracture-mechanics onset at
  GRAIN-BOUNDARY flaws (grain size = the model length); "three characteristic scales control the onset
  temperature, the dependence of surface temperature and penetration rate on heat flux, and the limiting
  rate at melting"; matches Wilkinson & Tester's rate AND surface T vs flux for Barre and Westerly, the
  difference attributed to grain size. ⚠ Closest published analogue to our chain; get the full text.
- **Lyu, Song, Li 2019 Geothermics 77:99** — analytical inversion of spallation temperature and Weibull
  heterogeneity from drilling data (the Rauenzahn family).

## B. The minimum flux: a THIN-LAYER BUCKLING / bifurcation condition (opposite sign to A)
- **Preston & White 1934** (mechanism); **Thirumalai 1969/1970** (ARMA-69-0705): "in good spalling rocks a
  thin heated layer with high thermal gradient forms" (quartzite ~1 mm, > 200 °C/mm); "inelastic behaviour
  inhibits thermal spalling"; **Germanovich 1997**: buckling of the compressed thin layer.
- **Dey & Kranz 1985/1987 "employed plate-buckling theory to BOUND THE HEAT FLUX REQUIRED to generate
  spalls from a granite surface"** (Walsh, Lomov & Roberts GRC 2011 §2, on disk) — the one published
  derivation of a minimum flux: the heated plate must be thin relative to its width to buckle.
- **Browning et al. 1965** (via Meier, Brkic et al. WGC 2015, on disk): "if the heat flux is too low the
  induced stress field is not strong enough to fracture the rock; if too high the heat diffuses faster than
  the spalls are removed, the profile flattens and the stresses are diminished ⇒ an optimal incident heat
  flux exists for every rock". The ETH group states spallation conditions as TWO numbers: heat flux ~1 MW/m²
  AND surface T ~500 °C (citing Rauenzahn).
- **Kant & Rudolf von Rohr 2016** (on disk): the measurement, h_min ≈ 500 W/m²K, q ≈ 0.3–0.5 MW/m².
- **LLNL grain-scale GEODYN sims (2011–2013, on disk):** fractures start at grain boundaries; spalls "one or
  more grain diameters in extent, thinner than a grain"; flat interfaces withstand "substantially higher heat
  fluxes" ⇒ "classical plate-buckling theory in a half-space may OVER-estimate the stress required";
  "the traditional buckling model is essentially a bifurcation event — small fluctuations have a large
  impact on onset"; "simple models based on idealized geometries and buckling modes are insufficient;
  Weibull models are empirical".
- A third route to a minimum flux: RATE-DEPENDENT INELASTICITY (Thirumalai's remark; our relief factor
  would become rate-dependent). Less constrained; not pursued.

## C. The same structure in other fields
- **Thermal-shock theory** (Kingery 1955; Manson 1966; Hasselman; Lu & Fleck 1998 Acta Mater 46:4755):
  σ_max = Eβ ΔT·ψ(Bi)/(1−ν), ψ → 1 at high Biot, **ψ ≈ 0.31·Bi at low Biot** (plate), Bi = h·L/k with L
  the body size ⇒ the stress under slow heating is ∝ q·L/k — our soaked-limit drive d_ref·q/k with the body
  size replaced by the spall's constraint length. In a half-space there is no L, which is why a surface-T
  criterion has no flux dependence; the rock supplies L through the flake.
- **Concrete explosive spalling** (NIST Phan; PMC8541664; S2214509521002199): criteria written on the
  GRADIENT (spalls at 14–21 °C/mm) and the heating RATE (HSC spalls at 1 °C/min, not at 0.1); mechanism
  mostly pore pressure, thermal stress secondary.
- **Laser rock spallation** (Xu & Reed ANL 2002; Xu & Graves): a power-density WINDOW — ≳ 0.8–1 kW/cm²
  (8–10 MW/m² at short exposure) to spall sandstone/shale, melting above.

## D. What it means for P2-2e
Nobody writes a "reference at a constraint depth"; the published minimum flux is Dey & Kranz's buckling
bound (plate thickness vs width), which our constraint depth proxies with one length (≈ half the buckle
wavelength of the thickest flake). The Weibull-volume models (A) give the steady-state T_s rise but no
minimum; our per-cell Weibull draw carries a weak form of (A). The thermal-shock ψ(Bi) is the formal
cousin. The one paper to obtain: Seagraves 2019.
