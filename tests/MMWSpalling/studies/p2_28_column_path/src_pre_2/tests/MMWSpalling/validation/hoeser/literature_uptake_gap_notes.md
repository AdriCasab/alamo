# Literature on the heat-uptake gap at a spalling face (advisor, 2026-10-26)

Question: why does a spalling granite face take in 3–6.5x less heat than a calorimeter under the same jet
(Höser 0 MPa bound 4.7–6.5x, P2-7), and what should the cut depend on? Three online searches (beam data,
flame-jet follow-ups, gas side). "read" = full text opened; "abstract" = abstract/snippet only.
Downloaded: scratchpad `gasside/` (Terekhov 2013, Halouane 2015) — session-local, re-fetch if needed.

## The gap is old, measured three times, and still open
- **Rauenzahn 1986 MIT PhD** (hdl 1721.1/14884, §7.3–7.5 pp. 426–447, read via OCR): copper calorimeter 3–6 MW/m²
  (peak 5.3–7.4), agrees with his CFD to 20 %; rate-implied flux ~4x lower; Browning field holes the same 3–4x.
  Ruled out: wall losses ≤ 25 %; a hotter surface gains < 2x (St rises too); judged T_s > 600 °C "improbable"
  (recovered chips show no phase change; Barre showed mineral glow). Proposed: partly ADHERED spalls held by
  "glue-like" FUSED BIOTITE trapping semi-stagnant gas over a shielded face; chip turbulence; particles absorbing
  heat; fractured rock. Verdict: the factor 3–4 "remains unexplained"; a smooth metal calorimeter on fractured
  rock is "unreliable".
- **Dey & Kranz 1985** via Dreesen & Bretz (NETL Phase I report ~2003, read): only **15 %** of available heat
  reached the rock at Pedernal (≈ 6.7x — same size as Höser). Efficiency 0.09–0.125 (mechanical Weibull) /
  0.27–0.47 (thermal Weibull). Barre flame-jet specific energy 8.7–18 kJ/cm³.
- **Wilkinson & Tester 1993 RMRE** (on disk): predicted St 3–5x field-implied (p. 34); thermal Weibull halves it.

## Rock side
- **Wilkinson 1989 MIT PhD** (hdl 1721.1/14031, §6.2 + Ch. 7, read): steady T_s > 2x onset above 1 MW/m²
  (ductile, stress-relieving yield ⇒ "overheating"); spalls leave the stagnation zone at **~13 Hz** (only
  release-frequency number found).
- **Kant 2017 ETH PhD** (doi 10.3929/ethz-b-000216571, read): Fig. 4.19 fits HÖSER's rates with gas T 1,200 °C
  (fitted, vs ~1,600 °C nozzle), Θ_SP 550 °C, Θ_M 1,100 °C, Δx 200 µm; melting-capped rate u_M = κ/Δx·ln(Θ_M/Θ_SP).
  Calorimeter held near rock T so h transfers. **Inverse-IR method for flux AT the rock surface** (Kant & von Rohr
  2016 IJHMT 99:1) — shown on metal only, never on spalling rock: the direct measurement of uptake nobody has made.
- **Gahan et al. 2001 SPE-71466** (abstract): laser minimum specific energy "just prior to the onset of melting" ⇒
  melting caps any overheating route.
- **Seagraves ARMA-2019-1799** (abstract only; paywalled): continuous thermally activated recession + GB-flaw onset;
  three scales set onset T, T_s(q) and rate(q), and a melting-limited rate; fits Wilkinson's Barre/Westerly.
- **Hu, Song, Liu, Li, Shen, Lyu 2018 JPSE 171:1390** (abstract): 280 spalls, aspect ratio 6.91.
  ⚠ our `spall_thickness/REFERENCES.md` calls this "Shen et al." — first author is HU.

## Gas side (sizes at Höser's conditions)
- **Blowing:** St/St0 = ln(1+2λB′)/(2λB′) (Kays & Crawford λ 0.5; Moyer & Rindal λ 0.4; via Padovan 2023
  arXiv:2310.07080 Eq. 36, read). 1 g/s of spalls over Ø40 ⇒ B′ ≈ 0.15 ⇒ −6 to −7 % (−15–18 % if stagnation-
  concentrated), and that counts solid spall mass as gas (true displacement B′ ~1e-5): an upper bound.
  Wilkinson 1989 modelled spall blowing −35–40 % (full-blowing limit St/St0 ≈ 0.12); Roberts, East & Pratt 1983
  (unseen): sand-laden hot boundary layer > 75 % loss at solids:gas ≈ 1. Höser's global ratio ≈ 0.07.
- **Particle-laden impinging jets mostly ENHANCE stagnation h** (Yoshida 1990 IJHMT 33:859 up to 2.7x; Shimizu 1980
  up to ~6x; abstracts) — but those particles arrive with the jet, not leave the wall.
- **Hot-wall gas properties:** variable c_p −30 % (Wilkinson 1989); Rothenfluh 2013 ETH (read): h depends strongly
  on target surface T for SCW jets. ⚠ Höser checked h independent of sensor T (180–550 °C, L > 9).
- **Cavity:** **Terekhov & Kalinina 2013** (doi 10.1299/jsmeicjwsf.2013.4._1089-1_, read): air jet d0 8.9 mm into a
  SHARP-edged spherical cavity Dc 46 mm (≈ Höser's geometry), depth/Dc 0.13–0.5, S/d0 0–10: h at r < Dc falls
  "almost twice" vs flat (toroidal vortex at the mouth); a rounded edge restores it. Terekhov 2009 IJHMT 52:2498
  (abstract): lower local h "almost fully compensated" by area at depth ≤ 0.5 Dc. Halouane 2015 (read): bottom acts
  flat for L ≤ 12d. ⚠ isothermal air; in a flame the trapped gas is hot.
- **Roughness** 0 to +6 % (Beitelmal 2000); **emissivity/catalysis** ≤ 10–12 % (Baukal & Gebhart 1997 ETFS 15:323).

## Beam data (no gas boundary layer)
- **DeepU — Słupski et al. 2025** Research Square doi 10.21203/rs.3.rs-7725822/v1 + EGC 2025 (read): 26–27 kW fibre
  laser, Ø50 mm, granite, 30 s; spall T ~550 °C; spalls 0.2–1 mm, ~0.5 mm thick; SE 6.35 kJ/cm³ ⇒ 3.5–4.6x on
  incident, 1.9–2.5x with absorptance 0.544 (Seo 2022 Sci Rep 12:680). ⚠ internally inconsistent (P·t/SE ⇒
  ~35 mm mean Ø vs stated ~80 mm; an 80 mm hole would mean ~0.9x). "Trapped spalls … partially melt … accumulate
  on the bottom" and stop penetration.
- **Reed et al. 2003 ANL/TD/TM03-02** (read): sedimentary; debris/exsolved gas "remains in the beam, parasitically
  absorbing energy". Graves & O'Brien 1998 SPE 49259 ~5x on incident (sandstone). Leong 2003 plume 10–20 % (snippet).
- **No MMW spallation energy balance exists** (Woskov & Michael 2011 PSFC/JA-11-24: 28 GHz, melt/vaporisation only;
  solid granite 1/e depth 7.1 cm ⇒ volumetric). A gap in the literature this project could fill.

## What the discriminating trends say (my reading)
| candidate | zero before 1st spall (Kant) | smaller cut at 47 MPa (P2-7 a1′) | rate rises with time | size at Höser |
|---|---|---|---|---|
| overheated / adhered spalls (rock) | yes | yes, if loaded spalls leave faster | neutral | ≤ ~2x before melting |
| cavity vortex (gas) | yes | NO — the 47 MPa bowl is deeper | NO — worsens with depth | ≤ ~2x local, offset by area |
| particle loading / blowing (gas) | yes | NO — more spalls at 47 MPa | neutral | tens of % |
| hot-wall gas properties (gas) | no (acts always) | neutral | neutral | ≤ 30 %, partly in her calibration |
⇒ several partial causes multiply to the observed 5x; only the rock-side spall behaviour matches the LOAD trend, so it
is the primary candidate; the cavity is the largest gas-side correction and is computable from the hole shape.
