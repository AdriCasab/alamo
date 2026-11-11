# The flake rule: what the candidates read, where they part, and what is underneath (2026-11-11)

Written after the user stopped a second quick switch of rule ("have a more comprehensive look ... there is something
deeper here"). Zero compute. Scripts and raw tables: `Claude_markdowns/scripts/flake_rule_study_1111/` (`study.py`,
`STUDY.txt`); the two earlier single-question checks are in `scripts/p2_32_plan/`.

## 1. Method, and what it was checked against

One 1-D explicit column (ρc_p(T), λ(T) as coded) supplies temperature profiles: Kant's six bracket legs to 200 s, the
six beam-column points to the solidus, and the flat fixture's flame on a base preheated to 293 / 450 / 550 K. Every
rule and every version of the rock law is evaluated on the SAME stored profiles.

Checks: C_L(851.21 K) = 12.727 and C_L(758.15 K) = 15.243 (both logged values, to the digit); R·σ(850 K) = 169.0 MPa
(the UCS); the onset criterion reaches 1 at 842.7 K (the code fires at 845.5–846.0 K); Kant's six logged spans within
1.5 %; the fixture's first firing at 2.37 s (code 2.38 s at 0.5 mm). Precision ≈ 1 %: it cannot resolve today's
0.15 %-wide span window.

Rules: TAN 0.4 × the tangent length kΔT/q; RISE the depth where the rise falls to e^-0.4 of the surface rise; ISO the
depth of a FIXED isotherm (664 K, where RISE sits at the onset); SIG the depth where the stress falls to 0.75; WB / WBs
Rauenzahn & Tester's Weibull median with and without the z² (cone) weight, m = 13. Each rule gets its own span from
Kant's bracket, as a re-derivation would give it.

Law versions (the stress σ = R·E(T)·[ε(T) − ε(T0)]/(1 − ν)): V0 as coded (Table 3 held flat above 600–650 °C);
V1 the expansion falls back after the quartz inversion (28.3 → 10 × 10⁻⁶/K by 650 °C); V2 the stiffness keeps falling
(E/E0 0.50 at 600 °C → 0.25 at 900 °C); V3 both. ⚠ V1–V3 are ILLUSTRATIVE brackets, not measurements: no Grimsel
expansion curve and no hot stiffness above 500 °C exist (literature report, 2026-11-02). They leave everything below
575 °C, hence the onset and R, unchanged.

## 2. Findings

**F1. Every profile rule is an isotherm rule.** With a local stress law, "where the stress has fallen to x %" is
"where the temperature has fallen to T*(T_s)". RISE, ISO and SIG differ only in how T* moves once the surface is hotter
than the onset. At the onset they coincide on EVERY profile shape (beam column at 846 K: 0.75 / 1.18 / 1.96 / 3.41 /
6.21 / 1.63 mm under all three). Kant's bracket, which is all first firings on exponential-like skins, fixes one
number for all of them and cannot tell them apart.

**F2. The rules disagree only about lines that were withdrawn at the onset and keep heating, and that disagreement is
made by the law above 600 °C.**

| | V0 (as coded) | V1 | V2 | V3 |
|---|---|---|---|---|
| R·σ at 600 / 727 / 1000 °C (MPa) | 170 / 243 / 400 | 167 / 195 / 250 | 170 / 192 / 200 | 167 / 154 / 125 |
| onset criterion satisfied for T_s in | 843 K → solidus | 843 K → solidus | 843 K → solidus | 843–957 K only |
| Kant's no-spall legs that would fire later (span on, 200 s) | SIG, WBs: KG140, KG150, KB160; RISE, ISO, TAN: none | none, any rule | none, any rule | none, any rule |
| beam, D 1 mm at 2 MW/m² | RISE spalls at 1126 K; SIG 946 K; WBs 926 K; ISO melts | melts, any rule | RISE spalls at 1126 K; others melt | melts, any rule |
| beam, D 2 mm at 4 MW/m² | SIG spalls at 986 K; others melt | melts | melts | melts |
| beam, D 0.5 mm at 2 MW/m² | spalls at the onset (WBs at 886 K) | spalls at the onset (WBs 986 K) | spalls (WBs 906 K) | spalls (WBs: marginal melt, bk 1.01) |
| beam, D ≥ 2 mm at 2 MW/m² | melts, any rule | melts | melts | melts |

As coded, the rock carries 400 MPa of elastic compression at 1000 °C (2.4 × its cold strength) and the onset criterion
stays satisfied up to the solidus. Every "hot spall" in the table comes from that tail. With any of the three
alternative tails all rules give the same Kant behaviour and the same beam outcome. My earlier "the stress rule is
refuted" was a statement about V0's tail as much as about the rule.

**F3. At the onset the law is a stress plateau, and the onset is set by the toughness table.** R·σ is 137 / 148 / 160 /
165 / 169 / 170 MPa at 700 / 750 / 800 / 820 / 846 / 873 K: +15 % over 120 K. Over the same range K_Ic falls from 0.76
to 0.60 of its cold value, which is what takes the score from 0.73 to 1.06. The log-slope of σ(T) is 0.36 at 846 K
(1 for a linear law). So a rule that asks "where has the stress fallen by x %" is ill-conditioned at the onset itself
(SIG's 0.36 → 0.40 → 0.43 → 0.38 ℓ between 800 and 900 K), and the onset level is the crossing of two estimated tables.

**F4. Today's two-cell rule is the outlier on flame decks too.** On the fixture's flame over a preheated base it reads
0.57 / 1.26 / 2.61 mm at T_b 293 / 450 / 550 K (2 mm cells) where TAN reads 0.57 throughout and RISE / SIG 0.57 / 0.61 /
0.67 mm; the 12 mm span withdraws above 0.93 mm. And in P2-31's own runs the buckling verdict flips with the cell:

| leg (P2-31 mesh test) | 2 mm | 1 mm | 0.5 mm |
|---|---|---|---|
| 45° tilt, lines withdrawn by B | 97.9 % | 45.7 % | 0 % |
| bore, lines withdrawn by B | 0 % | 0 % | 41.7 % (all on the wall) |
| Höser 10 s, lines withdrawn by B | 6.6 % | 0 % | 0.1 % |
| flat floor | 0 % | 0 % | 0 % |

The fixture runs at 65–83 % of the buckling limit (bk 0.83 / 0.70 / 0.65 on the flat floor), so a 10–20 % error in the
thickness reading decides whether a tilted or curved face spalls. ⚠ The profile length ℓ_p drifts with the cell as
well on these legs (TL45 2.29 / 2.08 / 1.94 mm), so a different rule alone does not cure this; a reading that converges
does. Hypothesis, not shown: P2-31's unconverged Höser hole shape and bore rate are this.

**F5. A measurement of the flake's base temperature exists and was never used.** Kant & von Rohr 2016 resolve single
spalls with a 1 kHz pyrometer: surface ≈ 600 °C just before, failure plane ≈ 550 °C just after (means, both rocks;
Δ ≈ 100 °C on a single spall, Fig. 9–11). The model's flake base is ≈ 180 K below its surface (0.4 ℓ). In ℓ units the
measurement is 0.09–0.19. Rauenzahn ("backside at ~80 % of the surface rise", 0.22 ℓ), Wilkinson (0.4 mm at 6e-4 m/s,
0.24 ℓ) and the Weibull median at m = 13 (0.07–0.26 ℓ) say the same. Only span / c_f is calibrated, so this does not
move the minimum flux; it moves the flake thickness, and it gives every rule a second observable to be scored on.

## 3. What this means

- The question "stress profile or temperature profile" was the wrong first question. Within the isotherm family the
  choice is nearly immaterial at the onset and on Kant; it matters only above ~900 K, where the model has no physics,
  only a table held flat.
- The quantity that decides the spall / melt boundary under a beam for absorption depths around 1–2 mm is what the
  rock does between 600 °C and the solidus: whether stress keeps building, whether the criterion stays on, how much
  strain the buckling coefficient sees. The 14 October direction listed exactly this ("softening toward the solidus")
  and it was never built. "No spalling past the solidus" (decided 11-11) is its crudest form.
- Robust across every rule and every tail tried: D 0.5 mm spalls at the onset (WBs marginal); D ≥ 2 mm at 2 MW/m² melts.
- The same thickness reading, near the buckling limit, is behind mesh-dependent removal on tilted and curved flame
  faces.

## 4. Open, and what would close it

1. The hot branch: stress, toughness and strain from 600 °C to the solidus. Evidence on disk: dilatometer curves spike
   near 570 °C and fall back; confined Westerly stays brittle to 900 °C and is ductile at 1000 °C; cooled specimens lose
   two thirds of their stiffness by 800 °C; the MIT theses' only cut-off is an upper, plastic one; Kant's operating
   range ends at "melting, ≈ 1000 °C". Kant's no-spall legs (surface well above the spall temperature, no spall) are
   the one flame constraint on it.
2. Which isotherm once the surface is hotter than the onset (RISE: a fraction of the rise; ISO: fixed). Kant cannot
   tell; a withdrawn-then-fired observation could.
3. c_f itself (Kant's failure-plane temperature says the flake is 2–4 × thinner than 0.4 ℓ).
4. Whether thickness tracks the heated layer at all at the thin end ("one to two grains", the aspect-ratio collapse at
   ~400 µm) — the assumption the buckling minimum flux rests on.

## 5. Addendum: the flake constant is not only roughness (`scripts/flake_rule_study_1111/cf_shape.py`, `CF_SHAPE.txt`)

Kant's bracket fixes span / c, so on a flame the minimum flux does not depend on c. Under a beam it does: on a
flat-topped profile an isotherm closer to the surface temperature sits proportionally deeper than on an exponential.
With the span re-derived for each c (12.0 / 6.0 / 3.0 mm at c 0.4 / 0.2 / 0.1):

| beam point | c 0.4 | c 0.2 | c 0.1 |
|---|---|---|---|
| D 0.5 mm, 2 MW/m²: bk at the onset | 0.80 (spalls) | 1.05 (marginal: spalls at 906 K only under RISE with the coded tail) | 1.44 (melts under every rule and tail) |
| D 1 mm, 2 MW/m²: bk at the onset | 1.26 | 1.68 | 2.32 |
| fixture's flame, flake at the first firing | 0.571 mm | 0.285 mm | 0.142 mm |
| flake base temperature under an 846 K surface | 664 K | 746 K | 793 K (Kant measured ≈ 823 K under ≈ 873 K) |

So the three open choices each move the beam boundary: the hot branch decides D ≈ 1 mm, the isotherm above the onset
decides the same case under the coded tail, and c decides D ≈ 0.5 mm. Across them the boundary at 2 MW/m² lies anywhere
from below 0.5 mm to between 1 and 2 mm.
