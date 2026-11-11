# Candidate tails for the rock law above 600 °C, and the offline regime map under the decided flake rule (2026-11-11)

Follows `2026-11-11_flake_rule_study.md`. The user's decisions (11-11): the law gets a tail above 600 °C built from how
other granites behave; the flake's base is a FIXED temperature (the one it has at the line's onset); the flake constant
is re-anchored on Kant & von Rohr 2016's mean temperatures; hot-tested data are preferred for the tail. Zero compute.
Scripts and raw output: `Claude_markdowns/scripts/regime_map_offline_1111/` (`thermal.py`, `map.py`, `MAP.txt`).

## 1. What other granites do above 600 °C

| property | measurement | source | as coded |
|---|---|---|---|
| expansion | total strain 0.87–1.08 % at 500 °C, 1.32–1.68 % at 700 °C, 1.58–1.91 % at 1000 °C (six of eight granites; two outliers) ⇒ only +0.14 to +0.29 % from 700 to 1000 °C, i.e. **5–10 × 10⁻⁶/K** | Alzahrani, Lasheen & Rashwan 2022, *Materials* 15:2041, Table 4 (dilatometer, 30–1000 °C, unconfined, includes crack opening) | 28.4 × 10⁻⁶/K held to 1000 °C, then the strain frozen |
| expansion, shape | sharp rise near 570 °C; a second rise 700–900 °C in four of the eight; Eibenstock: steps at 573 and 870 °C | same; Wang et al. 2020 *Acta Geotechnica* (abstract) | none |
| stiffness, tested hot, unconfined | E ≈ 29 / 29 / 24 / 21 / 17.5 / 8.5 / 8.5 GPa at 25 / 200 / 300 / 400 / 500 / 600 / 700 °C ⇒ E/E0 1.0 / 1.0 / 0.82 / 0.70 / 0.60 / **0.29 / 0.29**: a drop between 500 and 600 °C, then "no significant change" | Zhang, Zhang, Yang & Cai 2022, *Geofluids* 9426378, Fig. 6 and 8 (values READ FROM THE FIGURE, ±1 GPa); text: brittle-to-plastic transition between 600 and 700 °C | 0.65 / 0.58 / 0.50 / 0.50 at 500 / 550 / 600 / 650 °C, held |
| strength, tested hot, confined | Westerly, 130 MPa: peak stress 815 MPa at 600 °C → 719 MPa at 900 °C; brittle to 900 °C, ductile at 1000 °C | Violay et al. 2017, *Sci. Rep.* 7:7705 (via the 11-02 research note) | — |
| toughness, tested hot | data stop at 500 °C (50–76 % of the cold value) | Yang et al. 2022 (11-02 note) | 0.60 held from 600 °C |

⚠ The one hot stiffness dataset puts its drop between 500 and 600 °C, i.e. at the onset temperature itself. The tails
below start at 600 °C and leave the onset and R alone; that dataset also bears on the table below 600 °C (not acted on).

## 2. Three candidate tails

| | expansion above 600 °C | E/E0 above 600 °C | brittle limit | R·σ at 600 / 700 / 800 / 1000 °C (MPa) | onset criterion satisfied up to |
|---|---|---|---|---|---|
| coded | 28.4e-6 to 1000 °C | 0.50 | the solidus | 170 / 228 / 285 / 400 | the solidus |
| **T1** | 7e-6 from 650 °C | 0.50 | 1000 °C | 170 / 195 / 210 / 238 | the limit |
| T2 | as T1 | 0.30 from 650 °C (Zhang) | 1000 °C | 170 / 117 / 126 / 143 | 885 K |
| T3 | as T1 | linear to 0.10 at 1000 °C | 1000 °C | 170 / 156 / 126 / 48 | 960 K |

Proposed: **T1** — every number in it has a measurement behind it (the expansion; the brittle limit; "no significant
change" in hot stiffness from 600 to 700 °C). T2 and T3 are the brackets.

## 3. The offline map (P2-25's grid; a 1-D column, flaws off; onset 843 K)

Span re-derived on Kant's bracket for each constant: 12.02 / 5.73 / 2.73 mm at c 0.40 / 0.19 / 0.09 (flake base at
662 / 748 / 796 K under an 843 K surface). S = spalls at the onset; S@T = spalls later at surface temperature T;
M = melts first; U / C = undecided / cold by 600 s. Rows D ≥ 2 mm are M (or U / C at the lowest fluxes) in every case.

| absorbed flux (MW/m²) | 0.3 | 0.5 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|---|
| **today's form** (c 0.40, base at a fraction of the rise, coded law) — D 0.5 mm | M | M | S@980 | S | S | S |
| — D 1 mm | U | M | M | S@1120 | S | S |
| **decided form, c 0.19** — D 0.5 mm | M | M | M | M (bk 1.08) | S | S |
| — D 1 mm | M | M | M | M | M | M (bk 1.13) |
| **decided form, c 0.09** — D 0.5 mm | M | M | M | M | M (bk 1.15) | S (bk 0.92) |
| — D 1 mm | M | M | M | M | M | M |

For comparison, P2-25's simulation (3-D, 1.5625 mm cells, today's code) spalled in 10 cells: D 0.5 and 1 mm from 1 to
8 MW/m², D 2 mm at 4 and 8.

- **Under the fixed-base rule the tail changes no cell**: T1, T2, T3 and the coded law give identical outcomes (a line
  withdrawn at the onset only gets thicker). The tail matters for what is logged and for sanity, not for the map.
- **The constant decides the map**: 7 spall cells at 0.40, 2 at 0.19, 1 at 0.09. Reason: the span Kant's bracket
  allows a thin flame flake is short (2.7 mm ⇒ the thickest flake that can buckle is 0.21 mm), and under a beam the
  layer within ~50 K of the surface temperature is 0.19–0.9 mm thick at D 0.5 mm (flat-topped profile).
- Consistent with the published millimetre-wave record on granite (melt only). Limits: a 1-D column; Kant's two
  temperatures are means over spalls on Gotthard and Bethel, with overheated flakes in the "before" value (a smaller
  true difference would mean a smaller constant still); no Grimsel pair exists (Kant 2017 reports the failure plane
  only).
