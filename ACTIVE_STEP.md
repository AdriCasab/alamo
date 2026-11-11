# Active Step: none — P2-32 planning is PAUSED (the user, 2026-11-11)

Status: **no packet to implement.** P2-31 was accepted by review on 2026-11-11 and is archived (`docs/project/ARCHIVE_DONE.md` *"P2-31"*, the file verbatim). `/implement` has nothing to run until `/plan` writes P2-32.

## What P2-32 will be, once the open decision is taken

One flake-thickness rule and one surface temperature for flame and beam, then Kant's span re-derived once, the flat-floor mesh series, the beam column and three 3-D regime-boundary points re-run. Evidence that it is needed: `tests/MMWSpalling/studies/beam_column_1110/README.md` (under a beam the two-cell thermal length grows as 1/dz, so the spall / melt class flips with the cell).

## Decided (the user, 2026-11-11)

- No spalling once the surface cell passes the solidus (1,410.5 K).
- ONE rule on every face, flame and beam alike — not a beam-only branch. If the flame calibration moves, Kant's span is re-derived once.
- Under a beam the flake is measured from the surface; the peak-based depth is logged.
- Surface temperature: the face heat balance on every face; the exponential kept only as the sub-cell shape where heat enters through the face.

## OPEN — the reason for the pause, as it stands after the wider study (2026-11-11)

`Claude_markdowns/2026-11-11_flake_rule_study.md` (five rules × four versions of the hot rock law, on Kant's six legs,
the six beam-column points and a preheated flame profile; raw tables `Claude_markdowns/scripts/flake_rule_study_1111/`).

- Every profile rule is an isotherm rule; at the onset they coincide on every profile shape, and Kant's bracket
  cannot tell them apart.
- They part only on lines withdrawn at the onset that keep heating, and that difference is made by the law above
  600 °C, which is a table held flat (R·σ 400 MPa at 1000 °C, the criterion satisfied to the solidus). With any of
  three alternative tails every rule gives the same Kant behaviour and the same beam outcome. The earlier "stress rule
  refuted" check (`scripts/p2_32_plan/`) was as much a statement about that tail as about the rule.
- At the onset σ(T) is a plateau (+15 % over 750–873 K); the onset is set by the K_Ic table's fall.
- Today's two-cell rule is the outlier on flame decks: 2–4.6 × thicker on a preheated base, and in P2-31's runs the
  buckling verdict flips with the cell (45° tilt 98 / 46 / 0 % withdrawn at 2 / 1 / 0.5 mm; bore 0 / 0 / 42 %).
- Kant 2016 measured the failure-plane temperature (≈ 550 °C under a ≈ 600 °C surface): the flake is 2–4 × thinner
  than 0.4 ℓ. Unused so far; a second observable for any rule.

Decided by the user after the study (2026-11-11):
- The hot branch gets a TAIL: above ~600 °C the law is approximated from how other granites behave (no Grimsel data
  exists). Candidate curves with their sources are to be brought to the user before the packet.
- The flake constant is RE-ANCHORED on the data (Kant & von Rohr 2016: surface ≈ 600 °C before, failure plane ≈ 550 °C
  after a spall, means; Δ ≈ 100 °C on a single spall ⇒ c ≈ 0.09–0.19; Rauenzahn 0.22, Wilkinson 0.24). ⚠ Offline
  (`scripts/flake_rule_study_1111/CF_SHAPE.txt`): with the span re-derived, c 0.2 makes the beam's D 0.5 mm case at
  2 MW/m² marginal and c 0.1 makes it melt — the constant moves the beam boundary, not only the flake size.

Also decided (11-11, later): the flake's base is a FIXED temperature (the one at the line's onset); c anchored on
Kant 2016's MEANS (0.09; 0.19 the bracket); hot-tested data for the tail; the onset level stays.

Brought to the user (`Claude_markdowns/2026-11-11_tails_and_offline_map.md`): tail T1 proposed (expansion 7e-6/K above
650 °C from dilatometry; E/E0 held; brittle limit 1000 °C), T2 / T3 the brackets; the offline map: under the fixed-base
rule NO tail changes a cell; the constant decides it — 7 spall cells at c 0.40, 2 at 0.19, 1 at 0.09 (D 0.5 mm at
8 MW/m²). Awaiting the user's confirmation of the constant given that consequence, then the packets.

**ALL OPEN CHOICES NOW TAKEN (the user, 11-11): tail T1; production constant 0.19, bracket 0.09 (every beam result
reported at both); fixed base temperature; brittle limit 1000 °C; onset level unchanged.** Next: `/plan` writes P2-32
(the reading, the rule and the surface temperature at the present constant) and then P2-33 (the constant, the span, the
tail). Earlier wording, kept for the record: Tail T1 confirmed; The constant: the user asked for the difference between 0.09 and 0.19 before
confirming (50 ms pyrometer averaging biases the mean pair low by ~0.015; Wilkinson 0.24, the Weibull median ~0.22 and
Kant's single spall 0.19 cluster; on Höser's leg 0.09 / 0.19 / 0.40 give 80–120 / 170–250 / 360–530 µm flakes).

## Carried for the planner (P2-31 review)

- DR-3: `witness/REFERENCE.md`'s `## Current reference` block still leads with `face_p2_28.md5`; the live reference is `face_p2_31.md5`. `CLAUDE.md`'s gate example still reads `p2_30c` (a `CLAUDE.md` edit needs the user's approval).
- DR-1 (the phi half of the overlap rule has no RED-on-old witness), DR-8 (`retired_keys` K31's accepted side has no leg) — candidates for P2-32's unit clauses.
- Nothing is committed: P2-26 … P2-31 (HEAD `a237c13eb`).
