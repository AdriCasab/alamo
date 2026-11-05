# P2-9 offline geometry checks (planner, 2026-10-27) — evidence for rev 2

Run: `/Users/tzetze20/Desktop/code/.venv/bin/python <script>.py > <script>.txt` (numpy + scipy).
Binary voxel staircases of planes, spheres, a cube and a cylinder; no ALAMO run.

| script | measures | headline |
|---|---|---|
| `plane_line_vs_cdiff.py` | per-cell line area Σ 1/μ on tilted planes; ±1-cell normal (P2-8) vs Gaussian-smoothed normal | P2-8 normal +3 to +12 %; reproduces the implementer's a2 (vi) +8.6 / +5.2 / +10.0 % to within ~1.7 points |
| `line_area_shapes.py` | per-cell line area Σ 1/μ on curved/sharp shapes | **+4.1 % spherical cavity, +4.6 % bore**: the per-cell line area is plane-only and is NOT used for heat |
| `area_rules_s10_s15.py` | P2-8's per-face rule w = 1/\|n\|₁ vs a projected rule, normals cdiff / σ 1.0 / σ 1.5 (±3) | per-face rule + smoothed normal: planes ≤ ~1.8 %, spheres ≤ 1.2 % |
| `area_rules_s075_s10_trunc2.py` | the same, σ 0.75 / 1.0 truncated at ±2 cells (5×5×5) | **chosen: σ 1.0 ±2** — planes −1.8 to +0.1 %, spheres −1.2 to −0.2 %, bore −2.8 %, cube −5.0 % |
| `line_recession_sphere_exact_vs_smoothed.py`, `line_vs_face_sphere_by_angle.py` | one simultaneous firing's volume Σ 1/μ on spherical cavities, by angle | +5.5 / +2.4 / +1.5 % at R 12 / 24 / 36 cells even with the EXACT normal (falls as dx/R) |
| `a4_predictions.py` | the a4 per-tilt predictions at σ 1.0 ±2 | area error −1.1 / −0.1 / +0.4 / 0 / −1.6 / 0 % at 20 / 30 (=60) / 35 / 45 / 35-az45 / (1,1,1) |

⚠ **The one-firing volume excess does NOT bias the steady recession rate.** Heat in = flake
enthalpy out, and a line's flake carries enthalpy in proportion to the volume it removes, so the
volume rate is heat / (ρ c ΔT_flake) whatever the per-firing volume. The plane's (or sphere's)
mean normal speed error is therefore the AREA error alone. The excess changes how often cells
fire, not how fast the surface recedes. (Planner's first reading, corrected the same day.)
