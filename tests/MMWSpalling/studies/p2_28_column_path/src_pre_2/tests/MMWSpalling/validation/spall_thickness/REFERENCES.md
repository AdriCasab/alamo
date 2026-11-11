# Spall thickness — the literature this law is compared against

Supports P1g-3 `(a0)`, `(a2)`, `(a4)` and `(a6)`.

## ⚠⚠ Provenance of this note, stated first

These numbers were **extracted in the 2026-10-04 session** and are recorded with their
search/fetch status in `Claude_markdowns/2026-10-04.md` §1 (lines 5–40). **Three of the
four primary sources are paywalled and were not fetched**; what is recorded is the
journal-level citation plus the figure/table pointer and the extracted number, which is
what "landed" means for this step. Two PDFs *were* obtained that session (Kant SGW 2017;
Walsh/Lomov/Roberts LLNL 2012) and Rauenzahn's thesis OCR is on disk.

**No number in this note was produced by me reading a PDF in the current session.** Where
a pointer is a journal volume/page rather than a figure number, that is said explicitly
below rather than dressed up as a figure citation.

## (a0) THE CONVERSION — size is not thickness

⚠⚠ **Every field figure below is a SIZE** (a sieve diameter, a PSD bin, a granulometry
peak). **This law predicts a THICKNESS.** Comparing them directly is wrong by up to 6×,
and the error sits in the comparand where no later criterion can see it.

**The conversion used, and its source:**

| size | aspect ratio | source |
|---|---|---|
| ~200 µm | **1–2 : 1** | Walsh et al. 2014, **Fig. 4** — equant particles, tied to grain size |
| large spalls | **~6 : 1** | Walsh et al. 2014 — the commonly quoted value, for *large* spalls only |
| size-dependent fit | **6.91** recommended | Shen et al. 2018, from **280 measured spalls** |

`flake_thickness.py` interpolates log-linearly between the 200 µm anchor `(1–2):1` and a
`6:1` anchor at 1 mm, and carries a **range** rather than a point. ⚠ Walsh is recorded two
ways in the source log — "6:1 → 1:1 below **400 µm**" (§1) and "falling to 1–2:1 at
**~200 µm**" from Fig. 4 (§4) — so the transition width is itself uncertain; the Fig. 4
reading is the one used, being the more specific. ⚠ Shen's **6.91** is a *different
constant* from Walsh's 6, which is why a range is carried and not a single factor.

**Consequence, computed in `(a0)`:** a 200–250 µm size population is **0.1–0.25 mm
thick**, not 35–50 µm. This widens the implied modulus band upward (m ≈ 15–30) rather
than inverting any conclusion.

## The four primary sources

### 1. Höser 2016 — ETH Zürich PhD, *Flame-jet assisted drilling technology*
- **Rock:** Central Aare granite — the rock this model represents.
- **Extracted:** **average spall size 200 µm**.
- **Pointer:** the value Kant et al. 2017 (SGW) adopt as the surface-to-failure-plane
  distance, **Fig. 4 caption: Δx = 200 µm [22]**, with Θ_SP = 550 °C.
- **Status:** ETH thesis; not fetched. Recorded via the Kant 2017 SGW citation, whose PDF
  *was* obtained 2026-10-04.
- **As a THICKNESS:** 200 µm size at 1–2:1 ⇒ **100–200 µm**.

### 2. Rossi, Jamali, Saar, von Rohr 2020 — *J. Pet. Sci. Eng.* **190**:107005
- *Field test of a Combined Thermo-Mechanical Drilling technology. Mode I: thermal
  spallation drilling.*
- **Extracted:** spall granulometry shows **two populations, ≈250 µm and ≈1 mm**.
- **Pointer:** granulometry section; the same two populations are reported by the 2025
  DeepU laser-spallation paper for laser spalls.
- **Status:** paywalled, not fetched.
- **As THICKNESSES:** ≈250 µm ⇒ **98–148 µm**; ≈1 mm ⇒ **~167 µm**. ⚠⚠ **The two
  populations differ 4.0× in size but only ~1.36× in thickness**, because the larger
  spall carries the larger aspect ratio — so they are **not two independent thickness
  targets**, which is what weakens `(a4)` as a discriminator.

### 3. Shen, et al. 2018 — *J. Pet. Sci. Eng.* **171**:1390
- *Aspect ratio of spalls of granite in flame-jet spallation.*
- **Extracted:** **280 spalls measured**; size-dependent aspect ratio, **6.91
  recommended**; Weibull **homogeneity factor m = 10–17** recommended.
- **Status:** paywalled, not fetched.
- **Bearing:** independently brackets the modulus this law consumes. `input_drilling`
  uses **m = 20**, just above Shen's recommended band.

### 4. Wilkinson & Tester 1993 — *Rock Mech. Rock Eng.* **26**:29–62
- **Extracted:** Weibull **m = 20 mechanical** vs **8 (Barre) / 13 (Westerly)** from
  spallation; **representative spall 0.4 mm** at 6×10⁻⁴ m/s; `T_s(q)` relation.
- **Status:** paywalled, not fetched.
- **Bearing:** the source of the "0.1–1 mm flake band" that `spall.flake_h_max`'s own
  comment cites (`MMWSpalling.H:9821`). ⚠⚠ **Its band tops out at 1 mm, while the shipped
  cap is 1.5 mm** — the cap exceeds the band it claims to come from.

## Corroborating sources (not part of (a6)'s four)

| source | extracted | status |
|---|---|---|
| **Meier 2017** (in repo) | PSD **average < 100 µm, p90 < 300 µm**; fragments, a **lower bound** | `validation/meier/meier-2017-pilot-660kg-validation-reference.md:231-235` |
| **Kant & von Rohr 2016** | **80 % of spalls < 1 mm** (Central Aare, Ø85×150 mm cores) | not fetched |
| **Rauenzahn & Tester 1989**, IJRMMS **26**:381 | the Weibull thickness distribution; spalls "less than a grain diameter thick, several grains wide" | thesis OCR **on disk** |
| **Walsh, Lomov, Wideman, Potter 2014**, IJRMMS **70**:375 | measured thickness vs size (Sierra White); the aspect-ratio curve of `(a0)` | **on disk** |
| **Walsh/Lomov/Roberts LLNL 2012**, OSTI **1035600** | grain-scale simulation: spall thickness **< one grain**; model **over-predicts size** vs collected spalls | PDF **obtained** 2026-10-04 |
| **Kant et al. 2017**, JGR Solid Earth **122**:1805 | DOI **10.1002/2016JB013800**; LEFM theory, Spallability number | paywalled at Wiley |
| **Augustine 2009**, MIT thesis | "spalls **0.1–2 mm thick**, diameters 10–20×" | MIT DSpace blocks non-browser download |
| **Lyu, Song, Li 2019**, Geothermics **77**:99 | extracts spallation temperature **and m** from drilling data | — |

## What this note does and does not back

✅ **Backs:** `(a0)`'s conversion; the modulus band `m = 10–20` that `(a5)` tabulates;
the observation that `spall.flake_h_max = 1.5 mm` **exceeds** the 0.1–1 mm band its own
source comment cites.

⚠ **Does NOT back:** any claim that requires reading a figure I did not open. In
particular the exact shape of Walsh's Fig. 4 between 200 µm and 1 mm is an
**interpolation chosen here**, not a digitised curve, and `(a0)` carries a range for that
reason. ⚠ `w = 2.35` is **not** from this literature — it is the median of the
cone-weighted failure density under Rauenzahn's normalisation, derived in
`Claude_markdowns/2026-10-04.md` §3 (Rauenzahn Eq. 4.96–4.100, Walsh 2014 Eq. 3) and
stated in the source comment at the `flake_fraction` site.
