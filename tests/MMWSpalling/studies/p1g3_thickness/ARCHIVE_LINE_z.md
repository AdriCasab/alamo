# P1g-3 `(z)` — the archive line

Written as P2-0 deliverable 1(iii) / `(p3)`. Lift verbatim into
`docs/project/ARCHIVE_DONE.md`. Every number is from
`studies/p1g3_thickness/SCORE_p1g3.txt`, scored by `score_bore.py` on
`MEIER_LAT` vs `MEIER_CTRL2` — **one binary, `7908bd6b…`**, which is what makes
the pair comparable on a toolchain that is not byte-reproducible.

---

**P1g-3 — the flake thickness from the rock (2026-10-05).** The typed-in flake
thickness is **retired in favour of a derived one**: `h = w·ℓ_q/m`, where
`ℓ_q = κ(T_s−T₀)/q_net` is the flux-implied heated depth, `w = 2.35` is the
median of the cone-weighted failure density (Rauenzahn Eq. 4.96–4.100 / Walsh
2014 Eq. 3), and `m = weibull.m` is the Weibull modulus. Two keys,
`spall.flake_length = twocell|flux` and `spall.flake_fraction = cf|weibull`,
both **default off**; identity **0 of 2639**. `ℓ_q` is an **identity, not an
approximation** — `ℓ_eff·q_ratio ≡ κ(T_s−T₀)/q_real` — so it costs nothing to
compute. The law lands at **170 µm on the floor** and **742 µm on the wall**,
and it **retires `spall.flake_h_max` by making it unnecessary**: the 1.5 mm cap
governed **84.2 %** of wall firings in the control and **0.0 %** with the law
on, confirming P4. ⚠ The cap's own source comment cites a 0.1–1 mm band
(Wilkinson & Tester 1993), which **1.5 mm exceeded**; the key is still shipped
and still defaults to 1.5 mm.

**`m` remains a choice, and the choice is contradictory.** `h ∝ 1/m` exactly
(425.8 / 262.0 / 227.1 / 170.3 µm at m = 8 / 13 / 15 / 20), so the modulus sets
the thickness outright. Matching measured thicknesses wants **m ≳ 20**;
independent Weibull fits want **8–17** (Shen 2018 recommends 10–17 from 280
spalls; Wilkinson & Tester give 8 Barre / 13 Westerly from spallation and 20
mechanical). **m = 15 matches neither**, and the shipped Meier input reads
**20** while the fixtures read 15 — silently, until deliverable 4 moved the
`weibull.m` parse out of `if (weibull_enabled)`. Fitting `m` is deferred to
P1g-4.

**⚠⚠ THE RESULT IS A NEGATIVE, AND IT IS THE VALUABLE PART: the width does not
follow from the thickness.** Against Meier's measured Ø 84–96 mm, the control
runs **+10.7 mm** over the 96 mm top of that band and the fully-derived leg runs
**+10.2 mm** — **5 % of the excess closed**, mean band difference **−0.54 mm**
over five comparable depth bands, worst −2.06 mm. This is not a null from an
inert mechanism: the mechanism demonstrably bit (the cap went 84.2 % → 0.0 %,
wall thickness ×0.495, floor thickness ×0.227). It is a null because **removal
is energy-limited and the thickness only sets event granularity**. Both regimes
compensate almost exactly: vertical **thickness ×0.227 against firings ×5.03,
product 1.14**; lateral **×0.495 against ×1.96, product 0.97**. The beam's power
fixes the mass of rock removed; a thinner flake is simply fired proportionally
more often. Two independent regimes reaching the same product ≈ 1 is the
evidence, and it is stronger than P1g-2's earlier hint (~95 % of removal
surviving a 7× thinner flake) because here the thickness is **derived** rather
than imposed, at Meier scale, on the wall as well as the floor.

**What this hands the successor.** It **removes flake thickness as a candidate
explanation for the +10 mm bore**, which closes a line of enquiry rather than
merely failing to confirm one, and it points back at **D2ad's binding term — the
h-map's radial fall-off across 28–40 mm**. The 10–50 mm collar is untouched at
**+25.8 mm** and is a separate problem: it forms in the first seconds, before
either law bites. ⚠ The 300–400 mm band is excluded by the front rule (398 vs
368 mm deepest removal), and the 90–110 mm band is excluded by construction —
**never score the bore at 100 mm depth** (Meier p. 223 pause underream).

**⚠⚠ What P1g-3 did NOT deliver.** **`(b5)`, the MMW peak-column leg, was never
run, so the step has NO MMW answer** — the law's behaviour under a millimetre-
wave volumetric source, as opposed to a flame-jet surface flux, is **unmeasured**.
`mmw_peak.py` and `form_witness.py` were never written. Scaling the packet down
was not the implementer's call and is flagged as such: 14 deliverables and 18
criteria spanning offline analysis, two source sites, a new MMW leg and two full
Meier runs is plausibly two packets. One defect is carried open: the **`ell_q`
event-log column reads 0** on every event because `SharedSurfaceRecon` rebuilds
the struct from per-column arrays and a field with no backing array keeps its
default; adding `surf_ellq_col` did not fix it. The physics used the flux length
correctly (through `surf_hf_col`), but **the log does not carry the length the
thickness was built from**, so a successor cannot re-derive the bound from that
column. `(a6)` is partial: **3 of its 4 primary sources are paywalled and were
not fetched**, which `validation/spall_thickness/REFERENCES.md` states per
source rather than implying a read.
