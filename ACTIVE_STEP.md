# Active Step: P2-2b — the measured Grimsel properties and the relief factor (rev 2)

Status: ready for implementation — **rev 1 was BLOCKED at `(a1)`, correctly, and the fix is the user's: stop whitelisting hard-coded `rho`/`c_p` and convert every site to one style.** ord 0 and the rest of ord 1 PASSED and are carried. Then ord 2-5.

## Context

**2.27x too much rock = rock too CHEAP 1.33-1.54x x power OVER-DELIVERED 1.48-1.70x.** This step
takes the *cheap* half: connect the measured, temperature-dependent Grimsel laws to the places that
actually consume them, fit one relief factor on Kant with the flaw held, and test on four Meier legs
with nothing fitted.

**⚠⚠ MUCH OF THE MACHINERY EXISTS AND IS DEFAULT-OFF.** `UpdatePropertyFields`
(`MMWSpalling.H:5374-5440`) already writes `kappa_mf`, `beta_mf`, `E_phase_mf`, `mu_phase_mf`,
`rho_mf`, `Cp_mf` **per cell** from `Numeric::Material::Grimsel`, gated on
`microstructure.phaseN.law = grimsel` (default `constant` ⇒ byte-identical; `:310-328`, and its abort
is local to that key, inheriting nothing from `material.type` at `:134`). **What is missing is the
consumption:** `lambda` reaches conduction; `E`/`beta` do **not** reach the criterion (it reads
per-phase **scalars**, `Removal.H:398-405`); and `rho`/`c_p` are **plotted diagnostics only** —
conduction uses `phase_materials[p].effective_rhocp(T)`, and `phase_materials` is
`std::vector<Material::Constant>` built from the **scalars** (`:8074-8090`), which also builds the
enthalpy table.

**⚠⚠ `(a1)` WAS A PLANNER DEFECT TWICE OVER, AND THE FIX IS A RULE, NOT A LONGER LIST.** rev 1
enumerated the hard-coded `rho*c_p` sites by grepping for the **product**, found five, and then
barred the row at *"the five above and no sixth"*. Two sites read `rho` and `Cp` as **separate
scalars** and cannot match a product pattern — so the bar could only be satisfied by repeating the
planner's own search error. ⇒ **the whitelist is abolished.**

> **THE RULE (user's decision, 2026-10-08): no site in the rock path reads `rho` or `Cp` as a
> per-phase constant. Zero exceptions, so there is no list to keep correct.**
> Each site resolves one of exactly two ways — ⚠ the *mechanism* differs because the *context*
> does, and "consistent" does not mean "identical code":
> **(i) at a cell** ⇒ read the per-cell field at that cell;
> **(ii) building a per-phase object** (the enthalpy table) ⇒ build from the **law**, never the scalar.
> **THE STYLE: read `rho` and `c_p` SEPARATELY at the point of use; never store a pre-multiplied
> product.** A stored product cannot carry temperature dependence without being recomputed anyway,
> and the separate form is what the law provides and what the audit searches for.

⚠ The speed argument for keeping any of them fixed is worth nothing: **none is in the hot loop** —
`:1403` is owner-only per column, `:2027` a per-column branch, `:3174` per removal event,
`:8080`/`:8229` once at init. The crack integral's 129-calls-per-evaluation problem touches none of
them.

**The seven sites, derived against live source by `(a1)`'s regex rather than quoted:**

| # | site | what it feeds | resolves as |
|---|---|---|---|
| 1 | `Removal.H:1403` | pin reference → `pin_rhoCp` → `:3820`'s idle clock | (i), and the stored **product** retires |
| 2 | `MMWSpalling.H:3820` | the idle clock, consuming 1 | (i) |
| 3 | `Removal.H:3174` | flake energy estimate | (i) |
| 4 | `MMWSpalling.H:8080-8081` | `phase_materials` ⇒ the enthalpy table | **(ii)** |
| 5 | `MMWSpalling.H:8229-8230` | `InitializeMicrostructure`'s init fill of `rho_mf`/`Cp_mf` | ⚠⚠ **a second writer, but DO NOT simply delete it.** `UpdatePropertyFields` returns at `:5380` unless some phase sets `law = grimsel`, and skips non-Grimsel phases at `:5412` ⇒ under the shipped default the fields would be **never written (zero)** and the idle clock **divides by that**. ⇒ **ONE writer that covers EVERY phase under EVERY law** — at init under `constant`, per step under `grimsel` |
| 6 | `Removal.H:2027-2031` | the **vaporisation budget**, `L_tot = Cp_p*(T_vap_lo-T0_v)+L_m+L_v` | (i). ⚠ **inert here** (`vapor.enabled = 0` on both decks) ⇒ its witness is `(b6)`, not a Meier leg |
| — | `MMWSpalling.H:5041` | `St * wall_rho * wall_cp * v` | ⚠ **gas-side**, out of scope, named so it is not mistaken for rock |

⚠ Site 6 is the worst place for a constant `c_p`: it integrates one over the **entire span up to
vaporisation**, across which Grimsel's `c_p` rises ~1.5x.

**⚠⚠ THE SIGN OF THE PROPERTY CHANGE REVERSES, AND THE COLD END GETS CHEAPER.** Measured at ord 1:
Grimsel's `rho*c_p` is **0.9250** of the deck's at 293 K (`rho` 2650 < 2750 **and** `c_p(20 degC)`
758.3 < 790), crosses **1.0 at 325.7 K**, and reaches **1.4899 at 1200 K**. ⇒ **"~1.3x" is never a
scaling.** It is the integral over the span the removed rock actually spans, which is exactly why
`(d3)` reports the two factors apart. ⚠ A leg whose removed rock is cold gets **cheaper**.

**⚠⚠ THE FIT TARGET — THE USER'S DECISION, 2026-10-08, AND IT REPLACES A CONTRADICTION.** rev 2 said
`R` is fitted so Kant reproduces its **measured** onset **and** that Meier's surface temperature is
*unchanged within 5 K*. Both cannot hold: today's model fires Kant at **756 K** and Meier at a
**median 757 K**, against Kant's measured **823-873 K**. Fitting to the measurement therefore moves
the firing temperature up on **both** decks. "Pinned by construction" was right; "unchanged" was not —
it is pinned to a **new** value.

> **OPTION A IS TAKEN: `R` is fitted so Kant fires at 850 K.** The alternative (fit to today's 757 K)
> leaves Kant's measurement unused and makes `R` a pure renormalisation of the property change,
> defined by *"do not change the answer"* — circular, and unfalsifiable by construction.

⚠ **850 K is a deliberate UNDER-FIT, and its bias has a sign.** Kant's **823** (plane) and **873**
(surface) are **two different measurements, not an error bar**, and the criterion fires on a
reconstructed **surface** temperature, so 873 is the dimensionally matching figure. Fitting to the
midpoint therefore yields `R` slightly **too large** (too little relief) ⇒ the model stays slightly
**too easy** to spall. ⚠ The virtue is that **a result that closes at 850 K closes harder at 873 K**,
so a positive finding is robust to the choice rather than dependent on it.
⚠ The source also carries a **third** figure — **min 773 +/- 50 K** — whose upper edge *is* the 823 K
plane value, so the plane figure already sits near the low end of what Kant observed. Recorded here
so it is not later "discovered" and read as an omission.

⇒ **The firing temperature MOVES, by ~93 K, and that is the step's mechanism, not a side effect.**
Hotter rock leaves, hotter rock cost more to make, and the price rises. `(d1)` scores the **new**
value; `(d3)` carries a **temperature** factor beside the property one.

**⚠⚠ `K_Ic(T)` AND `R` ARE EXACTLY DEGENERATE ON KANT.** Both legs fire at **one** temperature, so a
toughness scale and a relief factor are indistinguishable there. ⇒ `(c2)` is scoreable **only** with
`kic.law = nasseri_table` and `kic.scale = 1.0` **pinned and stated in the bar**.

**⚠⚠ CONDUCTIVITY HALVING IS A MESH QUESTION THAT MAY INVALIDATE EVERY 2 mm RATE HERE.**
`lambda(600 degC)/lambda(20 degC) = 0.501`, so the flux skin `ell = k(T_s-T_0)/q` thins from
**~1.4 mm to ~0.70 mm** — **0.35 of a 2 mm cell**, and the surface reconstruction was built for
skins no thinner than half a cell. ⚠ The 1 mm ladder is **parked** (standing rule 12), so `(d8)`
records the sub-cell skin and marks every 2 mm rate **mesh-limited and provisional**. Pre-registered:
the 2 mm vs 1 mm rate gap **GROWS**. ⚠ Un-parking is the user's call and is **not** assumed.

**⚠ THE LATENT-HEAT HAZARD IS REAL IN CODE AND INERT ON THESE DECKS.** `Table::Build(const Material&)`
is polymorphic (`Table.H:30`) so it accepts a Grimsel object — but Grimsel parses only `T_ref`, so
`T_m`/`fus_width`/`L_m`/`T_vap_lo`/`T_vap_hi`/`L_v` sit at base defaults and the latent heat would
**silently vanish**. ⚠⚠ rev 1 called this *"the dangerous one"* and assigned it to `(b3)` — but
**neither deck sets a single phase-change scalar**, `effective_rhocp`'s boxcars need `fus_width > 0`
or `T_vap_hi > T_vap_lo`, so they contribute **exactly zero** on both forms and `(b3)` cannot see it.
⇒ the copy is done defensively and witnessed by **`(b6)`, a unit case that switches melting on**.
A gate aimed at an inert hazard passes for the wrong reason.

**⚠ WHAT `(b3)` IS REALLY ABOUT:** `(a2)` measured the two enthalpy curves diverging at equal
enthalpy by up to **-204.70 K**. That is the per-cell drift if conduction integrates one curve and
the `H->T` inversion reads the other. The ledger is the only row that can see it, and this class has
been caught twice by nothing else.

## Sources

`tests/MMWSpalling/validation/grimsel_properties/grimsel_granite_properties_reference.md` — 26 of 26
laws located by line at P2-2 ord 1. `theta` in degC. Consumed in `src/Numeric/Material/Grimsel.H`.

```
c_p(theta)    = 1370 - 178000/(theta+271)       [43]   Grimsel.H:122  cp.a/cp.b/cp.th0
lambda(theta) = 2000/(theta+563)                [89]   Grimsel.H:135  lambda.a/lambda.th0
kappa         = lambda/(rho c_p), NEVER set apart       [112]
rho           = 2650 +/- 50                     [141]  Grimsel.H:120  rho
E_0 = 53.3 GPa; E/E_0 = 1.0/0.85/0.65/0.50 at 200/400/500/600 degC   [194]  Grimsel.H:156
nu  = 0.30                                      [180,195]  Grimsel.H:106
alpha(theta) = 8.4/9.5/16.4/23.4/26.9/28.3 e-6 at 20/300/400/500/550/600 degC  [Table 3] ** [ESTIMATE] **
K_Ic,0 = 1.5 +/- 0.3; AT-TEMPERATURE 1.0/0.90/0.82/0.78/0.75 at 100-500 degC   [209,223]
UCS = 169 +/- 37 MPa                            [180]
Kant MEASURED onset: surface ~873 K, plane ~823 K, min 773 +/- 50 K, q 0.45-0.87 MW/m2  [251-253]
```

⚠⚠ **`alpha(theta)` IS THE ONE `[ESTIMATE]` AND IT IS LEAST CERTAIN WHERE SPALLING HAPPENS.** The ~3x
rise by 600 degC is the **quartz alpha-beta transition at 573 degC**, which the Meier and Kant firing
temperatures **straddle**; the reference gives **+/-30 %** above 300 degC, part of it opening
microcracks rather than loading the flake [172]. ⇒ `(c7)` fits `R` **twice** — `alpha` as given and
capped at 300 degC via `grimsel.alpha_cap_theta` — and **if `R` leaves the band between the two fits,
`(c2)` has no power and that is the finding.**
⚠ **`nu` = 0.300 and does not move** (source above; the deck's `E/(2*mu) - 1` reads 0.300
identically). For scale, 0.37 against 0.25 is a **19 %** difference in stress.
⚠ `R` is `spallation.sp.stress_relief`, default **1.0**. Carried: `sigma` **173.8 / 300.4 / 320.8
MPa** each by its own method; `R = UCS/sigma` **0.5625 / 0.5268**; implied flaw **~9 um**.

## Goal — deliverables in this order

1. **CARRIED: ord 0 and ord 1 are done.** `(a0)` PASS on the shipping binary (provenance
   self-asserting), `(a2)` and `(a3)` PASS. ⚠ Re-run `(a1)`'s regex **only**, against the rule below.
2. ⚠⚠ **CONVERT ALL SEVEN SITES TO ONE STYLE** — the rule in the context. Sites 1, 2, 3, 6 read the
   per-cell fields **separately**; site 4 builds from the **law**; site 5 is **deleted** in favour of
   one `UpdatePropertyFields` call at init. ⚠ The stored **products** `rev_rhoCp`/`pin_rhoCp` retire:
   carry `rho` and `c_p` and form the product at use. ⚠ This touches the pin/idle-clock path, which
   has its own history — `(b7)` is its witness.
3. **The tree gate ONCE, BEFORE the legs**, with `law = constant` everywhere — the shipped default.
   ⚠ `beam_void_closure` stays **rc = 1 by standing decision** (`ARCHIVE_DONE.md:16499`); carried,
   not caused, and **not** to be re-litigated.
4. **Strain and stiffness as dense uniform-grid tables**, with the crack integral evaluating both at
   the local temperature from the criterion's own reconstructed profile. ⚠⚠ **Strain is never
   interpolated in space.** `beta_mf` **stays** for the eigenstrain.
5. **The per-cell cache for conduction and the eigenstrain only**, registered **unplotted**
   (`RegisterNewFab(..., false)`) — the `true` on `kappa_mf`/`beta_mf`/`rho_mf`/`Cp_mf` (`:387-430`)
   is what puts them in the 34-component hashed cell plotfiles.
6. **Route `rho` and `c_p` into conduction**: build `phase_tables[p]` from the Grimsel object when
   `grimsel_law` is set, **after copying `T_m`/`fus_width`/`L_m`/`T_vap_lo`/`T_vap_hi`/`L_v` from the
   phase**, so the latent heat survives. ⚠ Conduction and the `H->T` inversion read the **same**
   curve.
7. **Extend `unit/lateral_flux_lock`** for `(b4)`: one strain source, one stiffness source, the
   tables agree with `Grimsel::sigma_thermal` to round-off, the cache unplotted, and ⚠⚠ **zero
   per-phase scalar reads of `rho` or `Cp` in the rock path, and zero stored pre-multiplied
   products** — the rule asserted as a predicate, not a count.
8. ⚠ **NEW unit case `(b6)`: switch melting ON** (`T_m`, `fus_width`, `L_m` set) and assert the
   rebuilt table carries the latent heat under both forms. This is the only witness the hazard has.
9. **KANT, short**, deck `sp_kant_onset_meiercfg`: `K0` baseline · `K1` + thermal set · `K2` + the
   mechanical set with the at-temperature `K_Ic`. ⚠ **Each leg states which criterion it runs**, so
   conduction separates cleanly from mechanics.
10. **FIT `R` ON `K2`, FLAW HELD** — twice, `alpha` as given and capped at 300 degC.
11. **FOUR MEIER LEGS LAST, NOTHING FITTED.** `CTRL6a` thermal only · `CTRL6b` full set · `CTRL7` =
    **CTRL5 minus the prescribed block plus the feet block, KEY BY KEY** from
    `input_feet_d2c:166-173` · `CTRL8` = CTRL7 + full set. ⚠⚠ The modes **refuse each other's keys**.
12. Price **per path**, five depth bands, per-15-degree sectors, each leg's **argv** beside its hash.
13. The archive line, and the P2-2d carry-forward.

## Criteria — run in order; **only a `gate` halts**

| # | ord | test | bar | kind | what a failure of this row means |
|---|---|---|---|---|---|
| (a0) | 0 | P2-2a `(c1)` on the shipping binary | **CARRIED: PASS** — R150/R500 dT **463.31 / 465.80**, \|diff\| 0.01 / 0.00, each leg stamping the shipping binary | — | — |
| (a1) | 1 | ⚠⚠ **re-barred as a RULE.** No per-phase scalar read of `rho` or `Cp` in the rock path, and no stored pre-multiplied product | **zero** matches for `phases[..].rho`/`.Cp`, `d_rho[`/`d_Cp[`, `m_rho`/`m_cp`, `rev_rhoCp`/`pin_rhoCp` outside the parse block and the law's own definition | **gate** | ⚠⚠ rev 1 barred this as *"five and no sixth"* and the bar could only pass by repeating the planner's product-grep error. A **count** is brittle; the predicate is not. A surviving scalar keeps part of the model at a fixed `rho*c_p` while the rest varies, and the books drift without saying where |
| (a2) | 1 | the two enthalpy curves and their inverse | **CARRIED: PASS** — `rho*c_p` ratio **0.9250** at 293 K, crossing 1.0 at **325.7 K**, **1.4899** at 1200 K; enthalpy ratio **1.298** at 1023 K; `T(H)` diverging to **-204.70 K** | — | — |
| (a3) | 1 | `(d2)`'s arithmetic, written out | **CARRIED: PASS** — form, both factors, the x4 quarter-domain factor and its cancellation | — | — |
| (b1) | 2 | tree identity with `law = constant` | **0 differing of 2639** | **gate** | the default is not byte-identical, so no property leg is attributable |
| (b2) | 2 | the tables vs `Grimsel::sigma_thermal` | agree to **<= 1e-12** relative over a 708-point sweep | **gate** | the `(b6)`-of-P2-2 pattern, which already caught a factor-of-two reading of the strain |
| (b3) | 2 | ⚠⚠ **conduction and the `H->T` inversion read ONE curve** | `ledger_err_solid` **<= 1e-6** (dimensionless) on `K1` | **gate** | `(a2)` measured the two curves diverging to **-204.70 K** at equal enthalpy; that is the per-cell drift if they disagree. Caught twice before by the ledger and by nothing else |
| (b4) | 2 | the lock, incl. `(a1)`'s predicate asserted in-binary | RED pre-change, green after | **gate** | a lock that passes on both binaries locks nothing |
| (b5) | 2 | the cost of the two hot paths | the crack integral's **and the `H->T` inversion's** share of runtime, **measured** before and after | report | the cache is *expected* cost-neutral; this is the number that confirms it rather than assuming it |
| (b6) | 2 | ⚠ **NEW. Melting switched ON** | the rebuilt table carries the latent heat under **both** forms, to **<= 1e-9** relative | **gate** | ⚠⚠ the latent-heat loss is **inert on both production decks** (`vapor.enabled = 0`, no phase-change scalar set), so `(b3)` cannot see it. Without this row the hazard ships unobserved |
| (b7) | 2 | ⚠ **NEW. The pin/idle-clock path** after the stored products retire | `t_cell` and the pin reference unchanged to **<= 1e-9** relative under `law = constant` | **gate** | that path has its own history (D2h/D2i's idle clock); retiring a stored product there must be identity-preserving on the shipped default before it is allowed to vary |
| (c1) | 3 | `K0`-`K2` onset level | reported per rung, each naming its criterion | report | which property moved the level |
| (c2) | 3 | `R` fitted on `K2`, flaw **held** | `R` in **0.41-0.69**, ⚠ with `kic.law = nasseri_table` and `kic.scale = 1.0` **PINNED and stated** | **gate** | outside the band the fit absorbs a property error. ⚠⚠ Unpinned, `K_Ic(T)` and `R` are **exactly degenerate** on Kant's single firing temperature and the band tests nothing |
| (c3) | 3 | the flaw **implied** by the fitted `R` | **5-20 um** | report | ⚠⚠ **this row cannot see the fit.** `R*sigma = UCS` ⇒ the implied flaw is `pi*(K_Ic/(3.526*UCS))^2` — no `sigma`, no `R`. A consistency check on **`K_Ic` and `UCS` only** |
| (c7) | 3 | the `alpha` sensitivity | `R` under **both** fits, with the span between them | **verdict** | ⚠⚠ if `R` leaves the band **between** the two fits, `(c2)` has no power. The quartz transition at 573 degC sits between the Meier and Kant firing temperatures |
| (d1) | 4 | `T_s` across the Meier floor, `CTRL6b` | ⚠ **re-barred.** Equals the **fitted Kant onset, 850 K, within 5 K** — NOT today's 757 K | report | ⚠⚠ rev 2 barred this at *unchanged*, contradicting its own fit target. `T_s` is pinned **to the fitted value**; a `T_s` still near 757 K means `R` did **not** reach the criterion, and a `T_s` far above 850 K means something other than the fit moved it |
| (d2) | 4 | **ENERGY LEAVING AS ROCK**, ratio to its own baseline. ⚠ **GATE on the FEET pair only; REPORT on the prescribed pair** | **0.95-1.10** on `CTRL8/CTRL7` — ⚠ equivalently, with the price up ~x1.47, a **volume fall of x0.64-0.74**; both factors reported beside it | **gate** | ⚠ a product is **not** immune to compensating errors in its factors — rev 5's "feedback-immune" was overstated. Stating it as a volume prediction first is what stops a volume change being read as confirmation of the price |
| (d3) | 4 | the price, **decomposed into property x temperature** | the **property** factor at the baseline's own temperatures (enthalpy ratio **1.19-1.23** at 757 K surface / 677 K flake mean) · the **temperature** factor from firing at 850 K instead of 757 K · and the realised price against **~1.37 GJ/m3** | report | ⚠⚠ rev 2's anchor of **1.298** was evaluated at **1023 K**, hard-coded as "the Meier firing temperature" with no source — a temperature the baseline never reaches. ⚠⚠ And rev 1's *"mean temperature of removed rock, x1.1-1.2"* was **never a factor this step produces**: the realised price already prices rock at the flake mean, and the flake-to-surface ratio is a **fixed 0.824** on a self-similar profile. It was a correction for comparing against literature prices quoted at the surface, and a prior revision mistook it for a physical effect ⇒ it is reported as a **diagnostic pre-registered at 0.82**, not as a factor |
| (d4) | 4 | volume, ratio to its own baseline per mode | reported, **not banded** | report | prescribed is damped by the standoff, feet is not |
| (d5) | 4 | `CTRL7/CTRL5` — descent mode alone | reported | report | ⚠⚠ the **standoff masking**: `feet` holds the nozzle ~75 mm off the axis against 179-224, map factor **0.45 -> 1.0**; D2ac's **1.98-2.23x** should appear |
| (d6) | 4 | `foot_stall_time`, and feet pads on a laterally-cut rim | ⚠ **numeric now:** completes to 490 s **and** cumulative `foot_stall_time` **< 10 %** of elapsed **and** mean ROP over the final 100 s **>= 1.0 m/h** | **gate** | ⚠⚠ **never run.** Feet jammed in the D2o era (ROP 1.49 -> 0.11). If it stalls, `CTRL8` is **not** launched |
| (d7) | 4 | bore, five bands, per sector; and the flake thickness | **+/- 3 mm**; thickness **0.2-0.3 mm** | report | the thickness is a **prediction** [289-291] — the one out-of-sample test |
| (d8) | 4 | ⚠⚠ **the skin against the mesh** | the firing-surface `ell` reported in **cells**; predicted **~0.35** at dz = 2 mm | **verdict** | below half a cell the surface reconstruction is outside its design range ⇒ **every 2 mm rate here is mesh-limited and provisional.** Pre-registered: the 2 mm vs 1 mm gap **GROWS**. ⚠ 1 mm is parked (rule 12) |
| (z) | 5 | the archive line and the P2-2d carry-forward | written | verdict | — |

## Pre-registered predictions

- ⚠⚠ **P1 is RESTORED, correctly derived this time.** rev 2 withdrew it entirely, which left `(d3)`
  a pure report — **and a report with no prediction cannot fail.** The predicted realised price is
  **~1.37 GJ/m3**, interpolated from the fit-target table (823 K -> 1.29, 873 K -> 1.43) at the
  chosen 850 K. That is **x1.47** on today's 0.930, which sits at the **top of the independently
  derived 1.33-1.54 "rock too cheap" band** — a genuine consistency check, since that band was
  measured from the volume/energy split and not from the laws.
  ⚠⚠ **Meier's price is 1.40, so the predicted residual is ~2 %. That is suspiciously good and is
  flagged, not celebrated.** What keeps it honest: Kant's onset is a published number we do not get
  to choose, `(c2)`'s band on `R` comes from the rock's own strength, and `(c7)` reports `R` under
  two expansion curves. If the residual lands near 2 % **and** `R` sits inside its band **and**
  `(c7)`'s two fits agree, the result stands; if it lands near 2 % with `R` outside the band, the fit
  absorbed an error and the agreement is an accident.
  ⚠ rev 1's arithmetic is withdrawn for the record: 0.930 -> ~1.25 needed **x1.344** while its own
  factors gave **x1.430-1.560**, unreachable from itself.
- ⚠⚠ **P2 is INVERTED.** rev 2 predicted `T_s` does not move. It **does** — by **~93 K**, from a
  median 757 K to the fitted 850 K — and that rise **is the step's mechanism**, not a side effect:
  hotter rock leaves, hotter rock cost more to make, and the price rises. ⚠ A `T_s` still near 757 K
  means `R` never reached the criterion.
- **P3 — `R` is ONE material number near 0.55**, implied flaw near 9 um.
- **P4 — the 2 mm vs 1 mm rate gap GROWS**, because `lambda(600)/lambda(20) = 0.501` thins the skin
  to ~0.35 cells at 2 mm. ⚠ Unfalsifiable here while 1 mm is parked; `(d8)` records it.
- **P5 — the feet legs look WORSE against Meier, and that is a finding.** `(d5)` makes two weeks of
  masked delivery error a number.
- **P6 — the cold end gets CHEAPER.** `rho*c_p` is **0.9250** of today's at 293 K and crosses 1.0 at
  **325.7 K**, so any leg whose removed rock is cold moves the other way. Anyone reading "measured
  properties make rock dearer" has the cold end inverted.
- **P7 — `(c7)` may dissolve `(c2)`.** If `R` under the two `alpha` fits brackets or exceeds
  0.41-0.69, the band was never a test. Writing that down first is the point.
- **P8 — `(b6)` and site 6 are inert on every leg this step runs**, by construction
  (`vapor.enabled = 0`, no phase-change scalar set). They are converted and witnessed **because the
  code is wrong**, not because this step can see it. ⚠ A prediction of inertness is normally a design
  smell; it is accepted here only because the witness is a dedicated unit case, not a production leg.

## What would make this inert — witness each

- **The law parsed but not reaching conduction** — `(a2)` offline, `(b3)` on a run.
- **A surviving per-phase scalar** — `(a1)`'s predicate, asserted in-binary by `(b4)`.
- **Latent heat silently lost** — `(b6)`; `(b3)` **cannot** see it on these decks.
- **The tables disagreeing with the reference implementation** — `(b2)`.
- **The cache registered plotted**, re-pinning the tree — `(b4)`.
- **The pin/idle-clock path moving** when a stored product retires — `(b7)`.
- **`R` absorbing a toughness error** — `(c2)` with `kic.*` pinned, and `(c7)`'s two fits.
- **A 2 mm result read as physical** when the skin is sub-cell — `(d8)`.
- **`feet` jamming, or jamming on a laterally-cut rim** — `(d6)`.

## Guardrails

- ⚠⚠ **NOTHING IS FITTED ON A MEIER LEG.** Kant calibrates, Meier tests. `R` is frozen before
  deliverable 11 and its value recorded.
- ⚠⚠ **THE SURFACE TEMPERATURE RISING ~93 K IS THE STEP'S MECHANISM, NOT A SIDE EFFECT** — rev 2's
  guardrail said the opposite ("a rate change explained by a hotter surface is wrong by
  construction") and that was exactly backwards once `R` is fitted to a measurement above today's
  onset. ⚠ What remains forbidden is a `T_s` that moves for any reason **other** than the fit: it is
  pinned to the fitted 850 K, not free. See `(d1)`.
- ⚠⚠ **STRAIN IS NEVER INTERPOLATED IN SPACE.** Tables on a uniform temperature grid, evaluated at
  the local temperature from the criterion's own reconstructed profile.
- ⚠⚠ **`beta_mf` STAYS** for the eigenstrain; the elastic solve is out of scope.
- ⚠ **`nu` does not move** (0.300, sourced above). ⚠ `kappa` is never set apart from
  `lambda/(rho c_p)` [112].
- ⚠ **A charged-volume total includes BOTH paths** (floor `h_applied` **and** wall
  `lat_bank_hf x 6.0226`); omitting the wall produced two withdrawn planner figures.
- ⚠ **Never sum a per-step `thermo` column**; `ledger_err_solid` is **dimensionless**.
- ⚠ **Score the bore from the void MASK**, five bands; take every gas-state and geometry number
  **from the leg under test**.
- ⚠ **Set `amr.thermo.plot_int` EXPLICITLY on every leg** — it defaults to **-1 = never**, so a leg
  can exit rc 0 with no `thermo.dat` and `(b3)`/`(d2)` would be unscoreable.
- ⚠ **Archive an output directory BEFORE launching into it**, never after; record each leg's **argv**
  beside its hash (rule 8); every long leg under `caffeinate -dis`. ⚠ **`pkill -f` is forbidden** —
  use `TaskStop`.
- ⚠ **A WAIT LOOP THAT GREPS FOR A PROCESS NAME MATCHES ITSELF.** `pgrep -f guarded_post_all` matches
  the waiter's own command line, so four gate-waiters polled forever in rev 1. Match on a PID or a
  sentinel file, never on a name the waiter also carries.
- ⚠ **When a harness can skip, assert it executed** — `run_sweep.py` returns `skip (done)` without
  `--force`, which turned a P2-2a measurement into a stale read.
- ⚠ **The manifest is MD5**, not SHA-256. ⚠ **An enumeration is derived against live source, never
  quoted** — that is how `(a1)` went stale.
- **Do not change** `side_face_*`, the wall law, the floor `h_expr`, or `V0`'s mesh-independence role.
- ⚠ **The flux exponent is P2-2d** — `size_volume` stays at `cell`. Do not sweep `m` here.

## Commands

```bash
# ord 1 -- OFFLINE, no build. (a0)/(a2)/(a3) are CARRIED; re-run (a1)'s predicate only.
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/studies/p2_2b_props/derive_enthalpy.py \
  > tests/MMWSpalling/studies/p2_2b_props/ENTHALPY.txt

# ord 2 -- the TREE GATE FIRST on the shipped default, then build, then the locks.
bash tests/MMWSpalling/studies/tree_sweep_0922/witness/guarded_post_all.sh p2_2b_pre
EIGEN=$PWD/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/lateral_flux_lock/test
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/latent_heat/test
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/pin_clock/test

# ord 3 -- KANT, short. Each leg names its criterion.
bash tests/MMWSpalling/studies/p2_2b_props/run_kant_ladder.sh   # K0 K1 K2
bash tests/MMWSpalling/studies/p2_2b_props/fit_relief.sh        # alpha as given AND capped

# ord 4 -- MEIER, last, nothing fitted. CTRL7 gates CTRL8 via (d6).
bash tests/MMWSpalling/studies/p2_2b_props/run_ctrl6a.sh
bash tests/MMWSpalling/studies/p2_2b_props/run_ctrl6b.sh
bash tests/MMWSpalling/studies/p2_2b_props/run_ctrl7.sh
bash tests/MMWSpalling/studies/p2_2b_props/run_ctrl8.sh

# ord 5 -- scoring.
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/studies/p2_2b_props/score_props.py \
  > tests/MMWSpalling/studies/p2_2b_props/SCORE_p2_2b.txt
```
