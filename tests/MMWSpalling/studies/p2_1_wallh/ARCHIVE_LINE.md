# P2-1 `(z)` — the archive line

Deliverable 11. Lift verbatim into `docs/project/ARCHIVE_DONE.md`. Numbers from
`WALL_H.txt`, `LOCK_prechange.log`/`LOCK_postchange.log`, `TREE_p2_1.log`,
`FIXTURE_p2_1.log`, `SCORE_p2_1.txt`.

---

**P2-1 — the wall heat-transfer coefficient, derived (2026-10-06).** The hole wall
was heated through one hand-typed scalar, `surface_patch.side_face_h = 169.3`,
applied flat to every wall face at every depth, and that scalar carried 54 % of
all power delivered to rock. It is now a law, behind
`surface_patch.side_face_h_law = typed | derived` (default `typed`).

**The law, and why it is one law and not two rivals.** The two earlier
derivations are the same flow at two heights. At the corner the gas has just
turned off the floor: the layer is **4.4 mm** thick and runs at **264 m/s** on a
free stream. Once it fills the wall-to-shield gap there is no free stream left
to entrain and the flow **is** the duct. So
`h_wall(z) = h_wj(z)` for `z < z_fill`, `h_duct` above — ⚠⚠ a **hard switch**,
never `max(h_wj, h_duct)`: `h_wj` is still **95.7** at z = 400 mm against the
duct's **57.5**, so under the `max` the duct would never be selected and the
deep wall would keep a model with nothing left to entrain.

**ord 1 ran offline, before any source change, and all three of its gates
passed.** `(a1)` corner **317.1 W/m²K** against the gate's 250–350 and the three
literature routes' 300 ± 50. `(a2)` `z_fill` **155.1 mm** against 150–220, with
the switch reported as a factor **2.60** step. `(a4)` all thirteen reachability
rows resolve, scored **by reading the source inside the script** so the row
re-runs rather than being asserted in prose.

**⚠⚠ THE SIGN CAME OUT THE OTHER WAY, AND THAT IS THE STEP'S RESULT.** `(a3)`
re-derived the cutting band against the defined datum and the law's
face-weighted mean over it is **242.2 W/m²K**, i.e. **1.43× the typed
constant**, not below it. The packet's P2 expected "close to or above" and its
P4 wording assumed a *fall*; the mechanism is unchanged but its **direction is
reversed**. The reason is scope, not physics: the `r ≥ 40 mm` gate admits two
populations under one name — the **crater floor** between the gate and the bore,
546 cells sitting in the radial wall jet at full strength (**mean 338.6**,
2.00× typed), and the **hole wall** proper, 1176 cells (**mean 203.6**, 1.20×).
P2-0's "the constant is 3.4× too high" was measured against a **duct** law
applied where the flow is a **jet**, and it is the duct comparison that was
wrong, not the constant's size.

**⚠⚠ THE HEIGHT DATUM IS LOAD-BEARING AND THE PACKET'S ONE-SENTENCE DEFINITION
DOES NOT DISAMBIGUATE IT.** Three readings were measured, not asserted:

| datum | band p50 | band mean h | /typed |
|---|---|---|---|
| nearest open neighbour's top (**rejected**) | 7.0 mm | 293.7 | 1.735 |
| **inward-ray corner (chosen)** | **41.0 mm** | **242.2** | **1.430** |
| axis front (**rejected**) | 89.0 mm | 176.5 | 1.043 |

The chosen rule is *the deepest column top on the ray from this column inward to
the axis, truncated at `side_face_h_r_min`*. It reproduces the packet's own
stated 30–65 mm band; the axis-front reading reproduces the **p50 90 mm** prior
packets quote, which is the check that these are the same three readings the
packet is talking about. ⚠ The nearest-neighbour reading is a **discretisation
artefact**: it measures the one-cell staircase tread, and the wall-jet layer is
4.4 mm at the corner growing to ~15 mm, so it spans several treads — a 2 mm
ledge every 10 mm is roughness, not a fresh corner. Resetting `z` at each tread
would also make the datum mesh-dependent, which is the D2j-0b failure mode.

**The resolution, and the lock that now forbids the alternative.** All four
roles resolve `h` through **one** function, `SideFaceH`, from **one** datum array
built **once per step** — building it per role would have made the four roles
disagree by a sub-step, which is the thing deliverable 4 exists to prevent.
⚠ The criterion's resolution had to move **inside its `k` loop**: its
column-level block has no `k`, so leaving it there while the kernel went
per-height would have been the **fourth** instance of the family
`Removal.H:589-594` pre-registers in source. `(b3)`: the extended lock is
**RED on the pre-change source (5 failures) and green after** — and the old `h`
row could not have caught this, because it is a *value* claim about four scalars
and a per-height/per-column split passes it unchanged.

**`(b1)` 0 of 2639.** `(b2)` upper half **p95/p05 = 2.254** on a derived leg
against the bar's 2×. ⚠ That span was **vacuous on the first attempt**: the
histogram was filled over every lateral face, including `r < 40 mm` where the
typed path uses the `h_expr` map, which spans 167→1358 on its own — so the row
would have passed with the law **switched off**. It is now restricted to the
law's own scope. `(b2)` lower half is scored as byte-identity between the
shipped fixture leg and one naming `side_face_h_law = typed`, plus zero added
thermo columns; ⚠ the packet's literal "`side_h_used` == 169.3 to 6 s.f. on
`G1_FLAT_SHIFT`" is **not scoreable** — that leg's coefficient is **400.0** by
argv, and the `side_h_used_*` columns are registered only under `derived`,
which is precisely what keeps `(b1)` clean. Registering them on a typed leg
would have traded `(b1)` for `(b2)`.

**⚠⚠ THE SMOKE RUN CAUGHT TWO DEFECTS A GREP CANNOT.** Both were silent and both
would have corrupted the leg: the azimuth probe sampled
`surface_patch.radius = 0.2 m` against a **0.12 m** domain, every sample fell
outside, and the fallback returned a **full** circle — so `mdot_full` came out
**4× too small** and the whole law with it. And before the hole clears the
shield the bore collapsed onto the `R_SHIELD` clamp, giving a **0.3 mm** gap,
`h_duct` **770** and `z_fill` **0**, i.e. every face on a duct branch that does
not exist. Both now report rather than default: the probe returns 0 and the law
declares itself off through **`wall_law_on`**, which reads 0 for the first ~40 s
and 1 thereafter. ⇒ **a degenerate input is a failure to be reported, not a
default to be guessed**, and an inert period that nothing shows is the inert-key
failure three keys in this line of packets already hit.

**ord 3: every gate passed and the one verdict row did not.** `(c1)`
**7.7847e-04** against the 5e-3 bar; `(c1w)` **1.591e-04** against the
comparison leg's own **1.538e-04**. `(c3)` worst band-mean **+0.73 mm** and
worst 15° sector **+1.22 mm**, both inside ±2 mm and **positive**, which is the
sign `(a3)` called offline before the leg ran. `(c5)` **2.496 m/h on both
legs** — identical to 4 s.f., and the depth at 150/200/250/300 s differs by
≤2 mm, so the law did **not** reach the floor risers and `(c4)` is attributable.
`(c4)` **2.352×** against the band's 2.15–2.23: outside, but **+3.48 %** on the
comparison leg, inside the ~5 % the packet's own arithmetic bound permits, so it
contradicts neither `(c3)` nor `(c5)`. ⚠ The band presumed a *fall* in wall
energy; it reverses with `(a3)`'s sign, exactly as P4's wording does. The mask
estimator reads CTRL4_FIXED at **2.886 L / 2.273×**, agreeing with the packet's
own stated figure to 3 s.f. — which is what makes the +3.48 % the law's.

**`(c2)` NOT MET — 1054.8 W against 1383.9 predicted, −23.8 % — and the miss is
in the BAR, which is shown rather than asserted.** The packet offers two
readings for a `(c2)` miss and **both are excluded by measurement**:
`ledger_err_solid` moves only **9.87e-05 → 1.012e-04**, so the ledger still
mirrors the kernel; the **kernel/march split holds at 1.0796 → 1.0743** where a
march left on the typed scalar would have driven it to ~1.179; and
`lat_fires` ×**1.279**, `lat_bank` ×**1.283**, `lat_shift` ×**1.535** on a
candidate population that moved only ×**1.027** — the criterion fired *more on
the same faces*, which is the law reaching it. The shortfall decomposes exactly:

> **h ×1.430 (predicted) × ΔT ×0.835 (not in the bar) × area ×0.9190 (not in the
> bar) = 1.0977**, against the measured wall-power ratio **1.0977**.

⚠⚠ **A hotter-driven wall runs hotter**, so `(T_gas − T_s)` shrinks and the power
rises far less than the coefficient: the wall FLUX rose ×1.1944 where the
coefficient rose ×1.43. The ord-1 bar said "dT held" in as many words, and it
also held the **wetted wall area** that the law itself changes (×0.919). ⇒ **an
energy bar for a Robin coefficient must carry the surface-temperature feedback,
or it predicts the coefficient and not the energy.** That is the generalisable
lesson and it is the fourth time in this programme that energy and geometry have
proved more weakly coupled than a hand bar assumed.

**⇒ P2-2 (the rock price, 1.33–1.54×) and P2-3 (the floor's level) still own the
2.27×.** This step's honest value is a derived coefficient, a closed book, and
`(c2)` as a prediction rather than a setting — plus the correction that the wall
coefficient is, if anything, **too low near the corner**, not 3.4× too high.

⚠ Carried, for the next packet to fix: `side_h_used_*` is accumulated over the
whole `r ≥ 40 mm` gate, which includes the **36 %** of wall faces above the
nozzle plane that `side_face_factor_above = 0.0` already holds at 0 W. Those sit
high on the wall on the duct branch, so they drag the reported p50 to **171.89**
where `(a3)`'s live-band mean is **242.2**. The energy rows are unaffected
(a 0 W face contributes nothing to `surf_side_share × ledger_P_robin`) and
`(b2)`'s span passes either way at **3.029**, so this was reported rather than
re-run: correcting it costs a rebuild, a second full tree gate and a second
32-minute Meier leg, and changes no scored row.

⚠ Also carried: the fixture's isotropy monotonicity row still fails, with
**identical values** (`['1.3','5.3','4.7']`, same line) in this run and in both
archived P1g-2c fixture logs — a standing red, not a regression.
