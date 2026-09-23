# D2r-0c item 10 — the campaign as specified cannot run

**The mechanism and its baseline require mutually exclusive configurations.**
The code (items 1–9) is implemented and verified; only the campaign is blocked.

## The contradiction, with evidence

| requirement | what it needs | verified |
|---|---|---|
| the axial pass is *"an axial extension of the jet enthalpy march"* | `jet_closure = enthalpy` — `JetEnthalpyMarch` is called only under `if (surface_patch.jet_enthalpy)` (`MMWSpalling.H:1774`) | the pass lives inside that function |
| criterion (a), the only scored criterion, is `M` vs the frozen `A_f02` pair | the **arbiter** (`d2j0_hotdisk/input_hotdisk`) | `A_f02.log` line 1 → `input_hotdisk` |
| but the arbiter has **no jet at all** | — | `d2j0_hotdisk/output/S2/metadata:219`: `surface_patch.jet_closure = none` |

So a leg cannot simultaneously be the arbiter (to be comparable to `A_f02`) and
carry a gas stream (to have anything to march up the annulus).

**Confirmed by inventory of `d2l_scan/output/`:**

```
A_f01 A_f02 A_f05 A_f10  + _1mm partners   ->  input_hotdisk   (arbiter, NO jet)
M_f01 M_f02 M_f05 M_f10  (no _1mm at all)  ->  input_drilling  (Meier, jet)
```

- The **mesh pairs exist only on the arbiter**, which cannot run the march.
- The **jet exists only on Meier**, which has **no frozen 1 mm partner**, so
  criterion (a) has no frozen baseline there.

## This fails loudly, not silently

`surface_patch.jet_annulus` is parsed only inside the `jet_closure = enthalpy`
block; outside it the key hits the existing allowlist and **aborts** with
`surface_patch.jet_annulus requires surface_patch.jet_closure = enthalpy`.
Setting it on an arbiter leg therefore stops the run rather than producing an
inert one. That is deliberate — "never silently fall back".

## What the arbiter's above-plane heating actually is

`A_f02` sets `surface_patch.h_expr = "if(s>0.0,700.0,0.001)"` — above the plane
`h = 0.001`, i.e. off. The arbiter expresses the above/below split **through the
h expression**, and `d2j0b_plane/run.py` already builds exactly that pair
(`WALL[case] = (h_up, T_up)` → `if(s>0.0,H_IN,h_up)` / `if(s>0.0,T_IN,T_up)`).
Above-plane heating in the arbiter is a **prescribed field**, not a marched one.

## The two ways forward — the planner's call, not the implementer's

1. **Keep the arbiter and the frozen baseline; drop the march.** Score
   above-plane heating as a *prescribed* `h_up`/`T_up` pair at the derived band
   (142 / 250 against `A_f02`'s `f = 0.2`), one key different from the frozen
   baseline. Directly answers the packet's title question and reuses
   `A_f02`/`A_f02_1mm` unchanged. **Cost: ~30 min.** But it tests prescribed
   above-plane heat, not the march, and a prescribed field is unbounded — the
   D2o-0 failure mode the march exists to avoid. It is bounded here only because
   the arbiter's wall is small.
2. **Keep the march; move to Meier.** Run `jet_annulus` legs on `input_drilling`
   against a freshly-run `M_f02` mesh pair. Tests the real mechanism with
   conservation intact and the shield in play. **Cost: the 1 mm Meier leg at
   250 s is the expensive item (~75 min extrapolating D2l's 27.9 min arbiter
   leg and the 3 h/600 s Meier figure), and `M_f02_1mm` must be run too, so the
   baseline is no longer frozen.**

Note that D2q-2c already contains a partial answer to the packet's question:
its **no-plane control converged to 3.7 / 4.5 / 2.0 %** where the with-plane case
stayed at 11 / 22 / 30 %, differing in nothing but what the above-plane rock
receives. That is evidence *for* option 1's hypothesis, obtained without the
march.
