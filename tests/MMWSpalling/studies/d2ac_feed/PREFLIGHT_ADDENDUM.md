# D2ac — PRE-FLIGHT ADDENDUM

`PREFLIGHT.md` (sha `6b4fc3cd…`) and `preflight.py` (`3be11f33…`) stay **frozen
and unmodified**. This file records an error in them, found by running, and the
forced key that follows. The rule is to report a frozen pre-flight's
disagreement, never to silently replace it.

## The error: `nozzle_collision_radius` under `prescribed`

**`PREFLIGHT.md` item 0(b) says the guard "defaults to 0.0 and cannot fire".
That is WRONG.** So is the packet's own version of the claim (`ACTIVE_STEP.md`
item 0(b) and criterion (b): "it cannot fire at any feed", assuming 0.002).

Measured at `MMWSpalling.H:7676-7677`:

```cpp
pp.query_default("nozzle_collision_radius", v.nozzle_collision_radius,
                 v.feet ? v.foot_r_inner : v.radius);
```

The default is **`foot_r_inner` under feet, but `surface_patch.radius` under
prescribed**. `A`'s argv never sets the key, and `surface_patch.radius = 0.2`,
so on a prescribed leg **every column in the domain is a collision column** and
the guard is maximally armed — the opposite of "cannot fire".

**What happened.** All four prescribed legs aborted at t = 120 s:

```
MMWSpalling.H:2046 (BuildFlameColumns) surface_patch nozzle collision:
stand-off s = -5.55112e-15 <= 0 at column (i, j) = (29, 0), time = 120
(z_nozzle = 0.5, z_face = 0.5)
```

Column (29, 0) is at **r ≈ 59 mm** — outside the ~85 mm hole, on rock the jet
never heats, still at its original z = 0.5. At 1.5 m/h the nozzle plane descends
50 mm in 120 s, reaching z = 0.5 exactly. So the plane meets unrecessed rock
beside the hole, which is **the normal and intended state** of a burner that has
entered its own hole. The code says so at `:7670-7674`:

> "Outside it a column may sit above the descending nozzle (unrecessed rock
> beside the hole); s <= 0 is then passed to the expressions, which must stay
> finite and > 0."

⚠ The packet's own "Refuting evidence on disk" flagged this — *"D2a2 ran a
prescribed 1.5 m/h feed and 4 of 9 runs hit the collision guard"* — and both the
packet and I read it as a drain-era artefact rather than as this guard's default.

## The forced key

```
surface_patch.nozzle_collision_radius=0.028
```

**Not chosen to make a number come out.** It is the code's **own feet-mode
default**, `foot_r_inner` = 28 mm = the Ø56 bore radius from the Meier drawing.
It keeps the guard doing its real job — aborting if the nozzle would enter rock
**inside the bore** — while letting rock beyond the burner body legitimately sit
above the descending plane.

**It changes no trajectory.** `rc2` appears at exactly two lines: its definition
(`:2023`) and the abort test (`:2044`). It feeds no flux, no recession, no
diagnostic. Verified by grep, stated here because "physics-neutral" is the kind
of claim that has to be checked rather than asserted.

## What this does NOT change

* **Criterion (b) is unaffected.** Pad contact is computed **offline** over
  28 ≤ r ≤ 40 mm against `z_n − 0.050`, exactly as `CRITERION.md` fixes it. The
  guard was never the test, and now it is demonstrably not even inert — but it
  still is not the test.
* **No scored metric moves**, because the key cannot move one.
* `PREFLIGHT.md` items (a), (c), (d), (e) stand unchanged. `nozzle_z0 = 0.550`
  stands: the plane passing the original surface at t = 120 s is the burner
  entering its own hole, not an error.

## Provenance

* `S0` (the identity leg) ran **before** this addendum and is **unaffected** —
  it is the feet-based `RB_D` argv with no shared block, rc = 0 in 10.9 min, on
  binary `ea833f2f…`. It is not re-run.
* The four prescribed legs are re-run with the key added. Same binary; **no
  rebuild**, so no identity re-proof is owed.
