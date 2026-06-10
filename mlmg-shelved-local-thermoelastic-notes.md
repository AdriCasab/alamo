# Shelving the elastic MLMG solve for brittle spallation

**Date:** 2026-05-30
**Status:** decision + first implementation step (Rossi experiment)
**Scope:** brittle thermal-spallation regime (end goal: MMW ablation)

## TL;DR

For the brittle-spallation criterion we do **not** need the heterogeneous
elastic MLMG solve. The stress that drives `Sp = K_I/K_Ic` can be reconstructed
algebraically from the temperature field via the 1-D-confined thermoelastic
relation

```
sigma_xx(z) = -E * beta * (T(z) - T_ref) / (1 - nu)        [phases[0]]
```

This is the Kant/von Rohr spallation formulation, it is already implemented and
validated in the code (it is the `prescribed_T` removal path that the Kant
onset test passes on), and it removes every blocker we hit:

- **No convergence cliff.** The late-run `MLMG failed` aborts (baseline died at
  t=22 s, mid at t=20.5 s once voids accumulated) disappear — there is no solve.
- **No resolution ceiling.** MLMG could not solve below ~0.4 mm cells (the
  0.20 mm run aborted on the first solve). The local form is a per-cell algebraic
  expression, so we can refine as fine as we like and finally do a real mesh
  convergence study.
- **Much faster.** The MLMG solve was the dominant per-step cost.

## Why this is physically appropriate (not just convenient)

Thermal spallation is driven by a **constrained heated skin**: a small hot spot
on a large cool body is laterally confined by the surrounding rock, so the
in-plane stress is the equibiaxial constrained value `-E*beta*dT/(1-nu)`. That
is exactly what the local form computes. The free-lateral FEM we were running
(`zlo_roller_321`) actually *relaxes* this stress (lateral expansion bleeds off
the confinement), so for the confined-drilling geometry the local form is
arguably **more** correct, not a downgrade. (This was already noted in
`input_pulsed_brittle`: "the free-lateral BC is relaxing the subsurface sigma_xx
below what a 1-D-confinement Kant flake needs.")

The Rossi input even confirms it backwards: to make MLMG fire, the test had to
replace free-lateral BCs with a **full lateral roller box** (`el.bc.type =
constant`, all 26 regions) precisely so the FEM would reproduce
`sigma_xx = -E*alpha*dT/(1-nu)` per cell — i.e. the MLMG was being tortured into
reproducing the local form anyway.

## What we give up

Committing to the 1-D-confinement assumption loses:

1. **Grain-scale stress heterogeneity** (quartz/feldspar elastic mismatch).
   Partly recoverable by evaluating the local form with per-cell phase `E, beta`;
   the Weibull flaw field already encodes the statistical weak-spots.
2. **Genuine stress redistribution around the evolving cavity.** Matters for
   cavity *shape* and post-onset evolution; second-order for *onset*.

If genuine elastic interaction is ever needed (grain-scale fracture, cavity
mechanics), MLMG comes back — but it would first need the convergence fix
(deeper coarsening / hypre / AMR-compatible operator). That is a separate,
larger workstream and is explicitly **out of scope** for the brittle-onset goal.

## Why this matters extra for the MMW-ablation end goal

Ablation is dominated by *thermal* removal (melt/vapour), not stress. Stress
drives only the **secondary** spallation channel (cool periphery, below-melt
cells). So sinking effort into fixing the elastic MLMG would be optimizing the
minor channel. Shelving it here is doubly justified: the spallation channel gets
a robust, validated, mesh-independent stress source for free, and the real
ablation effort (melt/vapour front, cavity MMW re-absorption, spall-vs-vapour
regime gating) is where the investment should go next.

## The architectural finding that made this clean

There are **three** consumers of the spall stress, and their selectors were
inconsistent:

| consumer | location | stress source (before) |
|---|---|---|
| h_spall removal depth scan | `Removal.H` `k_i_at_zcrack` | local thermoelastic when `surface_patch=prescribed_T` (gated by `sp_weibull_use_T_face`) |
| Sp_field surface diagnostic / **gate** | `MMWSpalling.H` `UpdateSpAfterMechanics` | `stress_mf` (MLMG), gated by the *different* flag `sp_prescribed_sigma` |
| cluster labeller (gates removal) | `UpdateSpClusters` | reads Sp_field (so inherits MLMG) |

So in Rossi the **gate** (`Sp_field >= 1`, from MLMG) and the **action**
(`h_spall`, from local thermoelastic) were using *different* stress models. They
roughly agreed only because the roller-box BC forced the MLMG result toward the
1-D-confined ideal. Making the Sp_field diagnostic use the same local
thermoelastic stress as the removal scan (a) removes the inconsistency and
(b) removes the only remaining MLMG dependency, so `el.type = disable` works.

## Implementation

New input flag (default preserves all existing behaviour byte-for-byte):

```
spall.stress_source = mechanics            # default: sample stress_mf (MLMG)
spall.stress_source = local_thermoelastic  # reconstruct sigma_xx from T(z)
```

- `UpdateSpAfterMechanics` gains a `local_thermoelastic` branch that reconstructs
  `sigma_at_depth` from `temp_mf` via `E*beta*(T-Tref)/(1-nu)` (phases[0],
  compression-positive), and skips all `stress_mf` access.
- `Removal.H`'s `sp_weibull_use_T_face` is extended to fire when
  `spall.stress_source = local_thermoelastic` even without `prescribed_T`
  (surface temperature then sourced from `temp_mf` top cell instead of
  `surface_patch.T_f`), so the same flag also decouples the MMW-beam case later.
- With the flag set, `el.type = disable` runs the full Sp / cluster / removal
  pipeline with no elastic solve.

## First experiment: Rossi without MLMG

`input_no_mlmg` variant: `spall.stress_source = local_thermoelastic`,
`el.type = disable`, everything else identical (keeps `prescribed_T` surface
ramp, so the removal scan's stress is unchanged — only the Sp_field gate moves
off MLMG). Expectation: matches the MLMG Rossi result closely (the roller-box
MLMG was reproducing the local form anyway), but runs faster, can't crash on
MLMG, and can be refined past the old 0.4 mm ceiling.
