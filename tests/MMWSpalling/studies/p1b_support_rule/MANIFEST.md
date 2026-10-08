# P1b -- MANIFEST

Binary, pinned for every scored run below. ⚠ NOTHING was rebuilt after the
first scored run: the toolchain is not byte-reproducible and there is no
post-hoc repair.

```
9e35d750c931f6ab3b5456131360f4a04339799e5205606505950f2efc3803c4  bin/mmwspalling-3d-g++
```

## Criterion (a) -- the frozen tree

| set | files | differing |
|---|---|---|
| `p1b.md5` vs `p1.md5` | 2639 | **0** |
| `p1b.md5` vs `p0b7.md5` | 2639 | **0** |

`beam_void_closure` rc = 1 by standing decision.

## Criterion (b) -- the unit test, 8 legs, one binary

`tests/MMWSpalling/unit/support_rule/`, sha256 of each leg's event log:

| leg | seed_void | T_init | T_flame | h_conv | support_rule | sha256 (24) |
|---|---|---|---|---|---|---|
| `TOWER` | `1 0 55 3 3 59 0 1 55 0 3 59` | 820 | 820 | 1.0 | 1 | `e7643aba1f8371046a5f7858` |
| `TOWER0` | `1 0 55 3 3 59 0 1 55 0 3 59` | 820 | 820 | 1.0 | 0 | `0a14c7924128448e5f43f41a` |
| `WALL` | `1 0 55 3 3 59` | 820 | 820 | 1.0 | 1 | `4c30865b1d07ede8d4e8739a` |
| `FLAT` | `-` | 820 | 1500 | 5000.0 | 1 | `1387601663260e35fc70580f` |
| `FLAT0` | `-` | 820 | 1500 | 5000.0 | 0 | `091b3cb4e0838d112af54965` |
| `SLAB` | `-` | 820 | 820 | 1.0 | 1 | `4c30865b1d07ede8d4e8739a` |
| `T745` | `1 0 55 3 3 59 0 1 55 0 3 59` | 745 | 745 | 1.0 | 1 | `4c30865b1d07ede8d4e8739a` |
| `T750` | `1 0 55 3 3 59 0 1 55 0 3 59` | 750 | 750 | 1.0 | 1 | `e2ee83b9cc23ec5a539625da` |

## Criterion (k) -- the exposure test, offline, no compute

Inputs, both read-only and imported by path:

| leg | file | sha256 (24) |
|---|---|---|
| rev 7 (healthy) | `p0b_one_surface_T/output_fh_rev7/FH_removal_events.csv` | `28ca953426c870d7854647a1` |
| P1 (stalled) | `p1_gradient_criterion/output_fh/FH_removal_events.csv` | `26d8996f26ae2c90c57989a4` |

Reduction: `exposure.py` (two axes, `support_min_sides` x `support_step_min`).
Full output in `EXPOSURE_SWEEP.txt`. Cut: `r < 40 mm`, `k_top` replayed from
the event log, quarter-domain axis neighbours skipped (their step is zero by
symmetry). Denominator: columns inside the burner radius that have fired at
least once by t.

## Deliverable 4 -- the Meier rerun

P1's `cmd.txt` **verbatim** plus `spall.support_rule=1`,
`spall.support_min_sides=2`, `spall.support_step_min=0.006`, and a new
`plot_file`. The exact argv is in `output_fh/cmd.txt`; `output_fh/run.sh`
launches it and asserts `caffeinate -dis` itself.
