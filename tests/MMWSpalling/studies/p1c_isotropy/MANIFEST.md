# P1c -- MANIFEST

⚠⚠ **NO REBUILD IN THIS STEP.** The binary is P1b's, unchanged, and no file under
`src/` was touched. The toolchain is not byte-reproducible, so rebuilding would
produce a different binary for no reason and forfeit P1b's one-binary discipline.

```
9e35d750c931f6ab3b5456131360f4a04339799e5205606505950f2efc3803c4  bin/mmwspalling-3d-g++
```

Criterion (a): this is byte-for-byte P1b's
`9e35d750c931f6ab3b5456131360f4a04339799e5205606505950f2efc3803c4`, and
`find src -newer bin/mmwspalling-3d-g++` returns empty.

⚠ `git diff --stat src/` is NOT the right check and the packet's wording invites a
mistake here: nothing in this repo is committed, so that command shows P1b's own
cumulative work against HEAD and would look "changed" however correct this step
is. The mtime test above is the one that answers "did source move since the
binary was built".

## Frozen inputs -- read-only, imported by path, never re-run

| leg | file | sha256 (24) |
|---|---|---|
| rev 7 (healthy) | `p0b_one_surface_T/output_fh_rev7/FH_removal_events.csv` | `28ca953426c870d7854647a1` |
| P1 (stalled) | `p1_gradient_criterion/output_fh/FH_removal_events.csv` | `26d8996f26ae2c90c57989a4` |
| P1b (min_sides=2) | `p1b_support_rule/output_fh/FH_removal_events.csv` | `5031129ca4e47710ed11f62b` |

The first two match the values the packet states (`28ca953426c870d7854647a1`,
`26d8996f26ae2c90c57989a4`), so the frozen studies have not moved.

## The one scored leg

P1b's `output_fh/cmd.txt` **verbatim** with exactly one key changed. Verified
mechanically, not asserted -- a token-level diff of the two argv files shows only:

```
< spall.support_min_sides=3        (P1c)
> spall.support_min_sides=2        (P1b)
< plot_file=.../p1c_isotropy/output_fh/FH
> plot_file=.../p1b_support_rule/output_fh/FH
```

`spall.support_rule=1`, `spall.support_step_min=0.006` and `stop_time=500.0` are
unchanged. The exact argv is in `output_fh/cmd.txt`; `output_fh/run.sh` launches
it and asserts `caffeinate -dis` itself, so it cannot be forgotten at the call
site.

| leg | output | event-log sha256 (24) | rows |
|---|---|---|---|
| P1c `min_sides=3` | `output_fh/` | `fff2513bb2a0c62ba95f33b6` | 113829 |

⚠⚠ **The leg ABORTED at t = 475.2 s** on a nozzle collision at column
(8, 10), r = 27.0 mm, stand-off s = -8.87e-14 m
(`MMWSpalling.H:2231`). It is scored as it stands; criterion (d) is the
primary gate and it FAILED, so per the packet no further long leg was
launched. ⚠ Criterion (f) is therefore **reported, not scored** --
witness 6 requires (d) to read 500.0 s before roundness is scorable, and
P1's own wall qualified only 4.5 % at 177 s.
