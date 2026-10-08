# P1f — MANIFEST

```
0bc96220d568c8fb7e456acadbb7b4ca7d24fc4e9d6f5ddea1140493f96ede2e  bin/mmwspalling-3d-g++
```

⚠ Deliverable 0 was built and gated separately, before the lateral code existed:
`BINARY_d0.sha256` = `ed7f050d64c24a4f33c86ed430a2144a7d13e03b8b7694356c4cf5ccf0adad7c`,
(a0) **0 of 2639** against `p1e.md5` and `p1b.md5` (`p1f_d0.md5`).

## Legs run

| leg | purpose | argv | outcome |
|---|---|---|---|
| `output_rev7sc/` | deliverable 7, the no-code ignition control | rev 7's own `cmd.txt` + `side_face_staircase=1` only; all five criterion/support keys asserted absent | **500.0 s, rc = 0**, 194,771 events, 398 mm |
| `unit/lateral_drive/` 8 legs | (c) | one input, geometry by `spall.seed_void` | see RESULTS |
| Kant `S150_2.0`, re-run | (a1) | the cheapest leg, on the d0 binary | **byte-identical** `5d3a2b9666bc5ad9b247443f` |

## NOT run, and why

**Leg A (the thesis leg) was not launched** — criterion (q), the run-order gate,
measured offline that **88.9 %** of P1e's viable faces lie inside the `r >= 40 mm`
scope cut, at a median radius of 25.8 mm. Leg A would have fired ~20 faces for a
reason with nothing to do with the thesis, so (f)/(g) would be unattributable.

**Leg B was not launched** — criterion (c) is a gate and it failed. The packet's
run order forbids launching a long leg after a gating failure, and the mode's
isotropy is the premise Leg B's results would rest on.

## Frozen inputs read

| leg | sha256 (24) |
|---|---|
| rev 7 | `28ca953426c870d7854647a1` |
| P1e | `b42716a39211e24fad5eb797` |
| rev7sc (this packet's control) | `d4d6ba671e8b38b81505247b` |
