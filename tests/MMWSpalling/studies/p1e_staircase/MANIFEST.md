# P1e — MANIFEST

One binary for every scored run; nothing rebuilt after the first scored run.

```
5913ff51b6a2b49d21238accff4e7d7941833f8d1bcb0164be95ac545ce4f127  bin/mmwspalling-3d-g++
```

⚠ This step DOES change source (four `AreaOfDelta` call sites, one new key), so
unlike P1c the binary is new and the laboratory anchors are re-run rather than
inherited. `find src -newer bin/mmwspalling-3d-g++` is empty from the first
scored run onward.

## Criterion (a) — the default did not leak

| set | files | differing |
|---|---|---|
| `p1e.md5` vs `p1b.md5` | 2639 | **0** |
| `p1e.md5` vs `p1.md5` | 2639 | **0** |

Run with `surface_patch.side_face_staircase` **absent**. `beam_void_closure`
rc = 1 by standing decision.

## The three parse guards

| command | aborts with |
|---|---|
| `side_face_staircase=1 side_face_step_height=0` (flux on) | requires `side_face_step_height = 1` |
| `side_face_staircase=1` with `side_face_flux` off | requires `side_face_flux = 1` ⚠ **added beyond the packet** — without this row the key is silently ignored, since AMReX only *warns* about an unused ParmParse variable |
| `side_face_staircase=2` | must be 0 or 1 |

## Criterion (d) — the unit test, 6 legs, one binary

`tests/MMWSpalling/unit/side_staircase/`. ⚠ The factor is exact, so the test
carries **no tolerance**; it compares the **printed token** of `side_sc_mean`
against `"%g" % want`, because `thermo.dat` is text written with `%g` (six
significant digits) and a bit-exact comparison against the full double fails at
2.2e-07 — the file's rounding, not the kernel's error.

| leg | seed_void | staircase | expected factor | sha256 of event log (24) |
|---|---|---|---|---|
| `WALL` | `1 0 55 3 3 59` | 1 | exactly 1 | `f7c3d6ceee436d0592dc4d36` |
| `WALL0` | `1 0 55 3 3 59` | 0 | (reference) | `f7c3d6ceee436d0592dc4d36` |
| `PILLAR` | `1 0 55 1 3 59 3 0 55 3 3 59` | 1 | exactly 1 | `f7c3d6ceee436d0592dc4d36` |
| `CORNER` | `1 0 55 3 3 59 0 1 55 0 3 59` | 1 | exactly 0.707107 | `f7c3d6ceee436d0592dc4d36` |
| `FLAT` | `-` | 1 | never applies | `f7c3d6ceee436d0592dc4d36` |
| `FLAT0` | `-` | 0 | (reference) | `f7c3d6ceee436d0592dc4d36` |

⚠ WALL/WALL0 and FLAT/FLAT0 identity is checked on the **shared `thermo.dat`
columns**, not on the event log: with the key off the four diagnostics are not
registered, and a short leg can produce an empty log — which would make the
comparison vacuous (the trap P1b's FLAT leg fell into). Non-vacuity is asserted
explicitly: `side_P_face` for WALL, `ledger_E_H` for FLAT (⚠ `side_P_face` is
zero on FLAT by construction and is the wrong probe there).

## Frozen inputs — read-only, imported by path

| leg | sha256 (24) |
|---|---|
| rev 7 (control) | `28ca953426c870d7854647a1` |
| P1 (stalled) | `26d8996f26ae2c90c57989a4` |
| P1b (sides=2) | `5031129ca4e47710ed11f62b` |
| P1c (sides=3, THE CONTROL) | `fff2513bb2a0c62ba95f33b6` |

⚠⚠ **P1c is the single-variable control for this step** — same
`spall.support_min_sides = 3`, staircase OFF. P1b differs by **two** keys, so a
P1b-only comparison cannot separate the correction from the support rule. Every
reported delta names which leg it is against.

## The scored leg

P1b's `cmd.txt` verbatim plus two keys and a new `plot_file`; verified by a
token-level diff, which shows only:

```
< surface_patch.side_face_staircase=1     (added)
< spall.support_min_sides=3               (P1e)
> spall.support_min_sides=2               (P1b)
< plot_file=.../p1e_staircase/output_fh/FH
> plot_file=.../p1b_support_rule/output_fh/FH
```

`output_fh/run.sh` launches it and asserts `caffeinate -dis` itself.

| leg | output | event-log sha256 (24) | rows |
|---|---|---|---|
| P1e `staircase=1, sides=3` | `output_fh/` | `b42716a39211e24fad5eb797` | 126132 |

Completed cleanly: STEP 31250, **t = 500.0 s**, AMReX finalized, zero abort
or collision lines (P1c aborted at 475.2 s).

## The laboratory anchors — re-RUN, not re-reduced

This step changes source, and neither anchor is in the frozen tree, so both were
re-run on the new binary and compared **byte for byte**.

| anchor | result |
|---|---|
| `kant/sp_kant_onset_meiercfg` `--check` | **PASS**, every value reproduced |
| ⚠ the same leg **re-run** on the new binary (`S150_2.0`) | event log **byte-identical**: `5d3a2b9666bc5ad9b247443f4fd013a81d64e40153d58b430fcef902f5170bb5` |
| `p0_steady_column/test`, all 24 legs re-run | reproduces P1 exactly (A08 1.51 %, A80 4.30 %, S 2.85 %); exits FAIL on the deliberately-unrelaxed (a) bar, as in P1 |
| ⚠ its legs vs P1 revision 5's archive | **byte-identical** on every leg checked: `S_2.0` `25d92ccc198406ab`, `C_2.0` `1ab9b17c0dd3ae65`, `A08_2.0` `47068feac77c7389`, `A80_2.0` `258c617b2349f7e6`, `B3_2.0` `e8f271b51d8edd49`, `S_0.5` `b2a53d0d4934cfac` |

⚠⚠ **`--check` ALONE COULD NOT HAVE ESTABLISHED THIS.** It re-reduces output written
by the OLD binary, so it cannot see a regression caused by a source change. The
re-run is what makes the anchor meaningful, and the byte-identity is what proves
the four-site change is inert on the laboratory — on a toolchain that is NOT
byte-reproducible.

⚠ `p0_steady_column/test` **rmtree's each leg directory**, and those outputs'
hashes are recorded in P1 revision 5's manifest. They were copied to
`output.p1rev5_archive/` BEFORE the re-run, because an overwrite would have made
that provenance permanently unverifiable. The archive is what the byte-identity
row above compares against.
