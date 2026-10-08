# `frozen_prechange/` — the pre-change reference `(b1)` is scored against

Full plotfile directories, not hashes. The manifest stores fingerprints only, so
rev 6 could tell that files had changed and never which **component** changed;
that is what blocked it. `fabcmp.py` and `score_restraint.py` compare against
this directory.

**Produced by the A-absent build** — changes B, C and D landed, change A reverted
— `BINARY_A_absent.sha256` = `2dfc6c74…`. The genuine pre-P2-2a plotfiles were
overwritten by rev 6's own gate run on 2026-10-07 at 16:36, with no archive, so
"freeze before editing source" arrived one revision too late. The A-absent build
replaces them and is **better than a snapshot**: it is regenerable, and it
isolates change A instead of the whole packet.

**Verified against the frozen `p2_1` pin: 363 of 369 manifest members match**
(`jet_enthalpy` 288/288, `beam_void_closure` 169/169, `dev2d` 64/64). ⇒ changes
B, C and D are **byte-inert**, which rev 6 asserted and never measured.

## ⚠⚠ Why `M-100_dz1` is NOT here

It was frozen, and it was **wrong** — 6 of its files were byte-identical to
**rev 6's output** rather than to the pin, because the freeze happened before I
diagnosed a stale read, and the leg was never re-frozen. Review row 13 caught it.

The stale read: `studies/s1_surface_resolution/run_sweep.py` **skips** any leg
carrying a `.done` marker unless `--force` is passed (`run_one` returns
`"skip (done)"`). The tree gate invokes it as `--smoke --force`; I had run plain
`--smoke`. So nothing re-ran and I compared rev 6's leftovers against the pin,
then concluded the leg was compile-sensitive. ⚠ The run-to-run "reproducibility"
check I ran to *validate* that conclusion was memoised the same way, so it
confirmed nothing.

**Measured properly, on the shipping binary through the gate (`--force`):
`M-100_dz1` matches the `p2_1` pin on 14 of 14 files.** So its pre-change
reference **is the pin** — the current output is byte-identical to it — and a
frozen copy would add nothing while a *mislabelled* one is an active trap for
the next reader. It is removed rather than re-frozen, deliberately.

⚠ `M-2000_dz1` **is** kept: it matches the pin, and it is the control that makes
the M-100 story legible (`alpha = 2000 /m` gives `exp(-500)`, which underflows to
exactly zero, against `alpha = 100 /m`'s `exp(-25)`, which does not).

## What is here

| leg | why |
|---|---|
| `B05 B1 B1e B1n B2 B3 B4 E Eon G Gbug` (+ `P1 P3`) | the strip-derived `jet_enthalpy` legs — the only files `(b1)` authorises to move, and the 26 that did |
| `QF` | a **quarter** leg that moved under rev 6 and must now read zero |
| `beam_void_closure` | moved under rev 6; must now read zero. ⚠ its own test is the carried standing `rc = 1` red |
| `dev2d` | the Meier pilot's 2-D dev leg; moved under rev 6, must now read zero |
| `M-2000_dz1` | the bit-stable absorption control |
