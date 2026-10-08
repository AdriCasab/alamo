# D2z-1 PRE-FLIGHT — step zero, before any scored run

Checked before `CRITERION.md` was hashed (`12db5609fec1652c…`) and before any
scored output existed on disk.

## ⚠ PROCESS NOTE — I launched, then killed, before hashing

I launched `W2`/`W1` **before** writing `CRITERION.md`, which inverts the
packet's "hash BEFORE the first run" rule. I killed both legs ~3 minutes in,
**deleted `output/` entirely**, wrote and hashed `CRITERION.md` with **zero**
outputs on disk, and relaunched. No scored output predates the hash and none was
inspected before it. Recorded here rather than quietly corrected.

While killing them I hit a second thing worth recording: **macOS `pgrep` ignores
`-a` and prints bare PIDs**, so `run.py --kill` read an empty command line,
matched nothing, and reported `killed 0 process(es)` while both legs kept
running. Fixed to ask `ps -o command=` per PID. A kill that reports success
without killing is worse than one that errors.

## Item 0 — name the call site, show each leg reaches it

The whole-cell override is one line, `MMWSpalling/Removal.H:1355`, now with no
key in front of it:

```cpp
if (reg_choice == 1 && !sp_weibull_presc_T_face)
    h_c = dz;
```

| condition | how each leg satisfies it |
|---|---|
| `reg_choice == 1` | `spall.enabled = 1`, `vapor.enabled = 0` (`input_hotdisk:83`) ⇒ the spall regime is the only one that can win |
| `!sp_weibull_presc_T_face` | `surface_patch.mode = convective_flame` (`input_hotdisk:127`) ⇒ false by its definition at `Removal.H:410-412` |

`spall.removal_events_csv = 1` (`:88`) supplies the event log (c)/(c0) read.
**No leg fails a row.** `RB_D` uses `input_drilling` (the Meier configuration),
also `convective_flame`.

## ⚠⚠ THE PACKET CONTRADICTS ITSELF ON GOAL 3 — resolved, with evidence

* **Goal 3** says: "**a loud parse-time message plus an on-disk witness, not an
  abort**", justified by a ⚠⚠ block: an abort would stop two sweep members from
  running, making criterion (a) and Goal 6 unmeetable.
* **The Guardrails, Goal 5 and P4** still say **abort** ("Step 20's retirement
  under whole-cell is an abort with a message naming the discarded scaling").

These cannot both hold. **Goal 3's reasoning is verified on disk and is taken as
authoritative.** Exactly four inputs set `spall.surface_normal = 1`:

```
tests/MMWSpalling/unit/beam_void_closure/input:90
tests/MMWSpalling/unit/seed_void/input:90
tests/MMWSpalling/validation/meier/sp_meier_pilot/input_2d_dev:106
tests/MMWSpalling/validation/meier/sp_meier_pilot/input_2d_dev_melt:136
```

⚠ **The packet names the wrong two as the in-sweep pair**, which only
strengthens its own point. Checked against `witness/hash.sh` and `witness/wit.sh`:

| input | in the sweep? |
|---|---|
| `unit/beam_void_closure` | **YES** — in both `hash.sh`'s directory list and `wit.sh`'s test list |
| `validation/.../input_2d_dev` (the `dev2d` leg) | **YES** — `hash.sh` hashes `sp_meier_pilot/output/dev2d` |
| `unit/seed_void` | **no** — absent from both |
| `input_2d_dev_melt` | **no** |

So the two in-sweep members are `beam_void_closure` and `dev2d`, not
`beam_void_closure` and `seed_void`. Either way an abort kills two sweep members
at parse. ⇒ **implemented as a message plus a conditional thermo column**, and
P4 is restated accordingly in `CRITERION.md`.

## Sources — every cited line verified

| what | cited | verified |
|---|---|---|
| the override to make unconditional | `Removal.H:1355` | ✓ |
| recession + flip (the phi tie) | `Removal.H:1441-1452` | ✓ `phi_new = phi_old - h_c`, flip on sign change, **no `cos(theta)` at this site** |
| `inc_full`/`inc_full_active` decls | `Removal.H:459-467` | ✓ |
| surface-normal cos scaling | `Removal.H:690-691` | ✓ and `sn_cth`/`sn_icth` appear **only** at `:689-992`, all inside the scan |
| `sn_cos_col` build, gated on the key | `Removal.H:256-261` | ✓ `if (spall_surface_normal && sp_weibull_spall)` |
| the scan's surface lerp (retained) | `Removal.H:822-843` | ✓ |
| criterion's slope awareness | `MMWSpalling.H:4665-4669` | ✓ `icth = sn_cos_p ? 1.0/sn_cos_p[...] : 1.0`, **outside the scan** |
| parse must precede the early return | `MMWSpalling.H:4145` | ✓ |
| `extensive=false` not rank-reduced | `Integrator.cpp:1141-1144` | ✓ |
| no `reference/` directory anywhere | count 0 | ✓ **0** |

**Finding 1 of the packet is exactly right**: `h_c = dz` discards the scan's
result including `h_col * sn_icth` (`:992`), so recession loses its slope
scaling, while the firing decision keeps `icth` at `MMWSpalling.H:4668`, which
never goes through the scan. The capability degrades; it does not die.

## The witness, and what it measures

`sn_recession_superseded` = mean `1/cos(theta)` over each step's spall-firing
columns — the factor the discarded scaling **would have** applied. `1.0` means
nothing was lost on that step. It is registered under exactly
`spall.surface_normal = 1`, so it lands in the four inputs above and nowhere
else.

On the regression's staircase it reads **1.00 to 6.67**: on a steep wall the old
path would have receded up to **6.7×** further per firing than whole-cell now
does. That is a real loss and it is now on disk in every affected run rather than
in a log line nobody diffs.

## Re-pin procedure (written BEFORE the sweep ran — Goal 6)

⚠ There is **no `run.py`** and no `--check`/`--repin`. The sweep runs the tests
against the binary at its hardcoded path and writes an md5 list; a re-pin is a
**declaration** about which list is the reference, not a script flag.

1. `bash witness/guarded_post_all.sh d2z1` → `witness/d2z1.md5` +
   `witness/d2z1_wit/`. The guard fingerprints the binary before and after and
   fails loudly if it moved mid-sweep.
2. **Audit** `diff witness/tree0922.md5 witness/d2z1.md5`, classifying every
   differing path against criterion (b)'s four bullets. Kept at
   `studies/d2z1_delete/tree_repin_audit.txt`.
3. **New reference: `witness/d2z1.md5`**, declared in `witness/REFERENCE.md`.
4. **`witness/tree0922.md5` preserved untouched** as the superseded pre-D2z-1
   reference. Nothing deleted or overwritten.

## Expected header changes, stated before the audit

* `firings_per_void` — **every** `thermo.dat` in the sweep (Goal 4,
  unconditional). This is deliberate and is the reason the sweep is re-pinned.
* `sn_recession_superseded` — `beam_void_closure` and `dev2d` only.

Any **other** header change refutes P5.

## Behaviour changes stated up front

* **Spall regime only** — vapour and bit removal untouched (the override is
  inside the spall branch).
* **`prescribed_T` excluded** ⇒ **Hu and Kant unchanged**. ⚠ This is
  conservatism, not a justified split: D2z-0 takeaway 7 found the same
  `dz/2 − a_f` landing under `prescribed_T` **with a real gradient**, pointing at
  the gate location rather than the lerp. Confirmed again here — the regression's
  `PRESC` leg reads mean `h/dz = 0.4941`.
* **The molten-surface stall survives** — a column is a spall candidate only when
  the scan returns a positive depth, and a molten surface returns zero.
* **Under `spall.melt_aware`, a solid crust over buried melt now leaves one cell
  per firing instead of the whole crust.** Should be inert: the crust case needs
  melt to form and the spall thermostat pins the surface near 850 K. Reported in
  `RESULTS.md` against the melt-aware regression rather than asserted.
* **Step 20 keeps its slope-aware criterion and loses its slope-aware
  recession.**

## Smoke / regression before the campaign

New `unit/removal_increment` (rewritten for D2z-1): **14/14 PASS**, including
`h_applied == dz` exactly on every firing with **no key set**, `PRESC` still at
`h/dz = 0.4941` (scan retained), `sn_recession_superseded` present on `SN` and
absent on `W`, the message found in `SN.log`, and `surface_normal` still
**running** rather than aborting.
