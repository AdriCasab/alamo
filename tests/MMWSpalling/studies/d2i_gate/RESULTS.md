# D2i results — the support-rule mesh gate with the idle clock off (2026-09-20/21)

**Verdict: FAIL.** With `pinned_idle_cycles` off, the 2 mm / 1 mm pair is **−13.9 %** over the
matched window 150–572 s, and the 100 s blocks run −2.1, −4.0, **−16.3**, **−36.1 %**. The gate
needed every one within 5 %. The idle clock was a real mechanism and removing it delays the
divergence by about 150 s, but **a second, slower mesh-dependent failure is underneath it**, and
it is not the clock: the idle-column count is **0 at both meshes in every block**.

D2g Stage 2 (the scored Meier runs) must not be planned on this rule. Nothing from D2b–D2d may
be quoted as a Meier result.

## 0. Criterion and binary

- `CRITERION.md`, sha256 `1c05225f80021d77814dfca4a4027edddd558efe0c6a5078bbd71201915ce2e8`,
  mtime **2026-09-20 19:09:58**, read-only, written **before** any run scored here (the output
  directory did not exist at that moment; first run directory 19:10).
- **Amendment on the record.** A first version (sha256 `051e93a5…`, 19:01:50) was written
  against the packet as it read then; the user revised `ACTIVE_STEP.md`, and the criterion was
  amended to match — the sensitivity run moved from 2 mm to 1 mm and the plotfile spacing from
  100 s to 50 s. **The deciding text (pair, ROP estimator, t_m, gap, conditions (a) and (b)) is
  byte-unchanged between the two versions.** The two runs started under the first version were
  killed at 164.8 s (2 mm) and 11.7 s (1 mm) and their output deleted; every run below starts
  from t = 0 after the amendment. `CRITERION.md` carries the same record.
- Binary `bin/mmwspalling-3d-g++` sha256 `6d047b50bdcde619952aff907d88f89ed0cae8f9d3ed645b43ba0ac2cd647005`,
  verified at the start and at the end. **No source edit, no build, no relink.**
- Harness `run.py` loads d2h → d2f → d2e → d2c by path. Nothing frozen was edited.

## 1. Identity

| run | rows (t ≤ 349.9 s) | vs d2h `I0` |
|---|---|---|
| `LI0_2mm` | 3645 | byte-identical |
| `LI0_1mm` | 3645 | byte-identical |

**PASS.** Both legs are exact extensions of the D2h diagnostic runs, so everything D2h reported
to 350 s carries over unchanged and the new information starts at 350 s. `I20_1mm` sets a
different key and has no identity partner.

## 2. Runs

| run | keys | status | t_end | wall | computing | disk |
|---|---|---|---|---|---|---|
| `LI0_2mm` | Q09 + `pinned_idle_cycles = 1e9` | bottom watchdog | 571.8 s | 15.7 min | 15.7 min | 2.4 GB |
| `LI0_1mm` | same, 1 mm | bottom watchdog | 591.1 s | 304.6 min | ≈ 188 min | 19 GB |
| `I20_1mm` | Q09 + `pinned_idle_cycles = 20`, 1 mm | stop cap | 349.9 s | 232.2 min | ≈ 115 min | 13 GB |

All three launched together as three four-rank jobs under `caffeinate -dis`. Suspended by
SIGSTOP 19:40:03–21:36:59 at the user's request and resumed by SIGCONT (`PAUSED.md`); that is a
suspension, not a restart, so the runs resume bit-identically, but it inflates `wall_s` by
≈ 117 min for the two 1 mm runs — hence the separate computing-time column. Plotfiles every
50 s at both meshes with the full D2h field list, 12 per run. `pinned_idle_cycles` is set per
input; the source default (2) is untouched.

## 3. The pair (t_m = 571.8 s)

Matched window 150–572 s: 2 mm **1.487** m/h, 1 mm **1.280** m/h, gap **−13.9 %** — outside 5 %,
so **(a) fails**.

| 100 s block | 2 mm ROP | 1 mm ROP | gap | D2g-1 reference (idle on) |
|---|---|---|---|---|
| 150–250 s | 1.491 | 1.461 | −2.1 % | −1.8 % |
| 250–350 s | 1.505 | 1.444 | −4.0 % | −12.5 % |
| 350–450 s | 1.498 | 1.254 | **−16.3 %** | −40.2 % |
| 450–550 s | 1.426 | 0.911 | **−36.1 %** | −64.8 % |

**(b) fails** on the last two blocks. Growth diagnostic: **−11.45 points per 100 s** (D2g-1:
−21.7). 50 s blocks, for continuity: +2.0, +0.1, −1.9, −0.2, −3.4, +0.0, −7.3, −13.1, −17.9,
−29.6, −39.4 %. The 2 mm leg is steady at 1.43–1.51 m/h across the whole run; the entire gap is
the 1 mm leg slowing down.

`d2i_blocks.png` overlays this on D2g-1's L09 and D2h's I0: LI0 sits on I0 to 350 s, as the
identity check requires, then leaves the band at the same place L09 did, on a shallower slope.

## 4. Per 100 s block (2 mm / 1 mm)

| block | rim spall 38–40 mm/min | annulus mean mm/min | 40–60 mm removal cm³/s | pinned share | idle cols | T_rec K | s_c mm | reach mm | P_face W | max abs ledger_err |
|---|---|---|---|---|---|---|---|---|---|---|
| 150–250 s | 24.8 / 24.3 | 28.3 / 27.7 | 1.77 / 1.67 | 0.209 / 0.205 | **0 / 0** | 1558 / 1560 | 136.4 / 136.0 | 62.0 / 61.4 | 1639 / 1630 | 3.4e-15 / 6.5e-15 |
| 250–350 s | 25.1 / 24.0 | 28.2 / 27.1 | 1.12 / 0.83 | 0.152 / 0.144 | **0 / 0** | 1618 / 1624 | 151.5 / 150.9 | 53.0 / 51.7 | 1332 / 1303 | 7.4e-15 / 7.4e-15 |
| 350–450 s | 24.5 / 20.8 | 27.7 / 24.9 | 0.72 / 0.32 | 0.127 / 0.116 | **0 / 0** | 1654 / 1667 | 161.5 / 163.5 | 48.6 / 46.6 | 1164 / 1102 | 9.0e-15 / 6.9e-15 |
| 450–550 s | 24.5 / **14.9** | 26.9 / 21.0 | 0.49 / **0.09** | 0.115 / 0.103 | **0 / 0** | 1677 / 1697 | 170.0 / **178.9** | 46.3 / 43.8 | 1056 / 959 | 9.1e-15 / 6.0e-15 |

`rim.py` reconstructs `foot_z` from the removal log to within 0.000 mm at both meshes, so the
support rule is being applied exactly as the analysis models it. The energy ledger closes to
≤ 9.1e-15 everywhere: this is not an accounting error.

### Silent fraction [%] by 4 mm band, and pad firings per column per 100 s

| run | t | 24–28 | 28–32 | 32–36 | 36–40 | 40–44 | 44–48 | 48–52 | fir/col |
|---|---|---|---|---|---|---|---|---|---|
| 2 mm | 250 s | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 46.9 |
| 2 mm | 350 s | 0 | 0 | 0 | 0 | 0 | 0 | 52 | 46.7 |
| 2 mm | 450 s | 0 | 0 | 0 | 0 | 0 | 28 | 100 | 45.9 |
| 2 mm | 550 s | 0 | 0 | 0 | 0 | 0 | 80 | 100 | 44.6 |
| 1 mm | 250 s | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 92.1 |
| 1 mm | 350 s | 0 | 0 | 0 | 0 | 1 | 25 | 99 | 90.3 |
| 1 mm | 450 s | 0 | 0 | 0 | 3 | 35 | 99 | 100 | 83.2 |
| 1 mm | 550 s | 0 | 1 | **11** | **31** | **80** | 100 | 100 | 70.0 |

t_cell = 7.65 s at 2 mm, 3.82 s at 1 mm. Silence in the 44–52 mm bands is the `s ≤ 0` rule above
the nozzle plane and does not count against the gate; it appears at **both** meshes.

## 5. Idle sensitivity at 1 mm — the clock is a threshold, not a dial

| `pinned_idle_cycles` | run | ROP 250–350 s | vs off | idle cols at 350 s | silent % 24→52 mm at 350 s |
|---|---|---|---|---|---|
| 2 | `d2g L09_1mm` | 1.322 m/h | −8.5 % | 593 | 3, 5, 12, 12, 46, 83, 100 |
| 20 | `I20_1mm` (new) | 1.444 m/h | +0.0 % | 0 | 0, 0, 0, 0, 1, 25, 99 |
| 1e9 (off) | `d2h I0_1mm` | 1.444 m/h | — | 0 | 0, 0, 0, 0, 1, 25, 99 |

`I20_1mm` is **byte-identical to `I0_1mm`** — `thermo.dat` and `_removal_events.csv` both — so at
1 mm a 20-cycle clock never trips once in 350 s. The setting does not trade rate against
stability: it is a threshold, and 2 lies below the interval between firings at 1 mm while 20
lies above it. **A finite but long clock does not cascade** (to 350 s), which is what needs to be
known if a physical idle timer is ever reintroduced. No value was chosen to change the rate.

## 6. Verdict and the fail-branch report

**FAIL**: (a) −13.9 %, (b) −16.3 % and −36.1 %. Identity PASS. Per `CRITERION.md`, stop.

**The new front, its position and timing.** A silent front still walks inward at 1 mm, starting
about 150 s later than with the clock on and moving more slowly. At 1 mm the 40–44 mm band goes
0 → 35 → 80 % silent between 350 and 550 s, 36–40 mm reaches 3 % at 450 s and 31 % at 550 s, and
32–36 mm reaches 11 % at 550 s. At 2 mm every band below 44 mm is 0 % silent for the whole run.
The first band to go is the one just outside the pads, and the pad band follows.

**What it is not.** Idle columns are 0 at both meshes in every block, so the clock is out of the
picture — this is a different mechanism with the same signature. The pinned share is nearly
equal between meshes (0.209/0.205 falling to 0.115/0.103), so it is not a pinning-rate
difference either. The ledger closes to 1e-15, so it is not accounting.

**The state at the divergence.** Between 350 and 550 s the 1 mm ring stops firing while the
centre keeps up: rim spall at 38–40 mm falls 24.3 → 14.9 mm/min while 2 mm holds 24.5; flank
removal over 40–60 mm collapses 1.67 → 0.09 cm³/s against 1.77 → 0.49; the feet end 23 mm behind
(195 vs 172 mm at 550 s) while the centre floor is only 9 mm apart (57 vs 48 mm); `s_c` runs away
to 178.9 vs 170.0 mm; `jet_r_reach` falls to 43.8 vs 45.8 mm; and `jet_P_face` ends 959 vs
1056 W. That is a closed loop — the ring lags, the stand-off grows, the jet reaches less far and
delivers less power to the rim, so the ring lags further — and it is the "no centre equilibrium
at 1 mm" lead carried over unresolved from D2g-1.

**Leads for the next packet, named and not tested here** (and not drawn from the final block,
which sits against the bottom stop):
1. **The `s ≤ 0` exclusion above the nozzle plane.** It is the one rule already known to silence
   bands at both meshes, and the front starts in the bands nearest it (40–44 mm) and moves in.
   Its interaction with a deepening hole is the first thing to look at.
2. **The stand-off feedback.** Whether a centre equilibrium exists at 1 mm at all, and whether
   the stagnation-decay closure `min(1, 5D/s_c)` is the term that makes it mesh-sensitive.
3. **The half-cell removal increment** (`h_col ≈ dz/2 − a_f`, S1): recession per firing scales
   with dz, so a band that fires at a marginal rate resolves differently at the two meshes.
4. **The pinned criterion at a lateral edge** — the D2h mechanism (a self-limiting face form
   plus lateral loss ∝ 1/dx) can park a column sub-critically without the clock, once the column
   stops being re-pinned for any other reason.

A diagnosis pair does not need to run to the bottom: the split is fully visible by 450 s, about
2.5 h at 1 mm, and the two legs are identical to 350 s.
