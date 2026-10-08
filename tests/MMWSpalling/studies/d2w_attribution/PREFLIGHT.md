# D2w pre-flight — written BEFORE any campaign leg launched

This is an **attribution** packet correcting a confound in **my own D2v result**.
D2v's `M1` changed the confined `h` map **and** switched on the below-plane
side-face path at `f = 0.2`, yet I reported "`M1 − M0` is the funnel fix". Both
mechanisms are live and both are large. The packet is right; this pre-flight
records what I verified before completing the 2×2.

---

## 0. Reachability ledger — every leg, every call site

| # | behaviour | call site | reached because | legs |
|---|---|---|---|---|
| 1 | confined map | `h_f(x, y, r, s, t)` under `jet_enthalpy`, `:1778` true branch | `A`/`M1` set `h_expr`; `M0`/`B` keep the stale one | `A`, `M1` |
| 2 | side-face term | `:1655-1656` inside `col_in_patch` → `nfx+nfy>0` | `B`/`M1` set `side_face_flux=1`; `A`/`M0` set **no side key at all** | `B`, `M1` |
| 3 | per-region factor | `:1653-1654` on `fs[col_c] > 0.0` | `B` sets `side_face_factor_above=0.0` explicitly | `B`, `M1` |

### The key diff, generated and compared against the packet's table

```
A: + surface_patch.h_expr="(if(s>0.0,((1449.0*..."     [751 chars, VERBATIM]
   - surface_patch.h_expr="1449.0*((0.0075/..."        [M0's stale map]
   (nothing else; 0 side keys)

B: + surface_patch.side_face_flux=1
   + surface_patch.side_face_factor=0.2
   + surface_patch.side_face_factor_above=0.0
   (h_expr unchanged from M0's)
```

**Both rows match the packet exactly.** The harness asserts the invariants at
job-build time and raises rather than launching: `A` must carry **0** side keys
(N1), `B` must carry **3** and keep `M0`'s `h_expr` (N2), and `A`'s map must be
**byte-equal** to `M1`'s token (N5).

⚠ **`A` sets no side key, rather than zeroing one.** `:6099` aborts on
`side_face_factor` without `side_face_flux`, so the correct "off" cell is the
key being **absent**.

### The map string is copied, never rebuilt (N5)

Taken verbatim from `d2v_hmap/output/M1.log`'s argv line — **751 characters**,
exactly one matching token. Cross-checked: `d2v_hmap/run.py`'s generator
reproduces that string **byte for byte**, and it **differs** from `M0`'s. So the
copy is provably the same map `M1` ran, and re-derivation (which would make `A`
incomparable) did not happen.

### Empirical demonstration — a 20 s smoke, not a grep

| leg | `side_P_face` | `side_P_face_above` | `patch_P_robin` max | `jet_P_face` max |
|---|---|---|---|---|
| **`A`** | **column ABSENT** | column ABSENT | 1660.97 W | 1660.97 W (= tops, no side term) |
| **`B`** | present, max **180.16 W** | **0.00 max AND 0.00 min** | 3067.77 W | 3067.77 W |

`A`'s side columns are **not registered at all**, which is the N1 witness. `B`'s
`side_P_face_above` is **exactly zero in both directions**, so it neither heats
nor **cools** its above-plane faces — N4 avoided. `ledger_err` ≤ 1.24e-14 on both.

⚠ **The 20 s smoke proves the below-plane side path and nothing above the plane.**
Above-plane columns do not exist until `nozzle_z` drops below the original
surface (t ≈ 120 s on `M0`), the trap D2v documented. `side_P_face_above` on `B`
is therefore re-checked **after 120 s** in the campaign, not taken from here.

**No leg fails to reach any site, so the rule does not stop this packet.**

---

## 1. Reference reproductions — all exact

| quantity | packet | measured |
|---|---|---|
| burner descent 150–400 s, `M0` | 1.50 m/h | **1.50** |
| burner descent, `M1` | 3.12 m/h | **3.12** |
| burner descent, `M2` | 3.10 m/h | **3.10** |
| `patch_P_robin` @ 400 s, `M1` | 1063.8 W | **1063.8** |
| `side_P_face` @ 400 s, `M1` | 624.8 W | **624.8** |
| `side_P_loss` @ 400 s, `M1` | 28.9 W | **28.9** |
| `jet_P_face` @ 400 s, `M1` | 1717.5 W | **1717.5** |
| **the identity** | `jet_P_face = tops + sides + loss` | **1717.5 vs 1717.5, diff +0.00 W** |
| `jet_P_face` mean 150–400, `M0` / `M1` | 1428 / 1760 W | **1427.5 / 1760.3** |
| (a) `M0` / `M1` | 37.2 / 3.8 mm | **37.23 / 3.82** |
| median radial step `M0` / `M1` | 6.9 / 1.3 mm | **6.9 / 1.3** |
| ring 50–76 mm recession `M0` / `M1` | 32.5 / 0.5 mm | **32.5 / 0.5** |
| Ø @ 10/25/50, `M0` | 144.3 / 129.1 / 114.3 | **identical** |
| Ø @ 10/25/50, `M1` | 100.4 / 95.4 / 91.7 | **identical** |

**No disagreements.** Every number this packet rests on reproduces with my own
estimator, including the three-way power identity to +0.00 W.

### ⚠ What my D2v write-up got wrong, in my own words

1. **The attribution was confounded.** I reported `M1 − M0` as "the funnel fix"
   when `M1` moved two levers. The packet specified that leg, but I had the data
   to see it and did not flag it.
2. **I never reported the rate.** D2v's packet said "score shape and volume,
   never rate", and I complied — but burner descent went **1.50 → 3.12 m/h**,
   i.e. **out of Meier's 1.3–1.6 band to twice it**, and that is a first-order
   problem with the configuration I presented as a success. A guardrail against
   *scoring* a quantity is not a licence to leave it **unreported**.
3. **I quoted `jet_P_face` as if it were the top-face power.** It is the
   **total** drawn from the stream; the march gathers side faces too. The
   correct decomposition is tops 1063.8 + sides 624.8 + losses 28.9.

---

## 2. Design readings

**The metric is imported, not re-implemented.** `d2v_hmap/analyze.py` is loaded
by path and its `profile` / `diameter` / `zface` / `MEIER` are used unchanged, so
all four cells share one metric. Verified by re-deriving `M0` and `M1`'s (a),
median step and ring recession through it — all exact.

**ROP is the slope of `nozzle_z` over 150–400 s, never pit depth.** The burner
rides the feet, not the deepest point.

**Legs will end at different times.** `M1` hit the frozen harness's bottom
watchdog at 401.4 s; a slower `A` may reach 500 s. (a) is scored at the common
time, and each leg's end time and depth are reported. D2v established that (a)
is window-insensitive (`M0` reads 37.2 at both 401.6 and 500 s).

**No `src/` change**, so the binary is unchanged from D2v's, whose identity
sweep was 2639/2639.
