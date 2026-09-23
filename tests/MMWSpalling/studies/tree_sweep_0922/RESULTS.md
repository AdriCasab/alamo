# Working-tree reference sweep, 2026-09-22 — PASS

**Purpose.** Certify that the current working-tree source is byte-identical to the
D2e reference with all new keys off. This was owed since 2026-09-22 00:01, when
another session added AMR surface-band refinement to `MMWSpalling.H` and
`MMWSpalling/Removal.H` and it touched D2k's own kernel line (`side_on && owner`,
:1369) — so **D2k's identity certified the 19:11 source, not what is on disk.**
It blocks any rebuild of `bin/mmwspalling-3d-g++`, AMR-related or not.

**Result: PASS. 2639/2639 byte-identical** (`witness/tree0922.md5` vs
`witness/ref_post_d2e.md5`, `diff` empty). All 11 witness tests rc = 0.

## Method

- **No rebuild.** The unit-test scripts hardcode `bin/mmwspalling-3d-g++`
  (e.g. `unit/robin_pinned/test:44`) with no override, so the binary was swapped
  rather than the scripts edited:
  - `bin/mmwspalling-3d-g++` (D2k campaign build, `62dea451…`) backed up to
    `bin/.mmwspalling-3d-g++.d2k-backup`;
  - `bin/mmwspalling-3d-g++-amr` (`2eeda1a5…`, built 00:01:56 from source stamped
    00:01:03, i.e. **the current tree**) copied into its place;
  - restore chained into the **same** command as the run, so no window existed in
    which it could be forgotten. **Verified restored: `62dea451…`.**
- **Frozen folders untouched.** `d2k_sideface/witness/` was copied here with
  paths rewritten rather than run in place (the `score.py` trap).
- Timings: `jet_d2c` 141 s, `dev2d` 44 s, then `robin_face` 5 s, `robin_pinned`
  16 s, `robin_jet` 21 s, `jet_enthalpy`, `robin_feet`, `beam_void_closure` 8 s,
  `scalar_flaw` 9 s, `spall_event` 1 s, `s1_smoke` 6 s. 13:54–14:03.

## What this certifies, and what it does not

- **Certifies:** with `surface.refine_depth` unset (default 0) and
  `surface_patch.side_face_flux` unset (default 0), the current source reproduces
  the D2e reference exactly. The AMR delta and the D2k delta are both key-off
  exact. **The tree is safe to commit and safe to rebuild.**
- **Does not certify:** anything with AMR *on*. `amr_band/` gated it on the
  hot-disk exclusion case only; the LI0 wiring identity, a long band gate and a
  two-level energy ledger are still owed before any scored AMR run. See ROADMAP.
