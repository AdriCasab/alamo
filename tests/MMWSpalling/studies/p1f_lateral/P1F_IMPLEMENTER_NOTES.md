# P1f revision 7 — implementer notes and takeaways (preserved)

⚠ These were appended to `ACTIVE_STEP.md` at ~12:13 on 2026-10-02 and were
**overwritten at 12:38:50** when the packet was rewritten to revision 8. They were
**not** carried into `docs/project/ARCHIVE_DONE.md` (last written 09:27). Preserved here
so the record is not lost. Scoring detail is in `RESULTS_LEGB.md`; the (c3)/(c2)
estimator history is in `C3_MESH_FINDING.txt`.

## Files changed in revision 7

- `src/Integrator/MMWSpalling.H` — deliverables 0, 7 (`T_face_out`, the three `lateral_*`
  guards, the Parse-time cell-form/face-form self-test).
- `src/Integrator/MMWSpalling/Removal.H` — the lateral host MFIter pass, regime 6, the
  unconditional ledger term (`:2476`, `:2697-2700`).
- `tests/MMWSpalling/unit/lateral_drive/test` — added `ell_field()` (field reduction of
  the final plotfile, unconditioned on firing); switched **(c2)** and **(c3)** onto it;
  `ell_last()` retained only to print beside, with a docstring recording why it is dead.
- `tests/MMWSpalling/studies/p1f_lateral/C3_MESH_FINDING.txt` — withdrawal banner,
  `CORRECTION`, `ADDENDUM`.
- `tests/MMWSpalling/studies/p1f_lateral/RESULTS_LEGB.md` — the Leg B scoring.

## Takeaways

1. ⚠⚠ **`lat_ell_med` is dead as an estimator.** `Removal.H:671-678` means it over
   **firing columns only**, and `lat_ell[c]` is itself a per-column **max over faces**
   (`ReduceRealMax` at `:664` after `:660-662` zeroes the non-winning columns). Doubly
   conditioned, and the conditioning *set* turns over step to step — measured **+54.8 %**
   and **+39.7 %** one-step jumps at **unchanged `lat_fires`**. It inflated (c3) from the
   field's +6.11 % to +33.56 %. The name is also wrong: it is a **mean**.
2. ⚠⚠ **(e) and (m2) are ONE defect.** The regime-6 ledger term is present, unconditional
   and correct in form. The tell that it is not a missing term: the error is **negative**
   and **does not grow with cumulative firings** — it peaks at −8.67e-3 at t = 200 s and
   *shrinks* while firings accumulate. What it tracks is the interior-cavity population,
   which the control does not have at all. ⇒ **P1h (collapse) is a prerequisite for P1g.**
3. ⚠⚠ **A passing scope gate can be the failure.** (h4) is green — 0 of 224,359 firings
   inside 40 mm — and that is *why* the mechanism runs away: the cut excludes the viable
   population (q) found inside it and leaves the mechanism unbounded outside it on soaked
   far-field wall with no descending front to extinguish it. p50 radius 75 mm, max
   167 mm, `riser_depth` max 374 mm. Widening the radius makes it worse.
4. **Volume is the diagnostic no criterion asked for.** 289,017 voided cells vs the
   control's 88,628 — **3.26×** — with `n_voided == 1` on all 224,359 regime-6 events, so
   it is real cells, not bookkeeping. A one-key diff that triples removed volume is a
   re-specification, and nothing in the table caught it: (h4) passed, (r) read a *better*
   isotropy, and only the energy invariant went red.
5. **Score the criterion and its removal rule separately.** The per-face criterion is
   corroborated by two independent instruments at the same order — (r)'s in-scope
   34× / 9.0× and the field's 8.1× — while the regime-6 removal rule fails on energy,
   cavities and volume. Keep the criterion; re-specify the removal.
6. ⚠ **Regime 1 is untouched** (−50 % … −607 %, diagonals firing more). The lateral mode
   does not cure the vertical mode's azimuthal bias; it adds a second population.
7. ⚠ Both physics-held (c3) legs are still-**rising** transients at t = 60 s (last-50 %
   spread 35 % axis / 50 % diagonal) and the diagonal's containment is 2.4 decay lengths,
   so +6.11 % is an order, not a converged figure. **Archive before re-running into a
   directory name** on a non-byte-reproducible toolchain.
