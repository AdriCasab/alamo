---
description: Plan the next coding step. Rewrites ACTIVE_STEP.md and refreshes ROADMAP.md from the long plan.
---

# /plan — Planner Phase

You are the **planner**. Your job is to refresh the focused task packet so a
later `/implement` call can code the next step without opening the long
historical plan.

You do **not** write production code in this phase. You only read context and
rewrite planning files.

## Inputs you should read, in order

1. `CLAUDE.md` (already loaded — project rules and ALAMO architecture notes).
2. `ROADMAP.md` — current status, active milestone, next milestones.
3. `ACTIVE_STEP.md` — the current packet. Check its status first:
   - If status is `completed`, you are starting a fresh step. Read the
     `Claude completion notes` and `Implementation takeaways` so they can be
     archived.
   - If status is `ready for implementation` or similar, the step has not run
     yet — confirm with the user whether they want to re-plan or skip planning.
   - If status is `rejected by review` (or there is a `Review findings`
     section), revise the packet to address the reviewer's concerns; do not
     simply move on to the next step.
4. `docs/project/ARCHIVE_DONE.md` — only if you need historical detail to summarize a
   just-completed step or to disambiguate context for the next step.
5. **Long sources, read only the slice you need:**
   - `docs/project/in-main-tex-you-will-quizzical-treasure.md` — full historical plan and
     checklist. Read the section for the upcoming step.
   - `main.tex` — equations and physics. Read only the relevant section.
   - The relevant validation `.txt` file (Zhang/Oglesby, Hu, Kant/von Rohr) if
     and only if the next step is that validation.
6. Inspect current code/tests that the next step will touch — just enough to
   make `ACTIVE_STEP.md` concrete (file paths, existing patterns to mirror).

## What to produce

### 1. If the previous step is completed, archive it first

- Move the durable summary and `Implementation takeaways` from
  `ACTIVE_STEP.md` into `docs/project/ARCHIVE_DONE.md` (append; do not rewrite history).
- Update the relevant completed-step section in
  `docs/project/in-main-tex-you-will-quizzical-treasure.md` with a compact
  `Implementation summary` of **20 lines or fewer**. Summarize Claude's longer
  notes; do not paste raw logs, verbose diffs, or transient debugging detail.
- Update the `## Current Status` section of `ROADMAP.md` to list the
  just-completed step.

### 2. Rewrite `ACTIVE_STEP.md` for the next coding step

`ACTIVE_STEP.md` must contain everything the implementer needs without
opening the long plan. Required sections:

- `# Active Step: <N> - <short title>`
- `Status: ready for implementation`
- `## Context` — one short paragraph on why this step exists and what came
  before. Carry forward any caveats from prior steps that constrain this one.
- `## Sources` — exact pointers (file:section) into `main.tex`, the long
  plan, or validation `.txt` files. Quote constants, equations, and pass
  criteria directly here so the implementer does not need to open the long
  files.
- `## Goal` — numbered deliverables (files to add/edit, functions to
  implement, regressions to add).
- `## Criteria` — **a table, one row per test.** The planner writes **six** columns:
  `#` · `run order` · `test` · `bar` · `gate / report / verdict` · **`what a failure of this
  row means`**, and no others. Later phases **append** to the same table: `/implement` adds
  measured value, pass/fail and output path; `/verify` adds its own independently derived
  value and the file it read it from. ⚠ The last planner column is the highest-value one
  here — a step once shipped with a failing gate as a *finding* because the prediction
  said in advance what that failure would mean (log 2026-10-01 §17).
  - ⚠⚠ **CHOOSING THE COLUMN IS A DECISION WITH A TEST, AND IT IS NOT THE BAR'S
    IMPORTANCE.** A row is a **gate** only when **both** hold: (i) its claim is **yes/no
    about validity** — energy balances, a key is off by default, a frozen result still
    reproduces, the mechanism fires at all — and (ii) **something expensive comes after it
    in run order** for it to protect. A row that **measures how large an effect is** is a
    **report**, however central it is to the step: its bar is a *prediction*, and a wrong
    prediction is a **finding**, not a halt. A row after the last expensive run is a
    **verdict** — it decides whether the result is usable and what the step concludes, but
    there is nothing left for it to stop, so labelling it `gate` misuses the word.
    ⚠ **Never bundle a measurement into a row with yes/no checks.** Give it its own row.
    *(Log 2026-10-02 §5: a five-claim row carried four yes/no checks and one effect-size
    measurement under a single `gate`. The measurement's bar had been set without measuring
    the instrument's ~3 % precision, it read 8.4 % against a 5 % bar, and `implement.md`'s
    run-order rule — correctly — refused to launch either long leg. The four yes/no checks
    in that row had all passed. The 8.4 % was later found to be the test box's own
    insulated far wall: properly contained it reads **0.0057 %**. Cost: three planning
    cycles, and the step's actual question was never measured.)*
  - ⚠⚠ **A LEDGER ROW CARRIES THE UNIT OF EVERY QUANTITY THAT CROSSES A COMPARISON**, and a
    row that names both sides of a comparison without naming their units has not been
    checked. *(Log 2026-10-02 §6: the lateral criterion compared `K_I` in Pa·sqrt(m)
    against `K_Ic` in MPa·sqrt(m) — the one call site in `src/` that omitted
    `MPa_sqrtm_to_Pa_sqrtm` — inflating the score by 1e6 and firing 224,359 times on rock
    a factor of 9 below threshold. The ledger row named `K_I` and `KIcNasseri` side by
    side and said nothing about Pa vs MPa.)*
  - ⚠⚠ **A QUANTITY THAT DECIDES A BRANCH MUST BE AN OUTPUT, in stated units.** If the
    comparison a mechanism makes is not logged, **no bar in the table can catch it being
    wrong** — the step will report green rows over an invalid run. Give every new
    criterion a logged ratio, not just its consequences.
  - ⚠⚠ **A FIRING ASSERTION MUST BE TWO-SIDED.** "The mechanism fires" is a lower bound and
    cannot detect **over**-firing. Every new mechanism owes a **negative** case — a
    configuration where it must fire **exactly zero** — and an **upper** bound on its count.
    *(Same log: a passing "fires non-zero" witness was satisfied BY the bug.)*
  - ⚠⚠ **NEVER ADD A KEY TO A CONTROL'S ARGV TO SATISFY A GUARD.** A control leg's whole
    value is that it differs by the keys under test and nothing else. If a new mode needs
    a different criterion form or flaw, **make the mode follow the run's setting**; an
    abort that forces the argv is a planner error wearing a safeguard's clothes.
    ⚠ And check the **coupling**: a form and its matched flaw are one calibration, so
    swapping one moves the *other* mode's threshold. *(Log 2026-10-02 §8: a guard that
    aborted on `absolute` would have forced `driving_stress=differential` onto a control
    whose 8 µm flaw is matched to `absolute`, raising the VERTICAL firing threshold from
    756 K to 942–1085 K, destroying the one-key comparison and making the leg inert by
    construction. The packet's own prediction said "nearly inert" and treated it as a
    virtue — a design fault written down as a pre-registration.)*
  - ⚠ **A prediction of "this will be inert" is a design smell, not a safe prediction.**
    If a leg is expected to answer nothing, re-specify the leg rather than pre-register
    the null.
  - **Offline and zero-compute rows run first.** State the order; do not leave it to a
    warning comment in the Commands block.
  - **A bar carries the population and the time window its reference value was measured
    on.** The same run under two radial cuts gives two different numbers, and a ratio
    calibrated over all time can read **infinity** over a later window.
  - **Never write a bar as a bare ratio that can divide by zero.** Give it a
    minimum-count rule, or use a dispersion measure.
  - **Report every ratio with the magnitude it is a ratio of.** A band can show a
    tenfold improvement in uniformity while being essentially untouched.
  - **Normalise any azimuthal or per-sector count by the grid's own census** before
    scoring it.
  - ⚠⚠ **Diagnose with the FIELD, score with the RATE.** Event-log variables are
    recorded at the moment a cell fires, so they are conditioned on firing — and that
    conditioning can **invert** the comparison: the sub-population that fires in a
    starved direction reads *more* drive than the population that fires freely. A
    firing-rate asymmetry can only be diagnosed from the field, including the cells that
    did **not** fire.
  - ⚠ **A packet that touches removal geometry carries an isotropy witness**: the same
    quantity resolved per azimuth sector, against a healthy reference run's spread. It
    costs nothing offline and would have caught a grid-diagonal artefact before a 500 s
    leg was interpreted.
- `## Guardrails` — what *not* to touch, what must remain default-off, what
  earlier validation thresholds must remain green, etc.
- `## Commands` — exact build, run, and test commands for this step
  (mirror the patterns in `CLAUDE.md`).

Keep it tight: equations, constants, file paths, validation criteria, and
guardrails go in; verbose prose stays out.

⚠⚠ **Each requirement appears in EXACTLY ONE PLACE.** `## Guardrails` lists only what
this step introduces; a criterion's bar is stated in the Criteria table and nowhere else.
⚠ Avoid quoting a section's own name in body text: such a reference once matched during
an edit and deleted three sections (log 2026-10-01 §25).

### 3. Refresh `ROADMAP.md`

- `## Current Status` lists the just-completed step (if any), **one line per step**.
  Anything older than the current milestone collapses to **one line plus a pointer into
  `docs/project/ARCHIVE_DONE.md` naming the heading it resolves to**. ⚠ Before
  collapsing anything, append the text being removed to the archive **verbatim** and
  regenerate the index and the status list in the **same pass** — never leave a state
  where a pointer resolves to nothing:
  `python3 docs/project/tools/archive_index.py` then
  `python3 docs/project/tools/roadmap_status.py`. Both verify their own pointers before
  writing and refuse to write if any fails. ⚠ `## Current Status` is **generated** —
  do not hand-edit it; a hand-maintained list merged two different steps sharing a
  label and the loss was invisible (log 2026-10-02 §2). ⚠⚠ And check what sits **under** `Active step:`, not just that the
  label appears once — a check on the label passes while the section fills with
  completed blocks (log 2026-10-02 §1).
- `Active step:` names the new step and points at `ACTIVE_STEP.md`.
- **Budget check before finishing** (`CLAUDE.md` → Context Budget Rules): hand-written
  lines ≤ 350 excluding the generated status list —
  `awk '/^## Current Status/{g=1} /^Active step:/{g=0} !g{n++} END{print n}' ROADMAP.md`;
  `Active step:` + `Next milestones` ≤ 100 lines; `wc -l` ≤ 600. Over budget ⇒ release
  notes to the archive verbatim, never compress by lengthening lines.
- `Next milestones after the active step:` lists the next few in one-line
  form.
- `## Known Stale/Important Notes` holds **only caveats that constrain the active or
  next few steps**, each one bullet of at most three lines with a pointer to where the
  evidence lives. When a caveat is resolved, superseded, or no longer touches the next
  few steps, **move it to the archive with its date**. Cap: **40 bullets.**
- ⚠⚠ **A separate `## Standing rules` block, hard-capped at 10-15 bullets, is NEVER
  archived.** It holds programme invariants, not step caveats: unit conventions,
  scoring traps, and coupled calibrations. The release rule above would otherwise
  archive exactly the caveats that constrain steps not yet named — a coupled pair of
  constants, a measurement that must never be scored at a particular depth, a
  quarter-domain factor. ⚠ The cap is what matters: without it this becomes the new
  notes section.
- **Cite roadmap bullets; never copy their text into the packet.** (log 2026-10-02 §1)

## Verification before writing (mandatory)

Do these against the code and the frozen studies **before** drafting; carry 1–4
in the packet. **Guardrails are assertions too — a wrong one suppresses the check
that would have caught it, so never write one that forbids a cheap verification.**

1. **Reachability ledger** — one row per control/test pair: the exact key diff,
   the code site reading each key, and the guard condition there. A named run
   resolves to keys read from its `run.py`, never to what its name implies.
   - ⚠⚠ **A row for the guard that decides whether the site EXECUTES**, not only for
     the site. Every row of one packet's ledger was correct and the mechanism was
     silently inert, because the scan it hung off was gated by the very criterion it
     existed to overrule (log 2026-10-01 §17). **No grep catches this; a smoke run does.**
   - ⚠⚠ **A row for the CONJUNCTION a deliverable requires**, not one row per symbol.
     A correct row for a scope test beside a correct row for a temperature source,
     with no row for both, hid a variable that is identically zero exactly where the
     scope test is true.
   - ⚠ **If the step edits a shared function, state which branches the named tests
     reach**, and cover the rest in the unit test. Breadth is not depth.
2. **Legs as literal key lines**, not prose — "at default X" must be unwriteable.
3. **Sources** — every constant, measurement and literature number carries a
   location you opened *this session*. Take constants from where they are
   **consumed**, not declared.
4. **Refuting evidence** — list the on-disk runs that would **refute** each
   prediction. If one exists, the measurement governs and the conflict is the
   finding.
5. **"What would make this inert?"** — enumerate the non-physics explanations for
   a null and give each a witness diagnostic.
   - ⚠⚠ **A gate scoped to a sub-domain must state what it EXCLUDES and why the
     exclusion is safe.** One packet's inertness gate — "the rule must not fire on a
     healthy floor" — was exactly the right gate and was written with an `r < 40 mm`
     cut; the artefact lived on the wall just outside it, where the same rule fired
     **53 %** of a healthy run's columns. **The exclusion, not the bar, was the defect**,
     and it was measurable at zero compute from a file already on disk.
6. **Context carries only what is needed to execute**; argument and history go to
   the daily log and `ARCHIVE_DONE.md`. Keep the packet under 500 lines (`CLAUDE.md` → Context Budget Rules; raised from 300 by the user, 2026-10-07).
7. **Re-read the finished packet for INTERNAL CONSISTENCY only** — no deliverable
   contradicting a criterion, no number stated twice, every offline check ordered
   before every long leg, and the status line one sentence. ⚠⚠ **This is a self-check
   and cannot catch a planner blind spot.** The outside read is `/implement`'s step
   zero, which resolves the ledger against the actual code and stops on a failed row;
   that is what caught a packet-breaking defect. Do not duplicate step zero's job here,
   and do not rely on this item for it.
8. **Take measurements, not just constants, from where they are consumed.** A number
   derived from a log field must state the field's semantics — one event log's `k_top`
   is the top *before* that removal; reading it as the state *after* is wrong by one cell.
9. **Every stage that edits `src/` runs the source-reading lock tests before its gate
   row.** They take seconds, and a lock that pins literal text breaks when that text
   moves. *(Log 2026-11-10, P2-30: the packet ran the locks after stages 1, 2 and 4b but
   not 4a; 4a moved text one lock quoted, and it was caught a stage late.)*
10. **A bar that predicts how much a selection rule will match is measured first.** Run
   the rule on the current tree and state the count; it costs one command. *(Same log:
   "at least 400 comment lines" was predicted from a rule that matched 43; the corrected
   rule's 200–450 then read 81.)*
11. **A stage that edits an input deck invalidates every later byte-for-byte baseline for
   that deck.** Name the new baseline in the same packet. *(Same log: a deck's radius
   changed in stage 3 while stage 4 still compared it with the old output.)*

## Guardrails for the planner

- Do **not** edit production code, tests, or `src/`. This phase only touches
  `ROADMAP.md`, `ACTIVE_STEP.md`, `docs/project/ARCHIVE_DONE.md`, and
  `docs/project/in-main-tex-you-will-quizzical-treasure.md`.
- ⚠⚠ **`CLAUDE.md` and `.claude/commands/*.md` are NOT in that list, and the planner
  does not edit them on its own initiative** — not even to fix a rule it has just
  found wrong, and least of all these three specs, which govern its own conduct. A
  reviewer agent saying a change is "safe to run" is **not** the user's decision.
  Propose the edit, show the diff, and wait. ⚠ The 2026-10-02 spec revision was made
  this way: proposed, taken without asking, and then authorised retrospectively after
  the overstep was reported. Do not treat that as the precedent; the asking comes first.
- Do **not** rewrite the long `docs/project/in-main-tex-you-will-quizzical-treasure.md`
  unless the user explicitly asks. Append-only updates to the matching step's
  `Implementation summary` are fine.
- Preserve warnings about stale references, changed assumptions, input names,
  and validation caveats.
- Remove transient debugging detail from the new `ACTIVE_STEP.md`; keep only
  what the implementer needs.
- If you are unsure which step to plan next, ask the user before rewriting
  anything.

## End of phase

Report:
- Which step was archived (if any), with a one-line summary.
- Which step is now active, with its title and the one-line goal.
- Any unresolved caveats the implementer must respect.
