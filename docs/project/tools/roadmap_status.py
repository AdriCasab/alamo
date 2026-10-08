#!/usr/bin/env python3
"""Regenerate ROADMAP.md's `## Current Status` list from ARCHIVE_DONE.md.

One bullet per archive step record, with the line it resolves to. Generated rather
than hand-maintained for one reason: a label-based dedupe of hand-written bullets
silently merged two DIFFERENT steps that share a label -- the melt-aware `P1` with the
gradient-criterion `P1`, and the two `Step 15b`s -- and the loss was invisible until
someone grepped for a word that had disappeared. Deriving from the archive makes that
class of error structurally impossible.

Steps with no archive record of their own are listed after, pointing at the verbatim
Current Status record.

Run:  python3 docs/project/tools/roadmap_status.py          (after any archive append)
Always run docs/project/tools/archive_index.py first -- these pointers must agree with it.
"""
import re, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
ROAD, ARCH = ROOT / "ROADMAP.md", ROOT / "docs/project/ARCHIVE_DONE.md"
RECORD = re.compile(r"^(Step\s+\d+[a-z]*\b|P\d+[a-z]*\b|D2[a-z0-9-]*\b|[EASC]\d+[a-z]*\b"
                    r"|Refactor\s+R\d\b|AMR\b|Item\s*\d+\b)")
ORPHANS = ROOT / "docs/project/tools/status_orphans.txt"   # steps with no record of their own


def main():
    A = ARCH.read_text().split("\n")
    recs = [(i + 1, l[3:].strip()) for i, l in enumerate(A)
            if l.startswith("## ") and RECORD.match(l[3:].strip())]
    if not recs:
        sys.exit("no step records found in the archive")
    try:
        vb = next(i + 1 for i, l in enumerate(A)
                  if l.startswith("## ROADMAP `## Current Status`"))
    except StopIteration:
        vb = None

    out = [f"- {h} — `ARCHIVE_DONE.md` L{ln}" for ln, h in recs]
    if ORPHANS.exists() and vb:
        for line in ORPHANS.read_text().split("\n"):
            if line.strip() and not line.startswith("#"):
                out.append(f"- {line.strip()} — `ARCHIVE_DONE.md` L{vb} "
                           f"(no record of its own; in the verbatim Current Status)")

    L = ROAD.read_text().split("\n")
    s = L.index("## Current Status")
    act = L.index("Active step:")
    hdr = ["",
           "Completed and passing — **one bullet per archive record**, generated so two different",
           "steps can never be merged under a shared label.",
           "⚠ Regenerate with `python3 docs/project/tools/roadmap_status.py` after any archive",
           "  append, and run `archive_index.py` first so the line numbers agree.", ""]
    ROAD.write_text("\n".join(L[:s + 1] + hdr + out + [""] + L[act:]))

    # verify every pointer lands on its heading
    A2 = ARCH.read_text().split("\n")
    bad = [ln for ln, h in recs if not A2[ln - 1].startswith("## ")]
    if bad:
        sys.exit(f"pointers do not resolve: {bad[:5]}")
    print(f"status list: {len(recs)} records + "
          f"{len(out) - len(recs)} orphans = {len(out)} bullets, all pointers verified")


if __name__ == "__main__":
    main()
