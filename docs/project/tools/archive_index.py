#!/usr/bin/env python3
"""Regenerate the index at the top of docs/project/ARCHIVE_DONE.md.

Two defects this script exists to prevent, both of which happened:
  1. The first index was built BEFORE being inserted, so every pointer was off by the
     index's own length. Line numbers are now computed at a fixed point and VERIFIED
     against the headings they name before anything is written.
  2. The second index called all 204 `## ` headings "completed-step records". 115 of
     them are sub-headings INSIDE records -- Guardrails, Commands, Review findings.
     Only step records are indexed now; sub-headings are findable once inside a record.

Run after any append:  python3 docs/project/tools/archive_index.py
"""
import re, sys, pathlib

P = pathlib.Path(__file__).resolve().parents[1] / "ARCHIVE_DONE.md"
BEGIN, END = "<!-- INDEX:BEGIN", "<!-- INDEX:END -->"
# A step record's heading starts with its programme label. Anything else is a
# sub-heading of the record it sits in.
RECORD = re.compile(
    r"^(Step\s+\d+[a-z]*\b|P\d+[a-z]*\b|D2[a-z0-9-]*\b|[EASC]\d+[a-z]*\b"
    r"|Refactor\s+R\d\b|AMR\b|Item\s*\d+\b|ROADMAP\b)")


def strip_index(lines):
    if lines and lines[0].startswith(BEGIN):
        for i, l in enumerate(lines):
            if l.startswith(END):
                return lines[i + 1:]
    return lines


def build(body, n_pad):
    recs, subs = [], 0
    for i, l in enumerate(body):
        if not l.startswith("## "):
            continue
        h = l[3:].strip()
        if RECORD.match(h):
            recs.append((i + 1 + n_pad, h))
        else:
            subs += 1
    steps = [r for r in recs if not r[1].startswith("ROADMAP")]
    snaps = len(recs) - len(steps)
    out = [f"{BEGIN}  regenerate: python3 docs/project/tools/archive_index.py -->",
           "# Archive index", "",
           f"{len(steps)} completed-step records"
           + (f" + {snaps} roadmap snapshots" if snaps else "")
           + ". Line numbers are exact as of the last",
           f"regeneration and are verified by the generator; re-run it after any append.",
           f"({subs} sub-headings inside records -- Goal, Guardrails, Commands, takeaways,",
           "review findings -- are deliberately not indexed; find them inside their record.)", ""]
    for ln, h in recs:
        out.append(f"- L{ln}: {h if len(h) <= 110 else h[:107] + '...'}")
    out += ["", END, ""]
    return out, recs


def main():
    body = strip_index(P.read_text().split("\n"))
    n = 0
    for _ in range(50):
        idx, recs = build(body, n)
        if len(idx) == n:
            break
        n = len(idx)
    else:
        sys.exit("index length did not converge")

    final = idx + body
    bad = [(ln, h) for ln, h in recs
           if ln > len(final) or not final[ln - 1].startswith("## ")
           or not final[ln - 1][3:].strip().startswith(h.split("...")[0][:40])]
    if bad:
        for b in bad[:10]:
            print("  BAD", b)
        sys.exit(f"{len(bad)} pointers do not resolve -- nothing written")

    P.write_text("\n".join(final))
    print(f"index: {len(recs)} step records indexed, {len(idx)} index lines, "
          f"{len(final)} lines total -- all pointers verified")


if __name__ == "__main__":
    main()
