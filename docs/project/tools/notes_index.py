#!/usr/bin/env python3
"""Index docs/project/NOTES_ARCHIVE.md by each note's lead phrase, with verified line numbers.
Run: python3 docs/project/tools/notes_index.py
"""
import re, sys, pathlib
P = pathlib.Path(__file__).resolve().parents[1] / "NOTES_ARCHIVE.md"
B, E = "<!-- INDEX:BEGIN -->", "<!-- INDEX:END -->"


def lead(line):
    t = re.sub(r"^[-*]\s*", "", line)
    m = re.match(r"\*\*(.+?)\*\*", t) or re.match(r"⚠+\s*\*\*(.+?)\*\*", t)
    t = m.group(1) if m else t
    t = re.sub(r"[`*⚠✅]", "", t).strip(" .—:")
    return t[:104]


def main():
    L = P.read_text().split("\n")
    b, e = L.index(B), L.index(E)
    body = L[e + 1:]
    n = 0
    for _ in range(50):
        idx = [B] + [f"- L{i + b + n + 3}: {lead(l)}"
                     for i, l in enumerate(body) if re.match(r"^[-*] ", l)] + [E]
        if len(idx) == n + 2:
            break
        n = len(idx) - 2
    final = L[:b] + idx + body
    bad = [l for l in idx if (m := re.match(r"- L(\d+):", l))
           and not re.match(r"^[-*] ", final[int(m.group(1)) - 1])]
    if bad:
        sys.exit(f"{len(bad)} pointers do not resolve -- nothing written")
    P.write_text("\n".join(final))
    print(f"notes index: {len(idx) - 2} notes, all pointers verified")


if __name__ == "__main__":
    main()
