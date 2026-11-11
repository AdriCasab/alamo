#!/usr/bin/env python3
"""P2-30 stage 2 (deliverable 4): list every comment block (a maximal run of whole-line // comments) in the class's source
whose text matches the history wording. The sort itself is by hand (S/COMMENT_SORT.txt). Writes S/HISTORY_CANDIDATES.txt.
Each candidate: an id, file, first-last line, the enclosing member, and the block's lines numbered from 1."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "unit"))
import _src  # noqa: E402

WORDING = re.compile(r"deleted|retired|until P2-|no longer|formerly|superseded|replaced by|used to")
MEMBER = re.compile(r"^    (?:static |inline |virtual |constexpr |const )*(?:[A-Za-z_][\w:<>,\* ]*?[\s\*&]+)?"
                    r"([A-Za-z_]\w*)\s*\([^;]*$|^    struct\s+(\w+)")


def blocks(lines):
    i, n = 0, len(lines)
    while i < n:
        if lines[i].lstrip().startswith("//"):
            j = i
            while j + 1 < n and lines[j + 1].lstrip().startswith("//"):
                j += 1
            yield i, j
            i = j + 1
        else:
            i += 1


def enclosing(lines, i):
    """the member a block sits in or stands above: the next member signature if the block is at member depth, else the last"""
    if lines[i].startswith("    //") and not lines[i].startswith("     "):
        for k in range(i + 1, min(i + 60, len(lines))):
            s = lines[k]
            if s.strip() == "" or s.lstrip().startswith("//") or s.strip().startswith("AMREX_GPU"):
                continue
            m = MEMBER.match(s)
            if m:
                return m.group(1) or m.group(2)
            break
    for k in range(i - 1, -1, -1):
        m = MEMBER.match(lines[k])
        if m and (m.group(1) or m.group(2)) not in ("if", "for", "while", "switch", "return"):
            return m.group(1) or m.group(2)
    return "(file head)"


def main():
    src = _src.tree()
    names = [_src.PARENT] + [_src.SLICE_PREFIX + n for n in _src.slices(src)]
    out, tot, tl, all_blocks, all_lines = [], 0, 0, 0, 0
    per = {}
    for rel in names:
        L = src[rel].split("\n")
        for a, b in blocks(L):
            all_blocks += 1
            all_lines += b - a + 1
            txt = "\n".join(L[a:b + 1])
            if not WORDING.search(txt):
                continue
            tot += 1
            tl += b - a + 1
            per[rel] = per.get(rel, 0) + 1
            cid = "C%03d" % tot
            out.append("%s  %s:%d-%d  member=%s  lines=%d" % (cid, os.path.basename(rel), a + 1, b + 1, enclosing(L, a), b - a + 1))
            for k, s in enumerate(L[a:b + 1], 1):
                out.append("  %3d| %s" % (k, s.strip()))
            out.append("")
    head = "history-worded comment blocks: %d of %d blocks; %d of %d comment lines (wording: %s)" % (
        tot, all_blocks, tl, all_lines, WORDING.pattern)
    open(os.path.join(HERE, "HISTORY_CANDIDATES.txt"), "w").write(head + "\n\n" + "\n".join(out) + "\n")
    print(head)
    for k, v in per.items():
        print("  %-40s %d" % (k, v))


if __name__ == "__main__":
    main()
