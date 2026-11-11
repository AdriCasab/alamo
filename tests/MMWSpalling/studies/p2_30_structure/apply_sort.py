#!/usr/bin/env python3
"""P2-30 stage 2 (deliverable 4): apply S/COMMENT_SORT.txt.
MOVE / SPLIT lines leave src/ and go VERBATIM to docs/project/CODE_HISTORY.md under their slice and member; when a whole
block leaves and both neighbours are blank, the second blank goes too (no double blank is left behind).
Reads the candidates from S/HISTORY_CANDIDATES.txt; refuses unless every file is still the src_pre_2 pre-image, so the
line numbers are the ones sorted. Also writes S/COMMENT_MOVED.txt, the moved lines one per row (<id> <file> <text>).
Usage: apply_sort.py check      only verify the sort is complete and print the counts and the flag list
       apply_sort.py apply      write src/ and CODE_HISTORY.md"""
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "unit"))
import _src  # noqa: E402

HIST = os.path.join(ROOT, "docs", "project", "CODE_HISTORY.md")
PRE = os.path.join(HERE, "src_pre_2", "Integrator")


def candidates():
    out, cur = {}, None
    for ln in open(os.path.join(HERE, "HISTORY_CANDIDATES.txt")):
        m = re.match(r"(C\d{3})  (\S+?):(\d+)-(\d+)  member=(\S+.*?)  lines=(\d+)", ln)
        if m:
            cur = m.group(1)
            out[cur] = dict(id=cur, file=m.group(2), a=int(m.group(3)), b=int(m.group(4)), member=m.group(5))
    return out


def sort_table():
    out = {}
    for ln in open(os.path.join(HERE, "COMMENT_SORT.txt")):
        if ln.startswith("#") or not ln.strip():
            continue
        head, reason = ln.rstrip("\n").split(" | ", 1)
        f = head.split()
        out[f[0]] = dict(verdict=f[1], lines=f[2] if len(f) > 2 else "", reason=reason.strip())
    return out


def rel_of(name):
    return _src.PARENT if name == "MMWSpalling.H" else _src.SLICE_PREFIX + name


def parse_lines(spec, n):
    out = set()
    for part in spec.split(","):
        a, _, b = part.partition("-")
        out |= set(range(int(a), int(b or a) + 1))
    if not out or max(out) > n:
        sys.exit("bad line spec %s (block of %d lines)" % (spec, n))
    return out


def plan():
    cand, tab = candidates(), sort_table()
    miss = sorted(set(cand) - set(tab))
    extra = sorted(set(tab) - set(cand))
    noreason = [k for k, v in tab.items() if not v["reason"]]
    bad = [k for k, v in tab.items() if v["verdict"] not in ("MOVE", "KEEP", "SPLIT")]
    if miss or extra or noreason or bad:
        sys.exit("sort incomplete: missing %s, extra %s, no reason %s, bad verdict %s" % (miss, extra, noreason, bad))
    src = _src.tree()
    moved = []      # (id, rel, line_no (1-based, in the slice), text)
    for cid, c in sorted(cand.items()):
        v = tab[cid]
        rel = rel_of(c["file"])
        L = src[rel].split("\n")
        n = c["b"] - c["a"] + 1
        if not all(L[i - 1].lstrip().startswith("//") for i in range(c["a"], c["b"] + 1)):
            sys.exit("%s: %s:%d-%d is not a comment block in the live file" % (cid, c["file"], c["a"], c["b"]))
        if v["verdict"] == "MOVE":
            which = set(range(1, n + 1))
        elif v["verdict"] == "SPLIT":
            which = parse_lines(v["lines"], n)
        else:
            which = set()
        for k in sorted(which):
            ln = c["a"] + k - 1
            moved.append((cid, rel, ln, L[ln - 1]))
    return cand, tab, moved


def flags(moved):
    """moved lines that name an identifier present in non-comment src/ beside must / never / NOT / the warning sign"""
    code = re.sub(r"//[^\n]*", "", _src.text() + "\n" + "\n".join(t for r, t in _src.tree().items() if r.startswith("Numeric")))
    ids = set(re.findall(r"[A-Za-z_]\w+", code))
    out = []
    for cid, rel, ln, t in moved:
        if not re.search(r"\bmust\b|\bnever\b|\bNOT\b|⚠", t):
            continue
        names = sorted(set(re.findall(r"[A-Za-z_]\w+", t)) & ids - {"must", "never", "NOT"})
        names = [n for n in names if n not in ("the", "a", "is", "of", "and", "it", "in", "to")]
        out.append((cid, rel, ln, t.strip(), names))
    return out


def main():
    cand, tab, moved = plan()
    per = collections.Counter(v["verdict"] for v in tab.values())
    print("candidates %d; MOVE %d, SPLIT %d, KEEP %d; lines moved %d of %d" % (
        len(cand), per["MOVE"], per["SPLIT"], per["KEEP"], len(moved), sum(c["b"] - c["a"] + 1 for c in cand.values())))
    fl = flags(moved)
    print("flag list (moved lines naming a live identifier beside must / never / NOT / the warning sign): %d" % len(fl))
    for cid, rel, ln, t, names in fl:
        print("  %s %s:%d  names %s\n      %s" % (cid, os.path.basename(rel), ln, names, t[:130]))
    if sys.argv[1] == "check":
        return 0
    # ------- apply --------
    src = _src.tree()
    for name in ["MMWSpalling.H"] + [os.path.join("MMWSpalling", os.path.basename(p)) for p in sorted(os.listdir(os.path.join(PRE, "MMWSpalling")))]:
        rel = "Integrator/" + name
        pre = open(os.path.join(PRE, name)).read()
        if src[rel] != pre:
            sys.exit("apply_sort: %s is not the src_pre_2 pre-image; refusing" % rel)
    by_file = collections.defaultdict(list)
    for cid, rel, ln, t in moved:
        by_file[rel].append(ln)
    for rel, lns in by_file.items():
        L = src[rel].split("\n")
        rm = set(lns)
        # a whole-block removal: if the lines on both sides are blank, drop the following blank as well
        for cid, c in cand.items():
            if rel_of(c["file"]) == rel and set(range(c["a"], c["b"] + 1)) <= rm:
                before = L[c["a"] - 2].strip() if c["a"] >= 2 else ""
                after = L[c["b"]].strip() if c["b"] < len(L) else ""
                if before == "" and after == "":
                    rm.add(c["b"] + 1)
        open(os.path.join(_src.SRC, rel), "w").write("\n".join(s for i, s in enumerate(L, 1) if i not in rm))
    # ---- the history file: under its slice and member, verbatim
    doc = ["# Code history: comments moved out of the MMWSpalling sources", "",
           "Moved **verbatim** by P2-30 (2026-11-10), per `tests/MMWSpalling/studies/p2_30_structure/COMMENT_SORT.txt`. "
           "Each block is listed under the slice it left and the member it stood at; the source no longer carries it. "
           "The packet that wrote a block is in its first words; the code it describes may be long gone.", ""]
    groups = collections.OrderedDict()
    for cid, rel, ln, t in moved:
        c = cand[cid]
        groups.setdefault(os.path.basename(rel), collections.OrderedDict()).setdefault(c["member"], collections.OrderedDict()).setdefault(cid, []).append(t)
    for f, members in groups.items():
        doc.append("## %s" % f)
        doc.append("")
        for m, blocks in members.items():
            doc.append("### %s" % m)
            doc.append("")
            for cid, lines in blocks.items():
                doc.append("```cpp")
                doc.extend(lines)
                doc.append("```")
                doc.append("")
    open(HIST, "w").write("\n".join(doc))
    with open(os.path.join(HERE, "COMMENT_MOVED.txt"), "w") as f:
        for cid, rel, ln, t in moved:
            f.write("%s\t%s\t%s\n" % (cid, os.path.basename(rel), t))
    print("applied: %d lines left src/; %s written" % (len(moved), os.path.relpath(HIST, ROOT)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
