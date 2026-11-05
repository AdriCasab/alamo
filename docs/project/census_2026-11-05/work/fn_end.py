#!/usr/bin/env python3
"""Census helper (2026-11-05): print [start, end] of a brace-delimited member
function given the line of its signature. Naive brace count (ignores braces in
comments/strings, adequate for this code)."""
import sys, re
path = sys.argv[1]
lines = open(path).read().split("\n")
for a in sys.argv[2:]:
    s = int(a)
    depth, started = 0, False
    for n in range(s - 1, len(lines)):
        t = re.sub(r"//.*", "", lines[n])
        t = re.sub(r'"([^"\\]|\\.)*"', '""', t)
        for ch in t:
            if ch == "{": depth += 1; started = True
            elif ch == "}": depth -= 1
        if started and depth == 0:
            print(f"{s}-{n+1}: {lines[s-1].strip()[:80]}")
            break
