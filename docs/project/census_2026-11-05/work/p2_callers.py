#!/usr/bin/env python3
"""For each depth-2 function in split_raw.tsv, list the functions (in any of the
three files) whose bodies call it. Body ranges from split_raw.tsv (MMWSpalling.H)
and a crude scan of Removal.H/FaceRemoval.H function heads."""
import re, sys, subprocess
base = 'src/Integrator/'
rows = [l.rstrip('\n').split('\t') for l in open('docs/project/census_2026-11-05/work/p2_split_raw.tsv')]
funcs = [(r[1], int(r[2]), int(r[3])) for r in rows if r[0] == '2']
main = open(base + 'MMWSpalling.H').read().split('\n')
def owner(ln):
    best = None
    for n, a, b in funcs:
        if a <= ln <= b: best = n
    return best or f'members@{ln}'
# included files: split with same script
def split_inc(path):
    out = subprocess.run(['python3', 'docs/project/census_2026-11-05/work/p2_split.py', path, '0'],
                         capture_output=True, text=True).stdout
    fs = []
    for l in out.splitlines():
        r = l.split('\t')
        if r[0] == '0': fs.append((r[1], int(r[2]), int(r[3])))
    return fs
inc = {p: (open(base + p).read().split('\n'), split_inc(base + p)) for p in
       ['MMWSpalling/Removal.H', 'MMWSpalling/FaceRemoval.H']}
names = sorted(set(n for n, a, b in funcs if not n.startswith(('struct', 'enum', 'MMWSpalling'))))
for nm in names:
    pat = re.compile(r'\b' + re.escape(nm) + r'\s*\(')
    hits = set()
    for i, line in enumerate(main, 1):
        s = line.split('//')[0]
        if pat.search(s):
            o = owner(i)
            if o != nm or not any(a == i for n, a, b in funcs if n == nm): hits.add(f'{o}')
    for p, (src, fs) in inc.items():
        for i, line in enumerate(src, 1):
            s = line.split('//')[0]
            if pat.search(s):
                o = next((n for n, a, b in fs if a <= i <= b), f'top@{i}')
                hits.add(f'{p.split("/")[1]}:{o}:{i}')
    hits.discard(nm)
    print(f'{nm}\t' + '; '.join(sorted(hits)))
