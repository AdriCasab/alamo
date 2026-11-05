#!/usr/bin/env python3
"""Group the non-function lines (members, includes, gaps) of MMWSpalling.H into runs of one
lineage, from p2_lineage_by_line.tsv; runs of >= 3 lines are printed as a markdown table."""
W = 'docs/project/census_2026-11-05/work/'
L = {}
for l in open(W + 'p2_lineage_by_line.tsv'):
    a = l.rstrip('\n').split('\t'); L[int(a[0])] = (a[1], a[2] == '1', a[3])
infn = set()
for l in open(W + 'p2_split_raw.tsv'):
    r = l.rstrip('\n').split('\t')
    if r[0] == '2' and not r[1].startswith('enum'):
        infn |= set(range(int(r[9]), int(r[3]) + 1))
N = max(L); runs = []; cur = None
for x in range(1, N + 1):
    if x in infn:
        if cur: runs.append(cur); cur = None
        continue
    lab = (L[x][0], L[x][1])
    if cur and cur[2] == lab and cur[1] == x - 1: cur[1] = x
    else:
        if cur: runs.append(cur)
        cur = [x, x, lab]
if cur: runs.append(cur)
print('| lines | count | lineage | code/comment |')
print('|---|---|---|---|')
for a, b, (lin, cl) in runs:
    if b - a + 1 < 3: continue
    code = sum(1 for x in range(a, b + 1) if L[x][2] == 'code')
    com = sum(1 for x in range(a, b + 1) if L[x][2] == 'comment')
    print(f'| {a}-{b} | {b-a+1} | {lin}{" CLOSED" if cl else ""} | {code}/{com} |')
