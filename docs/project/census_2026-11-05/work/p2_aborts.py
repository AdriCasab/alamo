#!/usr/bin/env python3
"""List every Util::Abort in a file with its guarding condition and message head."""
import sys, re
src = open(sys.argv[1]).read().split('\n')
for n, line in enumerate(src):
    if 'Util::Abort' not in line: continue
    # condition: walk back to the nearest line starting with if/else if/for (max 8)
    cond = ''
    for k in range(n, max(-1, n - 9), -1):
        s = src[k].strip()
        if re.match(r'(else\s+)?if\s*\(|.*\bif\s*\(|for\s*\(', s) and k != n or (k == n and re.match(r'.*\bif\s*\(', s)):
            # gather the full condition until a line that ends the condition
            j = k; buf = []
            while j <= n:
                buf.append(src[j].strip())
                if 'Util::Abort' in src[j]: break
                j += 1
            cond = ' '.join(buf)
            break
    msg = ' '.join(x.strip() for x in src[n:n+3])
    m = re.search(r'Util::Abort\(INFO,\s*(.*)', msg)
    print(f"{n+1}\t{cond[:260]}\t||\t{(m.group(1) if m else msg)[:200]}")
