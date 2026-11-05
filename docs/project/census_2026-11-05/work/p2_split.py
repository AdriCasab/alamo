#!/usr/bin/env python3
"""Brace-depth splitter for a C++ class-body file.
Emits every block opened at class-body depth (2) or inside a nested struct (3->4),
with signature start line, end line, comment/code/blank counts, and the
immediately-preceding comment block size.
Usage: split.py FILE [base_depth]"""
import sys, re
path = sys.argv[1]
base = int(sys.argv[2]) if len(sys.argv) > 2 else 2
src = open(path).read().split('\n')
# Strip comments/strings, produce per-line "code" text and a per-line flag of comment-only.
code_lines = []; is_comment = []; is_blank = []
in_block = False
for line in src:
    out = []; i = 0; n = len(line); had_code = False; had_comment = in_block
    while i < n:
        c = line[i]
        if in_block:
            j = line.find('*/', i)
            if j < 0: i = n; break
            in_block = False; i = j + 2; continue
        if line.startswith('//', i): had_comment = True; break
        if line.startswith('/*', i): in_block = True; had_comment = True; i += 2; continue
        if c == '"' or c == "'":
            q = c; i += 1
            while i < n and line[i] != q:
                if line[i] == '\\': i += 1
                i += 1
            i += 1; out.append('""'); had_code = True; continue
        if not c.isspace(): had_code = True
        out.append(c); i += 1
    s = ''.join(out)
    code_lines.append(s)
    is_blank.append(not had_code and not had_comment and line.strip() == '')
    is_comment.append(had_comment and not had_code)
depth = 0
items = []  # (open_line, close_line, depth_at_open)
stack = []
for ln, s in enumerate(code_lines):
    for c in s:
        if c == '{':
            stack.append((ln, depth)); depth += 1
        elif c == '}':
            depth -= 1
            o, d = stack.pop()
            if d in (base, base + 1):
                items.append((o, ln, d))
items.sort()
def sig_start(o, d):
    # walk back from open line to find the start of the declaration:
    # stop after a line ending with ; or } or { (at depth d) or a preprocessor/access label, skipping comments
    j = o
    # if '{' line itself contains the signature text, start there
    k = o - 1
    while k >= 0:
        s = code_lines[k].strip()
        if is_comment[k] or is_blank[k] or s == '':
            break
        if s.endswith(';') or s.endswith('}') or s.endswith('{') or s.endswith(':') and s in ('public:', 'private:', 'protected:') or s.startswith('#'):
            break
        j = k; k -= 1
    return j
def lead_comment(j):
    k = j - 1; c = 0
    while k >= 0 and (is_comment[k] or (is_blank[k] and k-1 >= 0 and is_comment[k-1])):
        if is_comment[k]: c += 1
        k -= 1
    return c, k + 1
for o, e, d in items:
    j = sig_start(o, d)
    name_txt = ' '.join(code_lines[x].strip() for x in range(j, o + 1))
    m = re.findall(r'([~A-Za-z_][A-Za-z0-9_:]*)\s*\(', name_txt)
    kw = {'if','for','while','switch','return','sizeof','decltype','static_cast','AMREX_GPU_HOST_DEVICE','AMREX_FORCE_INLINE'}
    nm = next((x for x in m if x not in kw), None)
    if nm is None:
        mm = re.search(r'(struct|class|enum(?: class)?)\s+(\w+)', name_txt)
        nm = (mm.group(1) + ' ' + mm.group(2)) if mm else name_txt[:40]
    nc = sum(is_comment[x] for x in range(j, e + 1))
    nb = sum(is_blank[x] for x in range(j, e + 1))
    ncode = (e - j + 1) - nc - nb
    lc, ls = lead_comment(j)
    print(f"{d}\t{nm}\t{j+1}\t{e+1}\t{e-j+1}\t{ncode}\t{nc}\t{nb}\t{lc}\t{ls+1}")
