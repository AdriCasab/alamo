#!/usr/bin/env python3
"""Assert the invariants of CLAUDE.md, the three phase specs and the roadmap head.

Every check here exists because the thing it checks was broken once, by hand, and found
only because someone asked twice. Run after any edit to those files.

  python3 docs/project/tools/spec_check.py
"""
import re, sys, pathlib

R = pathlib.Path(__file__).resolve().parents[3]
F = {n: (R / p).read_text() for n, p in {
    "plan": ".claude/commands/plan.md", "impl": ".claude/commands/implement.md",
    "ver": ".claude/commands/verify.md", "claude": "CLAUDE.md", "road": "ROADMAP.md"}.items()}
fail = []


def chk(cond, msg):
    if not cond:
        fail.append(msg)


# 1. exactly one spec CLAIMS the ledger-against-code read; the others must defer to it
chk("only outside read" in F["impl"], "implement.md no longer claims the outside read")
chk("Do not duplicate step zero" in F["plan"], "plan.md no longer defers to step zero")
chk("NOT here" in F["ver"], "verify.md no longer says where the pre-impl check lives")
for s in ("plan", "ver"):
    chk("only outside read" not in F[s], f"{s}.md now also claims the outside read")

# 2. the Criteria table's column ownership is stated once and consistently
chk("planner writes **six**" in F["plan"], "plan.md no longer says the planner writes six columns")
chk("no others" in F["plan"] and "append" in F["plan"],
    "plan.md must say six-and-no-others AND that later phases append")

# 3. structure: fences balanced, numbered runs contiguous, no path wrapped in a code span
for n in ("plan", "impl", "ver", "claude"):
    L = F[n].split("\n")
    chk(sum(1 for l in L if l.startswith("```")) % 2 == 0, f"{n}: odd number of code fences")
    chk(not any(re.search(r"`[^`]*/$", l) for l in L), f"{n}: a path is wrapped inside a code span")
    nums = [int(m.group(1)) for l in L if (m := re.match(r"^(\d+)\. ", l))]
    runs, prev = [], None
    for x in nums:
        if prev is None or x != prev + 1:
            runs.append([x])
        else:
            runs[-1].append(x)
        prev = x
    chk(all(r == list(range(1, len(r) + 1)) for r in runs), f"{n}: numbered list not contiguous from 1")

# 4. roadmap: standing rules capped, one Active step label, cap debt acknowledged
L = F["road"].split("\n")
s = L.index("## Standing rules")
e = next(i for i, l in enumerate(L) if i > s and l.startswith("## "))
n_rules = sum(1 for i in range(s, e) if re.match(r"^\d+\. ", L[i]))
chk(n_rules <= 15, f"standing rules: {n_rules} > cap of 15")
chk(F["road"].count("\nActive step:\n") == 1, "roadmap: not exactly one 'Active step:' label")
act = L.index("Active step:")
mil = next(i for i, l in enumerate(L) if l.startswith("Next milestones"))
n_act = sum(1 for i in range(act, mil) if L[i].startswith("- "))
chk(n_act == 1, f"'Active step:' holds {n_act} entries, must hold exactly 1 "
                "(this is the check the label-only version missed)")
chk("over its own cap, deliberately" in F["road"],
    "roadmap exceeds its cap without the dated exemption note")

# 5. the generators referenced by the specs must exist
for t in ("archive_index.py", "roadmap_status.py"):
    chk((R / "docs/project/tools" / t).exists(), f"missing generator docs/project/tools/{t}")
    chk(t in F["plan"], f"plan.md does not name {t}")

print(f"spec_check: {5} invariant groups, {len(fail)} failures")
for m in fail:
    print("  FAIL", m)
sys.exit(1 if fail else 0)
