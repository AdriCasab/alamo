# D2i paused

SIGSTOP at 2026-09-20 19:39:55 (suspend, not kill: the runs are not checkpointed, so a kill would lose them;
SIGCONT resumes bit-identically). `wall_s` in the .done files is inflated by the pause,
so computing time is quoted separately in RESULTS.md.

- LI0_1mm: pids [62989, 62992, 62997, 62999, 63002]
- I20_1mm: pids [62990, 62994, 62995, 62998, 63001]

Resumed 2026-09-20 21:36:32 by SIGCONT.
