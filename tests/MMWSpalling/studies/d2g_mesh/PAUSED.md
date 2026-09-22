# L09_1mm suspended (user request), 2026-09-19 19:44

- **Suspended with SIGSTOP** at t ≈ 103 s, on its way to the bottom stop at
  ~550 s. Not killed: it resumes bit-identically, so this is not a restart.
- **L09_2mm had already finished** (bottom watchdog).
- **PIDs** (the prterun launcher + 4 ranks): 41458 41460 41461 41463 41465 
- **Resume:** `kill -CONT $(cat tests/MMWSpalling/studies/d2g_mesh/paused_pids.txt)`
  (about 2.5 h left at the clean rate).
- **The pause counts in wall_s:** subtract it.
- **If the process dies** (reboot, sleep kill), rerun from t = 0 with
  `run.py --cases L09_1mm --force`.
- **Resumed** 2026-09-19 23:38 at t ≈ 107.808 s.
