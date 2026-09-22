# Q10_1mm and Q09_1mm suspended (user request), 2026-09-18 18:37

- **Suspended with SIGSTOP** at t ≈ 98 s of 300, at the user's request.
  Not killed: they resume bit-identically, so this is not a restart.
- **PIDs** (the two prterun launchers + 8 ranks): 59902 59903 59904 59905 59906 59907 59908 59909 59910 59911 
- **Resume:** `kill -CONT $(cat tests/MMWSpalling/studies/d2f_support/paused_pids.txt)`
- **Remaining:** about 200 s of simulated time each, roughly 70–90 min with
  both running.
- **The pause is inside the runs' wall time** (the `.done` wall_s): subtract
  it when quoting cost.
- **If the machine reboots or the processes die, the runs are lost.** Rerun
  them from t = 0 with `run.py --cases Q10_1mm Q09_1mm --jobs 2 --force`.
  Never restart them.
- **Resumed** with SIGCONT on 2026-09-19 10:33, at t ≈ 99 s. Paused from 18:35 to 10:33, about 16 h, which must be subtracted from the wall_s of both runs.
