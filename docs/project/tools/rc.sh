#!/bin/zsh
# Remote Control launcher for the ALAMO MMW project.
#
# Starts a Claude Code Remote Control SERVER rooted in this repo, so the Claude
# iPhone app (Code tab) or claude.ai/code can open sessions that run on this Mac
# with the full local toolchain (make, mpirun, the .venv python, /plan,
# /implement, /verify).
#
# Usage (from anywhere):
#   docs/project/tools/rc.sh                  # new server, session named "ALAMO"
#   docs/project/tools/rc.sh -c               # bring back the last server's session (<~4 h)
#   docs/project/tools/rc.sh --name "ALAMO P2-3"
#   any other `claude remote-control` flag is passed through (see --help).
#
# First run only:
#   1. answer `y` to "Enable Remote Control? (y/n)"
#   2. answer `y` to "Trust /Users/tzetze20/amr_tools/alamo? [y/N]"
#   3. in any claude session run /config and switch on
#      "Push when Claude decides" and "Push when actions required".
#
# While the server runs: SPACE shows a QR code for the phone, `w` toggles
# same-dir/worktree spawning (keep same-dir: bin/ and obj/ are built in place),
# Ctrl-C stops it. Stopped sessions come back with `rc.sh -c` for ~4 hours.
#
# Runs under `caffeinate -dis` (project rule: long runs must not lose the host
# to idle sleep). Lid-close sleep is NOT prevented: keep the lid open and the
# Mac on power while a phase runs from the phone.
set -eu

REPO="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO"

want_name=1
for a in "$@"; do
  case "$a" in
    --name|--name=*|-c|--continue|--session-id|--session-id=*) want_name=0 ;;
  esac
done
if (( want_name )); then
  set -- --name "ALAMO" "$@"
fi

echo "Remote Control server for $REPO"
echo "  flags: $*"
echo "  SPACE = QR code for the phone, Ctrl-C = stop (resume with rc.sh -c)"
exec caffeinate -dis claude remote-control "$@"
