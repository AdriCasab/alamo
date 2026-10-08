#!/bin/bash
# P1g-2b deliverable 3: CTRL3 -- banking OFF on the SHIPPING binary.
# The only thing that can attribute the Leg B -> CTRL2 bore difference, because
# it differs from CTRL2 by one key on one binary. See run_leg.sh for the argv
# provenance and the before/after hash guard.
exec bash "$(dirname "$0")/run_leg.sh" CTRL3
