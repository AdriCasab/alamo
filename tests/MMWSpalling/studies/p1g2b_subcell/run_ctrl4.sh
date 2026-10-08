#!/bin/bash
# P1g-2b deliverable 3: CTRL4 -- banking PLUS the thermal half.
# P1g-2's row shift measured at Meier scale for the FIRST time, two steps after
# it was implemented. Carries subcell=1 as well, because the thermal shift is
# reached only through `bank_on && sp_lateral_thermal` and would otherwise be
# inert. See run_leg.sh for the argv provenance and the before/after hash guard.
exec bash "$(dirname "$0")/run_leg.sh" CTRL4
