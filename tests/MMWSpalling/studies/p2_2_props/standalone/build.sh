#!/bin/sh
# Compile b6_dump against the REAL src/Numeric/Material/Grimsel.H.
# The stub dir comes FIRST so Set/Set.H and IO/ParmParse.H resolve to the
# offline stand-ins, while Material.H and Grimsel.H come from src/.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../../.." && pwd)
c++ -std=c++17 -O2 -I"$HERE" -I"$ROOT/src" \
    -o "$HERE/../b6_dump" "$HERE/../b6_dump.cc"
echo "built: $HERE/../b6_dump"
