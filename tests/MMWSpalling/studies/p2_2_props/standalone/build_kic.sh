#!/bin/sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../../.." && pwd)
c++ -std=c++17 -O2 -I"$HERE" -I"$ROOT/src" -o "$HERE/../kic_dump" "$HERE/../kic_dump.cc"
echo "built: $HERE/../kic_dump"
