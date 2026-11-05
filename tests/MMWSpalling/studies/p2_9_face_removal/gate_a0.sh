#!/bin/bash
# P2-9 rev 2 (a0): the build sees header edits. Writes GATE_a0.txt.
#   1. obj/obj-3d-g++/mmwspalling.cc.d lists FaceRemoval.H AND Grimsel.H
#   2. touching ONLY FaceRemoval.H, then make, recompiles mmwspalling.cc
#      (the object's mtime moves) -- and the binary hash is unchanged, since
#      the source content did not change.
ROOT=/Users/tzetze20/amr_tools/alamo
OUT=$ROOT/tests/MMWSpalling/studies/p2_9_face_removal/GATE_a0.txt
cd "$ROOT" || exit 1
D=obj/obj-3d-g++/mmwspalling.cc.d
O=obj/obj-3d-g++/mmwspalling.cc.o
{
echo "P2-9 rev 2 (a0) -- the build sees header edits"
echo "=================================================================="
f1=$(tr ' \\' '\n\n' < $D | grep -c "MMWSpalling/FaceRemoval.H")
f2=$(tr ' \\' '\n\n' < $D | grep -c "Numeric/Material/Grimsel.H")
echo "  $D lists FaceRemoval.H: $f1   Grimsel.H: $f2"
echo "  $D first target line: $(head -1 $D | cut -c1-120)"
m0=$(stat -f %m $O); h0=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-24)
sleep 1
touch src/Integrator/MMWSpalling/FaceRemoval.H
EIGEN=$PWD/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8 > /tmp/p29_a0_make.log 2>&1
rc=$?
m1=$(stat -f %m $O); h1=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-24)
rec=$(sed 's/\x1b\[[0-9;]*m//g' /tmp/p29_a0_make.log | grep -c "COMPILING.*src/mmwspalling.cc")
echo "  touch FaceRemoval.H only -> make rc $rc; 'COMPILING src/mmwspalling.cc' lines: $rec"
echo "  mmwspalling.cc.o mtime $m0 -> $m1 ($([ $m1 -gt $m0 ] && echo RECOMPILED || echo NOT recompiled))"
echo "  binary sha256 $h0 -> $h1"
if [ $f1 -ge 1 ] && [ $f2 -ge 1 ] && [ $m1 -gt $m0 ] && [ $rc -eq 0 ]; then echo "  a0 PASS"; else echo "  a0 FAIL"; fi
} | tee "$OUT"
