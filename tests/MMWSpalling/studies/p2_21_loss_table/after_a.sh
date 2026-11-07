#!/bin/bash
# waits for chain A (pid 39504) to exit, then runs chain B only if chain A logged its DONE line
cd /Users/tzetze20/amr_tools/alamo
while kill -0 39504 2>/dev/null; do sleep 20; done
grep -q "CHAIN A DONE" tests/MMWSpalling/studies/p2_21_loss_table/CHAIN.log && exec bash tests/MMWSpalling/studies/p2_21_loss_table/chain_b.sh
echo "$(date +%T) chain B NOT started: chain A did not finish clean" >> tests/MMWSpalling/studies/p2_21_loss_table/CHAIN.log
