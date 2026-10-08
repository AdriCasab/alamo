# P0b rev7 -- MANIFEST

Binary under test: `1c22018e9ad8e8bc93bf948b1324df3d5839600a3ec9e9b96eaebf64536408cf`

⚠ The toolchain is NOT byte-reproducible: rebuilding from identical
sources gives a different sha256. So the only way a leg's binary sha can
match the final build is for NOTHING to have been rebuilt after it ran.
Every row below was produced after the final `make` and before any other.

## (a) (c') (c'') (x) -- the 1-D column matrix, 18 legs

harness: `tests/MMWSpalling/validation/p0_steady_column/test`

| leg | output path | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|---|
| `B3_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/B3_0.5` | `f2c8d6a5462654df` | 25060 | ok |
| `B3_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/B3_1.0` | `67a57f22fee47d57` | 6104 | ok |
| `B3_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/B3_2.0` | `38f599542d477028` | 1476 | ok |
| `B5_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/B5_0.5` | `d45bb47b6c60e630` | 19688 | ok |
| `B5_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/B5_1.0` | `6073de8726cc4b11` | 4312 | ok |
| `B5_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/B5_2.0` | `94b15be54333fe7e` | 1084 | ok |
| `C_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/C_0.5` | `e83c21270164a437` | 25384 | ok |
| `C_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/C_1.0` | `72f916fc54e9d68e` | 6564 | ok |
| `C_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/C_2.0` | `40b74961aae27839` | 1772 | ok |
| `D_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/D_0.5` | `5c4ce6da85fea466` | 17408 | ⚠ **STALE — predates the final binary** |
| `D_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/D_1.0` | `3b64f99ee07a28df` | 2240 | ⚠ **STALE — predates the final binary** |
| `D_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/D_2.0` | `58538e1c3840155c` | 288 | ⚠ **STALE — predates the final binary** |
| `F_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/F_0.5` | `d14fd039692ad604` | 21280 | ok |
| `F_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/F_1.0` | `34ee3fd6893e7e10` | 4536 | ok |
| `F_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/F_2.0` | `a59afdaa5e50dbe7` | 940 | ok |
| `P_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/P_0.5` | `0ecd5588f6cffd41` | 23312 | ok |
| `P_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/P_1.0` | `6cc8ac3d22b79d21` | 5980 | ok |
| `P_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/P_2.0` | `79798032f1c2c473` | 1680 | ok |
| `S_0.5` | `tests/MMWSpalling/validation/p0_steady_column/output/S_0.5` | `9d5537518cefe306` | 21464 | ok |
| `S_1.0` | `tests/MMWSpalling/validation/p0_steady_column/output/S_1.0` | `c155b2b982bc7126` | 5020 | ok |
| `S_2.0` | `tests/MMWSpalling/validation/p0_steady_column/output/S_2.0` | `5668d23970adb65b` | 1220 | ok |

argv: see `tests/MMWSpalling/validation/p0_steady_column/test`

## (d) -- the Kant anchor, 15 legs

harness: `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/test`

| leg | output path | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|---|
| `C150_0.5` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/C150_0.5` | `cc940660a36446f8` | 11520 | ok |
| `C150_1.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/C150_1.0` | `ac2f8cb0aec82f93` | 1152 | ok |
| `C150_2.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/C150_2.0` | `3a962fd618d313f5` | 144 | ok |
| `C500_0.5` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/C500_0.5` | `41ea71cae94e94ca` | 94464 | ok |
| `C500_1.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/C500_1.0` | `17293d42305cd5eb` | 12096 | ok |
| `C500_2.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/C500_2.0` | `ecfb17e9d75991e7` | 1584 | ok |
| `R150_0.5` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/R150_0.5` | `ea0ab581728eada9` | 13824 | ok |
| `R150_1.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/R150_1.0` | `0d76066cc56a200e` | 2304 | ok |
| `R150_2.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/R150_2.0` | `e0b83f6db0f1bc59` | 432 | ok |
| `R500_0.5` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/R500_0.5` | `037cebde5962d1f4` | 108288 | ok |
| `R500_1.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/R500_1.0` | `f636c998a4281cb8` | 13248 | ok |
| `R500_2.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/R500_2.0` | `340fb2df916ccd67` | 1584 | ok |
| `W150_0.5` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/W150_0.5` | `f457f6098039490d` | 11444 | ok |
| `W150_1.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/W150_1.0` | `8b08402ac5e737f0` | 2026 | ok |
| `W150_2.0` | `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output/W150_2.0` | `c70823778beff27b` | 431 | ok |

argv: see `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/test`

## (f) (h) (i) d2 d6 -- the Meier FH rerun, closure FOLLOW

harness: `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7/cmd.txt`

| leg | output path | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|---|
| `FH` | `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7/FH` | `-` | - | ok |
| `FH` | `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7/FH_removal_events.csv` | `28ca953426c870d7` | 200836 | ok |

argv: see `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7/cmd.txt`

## d6 -- the Meier control, closure CELL

harness: `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7_cell/cmd.txt`

| leg | output path | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|---|
| `FH` | `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7_cell/FH` | `-` | - | ok |
| `FH` | `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7_cell/FH_removal_events.csv` | `c92125020742dc6a` | 123030 | ok |

argv: see `tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7_cell/cmd.txt`

