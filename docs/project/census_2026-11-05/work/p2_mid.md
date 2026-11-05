
## Table 2. Non-function blocks (members, includes), runs of 3 lines or more

These are the lines between functions: data members, nested enums, the two `#include`s of
`Removal.H`/`FaceRemoval.H` (`:6621-6622`), and the file header. With Table 1 they cover
all 9,375 lines (`work/lineage_by_line.tsv` has the per-line assignment). Members are still
interleaved with functions (B10): face-list members sit at `:3286-3309`, between
`FaceEntry` and `FrPathOffsets`.

