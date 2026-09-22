# Note for the D2f planner: the burner support rule (user, 2026-09-18)

Carried into the D2e completion notes. The user's correction after rereading
Meier Ch. 8 (printed pages):
- p. 209: the burner hangs on coil tubing kept in tension and progresses by
  gravity;
- p. 212: the feet keep the outlet unblocked and direct the drill;
- p. 219–220: the drawworks impose axial jolts "to keep the drill free";
- p. 221: the hook load shows whether the burner hangs or lies on the rock;
- p. 223: the clearance is set by the advance rate.

**The burner rests on the highest rock touching its feet** (r 28–40 mm).

- **The model's 0.9-quantile pad rule** (D2b) was a guard against one
  never-fired 2 mm column freezing the burner. It lets about 10 % of the
  columns stand inside the feet.
- **The `foot_body_clearance` clip** (D2c) removed that rock. It has no
  independent physical basis, and D2e P-a identifies it as the mesh seed.
- **D2f should use the physical rule:**
  - a per-pad maximum;
  - shielded rock under the feet, with no gas from above;
  - an explicit undercut or operator-push rule;
  - convergence shown with the 2 mm / 1 mm pair method.
- **Caveat:** the maximum is set by the single slowest column, so its mesh
  behaviour must be tested, not assumed.
