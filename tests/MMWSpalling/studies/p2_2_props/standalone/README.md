# `standalone/` — compile the REAL `Grimsel.H` offline

⚠⚠ WHY THIS EXISTS. `b6_dump.cc` used to `#include "g_standalone.H"`, a
240-line hand-edited **copy** of `src/Numeric/Material/Grimsel.H` that lived
only in a session scratchpad. Two things were wrong with that, and the second
is the dangerous one:

 1. the repo could not build its own `(b2)` witness — the header was not in it;
 2. **the copy drifts.** It was written before P2-2b's dense stiffness table,
    so running it after that change would have verified the OLD arithmetic and
    printed a pass that meant nothing. Standing rule 8 says a number is backed
    only if both the file it was read from and the code that reduced it are in
    the repository, and a scratchpad copy is neither.

⇒ these stubs supply ONLY what `Grimsel.H` needs from the wider tree
(`Set::Scalar`, a no-op `IO::ParmParse`), and the law itself is `#include`d
from `src/` so it CANNOT drift. If `Grimsel.H` stops compiling against these
stubs, that is a real interface change and the shim should be fixed, never
the copy forked again.

    sh build.sh && ./b6_dump | head
