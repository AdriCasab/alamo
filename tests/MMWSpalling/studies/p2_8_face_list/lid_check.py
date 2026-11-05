#!/usr/bin/env python3
"""P2-8 (a1) -- the lid under face-centre evaluation, at dz 2 / 1 / 0.5 mm.

Every leg's h_expr starts `if((r>0.02)and(s<0.0633),0.001,...)`. Under the face
list each face is heated at its OWN centre (F2), s_f = z_n - z_f, so:
  * the original-surface +z face (z_f = rock top)            must be LIDDED   (s_f < 0.0633)
  * the shallowest side-face centre (z_f = rock top - dz/2)  must NOT be      (s_f >= 0.0633)
with a margin >= 0.1 mm on every mesh. The base lid `s<0.0637` is shown too:
it TIES at 1 mm (side-face centre exactly 0.0637) and lids side faces at 0.5 mm.

Two evaluations per row: exact decimal (Fraction) for the margin, and IEEE
double built the way the binary builds it (z_face = PLO + (k+1)*dz, a lateral
face centre = PLO + (k+0.5)*dz, s = z_n - z_f) for the branch it actually takes.
"""
from fractions import Fraction as F

Z_N = "0.2132"          # surface_patch.nozzle_z0 on every leg
HI = "0.15"             # geometry.prob_hi z
LID_NEW, LID_OLD = "0.0633", "0.0637"
MESHES = (("L0c", 75), ("L0", 150), ("L0f", 300))   # amr.n_cell z on the 150 mm block


def main():
    print("P2-8 (a1) -- the lid at the face centre")
    print("=" * 96)
    ok = True
    for leg, nz in MESHES:
        dz_q = F(HI) / nz
        dz = float(HI) / nz
        k_top = nz - 1
        # exact
        s_top_q = F(Z_N) - (k_top + 1) * dz_q
        s_side_q = F(Z_N) - (F(2 * k_top + 1, 2)) * dz_q
        # double, as the binary forms it
        s_top = float(Z_N) - (0.0 + float(k_top + 1) * dz)
        s_side = float(Z_N) - (0.0 + (float(k_top) + 0.5) * dz)
        print(f"  {leg:4s} dz {dz*1e3:.2f} mm | +z face s {float(s_top_q):.6f} | shallowest side-face"
              f" centre s {float(s_side_q):.6f}")
        for tag, T in (("NEW s<0.0633", LID_NEW), ("old s<0.0637", LID_OLD)):
            Tq, Td = F(T), float(T)
            top_lid = s_top < Td
            side_lid = s_side < Td
            m_top = (Tq - s_top_q) * 1000          # mm, >0 means lidded with that margin
            m_side = (s_side_q - Tq) * 1000        # mm, >0 means NOT lidded with that margin
            good = top_lid and not side_lid and m_top >= F(1, 10) and m_side >= F(1, 10)
            print(f"       {tag}: +z lidded {str(top_lid):5s} (margin {float(m_top):+.3f} mm) |"
                  f" side lidded {str(side_lid):5s} (margin {float(m_side):+.3f} mm)"
                  f"  -> {'ok' if good else 'FAILS'}")
            if tag.startswith("NEW"):
                ok &= good
    print("=" * 96)
    print(f"  a1 {'PASS' if ok else 'FAIL'} -- with s<0.0633 the original-surface +z face is lidded and the"
          f" shallowest side face is heated, margin >= 0.1 mm on every mesh"
          if ok else "  a1 FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
