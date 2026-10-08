#!/usr/bin/env python3
"""P2-1: THE COMPOSITE WALL LAW -- one implementation, imported by everything.

⚠⚠ THIS FILE IS THE PYTHON TWIN OF `MMWSpalling.H::WallJetH()`. The offline
derivation (`derive_wall_h.py`), the scorer (`score_wallh.py`) and the C++
resolver must all produce the same number for the same inputs, and
`(b2)`/`(c2)` are only meaningful if they do. Change one, change both, and
re-run `derive_wall_h.py` -- its `CROSS-CHECK` block prints the C++ constants
it parsed out of the header and fails loudly when they drift apart.

The law, as the packet states it:

      h_wall(z) = h_wj(z)   for z <  z_fill        (wall jet, thin layer)
                = h_duct    for z >= z_fill        (fully-developed annulus)

  wall jet   v(r)     = 1.05 * u * D_e / r        radial wall-jet max velocity
             delta(r) = 0.09 * r                  layer thickness at the corner
             d delta/dz = 0.073                   plane-jet growth up the wall
             v(z)     = v_c * sqrt(delta_c/delta(z))   plane-jet momentum
             c_f      = 0.0315 * Re_delta^-0.182  WITH the 1/2 convention
             St*Pr^(2/3) = c_f/2                  Chilton-Colburn
             h_wj     = St * rho * c_p * v
  z_fill     delta(z_fill) == (r_bore - R_SHIELD), the wall-to-shield gap

⚠⚠ `u * D_e` IS NEVER SPLIT. Writing `u` alone would mean quoting a nozzle
velocity of ~1840 m/s, which is supersonic at 1900 K and is an artefact of
reading a momentum-equivalent diameter against a geometric exit area. The
radial wall jet depends on the jet's MOMENTUM, and the momentum-equivalent
pair satisfies rho_inf*u^2*D_e^2 = rho_j*U_j^2*D^2, so

      u * D_e = U_j * D * sqrt(rho_j/rho_inf) = 4*mdot / (pi * D * sqrt(rho_j*rho_inf))

which is convention-free and needs no velocity at all. Carried as one product.

⚠ THE HALF-FACTOR IS SETTLED BY MEASUREMENT, NOT BY READING A DEFINITION --
see the packet. `c_f/2` reproduces the measured friction-velocity band
u_tau/U_max = 0.050-0.064; the other reading needs 0.083, which no measurement
shows. `cf_half=False` is provided ONLY so the derivation can print what the
excluded convention would have given.
"""
import math

# ---------------------------------------------------------------------------
# THE LAW'S OWN CONSTANTS. Every one is quoted in the packet's `Sources`
# section; none is fitted here and none may be adjusted to move a criterion.
# ---------------------------------------------------------------------------
C_WJ_VEL   = 1.05      # v_max = C * u * D_e / r, radial wall jet
C_WJ_DELTA = 0.09      # delta = C * r at the corner
C_WJ_GROW  = 0.073     # d delta / dz, plane jet climbing the wall
C_CF_A     = 0.0315    # c_f = A * Re_delta^B
C_CF_B     = -0.182

P_ATM = 101325.0       # Pa; the model runs at atmospheric pressure
R_GAS = 287.0          # J/(kg K); combustion products ~ air to within a few %

# The burner's water-cooled shield is Oe67 mm => 33.5 mm radius (VT5 drawings).
# ⚠ NOT side_face_h_r_min = 40 mm: that gates WHICH faces take wall treatment,
# this is what the annulus GAP is measured against. Same constant, same
# provenance and same warning as p2_0_annulus/annulus_h.py:R_SHIELD.
R_SHIELD = 0.0335


def rho_gas(T):
    """Ideal-gas density at atmospheric pressure [kg/m3]."""
    return P_ATM / (R_GAS * T)


# ---------------------------------------------------------------------------
# THE HEIGHT DATUM -- deliverable 2. ONE definition, here and in C++.
# ---------------------------------------------------------------------------
# ⚠⚠ Python twin of `MMWSpalling.H::BuildWallCornerColumns`. The two must agree
# cell for cell; `derive_wall_h.py` prints the rule it used and the scorer
# re-derives it from the shipped run, so a drift shows up as a (c2) miss.
#
# THE RULE, stated once:
#
#   The corner beneath a wall face is the DEEPEST column top on the straight
#   ray from that face's column inward to the axis, truncated at the wall
#   band's own inner edge r = side_face_h_r_min.
#
#       z = (k - k_corner[col] - 0.5) * dz,  floored at dz/2
#
# ⚠ WHY NOT THE NEAREST OPEN NEIGHBOUR'S TOP. That reads the one-cell STAIRCASE
# tread, not the hole's floor: it gives a median height of 7 mm on the Meier
# wall where the packet's own floor profile puts the band at 30-65 mm. The
# wall-jet layer is 4.4 mm at the corner and ~15 mm by z = 150 mm, i.e. it spans
# several treads, so a 2 mm-wide ledge every 10 mm is ROUGHNESS and does not
# restart the jet. Resetting z at each tread would also make the datum
# mesh-dependent, which is the D2j-0b failure mode.
#
# ⚠ WHY NOT THE AXIS FRONT. That is the quantity prior packets quoted (p50
# 90 mm) and it ignores that the jet turns at the corner under THIS azimuth,
# not under the axis -- the crater floor rises ~50 mm from the axis to the band.
#
# ⚠ WHY TRUNCATE AT r_gate. Inside it the code is already in regime 7b (floor
# riser, the owning column's own profile); r_gate is therefore the existing,
# unchanged boundary between "floor" and "wall", and reusing it means the datum
# introduces no new radius.
def wall_corner_columns(ktop, rr, r_gate, missing):
    """k_corner[i][j] for every column in the band. `ktop` and the result are
    integer k indices; `missing` marks a column with no solid cell."""
    nx = len(ktop)
    ny = len(ktop[0])
    out = [[missing] * ny for _ in range(nx)]
    for i in range(nx):
        for j in range(ny):
            if rr[i][j] < r_gate or ktop[i][j] == missing:
                continue
            best = ktop[i][j]
            n = max(i, j)
            for s in range(1, n + 1):
                t = 1.0 - s / float(n + 1)
                ii = int(round((i + 0.5) * t - 0.5))
                jj = int(round((j + 0.5) * t - 0.5))
                if ii < 0 or jj < 0 or rr[ii][jj] < r_gate:
                    break
                if ktop[ii][jj] != missing and ktop[ii][jj] < best:
                    best = ktop[ii][jj]
            out[i][j] = best
    return out


def u_De(mdot_full, D_nozzle, T_nozzle, T_gas):
    """The momentum product u*D_e [m2/s]. See the module docstring."""
    rj, ri = rho_gas(T_nozzle), rho_gas(T_gas)
    return 4.0 * mdot_full / (math.pi * D_nozzle * math.sqrt(rj * ri))


def h_walljet(z, r, uDe, T_gas, cp, mu, Pr, cf_half=True):
    """The wall-jet branch at height `z` above the corner, at radius `r`.

    `z` is the height of the face above the floor AT THAT FACE'S OWN RADIUS --
    the datum the packet defines, not the height above the axis front.
    """
    if not (r > 0.0):
        return float("nan")
    rho = rho_gas(T_gas)
    nu = mu / rho
    v_c = C_WJ_VEL * uDe / r
    d_c = C_WJ_DELTA * r
    d_z = d_c + C_WJ_GROW * max(z, 0.0)
    v = v_c * math.sqrt(d_c / d_z)            # plane-jet momentum conservation
    Re = v * d_z / nu
    if not (Re > 0.0):
        return float("nan")
    cf = C_CF_A * Re ** C_CF_B
    St = (cf / (2.0 if cf_half else 1.0)) / Pr ** (2.0 / 3.0)
    return St * rho * cp * v


def z_fill(r):
    """Height at which the layer fills the wall-to-shield gap [m].

    Below it the thin-layer model has a free stream to entrain; at it the
    layer owns the whole gap and the flow IS the duct. Returns 0.0 when the
    gap is already thinner than the corner layer (no wall-jet zone at all).
    """
    gap = r - R_SHIELD
    d_c = C_WJ_DELTA * r
    if gap <= d_c:
        return 0.0
    return (gap - d_c) / C_WJ_GROW


def h_wall(z, r, uDe, T_gas, cp, mu, Pr, h_duct, cf_half=True):
    """The composite law. ⚠ A HARD SWITCH at z_fill, never max(h_wj, h_duct):
    under the `max` the duct is never selected and the deep wall keeps a
    wall-jet value from a model that has no free stream left to entrain."""
    if z >= z_fill(r):
        return h_duct
    return h_walljet(z, r, uDe, T_gas, cp, mu, Pr, cf_half)
