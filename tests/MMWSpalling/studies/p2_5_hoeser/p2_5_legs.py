"""P2-5 rev 2 -- the ONE leg table, imported by both run_p2_5.py and
derive_lab_map.py.

⚠ Shared on purpose. The harness needs mesh / stop_time / K / nozzle_z0 to build
the argv and the derivation needs shape / K / nozzle_z0 / dz to build the h_expr
string and the lid constant. Two copies of this table would be free to drift,
and a lid constant that disagrees with the mesh it is used at is exactly the
rev-1 defect (3).
"""

# Her geometry and flame, from ACTIVE_STEP.md Sources.
D_NOZ = 0.0079          # nozzle D1 [m]
S0 = 0.0632             # L = 8, the calibration standoff [m]
Z_ROCK = 0.10           # the deck's rock top [m]
T_NOZ = 1866.0          # her nozzle-exit temperature [K]
T_ENT = 293.15
H0 = 7100.0             # her h at L = 8, 38 kW [W/m2K]
R_DISC = 0.0045         # her Ø 9 mm sensing element's radius [m]
R_OPEN = 0.020          # the flange's Ø 40 mm opening, radius [m]
T_FACE_REF = 850.0      # the reference face temperature for the design-time q
LID_H = 0.001           # ⚠ NOT zero: BuildFlameColumns aborts on a
                        #   non-positive in-patch h (rev 1 takeaway 7).
FIT_L = (10, 11, 13, 14)

# sigma in units of D_NOZ; None = the plateau shape A.
SHAPES = {"A": None, "G10": 1.00, "G075": 0.75, "G05": 0.50, "G035": 0.35}

# leg -> (dz [m], stop_time [s], jet_core_length, nozzle_z0 [m], shape)
LEGS = {
    "SMOKE":    (0.002,   2.0, 30.0, 0.1632, "G05"),
    "A":        (0.001,  30.0, 30.0, 0.1632, "A"),
    "G05":      (0.001,  30.0, 30.0, 0.1632, "G05"),
    "G075":     (0.001,  30.0, 30.0, 0.1632, "G075"),
    "G10":      (0.001,  30.0, 30.0, 0.1632, "G10"),
    "G05K20":   (0.001,  30.0, 20.0, 0.1632, "G05"),
    "G05M2":    (0.002,  30.0, 30.0, 0.1632, "G05"),
    "G05T120":  (0.001, 120.0, 30.0, 0.1632, "G05"),
    "G05L13":   (0.001,  30.0, 30.0, 0.2027, "G05"),
    "G035":     (0.001,  30.0, 30.0, 0.1632, "G035"),
}
# The packet's run order (Goal 1's table order); G035 only on the c2 rule.
ORDER = ["SMOKE", "A", "G05", "G075", "G10", "G05K20", "G05M2", "G05T120",
         "G05L13", "G035"]


def s_lid(leg):
    """The lid lifts once a column has receded HALF a cell below the original
    surface: mesh-aware, and half a cell clear of any floating-point tie.
    Rev 1 defect (3): the packet's old constant was s0 + 2 mm, i.e. two cells at
    1 mm and an exact tie (0.1632 - 0.098) at 2 mm."""
    dz, _, _, z0, _ = LEGS[leg]
    return (z0 - Z_ROCK) + dz / 2.0


def s_leg(leg):
    """The leg's own nozzle-to-original-surface standoff."""
    return LEGS[leg][3] - Z_ROCK
