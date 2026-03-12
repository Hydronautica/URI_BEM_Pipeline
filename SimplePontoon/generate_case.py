#!/usr/bin/env python3
"""
Generate Nemoh BEM case for a twin-pontoon (catamaran) geometry.

Geometry: Two semicylindrical pontoons with axes along x-direction,
half-submerged (cylinder center at waterline z=0).

Pontoon 1: y = -0.125
Pontoon 2: y = +0.125  (mirror of pontoon 1 about y=0)

Each pontoon: R=0.026m radius, L=0.75m length, extending from x=-0.375 to x=+0.375.

Adapted from Salome CAD/mesh script for use with Nemoh BEM solver.
"""

import math
import os

# ============================================================
# INPUT PARAMETERS
# ============================================================
R     = 0.025       # pontoon radius (m)
L     = 0.75        # pontoon length (m)
y_c1  = -0.125      # y-center of pontoon 1
y_c2  =  0.125      # y-center of pontoon 2
z_c   = 0.0         # z-center of cylinder (at waterline)
rho   = 1000.0      # water density (kg/m^3)
g     = 9.81        # gravity (m/s^2)
z_G   = -0.005      # vertical CG position (m, below waterline)
depth = 1.0         # water depth (0 = infinite)

# Mesh resolution
N_circ = 16         # circumferential divisions (semicircle)
N_x    = 10         # axial divisions along pontoon length
N_r    = 3         # radial divisions on end caps
N_lid  = 4          # transverse divisions on waterplane lid (for irreg freq removal)
use_lid = True      # add lid panels at z=0 for irregular frequency removal

# Frequency range
n_freq    = 120      # number of frequencies
omega_min = 0.2     # min frequency (rad/s)
omega_max = 40.0    # max frequency (rad/s)

# Wave directions
n_dir   = 1
dir_min = 0.0
dir_max = 0.0

# ============================================================
# CASE DIRECTORY
# ============================================================
case_dir = os.path.dirname(os.path.abspath(__file__))
mech_dir = os.path.join(case_dir, "Mechanics")
results_dir = os.path.join(case_dir, "results")
os.makedirs(mech_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# ============================================================
# MESH GENERATION
# ============================================================

def generate_pontoon_mesh(y_c, R, L, z_c, N_circ, N_x, N_r, node_offset=0,
                          add_lid=False, N_lid=4):
    """
    Generate mesh nodes and panels for one semicylindrical pontoon.

    The submerged surface consists of:
    - Hull: semicylindrical surface below waterline
    - Two end caps: semicircular flat faces at x = -L/2 and x = +L/2
    - (Optional) Lid panels at z=0 for irregular frequency removal

    Panel normals point outward (into the fluid).
    Lid panel normals point upward (+z), into the fluid from inside the body.

    Returns:
        nodes: list of (x, y, z)
        panels: list of (n1, n2, n3, n4) with 1-based indices including offset
    """
    nodes = []
    panels = []

    x_start = -L / 2.0
    dx = L / N_x
    dtheta = math.pi / N_circ

    # --- Hull nodes: grid of (N_x+1) x (N_circ+1) ---
    hull_idx = []
    for i in range(N_x + 1):
        row = []
        x = x_start + i * dx
        for j in range(N_circ + 1):
            theta = math.pi + j * dtheta
            y = y_c + R * math.cos(theta)
            z = R * math.sin(theta)
            # Enforce exact z=0 at waterline edges
            if j == 0 or j == N_circ:
                z = 0.0
            nodes.append((x, y, z))
            row.append(len(nodes) + node_offset)
        hull_idx.append(row)

    # --- Hull panels ---
    # Ordering gives outward-pointing normals (verified analytically)
    for i in range(N_x):
        for j in range(N_circ):
            n1 = hull_idx[i][j]
            n2 = hull_idx[i][j + 1]
            n3 = hull_idx[i + 1][j + 1]
            n4 = hull_idx[i + 1][j]
            panels.append((n1, n2, n3, n4))

    # --- End cap mesh ---
    def make_end_cap(x_end, hull_i, normal_sign):
        """
        Generate end cap at given x position.
        hull_i: hull row index for boundary nodes
        normal_sign: -1 for x=-L/2 (normal in -x), +1 for x=+L/2 (normal in +x)
        """
        cap_nodes = {}  # (k, j) -> 1-based global index

        # Center node (at waterline, center of semicircle)
        nodes.append((x_end, y_c, 0.0))
        center_idx = len(nodes) + node_offset
        for j in range(N_circ + 1):
            cap_nodes[(0, j)] = center_idx

        # Interior ring nodes (k=1 to N_r-1)
        for k in range(1, N_r):
            r_k = k * R / N_r
            for j in range(N_circ + 1):
                theta = math.pi + j * dtheta
                y = y_c + r_k * math.cos(theta)
                z = r_k * math.sin(theta)
                if j == 0 or j == N_circ:
                    z = 0.0
                nodes.append((x_end, y, z))
                cap_nodes[(k, j)] = len(nodes) + node_offset

        # Outer ring (k=N_r) uses hull boundary nodes
        for j in range(N_circ + 1):
            cap_nodes[(N_r, j)] = hull_idx[hull_i][j]

        # End cap panels
        for k in range(1, N_r + 1):
            for j in range(N_circ):
                if normal_sign < 0:
                    # Normal in -x direction
                    p1 = cap_nodes[(k - 1, j)]
                    p2 = cap_nodes[(k - 1, j + 1)]
                    p3 = cap_nodes[(k, j + 1)]
                    p4 = cap_nodes[(k, j)]
                else:
                    # Normal in +x direction
                    p1 = cap_nodes[(k, j)]
                    p2 = cap_nodes[(k, j + 1)]
                    p3 = cap_nodes[(k - 1, j + 1)]
                    p4 = cap_nodes[(k - 1, j)]
                panels.append((p1, p2, p3, p4))

    # End cap at x = -L/2
    make_end_cap(x_start, 0, -1)

    # End cap at x = +L/2
    make_end_cap(x_start + L, N_x, +1)

    # --- Waterplane lid panels at z=0 (irregular frequency removal) ---
    if add_lid:
        # The lid is a rectangle at z=0 from x=-L/2 to x=+L/2, y=y_c-R to y_c+R.
        # Reuse hull waterline edge nodes at j=0 (y=y_c-R) and j=N_circ (y=y_c+R).
        # Create N_lid-1 interior rows of nodes between these edges.

        lid_idx = []  # (N_x+1) x (N_lid+1) grid of node indices

        for i in range(N_x + 1):
            row = []
            for m in range(N_lid + 1):
                if m == 0:
                    # Port waterline edge: hull node at j=0
                    row.append(hull_idx[i][0])
                elif m == N_lid:
                    # Starboard waterline edge: hull node at j=N_circ
                    row.append(hull_idx[i][N_circ])
                else:
                    # Interior lid node
                    x = x_start + i * dx
                    frac = m / N_lid
                    y_port = y_c - R  # j=0 side
                    y_stbd = y_c + R  # j=N_circ side
                    y = y_port + frac * (y_stbd - y_port)
                    nodes.append((x, y, 0.0))
                    row.append(len(nodes) + node_offset)
            lid_idx.append(row)

        # Lid panels - normal must point UP (+z, into the fluid from body interior)
        # For +z normal: use ordering that gives CCW when viewed from above
        for i in range(N_x):
            for m in range(N_lid):
                p1 = lid_idx[i][m]
                p2 = lid_idx[i + 1][m]
                p3 = lid_idx[i + 1][m + 1]
                p4 = lid_idx[i][m + 1]
                panels.append((p1, p2, p3, p4))

    return nodes, panels


# Generate mesh for both pontoons
nodes_1, panels_1 = generate_pontoon_mesh(
    y_c1, R, L, z_c, N_circ, N_x, N_r, node_offset=0,
    add_lid=use_lid, N_lid=N_lid
)
nodes_2, panels_2 = generate_pontoon_mesh(
    y_c2, R, L, z_c, N_circ, N_x, N_r, node_offset=len(nodes_1),
    add_lid=use_lid, N_lid=N_lid
)

all_nodes = nodes_1 + nodes_2
all_panels = panels_1 + panels_2
n_nodes = len(all_nodes)
n_panels = len(all_panels)

print(f"Mesh: {n_nodes} nodes, {n_panels} panels")

# ============================================================
# WRITE NEMOH MESH FILE (.dat)
# ============================================================
mesh_file = os.path.join(case_dir, "SimplePontoon.dat")

with open(mesh_file, 'w') as f:
    # Header: magic number=2, symmetry flag=0 (no symmetry)
    f.write("    2          0\n")

    # Nodes
    for idx, (x, y, z) in enumerate(all_nodes, start=1):
        f.write(f"    {idx:6d}     {x: 16.7E}     {y: 16.7E}     {z: 16.7E}\n")

    # End of nodes marker
    f.write(f"    {0:6d}     {0.0: 16.7E}     {0.0: 16.7E}     {0.0: 16.7E}\n")

    # Panels
    for n1, n2, n3, n4 in all_panels:
        f.write(f"    {n1:6d}     {n2:6d}     {n3:6d}     {n4:6d}\n")

    # End of panels marker
    f.write(f"    {0:6d}     {0:6d}     {0:6d}     {0:6d}\n")

print(f"Mesh file written: {mesh_file}")

# ============================================================
# HYDROSTATICS COMPUTATION
# ============================================================

def submerged_area(R, zc):
    if zc >= R:
        return 0.0
    if zc <= -R:
        return math.pi * R**2
    h = R - zc
    area = R**2 * math.acos((R - h) / R) - (R - h) * math.sqrt(max(2 * R * h - h * h, 0.0))
    return area


def waterplane_chord(R, zc):
    if abs(zc) >= R:
        return 0.0
    return 2.0 * math.sqrt(max(R * R - zc * zc, 0.0))


def submerged_centroid_z(R, zc, n=2000):
    z1 = zc - R
    z2 = min(0.0, zc + R)
    if z2 <= z1:
        return 0.0
    dz = (z2 - z1) / float(n)
    A = 0.0
    Mz = 0.0
    for i in range(n):
        zmid = z1 + (i + 0.5) * dz
        arg = R * R - (zmid - zc) ** 2
        if arg > 0.0:
            width = 2.0 * math.sqrt(arg)
            dA = width * dz
            A += dA
            Mz += zmid * dA
    if A <= 0.0:
        return 0.0
    return Mz / A


# Pontoon centers (centroid x=0 since pontoons span -L/2 to +L/2)
pontoon_centers = [(0.0, y_c1, z_c), (0.0, y_c2, z_c)]

# Compute hydrostatic stiffness matrix
K = [[0.0] * 6 for _ in range(6)]
Awp_total = 0.0
V_total = 0.0
zB_num = 0.0
Iwp_x = 0.0
Iwp_y = 0.0

for (xc, yc, zc_local) in pontoon_centers:
    A_sub = submerged_area(R, zc_local)
    V_i = A_sub * L
    c = waterplane_chord(R, zc_local)
    Awp_i = c * L
    zB_i = submerged_centroid_z(R, zc_local)

    V_total += V_i
    Awp_total += Awp_i
    zB_num += V_i * zB_i

    Iwp_x_local = L * c ** 3 / 12.0
    Iwp_y_local = c * L ** 3 / 12.0
    Iwp_x += Iwp_x_local + Awp_i * yc ** 2
    Iwp_y += Iwp_y_local + Awp_i * xc ** 2

z_B = zB_num / V_total if V_total > 0 else 0.0
BM_T = Iwp_x / V_total
BM_L = Iwp_y / V_total
GM_T = z_B + BM_T - z_G
GM_L = z_B + BM_L - z_G
W = rho * g * V_total

K[2][2] = rho * g * Awp_total  # Heave
K[3][3] = W * GM_T              # Roll
K[4][4] = W * GM_L              # Pitch

# ============================================================
# MASS MATRIX COMPUTATION
# ============================================================
M = [[0.0] * 6 for _ in range(6)]
total_mass = 0.0

for (xc, yc, zc_local) in pontoon_centers:
    A_sub = submerged_area(R, zc_local)
    V = A_sub * L
    m = rho * V
    total_mass += m

    Ixx = 0.5 * m * R ** 2
    Iyy = (1.0 / 12.0) * m * (3 * R ** 2 + L ** 2)
    Izz = Iyy

    # Parallel axis shift to global origin
    Ixx += m * (xc ** 2 + yc ** 2)
    Iyy += m * xc ** 2
    Izz += m * yc ** 2

    M[0][0] += m
    M[1][1] += m
    M[2][2] += m
    M[3][3] += Ixx
    M[4][4] += Iyy
    M[5][5] += Izz

print(f"\nHydrostatics:")
print(f"  Total displaced volume = {V_total:.6E} m^3")
print(f"  Total mass             = {total_mass:.6E} kg")
print(f"  Total waterplane area  = {Awp_total:.6E} m^2")
print(f"  Center of buoyancy z_B = {z_B:.6E} m")
print(f"  GM_T (roll)            = {GM_T:.6E} m")
print(f"  GM_L (pitch)           = {GM_L:.6E} m")
print(f"  K33 (heave)            = {K[2][2]:.6E} N/m")
print(f"  K44 (roll)             = {K[3][3]:.6E} N*m/rad")
print(f"  K55 (pitch)            = {K[4][4]:.6E} N*m/rad")

# ============================================================
# WRITE MECHANICS FILES
# ============================================================

def write_matrix(filepath, matrix):
    with open(filepath, 'w') as f:
        for row in matrix:
            f.write("  ".join(f"{v: .6E}" for v in row) + "\n")


write_matrix(os.path.join(mech_dir, "Inertia.dat"), M)
write_matrix(os.path.join(mech_dir, "Kh.dat"), K)
write_matrix(os.path.join(mech_dir, "Km.dat"), [[0.0] * 6 for _ in range(6)])
write_matrix(os.path.join(mech_dir, "Badd.dat"), [[0.0] * 6 for _ in range(6)])

print(f"\nMechanics files written to: {mech_dir}")

# ============================================================
# WRITE Nemoh.cal
# ============================================================
nemoh_file = os.path.join(case_dir, "Nemoh.cal")

with open(nemoh_file, 'w') as f:
    f.write("--- Environment ---\n")
    f.write(f"{rho:.1f}\t\t\t\t\t! RHO \t\t! KG/M**3 \t! Fluid specific volume \n")
    f.write(f"{g:.2f}\t\t\t\t\t! G\t\t\t! M/S**2\t! Gravity\n")
    f.write(f"{depth:.1f}\t\t\t\t\t! DEPTH\t\t\t! M\t\t! Water depth\n")
    f.write("0.\t0.\t\t\t\t\t! XEFF YEFF\t\t! M\t\t! Wave measurement point\n")
    f.write("--- Description of floating bodies ---\n")
    f.write("1\t\t\t\t\t\t! Number of bodies\n")
    f.write("--- Body 1 ---\n")
    f.write("SimplePontoon.dat\t\t! Name of mesh file\n")
    f.write(f"{n_nodes}\t{n_panels}\t\t\t\t! Number of points and number of panels\n")
    f.write("6\t\t\t\t\t\t! Number of degrees of freedom\n")
    f.write("1 1. 0. 0. 0. 0. 0.\t\t! Surge\n")
    f.write("1 0. 1. 0. 0. 0. 0.\t\t! Sway\n")
    f.write("1 0. 0. 1. 0. 0. 0.\t\t! Heave\n")
    f.write(f"2 1. 0. 0. 0. 0. {z_G:.6f}\t\t! Roll about CdG\n")
    f.write(f"2 0. 1. 0. 0. 0. {z_G:.6f}\t\t! Pitch about CdG\n")
    f.write(f"2 0. 0. 1. 0. 0. {z_G:.6f}\t\t! Yaw about CdG\n")
    f.write("6\t\t\t\t\t\t! Number of resulting generalised forces\n")
    f.write("1 1. 0. 0. 0. 0. 0.\t\t! Force in x direction\n")
    f.write("1 0. 1. 0. 0. 0. 0.\t\t! Force in y direction\n")
    f.write("1 0. 0. 1. 0. 0. 0.\t\t! Force in z direction\n")
    f.write(f"2 1. 0. 0. 0. 0. {z_G:.6f}\t\t! Moment about x at CdG\n")
    f.write(f"2 0. 1. 0. 0. 0. {z_G:.6f}\t\t! Moment about y at CdG\n")
    f.write(f"2 0. 0. 1. 0. 0. {z_G:.6f}\t\t! Moment about z at CdG\n")
    f.write("0\t\t\t\t\t\t! Number of lines of additional information\n")
    f.write("--- Load cases to be solved ---\n")
    f.write(f"1\t{n_freq}\t{omega_min:.6f}\t{omega_max:.6f}\t\t! Freq type 1=rad/s, Nfreq, Min, Max\n")
    f.write(f"{n_dir}\t{dir_min:.6f}\t{dir_max:.6f}\t\t\t! Ndir, Min, Max (degrees)\n")
    f.write("--- Post processing ---\n")
    f.write("1\t0.1\t10.\t\t\t\t! IRF (0=no), time step, duration\n")
    f.write("0\t\t\t\t\t\t! Show pressure\n")
    f.write("0\t0.\t180.\t\t\t\t! Kochin function\n")
    f.write("0\t10\t100.\t100.\t\t\t! Free surface elevation\n")
    f.write("0\t\t\t\t\t\t! RAO (0=no, 1=yes)\n")
    f.write("1\t\t\t\t\t\t! Output freq type (1=rad/s, 2=Hz, 3=s)\n")
    f.write("---QTF---\n")
    f.write("0\t\t\t\t\t\t! QTF flag (0=no, 1=yes)\n")

print(f"Nemoh.cal written: {nemoh_file}")

# ============================================================
# WRITE input_solver.txt
# ============================================================
solver_file = os.path.join(case_dir, "input_solver.txt")

with open(solver_file, 'w') as f:
    f.write("2\t\t\t\t! Gauss quadrature, N^2 nodes, specify N (1-4)\n")
    f.write("0.001\t\t\t! eps_zmin for panel z minimum\n")
    f.write("1\t\t\t\t! Linear solver: 0=Gauss, 1=LU, 2=GMRES\n")
    f.write("10 1e-5 1000\t! GMRES: restart, tolerance, max iterations\n")

print(f"input_solver.txt written: {solver_file}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("CASE SETUP COMPLETE")
print("=" * 60)
print(f"\nCase directory: {case_dir}")
print(f"Mesh file:      SimplePontoon.dat ({n_nodes} nodes, {n_panels} panels)")
print(f"\nMass matrix (diagonal):")
print(f"  m     = {M[0][0]:.6f} kg")
print(f"  Ixx   = {M[3][3]:.6f} kg*m^2")
print(f"  Iyy   = {M[4][4]:.6f} kg*m^2")
print(f"  Izz   = {M[5][5]:.6f} kg*m^2")
print(f"\nTo run Nemoh (from the SimplePontoon directory):")
print(f"  ../bin_windows/preProc.exe .")
print(f"  ../bin_windows/solver.exe .")
print(f"  ../bin_windows/postProc.exe .")
print(f"\nResults will appear in: {results_dir}")
