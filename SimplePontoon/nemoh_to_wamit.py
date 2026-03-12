#!/usr/bin/env python3
"""
Convert Nemoh BEM results to WAMIT format (.1, .3, .hst) for OpenFAST HydroDyn.

Reads:
  results/CM.dat           - Added mass A(omega) [6x6 per frequency]
  results/CA.dat           - Radiation damping B(omega) [6x6 per frequency]
  results/IRF.tec          - Impulse response functions (contains A_inf)
  results/Fe.dat           - Excitation forces (magnitude + phase)
  Mechanics/Kh.dat         - Hydrostatic stiffness [6x6]

Writes:
  Wamit_format/Buoy.1     - Added mass & radiation damping
  Wamit_format/Buoy.3     - Excitation forces
  Wamit_format/Buoy.hst   - Hydrostatic restoring matrix

Non-dimensionalization (WAMITULEN = 1):
  A_nd  = A / rho
  B_nd  = B / (rho * omega)
  F_nd  = F / (rho * g)
  C_nd  = C / (rho * g)

WAMIT .1/.3 files use wave PERIOD T = 2*pi/omega as the frequency parameter.
"""

import numpy as np
import os
import math

# ============================================================
# CONFIGURATION
# ============================================================
case_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(case_dir, 'results')
mechanics_dir = os.path.join(case_dir, 'Mechanics')
output_dir = os.path.join(case_dir, 'Wamit_format')

# Physical parameters (must match Nemoh.cal / generate_case.py)
rho = 1000.0       # fluid density [kg/m^3]
g   = 9.81         # gravity [m/s^2]

WAMITULEN = 1.0    # characteristic length for non-dimensionalization
n_dof = 6
heading = 0.0      # wave heading [degrees]
output_name = 'Buoy'

# ============================================================
# PARSERS
# ============================================================

def parse_matrix_dat(filepath):
    """Parse CM.dat or CA.dat.

    Format:
        Nb de frequency :   N
        omega_1
          A(1,1) A(1,2) ... A(1,6)
          A(2,1) ...
          ...
          A(6,1) ...
        omega_2
          ...

    Returns (freqs, matrices) where matrices[k] is 6x6 numpy array.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    n_freq = int(lines[0].split(':')[1].strip())
    freqs = np.zeros(n_freq)
    matrices = []

    idx = 1
    for k in range(n_freq):
        freqs[k] = float(lines[idx].strip())
        idx += 1
        mat = np.zeros((6, 6))
        for i in range(6):
            vals = lines[idx].split()
            for j in range(6):
                mat[i, j] = float(vals[j])
            idx += 1
        matrices.append(mat)

    return freqs, matrices


def parse_irf_ainf(filepath):
    """Extract infinite-frequency added mass A_inf [6x6] from IRF.tec.

    Each Zone (DoF i) has data lines with columns:
        time, A_inf(i,1), IRF(i,1), A_inf(i,2), IRF(i,2), ..., A_inf(i,6), IRF(i,6)
    The A_inf values are constant across all time steps.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    A_inf = np.zeros((6, 6))

    zone_starts = []
    for idx, line in enumerate(lines):
        if line.strip().startswith('Zone'):
            zone_starts.append(idx)

    for dof_i, start in enumerate(zone_starts):
        data = lines[start + 1].split()
        # columns: time, [A_inf(i,j), IRF(i,j)] for j=1..6
        for j in range(6):
            A_inf[dof_i, j] = float(data[1 + 2 * j])

    return A_inf


def parse_fe_dat(filepath):
    """Parse Fe.dat (excitation forces).

    Format:
        VARIABLES=...
        Zone t="Corps  1"
        Zone t="Excitation force - beta = ... deg",I=N,F=POINT
        omega  |F1| |F2| |F3| |F4| |F5| |F6|  ang1 ang2 ang3 ang4 ang5 ang6

    Returns (freqs, magnitudes[n_freq,6], phases_deg[n_freq,6]).
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    data_lines = []
    for line in lines:
        s = line.strip()
        if s and not s.startswith('VARIABLES') and not s.startswith('Zone'):
            data_lines.append(s)

    n_freq = len(data_lines)
    freqs = np.zeros(n_freq)
    magnitudes = np.zeros((n_freq, 6))
    phases = np.zeros((n_freq, 6))

    for k, line in enumerate(data_lines):
        vals = line.split()
        freqs[k] = float(vals[0])
        for i in range(6):
            magnitudes[k, i] = float(vals[1 + i])
            phases[k, i] = float(vals[7 + i])

    return freqs, magnitudes, phases


def parse_kh_dat(filepath):
    """Parse 6x6 hydrostatic stiffness matrix from Kh.dat."""
    mat = np.zeros((6, 6))
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i >= 6:
                break
            vals = line.split()
            for j in range(6):
                mat[i, j] = float(vals[j])
    return mat


# ============================================================
# READ ALL NEMOH DATA
# ============================================================
print("Reading Nemoh results...")

freqs_cm, A_matrices = parse_matrix_dat(os.path.join(results_dir, 'CM.dat'))
freqs_ca, B_matrices = parse_matrix_dat(os.path.join(results_dir, 'CA.dat'))
A_inf = parse_irf_ainf(os.path.join(results_dir, 'IRF.tec'))
Kh = parse_kh_dat(os.path.join(mechanics_dir, 'Kh.dat'))
freqs_fe, Fe_mag, Fe_phase_deg = parse_fe_dat(os.path.join(results_dir, 'Fe.dat'))

n_freq = len(freqs_cm)
print(f"  {n_freq} frequencies: {freqs_cm[0]:.4f} to {freqs_cm[-1]:.4f} rad/s")
print(f"  Period range: {2*math.pi/freqs_cm[-1]:.4f} to {2*math.pi/freqs_cm[0]:.4f} s")
print(f"  A_inf diagonal: [{', '.join(f'{v:.6e}' for v in np.diag(A_inf))}]")
print(f"  Kh   diagonal: [{', '.join(f'{v:.6e}' for v in np.diag(Kh))}]")

# ============================================================
# WRITE WAMIT FORMAT FILES
# ============================================================
os.makedirs(output_dir, exist_ok=True)

# --- .1 file: Added mass and radiation damping ---
print(f"\nWriting {output_name}.1 ...")
with open(os.path.join(output_dir, f'{output_name}.1'), 'w') as f:

    # Block 1: freq = -1 (infinite-frequency added mass), 4 columns
    for i in range(n_dof):
        for j in range(n_dof):
            A_nd = A_inf[i, j] / rho
            f.write(f" {-1.000000e+00:14.6E}  {i+1:4d}  {j+1:4d}  {A_nd:14.6E}\n")

    # Block 2: freq = 0 (zero-frequency added mass), 4 columns
    # Approximate using lowest available frequency
    A_0 = A_matrices[0]
    for i in range(n_dof):
        for j in range(n_dof):
            A_nd = A_0[i, j] / rho
            f.write(f"  {0.000000e+00:14.6E}  {i+1:4d}  {j+1:4d}  {A_nd:14.6E}\n")

    # Positive frequency blocks: 5 columns (T, I, J, A_nd, B_nd)
    # Listed in descending period (ascending omega) order
    for k in range(n_freq):
        omega = freqs_cm[k]
        T = 2.0 * math.pi / omega
        A = A_matrices[k]
        B = B_matrices[k]
        for i in range(n_dof):
            for j in range(n_dof):
                A_nd = A[i, j] / rho
                B_nd = B[i, j] / (rho * omega)
                f.write(f"  {T:14.6E}  {i+1:4d}  {j+1:4d}  {A_nd:14.6E} {B_nd:14.6E}\n")

print(f"  -> {36 + 36 + 36*n_freq} lines ({n_freq} freq + A_inf + A_0)")

# --- .3 file: Excitation forces ---
print(f"Writing {output_name}.3 ...")
with open(os.path.join(output_dir, f'{output_name}.3'), 'w') as f:
    for k in range(n_freq):
        omega = freqs_fe[k]
        T = 2.0 * math.pi / omega
        for i in range(n_dof):
            mag_nd = Fe_mag[k, i] / (rho * g)
            phase_deg = Fe_phase_deg[k, i]
            phase_rad = math.radians(phase_deg)
            re_nd = mag_nd * math.cos(phase_rad)
            im_nd = mag_nd * math.sin(phase_rad)
            f.write(f"  {T:14.6E}  {heading:14.6E}  {i+1:4d}  "
                    f"{mag_nd:14.6E}  {phase_deg:14.6E} {re_nd:14.6E}  {im_nd:14.6E}\n")

print(f"  -> {6*n_freq} lines ({n_freq} freq x 6 DOF)")

# --- .hst file: Hydrostatic stiffness ---
print(f"Writing {output_name}.hst ...")
with open(os.path.join(output_dir, f'{output_name}.hst'), 'w') as f:
    for i in range(n_dof):
        for j in range(n_dof):
            C_nd = Kh[i, j] / (rho * g)
            f.write(f"  {i+1:4d}  {j+1:4d}    {C_nd:14.6E}\n")

print(f"  -> 36 lines")

# ============================================================
# SUMMARY
# ============================================================
print(f"\nDone! Files written to: {output_dir}")
print(f"  {output_name}.1   - Added mass & radiation damping")
print(f"  {output_name}.3   - Excitation forces")
print(f"  {output_name}.hst - Hydrostatic restoring matrix")
print(f"\nTo use with OpenFAST HydroDyn:")
print(f"  1. Copy Buoy.1, Buoy.3, Buoy.hst to your OpenFAST model directory")
print(f"  2. Set in HydroDyn.dat: PotFile=\"Buoy\", WAMITULEN={WAMITULEN}")
