#!/usr/bin/env python3
"""
Remove irregular frequency artifacts from Nemoh BEM results using
Savitzky-Golay filter, then write cleaned WAMIT-format files.
"""

import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import os
import math

# ============================================================
# CONFIGURATION
# ============================================================
case_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(case_dir, 'results')
mechanics_dir = os.path.join(case_dir, 'Mechanics')
output_dir = os.path.join(case_dir, 'Wamit_format')

rho = 1000.0
g = 9.81
WAMITULEN = 1.0
n_dof = 6
heading = 0.0
output_name = 'Buoy'

# Savitzky-Golay parameters
SG_WINDOW = 15   # window length (odd integer)
SG_POLY   = 3    # polynomial order
CLIP_WINDOW = 7  # median filter half-window for outlier clipping
CLIP_SIGMA  = 4  # MAD multiplier for outlier threshold

# ============================================================
# PARSERS
# ============================================================

def parse_matrix_dat(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    n_freq = int(lines[0].split(':')[1].strip())
    freqs = np.zeros(n_freq)
    matrices = np.zeros((n_freq, 6, 6))
    idx = 1
    for k in range(n_freq):
        freqs[k] = float(lines[idx].strip())
        idx += 1
        for i in range(6):
            vals = lines[idx].split()
            for j in range(6):
                matrices[k, i, j] = float(vals[j])
            idx += 1
    return freqs, matrices


def parse_irf_ainf(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    A_inf = np.zeros((6, 6))
    zone_starts = [i for i, l in enumerate(lines) if l.strip().startswith('Zone')]
    for dof_i, start in enumerate(zone_starts):
        data = lines[start + 1].split()
        for j in range(6):
            A_inf[dof_i, j] = float(data[1 + 2 * j])
    return A_inf


def parse_fe_dat(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    data_lines = [l.strip() for l in lines
                  if l.strip() and not l.strip().startswith(('VARIABLES', 'Zone'))]
    n = len(data_lines)
    freqs = np.zeros(n)
    mag = np.zeros((n, 6))
    pha = np.zeros((n, 6))
    for k, line in enumerate(data_lines):
        v = line.split()
        freqs[k] = float(v[0])
        for i in range(6):
            mag[k, i] = float(v[1 + i])
            pha[k, i] = float(v[7 + i])
    return freqs, mag, pha


def parse_kh_dat(filepath):
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
# OUTLIER CLIPPING + SAVITZKY-GOLAY FILTER
# ============================================================

def clip_outliers_1d(y, half_win=CLIP_WINDOW, sigma=CLIP_SIGMA):
    """Replace outlier points with local median.
    Uses a rolling-median approach: any point deviating from its local
    median by more than sigma * local_MAD is replaced.
    """
    n = len(y)
    y_clip = y.copy()
    for k in range(n):
        lo = max(0, k - half_win)
        hi = min(n, k + half_win + 1)
        local = y[lo:hi]
        med = np.median(local)
        mad = np.median(np.abs(local - med))
        if mad < 1e-15:
            mad = np.std(local) + 1e-30
        if abs(y[k] - med) > sigma * mad:
            y_clip[k] = med
    return y_clip


def apply_savgol(freqs, A_mats, B_mats, Fe_mag, Fe_pha, window, poly):
    """Step 1: Clip extreme outliers using rolling median.
    Step 2: Apply Savitzky-Golay filter.
    Excitation forces are filtered in Re/Im space to avoid phase wrapping.
    """
    n_freq = len(freqs)
    A_f = A_mats.copy()
    B_f = B_mats.copy()
    Fm_f = Fe_mag.copy()
    Fp_f = Fe_pha.copy()

    n_clipped = 0
    for i in range(6):
        for j in range(6):
            # Clip outliers first
            a_clip = clip_outliers_1d(A_mats[:, i, j])
            b_clip = clip_outliers_1d(B_mats[:, i, j])
            n_clipped += np.sum(a_clip != A_mats[:, i, j])
            n_clipped += np.sum(b_clip != B_mats[:, i, j])
            # Then smooth
            A_f[:, i, j] = savgol_filter(a_clip, window, poly)
            B_f[:, i, j] = savgol_filter(b_clip, window, poly)

    for i in range(6):
        pha_rad = np.radians(Fe_pha[:, i])
        re = Fe_mag[:, i] * np.cos(pha_rad)
        im = Fe_mag[:, i] * np.sin(pha_rad)
        re_clip = clip_outliers_1d(re)
        im_clip = clip_outliers_1d(im)
        re_f = savgol_filter(re_clip, window, poly)
        im_f = savgol_filter(im_clip, window, poly)
        Fm_f[:, i] = np.sqrt(re_f**2 + im_f**2)
        Fp_f[:, i] = np.degrees(np.arctan2(im_f, re_f))

    print(f"  {n_clipped} outlier values clipped (A+B matrices)")
    return A_f, B_f, Fm_f, Fp_f


# ============================================================
# NEMOH .TEC WRITERS (for MATLAB time-domain script)
# ============================================================

def write_radiation_tec(filepath, freqs, A_mats, B_mats):
    """Write filtered RadiationCoefficients.tec in Nemoh format.

    Format: 6 zones (one per excitation DOF j).
    Each zone: omega, A(1,j), B(1,j), A(2,j), B(2,j), ..., A(6,j), B(6,j)
    """
    n_freq = len(freqs)
    with open(filepath, 'w') as f:
        # Header
        f.write('VARIABLES="w (rad/s)"  \n')
        for k in range(1, 7):
            f.write(f'"A   1   {k}" "B   1   {k}"\n')

        for j in range(6):  # excitation DOF (zone)
            f.write(f'Zone t="Motion of body    1 in DoF   {j+1}",I=  {n_freq:4d},F=POINT\n')
            for ki in range(n_freq):
                line = f"  {freqs[ki]:.7E}"
                for resp in range(6):  # response DOF
                    line += f"  {A_mats[ki, resp, j]:.7E}  {B_mats[ki, resp, j]:.7E}"
                f.write(line + "\n")


def write_excitation_tec(filepath, freqs, Fe_mag, Fe_pha_deg):
    """Write filtered ExcitationForce.tec in Nemoh format.

    Format: omega, |F1|, phase(F1) [rad], ..., |F6|, phase(F6) [rad]
    """
    n_freq = len(freqs)
    with open(filepath, 'w') as f:
        # Header
        f.write('VARIABLES="w (rad/s)"  \n')
        for k in range(1, 7):
            f.write(f'"abs(F   1   {k})" "angle(F   1   {k})"\n')

        f.write(f'Zone t="Excitation force - beta =   0.000 deg",I=  {n_freq:4d},F=POINT\n')
        for ki in range(n_freq):
            line = f"  {freqs[ki]:.7E}"
            for dof in range(6):
                phase_rad = math.radians(Fe_pha_deg[ki, dof])
                line += f"  {Fe_mag[ki, dof]:.7E}  {phase_rad:.7E}"
            f.write(line + "\n")


def write_cm_dat(filepath, freqs, A_mats):
    """Write filtered CM.dat (added mass) in Nemoh format."""
    n_freq = len(freqs)
    with open(filepath, 'w') as f:
        f.write(f"Nb de frequency :   {n_freq}\n")
        for k in range(n_freq):
            f.write(f" {freqs[k]:.4f}\n")
            for i in range(6):
                vals = "  ".join(f"{A_mats[k, i, j]:.6E}" for j in range(6))
                f.write(f"  {vals}\n")


def write_ca_dat(filepath, freqs, B_mats):
    """Write filtered CA.dat (damping) in Nemoh format."""
    n_freq = len(freqs)
    with open(filepath, 'w') as f:
        f.write(f"Nb de frequency :   {n_freq}\n")
        for k in range(n_freq):
            f.write(f" {freqs[k]:.4f}\n")
            for i in range(6):
                vals = "  ".join(f"{B_mats[k, i, j]:.6E}" for j in range(6))
                f.write(f"  {vals}\n")


# ============================================================
# WAMIT WRITER
# ============================================================

def write_wamit(out_dir, name, freqs, A_mats, B_mats, A_inf,
                Fe_mag, Fe_pha, Kh, rho, g, heading=0.0):
    os.makedirs(out_dir, exist_ok=True)
    n_freq = len(freqs)

    with open(os.path.join(out_dir, f'{name}.1'), 'w') as f:
        for i in range(6):
            for j in range(6):
                f.write(f" {-1.0:14.6E}  {i+1:4d}  {j+1:4d}  {A_inf[i,j]/rho:14.6E}\n")
        A_0 = A_mats[0]
        for i in range(6):
            for j in range(6):
                f.write(f"  {0.0:14.6E}  {i+1:4d}  {j+1:4d}  {A_0[i,j]/rho:14.6E}\n")
        for k in range(n_freq):
            w = freqs[k]
            T = 2.0 * math.pi / w
            for i in range(6):
                for j in range(6):
                    A_nd = A_mats[k, i, j] / rho
                    B_nd = B_mats[k, i, j] / (rho * w)
                    f.write(f"  {T:14.6E}  {i+1:4d}  {j+1:4d}  {A_nd:14.6E} {B_nd:14.6E}\n")

    with open(os.path.join(out_dir, f'{name}.3'), 'w') as f:
        for k in range(n_freq):
            w = freqs[k]
            T = 2.0 * math.pi / w
            for i in range(6):
                m_nd = Fe_mag[k, i] / (rho * g)
                p_deg = Fe_pha[k, i]
                p_rad = math.radians(p_deg)
                re_nd = m_nd * math.cos(p_rad)
                im_nd = m_nd * math.sin(p_rad)
                f.write(f"  {T:14.6E}  {heading:14.6E}  {i+1:4d}  "
                        f"{m_nd:14.6E}  {p_deg:14.6E} {re_nd:14.6E}  {im_nd:14.6E}\n")

    with open(os.path.join(out_dir, f'{name}.hst'), 'w') as f:
        for i in range(6):
            for j in range(6):
                f.write(f"  {i+1:4d}  {j+1:4d}    {Kh[i,j]/(rho*g):14.6E}\n")


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Irregular Frequency Removal (Savitzky-Golay)")
    print("=" * 60)

    # Read data
    print("\nReading Nemoh results...")
    freqs, A_raw = parse_matrix_dat(os.path.join(results_dir, 'CM.dat'))
    _, B_raw = parse_matrix_dat(os.path.join(results_dir, 'CA.dat'))
    A_inf = parse_irf_ainf(os.path.join(results_dir, 'IRF.tec'))
    Kh = parse_kh_dat(os.path.join(mechanics_dir, 'Kh.dat'))
    _, Fe_mag_raw, Fe_pha_raw = parse_fe_dat(os.path.join(results_dir, 'Fe.dat'))
    n_freq = len(freqs)
    print(f"  {n_freq} frequencies: {freqs[0]:.2f} to {freqs[-1]:.2f} rad/s")

    # Apply filter
    print(f"\nApplying Savitzky-Golay filter (window={SG_WINDOW}, poly={SG_POLY})...")
    A_sg, B_sg, Fm_sg, Fp_sg = apply_savgol(
        freqs, A_raw, B_raw, Fe_mag_raw, Fe_pha_raw, SG_WINDOW, SG_POLY)

    # Back up original Nemoh results and write filtered versions
    import shutil
    raw_backup = os.path.join(results_dir, 'raw')
    if not os.path.exists(raw_backup):
        os.makedirs(raw_backup)
        for fn in ['CM.dat', 'CA.dat', 'RadiationCoefficients.tec',
                    'ExcitationForce.tec', 'Fe.dat']:
            src = os.path.join(results_dir, fn)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(raw_backup, fn))
        print(f"\n  Original files backed up to: results/raw/")

    print("\nWriting filtered Nemoh results...")
    write_radiation_tec(os.path.join(results_dir, 'RadiationCoefficients.tec'),
                        freqs, A_sg, B_sg)
    write_excitation_tec(os.path.join(results_dir, 'ExcitationForce.tec'),
                         freqs, Fm_sg, Fp_sg)
    write_cm_dat(os.path.join(results_dir, 'CM.dat'), freqs, A_sg)
    write_ca_dat(os.path.join(results_dir, 'CA.dat'), freqs, B_sg)
    print("  -> RadiationCoefficients.tec, ExcitationForce.tec, CM.dat, CA.dat")

    # Write WAMIT files
    print("\nWriting WAMIT files...")
    write_wamit(output_dir, output_name,
                freqs, A_sg, B_sg, A_inf, Fm_sg, Fp_sg, Kh, rho, g, heading)
    print(f"  -> {output_dir}")

    # ========================================================
    # COMPARISON PLOTS
    # ========================================================
    print("\nGenerating comparison plots...")
    dof_labels = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']

    # Figure 1: Added mass diagonal (3x2)
    fig1, axes1 = plt.subplots(3, 2, figsize=(14, 10))
    fig1.suptitle('Added Mass - Raw vs Savitzky-Golay', fontsize=14, fontweight='bold')
    for d in range(6):
        ax = axes1.flat[d]
        ax.plot(freqs, A_raw[:, d, d], 'k.-', lw=0.8, ms=3, alpha=0.5, label='Raw')
        ax.plot(freqs, A_sg[:, d, d], 'b-', lw=2, label='SavGol')
        ax.set_title(f'A({d+1},{d+1}) - {dof_labels[d]}')
        ax.set_xlabel('w (rad/s)')
        ax.set_ylabel('Added mass')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(os.path.join(case_dir, 'irreg_freq_addedmass.png'), dpi=150)

    # Figure 2: Radiation damping diagonal (3x2)
    fig2, axes2 = plt.subplots(3, 2, figsize=(14, 10))
    fig2.suptitle('Radiation Damping - Raw vs Savitzky-Golay', fontsize=14, fontweight='bold')
    for d in range(6):
        ax = axes2.flat[d]
        ax.plot(freqs, B_raw[:, d, d], 'k.-', lw=0.8, ms=3, alpha=0.5, label='Raw')
        ax.plot(freqs, B_sg[:, d, d], 'b-', lw=2, label='SavGol')
        ax.set_title(f'B({d+1},{d+1}) - {dof_labels[d]}')
        ax.set_xlabel('w (rad/s)')
        ax.set_ylabel('Damping')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(case_dir, 'irreg_freq_damping.png'), dpi=150)

    # Figure 3: Excitation force magnitude (3x2)
    fig3, axes3 = plt.subplots(3, 2, figsize=(14, 10))
    fig3.suptitle('Excitation Force - Raw vs Savitzky-Golay', fontsize=14, fontweight='bold')
    for d in range(6):
        ax = axes3.flat[d]
        ax.plot(freqs, Fe_mag_raw[:, d], 'k.-', lw=0.8, ms=3, alpha=0.5, label='Raw')
        ax.plot(freqs, Fm_sg[:, d], 'b-', lw=2, label='SavGol')
        ax.set_title(f'|F{d+1}| - {dof_labels[d]}')
        ax.set_xlabel('w (rad/s)')
        ax.set_ylabel('Force')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(os.path.join(case_dir, 'irreg_freq_excitation.png'), dpi=150)

    # Figure 4: Heave detail (2x2)
    fig4, axes4 = plt.subplots(2, 2, figsize=(14, 8))
    fig4.suptitle('Heave (DOF 3) Detail', fontsize=14, fontweight='bold')

    ax = axes4[0, 0]
    ax.plot(freqs, A_raw[:, 2, 2], 'k.-', lw=0.8, ms=3, alpha=0.5, label='Raw')
    ax.plot(freqs, A_sg[:, 2, 2], 'b-', lw=2, label='SavGol')
    ax.set_title('A(3,3) - Full Range')
    ax.set_xlabel('w (rad/s)'); ax.legend(); ax.grid(True, alpha=0.3)

    # Zoomed near spike
    ax = axes4[0, 1]
    zoom = (freqs > 14) & (freqs < 24)
    ax.plot(freqs[zoom], A_raw[zoom, 2, 2], 'k.-', lw=0.8, ms=4, alpha=0.5, label='Raw')
    ax.plot(freqs[zoom], A_sg[zoom, 2, 2], 'b-', lw=2, label='SavGol')
    ax.set_title('A(3,3) - Zoomed Near Spike')
    ax.set_xlabel('w (rad/s)'); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes4[1, 0]
    ax.plot(freqs, B_raw[:, 2, 2], 'k.-', lw=0.8, ms=3, alpha=0.5, label='Raw')
    ax.plot(freqs, B_sg[:, 2, 2], 'b-', lw=2, label='SavGol')
    ax.set_title('B(3,3) - Full Range')
    ax.set_xlabel('w (rad/s)'); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes4[1, 1]
    ax.plot(freqs, Fe_mag_raw[:, 2], 'k.-', lw=0.8, ms=3, alpha=0.5, label='Raw')
    ax.plot(freqs, Fm_sg[:, 2], 'b-', lw=2, label='SavGol')
    ax.set_title('|F3| Heave Excitation')
    ax.set_xlabel('w (rad/s)'); ax.legend(); ax.grid(True, alpha=0.3)

    fig4.tight_layout()
    fig4.savefig(os.path.join(case_dir, 'irreg_freq_heave_detail.png'), dpi=150)

    plt.show()

    print("\nPlots saved to case directory:")
    print("  irreg_freq_addedmass.png")
    print("  irreg_freq_damping.png")
    print("  irreg_freq_excitation.png")
    print("  irreg_freq_heave_detail.png")
    print("\nDone!")
