"""
Trajectory and biophysical mechanics analysis:
RMSD, Radius of Gyration (Rg), Force-Extension profiles,
effective mechanical stiffness, and rupture peak detection.
"""

from typing import List, Dict, Tuple, Optional, Any
import math
import numpy as np
from .md_engine import SMDResult, SMDFrame


def calculate_radius_of_gyration(coords: np.ndarray) -> float:
    """
    Calculate the radius of gyration (Rg) for a set of atomic coordinates (N x 3).
    Rg = sqrt( (1/N) * sum_i ||r_i - r_cm||^2 )
    """
    if len(coords) == 0:
        return 0.0
    center_of_mass = np.mean(coords, axis=0)
    diff = coords - center_of_mass
    rg_sq = np.sum(diff ** 2) / len(coords)
    return float(np.sqrt(rg_sq))


def kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    """
    Compute optimal superposition RMSD between coordinate sets P and Q
    using the Kabsch algorithm. P and Q must have shape (N, 3).
    """
    if p.shape != q.shape:
        raise ValueError("Coordinate arrays must have the same shape.")

    n = p.shape[0]
    if n == 0:
        return 0.0

    # Center coordinates
    p_centered = p - np.mean(p, axis=0)
    q_centered = q - np.mean(q, axis=0)

    # Covariance matrix
    h = p_centered.T @ q_centered
    u, s, vt = np.linalg.svd(h)

    # Ensure right-handed coordinate system
    d = np.linalg.det(vt.T @ u.T)
    correction = np.eye(3)
    if d < 0.0:
        correction[2, 2] = -1.0

    rot = vt.T @ correction @ u.T
    p_rotated = (rot @ p_centered.T).T

    rmsd = np.sqrt(np.sum((p_rotated - q_centered) ** 2) / n)
    return float(rmsd)


def analyze_mechanics(smd_result: SMDResult) -> Dict[str, Any]:
    """
    Analyze the mechanical properties from an SMD pulling simulation:
    - Peak rupture force (F_max)
    - Plateau pulling force (mean force in unfolded regime)
    - Effective stiffness (k_eff in pN/nm) from initial linear response
    - Total mechanical work performed (kJ/mol and pN*nm)
    - Extension ratio (L_final / L_initial)
    """
    frames = smd_result.frames
    if not frames:
        return {}

    extensions_nm = [f.end_to_end_nm - smd_result.initial_distance_nm for f in frames]
    forces_pn = [f.force_pn for f in frames]
    times_ps = [f.time_ps for f in frames]

    peak_force = max(forces_pn)
    peak_idx = forces_pn.index(peak_force)
    extension_at_peak = extensions_nm[peak_idx]

    # Effective stiffness (slope in initial linear elastic regime)
    # Use first 25% of pulling or points before peak force
    fit_end = min(peak_idx + 1, max(3, len(forces_pn) // 4))
    if fit_end >= 3 and extensions_nm[fit_end - 1] > 0.01:
        x_sub = np.array(extensions_nm[:fit_end])
        f_sub = np.array(forces_pn[:fit_end])
        slope, _ = np.polyfit(x_sub, f_sub, 1)
        stiffness_pn_nm = max(0.0, float(slope))
    else:
        stiffness_pn_nm = 0.0

    # Total work conversion: 1 kJ/mol = 1.66054 pN*nm
    work_kj_mol = smd_result.total_work_kj_mol
    work_pn_nm = work_kj_mol * 1.660539

    ext_ratio = smd_result.final_distance_nm / max(0.001, smd_result.initial_distance_nm)

    return {
        "sequence": smd_result.sequence,
        "num_residues": smd_result.num_residues,
        "initial_distance_nm": smd_result.initial_distance_nm,
        "final_distance_nm": smd_result.final_distance_nm,
        "net_extension_nm": round(smd_result.final_distance_nm - smd_result.initial_distance_nm, 3),
        "extension_ratio": round(ext_ratio, 2),
        "peak_rupture_force_pn": round(peak_force, 2),
        "extension_at_peak_nm": round(extension_at_peak, 3),
        "effective_stiffness_pn_nm": round(stiffness_pn_nm, 2),
        "total_work_kj_mol": round(work_kj_mol, 2),
        "total_work_pn_nm": round(work_pn_nm, 2),
        "pull_velocity_nm_ps": smd_result.pull_velocity_nm_ps,
    }
