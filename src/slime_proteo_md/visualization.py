"""
Publication-quality biophysical and proteomics visualizations:
- Proteomics physicochemical overview & family comparisons
- Intrinsic disorder & charge profile along sequences
- Steered MD Force-Extension and mechanical tensile work profiles
"""

from typing import List, Dict, Any, Optional
import matplotlib
matplotlib.use("Agg")  # Non-interactive headless backend
import matplotlib.pyplot as plt
import numpy as np

from .md_engine import SMDResult
from .sequence import SlimeEntry
from .proteomics import analyze_protein_properties


def plot_proteomics_overview(
    entries: List[SlimeEntry], output_path: str
) -> None:
    """
    Generate 3-panel proteomics overview:
    1. Proline & Key AA composition by slime protein family (PH vs PM)
    2. Molecular Weight vs. Isoelectric Point (pI)
    3. Charge vs. Hydropathy (Uversky plot)
    """
    ph_props = [
        analyze_protein_properties(e.id, e.protein_sequence)
        for e in entries if e.family == "PH" and e.internal_stop_codons == 0 and len(e.protein_sequence) > 50
    ]
    pm_props = [
        analyze_protein_properties(e.id, e.protein_sequence)
        for e in entries if e.family == "PM" and e.internal_stop_codons == 0 and len(e.protein_sequence) > 50
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=300)

    # 1. Amino acid profile comparison (P, G, K, D, E, L, S)
    key_aas = ["P", "G", "K", "D", "E", "L", "S", "Y"]
    ph_means = [np.mean([p["composition"].get(aa, 0.0) for p in ph_props]) if ph_props else 0.0 for aa in key_aas]
    pm_means = [np.mean([p["composition"].get(aa, 0.0) for p in pm_props]) if pm_props else 0.0 for aa in key_aas]

    x = np.arange(len(key_aas))
    width = 0.35

    axes[0].bar(x - width / 2, ph_means, width, label="PH (High MW)", color="#2b5c8f", edgecolor="black", alpha=0.85)
    axes[0].bar(x + width / 2, pm_means, width, label="PM (Mid MW)", color="#e27c38", edgecolor="black", alpha=0.85)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(key_aas, fontweight="bold", fontsize=11)
    axes[0].set_ylabel("Mean Amino Acid %", fontsize=11)
    axes[0].set_title("A. Slime Protein Composition Profile", fontweight="bold", fontsize=12)
    axes[0].legend(frameon=True)
    axes[0].grid(axis="y", linestyle="--", alpha=0.5)

    # 2. MW vs pI
    if ph_props:
        axes[1].scatter(
            [p["pI"] for p in ph_props],
            [p["mw_kda"] for p in ph_props],
            color="#2b5c8f", s=80, alpha=0.8, edgecolors="black", label="PH Family"
        )
    if pm_props:
        axes[1].scatter(
            [p["pI"] for p in pm_props],
            [p["mw_kda"] for p in pm_props],
            color="#e27c38", s=80, alpha=0.8, edgecolors="black", label="PM Family"
        )
    axes[1].set_xlabel("Isoelectric Point (pI)", fontsize=11)
    axes[1].set_ylabel("Molecular Weight (kDa)", fontsize=11)
    axes[1].set_title("B. Slime Proteome 2D Separation (MW vs pI)", fontweight="bold", fontsize=12)
    axes[1].legend(frameon=True)
    axes[1].grid(True, linestyle="--", alpha=0.5)

    # 3. Charge-Hydropathy (Uversky Fold Space)
    # Folded vs Disordered boundary: <R> = 2.785 * <H> - 1.151
    h_vals = np.linspace(0.2, 0.6, 50)
    boundary_r = (2.785 * h_vals) - 1.151

    axes[2].plot(h_vals, boundary_r, "k--", linewidth=1.5, label="Uversky Boundary")
    if ph_props:
        ph_h = [(p["gravy"] + 4.5) / 9.0 for p in ph_props]
        ph_r = [abs(p["FPR"] - p["FNR"]) for p in ph_props]
        axes[2].scatter(ph_h, ph_r, color="#2b5c8f", s=80, edgecolors="black", alpha=0.8, label="PH")
    if pm_props:
        pm_h = [(p["gravy"] + 4.5) / 9.0 for p in pm_props]
        pm_r = [abs(p["FPR"] - p["FNR"]) for p in pm_props]
        axes[2].scatter(pm_h, pm_r, color="#e27c38", s=80, edgecolors="black", alpha=0.8, label="PM")

    axes[2].set_xlabel("Normalized Hydropathy <H>", fontsize=11)
    axes[2].set_ylabel("Normalized Net Charge <R>", fontsize=11)
    axes[2].set_title("C. Charge-Hydropathy Phase Space", fontweight="bold", fontsize=12)
    axes[2].legend(frameon=True, fontsize=9)
    axes[2].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_disorder_profile(
    profile_data: Dict[str, List[float]], protein_id: str, output_path: str
) -> None:
    """Plot sliding window FoldIndex intrinsic disorder and net charge profile."""
    positions = profile_data["positions"]
    foldindex = profile_data["foldindex"]
    charge = profile_data["charge"]

    fig, ax1 = plt.subplots(figsize=(10, 4.5), dpi=300)

    color_fold = "#1f77b4"
    ax1.set_xlabel("Residue Position", fontsize=11)
    ax1.set_ylabel("FoldIndex Score", color=color_fold, fontsize=11)
    ax1.plot(positions, foldindex, color=color_fold, linewidth=1.5, label="FoldIndex")
    ax1.axhline(0.0, color="red", linestyle="--", linewidth=1.2, label="Disorder Boundary (<0)")
    ax1.fill_between(positions, foldindex, 0.0, where=[f < 0 for f in foldindex], color="#aec7e8", alpha=0.5, label="Disordered Region")
    ax1.tick_params(axis="y", labelcolor=color_fold)

    ax2 = ax1.twinx()
    color_charge = "#2ca02c"
    ax2.set_ylabel("Net Charge (pH 7.4)", color=color_charge, fontsize=11)
    ax2.plot(positions, charge, color=color_charge, linewidth=1.2, linestyle=":", alpha=0.8, label="Local Net Charge")
    ax2.tick_params(axis="y", labelcolor=color_charge)

    plt.title(f"Intrinsic Disorder & Charge Architecture: {protein_id}", fontweight="bold", fontsize=12)
    fig.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_smd_mechanics(smd_result: SMDResult, output_path: str) -> None:
    """
    Generate 4-panel Steered Molecular Dynamics (SMD) mechanics figure:
    1. Force vs. Extension (F vs. Delta x) with peak force callout
    2. End-to-end distance vs. Simulation time
    3. Accumulated mechanical work (W vs. Delta x)
    4. Potential energy evolution during pulling
    """
    frames = smd_result.frames
    if not frames:
        return

    d0 = smd_result.initial_distance_nm
    extensions_nm = [f.end_to_end_nm - d0 for f in frames]
    forces_pn = [f.force_pn for f in frames]
    times_ps = [f.time_ps for f in frames]
    distances_nm = [f.end_to_end_nm for f in frames]
    works_kj = [f.work_kj_mol for f in frames]
    pot_energies = [f.potential_energy_kj_mol for f in frames]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=300)

    # 1. Force vs Extension
    axes[0, 0].plot(extensions_nm, forces_pn, color="#c0392b", linewidth=2.0)
    peak_f = smd_result.peak_force_pn
    peak_idx = forces_pn.index(max(forces_pn))
    peak_ext = extensions_nm[peak_idx]
    axes[0, 0].scatter([peak_ext], [max(forces_pn)], color="black", s=60, zorder=5)
    axes[0, 0].annotate(
        f"Rupture Peak: {max(forces_pn):.1f} pN",
        xy=(peak_ext, max(forces_pn)),
        xytext=(peak_ext + 0.1, max(forces_pn) * 0.9),
        arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=6),
        fontweight="bold",
        fontsize=9,
    )
    axes[0, 0].set_xlabel("Extension Δx (nm)", fontsize=11)
    axes[0, 0].set_ylabel("Tensile Pulling Force (pN)", fontsize=11)
    axes[0, 0].set_title(f"A. Force-Extension Curve ({smd_result.sequence})", fontweight="bold", fontsize=12)
    axes[0, 0].grid(True, linestyle="--", alpha=0.5)

    # 2. End-to-End Distance vs Time
    axes[0, 1].plot(times_ps, distances_nm, color="#2980b9", linewidth=2.0, label="Peptide End-to-End")
    target_d = [d0 + (smd_result.pull_velocity_nm_ps * t) for t in times_ps]
    axes[0, 1].plot(times_ps, target_d, "k--", linewidth=1.5, label="Spring Pull Target")
    axes[0, 1].set_xlabel("Simulation Time (ps)", fontsize=11)
    axes[0, 1].set_ylabel("Distance (nm)", fontsize=11)
    axes[0, 1].set_title("B. Molecular Extension Trajectory", fontweight="bold", fontsize=12)
    axes[0, 1].legend(frameon=True)
    axes[0, 1].grid(True, linestyle="--", alpha=0.5)

    # 3. Accumulated Work
    axes[1, 0].plot(extensions_nm, works_kj, color="#8e44ad", linewidth=2.0)
    axes[1, 0].set_xlabel("Extension Δx (nm)", fontsize=11)
    axes[1, 0].set_ylabel("Mechanical Work (kJ/mol)", fontsize=11)
    axes[1, 0].set_title(f"C. Tensile Work of Drawing (Total: {smd_result.total_work_kj_mol:.1f} kJ/mol)", fontweight="bold", fontsize=12)
    axes[1, 0].grid(True, linestyle="--", alpha=0.5)

    # 4. Potential Energy
    axes[1, 1].plot(times_ps, pot_energies, color="#27ae60", linewidth=1.8)
    axes[1, 1].set_xlabel("Simulation Time (ps)", fontsize=11)
    axes[1, 1].set_ylabel("Potential Energy (kJ/mol)", fontsize=11)
    axes[1, 1].set_title("D. Potential Energy Evolution", fontweight="bold", fontsize=12)
    axes[1, 1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
