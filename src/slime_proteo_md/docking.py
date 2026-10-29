"""Small-molecule ligand docking workflow using AutoDock Vina and Meeko.

This module provides an in silico small-molecule docking pipeline designed
for screening candidate ligand poses against an idealized or experimental
protein/peptide target structure.

Important scientific limitation:
Docking against an idealized or unrelaxed peptide model provides a demonstration
of computational docking algorithms and pose generation; it does NOT constitute
experimental binding affinity, target engagement, or biological validation.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import urllib.request
import stat
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def ensure_vina_binary(custom_path: Optional[str] = None) -> Optional[Path]:
    """Find an existing AutoDock Vina binary or download the official release."""
    if custom_path and Path(custom_path).is_file() and os.access(custom_path, os.X_OK):
        return Path(custom_path)

    # Check PATH
    which_vina = shutil.which("vina")
    if which_vina:
        return Path(which_vina)

    # Check common cache locations
    cache_dir = Path("/tmp/vina_bin")
    cached_bin = cache_dir / "vina"
    if cached_bin.is_file() and os.access(cached_bin, os.X_OK):
        return cached_bin

    # Attempt to download official Linux x86_64 binary
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        url = "https://github.com/forlilab/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64"
        urllib.request.urlretrieve(url, cached_bin)
        cached_bin.chmod(stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
        if cached_bin.is_file():
            return cached_bin
    except Exception:
        pass

    return None


def convert_pdb_to_pdbqt(pdb_path: Path | str, pdbqt_path: Path | str) -> Path:
    """Convert a standard PDB receptor file into PDBQT format for AutoDock Vina."""
    p_in = Path(pdb_path)
    p_out = Path(pdbqt_path)

    if not p_in.is_file():
        raise FileNotFoundError(f"Input PDB file not found: {p_in}")

    with p_in.open("r", encoding="utf-8") as fin, p_out.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                element = atom_name[0] if atom_name else "C"
                ad_type = element
                prefix = line[:54]
                fout.write(f"{prefix}  1.00  0.00    +0.000 {ad_type:>2}\n")
            elif line.startswith("TER") or line.startswith("END"):
                fout.write(line)

    return p_out


def prepare_ligand_pdbqt_from_smiles(
    smiles: str,
    output_pdbqt_path: Path | str,
    random_seed: int = 42,
) -> Tuple[Path, Dict[str, Any]]:
    """Generate 3D conformer and PDBQT representation for a ligand SMILES string."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import meeko

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    mol = Chem.AddHs(mol)
    embed_code = AllChem.EmbedMolecule(mol, randomSeed=random_seed)
    if embed_code != 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=random_seed)
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = meeko.MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if not mol_setups:
        raise RuntimeError("Meeko molecule preparation yielded no valid conformers.")

    pdbqt_string, is_ok, error_msg = meeko.PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        raise RuntimeError(f"Meeko PDBQT generation failed: {error_msg}")

    p_out = Path(output_pdbqt_path)
    p_out.parent.mkdir(parents=True, exist_ok=True)
    with p_out.open("w", encoding="utf-8") as f:
        f.write(pdbqt_string)

    properties = {
        "smiles": smiles,
        "num_atoms": mol.GetNumAtoms(),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
    }
    return p_out, properties


def compute_receptor_grid_box(pdb_path: Path | str, padding_angstrom: float = 8.0) -> Dict[str, float]:
    """Compute bounding box center and dimensions from a receptor PDB file."""
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []

    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue

    if not xs:
        return {
            "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
            "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
        }

    center_x = float(np.mean(xs))
    center_y = float(np.mean(ys))
    center_z = float(np.mean(zs))

    size_x = max(15.0, float(np.ptp(xs)) + padding_angstrom)
    size_y = max(15.0, float(np.ptp(ys)) + padding_angstrom)
    size_z = max(15.0, float(np.ptp(zs)) + padding_angstrom)

    return {
        "center_x": round(center_x, 3),
        "center_y": round(center_y, 3),
        "center_z": round(center_z, 3),
        "size_x": round(size_x, 3),
        "size_y": round(size_y, 3),
        "size_z": round(size_z, 3),
    }


def parse_vina_stdout(stdout_text: str) -> List[Dict[str, Any]]:
    """Parse pose modes and affinities from AutoDock Vina stdout."""
    poses: List[Dict[str, Any]] = []
    pattern = re.compile(r"^\s*(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", re.MULTILINE)
    for match in pattern.finditer(stdout_text):
        mode, affinity, rmsd_lb, rmsd_ub = match.groups()
        poses.append({
            "mode": int(mode),
            "affinity_kcal_mol": float(affinity),
            "rmsd_lb": float(rmsd_lb),
            "rmsd_ub": float(rmsd_ub),
        })
    return poses


def run_docking(
    *,
    receptor_pdb: Path | str,
    ligand_smiles: str,
    output_dir: Path | str,
    vina_binary: Optional[str] = None,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    seed: int = 42,
    grid_box: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """Execute end-to-end AutoDock Vina small-molecule docking."""
    out_p = Path(output_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    vina_bin = ensure_vina_binary(vina_binary)
    if not vina_bin:
        raise RuntimeError("AutoDock Vina binary could not be found or executed.")

    # 1. Prepare Receptor
    rec_pdb = Path(receptor_pdb)
    rec_pdbqt = out_p / "receptor.pdbqt"
    convert_pdb_to_pdbqt(rec_pdb, rec_pdbqt)

    # 2. Prepare Ligand
    lig_pdbqt = out_p / "ligand.pdbqt"
    lig_pdbqt, lig_meta = prepare_ligand_pdbqt_from_smiles(ligand_smiles, lig_pdbqt, random_seed=seed)

    # 3. Grid Box
    box = grid_box or compute_receptor_grid_box(rec_pdb)

    # 4. Run Vina
    docked_poses_pdbqt = out_p / "docked_poses.pdbqt"
    vina_cmd = [
        str(vina_bin),
        "--receptor", str(rec_pdbqt),
        "--ligand", str(lig_pdbqt),
        "--out", str(docked_poses_pdbqt),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed),
    ]

    p = subprocess.run(vina_cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Vina docking failed (code {p.returncode}):\n{p.stderr}")

    poses = parse_vina_stdout(p.stdout)
    best_affinity = poses[0]["affinity_kcal_mol"] if poses else None

    result = {
        "status": "completed",
        "receptor_file": str(rec_pdb),
        "ligand_smiles": ligand_smiles,
        "grid_box": box,
        "vina_parameters": {
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "energy_range": energy_range,
            "seed": seed,
        },
        "num_poses": len(poses),
        "best_affinity_kcal_mol": best_affinity,
        "poses": poses,
        "files": {
            "receptor_pdbqt": str(rec_pdbqt),
            "ligand_pdbqt": str(lig_pdbqt),
            "docked_poses_pdbqt": str(docked_poses_pdbqt),
        },
    }

    with (out_p / "docking_summary.json").open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    return result
