"""Tests for small-molecule AutoDock Vina and Meeko docking module."""

import os
import tempfile
from pathlib import Path

import pytest

from slime_proteo_md.docking import (
    compute_receptor_grid_box,
    convert_pdb_to_pdbqt,
    ensure_vina_binary,
    parse_vina_stdout,
    prepare_ligand_pdbqt_from_smiles,
    run_docking,
)
from slime_proteo_md.structure import generate_peptide_pdb
from slime_proteo_md.workflow import run_docking_for_slime_receptor


def test_ensure_vina_binary():
    bin_path = ensure_vina_binary()
    assert bin_path is not None
    assert bin_path.is_file()
    assert os.access(bin_path, os.X_OK)


def test_convert_pdb_to_pdbqt():
    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = Path(tmp) / "peptide.pdb"
        pdbqt_path = Path(tmp) / "peptide.pdbqt"

        generate_peptide_pdb("PGLP", output_filepath=str(pdb_path))
        out = convert_pdb_to_pdbqt(pdb_path, pdbqt_path)

        assert out.is_file()
        assert out.stat().st_size > 0
        content = out.read_text()
        assert "+0.000" in content


def test_prepare_ligand_pdbqt_from_smiles():
    with tempfile.TemporaryDirectory() as tmp:
        lig_path = Path(tmp) / "aspirin.pdbqt"
        out_path, meta = prepare_ligand_pdbqt_from_smiles("CC(=O)Oc1ccccc1C(=O)O", lig_path)

        assert out_path.is_file()
        assert meta["num_atoms"] > 0
        content = out_path.read_text()
        assert "REMARK" in content or "ROOT" in content


def test_compute_receptor_grid_box():
    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = Path(tmp) / "peptide.pdb"
        generate_peptide_pdb("PGLPG", output_filepath=str(pdb_path))

        box = compute_receptor_grid_box(pdb_path)
        assert "center_x" in box
        assert "size_x" in box
        assert box["size_x"] >= 15.0


def test_parse_vina_stdout():
    sample_stdout = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -2.150          0          0
   2       -1.920      2.100      3.200
"""
    poses = parse_vina_stdout(sample_stdout)
    assert len(poses) == 2
    assert poses[0]["mode"] == 1
    assert poses[0]["affinity_kcal_mol"] == -2.150


def test_run_docking_for_slime_receptor_smoke():
    with tempfile.TemporaryDirectory() as tmp:
        res = run_docking_for_slime_receptor(
            receptor_sequence="PSEPGYP",
            ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
            output_dir=tmp,
            exhaustiveness=1,
            num_modes=2,
            seed=42,
        )

        assert res["status"] == "completed"
        assert res["num_poses"] >= 1
        assert res["best_affinity_kcal_mol"] is not None
        assert Path(res["files"]["docked_poses_pdbqt"]).is_file()
        assert Path(res["files"]["affinity_plot_png"]).is_file()
        assert (Path(tmp) / "docking_run_manifest.json").is_file()
