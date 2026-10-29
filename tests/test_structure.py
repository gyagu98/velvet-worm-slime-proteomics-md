from slime_proteo_md.structure import build_backbone_coordinates, generate_peptide_pdb


def test_build_backbone_coordinates():
    pdb = build_backbone_coordinates("PGYPG")
    assert pdb.startswith("ATOM")
    assert "END" in pdb
    # 4 backbone atoms per residue plus OXT on the C-terminal residue
    assert pdb.count("ATOM") == 5 * 4 + 1


def test_build_ppii_conformation():
    pdb = build_backbone_coordinates("PPGPG", conformation="polyproline_ii")
    assert "ATOM" in pdb
    assert "TER" in pdb


def test_generate_peptide_pdb():
    pdb = generate_peptide_pdb("PGYPG", conformation="extended")
    assert "ATOM" in pdb
    assert " H " in pdb or " H1 " in pdb
    assert "END" in pdb


def test_run_short_smd_demo():
    from slime_proteo_md.md_engine import run_steered_md

    pdb = generate_peptide_pdb("PGYPG")
    result = run_steered_md(
        pdb_text=pdb,
        sequence="PGYPG",
        steps=20,
        report_interval=10,
        pull_velocity_nm_ps=0.01,
        platform_name="Reference",
    )
    assert len(result.frames) == 2
    assert result.total_steps == 20
    assert result.trajectory_pdb
