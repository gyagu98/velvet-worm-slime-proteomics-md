from pathlib import Path

from slime_proteo_md.workflow import analyze_sequence_dataset


def test_supplied_fasta_analysis(tmp_path):
    root = Path(__file__).resolve().parents[1]
    summary = analyze_sequence_dataset(
        [
            str(root / "data/Peripatus_slime_7_PH_no_gaps.fasta"),
            str(root / "data/Peripatus_slime_7_PM_no_gaps.fasta"),
        ],
        str(tmp_path),
    )
    assert summary["num_entries"] == 33
    assert summary["num_nucleotide_entries"] == 33
    assert summary["num_proteins_with_internal_stops"] == 1
    assert (tmp_path / "protein_properties.csv").exists()
    assert (tmp_path / "tryptic_peptides.csv").exists()
    assert (tmp_path / "proteomics_overview.png").exists()
    assert (tmp_path / "report.html").exists()
    assert "not experimental" in (tmp_path / "analysis_report.md").read_text().lower()
    html_report = (tmp_path / "report.html").read_text().lower()
    assert "provenance" in html_report
    assert "theoretical" in html_report
    assert "experimental" in html_report
    assert "in silico" in html_report


def test_build_report_from_directories_with_smd(tmp_path):
    import json

    from slime_proteo_md.report import build_report_from_directories

    root = Path(__file__).resolve().parents[1]
    analysis_dir = tmp_path / "analysis"
    analyze_sequence_dataset(
        [
            str(root / "data/Peripatus_slime_7_PH_no_gaps.fasta"),
            str(root / "data/Peripatus_slime_7_PM_no_gaps.fasta"),
        ],
        str(analysis_dir),
    )
    smd_dir = tmp_path / "smd"
    smd_dir.mkdir()
    (smd_dir / "smd_mechanics_summary.json").write_text(json.dumps({
        "sequence": "PSEPGYPGQPSEPGYP",
        "initial_distance_nm": 1.4,
        "final_distance_nm": 2.1,
        "net_extension_nm": 0.7,
        "peak_rupture_force_pn": 123.4,
        "total_work_kj_mol": 0.8,
        "total_work_pn_nm": 1.33,
        "total_steps": 1000,
        "total_time_ps": 2.0,
        "spring_k_kj_mol_nm2": 500.0,
        "extension_at_peak_nm": 0.4,
        "pull_velocity_nm_ps": 0.02,
    }))

    report_path = tmp_path / "combined.html"
    build_report_from_directories(str(analysis_dir), str(report_path), str(smd_dir))
    report = report_path.read_text()
    assert "123.4 pN" in report
    assert "PSEPGYPGQPSEPGYP" in report
    assert "SMD Peak Force" in report
    assert "No wet-lab mass spectrometry spectra" in report
