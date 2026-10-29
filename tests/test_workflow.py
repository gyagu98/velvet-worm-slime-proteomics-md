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
    assert "not experimental" in (tmp_path / "analysis_report.md").read_text().lower()
