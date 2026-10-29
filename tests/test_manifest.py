import json
from pathlib import Path
import pytest

from slime_proteo_md.manifest import compute_sha256, get_tool_versions, build_run_manifest


def test_compute_sha256(tmp_path: Path):
    test_file = tmp_path / "sample.txt"
    test_file.write_text("hello velvet worm slime", encoding="utf-8")
    sha = compute_sha256(test_file)
    assert len(sha) == 64
    assert sha == "5e9a65ade172431220d235e82e7d5825ca78883fd133efe845e1988fdfbe0ad5"

    missing = tmp_path / "does_not_exist.txt"
    assert compute_sha256(missing) == ""


def test_get_tool_versions():
    versions = get_tool_versions()
    assert "python" in versions
    assert "platform" in versions
    assert "biopython" in versions
    assert "pandas" in versions
    assert "openmm" in versions
    assert "click" in versions


def test_build_run_manifest(tmp_path: Path):
    f1 = tmp_path / "input1.fasta"
    f1.write_text(">seq1\nACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8")
    f2 = tmp_path / "input2.fasta"
    f2.write_text(">seq2\nAAAAAAAA\n", encoding="utf-8")

    manifest = build_run_manifest(
        command="analyze",
        parameters={"ph": 7.4, "missed_cleavages": 1},
        input_files={"fasta_1": f1, "fasta_2": f2},
        output_files={"summary": tmp_path / "missing_output.json"},
        extra_metadata={"user": "test_user"},
    )

    assert manifest["schema_version"] == "1.0"
    assert manifest["command"] == "analyze"
    assert "run_timestamp_utc" in manifest
    assert manifest["run_timestamp_utc"].endswith("Z")
    assert manifest["parameters"]["ph"] == 7.4

    # Inputs
    assert len(manifest["inputs"]) == 2
    in1 = next(item for item in manifest["inputs"] if item["label"] == "fasta_1")
    assert in1["filename"] == "input1.fasta"
    assert in1["exists"] is True
    assert in1["sha256"] is not None
    assert len(in1["sha256"]) == 64
    assert in1["size_bytes"] > 0

    # Outputs
    assert len(manifest["outputs"]) == 1
    out1 = manifest["outputs"][0]
    assert out1["exists"] is False
    assert out1["sha256"] is None

    # JSON serialization
    serialized = json.dumps(manifest)
    deserialized = json.loads(serialized)
    assert deserialized["command"] == "analyze"
