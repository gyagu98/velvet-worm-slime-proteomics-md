import math
from slime_proteo_md.proteomics import (
    calculate_molecular_weight,
    calculate_gravy,
    calculate_net_charge,
    calculate_isoelectric_point,
    calculate_aliphatic_index,
    calculate_polyampholyte_metrics,
    calculate_foldindex,
    sliding_window_profile,
    analyze_protein_properties,
)


def test_molecular_weight():
    # Ala-Ala: 2*71.0788 - 18.0153 (peptide bond) + 18.0153 terminal water
    assert math.isclose(calculate_molecular_weight("AA"), 160.1729, rel_tol=1e-4)
    assert calculate_molecular_weight("M") > 100


def test_gravy():
    assert calculate_gravy("IIIII") > 0
    assert calculate_gravy("DDDDD") < 0
    assert calculate_gravy("") == 0.0


def test_net_charge():
    assert calculate_net_charge("K", 7.4) > 0
    assert calculate_net_charge("D", 7.4) < 0
    assert abs(calculate_net_charge("", 7.4)) < 0.1


def test_isoelectric_point():
    assert 0 < calculate_isoelectric_point("K") < 14
    assert 0 < calculate_isoelectric_point("D") < 14


def test_aliphatic_index():
    assert calculate_aliphatic_index("VVV") > calculate_aliphatic_index("DDD")


def test_polyampholyte_metrics():
    metrics = calculate_polyampholyte_metrics("KKDDEE")
    assert metrics["FPR"] > 0
    assert metrics["FNR"] > 0
    assert metrics["FCR"] == metrics["FPR"] + metrics["FNR"]
    assert metrics["NCPR"] >= 0


def test_foldindex():
    assert calculate_foldindex("DDEEKR") < 0
    assert calculate_foldindex("ILVAFW") > -1.2


def test_sliding_window_profile():
    profile = sliding_window_profile("ACDEFGHIKLMNPQRSTVWY" * 2, window_size=10)
    assert len(profile["positions"]) > 0
    assert len(profile["positions"]) == len(profile["hydropathy"])
    assert len(profile["positions"]) == len(profile["foldindex"])


def test_analyze_protein_properties():
    props = analyze_protein_properties("test", "ACDEFGHIKLMNPQRSTVWY")
    assert props["id"] == "test"
    assert props["length"] == 20
    assert "mw_kda" in props
    assert "composition" in props
    assert len(props["composition"]) == 20
