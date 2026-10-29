from slime_proteo_md.digestion import (
    digest_sequence,
    generate_ms2_fragments,
    find_proteotypic_peptides,
    DigestedPeptide,
)


def test_trypsin_cleavage_after_kr():
    # AKR should cleave after K and R, except R-P
    peps = digest_sequence("AAAKRAAAP", missed_cleavages=0, min_length=1, max_length=100)
    seqs = [p.sequence for p in peps]
    assert "AAAK" in seqs
    assert "R" in seqs
    assert "AAAP" in seqs


def test_trypsin_no_cleavage_before_proline():
    peps = digest_sequence("AAAKPAAAR", missed_cleavages=0, min_length=1, max_length=100)
    seqs = [p.sequence for p in peps]
    assert "AAAKPAAAR" in seqs


def test_missed_cleavages():
    no_mc = digest_sequence("AAAKBBB", missed_cleavages=0, min_length=1, max_length=100)
    with_mc = digest_sequence("AAAKBBB", missed_cleavages=1, min_length=1, max_length=100)
    assert len(with_mc) >= len(no_mc)
    assert any(p.missed_cleavages == 1 for p in with_mc)


def test_peptide_mz():
    peps = digest_sequence("PEPTIDER", missed_cleavages=0, min_length=1, max_length=100)
    assert len(peps) == 1
    p = peps[0]
    assert p.mz_z1 > p.mz_z2 > p.mz_z3


def test_ms2_fragments():
    frags = generate_ms2_fragments("PEPTIDE", max_charge=1)
    b_frags = [f for f in frags if f.ion_type == "b" and f.loss == ""]
    y_frags = [f for f in frags if f.ion_type == "y" and f.loss == ""]
    assert len(b_frags) == 6
    assert len(y_frags) == 6
    assert b_frags[0].number == 1
    assert y_frags[0].number == 1


def test_proteotypic_peptides():
    pep1 = DigestedPeptide("A", "UNIQUEAA", 1, 8, 0, 8, 800, 801, 401, 268, 5)
    pep2 = DigestedPeptide("B", "UNIQUEAA", 1, 8, 0, 8, 800, 801, 401, 268, 5)
    pep3 = DigestedPeptide("B", "UNIQUEBB", 9, 16, 0, 8, 800, 801, 401, 268, 5)

    result = find_proteotypic_peptides({"A": [pep1], "B": [pep2, pep3]})
    assert not result["A"][0].is_unique
    assert not result["B"][0].is_unique
    assert result["B"][1].is_unique
