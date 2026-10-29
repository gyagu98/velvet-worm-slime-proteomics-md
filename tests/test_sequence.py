import os
import tempfile

from slime_proteo_md.sequence import (
    parse_fasta,
    is_nucleotide_sequence,
    translate_cds,
    parse_slime_entry,
)


def test_parse_fasta():
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write(">test1\nATGGCC\n>test2 description\nTTTAAA\n")
        path = f.name

    try:
        seqs = parse_fasta(path)
        assert len(seqs) == 2
        assert seqs["test1"] == "ATGGCC"
        assert seqs["test2 description"] == "TTTAAA"
    finally:
        os.remove(path)


def test_is_nucleotide():
    assert is_nucleotide_sequence("ATGC")
    assert is_nucleotide_sequence("atgc")
    assert is_nucleotide_sequence("ATGCN")
    assert not is_nucleotide_sequence("MKVL")


def test_translate_cds():
    assert translate_cds("ATGGCTTAA") == "MA*"
    assert translate_cds("ATGGCT") == "MA"
    assert translate_cds("ATGGC") == "M"  # trailing incomplete codon ignored


def test_parse_slime_entry_nt():
    entry = parse_slime_entry(
        "Eoperipatus_sp__PH1__PP123456_1", "ATGGCTTAA"
    )
    assert entry.is_nucleotide
    assert entry.protein_sequence == "MA"
    assert entry.family == "PH"
    assert entry.species == "Eoperipatus sp."
    assert entry.internal_stop_codons == 0


def test_parse_slime_entry_internal_stops():
    entry = parse_slime_entry(
        "Euperipatoides_rowelli__PH2__PP551270_1", "ATGTAGGCTTAA"
    )
    assert entry.is_nucleotide
    assert entry.protein_sequence == "M*A"
    assert entry.internal_stop_codons == 1


def test_parse_slime_entry_protein():
    entry = parse_slime_entry("some_protein", "MKVL")
    assert not entry.is_nucleotide
    assert entry.protein_sequence == "MKVL"
    assert entry.family == "Slime_Protein"
