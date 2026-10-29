"""
Sequence handling, FASTA parsing, translation, and validation
for velvet worm (Onychophora) slime transcripts and proteins.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


STANDARD_GENETIC_CODE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}


@dataclass
class SlimeEntry:
    header: str
    id: str
    species: str
    family: str  # "PH" (high MW proline-rich) or "PM" (mid MW) or "Other"
    raw_sequence: str
    is_nucleotide: bool
    protein_sequence: str
    cds_length: int
    protein_length: int
    stop_codons: int
    internal_stop_codons: int


def parse_fasta(filepath: str) -> Dict[str, str]:
    """Parse a FASTA file into a dictionary of header -> sequence."""
    sequences = {}
    current_header = None
    current_seq: List[str] = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = "".join(current_seq)

    return sequences


def is_nucleotide_sequence(seq: str) -> bool:
    """Check if sequence consists of nucleotide characters."""
    valid = set("ATGCN-atgcn")
    cleaned = seq.replace(" ", "").replace("\r", "").replace("\n", "")
    return set(cleaned).issubset(valid)


def translate_cds(nuc_seq: str, table: Optional[Dict[str, str]] = None) -> str:
    """Translate a nucleotide CDS string into an amino acid sequence."""
    if table is None:
        table = STANDARD_GENETIC_CODE
    nuc_clean = nuc_seq.replace("-", "").replace(" ", "").upper()
    end = len(nuc_clean) - (len(nuc_clean) % 3)
    protein = []
    for i in range(0, end, 3):
        codon = nuc_clean[i : i + 3]
        aa = table.get(codon, "X")
        protein.append(aa)
    return "".join(protein)


def parse_slime_entry(header: str, seq: str) -> SlimeEntry:
    """Classify and parse a velvet worm slime FASTA entry."""
    is_nuc = is_nucleotide_sequence(seq)
    if is_nuc:
        protein = translate_cds(seq)
        cds_len = len(seq)
    else:
        protein = seq.upper()
        cds_len = 0

    stops = protein.count("*")
    internal_stops = protein[:-1].count("*") if protein.endswith("*") else stops

    # Clean trailing stop codon for downstream physicochemical analysis
    clean_prot = protein[:-1] if protein.endswith("*") else protein

    # Determine species from the descriptive header tokens.
    header_lower = header.lower()
    if "epiperipatus_barbadensis" in header_lower:
        species = "Epiperipatus barbadensis"
    elif "euperipatoides_rowelli" in header_lower:
        species = "Euperipatoides rowelli"
    elif "eoperipatus" in header_lower:
        species = "Eoperipatus sp."
    elif "cuenca" in header_lower:
        species = "Cuenca Onychophora"
    else:
        species = "Onychophora sp."

    # Infer protein class from explicit annotations, not just compact PH/PM tokens.
    header_upper = header.upper()
    if "_PH" in header_upper or "HIGH_MW" in header_upper or "HIGH-MW" in header_upper:
        family = "PH"
    elif "_PM" in header_upper or "MID_MW" in header_upper or "MID-MW" in header_upper:
        family = "PM"
    else:
        family = "Slime_Protein"

    entry_id = header.split()[0].replace(">", "")

    return SlimeEntry(
        header=header,
        id=entry_id,
        species=species,
        family=family,
        raw_sequence=seq,
        is_nucleotide=is_nuc,
        protein_sequence=clean_prot,
        cds_length=cds_len,
        protein_length=len(clean_prot),
        stop_codons=stops,
        internal_stop_codons=internal_stops,
    )


def load_slime_dataset(filepaths: List[str]) -> List[SlimeEntry]:
    """Load multiple FASTA files and parse all slime entries."""
    entries = []
    for fp in filepaths:
        raw_dict = parse_fasta(fp)
        for h, s in raw_dict.items():
            entries.append(parse_slime_entry(h, s))
    return entries
