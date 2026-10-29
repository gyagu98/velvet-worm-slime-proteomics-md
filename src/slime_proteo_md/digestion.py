"""
In silico enzymatic digestion, peptide mass fingerprinting (PMF),
proteotypic peptide identification, and theoretical MS/MS fragmentation
(b- and y-ion series) for velvet worm slime proteomics.
"""

from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass
from .proteomics import MONOISOTOPIC_AA_MASSES, calculate_molecular_weight


PROTON_MASS = 1.007276
WATER_MONO = 18.01056


@dataclass
class DigestedPeptide:
    protein_id: str
    sequence: str
    start: int
    end: int
    missed_cleavages: int
    length: int
    monoisotopic_mw: float
    mz_z1: float
    mz_z2: float
    mz_z3: float
    hydrophobicity_index: float
    is_unique: bool = False


@dataclass
class FragmentIon:
    ion_type: str  # "b" or "y"
    number: int
    charge: int
    mz: float
    loss: str = ""  # "", "-H2O", "-NH3"


def digest_sequence(
    sequence: str,
    enzyme: str = "trypsin",
    missed_cleavages: int = 1,
    min_length: int = 6,
    max_length: int = 40,
    protein_id: str = "protein",
) -> List[DigestedPeptide]:
    """
    Perform in silico enzymatic digestion with configurable missed cleavages.
    Enzyme rules:
    - Trypsin: Cleaves after K or R, not before P.
    - Lys-C: Cleaves after K, not before P.
    - Chymotrypsin: Cleaves after F, Y, W, L, not before P.
    """
    seq = sequence.upper().replace("*", "")
    n = len(seq)
    if n == 0:
        return []

    # Find cleavage sites (1-based index of last residue before cut)
    cut_sites = [0]
    for i in range(n - 1):
        curr_aa = seq[i]
        next_aa = seq[i + 1]

        cleave = False
        if enzyme.lower() == "trypsin":
            if curr_aa in ("K", "R") and next_aa != "P":
                cleave = True
        elif enzyme.lower() == "lys-c":
            if curr_aa == "K" and next_aa != "P":
                cleave = True
        elif enzyme.lower() == "chymotrypsin":
            if curr_aa in ("F", "Y", "W", "L") and next_aa != "P":
                cleave = True
        else:
            # Default simple cut
            if curr_aa in ("K", "R"):
                cleave = True

        if cleave:
            cut_sites.append(i + 1)
    cut_sites.append(n)

    peptides = []
    num_segments = len(cut_sites) - 1

    for start_idx in range(num_segments):
        for mc in range(missed_cleavages + 1):
            end_idx = start_idx + 1 + mc
            if end_idx > num_segments:
                break

            start_pos = cut_sites[start_idx] + 1
            end_pos = cut_sites[end_idx]
            pep_seq = seq[start_pos - 1 : end_pos]

            if min_length <= len(pep_seq) <= max_length:
                mono_mw = calculate_molecular_weight(pep_seq, monoisotopic=True)
                mz_1 = mono_mw + PROTON_MASS
                mz_2 = (mono_mw + 2 * PROTON_MASS) / 2.0
                mz_3 = (mono_mw + 3 * PROTON_MASS) / 3.0

                # Approximate hydrophobicity index (Krokhin SSRCalc proxy)
                h_idx = estimate_retention_index(pep_seq)

                peptides.append(
                    DigestedPeptide(
                        protein_id=protein_id,
                        sequence=pep_seq,
                        start=start_pos,
                        end=end_pos,
                        missed_cleavages=mc,
                        length=len(pep_seq),
                        monoisotopic_mw=round(mono_mw, 4),
                        mz_z1=round(mz_1, 4),
                        mz_z2=round(mz_2, 4),
                        mz_z3=round(mz_3, 4),
                        hydrophobicity_index=round(h_idx, 2),
                    )
                )

    return peptides


def estimate_retention_index(sequence: str) -> float:
    """
    Krokhin-like empirical hydrophobicity index calculation for RP-HPLC retention prediction.
    """
    aa_weights = {
        "A": 0.42, "C": 0.84, "D": -0.83, "E": -0.44, "F": 2.52,
        "G": -0.40, "H": -0.50, "I": 2.45, "K": -1.60, "L": 2.40,
        "M": 1.25, "N": -0.85, "P": -0.15, "Q": -0.55, "R": -1.30,
        "S": -0.65, "T": -0.20, "V": 1.55, "W": 3.10, "Y": 1.50,
    }
    score = sum(aa_weights.get(aa, 0.0) for aa in sequence.upper())
    # Length adjustment
    score = score * (1.0 - 0.015 * len(sequence))
    return score


def generate_ms2_fragments(
    peptide: str, max_charge: int = 2
) -> List[FragmentIon]:
    """
    Generate theoretical b- and y-ion series for MS/MS fragmentation.
    b_i = sum(residues 0..i-1) + H+
    y_i = sum(residues n-i..n-1) + H2O + H+
    """
    seq = peptide.upper()
    n = len(seq)
    fragments = []

    # Calculate cumulative residue masses
    cum_b = 0.0
    for i in range(1, n):
        res = seq[i - 1]
        cum_b += MONOISOTOPIC_AA_MASSES.get(res, 110.0)

        # b-ion (charge 1)
        b_mz_1 = cum_b + PROTON_MASS
        fragments.append(FragmentIon(ion_type="b", number=i, charge=1, mz=round(b_mz_1, 4)))

        # Neutral losses for b-ion
        if "S" in seq[:i] or "T" in seq[:i] or "E" in seq[:i] or "D" in seq[:i]:
            fragments.append(FragmentIon(ion_type="b", number=i, charge=1, mz=round(b_mz_1 - WATER_MONO, 4), loss="-H2O"))
        if "R" in seq[:i] or "K" in seq[:i] or "Q" in seq[:i] or "N" in seq[:i]:
            fragments.append(FragmentIon(ion_type="b", number=i, charge=1, mz=round(b_mz_1 - 17.02655, 4), loss="-NH3"))

        if max_charge >= 2:
            b_mz_2 = (cum_b + 2 * PROTON_MASS) / 2.0
            fragments.append(FragmentIon(ion_type="b", number=i, charge=2, mz=round(b_mz_2, 4)))

    cum_y = WATER_MONO
    for j in range(1, n):
        res = seq[n - j]
        cum_y += MONOISOTOPIC_AA_MASSES.get(res, 110.0)

        # y-ion (charge 1)
        y_mz_1 = cum_y + PROTON_MASS
        fragments.append(FragmentIon(ion_type="y", number=j, charge=1, mz=round(y_mz_1, 4)))

        # Neutral losses for y-ion
        if "S" in seq[n - j :] or "T" in seq[n - j :] or "E" in seq[n - j :] or "D" in seq[n - j :]:
            fragments.append(FragmentIon(ion_type="y", number=j, charge=1, mz=round(y_mz_1 - WATER_MONO, 4), loss="-H2O"))
        if "R" in seq[n - j :] or "K" in seq[n - j :] or "Q" in seq[n - j :] or "N" in seq[n - j :]:
            fragments.append(FragmentIon(ion_type="y", number=j, charge=1, mz=round(y_mz_1 - 17.02655, 4), loss="-NH3"))

        if max_charge >= 2:
            y_mz_2 = (cum_y + 2 * PROTON_MASS) / 2.0
            fragments.append(FragmentIon(ion_type="y", number=j, charge=2, mz=round(y_mz_2, 4)))

    return fragments


def find_proteotypic_peptides(
    protein_peptides: Dict[str, List[DigestedPeptide]]
) -> Dict[str, List[DigestedPeptide]]:
    """
    Identify proteotypic (unique) peptides across a database of proteins.
    A peptide is proteotypic if its sequence appears in only one protein entry.
    """
    sequence_counts: Dict[str, int] = {}
    for pid, peps in protein_peptides.items():
        unique_in_prot = set(p.sequence for p in peps)
        for seq in unique_in_prot:
            sequence_counts[seq] = sequence_counts.get(seq, 0) + 1

    result = {}
    for pid, peps in protein_peptides.items():
        annotated = []
        for p in peps:
            is_unique = (sequence_counts.get(p.sequence, 0) == 1)
            p.is_unique = is_unique
            annotated.append(p)
        result[pid] = annotated

    return result


def export_skyline_transition_list(
    peptides: List[DigestedPeptide], output_path: str
) -> None:
    """
    Export targeted proteomics transition list in Skyline-compatible CSV format.
    """
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("ProteinName,PeptideSequence,PrecursorMz,PrecursorCharge,ProductMz,ProductCharge,FragmentIon,Loss\n")
        for pep in peptides:
            frags = generate_ms2_fragments(pep.sequence, max_charge=1)
            # Pick top 5 diagnostic y-ions (y3 through y7)
            y_frags = [fr for fr in frags if fr.ion_type == "y" and fr.loss == "" and 3 <= fr.number <= 8]
            for fr in y_frags:
                f.write(
                    f"{pep.protein_id},{pep.sequence},{pep.mz_z2:.4f},2,{fr.mz:.4f},{fr.charge},{fr.ion_type}{fr.number},\n"
                )
