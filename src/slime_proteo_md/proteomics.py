"""
Proteomics physicochemical characterization, charge-hydropathy profiling,
intrinsic disorder prediction (FoldIndex/Uversky), and repeat architecture
for velvet worm structural slime proteins.
"""

from typing import Dict, List, Tuple, Any
from collections import Counter
import math


# Monoisotopic and average residue masses (Da)
MONOISOTOPIC_AA_MASSES: Dict[str, float] = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276,  "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
}

AVERAGE_AA_MASSES: Dict[str, float] = {
    "A": 71.0788,  "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167,  "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}

# Kyte-Doolittle hydropathy scale
KYTE_DOOLITTLE: Dict[str, float] = {
    "A": 1.8,  "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8,  "K": -3.9, "M": 1.9,  "F": 2.8,  "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Bjellqvist pKa values
PKA_BJELLQVIST: Dict[str, float] = {
    "N_term": 7.5,
    "C_term": 3.55,
    "C": 9.0,
    "D": 4.05,
    "E": 4.45,
    "H": 5.98,
    "K": 10.0,
    "R": 12.0,
    "Y": 10.0,
}


def calculate_molecular_weight(sequence: str, monoisotopic: bool = False) -> float:
    """Calculate molecular weight in Daltons including terminal H2O (18.015 Da)."""
    table = MONOISOTOPIC_AA_MASSES if monoisotopic else AVERAGE_AA_MASSES
    water_mass = 18.01056 if monoisotopic else 18.0153
    total = water_mass
    for aa in sequence.upper():
        if aa in table:
            total += table[aa]
    return total


def calculate_gravy(sequence: str) -> float:
    """Calculate Grand Average of Hydropathy (GRAVY). Negative values indicate hydrophilic."""
    clean = [aa for aa in sequence.upper() if aa in KYTE_DOOLITTLE]
    if not clean:
        return 0.0
    return sum(KYTE_DOOLITTLE[aa] for aa in clean) / len(clean)


def calculate_net_charge(sequence: str, ph: float = 7.4) -> float:
    """Calculate net electrostatic charge of sequence at a given pH."""
    seq_upper = sequence.upper()
    if not seq_upper:
        return 0.0
    counts = Counter(seq_upper)

    # Positive contributions: N-term, K, R, H
    n_term = 1.0 / (1.0 + 10.0 ** (ph - PKA_BJELLQVIST["N_term"]))
    k_charge = counts.get("K", 0) / (1.0 + 10.0 ** (ph - PKA_BJELLQVIST["K"]))
    r_charge = counts.get("R", 0) / (1.0 + 10.0 ** (ph - PKA_BJELLQVIST["R"]))
    h_charge = counts.get("H", 0) / (1.0 + 10.0 ** (ph - PKA_BJELLQVIST["H"]))

    # Negative contributions: C-term, D, E, C, Y
    c_term = -1.0 / (1.0 + 10.0 ** (PKA_BJELLQVIST["C_term"] - ph))
    d_charge = -counts.get("D", 0) / (1.0 + 10.0 ** (PKA_BJELLQVIST["D"] - ph))
    e_charge = -counts.get("E", 0) / (1.0 + 10.0 ** (PKA_BJELLQVIST["E"] - ph))
    c_charge = -counts.get("C", 0) / (1.0 + 10.0 ** (PKA_BJELLQVIST["C"] - ph))
    y_charge = -counts.get("Y", 0) / (1.0 + 10.0 ** (PKA_BJELLQVIST["Y"] - ph))

    return n_term + k_charge + r_charge + h_charge + c_term + d_charge + e_charge + c_charge + y_charge


def calculate_isoelectric_point(sequence: str, precision: float = 0.01) -> float:
    """Calculate the isoelectric point (pI) using bisection on net charge."""
    low_ph = 0.0
    high_ph = 14.0

    while (high_ph - low_ph) > precision:
        mid_ph = (low_ph + high_ph) / 2.0
        charge = calculate_net_charge(sequence, mid_ph)
        if charge > 0:
            low_ph = mid_ph
        else:
            high_ph = mid_ph

    return round((low_ph + high_ph) / 2.0, 2)


def calculate_aliphatic_index(sequence: str) -> float:
    """Calculate aliphatic index (relative volume occupied by aliphatic side chains)."""
    counts = Counter(sequence.upper())
    total = len(sequence)
    if total == 0:
        return 0.0
    ala = counts.get("A", 0) / total * 100.0
    val = counts.get("V", 0) / total * 100.0
    ile = counts.get("I", 0) / total * 100.0
    leu = counts.get("L", 0) / total * 100.0
    return ala + (2.9 * val) + (3.9 * (ile + leu))


def calculate_polyampholyte_metrics(sequence: str) -> Dict[str, float]:
    """
    Calculate polyampholyte metrics (Pappu classification):
    - Fraction of Positive Residues (FPR: K, R)
    - Fraction of Negative Residues (FNR: D, E)
    - Fraction of Charged Residues (FCR = FPR + FNR)
    - Net Charge Per Residue (NCPR = |FPR - FNR|)
    """
    seq_upper = sequence.upper()
    total = len(seq_upper)
    if total == 0:
        return {"FPR": 0.0, "FNR": 0.0, "FCR": 0.0, "NCPR": 0.0}

    pos = seq_upper.count("K") + seq_upper.count("R")
    neg = seq_upper.count("D") + seq_upper.count("E")

    fpr = pos / total
    fnr = neg / total
    fcr = fpr + fnr
    ncpr = abs(fpr - fnr)

    return {
        "FPR": round(fpr, 4),
        "FNR": round(fnr, 4),
        "FCR": round(fcr, 4),
        "NCPR": round(ncpr, 4),
    }


def calculate_foldindex(sequence: str) -> float:
    """
    FoldIndex (Prilusky et al., Bioinformatics 2005):
    Predicts intrinsic disorder from normalized charge and hydropathy:
    FoldIndex = 2.785 * <H> - |<R>| - 1.151
    where <H> is normalized hydropathy [0, 1] and <R> is normalized net charge.
    Score < 0 signifies intrinsically disordered; Score >= 0 signifies folded.
    """
    seq_upper = sequence.upper()
    total = len(seq_upper)
    if total == 0:
        return 0.0

    # Normalized hydropathy: (KD + 4.5) / 9.0
    kd_sum = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq_upper)
    avg_h = (kd_sum / total + 4.5) / 9.0

    # Normalized net charge
    pos = seq_upper.count("K") + seq_upper.count("R")
    neg = seq_upper.count("D") + seq_upper.count("E")
    avg_r = (pos - neg) / total

    fold_score = (2.785 * avg_h) - abs(avg_r) - 1.151
    return round(fold_score, 4)


def sliding_window_profile(
    sequence: str, window_size: int = 21
) -> Dict[str, List[float]]:
    """Compute sliding window profiles for hydropathy, net charge, and FoldIndex."""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    if n < window_size:
        return {
            "positions": [n // 2],
            "hydropathy": [calculate_gravy(seq_upper)],
            "charge": [calculate_net_charge(seq_upper, 7.4)],
            "foldindex": [calculate_foldindex(seq_upper)],
        }

    positions = []
    hydropathy = []
    charge = []
    foldindex = []

    half = window_size // 2
    for i in range(n - window_size + 1):
        window = seq_upper[i : i + window_size]
        pos = i + half + 1
        positions.append(pos)
        hydropathy.append(calculate_gravy(window))
        charge.append(calculate_net_charge(window, 7.4))
        foldindex.append(calculate_foldindex(window))

    return {
        "positions": positions,
        "hydropathy": hydropathy,
        "charge": charge,
        "foldindex": foldindex,
    }


def analyze_protein_properties(entry_id: str, sequence: str, ph: float = 7.4) -> Dict[str, Any]:
    """Run full physicochemical and disorder characterization on a protein sequence."""
    mw_da = calculate_molecular_weight(sequence)
    mw_kda = mw_da / 1000.0
    pi = calculate_isoelectric_point(sequence)
    charge_74 = calculate_net_charge(sequence, ph)
    gravy = calculate_gravy(sequence)
    aliphatic = calculate_aliphatic_index(sequence)
    fold_idx = calculate_foldindex(sequence)
    polyamph = calculate_polyampholyte_metrics(sequence)

    # Amino acid composition
    counts = Counter(sequence.upper())
    total = len(sequence)
    composition = {aa: round(counts.get(aa, 0) / total * 100.0, 2) if total else 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY"}

    # Disordered fraction estimated by FoldIndex
    # A negative FoldIndex signifies predicted disorder
    is_disordered = fold_idx < 0.0

    return {
        "id": entry_id,
        "length": total,
        "mw_kda": round(mw_kda, 2),
        "pI": pi,
        "charge_ph": round(charge_74, 2),
        "ph": ph,
        "gravy": round(gravy, 3),
        "aliphatic_index": round(aliphatic, 2),
        "foldindex": fold_idx,
        "predicted_state": "Intrinsically Disordered" if is_disordered else "Globular/Folded",
        "FPR": polyamph["FPR"],
        "FNR": polyamph["FNR"],
        "FCR": polyamph["FCR"],
        "NCPR": polyamph["NCPR"],
        "composition": composition,
    }
