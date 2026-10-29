"""
PDB structure generation, backbone geometry, and all-atom coordinate construction
for velvet worm structural repeat peptides.
"""

from typing import Optional, Tuple
import io
import math
import numpy as np
from pdbfixer import PDBFixer
from openmm import app


THREE_LETTER_MAP = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def build_backbone_coordinates(
    sequence: str, conformation: str = "extended"
) -> str:
    """
    Construct idealized backbone coordinates (N, CA, C, O, OXT) for any sequence.
    Conformations:
    - 'extended': Beta-strand-like backbone (trans, extended spacing ~3.6 A / residue)
    - 'polyproline_ii': Helical polyproline II conformation (~3.1 A / residue)
    """
    seq = sequence.upper().replace("*", "")
    lines = []
    atom_id = 1

    spacing = 3.6 if conformation == "extended" else 3.12
    # Trans peptide bond parameters
    r_n_ca = 1.458
    r_ca_c = 1.525
    r_c_o = 1.229

    for i, aa in enumerate(seq):
        res_num = i + 1
        res_name = THREE_LETTER_MAP.get(aa, "ALA")
        x_base = i * spacing

        # Slight alternate y oscillation for zig-zag backbone
        y_osc = 0.5 if (i % 2 == 1) else -0.5
        if conformation == "polyproline_ii":
            # 3 residues per turn in PPII
            angle = i * (2.0 * math.pi / 3.0)
            y_base = 0.8 * math.sin(angle)
            z_base = 0.8 * math.cos(angle)
        else:
            y_base = y_osc
            z_base = 0.0

        n_x, n_y, n_z = x_base, y_base, z_base
        ca_x, ca_y, ca_z = x_base + r_n_ca, y_base, z_base
        c_x, c_y, c_z = x_base + r_n_ca + (r_ca_c * 0.5), y_base + 1.2, z_base
        o_x, o_y, o_z = c_x - 0.4, y_base + 1.2 + r_c_o, z_base

        atoms = [
            ("N", "N", n_x, n_y, n_z),
            ("CA", "C", ca_x, ca_y, ca_z),
            ("C", "C", c_x, c_y, c_z),
            ("O", "O", o_x, o_y, o_z),
        ]

        if i == len(seq) - 1:
            # C-terminal terminal oxygen
            oxt_x = c_x + 1.2
            oxt_y = y_base + 0.8
            oxt_z = z_base
            atoms.append(("OXT", "O", oxt_x, oxt_y, oxt_z))

        for aname, elem, x, y, z in atoms:
            line = f"ATOM  {atom_id:5d} {aname:<4s} {res_name:3s} A{res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {elem:>2s}"
            lines.append(line)
            atom_id += 1

    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def generate_peptide_pdb(
    sequence: str,
    conformation: str = "extended",
    ph: float = 7.4,
    output_filepath: Optional[str] = None,
) -> str:
    """
    Generate an all-atom, hydrogenated PDB structure for any peptide sequence
    using PDBFixer and OpenMM.
    """
    backbone_pdb = build_backbone_coordinates(sequence, conformation)
    fixer = PDBFixer(pdbfile=io.StringIO(backbone_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    out_io = io.StringIO()
    app.PDBFile.writeFile(fixer.topology, fixer.positions, out_io)
    pdb_text = out_io.getvalue()

    if output_filepath:
        with open(output_filepath, "w", encoding="utf-8") as f:
            f.write(pdb_text)

    return pdb_text
