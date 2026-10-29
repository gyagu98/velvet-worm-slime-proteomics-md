"""
Command-line interface for velvet worm slime proteomics and MD analysis.
"""

import argparse
import json
import os
import sys
from pathlib import Path

from .sequence import load_slime_dataset
from .workflow import analyze_sequence_dataset, run_smd_for_motif, run_docking_for_slime_receptor
from .report import build_report_from_directories
from .proteomics import analyze_protein_properties, sliding_window_profile


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="slime-proteo-md",
        description="Proteomics analysis, molecular dynamics, and docking of velvet worm slime proteins.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Analyze subcommand
    analyze_p = subparsers.add_parser(
        "analyze",
        help="Analyze FASTA nucleotide/protein sequences, theoretical tryptic peptides, and protein properties.",
    )
    analyze_p.add_argument(
        "--fasta", action="append", required=True,
        help="FASTA file path. Repeat to add multiple input files.",
    )
    analyze_p.add_argument("--out", required=True, help="Output directory.")
    analyze_p.add_argument("--ph", type=float, default=7.4, help="pH for charge calculations (default 7.4).")
    analyze_p.add_argument("--missed-cleavages", type=int, default=1, help="Maximum tryptic missed cleavages (default 1).")
    analyze_p.add_argument("--min-peptide-length", type=int, default=6, help="Minimum peptide length in amino acids (default 6).")
    analyze_p.add_argument("--max-peptide-length", type=int, default=40, help="Maximum peptide length in amino acids (default 40).")

    # SMD subcommand
    smd_p = subparsers.add_parser(
        "smd",
        help="Generate a peptide structure and run a short steered molecular dynamics (SMD) pulling demonstration.",
    )
    smd_p.add_argument("--sequence", required=True, help="Peptide sequence in one-letter amino-acid code.")
    smd_p.add_argument("--out", required=True, help="Output directory.")
    smd_p.add_argument("--conformation", choices=["extended", "polyproline_ii"], default="extended", help="Initial idealized peptide backbone conformation.")
    smd_p.add_argument("--steps", type=int, default=2000, help="Number of integration steps (default 2000, 4 ps at 2 fs).")
    smd_p.add_argument("--report-interval", type=int, default=20, help="Save trajectory/force output every N steps (default 20).")
    smd_p.add_argument("--pull-velocity", type=float, default=0.02, help="Pulling speed in nm/ps (default 0.02).")
    smd_p.add_argument("--spring-k", type=float, default=500.0, help="Pulling spring constant in kJ/mol/nm^2 (default 500).")
    smd_p.add_argument("--temperature", type=float, default=300.0, help="Temperature in Kelvin (default 300).")
    smd_p.add_argument("--platform", choices=["CPU", "Reference"], default="CPU", help="OpenMM platform.")

    # Dock subcommand
    dock_p = subparsers.add_parser(
        "dock",
        help="Run AutoDock Vina small-molecule ligand docking against a modeled slime peptide or input PDB receptor.",
    )
    dock_p.add_argument("--ligand-smiles", default="CC(=O)Oc1ccccc1C(=O)O", help="SMILES string of ligand molecule (default: Aspirin).")
    dock_p.add_argument("--sequence", default="PSEPGYPGQPSEPGYP", help="Target peptide receptor sequence to model (default: conserved slime proline motif).")
    dock_p.add_argument("--receptor-pdb", help="Optional path to existing receptor PDB structure file.")
    dock_p.add_argument("--out", required=True, help="Output directory for docking poses and scores.")
    dock_p.add_argument("--exhaustiveness", type=int, default=4, help="Vina search exhaustiveness (default 4).")
    dock_p.add_argument("--num-modes", type=int, default=5, help="Maximum number of docking poses to return (default 5).")
    dock_p.add_argument("--seed", type=int, default=42, help="Random seed for conformer and docking reproducibility (default 42).")

    # Report subcommand
    report_p = subparsers.add_parser(
        "report",
        help="Generate a self-contained HTML report from existing proteomics, SMD, and optional docking results.",
    )
    report_p.add_argument("--analysis-dir", required=True, help="Directory containing sequence analysis CSV/JSON outputs.")
    report_p.add_argument("--smd-dir", help="Optional directory containing SMD mechanics JSON and plot outputs.")
    report_p.add_argument("--docking-dir", help="Optional directory containing docking summary JSON and plot outputs.")
    report_p.add_argument("--out", required=True, help="Output HTML file path.")

    args = parser.parse_args(argv)

    if args.command == "analyze":
        if args.missed_cleavages < 0:
            parser.error("--missed-cleavages must be >= 0")
        if args.min_peptide_length < 1 or args.max_peptide_length < args.min_peptide_length:
            parser.error("peptide length bounds must satisfy 1 <= min <= max")
        if not 0 <= args.ph <= 14:
            parser.error("--ph must be between 0 and 14")

    elif args.command == "smd":
        if args.steps <= 0 or args.report_interval <= 0 or args.steps % args.report_interval:
            parser.error("--steps must be positive and divisible by --report-interval")
        if args.pull_velocity <= 0 or args.spring_k <= 0:
            parser.error("--pull-velocity and --spring-k must be positive")

    elif args.command == "dock":
        if args.exhaustiveness < 1:
            parser.error("--exhaustiveness must be >= 1")
        if args.num_modes < 1:
            parser.error("--num-modes must be >= 1")

    if args.command == "analyze":
        print("Starting velvet worm slime proteomics analysis...")
        summary = analyze_sequence_dataset(
            fasta_paths=args.fasta,
            output_dir=args.out,
            ph=args.ph,
            min_peptide_length=args.min_peptide_length,
            max_peptide_length=args.max_peptide_length,
            missed_cleavages=args.missed_cleavages,
            cli_parameters={
                "fasta": args.fasta,
                "out": args.out,
                "ph": args.ph,
                "missed_cleavages": args.missed_cleavages,
                "min_peptide_length": args.min_peptide_length,
                "max_peptide_length": args.max_peptide_length,
            },
        )
        print("\nAnalysis completed successfully!")
        print(f"  FASTA records: {summary['num_entries']}")
        print(f"  Theoretical tryptic peptides: {summary['num_tryptic_peptides']}")
        print(f"  Proteotypic peptide sequences: {summary['num_proteotypic_peptides']}")
        print(f"  Results written to: {args.out}")

    elif args.command == "smd":
        print("Starting steered molecular dynamics simulation...")
        print(f"  Sequence: {args.sequence} ({len(args.sequence)} residues)")
        print(f"  Steps: {args.steps} ({args.steps * 0.002:.2f} ps at 2 fs)")
        mechanics = run_smd_for_motif(
            sequence=args.sequence,
            output_dir=args.out,
            conformation=args.conformation,
            steps=args.steps,
            report_interval=args.report_interval,
            pull_velocity_nm_ps=args.pull_velocity,
            spring_k_kj_mol_nm2=args.spring_k,
            temperature_k=args.temperature,
            platform_name=args.platform,
        )
        print("\nSMD simulation completed successfully!")
        print(f"  Initial distance: {mechanics['initial_distance_nm']:.3f} nm")
        print(f"  Final distance: {mechanics['final_distance_nm']:.3f} nm")
        print(f"  Peak force: {mechanics['peak_rupture_force_pn']:.2f} pN")
        print(f"  Mechanical work: {mechanics['total_work_kj_mol']:.2f} kJ/mol")
        print(f"  Results written to: {args.out}")

    elif args.command == "dock":
        print("Starting small-molecule AutoDock Vina docking...")
        print(f"  Ligand SMILES: {args.ligand_smiles}")
        print(f"  Target receptor: {args.receptor_pdb or f'Modeled peptide ({args.sequence})'}")
        dock_res = run_docking_for_slime_receptor(
            receptor_sequence=args.sequence,
            receptor_pdb_path=args.receptor_pdb,
            ligand_smiles=args.ligand_smiles,
            output_dir=args.out,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            seed=args.seed,
        )
        print("\nDocking completed successfully!")
        print(f"  Generated poses: {dock_res['num_poses']}")
        if dock_res["best_affinity_kcal_mol"] is not None:
            print(f"  Best binding affinity: {dock_res['best_affinity_kcal_mol']:.2f} kcal/mol")
        print(f"  Docked poses PDBQT: {dock_res['files']['docked_poses_pdbqt']}")
        print(f"  Results written to: {args.out}")

    elif args.command == "report":
        build_report_from_directories(
            analysis_dir=args.analysis_dir,
            smd_dir=args.smd_dir,
            docking_dir=args.docking_dir,
            output_html_path=args.out,
        )
        print(f"HTML report written to: {args.out}")


if __name__ == "__main__":
    main()
