"""
Reproducible end-to-end workflow for velvet worm slime protein proteomics
and molecular dynamics (sequence analysis, theoretical MS, and SMD pulling).
"""

from typing import List, Dict, Any, Optional
import os
import csv
import json
from pathlib import Path

from .sequence import load_slime_dataset, SlimeEntry
from .proteomics import analyze_protein_properties, sliding_window_profile
from .digestion import digest_sequence, find_proteotypic_peptides, export_skyline_transition_list
from .structure import generate_peptide_pdb
from .md_engine import run_steered_md
from .analysis import analyze_mechanics
from .visualization import plot_proteomics_overview, plot_disorder_profile, plot_smd_mechanics
from .report import generate_html_report
from .manifest import build_run_manifest


def ensure_directory(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def analyze_sequence_dataset(
    fasta_paths: List[str], output_dir: str, ph: float = 7.4,
    min_peptide_length: int = 6, max_peptide_length: int = 40,
    missed_cleavages: int = 1,
    cli_parameters: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Run sequence/physicochemical analysis and in silico proteomics workflow."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    entries = load_slime_dataset(fasta_paths)
    if not entries:
        raise ValueError("No sequences found in FASTA input files.")

    # 1. Protein-level physicochemical table
    protein_rows = []
    all_peptides = {}
    protein_sequences = {}
    for entry in entries:
        props = analyze_protein_properties(entry.id, entry.protein_sequence, ph=ph)
        protein_sequences[entry.id] = entry.protein_sequence
        row = {
            "protein_id": entry.id,
            "species": entry.species,
            "family": entry.family,
            "sequence_type": "CDS (translated)" if entry.is_nucleotide else "protein",
            "cds_length_nt": entry.cds_length,
            "protein_length_aa": entry.protein_length,
            "internal_stop_codons": entry.internal_stop_codons,
            "mw_kda": props["mw_kda"],
            "pI": props["pI"],
            "charge_ph": props["charge_ph"],
            "analysis_ph": props["ph"],
            "gravy": props["gravy"],
            "proline_percent": props["composition"].get("P", 0.0),
            "glycine_percent": props["composition"].get("G", 0.0),
            "lysine_percent": props["composition"].get("K", 0.0),
            "aspartate_percent": props["composition"].get("D", 0.0),
            "glutamate_percent": props["composition"].get("E", 0.0),
            "foldindex": props["foldindex"],
            "predicted_state": props["predicted_state"],
            "FPR": props["FPR"],
            "FNR": props["FNR"],
            "FCR": props["FCR"],
            "NCPR": props["NCPR"],
        }
        protein_rows.append(row)

        # 2. In silico tryptic digest
        peps = [] if entry.internal_stop_codons else digest_sequence(
            entry.protein_sequence,
            enzyme="trypsin",
            missed_cleavages=missed_cleavages,
            min_length=min_peptide_length,
            max_length=max_peptide_length,
            protein_id=entry.id,
        )
        all_peptides[entry.id] = peps

    # 3. Identify proteotypic peptides across the full combined protein set
    annotated_peptides = find_proteotypic_peptides(all_peptides)

    # 4. Export protein summary CSV
    protein_csv_path = output_dir / "protein_properties.csv"
    if protein_rows:
        with open(protein_csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(protein_rows[0].keys()))
            writer.writeheader()
            writer.writerows(protein_rows)

    # 5. Export peptide-level table
    peptide_rows = []
    all_unique_peptides = []
    for pid, peps in annotated_peptides.items():
        entry = next(e for e in entries if e.id == pid)
        for pep in peps:
            row = {
                "protein_id": pep.protein_id,
                "species": entry.species,
                "family": entry.family,
                "peptide_sequence": pep.sequence,
                "start": pep.start,
                "end": pep.end,
                "length_aa": pep.length,
                "missed_cleavages": pep.missed_cleavages,
                "monoisotopic_mw": pep.monoisotopic_mw,
                "mz_z1": pep.mz_z1,
                "mz_z2": pep.mz_z2,
                "mz_z3": pep.mz_z3,
                "retention_index": pep.hydrophobicity_index,
                "is_proteotypic": pep.is_unique,
                "tryptic_termini": "Fully tryptic",
                "note": "Theoretical only: confirm by LC-MS/MS before biological interpretation",
            }
            peptide_rows.append(row)
            if pep.is_unique:
                all_unique_peptides.append(pep)

    peptide_csv_path = output_dir / "tryptic_peptides.csv"
    if peptide_rows:
        with open(peptide_csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(peptide_rows[0].keys()))
            writer.writeheader()
            writer.writerows(peptide_rows)

    # 6. Export Skyline transition list for unique peptides
    if all_unique_peptides:
        export_skyline_transition_list(
            all_unique_peptides,
            str(output_dir / "skyline_transition_list.csv"),
        )

    # 7. Overview visualization
    plot_proteomics_overview(entries, str(output_dir / "proteomics_overview.png"))

    # 8. Generate per-protein disorder profiles for representative entries
    representative_profiles = []
    for entry in entries:
        if len(entry.protein_sequence) < 100:
            continue
        profile = sliding_window_profile(entry.protein_sequence, window_size=21)
        safe_id = "".join(c for c in entry.id if c.isalnum() or c in "_-.")[:60]
        plot_disorder_profile(
            profile,
            entry.id,
            str(output_dir / f"{safe_id}_disorder_profile.png"),
        )
        representative_profiles.append(entry.id)
        if len(representative_profiles) >= 6:
            break

    summary = {
        "analysis_type": "velvet-worm-slime-proteomics-md",
        "num_entries": len(entries),
        "num_nucleotide_entries": sum(1 for e in entries if e.is_nucleotide),
        "num_protein_entries": sum(1 for e in entries if not e.is_nucleotide),
        "num_PH_entries": sum(1 for e in entries if e.family == "PH"),
        "num_PM_entries": sum(1 for e in entries if e.family == "PM"),
        "num_other_slime_entries": sum(1 for e in entries if e.family == "Slime_Protein"),
        "num_proteins_with_internal_stops": sum(1 for e in entries if e.internal_stop_codons > 0),
        "num_tryptic_peptides": len(peptide_rows),
        "num_proteotypic_peptides": sum(1 for p in peptide_rows if p["is_proteotypic"]),
        "fasta_inputs": fasta_paths,
        "ph": ph,
        "enzyme": "trypsin",
        "missed_cleavages": missed_cleavages,
        "min_peptide_length": min_peptide_length,
        "max_peptide_length": max_peptide_length,
        "representative_profiles": representative_profiles,
        "caveats": [
            "FASTAs are nucleotide CDS/transcript sequences and were translated in the +1 reading frame.",
            "One sequence labeled Euperipatoides_rowelli__PH2__PP551270_1 contains 41 internal stop codons in the +1 frame; its translation is flagged but generated values should not be interpreted until validated.",
            "Theoretical tryptic peptides and m/z values are hypotheses only; no experimental mass-spectrometry data were provided.",
            "No actual protein expression, abundance, PTM, peptide detection, or biological function is inferred.",
        ],
        "generated_outputs": [
            "protein_properties.csv",
            "tryptic_peptides.csv",
            "skyline_transition_list.csv" if all_unique_peptides else None,
            "proteomics_overview.png",
            *[f"{p}_disorder_profile.png" for p in representative_profiles],
        ],
    }

    with open(output_dir / "sequence_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # 9. Detailed Markdown report
    report_path = output_dir / "analysis_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Velvet Worm Slime Proteomics & MD Analysis\n\n")
        f.write("## Sequence input inventory\n\n")
        f.write(f"- Total FASTA records: {len(entries)}\n")
        f.write(f"- Nucleotide CDS/transcripts translated in +1 frame: {summary['num_nucleotide_entries']}\n")
        f.write(f"- Proteins with internal stop codons in +1 translation: {summary['num_proteins_with_internal_stops']}\n")
        f.write(f"- PH family labels: {summary['num_PH_entries']}\n")
        f.write(f"- PM family labels: {summary['num_PM_entries']}\n\n")
        f.write("## Proteomics-oriented analyses\n\n")
        f.write("- Molecular weight, predicted pI, net charge at pH 7.4, GRAVY, amino-acid composition, and FoldIndex-based disorder classification.\n")
        f.write("- In silico tryptic digestion with configurable missed cleavages.\n")
        f.write("- Theoretical monoisotopic precursor m/z and b/y fragment ions for use as a hypothesis generator.\n")
        f.write(f"- Theoretical tryptic peptides: {len(peptide_rows)}\n")
        f.write(f"- Peptide sequences unique among this FASTA collection: {summary['num_proteotypic_peptides']}\n\n")
        f.write("## Molecular dynamics\n\n")
        f.write("The MD demonstration builds short peptide structures from selected sequence motifs, performs energy minimization, and applies steered molecular dynamics (SMD) to generate a force-extension trace. This illustrates molecular-dynamics methods, but it does **not** simulate the native full-length slime protein, an assembled slime network, secretion, drying, or material behavior.\n\n")
        f.write("## Interpretation limits\n\n")
        f.write("- The provided FASTAs contain nucleotide sequences, not mass-spectrometry spectra. Therefore the project generates predicted/theoretical MS results only.\n")
        f.write("- Translation uses the standard genetic code and the +1 frame. Inspect sequence provenance and quality before biological conclusions.\n")
        f.write("- The PH/PM family assignments come from FASTA annotations; they are not independently validated.\n")
        f.write("- FoldIndex and Uversky classification are sequence-based heuristics, not direct measurements of structure or disorder.\n")
        f.write("- SMD force curves are highly dependent on initial structure, force field, spring constant, pulling speed, temperature, and coordinate construction. Treat the included values as teaching output, not experimental material parameters.\n")
        f.write("- The `Euperipatoides_rowelli__PH2__PP551270_1` record has 41 internal stops in +1 translation. It remains in generated tables and must be investigated before interpreting its protein-level values or using it for MD/proteomics validation.\n\n")
        f.write("## Reproducibility\n\n")
        f.write("Run the sequence analysis with:\n\n")
        f.write("```bash\nslime-proteo-md analyze --fasta data/Peripatus_slime_7_PH_no_gaps.fasta --fasta data/Peripatus_slime_7_PM_no_gaps.fasta --out results/\n```\n")
        f.write("Run a small SMD example with:\n\n")
        f.write("```bash\nslime-proteo-md smd --sequence PSEPGYPGQPSEPGYP --out results/smd_example/ --steps 1000 --pull-velocity 0.02\n```\n")
        f.write("\n**Workflow status: research scaffold.** The provided sequence files were analyzed computationally; no wet-lab proteomics data were supplied.\n")

    # 10. Record exact inputs, parameters, environment versions, and UTC run time.
    manifest = build_run_manifest(
        command="analyze",
        parameters=cli_parameters or {
            "fasta": [str(path) for path in fasta_paths],
            "out": str(output_dir),
            "ph": ph,
            "missed_cleavages": missed_cleavages,
            "min_peptide_length": min_peptide_length,
            "max_peptide_length": max_peptide_length,
        },
        input_files={f"fasta_{i + 1}": path for i, path in enumerate(fasta_paths)},
    )
    manifest_path = output_dir / "run_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    summary["generated_outputs"].extend(["run_manifest.json", "report.html"])
    summary["generated_outputs"] = list(dict.fromkeys(x for x in summary["generated_outputs"] if x))
    with open(output_dir / "sequence_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # 11. Self-contained HTML report (theoretical results only unless an external SMD run is supplied)
    generate_html_report(
        summary=summary,
        protein_rows=protein_rows,
        peptide_rows=peptide_rows,
        overview_img_path=str(output_dir / "proteomics_overview.png"),
        output_html_path=str(output_dir / "report.html"),
        run_manifest=manifest,
    )
    if "report.html" not in summary["generated_outputs"]:
        summary["generated_outputs"].append("report.html")
        with open(output_dir / "sequence_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    return summary


def run_smd_for_motif(
    sequence: str,
    output_dir: str,
    conformation: str = "extended",
    steps: int = 2000,
    report_interval: int = 20,
    pull_velocity_nm_ps: float = 0.02,
    spring_k_kj_mol_nm2: float = 500.0,
    temperature_k: float = 300.0,
    platform_name: str = "CPU",
) -> Dict[str, Any]:
    """Generate peptide coordinates, run SMD, and write trajectory/mechanics outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    clean_seq = sequence.upper().replace("*", "")
    if len(clean_seq) < 5:
        raise ValueError("SMD peptide must contain at least 5 amino acids.")
    if not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in clean_seq):
        raise ValueError("SMD sequence contains unsupported or ambiguous amino acids.")

    pdb_path = output_dir / "initial_peptide.pdb"
    pdb_text = generate_peptide_pdb(
        clean_seq,
        conformation=conformation,
        output_filepath=str(pdb_path),
    )

    result = run_steered_md(
        pdb_text=pdb_text,
        sequence=clean_seq,
        steps=steps,
        report_interval=report_interval,
        pull_velocity_nm_ps=pull_velocity_nm_ps,
        spring_k_kj_mol_nm2=spring_k_kj_mol_nm2,
        temperature_k=temperature_k,
        platform_name=platform_name,
    )

    # Write trajectory
    with open(output_dir / "smd_trajectory.pdb", "w", encoding="utf-8") as f:
        f.write(result.trajectory_pdb)

    # Write frame-level data
    frame_csv = output_dir / "smd_force_extension.csv"
    with open(frame_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(result.frames[0].__dict__.keys()) if result.frames else [])
        writer.writeheader()
        for frame in result.frames:
            writer.writerow(frame.__dict__)

    # Mechanics summary plus explicit simulation settings for report provenance.
    mechanics = analyze_mechanics(result)
    mechanics.update({
        "total_steps": result.total_steps,
        "total_time_ps": result.total_time_ps,
        "spring_k_kj_mol_nm2": spring_k_kj_mol_nm2,
        "temperature_k": temperature_k,
        "conformation": conformation,
        "platform": platform_name,
        "total_work_kcal_mol": round(mechanics.get("total_work_kj_mol", 0.0) / 4.184, 4),
    })
    with open(output_dir / "smd_mechanics_summary.json", "w", encoding="utf-8") as f:
        json.dump(mechanics, f, indent=2)

    smd_manifest = build_run_manifest(
        command="smd",
        parameters={
            "sequence": clean_seq,
            "out": str(output_dir),
            "conformation": conformation,
            "steps": steps,
            "report_interval": report_interval,
            "pull_velocity_nm_ps": pull_velocity_nm_ps,
            "spring_k_kj_mol_nm2": spring_k_kj_mol_nm2,
            "temperature_k": temperature_k,
            "platform": platform_name,
        },
        input_files={"initial_peptide_pdb": pdb_path},
    )
    with (output_dir / "smd_run_manifest.json").open("w", encoding="utf-8") as f:
        json.dump(smd_manifest, f, indent=2)

    # SMD figure
    plot_smd_mechanics(result, str(output_dir / "smd_mechanics.png"))

    # If the SMD directory is nested under a completed sequence analysis directory,
    # regenerate a combined report with both data products.
    analysis_dir = output_dir.parent
    if (analysis_dir / "protein_properties.csv").is_file():
        from .report import build_report_from_directories
        build_report_from_directories(
            analysis_dir=str(analysis_dir),
            smd_dir=str(output_dir),
            output_html_path=str(analysis_dir / "report.html"),
        )

    return mechanics
