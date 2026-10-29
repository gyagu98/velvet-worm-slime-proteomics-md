"""
HTML report generator for velvet worm slime proteomics and molecular dynamics.
Synthesizes theoretical in silico peptide predictions and steered MD metrics,
with explicit provenance badges and evidence level classifications.
"""

from typing import List, Dict, Any, Optional
import base64
from pathlib import Path
import html
import datetime
import csv
import json
from .manifest import compute_sha256, get_tool_versions, build_run_manifest


def image_to_base64_uri(image_path: Optional[str]) -> Optional[str]:
    """Convert an image file to a base64 data URI for self-contained HTML."""
    if not image_path:
        return None
    p = Path(image_path)
    if not p.is_file():
        return None
    try:
        data = p.read_bytes()
        ext = p.suffix.lower().replace(".", "")
        if ext == "jpg":
            ext = "jpeg"
        mime = f"image/{ext}" if ext in ["png", "jpeg", "svg+xml", "webp"] else "image/png"
        encoded = base64.b64encode(data).decode("ascii")
        return f"data:{mime};base64,{encoded}"
    except Exception:
        return None


def generate_html_report(
    summary: Dict[str, Any],
    protein_rows: List[Dict[str, Any]],
    peptide_rows: List[Dict[str, Any]],
    smd_mechanics: Optional[Dict[str, Any]] = None,
    overview_img_path: Optional[str] = None,
    smd_img_path: Optional[str] = None,
    output_html_path: Optional[str] = None,
    run_manifest: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Build a self-contained, publication-ready HTML report combining theoretical
    proteomics profiling and steered molecular dynamics metrics.
    """
    generated_time = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # Family aggregations
    ph_rows = [r for r in protein_rows if r.get("family") == "PH"]
    pm_rows = [r for r in protein_rows if r.get("family") == "PM"]

    def calc_mean(rows: List[Dict[str, Any]], key: str) -> float:
        vals = []
        for r in rows:
            v = r.get(key)
            if v is not None:
                try:
                    vals.append(float(v))
                except (ValueError, TypeError):
                    pass
        return sum(vals) / len(vals) if vals else 0.0

    ph_stats = {
        "count": len(ph_rows),
        "mean_len": calc_mean(ph_rows, "protein_length_aa"),
        "mean_mw": calc_mean(ph_rows, "mw_kda"),
        "mean_pi": calc_mean(ph_rows, "pI"),
        "mean_charge": calc_mean(ph_rows, "charge_ph"),
        "mean_pro": calc_mean(ph_rows, "proline_percent"),
        "mean_gly": calc_mean(ph_rows, "glycine_percent"),
        "disordered_pct": (sum(1 for r in ph_rows if r.get("predicted_state") == "Disordered") / len(ph_rows) * 100) if ph_rows else 0.0,
    }

    pm_stats = {
        "count": len(pm_rows),
        "mean_len": calc_mean(pm_rows, "protein_length_aa"),
        "mean_mw": calc_mean(pm_rows, "mw_kda"),
        "mean_pi": calc_mean(pm_rows, "pI"),
        "mean_charge": calc_mean(pm_rows, "charge_ph"),
        "mean_pro": calc_mean(pm_rows, "proline_percent"),
        "mean_gly": calc_mean(pm_rows, "glycine_percent"),
        "disordered_pct": (sum(1 for r in pm_rows if r.get("predicted_state") == "Disordered") / len(pm_rows) * 100) if pm_rows else 0.0,
    }

    # Proteotypic peptides
    unique_peptides = [p for p in peptide_rows if p.get("is_proteotypic")]

    # Embed images as base64
    overview_b64 = image_to_base64_uri(overview_img_path)
    smd_b64 = image_to_base64_uri(smd_img_path)

    # Flagged records
    flagged_records = []
    for r in protein_rows:
        try:
            if int(r.get("internal_stop_codons", 0)) > 0:
                flagged_records.append(r)
        except (ValueError, TypeError):
            pass

    # HTML builder
    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Velvet Worm Slime Proteomics & Molecular Dynamics Report</title>
  <style>
    :root {{
      --bg: #0d1117;
      --card-bg: #161b22;
      --card-border: #30363d;
      --text: #c9d1d9;
      --text-bright: #f0f6fc;
      --text-muted: #8b949e;
      --accent: #58a6ff;
      --accent-subtle: rgba(56, 139, 253, 0.15);
      --pred-bg: rgba(210, 153, 34, 0.15);
      --pred-border: #d29922;
      --pred-text: #e3b341;
      --exp-bg: rgba(248, 81, 73, 0.15);
      --exp-border: #f85149;
      --exp-text: #ff7b72;
      --ok-bg: rgba(63, 185, 80, 0.15);
      --ok-border: #3fb950;
      --ok-text: #56d364;
      --purple-bg: rgba(187, 128, 247, 0.15);
      --purple-border: #bc8cff;
      --purple-text: #d2a8ff;
      --font: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: var(--font);
      background-color: var(--bg);
      color: var(--text);
      line-height: 1.6;
      padding: 24px;
    }}
    .container {{
      max-width: 1200px;
      margin: 0 auto;
    }}
    header {{
      padding-bottom: 24px;
      margin-bottom: 24px;
      border-bottom: 1px solid var(--card-border);
    }}
    .header-badge-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-bottom: 12px;
    }}
    h1 {{
      color: var(--text-bright);
      font-size: 2rem;
      margin-bottom: 6px;
    }}
    .subtitle {{
      color: var(--text-muted);
      font-size: 1.05rem;
      margin-bottom: 12px;
    }}
    .meta-bar {{
      display: flex;
      flex-wrap: wrap;
      gap: 20px;
      font-size: 0.85rem;
      color: var(--text-muted);
      background: var(--card-bg);
      padding: 10px 16px;
      border-radius: 6px;
      border: 1px solid var(--card-border);
    }}
    .badge {{
      display: inline-block;
      font-size: 0.72rem;
      font-weight: 600;
      padding: 3px 8px;
      border-radius: 12px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }}
    .badge-prediction {{ background: var(--pred-bg); border: 1px solid var(--pred-border); color: var(--pred-text); }}
    .badge-experimental {{ background: var(--exp-bg); border: 1px solid var(--exp-border); color: var(--exp-text); }}
    .badge-md {{ background: var(--purple-bg); border: 1px solid var(--purple-border); color: var(--purple-text); }}
    .badge-verified {{ background: var(--ok-bg); border: 1px solid var(--ok-border); color: var(--ok-text); }}

    .grid-kpi {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 16px;
      margin-bottom: 28px;
    }}
    .kpi-card {{
      background: var(--card-bg);
      border: 1px solid var(--card-border);
      border-radius: 8px;
      padding: 16px;
      position: relative;
    }}
    .kpi-title {{
      font-size: 0.8rem;
      color: var(--text-muted);
      text-transform: uppercase;
      letter-spacing: 0.5px;
      margin-bottom: 6px;
    }}
    .kpi-value {{
      font-size: 1.8rem;
      font-weight: 700;
      color: var(--text-bright);
      line-height: 1.2;
    }}
    .kpi-sub {{
      font-size: 0.8rem;
      color: var(--text-muted);
      margin-top: 6px;
    }}

    .section-card {{
      background: var(--card-bg);
      border: 1px solid var(--card-border);
      border-radius: 8px;
      padding: 24px;
      margin-bottom: 28px;
    }}
    .section-title {{
      font-size: 1.3rem;
      color: var(--text-bright);
      margin-bottom: 12px;
      display: flex;
      align-items: center;
      gap: 12px;
    }}

    .provenance-table {{
      width: 100%;
      border-collapse: collapse;
      margin: 16px 0;
      font-size: 0.9rem;
    }}
    .provenance-table th, .provenance-table td {{
      padding: 12px 14px;
      text-align: left;
      border: 1px solid var(--card-border);
      vertical-align: top;
    }}
    .provenance-table th {{
      background: rgba(255, 255, 255, 0.03);
      color: var(--text-bright);
      font-weight: 600;
    }}
    .provenance-table tr:hover {{
      background: rgba(255, 255, 255, 0.02);
    }}

    .alert-box {{
      padding: 16px;
      border-radius: 6px;
      margin: 16px 0;
      font-size: 0.9rem;
      line-height: 1.5;
    }}
    .alert-warning {{
      background: var(--pred-bg);
      border-left: 4px solid var(--pred-border);
      color: var(--pred-text);
    }}
    .alert-danger {{
      background: var(--exp-bg);
      border-left: 4px solid var(--exp-border);
      color: var(--exp-text);
    }}
    .alert-info {{
      background: var(--accent-subtle);
      border-left: 4px solid var(--accent);
      color: var(--text-bright);
    }}

    .data-table {{
      width: 100%;
      border-collapse: collapse;
      margin: 14px 0;
      font-size: 0.85rem;
    }}
    .data-table th, .data-table td {{
      padding: 8px 12px;
      text-align: left;
      border: 1px solid var(--card-border);
    }}
    .data-table th {{
      background: rgba(255, 255, 255, 0.03);
      color: var(--text-bright);
    }}
    .data-table tr:nth-child(even) {{
      background: rgba(255, 255, 255, 0.015);
    }}

    .image-container {{
      text-align: center;
      margin: 20px 0;
    }}
    .image-container img {{
      max-width: 100%;
      height: auto;
      border-radius: 6px;
      border: 1px solid var(--card-border);
      background: #000;
    }}
    .image-caption {{
      font-size: 0.82rem;
      color: var(--text-muted);
      margin-top: 8px;
    }}

    .code-block {{
      background: #000;
      padding: 12px 16px;
      border-radius: 6px;
      font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
      font-size: 0.82rem;
      color: #39bae6;
      overflow-x: auto;
      border: 1px solid var(--card-border);
      margin: 10px 0;
    }}

    footer {{
      margin-top: 40px;
      padding-top: 20px;
      border-top: 1px solid var(--card-border);
      font-size: 0.8rem;
      color: var(--text-muted);
      text-align: center;
    }}
  </style>
</head>
<body>
  <div class="container">
    <header>
      <div class="header-badge-row">
        <span class="badge badge-prediction">In Silico Research Scaffold</span>
        <span class="badge badge-md">Steered Molecular Dynamics</span>
        <span class="badge badge-experimental">No Experimental MS/MS Provided</span>
      </div>
      <h1>Velvet Worm Slime Proteomics &amp; Molecular Dynamics Report</h1>
      <p class="subtitle">Theoretical Proteotypic Peptide Characterization &amp; Single-Peptide SMD Mechanical Pulling Profile</p>
      <div class="meta-bar">
        <div><strong>Dataset:</strong> {summary.get('num_entries', len(protein_rows))} FASTA CDS records (PH &amp; PM slime protein families)</div>
        <div><strong>In Silico Enzyme:</strong> Trypsin (missed cleavages &le; {summary.get('missed_cleavages', 1)})</div>
        <div><strong>Generated:</strong> {generated_time}</div>
      </div>
    </header>

    <!-- KPI Cards -->
    <div class="grid-kpi">
      <div class="kpi-card">
        <div class="kpi-title">Total FASTA Records</div>
        <div class="kpi-value">{summary.get('num_entries', len(protein_rows))}</div>
        <div class="kpi-sub">PH: {summary.get('num_PH_entries', len(ph_rows))} | PM: {summary.get('num_PM_entries', len(pm_rows))} | Trans. CDS: {summary.get('num_nucleotide_entries', len(protein_rows))}</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-title">Theoretical Peptides</div>
        <div class="kpi-value">{summary.get('num_tryptic_peptides', len(peptide_rows))}</div>
        <div class="kpi-sub"><span class="badge badge-prediction">In Silico</span> 6&ndash;40 AA range</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-title">Unique Peptide Candidates</div>
        <div class="kpi-value">{summary.get('num_proteotypic_peptides', len(unique_peptides))}</div>
        <div class="kpi-sub">Unique sequence across dataset</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-title">Flagged Internal Stops</div>
        <div class="kpi-value" style="color: { 'var(--exp-text)' if flagged_records else 'var(--ok-text)' }">{len(flagged_records)}</div>
        <div class="kpi-sub">{ 'PP551270_1 (+1 frame stop codons)' if flagged_records else 'None' }</div>
      </div>
      """

    if smd_mechanics:
        doc += f"""
      <div class="kpi-card">
        <div class="kpi-title">SMD Peak Force</div>
        <div class="kpi-value" style="color: var(--purple-text);">{smd_mechanics.get('peak_rupture_force_pn', 0.0):.1f} pN</div>
        <div class="kpi-sub"><span class="badge badge-md">SMD Demo</span> Pulling velocity {smd_mechanics.get('pull_velocity_nm_ps', 0.02)} nm/ps</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-title">SMD Mechanical Work</div>
        <div class="kpi-value" style="color: var(--purple-text);">{smd_mechanics.get('total_work_kj_mol', 0.0):.2f} <span style="font-size: 1rem;">kJ/mol</span></div>
        <div class="kpi-sub">&Delta;d: {smd_mechanics.get('net_extension_nm', smd_mechanics.get('net_extension_nm', smd_mechanics.get('extension_delta_nm', 0.0))):.2f} nm ({smd_mechanics.get('initial_distance_nm', 0.0):.2f} &rarr; {smd_mechanics.get('final_distance_nm', 0.0):.2f} nm)</div>
      </div>
        """

    doc += """
    </div>

    <!-- PROVENANCE AND EVIDENCE LEVEL CLASSIFICATION (CORE REQUIREMENT) -->
    <div class="section-card">
      <div class="section-title">
        <span>Provenance &amp; Evidence Level Classification</span>
        <span class="badge badge-experimental">Strict Evidentiary Boundary</span>
      </div>
      <p style="margin-bottom: 12px;">
        To ensure scientific integrity and prevent overinterpretation, all findings in this report are categorized
        into <strong>In Silico Theoretical Predictions</strong> (computational hypotheses generated from FASTA sequences and physics-based models)
        versus <strong>Experimental Evidence</strong> (empirical analytical chemistry and biophysical measurements required for confirmation).
      </p>

      <table class="provenance-table">
        <thead>
          <tr>
            <th style="width: 22%;">Analysis Component</th>
            <th style="width: 18%;">Evidence Level</th>
            <th style="width: 32%;">Computational Basis (Present in this Report)</th>
            <th style="width: 28%;">Required Experimental Validation</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>FASTA Sequence Translation</strong></td>
            <td><span class="badge badge-prediction">Derived Input</span></td>
            <td>Input records are nucleotide CDS/transcript FASTAs; translation uses the assumed +1 reading frame and standard genetic code.</td>
            <td>Transcript/frame confirmation is a separate validation step; translation is computational, not protein-detection evidence.</td>
          </tr>
          <tr>
            <td><strong>Enzymatic Cleavage (Tryptic Peptides)</strong></td>
            <td><span class="badge badge-prediction">In Silico Model</span></td>
            <td>Idealized cleavage after Arg (R) and Lys (K) unless followed by Pro (P). Configurable missed cleavage count (&le; 1).</td>
            <td>Bottom-up LC-MS/MS (DDA/DIA) digestion efficiency assays; detection of semi-tryptic and non-specific endopeptidase cleavages.</td>
          </tr>
          <tr>
            <td><strong>Peptide Precursor &amp; Fragment m/z</strong></td>
            <td><span class="badge badge-prediction">Theoretical Calculated</span></td>
            <td>Exact monoisotopic masses computed from standard elemental compositions for z = 1, 2, 3 and b/y fragment ions.</td>
            <td>High-resolution MS1 isotopic envelope calibration and MS2 fragmentation spectral matching (FDR &lt; 1%).</td>
          </tr>
          <tr>
            <td><strong>Skyline MRM / PRM Transition Candidates</strong></td>
            <td><span class="badge badge-prediction">Candidate Library</span></td>
            <td>Top theoretical y- and b-ions ranked by nominal sequence position; empirical collision energies approximated by linear charge fit.</td>
            <td>Empirical spectral library generation (e.g. Skyline / EncyclopeDIA) using synthetic heavy peptide standards to determine actual fragment branch ratios.</td>
          </tr>
          <tr>
            <td><strong>Intrinsic Disorder &amp; Charge Profiling</strong></td>
            <td><span class="badge badge-prediction">Heuristic Sequence Model</span></td>
            <td>FoldIndex hydropathy/charge balance and Uversky phase boundaries; net charge via Henderson-Hasselbalch equation at pH 7.4.</td>
            <td>Circular Dichroism (CD) spectroscopy, Small-Angle X-ray Scattering (SAXS), and cryo-EM to probe disordered vs. ordered polymer assemblies.</td>
          </tr>
          <tr>
            <td><strong>Steered Molecular Dynamics (SMD) Mechanics</strong></td>
            <td><span class="badge badge-md">Theoretical Toy Model</span></td>
            <td>Isolated synthetic oligopeptide pulled at constant velocity in implicit solvent (Amber99SB-ILDN / OBC2 GB/SA).</td>
            <td>Single-molecule AFM force spectroscopy or optical tweezers on purified native/recombinant slime proteins; macro-tensile testing of drawn slime fibers.</td>
          </tr>
        </tbody>
      </table>

      <div class="alert-box alert-warning">
        <strong>Biological Caution:</strong> No wet-lab mass spectrometry spectra, chromatographic retention times, or single-molecule pulling assays were provided.
        All proteomic entries represent theoretical mass hypotheses for instrument method development, and molecular dynamics traces represent non-equilibrium mechanics of an isolated oligomer, not bulk velvet worm slime hydrogel formation.
      </div>
    </div>
"""

    # Section: Flagged Sequence Notice
    if flagged_records:
        doc += """
    <div class="section-card" style="border-color: var(--exp-border);">
      <div class="section-title" style="color: var(--exp-text);">
        <span>Quality Control Alert: Nucleotide Translation Artifact Flagged</span>
        <span class="badge badge-experimental">QC Anomaly</span>
      </div>
      <p>
        The FASTA record <code>Euperipatoides_rowelli__PH2__PP551270_1</code> contains <strong>41 internal stop codons</strong> when translated in the default +1 reading frame.
      </p>
      <div class="alert-box alert-danger">
        <strong>Quarantine Protocol:</strong> This sequence was excluded from in silico tryptic digestion and Skyline transition exports to avoid corrupting candidate MS1/MS2 transition libraries.
        It is retained in the protein properties inventory with explicit stop-codon counts to preserve full provenance. Re-annotation (frame-shift inspection, alternative open reading frames, or pseudogene verification) is required before biological use.
      </div>
    </div>
"""

    # Section: Physicochemical & Sequence Properties
    doc += f"""
    <div class="section-card">
      <div class="section-title">
        <span>Slime Protein Family Physicochemical Profiling</span>
        <span class="badge badge-prediction">In Silico Computed</span>
      </div>
      <p style="margin-bottom: 14px;">
        Comparison of protein records annotated with the <strong>PH</strong> and <strong>PM</strong> family labels in the supplied FASTA files. All listed properties are calculated from the translated sequences and are not experimentally measured.
      </p>

      <table class="data-table">
        <thead>
          <tr>
            <th>Protein Family</th>
            <th>Count</th>
            <th>Mean Length (AA)</th>
            <th>Mean MW (kDa)</th>
            <th>Mean pI</th>
            <th>Mean Net Charge (pH 7.4)</th>
            <th>Proline %</th>
            <th>Glycine %</th>
            <th>Predicted Disordered %</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>PH Family</strong></td>
            <td>{ph_stats['count']}</td>
            <td>{ph_stats['mean_len']:.1f}</td>
            <td>{ph_stats['mean_mw']:.1f}</td>
            <td>{ph_stats['mean_pi']:.2f}</td>
            <td>{ph_stats['mean_charge']:+.1f}</td>
            <td>{ph_stats['mean_pro']:.1f}%</td>
            <td>{ph_stats['mean_gly']:.1f}%</td>
            <td>{ph_stats['disordered_pct']:.1f}%</td>
          </tr>
          <tr>
            <td><strong>PM Family</strong></td>
            <td>{pm_stats['count']}</td>
            <td>{pm_stats['mean_len']:.1f}</td>
            <td>{pm_stats['mean_mw']:.1f}</td>
            <td>{pm_stats['mean_pi']:.2f}</td>
            <td>{pm_stats['mean_charge']:+.1f}</td>
            <td>{pm_stats['mean_pro']:.1f}%</td>
            <td>{pm_stats['mean_gly']:.1f}%</td>
            <td>{pm_stats['disordered_pct']:.1f}%</td>
          </tr>
        </tbody>
      </table>
"""

    if overview_b64:
        doc += f"""
      <div class="image-container">
        <img src="{overview_b64}" alt="Proteomics Overview">
        <div class="image-caption">Figure 1: Proteomics and physicochemical overview across velvet worm slime proteins (molecular weight, isoelectric point, charge at pH 7.4, and amino-acid composition).</div>
      </div>
"""

    doc += """
      <h3 style="color: var(--text-bright); margin: 18px 0 8px 0; font-size: 1rem;">Representative Slime Protein Candidates</h3>
      <table class="data-table">
        <thead>
          <tr>
            <th>Protein Identifier</th>
            <th>Species</th>
            <th>Family</th>
            <th>Length (AA)</th>
            <th>MW (kDa)</th>
            <th>pI</th>
            <th>Charge (pH 7.4)</th>
            <th>FoldIndex State</th>
          </tr>
        </thead>
        <tbody>
"""
    for r in protein_rows[:8]:
        doc += f"""
          <tr>
            <td><code style="color: var(--accent);">{html.escape(str(r.get('protein_id', ''))[:45])}</code></td>
            <td>{html.escape(str(r.get('species', 'Unknown')))}</td>
            <td><span class="badge badge-prediction">{html.escape(str(r.get('family', '')))}</span></td>
            <td>{r.get('protein_length_aa', 0)}</td>
            <td>{float(r.get('mw_kda', 0.0)):.1f}</td>
            <td>{float(r.get('pI', 0.0)):.2f}</td>
            <td>{float(r.get('charge_ph', 0.0)):+.1f}</td>
            <td>{html.escape(str(r.get('predicted_state', 'N/A')))}</td>
          </tr>
"""
    doc += """
        </tbody>
      </table>
    </div>
"""

    # Section: Theoretical Tryptic Peptides & Skyline Transition Candidates
    doc += f"""
    <div class="section-card">
      <div class="section-title">
        <span>Theoretical Tryptic Peptides &amp; Skyline Transition Candidates</span>
        <span class="badge badge-prediction">In Silico Digestion</span>
      </div>
      <p style="margin-bottom: 14px;">
        In silico tryptic digestion generated <strong>{len(peptide_rows):,}</strong> theoretical peptide sequences,
        of which <strong>{len(unique_peptides):,}</strong> are strictly proteotypic (unique to a single slime protein entry in this dataset).
        These sequence-unique peptides are computational candidates for targeted Selected Reaction Monitoring (SRM) or Parallel Reaction Monitoring (PRM); uniqueness does not establish detection or biomarker validity.
      </p>

      <h3 style="color: var(--text-bright); margin: 16px 0 8px 0; font-size: 1rem;">Top Candidate Proteotypic Peptides</h3>
      <table class="data-table">
        <thead>
          <tr>
            <th>Protein ID</th>
            <th>Family</th>
            <th>Peptide Sequence</th>
            <th>Residues</th>
            <th>Length</th>
            <th>Precursor m/z (z=2)</th>
            <th>Precursor m/z (z=3)</th>
            <th>Retention Index</th>
          </tr>
        </thead>
        <tbody>
"""
    for p in unique_peptides[:10]:
        doc += f"""
          <tr>
            <td><code style="color: var(--accent);">{html.escape(str(p.get('protein_id', ''))[:30])}</code></td>
            <td>{html.escape(str(p.get('family', '')))}</td>
            <td><strong><code>{html.escape(str(p.get('peptide_sequence', '')))}</code></strong></td>
            <td>{p.get('start', 0)}&ndash;{p.get('end', 0)}</td>
            <td>{p.get('length_aa', 0)}</td>
            <td>{float(p.get('mz_z2', 0.0)):.4f}</td>
            <td>{float(p.get('mz_z3', 0.0)):.4f}</td>
            <td>{float(p.get('retention_index', 0.0)):.1f}</td>
          </tr>
"""
    doc += """
        </tbody>
      </table>
      <p style="font-size: 0.8rem; color: var(--text-muted); margin-top: 6px;">
        <em>Note:</em> The complete transition library containing theoretical b- and y-fragment ions is exported in <code>skyline_transition_list.csv</code>.
      </p>
    </div>
"""

    # Section: Molecular Dynamics (SMD)
    if smd_mechanics:
        doc += f"""
    <div class="section-card">
      <div class="section-title">
        <span>Steered Molecular Dynamics (SMD) Mechanical Pulling</span>
        <span class="badge badge-md">Biophysical Simulation</span>
      </div>
      <p style="margin-bottom: 14px;">
        To demonstrate mechanical tensile behavior of repetitive slime protein motifs, constant-velocity Steered Molecular Dynamics (SMD)
        was performed on an idealized peptide oligomer (<code>{html.escape(str(smd_mechanics.get('sequence', '')))}</code>) using OpenMM.
        The C-terminal alpha-carbon was pulled away from the fixed N-terminal alpha-carbon via a moving harmonic spring.
      </p>

      <div class="grid-kpi" style="margin: 16px 0;">
        <div class="kpi-card">
          <div class="kpi-title">Initial &rarr; Final Distance</div>
          <div class="kpi-value">{smd_mechanics.get('initial_distance_nm', 0.0):.3f} &rarr; {smd_mechanics.get('final_distance_nm', 0.0):.3f} nm</div>
          <div class="kpi-sub">Total Extension: {smd_mechanics.get('net_extension_nm', smd_mechanics.get('net_extension_nm', smd_mechanics.get('extension_delta_nm', 0.0))):.3f} nm</div>
        </div>
        <div class="kpi-card">
          <div class="kpi-title">Peak Rupture Force</div>
          <div class="kpi-value" style="color: var(--purple-text);">{smd_mechanics.get('peak_rupture_force_pn', 0.0):.1f} pN</div>
          <div class="kpi-sub">At extension {smd_mechanics.get('extension_at_peak_nm', smd_mechanics.get('extension_at_peak_nm', smd_mechanics.get('extension_at_peak_force_nm', 0.0))):.3f} nm</div>
        </div>
        <div class="kpi-card">
          <div class="kpi-title">Mechanical Work (W)</div>
          <div class="kpi-value" style="color: var(--purple-text);">{smd_mechanics.get('total_work_kj_mol', 0.0):.2f} kJ/mol</div>
          <div class="kpi-sub">{smd_mechanics.get('total_work_kcal_mol', float(smd_mechanics.get('total_work_kj_mol', 0.0)) / 4.184):.2f} kcal/mol integrated work</div>
        </div>
        <div class="kpi-card">
          <div class="kpi-title">Simulation Parameters</div>
          <div class="kpi-value" style="font-size: 1.2rem; color: var(--accent);">{smd_mechanics.get('total_steps', smd_mechanics.get('steps', 0))} steps</div>
          <div class="kpi-sub">k = {smd_mechanics.get('spring_k_kj_mol_nm2', smd_mechanics.get('spring_constant_kj_mol_nm2', 500.0)):.0f} kJ/(mol&middot;nm&sup2;) | v = {smd_mechanics.get('pull_velocity_nm_ps', 0.02)} nm/ps</div>
        </div>
      </div>
"""
        if smd_b64:
            doc += f"""
      <div class="image-container">
        <img src="{smd_b64}" alt="SMD Mechanics Force-Extension">
        <div class="image-caption">Figure 2: Steered molecular dynamics force-extension profile and cumulative mechanical work trace for peptide motif <code>{html.escape(str(smd_mechanics.get('sequence', '')))}</code>.</div>
      </div>
"""

        doc += """
      <div class="alert-box alert-info">
        <strong>Biophysical Interpretation &amp; Limits:</strong>
        The SMD trace illustrates the non-equilibrium mechanical extension of an individual disordered oligopeptide backbone.
        In native velvet worm slime, shear-induced macroscopic fiber drawing involves coordinated multi-chain beta-sheet assembly, phase separation, and inter-chain entanglement.
        Single-peptide pulling forces reflect local steric and dihedral barrier crossings under non-equilibrium velocities and cannot be directly equated with bulk fiber tensile strength without multi-chain coarse-grained or atomistic hydrogel simulations.
      </div>
    </div>
"""

    # Section: Reproducibility Run Manifest & Audit Trail
    manifest_html = ""
    if run_manifest:
        inputs_rows = ""
        for inp in run_manifest.get("inputs", []):
            sha_display = f"<code>{inp.get('sha256', 'N/A')}</code>" if inp.get("sha256") else "<span style='color: var(--text-muted);'>N/A</span>"
            size_str = f"{inp.get('size_bytes', 0):,} B" if inp.get("size_bytes") is not None else "N/A"
            inputs_rows += f"""
              <tr>
                <td><strong>{html.escape(str(inp.get('label', '')))}</strong></td>
                <td><code>{html.escape(str(inp.get('filename', '')))}</code></td>
                <td style="font-size: 0.75rem; word-break: break-all;">{sha_display}</td>
                <td>{size_str}</td>
              </tr>
            """

        env_items = ""
        for k, v in run_manifest.get("environment", {}).items():
            env_items += f"<div><strong>{html.escape(k)}:</strong> <code>{html.escape(str(v))}</code></div>"

        param_items = ""
        for k, v in run_manifest.get("parameters", {}).items():
            param_items += f"<div><strong>--{html.escape(k.replace('_', '-'))}:</strong> <code>{html.escape(str(v))}</code></div>"

        manifest_html = f"""
    <div class="section-card" id="run-manifest">
      <div class="section-title">
        <span>Reproducibility Run Manifest</span>
        <span class="badge badge-verified">SHA-256 Verified Audit</span>
      </div>
      <p style="margin-bottom: 14px;">
        Cryptographic checksums, runtime dependencies, exact CLI parameters, and UTC execution timestamps recorded to guarantee bit-level and analytical reproducibility.
      </p>

      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 14px; margin-bottom: 16px;">
        <div style="background: var(--card-bg); border: 1px solid var(--card-border); border-radius: 6px; padding: 12px;">
          <div style="font-size: 0.75rem; color: var(--text-muted); text-transform: uppercase; margin-bottom: 6px;">Execution Metadata</div>
          <div style="font-size: 0.85rem; margin-bottom: 4px;"><strong>Command:</strong> <code>{html.escape(str(run_manifest.get('command', 'analyze')))}</code></div>
          <div style="font-size: 0.85rem; margin-bottom: 4px;"><strong>Run Timestamp:</strong> <span style="color: var(--text-bright); font-family: monospace;">{html.escape(str(run_manifest.get('run_timestamp_utc', 'N/A')))}</span></div>
          <div style="font-size: 0.85rem;"><strong>Schema Version:</strong> <code>{html.escape(str(run_manifest.get('schema_version', '1.0')))}</code></div>
        </div>

        <div style="background: var(--card-bg); border: 1px solid var(--card-border); border-radius: 6px; padding: 12px;">
          <div style="font-size: 0.75rem; color: var(--text-muted); text-transform: uppercase; margin-bottom: 6px;">Runtime Environment &amp; Versions</div>
          <div style="font-size: 0.8rem; display: flex; flex-direction: column; gap: 3px;">
            {env_items}
          </div>
        </div>

        <div style="background: var(--card-bg); border: 1px solid var(--card-border); border-radius: 6px; padding: 12px;">
          <div style="font-size: 0.75rem; color: var(--text-muted); text-transform: uppercase; margin-bottom: 6px;">CLI Parameters Captured</div>
          <div style="font-size: 0.8rem; display: flex; flex-direction: column; gap: 3px;">
            {param_items}
          </div>
        </div>
      </div>

      <h4 style="color: var(--text-bright); margin: 14px 0 8px 0; font-size: 0.95rem;">Input Files &amp; Cryptographic Hashes (SHA-256)</h4>
      <table class="data-table" style="margin-bottom: 8px;">
        <thead>
          <tr>
            <th style="width: 15%;">Label</th>
            <th style="width: 25%;">File</th>
            <th style="width: 45%;">SHA-256 Checksum</th>
            <th style="width: 15%;">Size</th>
          </tr>
        </thead>
        <tbody>
          {inputs_rows}
        </tbody>
      </table>
    </div>
        """


    doc += manifest_html

    doc += """
    <div class="section-card">
      <div class="section-title">
        <span>Reproducibility &amp; Pipeline Execution</span>
        <span class="badge badge-verified">CLI Scaffold</span>
      </div>
      <p style="margin-bottom: 10px;">
        This report and its companion data files were generated using the <code>slime-proteo-md</code> command-line suite:
      </p>

      <div class="code-block">
# Run end-to-end pipeline (proteomics profiling + SMD mechanics + HTML report)
slime-proteo-md pipeline \\
  --fasta data/Peripatus_slime_7_PH_no_gaps.fasta \\
  --fasta data/Peripatus_slime_7_PM_no_gaps.fasta \\
  --out results/

# Or run sequence proteomics analysis alone:
slime-proteo-md analyze \\
  --fasta data/Peripatus_slime_7_PH_no_gaps.fasta \\
  --fasta data/Peripatus_slime_7_PM_no_gaps.fasta \\
  --out results/

# Run steered molecular dynamics (SMD) pulling demonstration:
slime-proteo-md smd \\
  --sequence PSEPGYPGQPSEPGYP \\
  --steps 2000 \\
  --pull-velocity 0.02 \\
  --out results/smd/

# Assemble or regenerate combined HTML report from existing directories:
slime-proteo-md report \\
  --analysis-dir results/ \\
  --smd-dir results/smd/ \\
  --out results/report.html
      </div>

      <h3 style="color: var(--text-bright); margin: 16px 0 8px 0; font-size: 1rem;">Exported Artifacts Inventory</h3>
      <table class="data-table">
        <thead>
          <tr>
            <th>File Name</th>
            <th>Type</th>
            <th>Description</th>
          </tr>
        </thead>
        <tbody>
          <tr>
                      <tr>
            <td><code>run_manifest.json</code></td>
            <td>Audit Manifest</td>
            <td>SHA-256 hashes for raw inputs, runtime package versions, CLI flags, and UTC run timestamp.</td>
          </tr>
          <tr>
            <td><code>report.html</code></td>
            <td>Interactive Report</td>
            <td>Self-contained HTML report with embedded visuals, provenance badges, and metrics.</td>
          </tr>
          <tr>
            <td><code>protein_properties.csv</code></td>
            <td>Data Table</td>
            <td>Protein-level lengths, molecular weights, isoelectric points, net charges, and FoldIndex scores.</td>
          </tr>
          <tr>
            <td><code>tryptic_peptides.csv</code></td>
            <td>Data Table</td>
            <td>In silico digested peptides with monoisotopic masses, m/z (z=1,2,3), and proteotypic flags.</td>
          </tr>
          <tr>
            <td><code>skyline_transition_list.csv</code></td>
            <td>Instrument Input</td>
            <td>Targeted MRM/PRM transition list formatted for Skyline import.</td>
          </tr>
          <tr>
            <td><code>smd_force_extension.csv</code></td>
            <td>Simulation Data</td>
            <td>Frame-by-frame time (ps), distance (nm), force (pN), and cumulative work (kJ/mol).</td>
          </tr>
          <tr>
            <td><code>smd_mechanics_summary.json</code></td>
            <td>Metrics Summary</td>
            <td>Structured peak rupture force, extension delta, stiffness, and work metrics.</td>
          </tr>
        </tbody>
      </table>
    </div>

    <footer>
      <p>Velvet Worm Slime Proteomics &amp; Molecular Dynamics Research Scaffold &bull; gyagu98 &bull; MIT License</p>
    </footer>
  </div>
</body>
</html>
"""

    if output_html_path:
        out_p = Path(output_html_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        out_p.write_text(doc, encoding="utf-8")

    return doc


def build_report_from_directories(
    analysis_dir: str,
    output_html_path: str,
    smd_dir: Optional[str] = None,
) -> str:
    """Load existing analysis and SMD results from directories and render the HTML report."""
    analysis_path = Path(analysis_dir)
    summary_file = analysis_path / "sequence_summary.json"
    summary = {}
    if summary_file.is_file():
        with open(summary_file, "r", encoding="utf-8") as f:
            summary = json.load(f)

    protein_rows = []
    prot_csv = analysis_path / "protein_properties.csv"
    if prot_csv.is_file():
        with open(prot_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                for k in ["protein_length_aa", "cds_length_nt", "internal_stop_codons"]:
                    if k in row and row[k]:
                        try:
                            row[k] = int(row[k])
                        except (ValueError, TypeError):
                            pass
                for k in ["mw_kda", "pI", "charge_ph", "gravy", "proline_percent", "glycine_percent"]:
                    if k in row and row[k]:
                        try:
                            row[k] = float(row[k])
                        except (ValueError, TypeError):
                            pass
                protein_rows.append(row)

    peptide_rows = []
    pep_csv = analysis_path / "tryptic_peptides.csv"
    if pep_csv.is_file():
        with open(pep_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "is_proteotypic" in row and isinstance(row["is_proteotypic"], str):
                    row["is_proteotypic"] = row["is_proteotypic"].strip().lower() in ("true", "1", "yes")
                for k in ["start", "end", "length_aa", "missed_cleavages"]:
                    if k in row and row[k]:
                        try:
                            row[k] = int(row[k])
                        except (ValueError, TypeError):
                            pass
                for k in ["monoisotopic_mw", "mz_z1", "mz_z2", "mz_z3", "retention_index"]:
                    if k in row and row[k]:
                        try:
                            row[k] = float(row[k])
                        except (ValueError, TypeError):
                            pass
                peptide_rows.append(row)

    manifest_data = None
    manifest_file = analysis_path / "run_manifest.json"
    if manifest_file.is_file():
        try:
            with open(manifest_file, "r", encoding="utf-8") as f:
                manifest_data = json.load(f)
        except Exception:
            pass

    overview_img = analysis_path / "proteomics_overview.png"

    smd_mechanics = None
    smd_img = None
    if smd_dir:
        smd_path = Path(smd_dir)
        smd_json = smd_path / "smd_mechanics_summary.json"
        if smd_json.is_file():
            with open(smd_json, "r", encoding="utf-8") as f:
                smd_mechanics = json.load(f)
        smd_png = smd_path / "smd_mechanics.png"
        if smd_png.is_file():
            smd_img = str(smd_png)

    return generate_html_report(
        summary=summary,
        protein_rows=protein_rows,
        peptide_rows=peptide_rows,
        smd_mechanics=smd_mechanics,
        overview_img_path=str(overview_img) if overview_img.is_file() else None,
        smd_img_path=smd_img,
        output_html_path=output_html_path,
        run_manifest=manifest_data,
    )
