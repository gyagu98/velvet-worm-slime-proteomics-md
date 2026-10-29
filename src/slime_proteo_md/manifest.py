"""Reproducibility run manifest generation and auditing for slime proteomics & MD."""

from __future__ import annotations

import hashlib
import platform
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


def compute_sha256(path: Path | str) -> str:
    """Compute the SHA-256 hex digest for a file. Returns empty string if missing."""
    p = Path(path)
    if not p.is_file():
        return ""
    h = hashlib.sha256()
    with p.open("rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def get_tool_versions(extra_packages: Optional[Iterable[str]] = None) -> Dict[str, str]:
    """Capture runtime versions of Python, OS, platform, and key libraries."""
    packages = [
        "slime-proteo-md",
        "openmm",
        "biopython",
        "pandas",
        "numpy",
        "matplotlib",
        "click",
        "pytest",
    ]
    if extra_packages:
        packages.extend(extra_packages)

    versions: Dict[str, str] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "architecture": " ".join(platform.architecture()),
    }
    for pkg in packages:
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            versions[pkg] = "not-installed"
        except Exception:
            versions[pkg] = "unknown"

    return versions


def build_run_manifest(
    *,
    command: str,
    parameters: Dict[str, Any],
    input_files: Dict[str, Path | str],
    output_files: Optional[Dict[str, Path | str]] = None,
    extra_metadata: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build a serializable reproducibility run manifest."""
    timestamp_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    resolved_inputs: List[Dict[str, Any]] = []
    for label, path in input_files.items():
        p = Path(path)
        sha = compute_sha256(p)
        resolved_inputs.append({
            "label": label,
            "filename": p.name,
            "path": str(p),
            "exists": p.is_file(),
            "size_bytes": p.stat().st_size if p.is_file() else None,
            "sha256": sha or None,
        })

    resolved_outputs: List[Dict[str, Any]] = []
    if output_files:
        for label, path in output_files.items():
            p = Path(path)
            sha = compute_sha256(p) if p.is_file() else None
            resolved_outputs.append({
                "label": label,
                "filename": p.name,
                "path": str(p),
                "exists": p.is_file(),
                "size_bytes": p.stat().st_size if p.is_file() else None,
                "sha256": sha,
            })

    manifest: Dict[str, Any] = {
        "schema_version": "1.0",
        "command": command,
        "run_timestamp_utc": timestamp_utc,
        "parameters": parameters,
        "environment": get_tool_versions(),
        "inputs": resolved_inputs,
        "outputs": resolved_outputs,
    }
    if extra_metadata:
        manifest["extra_metadata"] = extra_metadata

    return manifest
