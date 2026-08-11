#!/usr/bin/env python3
"""Thin wrapper around Pensieve's real PAML marginal-ASR parser.

There is deliberately no independent node-count logic here.  The validator and
03_integrate_asr_evidence.py use the same parser and the same source of truth:
PAML's own ``Nodes X to Y are ancestral`` declaration in ``rst``.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path


def load_integrator(script_dir: Path):
    path = script_dir / "03_integrate_asr_evidence.py"
    spec = importlib.util.spec_from_file_location("pensieve_integrator", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not load PAML parser: {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_status(path: Path, row: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main():
    ap = argparse.ArgumentParser(description="Validate codeml marginal ASR using the rst-declared ancestral-node range.")
    ap.add_argument("--rst", required=True)
    ap.add_argument("--phylip", required=True)
    ap.add_argument("--status-out", required=True)
    args = ap.parse_args()

    rst = Path(args.rst)
    phylip = Path(args.phylip)
    status_out = Path(args.status_out)
    mod = load_integrator(Path(__file__).resolve().parent)

    default_row = {
        "marginal_asr_valid": False,
        "reason": "validation_not_started",
        "tip_count": "NA",
        "declared_ancestral_start": "NA",
        "declared_ancestral_end": "NA",
        "expected_internal_sequences": "NA",
        "observed_internal_sequences": "NA",
        "alignment_length": "NA",
        "joint_section_present": False,
        "overall_accuracy_marker_present": False,
        "validation_basis": "rst_line:Nodes_X_to_Y_are_ancestral",
    }

    try:
        row, _marginal, _diagnostics, declaration = mod.validate_paml_marginal_output(rst, phylip)
        write_status(status_out, row)
        print(
            f"Validated complete marginal ASR from PAML declaration: nodes "
            f"{declaration['start']}..{declaration['end']} ({declaration['count']} sequences); "
            f"alignment_length={row['alignment_length']}"
        )
        return 0
    except (SystemExit, Exception) as exc:
        row = dict(default_row)
        row["reason"] = str(exc) or exc.__class__.__name__
        try:
            if rst.exists() and rst.stat().st_size:
                text = rst.read_text(errors="replace")
                row["joint_section_present"] = "Joint reconstruction of ancestral sequences" in text
                row["overall_accuracy_marker_present"] = "Overall accuracy of the" in text
                try:
                    declaration = mod.parse_declared_ancestral_node_range(rst)
                    row["declared_ancestral_start"] = declaration["start"]
                    row["declared_ancestral_end"] = declaration["end"]
                    row["expected_internal_sequences"] = declaration["count"]
                except BaseException:
                    pass
            if phylip.exists() and phylip.stat().st_size:
                tips, length = mod.read_phylip_dimensions(phylip)
                row["tip_count"] = tips
                row["alignment_length"] = length
        except Exception:
            pass
        write_status(status_out, row)
        print(f"Marginal ASR validation failed: {row['reason']}", flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
