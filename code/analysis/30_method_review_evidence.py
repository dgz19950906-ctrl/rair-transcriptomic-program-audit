#!/usr/bin/env python3
"""Deterministic evidence extraction for the RAIR methodological self-review."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path("/home/dony/ThyroidCancer_Project/rair_audit")
M = ROOT / "manuscript/two_null_v2/MANUSCRIPT_FULL_DRAFT.two_null_v2.md"
REG = ROOT / "covariance_null_v2/inputs/frozen_programs_all9.tsv"
BOOT = ROOT / "results/08_three_atlas_bootstrap_effects.tsv"
OUT = ROOT / "pre_submission_audit_v1/method_review/evidence_v1"


def tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    text = M.read_text(encoding="utf-8")
    registry = tsv(REG)
    programs = {}
    for row in registry:
        programs.setdefault(row["signature_id"], set()).add(row["gene"])

    marker_patterns = {
        "T_NK": r"T/NK cells \((.*?)\); B cells",
        "B": r"B cells \((.*?)\); thyroid cells",
        "Thyroid": r"thyroid cells \((.*?)\); myeloid cells",
        "Myeloid": r"myeloid cells \((.*?)\); fibroblasts",
        "Fibroblast": r"fibroblasts \((.*?)\); and endothelial cells",
        "Endothelial": r"endothelial cells \((.*?)\)\. Marker expression",
    }
    markers = {}
    for name, pat in marker_patterns.items():
        match = re.search(pat, text, flags=re.S | re.I)
        if not match:
            raise RuntimeError(f"Could not extract {name} marker panel")
        markers[name] = set(re.findall(r"`([^`]+)`", match.group(1)))

    target = {
        "TDS_16": "Thyroid",
        "IODIDE_HANDLING_11": "Thyroid",
        "CONDELLO_2025_SIX": "T_NK",
    }
    overlap_rows = []
    for sid, compartment in target.items():
        overlap = sorted(programs[sid] & markers[compartment])
        overlap_rows.append({
            "signature_id": sid,
            "target_compartment": compartment,
            "program_gene_n": len(programs[sid]),
            "marker_gene_n": len(markers[compartment]),
            "overlap_n": len(overlap),
            "program_overlap_fraction": len(overlap) / len(programs[sid]),
            "marker_overlap_fraction": len(overlap) / len(markers[compartment]),
            "overlap_genes": ";".join(overlap),
        })
    with (OUT / "PROGRAM_MARKER_OVERLAP.tsv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(overlap_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(overlap_rows)

    boot = [r for r in tsv(BOOT) if r["signature_id"] in target]
    donor_rows = [{
        "dataset": r["dataset"], "signature_id": r["signature_id"],
        "target_compartment": r["target_compartment"], "n_donors": int(r["n_donors"]),
        "estimate": float(r["estimate"]), "conf_low": float(r["conf_low"]),
        "conf_high": float(r["conf_high"]),
    } for r in boot]
    with (OUT / "FIG5_PRIMARY_DONOR_PRECISION.tsv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(donor_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(donor_rows)

    literature_section = text.split("## Systematic evidence and reproducibility audit", 1)[1].split("## Public transcriptomic datasets", 1)[0]
    reviewer_tokens = re.findall(r"\b(?:reviewer|screened independently|dual screening|adjudicat\w*|consensus)\b", literature_section, flags=re.I)
    placeholders = []
    for i, line in enumerate(text.splitlines(), 1):
        if re.search(r"\[To be completed|to be added|repository name|persistent identifier|\bTBD\b|\bTODO\b|\[VERIFY", line, flags=re.I):
            placeholders.append({"line": i, "text": line})

    evidence = {
        "program_marker_overlap": overlap_rows,
        "gse191288_primary_donors": [r for r in donor_rows if r["dataset"] == "GSE191288"],
        "literature_audit_reviewer_protocol_tokens": reviewer_tokens,
        "literature_audit_reviewer_protocol_reported": bool(reviewer_tokens),
        "ethics_heading_count": len(re.findall(r"^## Ethics statement\s*$", text, flags=re.M)),
        "submission_placeholders": placeholders,
        "marker_overlap_interpretation": {
            "TDS_16": "annotation-program non-independence requires leave-program-genes-out sensitivity",
            "IODIDE_HANDLING_11": "annotation-program non-independence requires leave-program-genes-out sensitivity",
            "CONDELLO_2025_SIX": "no direct overlap with T/NK marker panel",
        },
        "bootstrap_precision_note": "GSE191288 has n=3 paired donors for all three primary programs; percentile bootstrap intervals are descriptive and have limited support.",
    }
    (OUT / "METHOD_REVIEW_EVIDENCE.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
