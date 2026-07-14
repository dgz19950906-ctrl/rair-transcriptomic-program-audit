#!/usr/bin/env python3
"""Write SHA-256 hashes for the frozen systematic-search evidence set."""

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "SYSTEMATIC_LITERATURE_SEARCH_FREEZE.md",
    ROOT / "LITERATURE_AUDIT_RESULTS_DRAFT.md",
    ROOT / "literature_search_freeze_template.tsv",
    ROOT / "literature_exclusion_template.tsv",
    ROOT / "literature_audit_funnel_template.tsv",
    ROOT / "literature_search" / "pubmed_query.txt",
    ROOT / "literature_search" / "pubmed_search_metadata.json",
    ROOT / "literature_search" / "pubmed_search_records.tsv",
    ROOT / "literature_search" / "europe_pmc_search_records.tsv",
    ROOT / "literature_search" / "semantic_scholar_search_records.tsv",
    ROOT / "literature_search" / "geo_search_records.tsv",
    ROOT / "literature_search" / "supplementary_search_metadata.json",
    ROOT / "literature_search" / "master_screening.tsv",
    ROOT / "literature_search" / "pubmed_screening.tsv",
    ROOT / "literature_search" / "literature_tool_audit.tsv",
    ROOT / "literature_search" / "literature_flow_counts.json",
    ROOT / "literature_search" / "fulltext_xml_manifest.tsv",
    ROOT / "references" / "library.bib",
    ROOT / "references" / "snowball_candidates.bib",
]

lines = ["sha256\tbytes\tpath"]
for path in FILES:
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    lines.append(f"{digest}\t{len(content)}\t{path.relative_to(ROOT)}")

out = ROOT / "literature_search" / "FROZEN_SHA256.tsv"
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"Wrote {len(FILES)} hashes to {out}")
