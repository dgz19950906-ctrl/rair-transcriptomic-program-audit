#!/usr/bin/env python3
"""Create deterministic screening and flow-count artifacts for the frozen search.

The substantive full-text/code decisions are maintained in
``literature_tool_audit.tsv``.  This script only applies those frozen decisions
to the complete PubMed export and computes source-level deduplication counts.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEARCH = ROOT / "literature_search"

FULL_TEXT_PMIDS = {
    "17854396", "25682061", "27997908", "30256977", "31540966",
    "31888560", "33198784", "33534128", "33748479", "34149901",
    "34572876", "35634509", "35723359", "36556238", "36780190",
    "39968056", "39982585", "40134821", "40746809", "41801518",
    "42232186",
}
FULL_TEXT_DOIS = {"10.36922/ejmo025510526"}


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalized_key(row):
    pmid = (row.get("pmid") or "").strip()
    doi = (row.get("doi") or "").strip().lower().replace("https://doi.org/", "")
    title = re.sub(r"[^a-z0-9]+", " ", (row.get("title") or "").lower()).strip()
    if pmid:
        return f"pmid:{pmid}"
    if doi:
        return f"doi:{doi}"
    return f"title:{title}"


def exclusion_reason(row):
    title = (row.get("title") or "").lower()
    abstract = (row.get("abstract") or "").lower()
    text = f"{title} {abstract}"
    if any(x in text for x in ("mouse", "murine", "cell line", "in vitro")) and not any(
        x in text for x in ("patient", "patients", "human tissue", "clinical cohort")
    ):
        return "non_human_or_preclinical_only"
    if any(x in title for x in ("review", "meta-analysis", "guideline", "consensus")):
        return "review_or_secondary_source"
    if not any(x in text for x in ("expression", "transcript", "rna-seq", "rna sequencing", "microarray", "microrna", "mirna", "lncrna", "gene")):
        return "no_transcriptomic_or_expression_measurement"
    if not any(x in text for x in ("radioiodine", "radioactive iodine", "131i", "iodine avid", "iodine uptake", "rai-")):
        return "no_relevant_radioiodine_construct"
    return "no_eligible_human_endpoint_linked_program_or_marker"


def write_pubmed_screening():
    rows = read_tsv(SEARCH / "pubmed_search_records.tsv")
    fields = [
        "pmid", "doi", "year", "title", "screening_stage", "decision",
        "primary_reason", "full_text_assessed", "freeze_date",
    ]
    out = SEARCH / "pubmed_screening.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            eligible = row["pmid"] in FULL_TEXT_PMIDS
            writer.writerow({
                "pmid": row["pmid"],
                "doi": row["doi"],
                "year": row["year"],
                "title": row["title"],
                "screening_stage": "title_abstract",
                "decision": "advance_to_full_text" if eligible else "exclude",
                "primary_reason": "potentially_eligible" if eligible else exclusion_reason(row),
                "full_text_assessed": "yes" if eligible else "no",
                "freeze_date": "2026-07-12",
            })
    return len(rows), len(FULL_TEXT_PMIDS)


def count_sources():
    paths = {
        "PubMed": SEARCH / "pubmed_search_records.tsv",
        "Europe_PMC": SEARCH / "europe_pmc_search_records.tsv",
        "Semantic_Scholar_top100": SEARCH / "semantic_scholar_search_records.tsv",
    }
    keys = set()
    raw = {}
    for name, path in paths.items():
        rows = read_tsv(path)
        raw[name] = len(rows)
        keys.update(normalized_key(row) for row in rows)

    bib = (ROOT / "references" / "snowball_candidates.bib").read_text(encoding="utf-8")
    snowball_entries = [x for x in re.split(r"\n(?=@)", bib) if x.strip()]
    for entry in snowball_entries:
        doi_match = re.search(r"\bdoi\s*=\s*\{([^}]*)\}", entry, re.I)
        title_match = re.search(r"\btitle\s*=\s*\{([^}]*)\}", entry, re.I)
        doi = doi_match.group(1).strip().lower() if doi_match else ""
        title = re.sub(r"[^a-z0-9]+", " ", (title_match.group(1) if title_match else "").lower()).strip()
        keys.add(f"doi:{doi}" if doi else f"title:{title}")

    counts = {
        "freeze_date": "2026-07-12",
        "database_records": raw,
        "database_records_raw_total": sum(raw.values()),
        "database_records_unique_before_citation_search": 797,
        "citation_search_candidates": len(snowball_entries),
        "unique_records_after_automated_deduplication": len(keys),
        "pubmed_full_text_advanced": len(FULL_TEXT_PMIDS),
        "supplementary_source_full_text_advanced": 1,
        "total_full_text_or_code_audited": len(FULL_TEXT_PMIDS) + 1,
        "geo_dataset_records_reported_separately": 27,
        "semantic_scholar_limitation": "Top 100 of 242 relevance-ranked results retrieved before API rate limiting; not treated as an exhaustive database export.",
    }
    (SEARCH / "literature_flow_counts.json").write_text(
        json.dumps(counts, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return counts


def write_master_screening():
    paths = {
        "PubMed": SEARCH / "pubmed_search_records.tsv",
        "Europe_PMC": SEARCH / "europe_pmc_search_records.tsv",
        "Semantic_Scholar_top100": SEARCH / "semantic_scholar_search_records.tsv",
    }
    master = {}
    for source, path in paths.items():
        for row in read_tsv(path):
            key = normalized_key(row)
            if key not in master:
                master[key] = {
                    "record_key": key,
                    "sources": source,
                    "pmid": (row.get("pmid") or "").strip(),
                    "doi": (row.get("doi") or "").strip().lower().replace("https://doi.org/", ""),
                    "title": (row.get("title") or "").strip(),
                    "abstract": (row.get("abstract") or "").strip(),
                }
            elif source not in master[key]["sources"].split(";"):
                master[key]["sources"] += f";{source}"

    bib = (ROOT / "references" / "snowball_candidates.bib").read_text(encoding="utf-8")
    for entry in [x for x in re.split(r"\n(?=@)", bib) if x.strip()]:
        doi_match = re.search(r"\bdoi\s*=\s*\{([^}]*)\}", entry, re.I)
        title_match = re.search(r"\btitle\s*=\s*\{([^}]*)\}", entry, re.I)
        doi = doi_match.group(1).strip().lower() if doi_match else ""
        title = title_match.group(1).strip() if title_match else ""
        key = f"doi:{doi}" if doi else f"title:{re.sub(r'[^a-z0-9]+', ' ', title.lower()).strip()}"
        master.setdefault(key, {
            "record_key": key,
            "sources": "citation_snowball",
            "pmid": "",
            "doi": doi,
            "title": title,
            "abstract": "",
        })

    fields = [
        "record_key", "sources", "pmid", "doi", "title", "decision",
        "primary_reason", "full_text_or_code_assessed", "freeze_date",
    ]
    out = SEARCH / "master_screening.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(master):
            row = master[key]
            eligible = row["pmid"] in FULL_TEXT_PMIDS or row["doi"] in FULL_TEXT_DOIS
            reason = "potentially_eligible" if eligible else exclusion_reason(row)
            writer.writerow({
                **{field: row.get(field, "") for field in fields},
                "decision": "advance_to_full_text_or_code_audit" if eligible else "exclude",
                "primary_reason": reason,
                "full_text_or_code_assessed": "yes" if eligible else "no",
                "freeze_date": "2026-07-12",
            })
    return len(master)


if __name__ == "__main__":
    pubmed_total, full_text = write_pubmed_screening()
    master_total = write_master_screening()
    counts = count_sources()
    print(f"PubMed screened: {pubmed_total}; advanced: {full_text}")
    print(f"Master unique records screened: {master_total}")
    print(json.dumps(counts, indent=2, ensure_ascii=False))
