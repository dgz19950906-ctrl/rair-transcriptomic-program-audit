#!/usr/bin/env python3
"""Freeze supplementary discovery searches for the RAIR transcriptomic audit."""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


SEARCH_DATE = "2026-07-12"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "literature_search"
UA = {"User-Agent": "RAIR-transcriptomic-audit/1.0 (systematic-search-freeze)"}

# This is the exact relevance query that produced the archived 100-record export.
# The endpoint reported 242 matches; only the first 100 were retained before API
# rate limiting, so this source is explicitly a supplementary sensitivity search.
SEMANTIC_QUERY = "radioiodine refractory thyroid cancer transcriptomic gene expression signature"
GEO_QUERY = '(thyroid[All Fields]) AND (radioiodine[All Fields] OR "radioactive iodine"[All Fields] OR RAIR[All Fields] OR "iodine avidity"[All Fields]) AND "Homo sapiens"[Organism]'
EUROPE_PMC_QUERY = '((TITLE_ABS:"radioiodine" OR TITLE_ABS:"radioactive iodine" OR TITLE_ABS:"iodine avidity" OR TITLE_ABS:"RAIR") AND (TITLE_ABS:"thyroid cancer" OR TITLE_ABS:"thyroid carcinoma") AND (TITLE_ABS:transcriptom* OR TITLE_ABS:"gene expression" OR TITLE_ABS:signature OR TITLE_ABS:"RNA-seq" OR TITLE_ABS:microarray)) AND FIRST_PDATE:[1900-01-01 TO 2026-07-12]'


def get_json(url: str) -> dict:
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=90) as response:
        payload = json.loads(response.read().decode("utf-8"))
    time.sleep(0.5)
    return payload


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def semantic_scholar() -> dict:
    query_params = {
        "query": SEMANTIC_QUERY,
        "limit": 100,
        "fields": "title,year,venue,authors,abstract,externalIds,url,publicationTypes",
    }
    params = urllib.parse.urlencode(query_params)
    payload = get_json(f"https://api.semanticscholar.org/graph/v1/paper/search?{params}")
    reported_total = payload.get("total")
    rows = []
    for item in payload.get("data", []):
        external = item.get("externalIds") or {}
        rows.append({
            "source": "Semantic Scholar",
            "paper_id": item.get("paperId", ""),
            "pmid": external.get("PubMed", ""),
            "doi": external.get("DOI", ""),
            "year": item.get("year", ""),
            "title": item.get("title", ""),
            "venue": item.get("venue", ""),
            "authors": "; ".join(x.get("name", "") for x in item.get("authors", [])),
            "abstract": item.get("abstract") or "",
            "url": item.get("url", ""),
            "screening_decision": "pending_manual_screen",
            "screening_reason": "",
        })
    fields = list(rows[0]) if rows else ["source", "paper_id", "pmid", "doi", "year", "title", "venue", "authors", "abstract", "url", "screening_decision", "screening_reason"]
    write_tsv(OUT / "semantic_scholar_search_records.tsv", rows, fields)
    return {
        "database": "Semantic Scholar",
        "query": SEMANTIC_QUERY,
        "reported_total": reported_total,
        "retrieved": len(rows),
        "limitation": "Top 100 relevance-ranked results only; API rate limiting prevented exhaustive retrieval.",
    }


def geo_datasets() -> dict:
    params = urllib.parse.urlencode({"db": "gds", "term": GEO_QUERY, "retmax": 500, "retmode": "json"})
    search = get_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}")
    result = search["esearchresult"]
    ids = result.get("idlist", [])
    rows = []
    if ids:
        summary_params = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
        summary = get_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{summary_params}")
        for uid in summary.get("result", {}).get("uids", []):
            item = summary["result"].get(uid, {})
            rows.append({
                "source": "NCBI GEO DataSets",
                "uid": uid,
                "accession": item.get("accession", ""),
                "gds_type": item.get("gdstype", ""),
                "organism": item.get("taxon", ""),
                "sample_count": item.get("n_samples", ""),
                "title": item.get("title", ""),
                "summary": item.get("summary", ""),
                "pubmed_ids": ";".join(str(x) for x in item.get("pubmedids", [])),
                "screening_decision": "pending_manual_screen",
                "screening_reason": "",
            })
    fields = list(rows[0]) if rows else ["source", "uid", "accession", "gds_type", "organism", "sample_count", "title", "summary", "pubmed_ids", "screening_decision", "screening_reason"]
    write_tsv(OUT / "geo_search_records.tsv", rows, fields)
    return {"database": "NCBI GEO DataSets", "query": GEO_QUERY, "reported_total": int(result.get("count", 0)), "retrieved": len(rows)}


def europe_pmc() -> dict:
    cursor = "*"
    rows = []
    hit_count = None
    while True:
        params = urllib.parse.urlencode({
            "query": EUROPE_PMC_QUERY,
            "format": "json",
            "pageSize": 1000,
            "cursorMark": cursor,
            "resultType": "core",
        })
        payload = get_json(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{params}")
        hit_count = payload.get("hitCount", hit_count)
        items = payload.get("resultList", {}).get("result", [])
        for item in items:
            rows.append({
                "source": "Europe PMC",
                "source_id": item.get("id", ""),
                "pmid": item.get("pmid", ""),
                "pmcid": item.get("pmcid", ""),
                "doi": item.get("doi", ""),
                "year": item.get("pubYear", ""),
                "title": item.get("title", ""),
                "journal": item.get("journalTitle", ""),
                "authors": item.get("authorString", ""),
                "abstract": item.get("abstractText", ""),
                "publication_types": "; ".join(item.get("pubTypeList", {}).get("pubType", [])),
                "screening_decision": "pending_manual_screen",
                "screening_reason": "",
            })
        next_cursor = payload.get("nextCursorMark")
        if not items or not next_cursor or next_cursor == cursor or len(rows) >= int(hit_count or 0):
            break
        cursor = next_cursor
    fields = list(rows[0]) if rows else ["source", "source_id", "pmid", "pmcid", "doi", "year", "title", "journal", "authors", "abstract", "publication_types", "screening_decision", "screening_reason"]
    write_tsv(OUT / "europe_pmc_search_records.tsv", rows, fields)
    return {"database": "Europe PMC", "query": EUROPE_PMC_QUERY, "reported_total": hit_count, "retrieved": len(rows)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reports = []
    errors = []
    for name, function in (("semantic_scholar", semantic_scholar), ("geo", geo_datasets), ("europe_pmc", europe_pmc)):
        try:
            reports.append(function())
        except Exception as exc:
            errors.append({"source": name, "error_type": type(exc).__name__, "message": str(exc)})
    metadata = {"search_date": SEARCH_DATE, "reports": reports, "errors": errors}
    (OUT / "supplementary_search_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False))


if __name__ == "__main__":
    main()
