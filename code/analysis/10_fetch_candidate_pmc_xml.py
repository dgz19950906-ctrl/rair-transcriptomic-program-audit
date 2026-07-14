#!/usr/bin/env python3
"""Fetch open PMC XML for prespecified full-text audit candidates."""

from __future__ import annotations

import csv
import hashlib
import time
import urllib.parse
import urllib.request
from pathlib import Path


SEARCH_DATE = "2026-07-12"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "literature_search" / "fulltext_xml"
UA = {"User-Agent": "RAIR-transcriptomic-audit/1.0 (full-text-audit)"}

PMCIDS = {
    "42232186": "PMC13223037",
    "41801518": "PMC12971732",
    "40746809": "PMC12310417",
    "39982585": "PMC11845550",
    "39968056": "PMC11834997",
    "40134821": "PMC11934929",
    "36780190": "PMC10106408",
    "36556238": "PMC9788488",
    "35634509": "PMC9132198",
    "34572876": "PMC8468667",
    "34149901": "PMC8200939",
    "33748479": "PMC7970325",
    "33534128": "PMC8213564",
    "33198784": "PMC7667839",
    "31888560": "PMC6937781",
    "30256977": "PMC6435099",
    "35723359": "PMC9164071",
}


def fetch(pmcid: str) -> bytes:
    params = urllib.parse.urlencode({
        "db": "pmc",
        "id": pmcid.removeprefix("PMC"),
        "rettype": "xml",
        "retmode": "xml",
        "tool": "RAIR-transcriptomic-audit",
    })
    request = urllib.request.Request(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}",
        headers=UA,
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        payload = response.read()
    time.sleep(0.4)
    return payload


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for pmid, pmcid in PMCIDS.items():
        payload = fetch(pmcid)
        if not payload.lstrip().startswith(b"<?xml") and b"<article" not in payload[:1000]:
            raise RuntimeError(f"Unexpected PMC payload for {pmcid}")
        path = OUT / f"{pmcid}.xml"
        path.write_bytes(payload)
        rows.append({
            "pmid": pmid,
            "pmcid": pmcid,
            "retrieved_on": SEARCH_DATE,
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "local_path": str(path.relative_to(ROOT)),
        })
    with (OUT.parent / "fulltext_xml_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Fetched {len(rows)} PMC full texts")


if __name__ == "__main__":
    main()
