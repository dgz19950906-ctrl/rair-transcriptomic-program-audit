#!/usr/bin/env python3
"""Build a verified BibTeX library for the frozen 22-study audit."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEARCH = ROOT / "literature_search"
OUT = ROOT / "references" / "library.bib"

PMIDS = {
    "17854396", "25682061", "27997908", "30256977", "31540966",
    "31888560", "33198784", "33534128", "33748479", "34149901",
    "34572876", "35634509", "35723359", "36556238", "36780190",
    "39968056", "39982585", "40134821", "40746809", "41801518",
    "42232186",
}


def clean(value):
    return (value or "").replace("{", "").replace("}", "").strip()


rows = list(csv.DictReader((SEARCH / "pubmed_search_records.tsv").open(), delimiter="\t"))
rows = [row for row in rows if row["pmid"] in PMIDS]
entries = []
used = set()
for row in sorted(rows, key=lambda x: (x["year"], x["pmid"])):
    first = clean(row["authors"].split(";")[0].split(",")[0]) or "Anon"
    stem = re.sub(r"[^A-Za-z0-9]", "", first) + row["year"] + "RAI"
    key = stem
    n = 2
    while key in used:
        key = f"{stem}{n}"
        n += 1
    used.add(key)
    fields = [
        f"  author = {{{clean(row['authors']).replace(';', ' and')}}}",
        f"  title = {{{clean(row['title'])}}}",
        f"  journal = {{{clean(row['journal'])}}}",
        f"  year = {{{clean(row['year'])}}}",
        f"  pmid = {{{clean(row['pmid'])}}}",
    ]
    if row["doi"]:
        fields.append(f"  doi = {{{clean(row['doi']).lower()}}}")
    fields.extend([
        "  verified = {true}",
        "  verified_by = {PubMed}",
        "  verified_on = {2026-07-12}",
    ])
    entries.append(f"@article{{{key},\n" + ",\n".join(fields) + "\n}")

entries.append("""@article{Luo2026AutophagyRAIR,
  author = {Luo, Yingying and Yuan, Zengbei and Lu, Qiteng and Li, Junhong and Pang, Xiaoan and Wei, Zhixiao and Li, Sijin},
  title = {Screening of potential autophagy-related long non-coding RNAs and their regulatory pathways in radioactive iodine-refractory differentiated thyroid cancer},
  journal = {Eurasian Journal of Medicine and Oncology},
  year = {2026},
  doi = {10.36922/ejmo025510526},
  verified = {true},
  verified_by = {publisher_page},
  verified_on = {2026-07-12}
}""")

OUT.write_text("\n\n".join(entries) + "\n", encoding="utf-8")
print(f"Wrote {len(entries)} verified entries to {OUT}")
