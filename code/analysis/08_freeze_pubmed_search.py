#!/usr/bin/env python3
"""Freeze the prespecified PubMed search and export auditable records."""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


SEARCH_DATE = "2026-07-12"
QUERY = r'''((thyroid[Title/Abstract] AND (cancer*[Title/Abstract] OR carcinoma*[Title/Abstract] OR tumor*[Title/Abstract] OR tumour*[Title/Abstract] OR neoplasm*[Title/Abstract])) AND (radioiodine[Title/Abstract] OR radioiodide[Title/Abstract] OR "radioactive iodine"[Title/Abstract] OR "iodine-131"[Title/Abstract] OR 131I[Title/Abstract] OR "iodine avidity"[Title/Abstract] OR "iodine uptake"[Title/Abstract] OR "iodine content"[Title/Abstract] OR RAI[Title/Abstract]) AND (refractor*[Title/Abstract] OR resist*[Title/Abstract] OR nonavid[Title/Abstract] OR non-avid[Title/Abstract] OR avid*[Title/Abstract] OR uptake[Title/Abstract] OR response[Title/Abstract] OR ablation[Title/Abstract] OR sensitive[Title/Abstract] OR progression[Title/Abstract] OR recurrence[Title/Abstract]) AND (transcriptom*[Title/Abstract] OR "gene expression"[Title/Abstract] OR RNA-seq[Title/Abstract] OR "RNA sequencing"[Title/Abstract] OR microarray[Title/Abstract] OR miRNA[Title/Abstract] OR microRNA[Title/Abstract] OR signature[Title/Abstract] OR "gene set"[Title/Abstract] OR module[Title/Abstract] OR profile[Title/Abstract] OR score[Title/Abstract] OR biomarker*[Title/Abstract])) AND ("1900/01/01"[Date - Publication] : "2026/07/12"[Date - Publication])'''

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "RAIR-transcriptomic-audit/1.0 (systematic-search-freeze)"}
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "literature_search"


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=90) as response:
        payload = response.read()
    time.sleep(0.4)
    return payload


def clean_text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return " ".join("".join(node.itertext()).split())


def article_year(article: ET.Element) -> str:
    for path in (
        ".//JournalIssue/PubDate/Year",
        ".//ArticleDate/Year",
        ".//PubMedPubDate[@PubStatus='pubmed']/Year",
        ".//PubMedPubDate[@PubStatus='entrez']/Year",
    ):
        value = article.findtext(path, default="").strip()
        if value:
            return value
    medline = article.findtext(".//JournalIssue/PubDate/MedlineDate", default="")
    match = re.search(r"(?:19|20)\d{2}", medline)
    return match.group(0) if match else ""


def parse_article(node: ET.Element) -> dict[str, str]:
    article = node.find("MedlineCitation/Article")
    citation = node.find("MedlineCitation")
    pmid = citation.findtext("PMID", default="") if citation is not None else ""
    title = clean_text(article.find("ArticleTitle")) if article is not None else ""
    abstract = " ".join(
        clean_text(x) for x in article.findall("Abstract/AbstractText")
    ) if article is not None else ""
    journal = clean_text(article.find("Journal/Title")) if article is not None else ""
    authors = []
    if article is not None:
        for author in article.findall("AuthorList/Author"):
            collective = author.findtext("CollectiveName", default="").strip()
            if collective:
                authors.append(collective)
                continue
            last = author.findtext("LastName", default="").strip()
            fore = author.findtext("ForeName", default="").strip()
            name = ", ".join(x for x in (last, fore) if x)
            if name:
                authors.append(name)
    doi = ""
    for identifier in node.findall("PubmedData/ArticleIdList/ArticleId"):
        if identifier.attrib.get("IdType") == "doi":
            doi = clean_text(identifier)
            break
    publication_types = "; ".join(
        clean_text(x) for x in article.findall("PublicationTypeList/PublicationType")
    ) if article is not None else ""
    language = "; ".join(
        clean_text(x) for x in article.findall("Language")
    ) if article is not None else ""
    haystack = f"{title} {abstract}".lower()
    screen_terms = {
        "transcriptomic_or_expression": bool(re.search(r"transcriptom|rna[- ]?seq|gene expression|microarray", haystack)),
        "rai_endpoint": bool(re.search(r"radioiod|radioactive iodine|iodine avid|iodine uptake|iodine content|rai[- ]?refract", haystack)),
        "tool_language": bool(re.search(r"signature|score|model|module|profile|gene set|classifier|predict", haystack)),
        "human_tissue_or_clinical": bool(re.search(r"patient|tumou?r tissue|cohort|clinical|metasta", haystack)),
    }
    return {
        "record_type": "journal_article",
        "pmid": pmid,
        "doi": doi,
        "year": article_year(node),
        "title": title,
        "journal": journal,
        "authors": "; ".join(authors),
        "abstract": abstract,
        "publication_types": publication_types,
        "language": language,
        **{key: "yes" if value else "no" for key, value in screen_terms.items()},
        "screening_decision": "pending_manual_screen",
        "screening_reason": "",
    }


def parse_book_article(node: ET.Element) -> dict[str, str]:
    document = node.find("BookDocument")
    pmid = document.findtext("PMID", default="") if document is not None else ""
    title = clean_text(document.find("ArticleTitle")) if document is not None else ""
    abstract = " ".join(
        clean_text(x) for x in document.findall("Abstract/AbstractText")
    ) if document is not None else ""
    journal = clean_text(document.find("Book/BookTitle")) if document is not None else ""
    year = document.findtext("Book/PubDate/Year", default="") if document is not None else ""
    authors = []
    if document is not None:
        for author in document.findall("AuthorList/Author"):
            last = author.findtext("LastName", default="").strip()
            fore = author.findtext("ForeName", default="").strip()
            name = ", ".join(x for x in (last, fore) if x)
            if name:
                authors.append(name)
    publication_types = "; ".join(
        clean_text(x) for x in document.findall("PublicationType")
    ) if document is not None else ""
    language = "; ".join(
        clean_text(x) for x in document.findall("Language")
    ) if document is not None else ""
    haystack = f"{title} {abstract}".lower()
    return {
        "record_type": "book_or_chapter",
        "pmid": pmid,
        "doi": "",
        "year": year,
        "title": title,
        "journal": journal,
        "authors": "; ".join(authors),
        "abstract": abstract,
        "publication_types": publication_types,
        "language": language,
        "transcriptomic_or_expression": "yes" if re.search(r"transcriptom|rna[- ]?seq|gene expression|microarray", haystack) else "no",
        "rai_endpoint": "yes" if re.search(r"radioiod|radioactive iodine|iodine avid|iodine uptake|iodine content|rai[- ]?refract", haystack) else "no",
        "tool_language": "yes" if re.search(r"signature|score|model|module|profile|gene set|classifier|predict", haystack) else "no",
        "human_tissue_or_clinical": "yes" if re.search(r"patient|tumou?r tissue|cohort|clinical|metasta", haystack) else "no",
        "screening_decision": "exclude",
        "screening_reason": "not_original_research_book_or_chapter",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    params = urllib.parse.urlencode({
        "db": "pubmed",
        "term": QUERY,
        "retmax": 1000,
        "retmode": "json",
        "sort": "pub date",
        "tool": "RAIR-transcriptomic-audit",
    })
    search_payload = json.loads(get(f"{BASE}/esearch.fcgi?{params}"))
    result = search_payload["esearchresult"]
    ids = result["idlist"]
    records = []
    for start in range(0, len(ids), 150):
        batch = ids[start:start + 150]
        fetch_params = urllib.parse.urlencode({
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "xml",
            "retmode": "xml",
            "tool": "RAIR-transcriptomic-audit",
        })
        root = ET.fromstring(get(f"{BASE}/efetch.fcgi?{fetch_params}"))
        records.extend(parse_article(node) for node in root.findall("PubmedArticle"))
        records.extend(parse_book_article(node) for node in root.findall("PubmedBookArticle"))
    order = {pmid: index for index, pmid in enumerate(ids)}
    records.sort(key=lambda record: order.get(record["pmid"], len(order)))
    if len(records) != len(ids):
        raise RuntimeError(f"Fetched {len(records)} articles for {len(ids)} PMIDs")

    metadata = {
        "database": "PubMed",
        "search_date": SEARCH_DATE,
        "coverage": "database inception through 2026-07-12",
        "query_version": "v2_high_sensitivity_after_supplementary-source_recall_check",
        "language_filter": "none at retrieval; English assessed during screening",
        "article_type_filter": "none at retrieval",
        "query": QUERY,
        "record_count": int(result["count"]),
        "returned_count": len(ids),
        "pmids": ids,
        "retrieval_tool": "NCBI E-utilities esearch+efetch",
    }
    (OUT / "pubmed_search_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (OUT / "pubmed_query.txt").write_text(QUERY + "\n", encoding="utf-8")
    fields = list(records[0])
    with (OUT / "pubmed_search_records.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps({"count": len(records), "out": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
