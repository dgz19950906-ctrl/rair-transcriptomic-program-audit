#!/usr/bin/env python3
"""Extract section-aware plain text from downloaded PMC XML for audit searches."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "literature_search" / "fulltext_xml"
OUT = ROOT / "literature_search" / "fulltext_text"


def text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path in sorted(SOURCE.glob("PMC*.xml")):
        root = ET.parse(path).getroot()
        article = root.find(".//article")
        if article is None:
            continue
        lines = []
        title = text(article.find(".//article-title"))
        if title:
            lines.extend((f"# {title}", ""))
        abstract = article.find(".//abstract")
        if abstract is not None:
            lines.extend(("## ABSTRACT", text(abstract), ""))
        for section in article.findall(".//body//sec"):
            heading = text(section.find("title")) or "UNTITLED SECTION"
            lines.append(f"## {heading}")
            for paragraph in section.findall("./p"):
                value = text(paragraph)
                if value:
                    lines.append(value)
            for table in section.findall("./table-wrap"):
                label = text(table.find("label"))
                caption = text(table.find("caption"))
                body = text(table.find("table"))
                lines.append(f"[TABLE {label}] {caption} {body}".strip())
            lines.append("")
        for section in article.findall(".//back/sec"):
            heading = text(section.find("title")) or "BACK MATTER"
            lines.append(f"## {heading}")
            for paragraph in section.findall("./p"):
                value = text(paragraph)
                if value:
                    lines.append(value)
            lines.append("")
        (OUT / f"{path.stem}.txt").write_text("\n".join(lines), encoding="utf-8")
    print(f"Extracted {len(list(OUT.glob('PMC*.txt')))} full texts")


if __name__ == "__main__":
    main()
