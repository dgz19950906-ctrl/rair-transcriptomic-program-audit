#!/usr/bin/env python3
"""Freeze first-appearance Vancouver numbering for the RAIR manuscript.

This does not hand-format the bibliography. It provides a deterministic
citekey-to-number map for a later CSL render and verifies that every mapped key
has a fully verified reference record.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path("/home/dony/ThyroidCancer_Project/rair_audit")
MANUSCRIPT = ROOT / "manuscript/two_null_v2/MANUSCRIPT_FULL_DRAFT.two_null_v2.md"
AUDIT = ROOT / "pre_submission_audit_v1/reference_gate/qc/reference_audit.json"
OUT = ROOT / "pre_submission_audit_v1/reference_gate/vancouver_map_v1"
CITE_RE = re.compile(r"(?<![A-Za-z0-9_])-?@([A-Za-z][\w:.\-/+]*)")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    text = MANUSCRIPT.read_text(encoding="utf-8")
    text_no_code = re.sub(r"```.*?```", "", text, flags=re.S)
    text_no_code = re.sub(r"`[^`\n]+`", "", text_no_code)
    ordered: list[str] = []
    for key in CITE_RE.findall(text_no_code):
        if key not in ordered:
            ordered.append(key)

    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    by_key = {r["ref_id"]: r for r in audit["records"]}
    rows = []
    for n, key in enumerate(ordered, 1):
        r = by_key.get(key)
        rows.append({
            "vancouver_number": n,
            "citekey": key,
            "verification_status": r["status"] if r else "MISSING",
            "pmid": r["pmid"] if r else "",
            "doi": r["doi"] if r else "",
            "title": r["title_guess"] if r else "",
        })

    with (OUT / "VANCOUVER_NUMBER_MAP.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    safe = (
        len(rows) == 27
        and all(r["verification_status"] == "OK" for r in rows)
        and audit.get("submission_safe") is True
        and audit.get("fully_verified") is True
        and not audit.get("duplicate_findings")
    )
    report = [
        "# Vancouver citation-order gate",
        "",
        f"- Unique in-text citekeys: {len(rows)}",
        f"- Fully verified mappings: {sum(r['verification_status'] == 'OK' for r in rows)}",
        f"- Undefined mappings: {sum(r['verification_status'] == 'MISSING' for r in rows)}",
        f"- Duplicate PMID/DOI findings: {len(audit.get('duplicate_findings', []))}",
        f"- Number-map safe: {'yes' if safe else 'no'}",
        "",
        "The map follows order of first appearance and is suitable for a Vancouver numeric CSL render. It is not a hand-formatted reference list.",
        "The server Pandoc version is 2.9.2.1 and lacks both the built-in --citeproc option and the legacy pandoc-citeproc executable; final CSL rendering therefore remains a submission-package build step unless a portable newer Pandoc is authorized.",
        "",
        f"Manuscript SHA-256: `{sha256(MANUSCRIPT)}`",
        f"Reference audit SHA-256: `{sha256(AUDIT)}`",
    ]
    (OUT / "VANCOUVER_STYLE_GATE.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (OUT / "VANCOUVER_MAP_AUDIT.json").write_text(json.dumps({
        "manuscript_sha256": sha256(MANUSCRIPT),
        "reference_audit_sha256": sha256(AUDIT),
        "n_citekeys": len(rows),
        "all_verified": all(r["verification_status"] == "OK" for r in rows),
        "submission_safe_number_map": safe,
        "final_csl_render_complete": False,
        "render_blocker": "server_pandoc_2.9.2.1_without_citeproc",
    }, indent=2), encoding="utf-8")
    print(json.dumps({"n": len(rows), "safe": safe, "render_complete": False}, indent=2))


if __name__ == "__main__":
    main()
