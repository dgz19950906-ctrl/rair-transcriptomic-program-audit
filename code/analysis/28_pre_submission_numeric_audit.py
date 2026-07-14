#!/usr/bin/env python3
"""Pre-submission numerical reconciliation for RAIR audit Figures 3--5.

This script is intentionally stdlib-only and is designed to run on Xiyou Cloud.
It compares frozen analysis tables, figure source tables, manuscript claims, and
figure provenance. It does not recompute biological or statistical results.
"""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path("/home/dony/ThyroidCancer_Project/rair_audit")
OUT = ROOT / "pre_submission_audit_v1" / "numeric_gate_v3"
MANUSCRIPT = ROOT / "manuscript/two_null_v2/MANUSCRIPT_FULL_DRAFT.two_null_v2.md"

FIG3_CLIN = ROOT / "covariance_null_v2/clinical_label_layer_v1/tables/clinical_endpoint_two_null_results.tsv"
FIG3_ALL = ROOT / "covariance_null_v2/figure3_two_null_v4/tables/Figure3_all18_results.tsv"
FIG3_B = ROOT / "covariance_null_v2/figure3_two_null_v4/tables/Figure3_panel_b_two_null_bh_source.tsv"
FIG3_PDF = ROOT / "figures/final/Figure3_lopo_two_null_calibration.pdf"
FIG3_CANON_PDF = ROOT / "covariance_null_v2/figure3_two_null_v4/figures/Figure3_lopo_two_null_calibration.pdf"
FIG3_SVG = ROOT / "figures/final/Figure3_lopo_two_null_calibration.svg"

FIG2_SOURCE = ROOT / "results/10_figure2_two_axis_source.tsv"

FIG4_CLIN = ROOT / "covariance_null_gse299988_v1/clinical_label_layer_v1/tables/GSE299988_two_null_challenge.tsv"
FIG4_B = ROOT / "figures/figure4_two_null_v3/Figure4_panel_b_effects_source.tsv"
FIG4_C = ROOT / "figures/figure4_two_null_v3/Figure4_panel_c_two_null_q_source.tsv"
FIG4_D = ROOT / "figures/figure4_two_null_v3/Figure4_panel_d_shift_source.tsv"
FIG4_SVG = ROOT / "figures/figure4_two_null_v3/Figure4_GSE299988_two_null_challenge_v3.svg"

FIG5_BOOT = ROOT / "results/08_three_atlas_bootstrap_effects.tsv"
FIG5_GATE = ROOT / "results/07_three_atlas_AND_gate.tsv"
FIG5_LODO = ROOT / "results/04_lodo_summary.tsv"
FIG5_SCRIPT = ROOT / "scripts/build_figure5.R"
FIG5_PDF = ROOT / "figures/final/Figure5_cellular_compartment_validation.pdf"
FIG5_FREEZE = ROOT / "manuscript/FIGURE5_FREEZE.md"


rows: list[dict[str, str]] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def f(value: str | float) -> float:
    return float(value)


def equal(a: str | float, b: str | float, tol: float = 1e-12) -> bool:
    try:
        aa, bb = f(a), f(b)
        if math.isnan(aa) and math.isnan(bb):
            return True
        return abs(aa - bb) <= tol
    except (TypeError, ValueError):
        return str(a) == str(b)


def add(gate: str, item: str, source: str, target: str, status: bool,
        severity: str = "P0", note: str = "") -> None:
    rows.append({
        "gate": gate,
        "item": item,
        "source": source,
        "target": target,
        "status": "PASS" if status else "FAIL",
        "severity_if_fail": severity,
        "note": note,
    })


def index(data: list[dict[str, str]], *keys: str) -> dict[tuple[str, ...], dict[str, str]]:
    return {tuple(row[k] for k in keys): row for row in data}


def require_text(label: str, expected: str, manuscript: str, note: str = "") -> None:
    add("manuscript", label, expected, "MANUSCRIPT_FULL_DRAFT", expected in manuscript,
        "P0", note)


def svg_text(path: Path) -> str:
    raw = path.read_text(encoding="utf-8", errors="replace")
    raw = re.sub(r"<[^>]+>", " ", raw)
    return re.sub(r"\s+", " ", html.unescape(raw))


def read_fig4_panel_d(path: Path, record_schema: bool = True) -> list[dict[str, str]]:
    """Read the currently frozen panel-d TSV while documenting its newline defect.

    The display labels contain literal newlines and were written without field
    quoting, splitting every intended row into two physical lines. Reconstructing
    the rows here permits value auditing but does not waive the packaging failure.
    """
    physical = path.read_text(encoding="utf-8").splitlines()
    header = physical[0].split("\t")
    normal = [line.split("\t") for line in physical[1:]]
    rectangular = all(len(parts) == len(header) for parts in normal)
    if record_schema:
        add("Fig4", "panel-d Source Data is rectangular TSV", "6 rows x 6 columns",
            f"{len(normal)} physical rows; widths={sorted(set(map(len, normal)))}",
            rectangular, "P0", "Literal display newlines split each intended data row")
    if rectangular:
        return [dict(zip(header, parts)) for parts in normal]
    rebuilt: list[dict[str, str]] = []
    for i in range(0, len(normal), 2):
        first, second = normal[i], normal[i + 1]
        if len(first) != 3 or len(second) != 4:
            raise ValueError(f"Unexpected malformed panel-d lines at {i + 2}")
        rebuilt.append({
            "signature_id": first[0], "program": first[1],
            "dataset": first[2] + "\n" + second[0],
            "estimate": second[1], "lower": second[2], "upper": second[3],
        })
    return rebuilt


def audit_figure3(manuscript: str) -> None:
    clinical = read_tsv(FIG3_CLIN)
    all18 = read_tsv(FIG3_ALL)
    panel = read_tsv(FIG3_B)
    add("Fig3", "18 program-endpoint rows", "clinical label layer", "Figure3_all18_results.tsv",
        len(clinical) == len(all18) == 18, "P0", f"n={len(clinical)}/{len(all18)}")

    ci = index(clinical, "signature_id", "contrast")
    ai = index(all18, "signature_id", "contrast")
    columns = [
        "raw_hedges_g", "adverse_aligned_hedges_g", "bootstrap_ci_low", "bootstrap_ci_high",
        "exact_label_p", "covariance_program_p", "exact_label_q_bh18", "covariance_program_q_bh18",
    ]
    for key, src in ci.items():
        dst = ai.get(key)
        add("Fig3", f"row present {key}", "clinical label layer", "Figure3 all18", dst is not None)
        if dst:
            for col in columns:
                add("Fig3", f"{key} {col}", src[col], dst[col], equal(src[col], dst[col]))

    pi = index(panel, "signature_id", "contrast", "null_type")
    for key, src in ci.items():
        sid, contrast = key
        for null_type, pcol, qcol in [
            ("Exact patient-label null", "exact_label_p", "exact_label_q_bh18"),
            ("Covariance-matched program null", "covariance_program_p", "covariance_program_q_bh18"),
        ]:
            dst = pi.get((sid, contrast, null_type))
            add("Fig3", f"panel-b row {(sid, contrast,null_type)}", "clinical label layer", "panel b", dst is not None)
            if dst:
                add("Fig3", f"panel-b p {(sid,contrast,null_type)}", src[pcol], dst["p_value"], equal(src[pcol], dst["p_value"]))
                add("Fig3", f"panel-b q {(sid,contrast,null_type)}", src[qcol], dst["q_value"], equal(src[qcol], dst["q_value"]))

    add("Fig3", "final PDF equals v4 canonical PDF", sha256(FIG3_CANON_PDF), sha256(FIG3_PDF),
        sha256(FIG3_CANON_PDF) == sha256(FIG3_PDF))

    # Every panel-b q value is printed to three decimals in the SVG. Verify multiplicity.
    observed = Counter(re.findall(r"(?<!\d)(?:0|1)\.\d{3}(?!\d)", svg_text(FIG3_SVG)))
    expected = Counter(f"{f(r['q_value']):.3f}" for r in panel)
    for value, count in sorted(expected.items()):
        add("Fig3", f"SVG printed q={value} multiplicity", str(count), str(observed[value]),
            observed[value] >= count, "P0")

    min_exact_q = min(f(r["exact_label_q_bh18"]) for r in clinical)
    primary = [r for r in clinical if r["signature_id"] in {"TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX"}]
    require_text("Fig3 minimum exact q", f"minimum q = {min_exact_q:.3f}", manuscript)
    exact_range = (min(f(r["exact_label_q_bh18"]) for r in primary), max(f(r["exact_label_q_bh18"]) for r in primary))
    cov_range = (min(f(r["covariance_program_q_bh18"]) for r in primary), max(f(r["covariance_program_q_bh18"]) for r in primary))
    require_text("Fig3 primary exact q range", f"exact-label q values ranged from {exact_range[0]:.3f} to {exact_range[1]:.3f}", manuscript)
    require_text("Fig3 primary covariance q range", f"program-identity q values from {cov_range[0]:.3f} to {cov_range[1]:.3f}", manuscript)

    lookups = index(clinical, "program_label", "contrast")
    for program, contrast, label in [
        ("Hypoxia", "uptake_failure", "Hypoxia uptake"),
        ("Hypoxia", "response_failure_with_uptake", "Hypoxia response"),
        ("Inflammatory response", "response_failure_with_uptake", "Inflammatory response response"),
    ]:
        r = lookups[(program, contrast)]
        token = f"Pcov = {f(r['covariance_program_p']):.3f}, q = {f(r['covariance_program_q_bh18']):.3f}"
        require_text(label, token, manuscript)


def audit_cross_figure_ci() -> None:
    clinical = index(read_tsv(FIG3_CLIN), "signature_id", "contrast")
    fig2 = index(read_tsv(FIG2_SOURCE), "signature_id")
    fig4d = read_fig4_panel_d(FIG4_D)
    for r in fig4d:
        if r["dataset"].replace("\n", " ").startswith("GSE151179"):
            sid = r["signature_id"]
            c = clinical[(sid, "uptake_failure")]
            f2 = fig2[(sid,)]
            add("cross-figure", f"{sid} point estimate Fig2/Fig3/Fig4",
                c["adverse_aligned_hedges_g"], r["estimate"], equal(c["adverse_aligned_hedges_g"], r["estimate"]))
            add("cross-figure", f"{sid} Fig4 CI matches Figure2 source",
                f"{f2['ci_low_uptake_failure']}..{f2['ci_high_uptake_failure']}",
                f"{r['lower']}..{r['upper']}",
                equal(f2["ci_low_uptake_failure"], r["lower"]) and equal(f2["ci_high_uptake_failure"], r["upper"]))
            add("cross-figure", f"{sid} Fig4 CI matches canonical Figure3 clinical layer",
                f"{c['bootstrap_ci_low']}..{c['bootstrap_ci_high']}",
                f"{r['lower']}..{r['upper']}",
                equal(c["bootstrap_ci_low"], r["lower"]) and equal(c["bootstrap_ci_high"], r["upper"]),
                "P0", "Same estimand has different frozen bootstrap intervals across figures")


def audit_figure4(manuscript: str) -> None:
    clinical = index(read_tsv(FIG4_CLIN), "signature_id")
    panel_b = index(read_tsv(FIG4_B), "signature_id")
    panel_c = read_tsv(FIG4_C)
    panel_d = read_fig4_panel_d(FIG4_D, record_schema=False)
    add("Fig4", "nine challenge rows", "clinical label layer", "panel b", len(clinical) == len(panel_b) == 9)
    for key, src in clinical.items():
        dst = panel_b.get(key)
        add("Fig4", f"panel-b row {key[0]}", "clinical label layer", "panel b", dst is not None)
        if dst:
            for scol, dcol in [
                ("adverse_aligned_hedges_g", "adverse_aligned_hedges_g"),
                ("frozen_bootstrap_ci_low", "adverse_aligned_bootstrap_g_ci_low"),
                ("frozen_bootstrap_ci_high", "adverse_aligned_bootstrap_g_ci_high"),
                ("exact_label_p", "exact_label_p"),
                ("exact_label_q_bh9", "exact_label_q_bh9"),
                ("covariance_program_p", "covariance_program_p"),
                ("covariance_program_q_bh9", "covariance_program_q_bh9"),
            ]:
                add("Fig4", f"{key[0]} {scol}", src[scol], dst[dcol], equal(src[scol], dst[dcol]))

    # Panel c contains two rows per program, one per null.
    for r in panel_c:
        src = clinical[(r["signature_id"],)]
        if "Exact" in r["null_model"]:
            expected = src["exact_label_q_bh9"]
        else:
            expected = src["covariance_program_q_bh9"]
        add("Fig4", f"panel-c q {r['signature_id']} {r['null_model']}", expected, r["q"], equal(expected, r["q"]))

    # GSE299988 rows in panel d must match the challenge layer.
    for r in panel_d:
        if r["dataset"].replace("\n", " ").startswith("GSE299988"):
            src = clinical[(r["signature_id"],)]
            ok = (equal(src["adverse_aligned_hedges_g"], r["estimate"]) and
                  equal(src["frozen_bootstrap_ci_low"], r["lower"]) and
                  equal(src["frozen_bootstrap_ci_high"], r["upper"]))
            add("Fig4", f"panel-d challenge {r['signature_id']}",
                f"{src['adverse_aligned_hedges_g']} [{src['frozen_bootstrap_ci_low']},{src['frozen_bootstrap_ci_high']}]",
                f"{r['estimate']} [{r['lower']},{r['upper']}]", ok)

    observed = Counter(re.findall(r"(?<!\d)(?:0|1)\.\d{3}(?!\d)", svg_text(FIG4_SVG)))
    expected = Counter(f"{f(r['q']):.3f}" for r in panel_c)
    for value, count in sorted(expected.items()):
        add("Fig4", f"SVG printed q={value} multiplicity", str(count), str(observed[value]),
            observed[value] >= count)

    for sid, label in [("TDS_16", "TDS-16"), ("IODIDE_HANDLING_11", "iodide handling-11"), ("CONDELLO_2025_SIX", "Condello-6")]:
        r = clinical[(sid,)]
        token = f"{f(r['adverse_aligned_hedges_g']):.2f} ({f(r['frozen_bootstrap_ci_low']):.2f} to {f(r['frozen_bootstrap_ci_high']):.2f})"
        # Manuscript uses brackets and optional '95% bootstrap CI' before the interval; assert the three values independently.
        for value_name, value in [("g", r["adverse_aligned_hedges_g"]), ("CI low", r["frozen_bootstrap_ci_low"]), ("CI high", r["frozen_bootstrap_ci_high"])]:
            require_text(f"Fig4 {label} {value_name}", f"{f(value):.2f}", manuscript, token)


def audit_figure5(manuscript: str) -> None:
    boot = read_tsv(FIG5_BOOT)
    gate = index(read_tsv(FIG5_GATE), "signature_id")
    lodo = index(read_tsv(FIG5_LODO), "signature_id")
    primary_ids = ["TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX"]
    primary = [r for r in boot if r["signature_id"] in primary_ids]
    add("Fig5", "three programs x three atlases", "08 source", "expected grid", len(primary) == 9)
    for sid in primary_ids:
        g = gate[(sid,)]
        add("Fig5", f"{sid} discovery donor gate", "TRUE", g["GSE184362_donor_gate"], g["GSE184362_donor_gate"] == "TRUE")
        add("Fig5", f"{sid} three-atlas direction", "TRUE", g["three_atlas_same_direction"], g["three_atlas_same_direction"] == "TRUE")
        add("Fig5", f"{sid} AND gate", "TRUE", g["three_atlas_AND_gate"], g["three_atlas_AND_gate"] == "TRUE")
        ld = lodo[(sid,)]
        add("Fig5", f"{sid} LODO direction count", "at least 9/11",
            f"{ld['same_direction_estimable_iterations']}/{ld['direction_consistency_denominator']}",
            int(ld["same_direction_estimable_iterations"]) >= 9 and int(ld["direction_consistency_denominator"]) == 11)

    # Results must report every primary atlas estimate and both CI bounds at 2 decimals.
    for r in primary:
        for name, col in [("estimate", "estimate"), ("CI low", "conf_low"), ("CI high", "conf_high")]:
            require_text(f"Fig5 {r['program']} {r['dataset']} {name}", f"{f(r[col]):.2f}", manuscript)

    # Provenance: source table is written from boot_summary and the same object is plotted.
    script = FIG5_SCRIPT.read_text(encoding="utf-8")
    add("Fig5", "figure script writes 08 source", "write_tsv(boot_summary", "build_figure5.R",
        "write_tsv(boot_summary" in script)
    add("Fig5", "figure script plots boot_summary", "boot_summary used in forest", "build_figure5.R",
        script.count("boot_summary") >= 3, "P0",
        "PDF does not print exact values; provenance links plotted object to the frozen source table")
    freeze = FIG5_FREEZE.read_text(encoding="utf-8")
    add("Fig5", "final PDF matches freeze hash", re.search(r"Figure5_cellular_compartment_validation\.pdf.*?`([0-9a-f]{64})`", freeze).group(1),
        sha256(FIG5_PDF),
        re.search(r"Figure5_cellular_compartment_validation\.pdf.*?`([0-9a-f]{64})`", freeze).group(1) == sha256(FIG5_PDF))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    audit_figure3(manuscript)
    audit_cross_figure_ci()
    audit_figure4(manuscript)
    audit_figure5(manuscript)

    tsv_path = OUT / "NUMERIC_RECONCILIATION_FIG3_FIG5.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    failures = [r for r in rows if r["status"] == "FAIL"]
    payload = {
        "schema_version": 1,
        "manuscript_sha256": sha256(MANUSCRIPT),
        "checks": len(rows),
        "pass": len(rows) - len(failures),
        "fail": len(failures),
        "submission_safe": not failures,
        "failures": failures,
        "artifacts": {str(p): sha256(p) for p in [FIG3_ALL, FIG3_B, FIG3_PDF, FIG4_B, FIG4_C, FIG4_D, FIG4_SVG, FIG5_BOOT, FIG5_GATE, FIG5_PDF]},
    }
    (OUT / "NUMERIC_AUDIT.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "# Numerical pre-submission audit: Figures 3--5",
        "",
        f"- Checks: {len(rows)}",
        f"- Passed: {len(rows) - len(failures)}",
        f"- Failed: {len(failures)}",
        f"- Submission-safe: {'yes' if not failures else 'no'}",
        "",
        "## Failures",
        "",
    ]
    if not failures:
        lines.append("None.")
    else:
        lines += ["| Gate | Item | Source | Target | Note |", "|---|---|---:|---:|---|"]
        for r in failures:
            lines.append(f"| {r['gate']} | {r['item']} | {r['source']} | {r['target']} | {r['note']} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "A failure denotes a transcription, cross-artifact, or provenance discrepancy. It does not by itself imply that the biological conclusion changes.",
    ]
    (OUT / "NUMERIC_AUDIT_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({k: payload[k] for k in ["checks", "pass", "fail", "submission_safe"]}, indent=2))
    for row in failures:
        print("FAIL", row["gate"], row["item"], row["source"], row["target"])
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
