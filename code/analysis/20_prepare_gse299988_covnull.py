#!/usr/bin/env python3
"""Rebuild the GSE299988 expression matrix for label-blind covariance nulls.

All numerical processing is intended to run on Xiyou Cloud. The script reads
only official GEO source files already downloaded to the remote project and
writes a tumor-only gene matrix whose sample columns are GEO accessions.
"""

import argparse
import csv
import gzip
import hashlib
import json
import statistics
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile


NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_platform_and_samples(soft_path: Path):
    mapping = {}
    samples = []
    in_platform = False
    header = None
    pending_title = None
    with gzip.open(soft_path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("!Sample_title = "):
                pending_title = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_geo_accession = "):
                if pending_title is None:
                    raise RuntimeError("Sample accession encountered without title")
                samples.append((line.split(" = ", 1)[1], pending_title))
                pending_title = None
            if line == "!platform_table_begin":
                in_platform = True
                continue
            if line == "!platform_table_end":
                in_platform = False
                continue
            if not in_platform:
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                idx_id = header.index("ID")
                idx_symbol = header.index("GENE_SYMBOL")
                idx_control = header.index("CONTROL_TYPE")
                continue
            if len(fields) <= max(idx_id, idx_symbol, idx_control):
                continue
            probe = fields[idx_id]
            symbol = fields[idx_symbol].strip()
            control = fields[idx_control].upper()
            if control == "FALSE" and symbol and "///" not in symbol:
                mapping[probe] = symbol
    if header is None:
        raise RuntimeError("Platform table was not found in family SOFT")
    return mapping, samples


def read_xlsx(xlsx_path: Path):
    with ZipFile(xlsx_path) as zf:
        ss_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        shared = [
            "".join(t.text or "" for t in si.findall(".//m:t", NS))
            for si in ss_root.findall("m:si", NS)
        ]
        root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
    rows = root.findall(".//m:sheetData/m:row", NS)
    parsed = []
    for row in rows:
        values = []
        for cell in row.findall("m:c", NS):
            val_node = cell.find("m:v", NS)
            val = "" if val_node is None else val_node.text
            if cell.attrib.get("t") == "s" and val:
                val = shared[int(val)]
            values.append(val)
        parsed.append(values)
    return parsed[0], parsed[1:]


def workbook_name(accession: str, title: str) -> str:
    if "_tumor_" in title and title.endswith(tuple(f"RAIA{i}" for i in range(1, 6))):
        return "APTC" + title.rsplit("RAIA", 1)[1]
    if "_tumor_" in title and title.endswith(tuple(f"RAIR{i}" for i in range(1, 6))):
        return "NaPTC" + title.rsplit("RAIR", 1)[1]
    normal_map = {
        "GSM9051410": "Normal1",
        "GSM9051411": "Normal2",
        "GSM9051412": "Normal3",
        "GSM9051413": "Normal4",
    }
    if accession in normal_map:
        return normal_map[accession]
    raise RuntimeError(f"Unrecognized sample title: {accession} {title}")


def iqr(values):
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[2] - q[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    args = parser.parse_args()
    base = Path(args.base).resolve()
    raw = base / "raw"
    inputs = base / "inputs"
    manifests = base / "manifests"
    inputs.mkdir(parents=True, exist_ok=True)
    manifests.mkdir(parents=True, exist_ok=True)

    xlsx = raw / "GSE299988_Processed_data_files.xlsx"
    soft = raw / "GSE299988_family.soft.gz"
    out_all = inputs / "GSE299988_all13_gene_expression.tsv.gz"
    out_tumor = inputs / "GSE299988_tumor_gene_expression.tsv.gz"
    out_map = manifests / "GSE299988_sample_accession_map.tsv"
    out_qc = manifests / "GSE299988_matrix_rebuild_qc.json"
    out_hash = manifests / "GSE299988_matrix_rebuild_SHA256.tsv"
    for path in (out_all, out_tumor, out_map, out_qc, out_hash):
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite: {path}")

    platform, geo_samples = read_platform_and_samples(soft)
    header, rows = read_xlsx(xlsx)
    if header[0] != "Probesets":
        raise RuntimeError("Unexpected workbook identifier column")
    observed = header[1:]

    sample_rows = []
    by_workbook = {}
    for accession, title in geo_samples:
        wb = workbook_name(accession, title)
        role = "tumor" if "_tumor_" in title else "adjacent_normal"
        sample_rows.append((wb, accession, role, title, wb in observed))
        by_workbook[wb] = accession
    if set(observed) != {x[0] for x in sample_rows if x[4]}:
        raise RuntimeError("Workbook-to-GEO sample mapping is not bijective")
    if "Normal2" in observed or len(observed) != 13:
        raise RuntimeError("Expected 13 observed columns with Normal2 absent")

    by_gene = defaultdict(list)
    mapped_rows = 0
    malformed_rows = 0
    for row in rows:
        if len(row) != len(header):
            malformed_rows += 1
            continue
        symbol = platform.get(row[0])
        if not symbol:
            continue
        values = [float(x) for x in row[1:]]
        by_gene[symbol].append((iqr(values), row[0], values))
        mapped_rows += 1
    selected = {gene: max(candidates, key=lambda x: (x[0], x[1])) for gene, candidates in by_gene.items()}

    all_accessions = [by_workbook[x] for x in observed]
    tumor_positions = [i for i, name in enumerate(observed) if name.startswith(("APTC", "NaPTC"))]
    tumor_accessions = [all_accessions[i] for i in tumor_positions]
    if len(tumor_accessions) != 10 or not all(x.startswith("GSM") for x in tumor_accessions):
        raise RuntimeError("Tumor-only GSM accession gate failed")

    def write_matrix(path, accessions, positions):
        with gzip.open(path, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["symbol"] + accessions)
            for gene in sorted(selected):
                values = selected[gene][2]
                writer.writerow([gene] + [values[i] for i in positions])

    write_matrix(out_all, all_accessions, list(range(len(observed))))
    write_matrix(out_tumor, tumor_accessions, tumor_positions)

    with out_map.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["workbook_column", "geo_accession", "sample_role", "geo_title", "observed_in_workbook"])
        writer.writerows(sample_rows)

    qc = {
        "official_source_sha256": {xlsx.name: sha256(xlsx), soft.name: sha256(soft)},
        "expression_probe_rows": len(rows),
        "malformed_expression_rows": malformed_rows,
        "observed_sample_columns": len(observed),
        "reported_geo_samples": len(geo_samples),
        "missing_expected_workbook_column": "Normal2",
        "eligible_platform_probe_mappings": len(platform),
        "mapped_expression_rows": mapped_rows,
        "collapsed_gene_symbols": len(selected),
        "tumor_samples": len(tumor_accessions),
        "collapse_rule": "largest inclusive IQR across all 13 observed samples; ties by lexicographically largest probe ID",
        "null_generator_input": str(out_tumor),
        "null_generator_columns_are_geo_accessions_only": True,
        "clinical_endpoint_labels_in_null_generator_input": False,
    }
    expected = {
        "expression_probe_rows": 58341,
        "observed_sample_columns": 13,
        "reported_geo_samples": 14,
        "mapped_expression_rows": 48862,
        "collapsed_gene_symbols": 34729,
        "tumor_samples": 10,
    }
    deviations = {k: (qc[k], v) for k, v in expected.items() if qc[k] != v}
    qc["expected_value_gate_passed"] = not deviations
    qc["expected_value_deviations"] = deviations
    out_qc.write_text(json.dumps(qc, indent=2) + "\n")

    output_files = [out_all, out_tumor, out_map, out_qc]
    with out_hash.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "sha256"])
        for path in output_files:
            writer.writerow([str(path), sha256(path)])

    print(json.dumps(qc, indent=2))
    if deviations:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
