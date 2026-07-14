#!/usr/bin/env python3
"""Prepare the stage-confounded GSE299988 challenge matrix reproducibly."""

import csv
import gzip
import json
import statistics
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "phase1_cross_definition"
XLSX = BASE / "raw" / "GSE299988_Processed_data_files.xlsx"
GPL = BASE / "raw" / "GPL21185_family.soft.gz"
OUT_EXPR = BASE / "processed" / "GSE299988_gene_expression.tsv.gz"
OUT_META = BASE / "processed" / "GSE299988_samples.tsv"
OUT_QC = BASE / "qc" / "GSE299988_preparation_qc.json"

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def read_platform():
    mapping = {}
    in_table = False
    header = None
    with gzip.open(GPL, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line == "!platform_table_begin":
                in_table = True
                continue
            if line == "!platform_table_end":
                break
            if not in_table:
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
            probe, symbol, control = fields[idx_id], fields[idx_symbol].strip(), fields[idx_control]
            if control.upper() == "FALSE" and symbol and "///" not in symbol:
                mapping[probe] = symbol
    return mapping


def read_xlsx():
    with ZipFile(XLSX) as zf:
        ss_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        shared = ["".join(t.text or "" for t in si.findall(".//m:t", NS))
                  for si in ss_root.findall("m:si", NS)]
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


def iqr(values):
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[2] - q[0]


def main():
    mapping = read_platform()
    header, rows = read_xlsx()
    samples = header[1:]
    by_gene = defaultdict(list)
    mapped_rows = 0
    for row in rows:
        if len(row) != len(header):
            continue
        symbol = mapping.get(row[0])
        if not symbol:
            continue
        values = [float(x) for x in row[1:]]
        by_gene[symbol].append((iqr(values), row[0], values))
        mapped_rows += 1

    selected = {}
    for symbol, candidates in by_gene.items():
        selected[symbol] = max(candidates, key=lambda x: (x[0], x[1]))

    OUT_EXPR.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT_EXPR, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene"] + samples)
        for symbol in sorted(selected):
            writer.writerow([symbol] + selected[symbol][2])

    meta_rows = []
    for sample in samples:
        if sample.startswith("APTC"):
            group, rai, ln = "RAI_avid_LN_negative", "avid", "negative"
        elif sample.startswith("NaPTC"):
            group, rai, ln = "RAI_nonavid_LN_positive", "nonavid", "positive"
        else:
            group, rai, ln = "adjacent_normal", "not_applicable", "not_applicable"
        meta_rows.append([sample, group, rai, ln])
    with open(OUT_META, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "analysis_group", "rai_status", "ln_status"])
        writer.writerows(meta_rows)

    qc = {
        "xlsx_dimension": "A1:N58342",
        "matrix_samples_observed": len(samples),
        "series_samples_reported": 14,
        "normal_samples_observed": sum(x.startswith("Normal") for x in samples),
        "normal_samples_reported": 4,
        "missing_expected_column": "Normal2",
        "platform_probe_mappings": len(mapping),
        "expression_probe_rows": len(rows),
        "mapped_expression_rows": mapped_rows,
        "collapsed_gene_symbols": len(selected),
        "collapse_rule": "largest IQR across all 13 matrix samples; ties by probe ID",
        "primary_challenge": "5 RAI-nonavid/LN-positive vs 5 RAI-avid/LN-negative tumors",
        "identifiability": "RAI status and LN status are perfectly co-linear",
    }
    OUT_QC.parent.mkdir(parents=True, exist_ok=True)
    OUT_QC.write_text(json.dumps(qc, indent=2) + "\n")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
