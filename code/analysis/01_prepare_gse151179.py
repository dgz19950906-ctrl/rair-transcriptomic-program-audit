#!/usr/bin/env python3
"""Prepare a gene-level GSE151179 matrix for the locked Phase-1 analysis.

Only complete, locally cached inputs are used. Probe-to-gene collapse is performed
without outcome labels: probes with ambiguous mappings are removed and, for each
gene, the probe with the largest IQR across the 17 eligible pre-RAI primary tumors
is retained.
"""

from __future__ import annotations

import csv
import gzip
import json
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SOFT = ROOT / "phase0_rai_audit/raw_metadata/GSE151179_samples_full.soft"
ANNOT = ROOT / "phase0_rai_audit/sample_annotation_audit.tsv"
SQLITE = ROOT / (
    "phase1_cross_definition/raw/clariomshumantranscriptcluster.db/inst/extdata/"
    "clariomshumantranscriptcluster.sqlite"
)
GENE_INFO = ROOT / "phase1_cross_definition/raw/Homo_sapiens.gene_info.gz"
OUT = ROOT / "phase1_cross_definition/processed"
QC = ROOT / "phase1_cross_definition/qc"


def parse_soft(path: Path) -> pd.DataFrame:
    sample_id = None
    in_table = False
    probe_order: list[str] | None = None
    values_by_sample: dict[str, list[float]] = {}
    current_probes: list[str] = []
    current_values: list[float] = []

    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                sample_id = line.split(" = ", 1)[1]
            elif line == "!sample_table_begin":
                if sample_id is None:
                    raise ValueError("Sample table found before sample accession")
                in_table = True
                current_probes, current_values = [], []
            elif line == "!sample_table_end":
                if sample_id is None:
                    raise ValueError("Sample table ended without sample accession")
                if probe_order is None:
                    probe_order = current_probes
                elif current_probes != probe_order:
                    raise ValueError(f"Probe order mismatch in {sample_id}")
                values_by_sample[sample_id] = current_values
                in_table = False
            elif in_table:
                if line == "ID_REF\tVALUE" or not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 2:
                    raise ValueError(f"Malformed expression row in {sample_id}: {line}")
                current_probes.append(fields[0])
                current_values.append(float(fields[1]))

    if in_table:
        raise ValueError("SOFT file ended inside a sample table")
    if probe_order is None or len(values_by_sample) != 52:
        raise ValueError(
            f"Expected 52 samples; found {len(values_by_sample)} (file may be incomplete)"
        )
    if len(probe_order) != 27189:
        raise ValueError(
            f"Expected 27,189 probes; found {len(probe_order)} (file may be incomplete)"
        )
    return pd.DataFrame(values_by_sample, index=probe_order, dtype=float)


def load_probe_map(path: Path) -> pd.DataFrame:
    with sqlite3.connect(path) as con:
        mapping = pd.read_sql_query(
            "SELECT probe_id, gene_id, is_multiple FROM probes", con
        )
    mapping["gene_id"] = mapping["gene_id"].astype("string")
    return mapping


def load_gene_symbols(path: Path) -> dict[str, str]:
    symbols: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["#tax_id"] == "9606":
                symbols[row["GeneID"]] = row["Symbol"]
    return symbols


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)

    probe_expr = parse_soft(SOFT)
    samples = pd.read_csv(ANNOT, sep="\t", dtype=str).fillna("")
    samples = samples[samples["geo_accession"].isin(probe_expr.columns)].copy()
    samples["eligible_primary_pre_rai"] = (
        samples["tissue_type"].str.casefold().eq("primary tumor")
        & samples["collection_before_after_rai"].str.casefold().eq("before")
    )
    samples["analysis_group"] = "excluded"
    eligible = samples["eligible_primary_pre_rai"]
    samples.loc[
        eligible
        & samples["patient_rai_response"].eq("Avid")
        & samples["disease"].eq("Remission"),
        "analysis_group",
    ] = "RAI_avid_remission"
    samples.loc[
        eligible
        & samples["patient_rai_response"].eq("Refractory")
        & samples["rai_uptake_at_metastatic_site"].eq("Yes")
        & samples["disease"].eq("Persistence"),
        "analysis_group",
    ] = "RAI_avid_persistent"
    samples.loc[
        eligible
        & samples["patient_rai_response"].eq("Refractory")
        & samples["rai_uptake_at_metastatic_site"].eq("No")
        & samples["disease"].eq("Persistence"),
        "analysis_group",
    ] = "RAI_nonavid_persistent"

    analysis_samples = samples.loc[
        samples["analysis_group"].ne("excluded"), "geo_accession"
    ].tolist()
    if len(analysis_samples) != 17:
        raise ValueError(f"Expected 17 eligible samples, found {len(analysis_samples)}")
    if samples.loc[samples["analysis_group"].ne("excluded"), "patient_id"].duplicated().any():
        raise ValueError("Eligible matrix contains repeated patients")

    probe_map = load_probe_map(SQLITE)
    in_platform = probe_map[probe_map["probe_id"].isin(probe_expr.index)].copy()
    counts_per_probe = in_platform.groupby("probe_id")["gene_id"].nunique(dropna=True)
    good_probe_ids = counts_per_probe[counts_per_probe.eq(1)].index
    usable = in_platform[
        in_platform["probe_id"].isin(good_probe_ids)
        & in_platform["gene_id"].notna()
        & in_platform["is_multiple"].eq(0)
    ].drop_duplicates(["probe_id", "gene_id"])

    eligible_expr = probe_expr.loc[usable["probe_id"], analysis_samples]
    q75 = eligible_expr.quantile(0.75, axis=1)
    q25 = eligible_expr.quantile(0.25, axis=1)
    iqr = (q75 - q25).rename("iqr")
    usable = usable.set_index("probe_id").join(iqr).reset_index()
    usable = usable.sort_values(
        ["gene_id", "iqr", "probe_id"], ascending=[True, False, True]
    )
    selected = usable.drop_duplicates("gene_id", keep="first").copy()

    symbols = load_gene_symbols(GENE_INFO)
    selected["symbol"] = selected["gene_id"].map(symbols)
    selected = selected[selected["symbol"].notna() & selected["symbol"].ne("-")].copy()
    # Current official symbols are expected to be unique; retain highest-IQR row if not.
    selected = selected.sort_values(
        ["symbol", "iqr", "probe_id"], ascending=[True, False, True]
    ).drop_duplicates("symbol", keep="first")

    gene_expr = probe_expr.loc[selected["probe_id"], analysis_samples].copy()
    gene_expr.index = selected["symbol"].tolist()
    gene_expr.index.name = "symbol"

    group_order = {
        "RAI_avid_remission": 0,
        "RAI_avid_persistent": 1,
        "RAI_nonavid_persistent": 2,
    }
    sample_sheet = samples[samples["geo_accession"].isin(analysis_samples)].copy()
    sample_sheet["group_order"] = sample_sheet["analysis_group"].map(group_order)
    sample_sheet = sample_sheet.sort_values(["group_order", "patient_id"])
    gene_expr = gene_expr.loc[:, sample_sheet["geo_accession"]]

    gene_expr.to_csv(OUT / "GSE151179_primary_preRAI_gene_expression.tsv.gz", sep="\t")
    sample_sheet.to_csv(OUT / "GSE151179_primary_preRAI_samples.tsv", sep="\t", index=False)
    selected[["probe_id", "gene_id", "symbol", "iqr"]].to_csv(
        OUT / "GSE151179_probe_to_gene_selected.tsv.gz", sep="\t", index=False
    )

    group_counts = sample_sheet["analysis_group"].value_counts().to_dict()
    qc = {
        "source_soft": str(SOFT.relative_to(ROOT)),
        "n_samples_in_soft": int(probe_expr.shape[1]),
        "n_probes_in_soft": int(probe_expr.shape[0]),
        "n_eligible_primary_pre_rai": len(analysis_samples),
        "n_unique_eligible_patients": int(sample_sheet["patient_id"].nunique()),
        "group_counts": group_counts,
        "n_platform_probes_with_annotation_rows": int(in_platform["probe_id"].nunique()),
        "n_unambiguous_probe_gene_rows": int(len(usable)),
        "n_gene_symbols_after_iqr_collapse": int(gene_expr.shape[0]),
        "collapse_rule": "highest probe IQR across all 17 eligible samples; no outcome labels used",
    }
    with (QC / "GSE151179_preparation_qc.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(json.dumps(qc, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
