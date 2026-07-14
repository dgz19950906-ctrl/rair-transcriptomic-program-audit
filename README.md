# RAIR transcriptomic program portability audit

Version 1.0.0 release candidate accompanying the manuscript:

> Auditing transcriptomic program portability: a methodological study of reproducibility, endpoint alignment and cohort-structure sensitivity in radioiodine-refractory thyroid cancer

## Contents

- `code/analysis/`: frozen bulk-transcriptomic, null-model, literature-audit and numerical-QC scripts.
- `code/single_cell/`: frozen cloud-execution scripts for the three-atlas single-cell analysis.
- `code/figures/`: canonical builders for the main and sensitivity figures.
- `config/`: locked program definitions, clinical challenge labels and the single-cell execution contract.
- `data/source_data/`: source tables underlying Figures 1–6 and Supplementary Figures S5–S6.
- `data/derived_tables/`: additional derived tables and sensitivity outputs.
- `literature_audit/`: frozen queries, search metadata, record-level screening decisions, flow counts and the 22-study reproducibility audit. Raw database exports containing third-party abstracts are not redistributed; they can be regenerated from the supplied queries and metadata.
- `qc/`: numerical, citation-order and integrity audit records.

No new primary sequencing data are redistributed. The transcriptomic inputs remain available from NCBI Gene Expression Omnibus under GSE151179, GSE299988, GSE184362, GSE191288 and GSE281736.

## Reproduction notes

The original analyses ran on Ubuntu 22.04.5 LTS with R 4.5.1. Principal package versions are recorded in `software_versions.tsv` and in the manuscript Methods. Some frozen scripts retain the original project root (`/home/dony/ThyroidCancer_Project/rair_audit`) to preserve the exact executed workflow; users should replace this root or expose an equivalent directory structure before rerunning.

The intended high-level order is:

1. Recreate the public-data inputs using the GEO accessions and preparation scripts.
2. Score the locked programs and run the two endpoint contrasts.
3. Generate exact patient-label nulls and label-blind covariance-matched program-identity nulls.
4. Run the GSE299988 collinear challenge independently.
5. Run donor-pseudobulk single-cell analyses and the prespecified cross-atlas AND gate.
6. Build figures only from frozen source tables and execute the numerical QC scripts.

Large raw and processed single-cell objects are intentionally excluded from this release. They can be reconstructed from the public GEO records using the supplied scripts.

## Licensing

Original software in `code/` is released under the MIT License (`LICENSE`). Original source and derived tables created for this audit are released under CC BY 4.0 (`DATA_LICENSE.txt`). Third-party datasets, article metadata and externally defined gene sets remain subject to their original providers' terms.

## Citation

Until the accompanying article and repository DOI are available, cite this release using `CITATION.cff`. The final Zenodo DOI and article DOI should be added to a tagged release without rewriting version 1.0.0.

Repository: https://github.com/dgz19950906-ctrl/rair-transcriptomic-program-audit

Version 1.0.0: https://github.com/dgz19950906-ctrl/rair-transcriptomic-program-audit/releases/tag/v1.0.0

## Contact

- Guozhang Dong: dgz_1995@foxmail.com
- Yanbin Zhao (corresponding author): 360348468@qq.com
