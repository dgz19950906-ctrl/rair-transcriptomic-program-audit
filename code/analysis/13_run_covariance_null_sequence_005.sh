#!/usr/bin/env bash
set -euo pipefail

base=/home/dony/ThyroidCancer_Project/rair_audit/covariance_null_v2
expr="$base/inputs/GSE151179_primary_preRAI_gene_expression.tsv.gz"
registry="$base/inputs/frozen_programs_all9.tsv"
generator="$base/scripts/11_covariance_null_generator_general.R"
validator="$base/scripts/12_validate_covariance_null_general.R"
workers=8

programs=(
  CONDELLO_2025_SIX
  TDS_16
  IODIDE_HANDLING_11
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
  HALLMARK_HYPOXIA
  HALLMARK_ANGIOGENESIS
  HALLMARK_G2M_CHECKPOINT
  HALLMARK_E2F_TARGETS
  HALLMARK_INFLAMMATORY_RESPONSE
)

progress="$base/logs/sequence_0.05_progress.tsv"
if [[ -e "$progress" ]]; then
  echo "Progress file already exists; refusing to overwrite" >&2
  exit 2
fi
printf 'program\tstage\ttimestamp_utc\n' > "$progress"

for program in "${programs[@]}"; do
  result="$base/results/GSE151179/$program/tolerance_0.05"
  validation="$base/validation/$program/tolerance_0.05"
  mkdir -p "$result" "$validation"
  if [[ -e "$result/tier_manifest.tsv" || -e "$validation/validation_report.tsv" ]]; then
    echo "Existing sentinel for $program; refusing to overwrite" >&2
    exit 3
  fi

  printf '%s\tgenerator_started\t%s\n' "$program" "$(date -u +%FT%TZ)" >> "$progress"
  Rscript "$generator" \
    --expr "$expr" \
    --signatures "$registry" \
    --out "$result" \
    --program "$program" \
    --tolerance 0.05 \
    --workers "$workers" \
    > "$base/logs/${program}_tolerance_0.05_generator.log" 2>&1
  printf '%s\tgenerator_finished\t%s\n' "$program" "$(date -u +%FT%TZ)" >> "$progress"

  Rscript "$validator" \
    --result "$result" \
    --signatures "$registry" \
    --generator "$generator" \
    --out "$validation" \
    --program "$program" \
    > "$base/logs/${program}_tolerance_0.05_validation.log" 2>&1
  printf '%s\tvalidation_passed\t%s\n' "$program" "$(date -u +%FT%TZ)" >> "$progress"
done

printf 'ALL_PROGRAMS\tsequence_complete\t%s\n' "$(date -u +%FT%TZ)" >> "$progress"
