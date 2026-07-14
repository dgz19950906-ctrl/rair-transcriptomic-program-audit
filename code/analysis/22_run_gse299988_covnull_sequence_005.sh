#!/usr/bin/env bash
set -euo pipefail

BASE=/home/dony/ThyroidCancer_Project/rair_audit/covariance_null_gse299988_v1
EXPR="$BASE/inputs/GSE299988_tumor_gene_expression.tsv.gz"
REGISTRY="$BASE/inputs/frozen_programs_all9.tsv"
GENERATOR="$BASE/scripts/21_covariance_null_generator_gse299988.R"
VALIDATOR="$BASE/scripts/12_validate_covariance_null_general.R"
WORKERS=8

PROGRAMS=(
  TDS_16
  IODIDE_HANDLING_11
  CONDELLO_2025_SIX
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
  HALLMARK_HYPOXIA
  HALLMARK_ANGIOGENESIS
  HALLMARK_G2M_CHECKPOINT
  HALLMARK_E2F_TARGETS
  HALLMARK_INFLAMMATORY_RESPONSE
)

PROGRESS="$BASE/logs/sequence_0.05_progress.tsv"
if [[ -e "$PROGRESS" ]]; then
  echo "Refusing to overwrite an existing progress log" >&2
  exit 17
fi
printf 'timestamp\tprogram\tevent\tstatus\n' > "$PROGRESS"

for program in "${PROGRAMS[@]}"; do
  result="$BASE/results/GSE299988/$program/tolerance_0.05"
  validation="$BASE/validation/$program/tolerance_0.05"
  if [[ -e "$result" || -e "$validation" ]]; then
    printf '%s\t%s\tpreflight\texisting_output_refused\n' "$(date -Is)" "$program" >> "$PROGRESS"
    exit 18
  fi
  printf '%s\t%s\tgenerator_started\trunning\n' "$(date -Is)" "$program" >> "$PROGRESS"
  Rscript "$GENERATOR" \
    --expr "$EXPR" \
    --signatures "$REGISTRY" \
    --out "$result" \
    --program "$program" \
    --tolerance 0.05 \
    --workers "$WORKERS" \
    > "$BASE/logs/${program}_generator.log" 2>&1
  printf '%s\t%s\tgenerator_completed\tpassed\n' "$(date -Is)" "$program" >> "$PROGRESS"

  Rscript "$VALIDATOR" \
    --result "$result" \
    --signatures "$REGISTRY" \
    --generator "$GENERATOR" \
    --out "$validation" \
    --program "$program" \
    > "$BASE/logs/${program}_validation.log" 2>&1
  printf '%s\t%s\tvalidation_completed\tpassed\n' "$(date -Is)" "$program" >> "$PROGRESS"
done

printf '%s\tALL\tsequence_completed\tpassed\n' "$(date -Is)" >> "$PROGRESS"
touch "$BASE/logs/sequence_0.05.complete"
