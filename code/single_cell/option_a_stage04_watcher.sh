#!/usr/bin/env bash
set -u

AUDIT=/home/dony/ThyroidCancer_Project/rair_audit
STAGE03_OUTPUT=${AUDIT}/checkpoints/03_mnn_clustered_counts.rds
STAGE03_LOG=${AUDIT}/logs/03_fastmnn_author_compartments.log
WATCH_LOG=${AUDIT}/logs/04_watcher.log

date '+[%Y-%m-%d %H:%M:%S %Z] Stage 04 watcher started' >> "${WATCH_LOG}"

while true; do
  if [[ -s "${STAGE03_OUTPUT}" ]] && grep -q 'Stage 03 complete' "${STAGE03_LOG}"; then
    date '+[%Y-%m-%d %H:%M:%S %Z] Stage 03 completion detected; launching Stage 04' >> "${WATCH_LOG}"
    cd "${AUDIT}" || exit 10
    Rscript scripts/04_pseudobulk_lodo_programs.R > logs/04_pseudobulk_lodo_programs.log 2>&1
    status=$?
    date "+[%Y-%m-%d %H:%M:%S %Z] Stage 04 exited with status ${status}" >> "${WATCH_LOG}"
    exit "${status}"
  fi

  if ! pgrep -f 'R.*02_author_qc_doublets.R|R.*03_fastmnn_author_compartments.R|option_a_stage03_watcher.sh' >/dev/null; then
    date '+[%Y-%m-%d %H:%M:%S %Z] Upstream pipeline is not running and no valid Stage 03 checkpoint exists' >> "${WATCH_LOG}"
    exit 20
  fi

  sleep 60
done
