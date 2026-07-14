#!/usr/bin/env bash
set -u

AUDIT=/home/dony/ThyroidCancer_Project/rair_audit
STAGE02_OUTPUT=${AUDIT}/checkpoints/02_author_qc_singlets.rds
STAGE02_LOG=${AUDIT}/logs/02_author_qc_doublets.log
WATCH_LOG=${AUDIT}/logs/03_watcher.log

date '+[%Y-%m-%d %H:%M:%S %Z] Stage 03 watcher started' >> "${WATCH_LOG}"

while true; do
  if [[ -s "${STAGE02_OUTPUT}" ]] && grep -q 'Stage 02 complete' "${STAGE02_LOG}"; then
    date '+[%Y-%m-%d %H:%M:%S %Z] Stage 02 completion detected; launching Stage 03' >> "${WATCH_LOG}"
    cd "${AUDIT}" || exit 10
    Rscript scripts/03_fastmnn_author_compartments.R > logs/03_fastmnn_author_compartments.log 2>&1
    status=$?
    date "+[%Y-%m-%d %H:%M:%S %Z] Stage 03 exited with status ${status}" >> "${WATCH_LOG}"
    exit "${status}"
  fi

  if ! pgrep -f 'R.*02_author_qc_doublets.R' >/dev/null; then
    date '+[%Y-%m-%d %H:%M:%S %Z] Stage 02 is not running and no valid completion checkpoint exists' >> "${WATCH_LOG}"
    exit 20
  fi

  sleep 60
done
