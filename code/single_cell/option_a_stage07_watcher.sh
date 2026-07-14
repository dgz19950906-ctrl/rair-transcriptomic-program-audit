#!/usr/bin/env bash
set -u

AUDIT=/home/dony/ThyroidCancer_Project/rair_audit
STAGE06_COMPLETE=${AUDIT}/checkpoints/06_COMPLETE
STAGE06_LOG=${AUDIT}/logs/06_prepare_independent_atlases.log
WATCH_LOG=${AUDIT}/logs/07_watcher.log

date '+[%Y-%m-%d %H:%M:%S %Z] Stage 07 watcher started' >> "${WATCH_LOG}"
while true; do
  if [[ -s "${STAGE06_COMPLETE}" ]] && grep -q 'Stage 06 complete' "${STAGE06_LOG}"; then
    date '+[%Y-%m-%d %H:%M:%S %Z] Stage 06 completion detected; launching Stage 07' >> "${WATCH_LOG}"
    cd "${AUDIT}" || exit 10
    Rscript scripts/07_cross_atlas_and_gate.R > logs/07_cross_atlas_and_gate.log 2>&1
    status=$?
    date "+[%Y-%m-%d %H:%M:%S %Z] Stage 07 exited with status ${status}" >> "${WATCH_LOG}"
    exit "${status}"
  fi
  if ! pgrep -f 'R.*06_prepare_independent_atlases.R' >/dev/null; then
    date '+[%Y-%m-%d %H:%M:%S %Z] Stage 06 is not running and no valid completion checkpoint exists' >> "${WATCH_LOG}"
    exit 20
  fi
  sleep 60
done
