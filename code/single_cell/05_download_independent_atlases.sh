#!/usr/bin/env bash
set -euo pipefail

AUDIT=/home/dony/ThyroidCancer_Project/rair_audit
DATA=${AUDIT}/external_data
MANIFEST=${AUDIT}/manifests/05_external_atlas_sha256.tsv
mkdir -p "${DATA}/GSE191288/raw" "${DATA}/GSE281736/raw"

GSE191_URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE191nnn/GSE191288/suppl/GSE191288_RAW.tar
GSE281_URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281736/suppl/GSE281736_RAW.tar

echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] Downloading GSE191288 from official GEO" 
wget -4 -c --timeout=60 --read-timeout=60 --tries=20 --waitretry=30 \
  -O "${DATA}/GSE191288/GSE191288_RAW.tar" "${GSE191_URL}"
tar -tf "${DATA}/GSE191288/GSE191288_RAW.tar" > "${DATA}/GSE191288/archive_contents.txt"
tar -xf "${DATA}/GSE191288/GSE191288_RAW.tar" -C "${DATA}/GSE191288/raw"

echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] Downloading GSE281736 from official GEO"
wget -4 -c --timeout=60 --read-timeout=60 --tries=20 --waitretry=30 \
  -O "${DATA}/GSE281736/GSE281736_RAW.tar" "${GSE281_URL}"
tar -tf "${DATA}/GSE281736/GSE281736_RAW.tar" > "${DATA}/GSE281736/archive_contents.txt"
tar -xf "${DATA}/GSE281736/GSE281736_RAW.tar" -C "${DATA}/GSE281736/raw"

{
  printf 'sha256\tbytes\tpath\n'
  for file in "${DATA}/GSE191288/GSE191288_RAW.tar" "${DATA}/GSE281736/GSE281736_RAW.tar"; do
    hash=$(sha256sum "${file}" | awk '{print $1}')
    bytes=$(stat -c '%s' "${file}")
    printf '%s\t%s\t%s\n' "${hash}" "${bytes}" "${file}"
  done
} > "${MANIFEST}"

echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] Independent-atlas downloads and extraction complete"
