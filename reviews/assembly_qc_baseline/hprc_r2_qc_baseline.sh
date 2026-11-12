#!/usr/bin/env bash
# Reproduce HPRC release-2 QC baselines (T2T counts, gaps, HMM-Flagger HiFi/ONT, NucFlag) from public index tables.
# Usage: hprc_r2_qc_baseline.sh <out_dir> [repo_raw_base_url] [n_io_workers]
set -euo pipefail
OUT=${1:?out_dir}; BASE=${2:-https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main}; NW=${3:-4}
d=$(mktemp -d "$PWD/scratch.XXXXXX"); trap 'rm -rf "$d"' EXIT
python3 "$(dirname "$0")/hprc_r2_qc_baseline.py" "$BASE" "$d" "$OUT" "$NW"
