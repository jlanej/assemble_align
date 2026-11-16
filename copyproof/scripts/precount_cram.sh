#!/usr/bin/env bash
# Count the k-mers of CRAM files as they finish downloading, into the databases a later
# `copyproof run --out OUT` will reuse.  precount_cram.sh OUT REF THREADS LABEL=CRAM ...
set -euo pipefail
OUT=${1:?}; REF=${2:?}; T=${3:?}; shift 3
here=$(cd "$(dirname "$0")/.." && pwd); CPK=$here/engine/target/release/cpk
mkdir -p "$OUT"
for spec in "$@"; do
  label=${spec%%=*}; cram=${spec#*=}
  [ -s "$OUT/reads.$label.cpkdb" ] && continue
  until [ -s "$cram" ]; do sleep 20; done
  echo "[$(date '+%T')] counting $label"
  bash -o pipefail -c "samtools fastq -@ 4 -F 0x900 --reference '$REF' '$cram' 2>/dev/null | '$CPK' count --out '$OUT/reads.$label.tmp' --threads $T -" \
    2> "$OUT/reads.$label.count.log"
  for ext in .cpkdb .hist.tsv .json; do mv "$OUT/reads.$label.tmp$ext" "$OUT/reads.$label$ext"; done
  echo "[$(date '+%T')] done $label: $(grep 'count: wrote' "$OUT/reads.$label.count.log" | cut -c1-200)"
done
