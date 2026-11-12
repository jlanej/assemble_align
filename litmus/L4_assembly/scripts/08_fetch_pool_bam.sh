#!/usr/bin/env bash
# Fetch the candidate read pool for locus recruitment from a remote GRCh38-aligned long-read BAM:
# primary+supplementary records overlapping pool.bed (locus windows +/- 1-2 Mb) and all unplaced unmapped reads.
# Usage: 08_fetch_pool_bam.sh <BAM_URL> <pool.bed> <OUT_PREFIX>   -> OUT_PREFIX.pool.bam, OUT_PREFIX.unmapped.bam
set -euo pipefail
url=$1; bed=$2; out=$3
tmp=$(mktemp -d "${TMPDIR:-$PWD}/scratch.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
curl -sSfL --retry 5 -o "$tmp/idx.bai" "$url.bai"
s=$(date +%s)
samtools view -b -F 0x100 -M -L "$bed" -X "$url" "$tmp/idx.bai" -o "$out.pool.bam"
samtools view -b -X "$url" "$tmp/idx.bai" '*' -o "$out.unmapped.bam"
samtools index "$out.pool.bam"
echo "url=$url seconds=$(( $(date +%s)-s )) pool_records=$(samtools view -c "$out.pool.bam") unmapped=$(samtools view -c "$out.unmapped.bam")" | tee "$out.pool.log"
