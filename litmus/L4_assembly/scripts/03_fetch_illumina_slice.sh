#!/usr/bin/env bash
# Slice a remote NYGC 1kGP high-coverage CRAM (GRCh38_full_analysis_set_plus_decoy_hla) for one region and
# write primary, non-duplicate, QC-pass reads as FASTQ (for yak k-mer counting).
# Usage: 03_fetch_illumina_slice.sh <CRAM_URL> <CRAI_URL> <REF_FASTA_with_region_contig> <REGION> <OUT_PREFIX>
# Note: decoding needs only the reference contig(s) of the sliced region (UCSC hg38 chr20 == GRCh38 chr20 up to case).
set -euo pipefail
cram=$1; crai=$2; ref=$3; reg=$4; out=$5
mkdir -p "$(dirname "$out")"
tmp=$(mktemp -d "${TMPDIR:-$PWD}/scratch.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
curl -sSfL --retry 5 -o "$tmp/idx.crai" "$crai"
start=$(date +%s)
samtools view -b -F 0xF00 -T "$ref" -X "$cram" "$tmp/idx.crai" "$reg" -o "$tmp/s.bam"
samtools fastq -n "$tmp/s.bam" 2>/dev/null | gzip -1 > "$out.fastq.gz"
{ echo "cram=$cram crai=$crai region=$reg seconds=$(( $(date +%s)-start )) crai_md5=$(md5 -q "$tmp/idx.crai")";
  seqkit stats -a -T "$out.fastq.gz"; } | tee "$out.log"
