#!/usr/bin/env bash
# Slice PRIMARY alignments of a remote, indexed long-read BAM for one or more regions and
# write the read set as FASTQ (reads de-duplicated by name).  Region-sliceable recruitment
# ("recruitment by alignment"): reads are assigned to a region by the aligner used to make the BAM.
# Usage: 01_fetch_longread_slice.sh <BAM_URL> <REGION[,REGION...] | regions.bed> <OUT_PREFIX>
#   e.g. 01_fetch_longread_slice.sh https://.../NA19240_aligned_CHM13v1Y_winnowmap.sorted.bam chr20 out/NA19240.hifi.chr20
# Outputs: OUT_PREFIX.fastq.gz, OUT_PREFIX.log (timing, record counts, index md5)
set -euo pipefail
url=$1; reg=$2; out=$3
mkdir -p "$(dirname "$out")"
tmp=$(mktemp -d "${TMPDIR:-$PWD}/scratch.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
log="$out.log"; : > "$log"
echo "url=$url region=$reg date=$(date -u +%FT%TZ) samtools=$(samtools --version | head -1)" >> "$log"
curl -sSfL --retry 5 -o "$tmp/idx.bai" "$url.bai"
echo "bai_md5=$(md5 -q "$tmp/idx.bai" 2>/dev/null || md5sum "$tmp/idx.bai" | cut -d' ' -f1)" >> "$log"
start=$(date +%s)
if [[ -f "$reg" ]]; then
  # multi-region iterator: each record is returned once even if it overlaps several BED intervals
  samtools view -b -F 0x900 -M -L "$reg" -X "$url" "$tmp/idx.bai" -o "$tmp/slice.bam"
else
  samtools view -b -F 0x900 -M -X "$url" "$tmp/idx.bai" ${reg//,/ } -o "$tmp/slice.bam"
fi
echo "slice_seconds=$(( $(date +%s)-start )) slice_bytes=$(wc -c < "$tmp/slice.bam")" >> "$log"
samtools fastq -n "$tmp/slice.bam" 2>/dev/null | seqkit rmdup -n 2>/dev/null | gzip -1 > "$out.fastq.gz"
seqkit stats -a -T "$out.fastq.gz" | sed 's/^/stats\t/' >> "$log"
cat "$log"
