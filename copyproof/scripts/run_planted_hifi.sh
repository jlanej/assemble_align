#!/usr/bin/env bash
# Rung 2 with long reads. No long reads are aligned to the planted genome, so an alignment is made
# from real ones: the HiFi reads that the public alignment places near each edit of the unedited
# genome are fetched by region and aligned again to the planted genome. Then the pipeline runs on
# the planted genome with --long-read-bam, reusing the scan, the windows, the unsupported runs and
# the recruited reads of the short-read run, and the audit is scored. It runs twice: with the
# alignment as minimap2 writes it, and with its secondary records removed.
#   run_planted_hifi.sh DATA [THREADS]
# Needs scripts/run_hg002_truth.sh to have run (it makes DATA/work/plant/run35), samtools and
# minimap2. JOBS sets the number of regions fetched at a time (6). BEFORE_HEAVY, if set, is a
# command run before each step that reads the genome or a track from DATA (on a disk shared with
# other jobs: a script that waits for its turn).
set -euo pipefail
D=${1:?data root}; T=${2:-10}
here=$(cd "$(dirname "$0")/.." && pwd)
PY=$here/.venv/bin/python
X=$here/experiments/2026-11-13_planted_errors
P=$D/work/plant; H=$P/hifi
# HiFi reads (Revio) aligned to hg002v1.1 by the Q100 project: resources.tsv, hg002_q100_hifi_bam
BAM=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/hifi_revio_pbmay24/hg002v1.1_hifi_revio_pbmay24.bam
FLANK=50000
heavy() { [ -z "${BEFORE_HEAVY:-}" ] || "$BEFORE_HEAVY"; }
mkdir -p "$H"

# 1. the stretches to fetch, in the coordinates of the unedited genome: every edit and its
#    homologous stretch, with FLANK bases on both sides
"$PY" "$X/long_read_audit.py" regions --edits "$P/plant1.edits.tsv" --fai "$D/hg002/asm/hg002v1.1.fasta.gz.fai" \
  --flank "$FLANK" --out "$H/regions.bed"

# 2. every alignment record of the public file in those stretches, one request per stretch and
#    JOBS at a time (the stretches lie megabases apart, so no record is fetched twice). samtools
#    keeps the index in the working directory; a cached copy is used when there is one. A part
#    that failed is fetched again when the script is run again.
if [ ! -s "$H/q100.bam" ]; then
  mkdir -p "$H/parts"
  [ -s "$H/$(basename "$BAM").bai" ] || cp "$D/work/hg002_v1.1/bam_cache/$(basename "$BAM").bai" "$H/" 2>/dev/null \
    || (cd "$H" && samtools idxstats "$BAM" > /dev/null)
  awk '{printf "%04d %s:%d-%d\n", NR, $1, $2 + 1, $3}' "$H/regions.bed" | (cd "$H" && xargs -P "${JOBS:-6}" -L 1 sh -c \
    '[ -s "parts/$1.bam" ] || { samtools view -b -o "parts/$1.tmp" "$0" "$2" 2>> fetch.log && mv "parts/$1.tmp" "parts/$1.bam"; }' "$BAM")
  if [ "$(ls "$H/parts" | grep -c '\.bam$')" -ne "$(wc -l < "$H/regions.bed")" ] || ! samtools quickcheck "$H"/parts/*.bam; then
    echo "the fetch is incomplete (see $H/fetch.log); run again" >&2; exit 1
  fi
  samtools cat -o "$H/q100.tmp.bam" "$H"/parts/*.bam
  mv "$H/q100.tmp.bam" "$H/q100.bam"
  rm -r "$H/parts"
fi

# 3. the reads themselves (a primary record holds the whole read), aligned to the planted genome;
#    then the same alignment without its secondary records
[ -s "$H/reads.fastq.gz" ] || { samtools fastq -F 0x900 "$H/q100.bam" 2> "$H/fastq.log" | gzip -1 > "$H/reads.tmp.gz" && mv "$H/reads.tmp.gz" "$H/reads.fastq.gz"; }
if [ ! -s "$H/plant1.hifi.bam.bai" ]; then
  heavy
  minimap2 -ax map-hifi -t "$T" "$P/plant1.fa" "$H/reads.fastq.gz" 2> "$H/minimap2.log" \
    | samtools sort -@ 4 -m 1G -T "$H/sort" -o "$H/plant1.hifi.bam" -
  samtools index "$H/plant1.hifi.bam"
fi
if [ ! -s "$H/plant1.hifi.no_secondary.bam.bai" ]; then
  samtools view -b -F 0x100 -o "$H/plant1.hifi.no_secondary.bam" "$H/plant1.hifi.bam"
  samtools index "$H/plant1.hifi.no_secondary.bam"
fi
"$PY" "$X/long_read_audit.py" reads --regions "$H/regions.bed" --fetched "$H/q100.bam" --aligned "$H/plant1.hifi.bam" \
  --out "$P/run35_hifi.reads.tsv" > "$P/run35_hifi.reads.log" 2>&1

# 4. the pipeline on the planted genome with each alignment, and 5. the scores: the audit at every
#    planted false adjacency, and the scoring of rung 2 again. The scan, the windows, the list of
#    unsupported runs and the recruited reads are links into the short-read run. The located
#    copies are copied: that stage writes its two files again when the segments have changed, and
#    through links it would leave the short-read run with one new file and one old.
for name in run35_hifi run35_hifi_no_secondary; do
  R=$P/$name
  aln=$H/plant1.hifi.bam; [ "$name" = run35_hifi ] || aln=$H/plant1.hifi.no_secondary.bam
  mkdir -p "$R"
  for f in "$P"/run35/scan.* "$P"/run35/windows.w1000.* "$P/run35/unsupported_clusters.tsv" "$P/run35/missing"; do
    ln -sfn "$f" "$R/$(basename "$f")"
  done
  [ -s "$R/copies.hits.tsv" ] || cp "$P/run35/copies.regions.bed" "$P/run35/copies.hits.tsv" "$R/"
  heavy
  "$PY" "$here/scripts/run_sample.py" --out "$R" --asm mat="$P/plant1.fa" --group-map "$P/plant1.groups.tsv" \
    --db sample="$D/kmerdb/HG002.novaseq.35x.cpkdb" --db truth="$D/kmerdb/hg002v1.1.truth.cpkdb" \
    --reads "$D/hg002/reads/HG002.novaseq.pcr-free.30x.R1.fastq.gz" "$D/hg002/reads/HG002.novaseq.pcr-free.30x.R2.fastq.gz" \
    --long-read-bam "$aln" --threads "$T" > "$R.log" 2>&1
  "$PY" "$X/long_read_audit.py" score --edits "$P/plant1.edits.tsv" --expected "$P/plant1.expected.bed" --regions "$H/regions.bed" \
    --run "$R" --short-run "$P/run35" --truth-run "$D/work/hg002_v1.1/run35_hifi" --out "$R.audit" > "$R.audit.log" 2>&1
  "$PY" "$X/long_read_audit.py" records --edits "$P/plant1.edits.tsv" --expected "$P/plant1.expected.bed" --bam "$aln" \
    --out "$R.audit" > "$R.records.log" 2>&1
  "$PY" "$here/scripts/eval_planted.py" --expected "$P/plant1.expected.bed" --edits "$P/plant1.edits.tsv" --run "$R" \
    --baseline "$D/work/hg002_v1.1/run35" --out "$R.eval" > "$R.eval.log" 2>&1
done
echo "planted genome with long reads finished"
