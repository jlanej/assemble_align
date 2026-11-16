#!/usr/bin/env bash
# HG002 with both parents' reads (HG004 mother, HG003 father): the truth genome and the planted
# genome again, now with the hapmer channel (S3). On the truth genome every phase switch is a
# false flag; on the planted genome the haplotype swaps are the positives.
#   run_hg002_trio.sh DATA [THREADS]
set -euo pipefail
D=${1:?data root}; T=${2:-8}
here=$(cd "$(dirname "$0")/.." && pwd)
PY=$here/.venv/bin/python; CPK=$here/engine/target/release/cpk
R=$D/hg002/reads; K=$D/kmerdb
for s in HG003 HG004; do
  [ -s "$K/$s.novaseq.cpkdb" ] || "$CPK" count --out "$K/$s.novaseq" --threads "$T" "$R/$s.novaseq.pcr-free.30x.R1.fastq.gz" "$R/$s.novaseq.pcr-free.30x.R2.fastq.gz"
done
dbs=(--db sample="$K/HG002.novaseq.35x.cpkdb" --db mother="$K/HG004.novaseq.cpkdb" --db father="$K/HG003.novaseq.cpkdb" --db truth="$K/hg002v1.1.truth.cpkdb")
reads=(--reads "$R/HG002.novaseq.pcr-free.30x.R1.fastq.gz" "$R/HG002.novaseq.pcr-free.30x.R2.fastq.gz")
W=$D/work/hg002_v1.1; P=$D/work/plant
# the recruited reads do not depend on the parents: reuse them
mkdir -p "$W/run35_trio" "$P/run35_trio"
[ -e "$W/run35_trio/missing" ] || ln -s "$W/run35/missing" "$W/run35_trio/missing"
[ -e "$P/run35_trio/missing" ] || ln -s "$P/run35/missing" "$P/run35_trio/missing"
"$PY" "$here/scripts/run_sample.py" --out "$W/run35_trio" --asm mat="$D/hg002/asm/hg002v1.1.fasta.gz" --group-map "$W/groups.tsv" \
  "${dbs[@]}" "${reads[@]}" --threads "$T" > "$W/run35_trio.log" 2>&1
"$PY" "$here/scripts/run_sample.py" --out "$P/run35_trio" --asm mat="$P/plant1.fa" --group-map "$P/plant1.groups.tsv" \
  "${dbs[@]}" "${reads[@]}" --threads "$T" > "$P/run35_trio.log" 2>&1
"$PY" "$here/scripts/eval_planted.py" --expected "$P/plant1.expected.bed" --edits "$P/plant1.edits.tsv" --run "$P/run35_trio" \
  --baseline "$W/run35_trio" --out "$P/run35_trio.eval" > "$P/run35_trio.eval.log" 2>&1
echo "hg002 trio runs finished"
