#!/usr/bin/env bash
# HG002: count the child's reads once, then run copyproof on the truth genome (rung 3) and on the
# truth genome with planted errors (rung 2). Large files live under DATA (outside git).
#   run_hg002_truth.sh DATA [THREADS]
set -euo pipefail
D=${1:?data root}; T=${2:-10}
here=$(cd "$(dirname "$0")/.." && pwd)
PY=$here/.venv/bin/python; CPK=$here/engine/target/release/cpk
R1=$D/hg002/reads/HG002.novaseq.pcr-free.30x.R1.fastq.gz
R2=$D/hg002/reads/HG002.novaseq.pcr-free.30x.R2.fastq.gz
DB=$D/kmerdb/HG002.novaseq.35x
TRUTH=$D/kmerdb/hg002v1.1.truth
until [ -s "$R2" ]; do sleep 10; done
mkdir -p "$D/kmerdb"
[ -s "$TRUTH.cpkdb" ] || "$CPK" count --exact true --out "$TRUTH" --threads "$T" "$D/hg002/asm/hg002v1.1.fasta.gz"
[ -s "$DB.cpkdb" ] || /usr/bin/time -l "$CPK" count --out "$DB" --threads "$T" "$R1" "$R2" 2> "$DB.count.log"
W=$D/work/hg002_v1.1
"$PY" "$here/scripts/run_sample.py" --out "$W/run35" --asm mat="$D/hg002/asm/hg002v1.1.fasta.gz" --group-map "$W/groups.tsv" \
  --db sample="$DB.cpkdb" --db truth="$TRUTH.cpkdb" --reads "$R1" "$R2" --threads "$T" > "$W/run35.log" 2>&1
"$PY" "$here/scripts/eval_truth.py" --run "$W/run35" --out "$W/run35.eval" > "$W/run35.eval.log" 2>&1
P=$D/work/plant
if [ -s "$P/plant1.fa" ]; then
  "$PY" "$here/scripts/run_sample.py" --out "$P/run35" --asm mat="$P/plant1.fa" --group-map "$P/plant1.groups.tsv" \
    --db sample="$DB.cpkdb" --db truth="$TRUTH.cpkdb" --reads "$R1" "$R2" --threads "$T" > "$P/run35.log" 2>&1
  "$PY" "$here/scripts/eval_planted.py" --expected "$P/plant1.expected.bed" --edits "$P/plant1.edits.tsv" --run "$P/run35" \
    --baseline "$W/run35" --out "$P/run35.eval" \
    > "$P/run35.eval.log" 2>&1
fi
echo "hg002 truth runs finished"
