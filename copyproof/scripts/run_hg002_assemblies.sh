#!/usr/bin/env bash
# Rung 4: run copyproof on HG002 assemblies with the child's 35x reads and score every run against
# the truth genome's k-mers. Assemblies are processed as their FASTA appears (they may still be
# downloading); finished runs are skipped unless FORCE=1.
#   [FORCE=1] run_hg002_assemblies.sh DATA LIST.tsv [THREADS] [NAME ...]
set -uo pipefail
D=${1:?data root}; LIST=${2:?assemblies.tsv}; T=${3:-10}; shift 3 || true
here=$(cd "$(dirname "$0")/.." && pwd)
PY=$here/.venv/bin/python
DB=$D/kmerdb/HG002.novaseq.35x.cpkdb; TRUTH=$D/kmerdb/hg002v1.1.truth.cpkdb
R1=$D/hg002/reads/HG002.novaseq.pcr-free.30x.R1.fastq.gz; R2=$D/hg002/reads/HG002.novaseq.pcr-free.30x.R2.fastq.gz
W=$D/work/hg002_assemblies; mkdir -p "$W"
grep -v '^#' "$LIST" | while IFS=$'\t' read -r name fasta groups what; do
  [ -n "$name" ] || continue
  if [ $# -gt 0 ] && ! printf '%s\n' "$@" | grep -qx "$name"; then continue; fi
  fa=$D/hg002/gqc/$name/$fasta
  [ -s "$fa" ] || { echo "skip $name: $fa is not there yet"; continue; }
  # FORCE=1 calls and scores again (the scan, the windows and the truth labels are reused)
  [ -z "${FORCE:-}" ] && [ -s "$W/$name.eval.truth_eval.json" ] && { echo "have $name"; continue; }
  args=()
  if [ -n "$groups" ]; then IFS=';' read -ra gs <<< "$groups"; for g in "${gs[@]}"; do args+=(--group-regex "$g"); done; fi
  # sequences that match no pattern (EBV, mitochondrion) go to "other"
  first=other; [ -n "$groups" ] || first=asm
  reads=(); [ -n "${WITH_READS:-}" ] && reads=(--reads "$R1" "$R2")
  [ -n "${WITH_READS:-}" ] || reads=(--no-missing)
  echo "[$(date '+%T')] run $name"
  # ${x[@]+"${x[@]}"}: an empty array must expand to nothing, also under bash 3.2 with set -u
  "$PY" "$here/scripts/run_sample.py" --out "$W/$name" --asm "$first=$fa" ${args[@]+"${args[@]}"} --db sample="$DB" \
     --db truth="$TRUTH" ${reads[@]+"${reads[@]}"} --threads "$T" > "$W/$name.log" 2>&1 \
     || { echo "FAILED $name (see $W/$name.log)"; continue; }
  "$PY" "$here/scripts/eval_truth.py" --run "$W/$name" --out "$W/$name.eval" > "$W/$name.eval.log" 2>&1 || echo "eval FAILED $name"
  # the positional track (about 10 GB) can be rebuilt from the FASTA and the databases
  [ -z "${DROP_TRACKS:-}" ] || rm -f "$W/$name/scan.track.bin"
done
echo "[$(date '+%T')] assemblies pass finished"
