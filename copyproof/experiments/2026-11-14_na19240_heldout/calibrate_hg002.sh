#!/usr/bin/env bash
# The proxy of the held-out test, measured where the truth is known (NOTEBOOK entry p).
# HG002's HPRC year-2 assembly is run as NA19240 will be (the child's reads, both parents, no long
# reads) and its issues are judged twice: by two other HG002 assemblies, as NA19240's will be, and
# by the truth genome. The primary pair stands in for the two HGSVC3 assemblies: a Verkko assembly
# (v2_trio) and a hifiasm assembly from other reads (hprc_year1). The second pair
# (verkko_asm_thic, hprc_year2_v2_thic) shows how much the proxy moves with the choice of
# assemblies; it sets no prediction.
#   calibrate_hg002.sh DATA [THREADS]
set -euo pipefail
D=${1:?data root}; T=${2:-10}
here=$(cd "$(dirname "$0")/../.." && pwd)
PY=$here/.venv/bin/python; CPK=$here/engine/target/release/cpk
K=$D/kmerdb; G=$D/hg002/gqc; R=$D/hg002/reads; W=$D/work/rung6/hg002_calibration
mkdir -p "$W/hprc_year2"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

frozen=$(cat "$(dirname "$0")/FROZEN_COMMIT")
git -C "$here" diff --quiet "$frozen" -- src engine/src engine/Cargo.toml engine/Cargo.lock scripts/run_sample.py \
  && [ -z "$(git -C "$here" status --porcelain -- src engine/src)" ] \
  || { log "the code differs from the frozen commit $frozen: not running"; exit 1; }
log "code is at the frozen commit $frozen"

# exact k-mer counts of the other assemblies (every sighting counted, no adapter trimming)
exact() {  # name, fasta
  [ -s "$K/hg002.$1.exact.cpkdb" ] && return 0
  log "counting the k-mers of $1"
  "$CPK" count --exact true --out "$K/hg002.$1.exact.tmp" --threads "$T" "$2" 2> "$K/hg002.$1.exact.count.log"
  for ext in .cpkdb .hist.tsv .json; do mv "$K/hg002.$1.exact.tmp$ext" "$K/hg002.$1.exact$ext"; done
}
exact v2_trio "$G/v2_trio/v2_trio.fasta.gz"
exact hprc_year1 "$G/hprc_year1/HG002.f1_assembly_v2_genbank.fa.gz"
exact verkko_asm_thic "$G/verkko_asm_thic/verkko_asm_thic.fasta.gz"
exact hprc_year2_v2_thic "$G/hprc_year2_v2_thic/hprc_year2_v2_thic.fasta.gz"

# the reads recruited for the missing-sequence channel do not depend on the other read sets:
# start from a copy of those of the rung-4 run of the same assembly (re-made if they do not fit)
old=$D/work/hg002_assemblies/hprc_year2/missing
[ -e "$W/hprc_year2/missing" ] || { [ -d "$old" ] && cp -R "$old" "$W/hprc_year2/missing"; }

log "running hprc_year2 with the parents and the other assemblies' k-mers"
"$PY" "$here/scripts/run_sample.py" --out "$W/hprc_year2" --asm other="$G/hprc_year2/HG002.polished.fasta.gz" \
  --group-regex 'hap1=^h1tg' --group-regex 'hap2=^h2tg' \
  --db sample="$K/HG002.novaseq.35x.cpkdb" --db mother="$K/HG004.novaseq.cpkdb" --db father="$K/HG003.novaseq.cpkdb" \
  --db truth="$K/hg002v1.1.truth.cpkdb" --db v2_trio="$K/hg002.v2_trio.exact.cpkdb" --db hprc_year1="$K/hg002.hprc_year1.exact.cpkdb" \
  --db verkko_asm_thic="$K/hg002.verkko_asm_thic.exact.cpkdb" --db hprc_year2_v2_thic="$K/hg002.hprc_year2_v2_thic.exact.cpkdb" \
  --reads "$R/HG002.novaseq.pcr-free.30x.R1.fastq.gz" "$R/HG002.novaseq.pcr-free.30x.R2.fastq.gz" \
  --threads "$T" > "$W/hprc_year2.log" 2>&1
log "scoring"
"$PY" "$here/scripts/eval_assemblies.py" --run "$W/hprc_year2" --refs v2_trio hprc_year1 --truth-label truth \
  --out "$W/hprc_year2.eval" > "$W/hprc_year2.eval.log" 2>&1
"$PY" "$here/scripts/eval_assemblies.py" --run "$W/hprc_year2" --refs verkko_asm_thic hprc_year2_v2_thic --truth-label truth \
  --out "$W/hprc_year2.eval_pair2" > "$W/hprc_year2.eval_pair2.log" 2>&1
"$PY" "$here/scripts/eval_truth.py" --run "$W/hprc_year2" --out "$W/hprc_year2.truth" > "$W/hprc_year2.truth.log" 2>&1
log "calibration finished"
