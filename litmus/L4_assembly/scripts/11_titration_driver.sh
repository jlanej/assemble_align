#!/usr/bin/env bash
# Coverage titration: nested subsamples -> hifiasm (hifi | hifi_ul | ont) -> evaluation vs truth haplotypes.
# Usage: 11_titration_driver.sh <SAMPLE> <MODE> <READS.fq.gz> <GENOME_SIZE> <COVS e.g. full,24,16,12,8> <WORKDIR> <THREADS> [UL.fq.gz]
#   WORKDIR must contain truth/<SAMPLE>.{pat,mat}.chr20.fa and yak/<SAMPLE>/{child,pat,mat}.yak
set -euo pipefail
export OMP_NUM_THREADS=3 NUMBA_NUM_THREADS=3 MKL_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
S=$1; mode=$2; reads=$3; G=$4; covs=$5; W=$6; t=$7; ul=${8:-}
here=$(cd "$(dirname "$0")" && pwd)
sub="$W/sub/$S.$mode"; mkdir -p "$W/sub" "$W/asm" "$W/eval"
need=""; for c in ${covs//,/ }; do [[ -s "$W/asm/$S.$mode.${c}x/hap1.fa" || -f "$sub.${c}x.fastq.gz" ]] || need="$need,$c"; done
[[ -z "$need" ]] || python3 "$here/05_nested_subsample.py" "$reads" "$G" "$sub" "${need#,}"
for c in ${covs//,/ }; do
  fq="$sub.${c}x.fastq.gz"; [[ -f "$fq" ]] || { echo "skip $c (not reachable)"; continue; }
  od="$W/asm/$S.$mode.${c}x"
  if [[ ! -s "$od/hap1.fa" ]]; then
    if [[ "$mode" == hifi_ul ]]; then "$here/06_run_hifiasm.sh" hifi_ul "$od" "$t" "$fq" "$ul"
    else "$here/06_run_hifiasm.sh" "$mode" "$od" "$t" "$fq"; fi
  fi
  [[ -s "$od/eval.jsonl" ]] || python3 "$here/07_eval_asm.py" --asm "$od" --pat "$W/truth/$S.pat.chr20.fa" --mat "$W/truth/$S.mat.chr20.fa" \
     --child "$W/yak/$S/child.yak" --patyak "$W/yak/$S/pat.yak" --matyak "$W/yak/$S/mat.yak" \
     --label "$S|$mode|${c}x" --out "$W/eval/titration_raw.tsv" --threads "$t" > "$od/eval.jsonl"
  rm -f "$fq"   # subsample no longer needed once assembled+evaluated
  echo "$(date -u +%FT%TZ) done $S $mode ${c}x $(cat "$od/time.txt" | tr '\n' ' ')"
done
