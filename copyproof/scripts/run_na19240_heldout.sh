#!/usr/bin/env bash
# Rung 6, the held-out test: NA19240 in HPRC release 2, with every threshold fixed beforehand
# (NOTEBOOK entry p; record in experiments/2026-11-14_na19240_heldout/).
#
# The trio's reads come from the NYGC 30x CRAMs and are counted once. The two HGSVC3 assemblies of
# the same person ride along as extra k-mer sets: the scoring reads them, the calls do not.
# Two runs of the unedited release: `short` (the child's reads and both parents: the run that is
# scored against the prediction, like the HG002 calibration) and `long` (the same with HPRC's
# HiFi-to-assembly alignment, read remotely).
# A run added to the test (entry p, "a run added by the user"): the hapmer audit on HGSVC3's Verkko
# assembly, which was phased with Strand-seq and not with the parents; HPRC's two trio-phased
# haplotypes ride along as labels, and every call is held against them.
#
#   run_na19240_heldout.sh DATA [THREADS] [STAGE ...]
# Stages, in this order when none is named: count refs short score_short long score_long
# The added run, named explicitly: verkko_labels verkko score_verkko
set -euo pipefail
D=${1:?data root}; T=${2:-10}; shift 2 2>/dev/null || shift $#
stages=${*:-count refs short score_short long score_long}
here=$(cd "$(dirname "$0")/.." && pwd)
PY=$here/.venv/bin/python; CPK=$here/engine/target/release/cpk
E=$here/experiments/2026-11-14_na19240_heldout
K=$D/kmerdb; W=$D/work/rung6; mkdir -p "$W" "$K"
REF=$D/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
A=$D/hprc_r2/NA19240; V=$D/hgsvc3/NA19240/verkko; H=$D/hgsvc3/NA19240/hifiasm
BAM=https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/NA19240/hprc_r2/assembly_qc/read_alignments/hifi/NA19240.hifi_minimap2_2.28.corrected.bam
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

# The code under test is the commit named in the pre-registration: refuse to run anything else.
frozen=$(cat "$E/FROZEN_COMMIT")
git -C "$here" diff --quiet "$frozen" -- src engine/src engine/Cargo.toml engine/Cargo.lock scripts/run_sample.py \
  && [ -z "$(git -C "$here" status --porcelain -- src engine/src)" ] \
  || { log "the code differs from the frozen commit $frozen: not running"; exit 1; }
log "code is at the frozen commit $frozen; stages: $stages"

count_cram() {  # sample: the command `copyproof run` uses for a CRAM, into a database that every run shares
  local out=$K/$1.nygc
  [ -s "$out.cpkdb" ] && { log "have $1"; return 0; }
  # a download is renamed into place only when it is complete: wait for it
  until [ -s "$D/nygc/$1.final.cram" ]; do sleep 30; done
  log "counting $1"
  bash -o pipefail -c "samtools fastq -@ 4 -F 0x900 --reference '$REF' '$D/nygc/$1.final.cram' 2>/dev/null | '$CPK' count --out '$out.tmp' --k 31 --sample 4 --threads $T --genome-size 3.1g -" 2> "$out.count.log"
  for ext in .cpkdb .hist.tsv .json; do mv "$out.tmp$ext" "$out$ext"; done
  log "done $1: $(grep 'count: ' "$out.count.log" | tail -1 | cut -c1-200)"
}

exact() {  # name, FASTA ...: every k-mer sighting of another assembly counted, no adapter trimming
  local out=$K/NA19240.$1.exact; shift
  [ -s "$out.cpkdb" ] && return 0
  log "counting the k-mers of $(basename "$out")"
  "$CPK" count --exact true --out "$out.tmp" --threads "$T" "$@" 2> "$out.count.log"
  for ext in .cpkdb .hist.tsv .json; do mv "$out.tmp$ext" "$out$ext"; done
}

asm=(--asm pat="$A/NA19240_pat_hprc_r2_v1.0.1.fa.gz" --asm mat="$A/NA19240_mat_hprc_r2_v1.0.1.fa.gz")
dbs=(--db sample="$K/NA19240.nygc.cpkdb" --db mother="$K/NA19238.nygc.cpkdb" --db father="$K/NA19239.nygc.cpkdb"
     --db verkko="$K/NA19240.hgsvc3_verkko.exact.cpkdb" --db hifiasm="$K/NA19240.hgsvc3_hifiasm.exact.cpkdb")
reads=(--reads "$D/nygc/NA19240.final.cram" --reference "$REF")

for stage in $stages; do
  case $stage in
    count)
      for s in NA19240 NA19238 NA19239; do count_cram "$s"; done ;;
    refs)
      # a release is every file shipped: both haplotypes, unassigned and disconnected sequence, rDNA, EBV, mitochondrion
      exact hgsvc3_verkko "$V/NA19240.vrk-ps-sseq.asm-hap1.fasta.gz" "$V/NA19240.vrk-ps-sseq.asm-hap2.fasta.gz" \
        "$V/NA19240.vrk-ps-sseq.asm-unassigned.fasta.gz" "$V/NA19240.vrk-ps-sseq.asm-disconnected.fasta.gz" \
        "$V/NA19240.vrk-ps-sseq.asm-rdna.fasta.gz" "$V/NA19240.ps-sseq.exemplar-ebv.fasta.gz" "$V/NA19240.ps-sseq.exemplar-mito.fasta.gz"
      exact hgsvc3_hifiasm "$H/NA19240-hifiasm-hap1.fasta.gz" "$H/NA19240-hifiasm-hap2.fasta.gz" ;;
    short)
      log "run: short reads of the trio"
      "$PY" "$here/scripts/run_sample.py" --out "$W/NA19240_short" "${asm[@]}" "${dbs[@]}" "${reads[@]}" \
        --threads "$T" > "$W/NA19240_short.log" 2>&1 ;;
    score_short)
      "$PY" "$here/scripts/eval_assemblies.py" --run "$W/NA19240_short" --refs verkko hifiasm \
        --out "$W/NA19240_short.eval" > "$W/NA19240_short.eval.log" 2>&1 ;;
    long)
      # the scan, the windows and the recruited reads do not depend on the long reads: share them
      mkdir -p "$W/NA19240_long"
      for f in "$W"/NA19240_short/scan.* "$W"/NA19240_short/windows.* "$W/NA19240_short/unsupported_clusters.tsv" "$W/NA19240_short/missing"; do
        [ -e "$W/NA19240_long/$(basename "$f")" ] || ln -s "$f" "$W/NA19240_long/$(basename "$f")"
      done
      log "run: with HPRC's HiFi alignment"
      "$PY" "$here/scripts/run_sample.py" --out "$W/NA19240_long" "${asm[@]}" "${dbs[@]}" "${reads[@]}" \
        --long-read-bam "$BAM" --cards 150 --threads "$T" > "$W/NA19240_long.log" 2>&1 ;;
    score_long)
      "$PY" "$here/scripts/eval_assemblies.py" --run "$W/NA19240_long" --refs verkko hifiasm \
        --out "$W/NA19240_long.eval" > "$W/NA19240_long.eval.log" 2>&1 ;;
    verkko_labels)
      exact hprc_r2_pat "$A/NA19240_pat_hprc_r2_v1.0.1.fa.gz"
      exact hprc_r2_mat "$A/NA19240_mat_hprc_r2_v1.0.1.fa.gz" ;;
    verkko)
      log "run: the Verkko assembly, hapmer audit"
      "$PY" "$here/scripts/run_sample.py" --out "$W/NA19240_verkko" \
        --asm hap1="$V/NA19240.vrk-ps-sseq.asm-hap1.fasta.gz" --asm hap2="$V/NA19240.vrk-ps-sseq.asm-hap2.fasta.gz" \
        --asm other="$V/NA19240.vrk-ps-sseq.asm-unassigned.fasta.gz" --asm other="$V/NA19240.vrk-ps-sseq.asm-disconnected.fasta.gz" \
        --asm other="$V/NA19240.vrk-ps-sseq.asm-rdna.fasta.gz" --asm other="$V/NA19240.ps-sseq.exemplar-ebv.fasta.gz" \
        --asm other="$V/NA19240.ps-sseq.exemplar-mito.fasta.gz" \
        --db sample="$K/NA19240.nygc.cpkdb" --db mother="$K/NA19238.nygc.cpkdb" --db father="$K/NA19239.nygc.cpkdb" \
        --db hprc_pat="$K/NA19240.hprc_r2_pat.exact.cpkdb" --db hprc_mat="$K/NA19240.hprc_r2_mat.exact.cpkdb" \
        --no-missing --threads "$T" > "$W/NA19240_verkko.log" 2>&1 ;;
    score_verkko)
      "$PY" "$E/hapmers_against_hprc.py" "$W/NA19240_verkko" "$W/NA19240_verkko.hapmers" \
        --labels mother=hprc_mat father=hprc_pat > "$W/NA19240_verkko.hapmers.log" 2>&1 ;;
    *) log "unknown stage $stage"; exit 1 ;;
  esac
  log "stage $stage finished"
done
