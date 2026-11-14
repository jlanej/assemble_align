#!/usr/bin/env bash
# Score 3 of the held-out test: errors planted in NA19240's release with the design of rung 2
# (NOTEBOOK entry p). The stretch homologous to each edit comes from an alignment of the two
# haplotypes. Edits stay away from satellites, from HPRC's own Flagger and NucFlag calls, and from
# the copy-number mask of the unedited run. The planted genome is run like the unedited one
# (the child's reads and both parents, no long reads) and scored against it.
#   plant.sh DATA [THREADS] [STAGE ...]
# Stages, in this order when none is named: align chain plant run score
set -euo pipefail
D=${1:?data root}; T=${2:-10}; shift 2 2>/dev/null || shift $#
stages=${*:-align chain plant run score}
here=$(cd "$(dirname "$0")/../.." && pwd)
PY=$here/.venv/bin/python
K=$D/kmerdb; W=$D/work/rung6; P=$W/plant; mkdir -p "$P"
A=$D/hprc_r2/NA19240; REF=$D/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
pat=$A/NA19240_pat_hprc_r2_v1.0.1.fa.gz; mat=$A/NA19240_mat_hprc_r2_v1.0.1.fa.gz
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

frozen=$(cat "$(dirname "$0")/FROZEN_COMMIT")
git -C "$here" diff --quiet "$frozen" -- src engine/src engine/Cargo.toml engine/Cargo.lock scripts/run_sample.py scripts/eval_planted.py \
  && [ -z "$(git -C "$here" status --porcelain -- src engine/src)" ] \
  || { log "the code differs from the frozen commit $frozen: not running"; exit 1; }
log "code is at the frozen commit $frozen; stages: $stages"

for stage in $stages; do
  case $stage in
    align)
      if [ ! -s "$P/pat_mat.paf" ]; then
        log "aligning the maternal haplotype to the paternal one"
        minimap2 -cx asm5 -t "$T" "$pat" "$mat" > "$P/pat_mat.paf.tmp" 2> "$P/pat_mat.minimap2.log"
        mv "$P/pat_mat.paf.tmp" "$P/pat_mat.paf"
      fi ;;
    chain)
      "$PY" "$here/scripts/paf_to_chain.py" "$P/pat_mat.paf" > "$P/hap_to_hap.chain" 2> "$P/hap_to_hap.chain.log"
      cat "$P/hap_to_hap.chain.log" ;;
    plant)
      cat "$A"/annot/*.segdups_eichler.bed > "$P/segdups.bed"
      "$PY" "$here/scripts/plant_errors.py" --fasta "$pat" "$mat" --chain "$P/hap_to_hap.chain" --sd "$P/segdups.bed" \
        --avoid "$A"/annot/*.cenSat.bed "$A"/qc/*.bed --assessable "$W/NA19240_short/cannot_assess.bed" \
        --min-length 12000000 --group 'pat=#1#' --group 'mat=#2#' --out "$P/plant1" > "$P/plant1.log" 2>&1
      tail -2 "$P/plant1.log" ;;
    run)
      log "run: the planted genome, short reads of the trio"
      "$PY" "$here/scripts/run_sample.py" --out "$P/run" --asm mat="$P/plant1.fa" --group-map "$P/plant1.groups.tsv" \
        --db sample="$K/NA19240.nygc.cpkdb" --db mother="$K/NA19238.nygc.cpkdb" --db father="$K/NA19239.nygc.cpkdb" \
        --reads "$D/nygc/NA19240.final.cram" --reference "$REF" --threads "$T" > "$P/run.log" 2>&1 ;;
    score)
      "$PY" "$here/scripts/eval_planted.py" --expected "$P/plant1.expected.bed" --edits "$P/plant1.edits.tsv" --run "$P/run" \
        --baseline "$W/NA19240_short" --out "$P/run.eval" > "$P/run.eval.log" 2>&1 ;;
    *) log "unknown stage $stage"; exit 1 ;;
  esac
  log "stage $stage finished"
done
