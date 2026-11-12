#!/usr/bin/env bash
# Sequential L4 compute queue (never more than THREADS concurrent threads). Re-entrant: finished steps are skipped.
# Usage: 13_master_queue.sh <WORKDIR> <SCRIPTS_DIR> [THREADS=3]
set -uo pipefail
W=$(cd "$1" && pwd); P=$(cd "$2" && pwd); t=${3:-3}
export OMP_NUM_THREADS=$t NUMBA_NUM_THREADS=$t MKL_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t VECLIB_MAXIMUM_THREADS=$t
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$W/queue.log"; }
G_NA19240=66977342; G_HG01258=66360199
HIFI_NA19240=$W/hifi/NA19240.hifi.fullx.fastq.gz; HIFI_HG01258=$W/hifi/HG01258.hifi.chm13chr20.fastq.gz
ONT_HG01258=$W/ont/HG01258.ontR10.chm13chr20.fastq.gz; UL_HG01258=$W/ont/HG01258.ontR10.ge20kb.fastq.gz
# 0. yak k-mer DBs (chr20 Illumina; child, father, mother)
for S in NA19240 HG01258; do
  if [[ ! -s $W/yak/$S/mat.yak ]]; then
    if [[ $S == NA19240 ]]; then f=NA19239; m=NA19238; else f=HG01256; m=HG01257; fi
    log "yak $S"; "$P/10_yak_dbs.sh" "$W/yak/$S" "$t" "$W/illumina/$S.ilmn.chr20.fastq.gz" "$W/illumina/$f.ilmn.chr20.fastq.gz" "$W/illumina/$m.ilmn.chr20.fastq.gz" >> "$W/queue.log" 2>&1
  fi
done
# 1. truth controls
for S in NA19240 HG01258; do
  od=$W/asm/$S.truth.-; [[ -s $od/eval.jsonl ]] && continue; log "eval truth $S"
  python3 "$P/07_eval_asm.py" --asm "$od" --pat "$W/truth/$S.pat.chr20.fa" --mat "$W/truth/$S.mat.chr20.fa" --child "$W/yak/$S/child.yak" \
    --patyak "$W/yak/$S/pat.yak" --matyak "$W/yak/$S/mat.yak" --label "$S|truth|-x" --out "$W/eval/titration_raw.tsv" --threads "$t" > "$od/eval.jsonl" 2>> "$W/queue.log"
done
run() { log "titration $*"; "$P/11_titration_driver.sh" "$@" >> "$W/queue.log" 2>&1 || log "FAILED $*"; }
# 2. full-coverage reference points for each mode (HG01258), then NA19240 HiFi full
run HG01258 hifi    "$HIFI_HG01258" $G_HG01258 full "$W" "$t"
run HG01258 ont     "$ONT_HG01258"  $G_HG01258 full "$W" "$t"
run HG01258 hifi_ul "$HIFI_HG01258" $G_HG01258 full "$W" "$t" "$UL_HG01258"
# 3. hard loci (HiFi): reference vs population vs oracle recruitment
for S in NA19240 HG01258; do
  [[ -s $W/loci_out/$S.loci.tsv ]] && continue; log "loci $S"
  python3 "$P/09_loci_pipeline.py" --sample $S --pool "$W/loci/$S.hifi.pool.bam" --unmapped "$W/loci/$S.hifi.unmapped.bam" \
    --loci "$W/../meta/loci_grch38.tsv" --proj "$W/loci/truth_projections.tsv" --truthdir "$W/loci/haps_truth,$W/loci/haps_truth_smn" \
    --popdir "$W/loci/haps_panel" --grch38dir "$W/loci/grch38" --gaps "$W/truth/gaps/${S}_pat.gaps.bed,$W/truth/gaps/${S}_mat.gaps.bed" \
    --out "$W/loci_out" --threads "$t" >> "$W/loci_out.$S.log" 2>&1 || log "FAILED loci $S"
done
# 4. lower coverages
run HG01258 hifi    "$HIFI_HG01258" $G_HG01258 24,16,12,8 "$W" "$t"
run HG01258 hifi_ul "$HIFI_HG01258" $G_HG01258 8,16,24,12 "$W" "$t" "$UL_HG01258"
run HG01258 ont     "$ONT_HG01258"  $G_HG01258 16,24,8,12 "$W" "$t"
# 5. NA19240 HiFi replicate
run NA19240 hifi    "$HIFI_NA19240" $G_NA19240 full,8,16,24,12 "$W" "$t"
log "queue finished"
