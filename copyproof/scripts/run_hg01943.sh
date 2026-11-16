#!/usr/bin/env bash
# Rung 1: HG01943, HPRC release 2. Child and parents from the NYGC 30x CRAMs, long reads through
# HPRC's published HiFi-to-assembly BAM (read remotely, by index), then the R2 regression scorer.
#   run_hg01943.sh DATA [THREADS]
set -euo pipefail
D=${1:?data root}; T=${2:-10}
here=$(cd "$(dirname "$0")/.." && pwd)
OUT=$D/work/hprc_r2/HG01943; mkdir -p "$OUT"
REF=$D/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
A=$D/hprc_r2/HG01943
BAM=https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/HG01943/hprc_r2/assembly_qc/read_alignments/hifi/HG01943.hifi_minimap2_2.28.corrected.bam
for f in "$D/nygc/HG01943.final.cram" "$D/nygc/HG01942.final.cram" "$D/nygc/HG01941.final.cram"; do
  until [ -s "$f" ]; do sleep 20; done
done
"$here/.venv/bin/copyproof" run --out "$OUT" \
  --asm pat="$A/HG01943_pat_hprc_r2_v1.0.1.fa.gz" --asm mat="$A/HG01943_mat_hprc_r2_v1.0.1.fa.gz" \
  --reads "$D/nygc/HG01943.final.cram" --reference "$REF" \
  --parent mother="$D/nygc/HG01942.final.cram" --parent father="$D/nygc/HG01941.final.cram" \
  --long-read-bam "$BAM" --threads "$T" --cards 150 > "$OUT.log" 2>&1
python3 "$here/../litmus/cases_HG01943_HG01103/score_flagged_cases.py" qc-hg01943 --bed "$OUT/issues.bed" --out "$OUT/score_R2.json" \
  > "$OUT.score.log" 2>&1 || true
echo "HG01943 finished"
