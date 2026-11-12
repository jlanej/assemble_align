#!/usr/bin/env bash
# L2/L3 reference-headroom pilot (NA19240, chr20): end-to-end driver.
# usage: run_all.sh OUTDIR [threads=4] [keep_scratch=0]
# All large intermediates (references, CRAM, reads, BAMs, indexes; ~120 GB peak) live in a
# mktemp scratch directory under $PWD that is deleted on exit unless keep_scratch=1.
# Only small result tables/figures are copied to OUTDIR.
# Tools (2026-11-11 run): minimap2 2.31-r1302, samtools/bcftools/htslib 1.25, seqkit, pysam 0.24.1,
# pytrf 1.5.0, numpy 2.4.6, pandas 3.1.1, matplotlib (see ../PROVENANCE.tsv).
# Resource budget used: <=4 threads, ~14 GB peak RSS (minimap2 sr index build 13.9 GB).
set -euo pipefail
OUT=$(mkdir -p "$1" && cd "$1" && pwd); T=${2:-4}; KEEP=${3:-0}
SD=$(cd "$(dirname "$0")" && pwd)
S=$(mktemp -d "$PWD/scratch.l2l3.XXXXXX")
[ "$KEEP" = 1 ] || trap 'rm -rf "$S"' EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $S/logs $S/map $S/lr $S/pers $S/eval
bash $SD/00_fetch.sh $S                                          # references, assembly, panel, CRAM, long reads
bash $SD/01b_truth_and_strata.sh $S $T                           # truth het SNVs, strata, site tables (GRCh38/CHM13/haps)
bash $SD/02_illumina_slice.sh $S/ilmn/NA19240.final.cram $S/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa $S/ilmn $T
# chr20-union prefilter index and long-read candidates
samtools faidx $S/ref/CHM13.fa chr20 | sed 's/^>chr20/>CHM13_chr20/' > $S/lr/P.fa
sed 's/^>chr20/>GRCh38_chr20/' $S/truth/GRCh38.chr20.fa >> $S/lr/P.fa
cat $S/truth/hap1.chr20.fa $S/truth/hap2.chr20.fa >> $S/lr/P.fa
minimap2 -t $T -x map-hifi -d $S/lr/P.lr.mmi $S/lr/P.fa
bash $SD/03_lr_prefilter.sh $S/lr/P.lr.mmi $S/reads $S/lr $T 0.8
# GRCh38 mapping first (needed for the two-pass personalization)
bash $SD/05_map_reference.sh GRCh38 $S/ref/GRCh38.fa $S/map $S/ilmn $S/lr $T 0
samtools sort -@ $T -m 768M -o $S/map/ilmn.GRCh38.sorted.bam $S/map/ilmn.GRCh38.bam && samtools index $S/map/ilmn.GRCh38.sorted.bam
bash $SD/04_two_pass_personalize.sh panel $S/ref/GRCh38.fa $S/vcf/panel.chr20.vcf.gz $S/pers
bash $SD/04_two_pass_personalize.sh twopass $S/ref/GRCh38.fa $S/map/ilmn.GRCh38.sorted.bam $S/pers $T
python $SD/l2l3_sites.py --truth $S/truth/truth.chr20.het_snv.tsv.gz --strata $S/strata --outdir $S/sites \
   --chain pers1KG_1=$S/pers/pers1KG_1.chain --chain pers1KG_2=$S/pers/pers1KG_2.chain \
   --chain pers2p_A=$S/pers/pers2p_A.chain --chain pers2p_R=$S/pers/pers2p_R.chain
# all references, sequentially (one minimap2 at a time)
bash $SD/07_run_refs.sh $S $SD $T GRCh38 CHM13 oracle pers1KG pers2p
# summaries, long-read error profiles, figures
bash $SD/08_summarize.sh $S $SD $OUT
