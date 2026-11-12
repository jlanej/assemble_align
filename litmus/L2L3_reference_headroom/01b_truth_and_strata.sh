#!/usr/bin/env bash
# L2/L3 pilot: assembly-based chr20 truth (dipcall-like) + GRCh38 strata + per-reference site tables.
# usage: 01b_truth_and_strata.sh WORKDIR [threads]
# hap1/hap2 chr20 = NA19240#1#CM099588.1 / NA19240#2#CM099626.1 (NCBI assembly reports of
# GCA_018503265.2 / GCA_018503275.2, 'assembled-molecule 20').
set -euo pipefail
S=$1; T=${2:-4}; SD=$(cd "$(dirname "$0")" && pwd)
mkdir -p $S/truth $S/strata $S/sites $S/map
cd $S/truth
samtools faidx $S/ref/GRCh38.fa chr20 > GRCh38.chr20.fa
samtools faidx $S/asm/hap1.fa "NA19240#1#CM099588.1" > hap1.chr20.fa
samtools faidx $S/asm/hap2.fa "NA19240#2#CM099626.1" > hap2.chr20.fa
for R in GRCh38 CHM13; do
  [ -s $S/map/$R.lr.mmi ] || minimap2 -t $T -x map-hifi -d $S/map/$R.lr.mmi $S/ref/$R.fa   # k19/w19 = asm5 index
  for h in hap1 hap2; do
    minimap2 -t $T -c -x asm5 --cs -z200000,10000 -r2k --secondary=no $S/map/$R.lr.mmi $h.chr20.fa 2> $h.$R.log | gzip > $h.$R.paf.gz
  done
done
python $SD/l2l3_truth.py --ref $S/ref/GRCh38.fa --chrom chr20 --h1 hap1.GRCh38.paf.gz --h2 hap2.GRCh38.paf.gz \
    --h1c hap1.CHM13.paf.gz --h2c hap2.CHM13.paf.gz --out truth.chr20
python $SD/l2l3_strata.py --fasta $S/ref/GRCh38.fa --sd $S/ann/hg38.genomicSuperDups.txt.gz \
    --srep $S/ann/hg38.simpleRepeat.txt.gz --cen $S/ann/hg38.centromeres.txt.gz --outdir $S/strata
python $SD/l2l3_sites.py --truth truth.chr20.het_snv.tsv.gz --strata $S/strata --outdir $S/sites \
    $(for c in pers1KG_1 pers1KG_2 pers2p_A pers2p_R; do [ -s $S/pers/$c.chain ] && echo "--chain $c=$S/pers/$c.chain"; done)
