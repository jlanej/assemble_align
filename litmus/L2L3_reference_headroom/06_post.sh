#!/usr/bin/env bash
# L2/L3 pilot: per-read metrics, best-of-two selection, coordinate sort and allele counts.
# usage (haploid): 06_post.sh hap  TECH LABEL MAPDIR REF.fa SITES.tsv [threads]
#        (diploid): 06_post.sh dip  TECH LABEL MAPDIR H1LABEL H2LABEL H1.fa H2.fa H1SITES H2SITES [threads]
# TECH: ilmn | hifi | ont   (pileup base-quality floor: 20 for ilmn, 10 for long reads)
set -euo pipefail
SD=$(cd "$(dirname "$0")" && pwd)
mode=$1; tech=$2; L=$3; M=$4; shift 4
bq=20; [ "$tech" = ilmn ] || bq=10
if [ "$mode" = hap ]; then
  REF=$1; SITES=$2; T=${3:-4}
  python $SD/l2l3_reads.py --bam $M/$tech.$L.bam --out $M/$tech.$L.npz
  [ -s $M/$tech.$L.sorted.bam.bai ] || { samtools sort -@ $T -m 768M -o $M/$tech.$L.sorted.bam $M/$tech.$L.bam && samtools index $M/$tech.$L.sorted.bam; }
  python $SD/l2l3_pileup.py --bam $M/$tech.$L.sorted.bam --fasta $REF --sites $SITES --minbq $bq --out $M/$tech.$L.ab.tsv.gz
else
  H1=$1; H2=$2; F1=$3; F2=$4; S1=$5; S2=$6; T=${7:-4}
  python $SD/l2l3_reads.py --bam $M/$tech.$H1.bam --bam2 $M/$tech.$H2.bam --out $M/$tech.$L.npz \
      --chosen1 $M/$tech.$L.c1.bam --chosen2 $M/$tech.$L.c2.bam
  for k in 1 2; do
    samtools sort -@ $T -m 768M -o $M/$tech.$L.c$k.sorted.bam $M/$tech.$L.c$k.bam && samtools index $M/$tech.$L.c$k.sorted.bam
    rm -f $M/$tech.$L.c$k.bam
  done
  python $SD/l2l3_pileup.py --bam $M/$tech.$L.c1.sorted.bam --fasta $F1 --sites $S1 --minbq $bq --out $M/$tech.$L.ab1.tsv.gz
  python $SD/l2l3_pileup.py --bam $M/$tech.$L.c2.sorted.bam --fasta $F2 --sites $S2 --minbq $bq --out $M/$tech.$L.ab2.tsv.gz
fi
