#!/usr/bin/env bash
# Stream the first TARGET bases of each long-read file (HPRC release 2, NA19240) into bgzipped FASTQ.
# usage: 01_stream_long_reads.sh OUTDIR TARGET_BP_PER_FILE
# The six streams run concurrently (network-bound; each uses ~1 core for decompression).
# Reads are the FIRST TARGET bases of each file (file order), not a uniform random sample.
set -u
OUT=$1; T=$2; mkdir -p $OUT
B=https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA19240/raw_data
cut_fq(){ awk -v T=$1 'NR%4==2{b+=length($0)} {print} NR%4==0 && b>=T {exit}'; }
i=0
for f in NA19240.m64043_200128_181438.dc.q20.fastq.gz NA19240.m64043_200201_000449.dc.q20.fastq.gz NA19240.m64043_200202_062937.dc.q20.fastq.gz; do
  i=$((i+1)); ( curl -sS -L "$B/PacBio_HiFi/deepconsensus/v1pt2/$f" | gzip -dc 2>/dev/null | cut_fq $T | bgzip -@1 > $OUT/hifi_$i.fq.gz; echo "hifi_$i done $f" ) &
done
i=0
for f in 03_14_23_R1041_NA19240_1_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam 03_14_23_R1041_NA19240_2_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam 03_14_23_R1041_NA19240_3_dorado0.7.2_sup4.1.0_5mCG_5hmCG.bam; do
  i=$((i+1)); ( samtools fastq -n "$B/nanopore/dorado0.7.2_sup4.1.0_5mCG_5hmCG/$f" 2>/dev/null | cut_fq $T | bgzip -@1 > $OUT/ont_$i.fq.gz; echo "ont_$i done $f" ) &
done
wait
for x in $OUT/*.fq.gz; do printf "%s\t" $(basename $x); seqkit stats -T $x | tail -1 | cut -f4-8; done
