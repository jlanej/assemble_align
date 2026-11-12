#!/usr/bin/env bash
# L2/L3 pilot: extract all read pairs with >=1 primary, non-duplicate, non-QC-fail record on
# chr20 from the NYGC 30x CRAM (NA19240), including mates mapped elsewhere and unmapped
# mates, in ONE streaming pass (faster than samtools view --fetch-pairs for scattered mates).
# usage: 02_illumina_slice.sh CRAM(local or URL) GRCh38_full_analysis_set_plus_decoy_hla.fa OUTDIR [threads]
set -euo pipefail
CRAM=$1; D=$2; O=$3; T=${4:-4}
mkdir -p $O
samtools view -@ $((T-1)) -T $D -F 0xF00 -e 'rname=="chr20" || rnext=="chr20"' -b -o $O/chr20.pairs.bam $CRAM
samtools flagstat -@ $((T-1)) $O/chr20.pairs.bam > $O/chr20.pairs.flagstat
samtools collate -@ $((T-1)) -O -u $O/chr20.pairs.bam $O/collate.tmp | \
  samtools fastq -n -F 0x900 -1 $O/R1.fq.gz -2 $O/R2.fq.gz -0 $O/other.fq.gz -s $O/single.fq.gz - 2> $O/fastq.log
seqkit stats -T $O/R1.fq.gz $O/R2.fq.gz $O/other.fq.gz $O/single.fq.gz > $O/fastq.stats
