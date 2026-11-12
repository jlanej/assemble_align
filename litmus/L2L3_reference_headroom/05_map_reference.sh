#!/usr/bin/env bash
# L2/L3 pilot: map the three NA19240 read sets to ONE whole-genome reference with matched
# minimap2 settings; output unsorted BAMs in input order (needed for per-read best-of-two).
# usage: 05_map_reference.sh LABEL REF.fa OUTDIR ILMN_DIR LR_DIR [threads] [keep_index]
#   ILMN_DIR: R1.fq.gz R2.fq.gz (chr20 slice pairs);  LR_DIR: hifi.cand.fq.gz ont.cand.fq.gz
# Indexes are built once per reference (-d) and deleted afterwards unless keep_index=1.
#   sr  index: k21 w11 (minimap2 -x sr)       -> Illumina
#   lr  index: k19 w19 (shared by map-hifi, lr:hq and asm5) -> HiFi (map-hifi), ONT R10 (lr:hq)
set -euo pipefail
L=$1; REF=$2; O=$3; I=$4; LR=$5; T=${6:-4}; KEEP=${7:-0}
mkdir -p $O
log=$O/$L.map.log; : > $log
if [ ! -s $O/$L.sr.mmi ]; then minimap2 -t $T -x sr -d $O/$L.sr.mmi $REF 2>>$log; fi
minimap2 -t $T -ax sr $O/$L.sr.mmi $I/R1.fq.gz $I/R2.fq.gz 2>>$log | samtools view -@1 -b -1 -o $O/ilmn.$L.bam -
[ "$KEEP" = 1 ] || rm -f $O/$L.sr.mmi
if [ ! -s $O/$L.lr.mmi ]; then minimap2 -t $T -x map-hifi -d $O/$L.lr.mmi $REF 2>>$log; fi
minimap2 -t $T -K 200M -ax map-hifi --eqx -Y --secondary=no $O/$L.lr.mmi $LR/hifi.cand.fq.gz 2>>$log | samtools view -@1 -b -1 -o $O/hifi.$L.bam -
minimap2 -t $T -K 200M -ax lr:hq   --eqx -Y --secondary=no $O/$L.lr.mmi $LR/ont.cand.fq.gz  2>>$log | samtools view -@1 -b -1 -o $O/ont.$L.bam -
[ "$KEEP" = 1 ] || rm -f $O/$L.lr.mmi
grep -E "Real time|CMD" $log | sed 's/^/  /'
