#!/usr/bin/env bash
# Long-read chr20 candidate prefilter (L2/L3 pilot).
# usage: 03_lr_prefilter.sh P.lr.mmi READS_DIR OUT_DIR [threads] [min_best_hit_qcov]
#   P.lr.mmi : minimap2 (k19,w19) index of {GRCh38 chr20, CHM13 chr20, hap1 chr20, hap2 chr20}
#   READS_DIR: hifi_{1,2,3}.fq.gz and ont_{1,2,3}.fq.gz (streamed subsets)
# A read is a candidate if its single best approximate chain (no base-level alignment)
# on any of the four chr20 sequences covers >= min_best_hit_qcov (default 0.8) of the read.
set -euo pipefail
IDX=$1; R=$2; O=$3; T=${4:-4}; MINCOV=${5:-0.8}
mkdir -p $O
for tech in hifi ont; do
  preset=map-hifi; [ $tech = ont ] && preset=lr:hq
  for i in 1 2 3; do
    minimap2 -t $T -x $preset --secondary=no $IDX $R/${tech}_$i.fq.gz 2> $O/${tech}_$i.P.log | \
      awk -v m=$MINCOV 'BEGIN{OFS="\t"} {f=($4-$3)/$2; if(f>b[$1]) b[$1]=f} END{for(r in b) if(b[r]>=m) print r, b[r]}' > $O/${tech}_$i.cand.tsv
    grep "Real time" $O/${tech}_$i.P.log || true
  done
  cat $O/${tech}_?.cand.tsv | cut -f1 > $O/$tech.cand.ids
  for i in 1 2 3; do seqkit grep -j 2 -f $O/$tech.cand.ids $R/${tech}_$i.fq.gz; done | bgzip -@2 > $O/$tech.cand.fq.gz
  echo "$tech candidates: $(wc -l < $O/$tech.cand.ids) reads"; seqkit stats -T $O/$tech.cand.fq.gz | tail -1
done
