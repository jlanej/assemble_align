#!/usr/bin/env bash
# L5: map simulated reads with bwa mem 0.7.19 and minimap2 2.31 presets; SAM output (+ peak RSS log).
# usage: l5_map.sh REFNAME REF_FASTA IDXDIR SIMDIR OUTDIR [THREADS=2] [RUNS="bwa_sr mm2_sr mm2_hifi mm2_lrhq mm2_ont"]
#   IDXDIR must hold REFNAME.k21w11.mmi, REFNAME.k19w19.mmi, REFNAME.k15w10.mmi and a bwa index REFNAME.*
set -euo pipefail
REF=$1; FA=$2; IDX=$3; SIM=$4; OUT=$5; T=${6:-2}; RUNS=${7:-"bwa_sr mm2_sr mm2_hifi mm2_lrhq mm2_ont"}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); mkdir -p "$OUT"; LOG=$OUT/rss.log
for run in $RUNS; do
  o=$OUT/$REF.$run.sam; [ -s "$o" ] && { echo "skip $o"; continue; }
  case $run in
    bwa_sr)   cmd="bwa mem -t $T -K 10000000 $IDX/$REF $SIM/ill_R1.fq.gz $SIM/ill_R2.fq.gz" ;;
    mm2_sr)   cmd="minimap2 -t $T -K 100M -ax sr --secondary=yes $IDX/$REF.k21w11.mmi $SIM/ill_R1.fq.gz $SIM/ill_R2.fq.gz" ;;
    mm2_hifi) cmd="minimap2 -t $T -K 100M -ax map-hifi $IDX/$REF.k19w19.mmi $SIM/hifi.fa.gz" ;;
    mm2_lrhq) cmd="minimap2 -t $T -K 100M -ax lr:hq $IDX/$REF.k19w19.mmi $SIM/ont.fa.gz" ;;
    mm2_ont)  cmd="minimap2 -t $T -K 100M -ax map-ont $IDX/$REF.k15w10.mmi $SIM/ont.fa.gz" ;;
  esac
  python3 "$HERE/l5_runrss.py" "$LOG" "$cmd > $o.part 2> $OUT/$REF.$run.log" && mv "$o.part" "$o"
done
