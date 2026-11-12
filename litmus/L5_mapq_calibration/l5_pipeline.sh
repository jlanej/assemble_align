#!/usr/bin/env bash
# L5 driver: minimap2 index builds, mapping runs and feature extraction, strictly one memory-heavy step
# at a time (index builds peak at 14-22 GB RSS for a single-part human index; mapping <= ~11 GB).
# usage: l5_pipeline.sh WORKDIR [THREADS=2] [REFS="chm13 grch38"]
#   WORKDIR must contain data/ (l5_fetch_data.sh), sim/ (l5_simulate.py) and idx/ with bwa indexes
#   named <ref>.{amb,ann,bwt,pac,sa}.
set -euo pipefail
W=$1; T=${2:-2}; REFS=${3:-"chm13 grch38"}
S=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$W"; mkdir -p aln feat
fa(){ case $1 in grch38) echo data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz;; chm13) echo data/chm13v2.0_maskedY_rCRS.fa.gz;; esac; }
sd(){ case $1 in grch38) echo data/hg38.genomicSuperDups.txt.gz;; chm13) echo data/chm13v2.0_SD.full.bed;; esac; }
truth(){ case $1 in bwa_sr|mm2_sr) echo ill;; mm2_hifi) echo hifi;; mm2_lrhq|mm2_ont) echo ont;; esac; }
build(){ [ -s idx/$1.$2.mmi ] || python3 "$S/l5_runrss.py" aln/rss.log "minimap2 -t $T -x $3 -d idx/$1.$2.mmi $(fa $1) 2> idx/$1.$2.log"; }
mapx(){ ref=$1; shift; for run in "$@"; do
          bash "$S/l5_map.sh" $ref $(fa $ref) idx sim aln $T "$run"
          [ -s feat/$ref.$run.feat.tsv.gz ] || python3 "$S/l5_runrss.py" aln/rss.log \
            "python3 $S/l5_extract.py --bam aln/$ref.$run.sam --truth sim/$(truth $run).truth.tsv.gz --ref $ref --sd $(sd $ref) --run $ref.$run --out feat/$ref.$run.feat.tsv.gz 2>> feat/extract.log"
        done; }
for ref in $REFS; do
  build $ref k19w19 map-hifi; mapx $ref mm2_hifi mm2_lrhq
  mapx $ref bwa_sr
  build $ref k21w11 sr;       mapx $ref mm2_sr
  build $ref k15w10 map-ont;  mapx $ref mm2_ont
done
echo pipeline-done
