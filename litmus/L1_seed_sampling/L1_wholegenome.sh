#!/usr/bin/env bash
# L1_wholegenome.sh -- whole-genome CHM13v2.0 index: stock minimap2 2.31 vs 2.31 + minimap2_sampling.patch (--mod-oc).
# Measures .mmi size, index build CPU/RSS (split index -I to bound RAM) and, for map-hifi, genome-wide mapping
# accuracy (paftools.js mapeval) of the chr20-simulated HiFi reads with --split-prefix.
# usage: L1_wholegenome.sh chm13v2.0.fa.gz MINIMAP2_STOCK MINIMAP2_MODOC hifi.fq PAFTOOLS_JS OUT_DIR [THREADS] ["map-hifi map-ont"]
# Peak RSS measured on CHM13v2.0 (1 thread): map-hifi -I 1G stock 7.7 GB / mod-oc 5.2 GB.
set -euo pipefail
REF=$1; MS=$2; MM=$3; FQ=$4; PT=$5; mkdir -p "$6"; OUT=$(cd "$6" && pwd); T=${7:-1}; PRESETS=${8:-"map-hifi map-ont"}
scr=$(mktemp -d "$PWD/scratch.XXXXXX"); trap 'rm -rf "$scr"' EXIT
rt() { grep "Real time" "$1" | sed -E 's/.*Real time: ([0-9.]+) sec; CPU: ([0-9.]+) sec; Peak RSS: ([0-9.]+) GB.*/\1\t\2\t\3/'; }
[ -s "$OUT/wholegenome.tsv" ] || printf "preset\tvariant\tsplit\tidx_real_s\tidx_cpu_s\tidx_rss_gb\tmmi_bytes\tmap_real_s\tmap_cpu_s\tmap_rss_gb\tmapped\ttop_mapq\tmapped_top\terr_top\tmapped_q1\terr_q1\tmapped_all\terr_all\n" > "$OUT/wholegenome.tsv"
for preset in $PRESETS; do
  split=1G; [ $preset = map-ont ] && split=500M
  for v in stock modoc; do
    B=$MS; extra=""; [ $v = modoc ] && B=$MM && extra="--mod-oc"
    "$B" -x $preset $extra -I $split -t $T -d "$scr/wg.mmi" "$REF" 2> "$scr/idx.log"
    il=$(rt "$scr/idx.log"); sz=$(wc -c < "$scr/wg.mmi" | tr -d ' ')
    ml="NA	NA	NA"; mp=NA; qt="NA	NA	NA"; q1="NA	NA"; qa="NA	NA"
    if [ $preset = map-hifi ]; then
      "$B" -x $preset -t $T --split-prefix "$scr/sp_$v" "$scr/wg.mmi" "$FQ" > "$scr/wg.paf" 2> "$scr/map.log"
      ml=$(rt "$scr/map.log"); mp=$(cut -f1 "$scr/wg.paf" | sort -u | wc -l | tr -d ' ')
      k8 "$PT" mapeval "$scr/wg.paf" > "$OUT/wg_mapeval_${v}_${preset}.txt"
      qt=$(awk '$1=="Q" {print $2"\t"$6"\t"$5; exit}' "$OUT/wg_mapeval_${v}_${preset}.txt")
      q1=$(awk '$1=="Q" && $2>=1 {l=$6"\t"$5} END{print l}' "$OUT/wg_mapeval_${v}_${preset}.txt")
      qa=$(awk '$1=="Q" {l=$6"\t"$5} END{print l}' "$OUT/wg_mapeval_${v}_${preset}.txt")
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $preset $v $split "$il" $sz "$ml" $mp "$qt" "$q1" "$qa" >> "$OUT/wholegenome.tsv"
    rm -f "$scr/wg.mmi" "$scr/wg.paf"
  done
done
cat "$OUT/wholegenome.tsv"
