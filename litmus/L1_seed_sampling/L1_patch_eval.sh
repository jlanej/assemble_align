#!/usr/bin/env bash
# L1_patch_eval.sh -- stock minimap2 v2.31 vs the same source + minimap2_sampling.patch (--mod-oc).
# Builds both on macOS arm64 (make arm_neon=1 aarch64=1), indexes a reference with each preset,
# maps simulated reads with known origin and scores them with paftools.js mapeval.
#
# usage: L1_patch_eval.sh REF.fa FASTQ_DIR PATCH OUT_DIR [THREADS]
#   FASTQ_DIR holds ill.fq, hifi.fq, ont.fq written by L1_make_fastq.py (mapeval-style read names)
# outputs: OUT_DIR/patch_eval.tsv (+ mapeval_<variant>_<preset>.txt)
# requires: curl, make, cc, k8 (for paftools.js), md5/md5sum
set -euo pipefail
REF=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); FQ=$(cd "$2" && pwd); PATCH=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
mkdir -p "$4"; OUT=$(cd "$4" && pwd); T=${5:-1}
URL=https://github.com/lh3/minimap2/archive/refs/tags/v2.31.tar.gz
scr=$(mktemp -d "$PWD/scratch.XXXXXX"); trap 'rm -rf "$scr"' EXIT
curl -sSL -o "$scr/mm2.tgz" "$URL"
( md5 -q "$scr/mm2.tgz" 2>/dev/null || md5sum "$scr/mm2.tgz" | cut -d' ' -f1 ) > "$OUT/minimap2-2.31.tar.gz.md5"
for v in stock modoc; do
  mkdir -p "$scr/$v"; tar xzf "$scr/mm2.tgz" -C "$scr/$v" --strip-components 1 --exclude='*/.gitmodules' 2>/dev/null || true
done
( cd "$scr/modoc" && patch -p1 -s < "$PATCH" )
ARCH=$(uname -m)
for v in stock modoc; do
  if [ "$ARCH" = arm64 ] || [ "$ARCH" = aarch64 ]; then ( cd "$scr/$v" && make -s arm_neon=1 aarch64=1 >/dev/null 2>&1 )
  else ( cd "$scr/$v" && make -s >/dev/null 2>&1 ); fi
done
PAFTOOLS="k8 $scr/stock/misc/paftools.js"
rt() { grep "Real time" "$1" | sed -E 's/.*Real time: ([0-9.]+) sec; CPU: ([0-9.]+) sec; Peak RSS: ([0-9.]+) GB.*/\1\t\2\t\3/'; }
printf "variant\tpreset\treads\tidx_real_s\tidx_cpu_s\tidx_rss_gb\tmmi_bytes\tmap_real_s\tmap_cpu_s\tmap_rss_gb\tmapped\ttop_mapq\tmapped_top\terr_top\tmapped_q1\terr_q1\tmapped_all\terr_all\n" > "$OUT/patch_eval.tsv"
for preset in sr map-ont map-hifi; do
  case $preset in sr) fq=$FQ/ill.fq;; map-ont) fq=$FQ/ont.fq;; map-hifi) fq=$FQ/hifi.fq;; esac
  nreads=$(awk 'NR%4==1' "$fq" | wc -l | tr -d ' ')
  for v in stock modoc; do
    extra=""; [ $v = modoc ] && extra="--mod-oc"
    "$scr/$v/minimap2" -x $preset $extra -t $T -d "$scr/$v.$preset.mmi" "$REF" 2> "$scr/idx.log"
    il=$(rt "$scr/idx.log"); sz=$(wc -c < "$scr/$v.$preset.mmi" | tr -d ' ')
    "$scr/$v/minimap2" -x $preset -t $T "$scr/$v.$preset.mmi" "$fq" > "$scr/$v.$preset.paf" 2> "$scr/map.log"
    ml=$(rt "$scr/map.log")
    mapped=$(cut -f1 "$scr/$v.$preset.paf" | sort -u | wc -l | tr -d ' ')
    $PAFTOOLS mapeval "$scr/$v.$preset.paf" > "$OUT/mapeval_${v}_${preset}.txt"
    # mapeval lines: Q <mapq_bin> <n_in_bin> <err_in_bin> <cumulative_err_rate> <cumulative_mapped>
    qt=$(awk '$1=="Q" {print $2"\t"$6"\t"$5; exit}' "$OUT/mapeval_${v}_${preset}.txt")
    q1=$(awk '$1=="Q" && $2>=1 {l=$6"\t"$5} END{print l}' "$OUT/mapeval_${v}_${preset}.txt")
    qa=$(awk '$1=="Q" {l=$6"\t"$5} END{print l}' "$OUT/mapeval_${v}_${preset}.txt")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $v $preset $nreads "$il" $sz "$ml" $mapped "${qt:-NA	NA	NA}" "${q1:-NA	NA}" "${qa:-NA	NA}" >> "$OUT/patch_eval.tsv"
    rm -f "$scr/$v.$preset.mmi"
  done
done
"$scr/stock/minimap2" --version > "$OUT/minimap2_version.txt"
cat "$OUT/patch_eval.tsv"
