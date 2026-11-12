#!/usr/bin/env bash
# L1_fetch_data.sh -- fetch all public inputs of litmus pilot L1 into DATA_DIR and record md5s.
# usage: L1_fetch_data.sh DATA_DIR [THREADS] [--keep-genome]
#   chm13_chr20.fa, hg38_chr20.fa, NA19240_{pat,mat}_chr20.fa (HPRC release 2, remote samtools faidx),
#   NA19240_{pat,mat}_vs_chm13_chr20.paf (minimap2 -x asm5 -c), GreedyMini orders (gm/*.txt),
#   minimap2-2.31.tar.gz; with --keep-genome also chm13v2.0.fa.gz(+.fai/.gzi) for L1_wholegenome.sh.
# requires: curl, samtools (>=1.10, with libcurl), minimap2 2.31, md5 or md5sum
set -euo pipefail
mkdir -p "$1"; D=$(cd "$1" && pwd); T=${2:-2}; KEEP=${3:-}
scr=$(mktemp -d "$PWD/scratch.XXXXXX"); trap 'rm -rf "$scr"' EXIT
md5f() { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }
B=https://s3-us-west-2.amazonaws.com/human-pangenomics
: > "$D/MD5SUMS.txt"
# CHM13v2.0 (bgzipped, with .gzi): download, index, extract chr20
curl -sSL -o "$scr/chm13v2.0.fa.gz" $B/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
curl -sSL -o "$scr/chm13v2.0.fa.gz.gzi" $B/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.gzi
echo "$(md5f "$scr/chm13v2.0.fa.gz")  chm13v2.0.fa.gz  $B/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" >> "$D/MD5SUMS.txt"
samtools faidx "$scr/chm13v2.0.fa.gz"
samtools faidx "$scr/chm13v2.0.fa.gz" chr20 > "$D/chm13_chr20.fa"
if [ "$KEEP" = "--keep-genome" ]; then mv "$scr"/chm13v2.0.fa.gz* "$D/"; fi
# GRCh38 chr20 (UCSC)
curl -sSL -o "$scr/hg38_chr20.fa.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
echo "$(md5f "$scr/hg38_chr20.fa.gz")  hg38_chr20.fa.gz  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz" >> "$D/MD5SUMS.txt"
gzip -dc "$scr/hg38_chr20.fa.gz" > "$D/hg38_chr20.fa"
# NA19240 HPRC release-2 haplotypes: chr20 = CM099588.1 (pat, GCA_018503265.2) / CM099626.1 (mat, GCA_018503275.2)
# (chromosome assignment from the NCBI assembly reports; see L1_sampling_schemes.md)
A=$B/working/HPRC_PLUS/NA19240/assemblies/release2
( cd "$scr" && samtools faidx $A/NA19240_pat_hprc_r2_v1.0.1.fa.gz "NA19240#1#CM099588.1" > "$D/NA19240_pat_chr20.fa" )
( cd "$scr" && samtools faidx $A/NA19240_mat_hprc_r2_v1.0.1.fa.gz "NA19240#2#CM099626.1" > "$D/NA19240_mat_chr20.fa" )
for f in chm13_chr20.fa hg38_chr20.fa NA19240_pat_chr20.fa NA19240_mat_chr20.fa; do
  samtools faidx "$D/$f"; echo "$(md5f "$D/$f")  $f" >> "$D/MD5SUMS.txt"; done
# liftover alignments
for h in pat mat; do
  minimap2 -x asm5 -c --secondary=no -t "$T" "$D/chm13_chr20.fa" "$D/NA19240_${h}_chr20.fa" > "$D/NA19240_${h}_vs_chm13_chr20.paf" 2> "$scr/mm2.log"
done
# GreedyMini binary orders (OrensteinLab/Minimizers, lowest-density/binary/individual-orders/ranks)
mkdir -p "$D/gm"
for o in k15-w10 k15-w19 k16-w11 k15-w5; do
  curl -sSL -o "$D/gm/$o.txt" "https://raw.githubusercontent.com/OrensteinLab/Minimizers/main/lowest-density/binary/individual-orders/ranks/$o.txt"
  echo "$(md5f "$D/gm/$o.txt")  gm/$o.txt" >> "$D/MD5SUMS.txt"; done
# minimap2 2.31 source (for the sketch.c harness and the patch)
curl -sSL -o "$D/minimap2-2.31.tar.gz" https://github.com/lh3/minimap2/archive/refs/tags/v2.31.tar.gz
echo "$(md5f "$D/minimap2-2.31.tar.gz")  minimap2-2.31.tar.gz" >> "$D/MD5SUMS.txt"
cat "$D/MD5SUMS.txt"
