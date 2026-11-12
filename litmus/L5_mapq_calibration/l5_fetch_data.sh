#!/usr/bin/env bash
# L5 MAPQ-calibration pilot: fetch public inputs (NA19240 HPRC r2 haplotypes, GRCh38/CHM13 analysis sets,
# NCBI prebuilt bwa index for GRCh38, impg wfmash PAFs filtered to NA19240, SD annotations).
# Usage: l5_fetch_data.sh OUTDIR
# Streams the 6.5 GB / 5.6 GB impg PAFs through a temp pipe; only NA19240 rows are kept.
set -euo pipefail
OUT=${1:?usage: l5_fetch_data.sh OUTDIR}; mkdir -p "$OUT"; cd "$OUT"
S3=https://s3-us-west-2.amazonaws.com/human-pangenomics
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
G38=GCA_000001405.15_GRCh38_no_alt_analysis_set
get(){ if [ ! -s "$2" ]; then
  curl -sSfL --retry 3 -o "$2.part" "$1"
  mv -f "$2.part" "$2"
fi; }
# 1. references
get $S3/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz chm13v2.0_maskedY_rCRS.fa.gz
get $S3/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz.gzi chm13v2.0_maskedY_rCRS.fa.gz.gzi
get $NCBI/$G38.fna.gz $G38.fna.gz
get $NCBI/md5checksums.txt ncbi_md5checksums.txt
# 2. NA19240 HPRC release-2 haplotypes (hap1 = paternal, hap2 = maternal; hifiasm 0.19.9 trio)
for h in pat mat; do
  for ext in "" .fai .gzi .md5; do
    get $S3/working/HPRC_PLUS/NA19240/assemblies/release2/NA19240_${h}_hprc_r2_v1.0.1.fa.gz$ext NA19240_${h}_hprc_r2_v1.0.1.fa.gz$ext
  done
done
# 3. SD annotations
get https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz hg38.genomicSuperDups.txt.gz
get $S3/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed chm13v2.0_SD.bed
# 4. impg all-vs-1 PAFs (wfmash; query = NA19240#1|2, target = reference), filtered on the fly
P=$S3/pangenomes/freeze/release2/impg/pafs/all-vs-1
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# fast path: BGZF range extraction (rows are grouped by query sample); fallback: full stream
[ -s NA19240.vs_grch38.paf.gz ] || python3 "$HERE/l5_fetch_paf_range.py" $P/$G38.PanSN.merged.paf.gz NA19240 NA19240.vs_grch38.paf.gz || true
[ -s NA19240.vs_chm13.paf.gz ] || python3 "$HERE/l5_fetch_paf_range.py" $P/chm13v2.0_maskedY_rCRS.fa.PanSN.merged.paf.gz NA19240 NA19240.vs_chm13.paf.gz || true
[ -s NA19240.vs_grch38.paf.gz ] || { curl -sSfL --retry 3 $P/$G38.PanSN.merged.paf.gz | gzip -dc | LC_ALL=C grep '^NA19240#' | gzip -c > NA19240.vs_grch38.paf.gz.part && mv NA19240.vs_grch38.paf.gz.part NA19240.vs_grch38.paf.gz; }
[ -s NA19240.vs_chm13.paf.gz ] || { curl -sSfL --retry 3 $P/chm13v2.0_maskedY_rCRS.fa.PanSN.merged.paf.gz | gzip -dc | LC_ALL=C grep '^NA19240#' | gzip -c > NA19240.vs_chm13.paf.gz.part && mv NA19240.vs_chm13.paf.gz.part NA19240.vs_chm13.paf.gz; }
# 5. NCBI prebuilt bwa index for GRCh38 no-alt analysis set (3.4 GB)
get $NCBI/$G38.fna.bwa_index.tar.gz $G38.fna.bwa_index.tar.gz
# 6. md5 checks (NCBI-provided and HPRC-provided)
{ for f in $G38.fna.gz $G38.fna.bwa_index.tar.gz; do
    exp=$(grep -E "/$f\$" ncbi_md5checksums.txt | awk '{print $1}'); got=$(md5 -q $f 2>/dev/null || md5sum $f | cut -d' ' -f1)
    echo -e "$f\t$exp\t$got\t$([ "$exp" = "$got" ] && echo OK || echo MISMATCH)"; done
  for h in pat mat; do f=NA19240_${h}_hprc_r2_v1.0.1.fa.gz; exp=$(awk '{print $1}' $f.md5); got=$(md5 -q $f 2>/dev/null || md5sum $f | cut -d' ' -f1)
    echo -e "$f\t$exp\t$got\t$([ "$exp" = "$got" ] && echo OK || echo MISMATCH)"; done
  for f in chm13v2.0_maskedY_rCRS.fa.gz hg38.genomicSuperDups.txt.gz chm13v2.0_SD.bed NA19240.vs_grch38.paf.gz NA19240.vs_chm13.paf.gz; do
    echo -e "$f\tNA\t$(md5 -q $f 2>/dev/null || md5sum $f | cut -d' ' -f1)\tcomputed"; done; } > md5_report.tsv
echo done
