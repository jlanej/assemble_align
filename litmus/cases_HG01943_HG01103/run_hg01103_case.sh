#!/usr/bin/env bash
# HG01103 flagged case (sub-clonal chromothripsis of the paternal chr2 + chr14 in the lymphoblastoid cell line; 16 junctions).
# Re-tests the event and the 1KGP-PacBio hifiasm assemblies from public data; writes small tables to OUTDIR.
# usage: run_hg01103_case.sh OUTDIR [THREADS=2] [ASM_SCAN=1] [ASM_BAMS=0]
#   ASM_SCAN=1 streams the trio and non-trio hap1/hap2 FASTAs (4 x ~0.85 GB) and scores junction presence (~25 min).
#   ASM_BAMS=1 downloads the unindexed ASM_to_REF trio BAMs (2 x ~1.47 GB) and classifies contigs at every junction end.
# needs on PATH: samtools>=1.20 minimap2 python3 (pysam, numpy, pyahocorasick); bash 3.2 compatible; scratch removed on exit.
set -euo pipefail
OUT=${1:?usage: run_hg01103_case.sh OUTDIR [THREADS] [ASM_SCAN] [ASM_BAMS]}; T=${2:-2}; SCAN=${3:-1}; ABAM=${4:-0}
HERE=$(cd "$(dirname "$0")" && pwd); mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
d=$(mktemp -d "${TMPDIR:-/tmp}/hg01103.XXXXXX"); trap 'rm -rf "$d"' EXIT; cd "$d"
REF=https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
CR=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage
CHILD=$CR/additional_698_related/data/ERR3988840/HG01103.final.cram; FATHER=$CR/data/ERR3241820/HG01101.final.cram; MOTHER=$CR/data/ERR3241821/HG01102.final.cram
B=https://1000g-ont.s3.amazonaws.com/1KGP_PacBio_WGS     # use the virtual-hosted bucket name; path-style s3.amazonaws.com/1000g-ont may be blocked
HIFI=$B/Alignment/GRCh38/HG01103-W-0.GRCh38_no_alt.longphase.haplotagged.bam
TRIO1=$B/Assembly/Trio_Assembly/HG01103-W-0/HG01103-W-0.hifiasm.dip.hap1.p_ctg.gfa.fasta.gz; TRIO2=${TRIO1/hap1/hap2}
REG1=$B/Assembly/Regular_Assembly/HG01103-W-0/HG01103-W-0.hap1.p_ctg.fa.gz; REG2=${REG1/hap1/hap2}
ABAM1=$B/ASM_to_REF_BAM/trio_asm/GRCh38/HG01103-W-0.hap1.bam; ABAM2=${ABAM1/hap1/hap2}
ONT=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38
JB="$HERE/data/hg01103_junctions.bedpe"; SB="$HERE/data/hg01103_segments.bed"
{ echo "date	$(date -u +%FT%TZ)"; for x in samtools minimap2; do echo "$x	$($x --version 2>&1 | head -1)"; done; echo "python	$(python3 --version)"
  for u in REF CHILD FATHER MOTHER HIFI TRIO1 TRIO2 REG1 REG2 ABAM1 ABAM2 ONT; do eval echo "$u	\$$u"; done; } > "$OUT/versions_and_sources.tsv"
# ---- 1. copy number from CRAM indexes (segments >= 100 kb) -------------------------------------------------------
python3 "$HERE/crai_cn.py" --child $CHILD.crai --father $FATHER.crai --mother $MOTHER.crai --regions "$SB" --bed --out "$OUT/crai_cn_segments.tsv"
# ---- 2. reference chr2 + chr14 for CRAM decoding ------------------------------------------------------------------
samtools faidx $REF chr2 chr14 > ref.fa && samtools faidx ref.fa
# ---- 3. HiFi: junction reads, read-derived junction 31-mers (dropping any that occur in GRCh38 chr2/chr14) ---------
grep -v '^#' "$JB" | awk 'BEGIN{OFS="\t"}{print $1,($3>3000?$3-3000:0),$3+3000; print $4,($6>3000?$6-3000:0),$6+3000}' | sort -k1,1 -k2,2n > ends.bed
samtools view -b -L ends.bed -M -o hifi_ends.bam $HIFI && samtools index hifi_ends.bam; rm -f HG01103-W-0.GRCh38_no_alt.longphase.haplotagged.bam.bai
python3 "$HERE/junction_from_reads.py" hifi_ends.bam "$JB" "$OUT/hifi_junction_reads.tsv" kmers_raw.tsv 31
python3 "$HERE/asm_kmer_search.py" kmers_raw.tsv ref_hits.json ref.fa > /dev/null
python3 - <<'PY'
import json; bad={k for k,v in json.load(open('ref_hits.json')).items() if v}
open('kmers.tsv','w').writelines(l for l in open('kmers_raw.tsv') if l.split()[0] not in bad)
PY
cp kmers.tsv "$OUT/hg01103_junction_31mers.tsv"
# ---- 4. short reads of the trio and nanopore reads of the parents at every junction -------------------------------
for x in "child $CHILD" "father $FATHER" "mother $MOTHER"; do set -- $x; python3 "$HERE/sr_junctions_bedpe.py" $2 ref.fa "$JB" $1; done > "$OUT/sr_junctions.jsonl" 2> sr.log
for s in HG01101 HG01102; do python3 "$HERE/lr_junctions_bedpe.py" $ONT/$s.hg38.cram ref.fa "$JB" $s; done > "$OUT/ont_parents_junctions.jsonl" 2> ont.log
rm -f ./*.crai
# ---- 5. haplotagged HiFi depth in kept vs lost segments (coverage signal a read-to-assembly QC would see) -----------
python3 - "$HIFI" > "$OUT/hifi_haplotag_depth.tsv" <<'PY'
import pysam, sys
B = pysam.AlignmentFile(sys.argv[1])
W = [('kept','chr2',120000000,120200000),('kept','chr2',190000000,190200000),('kept','chr14',60000000,60200000),
     ('lost','chr2',145000000,145200000),('lost','chr2',160000000,160200000),('lost','chr14',38000000,38200000)]
print('class\tregion\tHP1_depth\tHP2_depth\tuntagged_depth\tminor_over_major')
for lab, c, s, e in W:
    d = {'1': 0, '2': 0, 'u': 0}
    for r in B.fetch(c, s, e):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
        ov = min(e, r.reference_end) - max(s, r.reference_start)
        if ov > 0: d[str(r.get_tag('HP')) if r.has_tag('HP') else 'u'] += ov
    L = e - s; a, b = sorted([d['1'] / L, d['2'] / L])
    print(f"{lab}\t{c}:{s}-{e}\t{d['1']/L:.1f}\t{d['2']/L:.1f}\t{d['u']/L:.1f}\t{a/b:.2f}")
PY
rm -f HG01103-W-0.GRCh38_no_alt.longphase.haplotagged.bam.bai
# ---- 6. junction presence and context in the assemblies -----------------------------------------------------------
if [ "$SCAN" = 1 ]; then
  python3 "$HERE/score_flagged_cases.py" asm-hg01103 --hap1 $TRIO1 --hap2 $TRIO2 --junction-kmers kmers.tsv --out "$OUT/score_R3_trio_hifiasm.json" > /dev/null
  python3 "$HERE/score_flagged_cases.py" asm-hg01103 --hap1 $REG1 --hap2 $REG2 --junction-kmers kmers.tsv --out "$OUT/score_R3_nontrio_hifiasm.json" > /dev/null
fi
# ---- 7. the published assembly-to-reference alignments (what an alignment-only review sees) ----------------------
if [ "$ABAM" = 1 ]; then
  for h in hap1 hap2; do curl -sS --retry 3 -o trio_$h.bam "$(eval echo \$ABAM${h#hap})"; samtools index trio_$h.bam; md5sum trio_$h.bam 2>/dev/null || md5 trio_$h.bam; done >> "$OUT/versions_and_sources.tsv"
  python3 "$HERE/asm_junction_status.py" "$JB" "$SB" "$OUT/asm_bam_junction_status.tsv" trio_hap1.bam trio_hap2.bam > /dev/null
  for h in hap1 hap2; do echo "$h contigs_with_records $(samtools view trio_$h.bam | cut -f1 | sort -u | wc -l)"; done > "$OUT/asm_bam_contig_counts.txt"
fi
echo "done: $OUT"
