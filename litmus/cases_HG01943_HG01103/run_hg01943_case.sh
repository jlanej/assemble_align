#!/usr/bin/env bash
# HG01943 flagged case (chr2:46,970,725-48,057,041 maternal de novo tandem duplication; HPRC release-2 assembly).
# Re-tests the call and the assembly from public, region-sliced data and writes small tables to OUTDIR.
# usage: run_hg01943_case.sh OUTDIR [THREADS=2] [WHOLE_ASM_SCAN=0|1]
#   WHOLE_ASM_SCAN=1 streams both HPRC haplotype FASTAs (1.77 GB) to search the junction k-mers genome-wide (~10 min).
# needs on PATH: samtools>=1.20 bcftools minimap2 paftools.js(k8) meryl meryl-lookup yak hifiasm python3 (pysam, numpy, pyahocorasick)
# bash 3.2 compatible; all large intermediates live in a mktemp scratch removed on exit.
set -euo pipefail
OUT=${1:?usage: run_hg01943_case.sh OUTDIR [THREADS] [WHOLE_ASM_SCAN]}; T=${2:-2}; WHOLE=${3:-0}
HERE=$(cd "$(dirname "$0")" && pwd); mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
d=$(mktemp -d "${TMPDIR:-/tmp}/hg01943.XXXXXX"); trap 'rm -rf "$d"' EXIT; cd "$d"
export OMP_NUM_THREADS=$T
# ---- public sources ----------------------------------------------------------------------------------------------
REF=https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
CR=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage
CHILD=$CR/additional_698_related/data/ERR3988955/HG01943.final.cram; FATHER=$CR/data/ERR3242006/HG01941.final.cram; MOTHER=$CR/data/ERR3242007/HG01942.final.cram
HIFI=https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/samples/HG01943/alignments/GRCh38_no_alt/hifi_reads/HG01943.GRCh38_no_alt.sorted.bam
ASM=https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01943/assemblies/release2
WM=$ASM/alignments/winnowmap/GRCh38_no_alt
QC=https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions
FL=$QC/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/HG01943/hprc_r2/assembly_qc/hmm_flagger
NF=$QC/5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA--HPRC_R2_QC_NUCFLAG/HG01943/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi
WIN=chr2:46400001-48650000; LWIN=chr2:46790001-48240000
{ echo "date	$(date -u +%FT%TZ)"; for x in samtools bcftools minimap2 hifiasm; do echo "$x	$($x --version 2>&1 | head -1)"; done
  echo "meryl	$(meryl --version 2>&1 | head -1)"; echo "yak	$(yak version 2>&1 | head -1)"; echo "python	$(python3 --version)"
  for u in REF CHILD FATHER MOTHER HIFI ASM WM FL NF; do eval echo "$u	\$$u"; done; } > "$OUT/versions_and_sources.tsv"
curl -sS "$ASM/HG01943_pat_hprc_r2_v1.0.1.fa.gz.md5" >> "$OUT/versions_and_sources.tsv" || true
curl -sS "$ASM/HG01943_mat_hprc_r2_v1.0.1.fa.gz.md5" >> "$OUT/versions_and_sources.tsv" || true
# ---- 1. copy number from CRAM indexes ---------------------------------------------------------------------------
printf 'left_flank\tchr2\t45900000\t46900000\nsegment\tchr2\t46970725\t48057041\nright_flank\tchr2\t48130000\t49130000\n' > regions.tsv
python3 "$HERE/crai_cn.py" --child $CHILD.crai --father $FATHER.crai --mother $MOTHER.crai --regions regions.tsv --norm chr2:45900000-46900000 --out "$OUT/crai_cn.tsv"
# ---- 2. short reads over the window (CRAM decode needs GRCh38 chr2 of the analysis set) -------------------------
samtools faidx $REF chr2 > chr2.fa && samtools faidx chr2.fa
for s in child father mother; do eval u=\$$(echo $s | tr a-z A-Z); samtools view -@$T -T chr2.fa -b -o sr_$s.bam "$u" $WIN; samtools index sr_$s.bam; done
rm -f ./*.crai
for s in child father mother; do python3 "$HERE/sr_junction.py" sr_$s.bam chr2 48057041 L chr2 46970725 R - 1000 10 20 | sed "s/^{/{\"sample\": \"$s\", /"; done > "$OUT/sr_junction.jsonl"
# ---- 3. HiFi junction reads, trio SNVs, parent of origin ---------------------------------------------------------
samtools view -@$T -b -o hifi.bam $HIFI $LWIN && samtools index hifi.bam; rm -f ./*.bai.tmp
python3 "$HERE/hifi_junction.py" hifi.bam chr2 48057041 L chr2 46970725 R 50 20 1000 junction_reads.txt > "$OUT/hifi_junction.json"
python3 "$HERE/junction_from_reads.py" hifi.bam "$HERE/data/hg01943_junction.bedpe" "$OUT/hifi_junction_consensus.tsv" "$OUT/hg01943_junction_31mers.readderived.tsv" 31
bcftools mpileup --threads $T -f chr2.fa -r $WIN -q 20 -Q 20 -a AD,DP -Ou sr_child.bam sr_father.bam sr_mother.bam 2>/dev/null | bcftools call --threads $T -mv -V indels -Ob -o trio.bcf
bcftools query -i 'QUAL>=30' -f '%POS\t%REF\t%ALT[\t%GT]\n' trio.bcf > trio.tsv
# ---- 4. the assembly: HPRC alignments, region sequence, pieces ---------------------------------------------------
for h in hap1 hap2; do samtools view -b -o asm_$h.bam $WM/HG01943.GRCh38_no_alt.$h.sorted.bam $WIN; samtools index asm_$h.bam; samtools view -h asm_$h.bam | paftools.js sam2paf - > asm_$h.paf; done
rm -f ./*.bai.tmp; cat asm_hap1.paf asm_hap2.paf | cut -f1-12 > "$OUT/asm_pieces_vs_grch38.paf"
python3 "$HERE/trio_phase.py" trio.tsv sr_child.bam hifi.bam junction_reads.txt asm_hap1.bam,asm_hap2.bam chr2 46970725 48057041 > "$OUT/trio_phase.json"
samtools faidx $ASM/HG01943_pat_hprc_r2_v1.0.1.fa.gz 'HG01943#1#CM089354.1:46450001-49420000' | sed '1s/.*/>hap1/' > asm_local.fa
samtools faidx $ASM/HG01943_mat_hprc_r2_v1.0.1.fa.gz 'HG01943#2#CM089362.1:46420001-48680000' | sed '1s/.*/>hap2/' >> asm_local.fa
samtools faidx $ASM/HG01943_mat_hprc_r2_v1.0.1.fa.gz 'HG01943#2#JBHIIL010000075.1' | sed '1s/.*/>U/' >> asm_local.fa
samtools faidx asm_local.fa
# ---- 5. reads given back to the assembly (Flagger/NucFlag-style signals) -----------------------------------------
samtools fastq -F 0x900 hifi.bam 2>/dev/null > hifi.fq
python3 "$HERE/read_parent.py" trio.tsv hifi.bam chr2 read_parent.tsv
minimap2 -t $T -ax map-hifi --eqx -Y asm_local.fa hifi.fq 2>/dev/null | samtools sort -@1 -o hifi_on_asm.bam - && samtools index hifi_on_asm.bam
python3 "$HERE/asm_read_qc.py" hifi_on_asm.bam asm_local.fa read_parent.tsv "$OUT/asmqc" 5000 500 > /dev/null
# ---- 6. k-mers: copy-number track, inherited hapmers per piece, junction k-mers ----------------------------------
for s in child father mother; do samtools fastq -F 0xD00 sr_$s.bam 2>/dev/null > $s.fq; meryl count k=21 threads=$T memory=3 $s.fq output $s.meryl >/dev/null 2>&1; done
meryl count k=21 threads=$T memory=2 asm_local.fa output asm.meryl >/dev/null 2>&1
for s in mother father child; do meryl greater-than 4 $s.meryl output $s.solid.meryl >/dev/null 2>&1; done
meryl difference mother.solid.meryl father.meryl output mat_only.meryl >/dev/null 2>&1; meryl difference father.solid.meryl mother.meryl output pat_only.meryl >/dev/null 2>&1
meryl intersect mat_only.meryl child.solid.meryl output mat.inh.meryl >/dev/null 2>&1; meryl intersect pat_only.meryl child.solid.meryl output pat.inh.meryl >/dev/null 2>&1
samtools faidx chr2.fa $WIN | sed '1s/.*/>win/' > win.fa
for s in child father mother asm; do meryl-lookup -wig-count -sequence win.fa -output win.$s.wig -mers $s.meryl >/dev/null 2>&1; done
python3 "$HERE/kmer_cn_track.py" --offset 46400000 --haploid-depth 15 --out "$OUT/kmer_cn_bins.tsv" win.child.wig win.father.wig win.mother.wig win.asm.wig
python3 - <<'PY'
import re
s={};n=None
for l in open('asm_local.fa'):
    if l[0]=='>': n=l[1:].strip(); s[n]=[]
    else: s[n].append(l.strip())
s={k:''.join(v) for k,v in s.items()}
with open('pieces.fa','w') as f:
    for k,v in s.items():
        for i,p in enumerate(re.split('N{10,}',v)): f.write(f'>{k}_piece{i}\n{p}\n')
PY
meryl-lookup -existence -sequence pieces.fa -output "$OUT/hapmers_per_piece.tsv" -mers mat.inh.meryl pat.inh.meryl >/dev/null 2>&1
# ---- 7. junction k-mers in the whole assembly (optional) ---------------------------------------------------------
if [ "$WHOLE" = 1 ]; then python3 "$HERE/asm_kmer_search.py" "$OUT/hg01943_junction_31mers.readderived.tsv" "$OUT/junction_kmers_whole_assembly.json" $ASM/HG01943_pat_hprc_r2_v1.0.1.fa.gz $ASM/HG01943_mat_hprc_r2_v1.0.1.fa.gz > "$OUT/junction_kmers_whole_assembly.log"; fi
# ---- 8. local reassembly and regression scores -------------------------------------------------------------------
yak count -k31 -b37 -t$T -o pat.yak father.fq 2>/dev/null; yak count -k31 -b37 -t$T -o mat.yak mother.fq 2>/dev/null
hifiasm -o trio -t $T -1 pat.yak -2 mat.yak hifi.fq > hifiasm.log 2>&1
for h in hap1 hap2; do awk '$1=="S"{print ">"$2"\n"$3}' trio.dip.$h.p_ctg.gfa > local_$h.fa; done
samtools faidx asm_local.fa hap1 > hprc_pat.fa; samtools faidx asm_local.fa hap2 U > hprc_mat.fa
for x in "hprc_r2 hprc_pat.fa hprc_mat.fa" "local_hifiasm_trio local_hap1.fa local_hap2.fa"; do set -- $x
  python3 "$HERE/score_flagged_cases.py" asm-hg01943 --pat $2 --mat $3 --junction-kmers "$OUT/hg01943_junction_31mers.readderived.tsv" --ref-window win.fa \
     --hapmers-mat mat.inh.meryl --hapmers-pat pat.inh.meryl --threads $T --out "$OUT/score_R1_$1.json" > /dev/null; done
for h in hap1 hap2; do curl -sS -o fl_hifi_cons_$h.bed "$FL/v1.2.0_hifi/HG01943.hifi_minimap2_2.28.corrected.downsample_1.0.hmm_flagger_v1.2.0_prediction.conservative.canonical.no_Hap.$h.bed"
  curl -sS -o fl_hifi_full_$h.bed "$FL/v1.2.0_hifi/HG01943.hifi_minimap2_2.28.corrected.downsample_1.0.hmm_flagger_v1.2.0_prediction.canonical.no_Hap.$h.bed"; done
for h in pat mat; do curl -sS -o nucflag_$h.bed "$NF/HG01943_${h}_hprc_r2_v1.0.1.nucflag.bed"; done
python3 "$HERE/score_flagged_cases.py" qc-hg01943 --bed fl_hifi_cons_hap1.bed,fl_hifi_cons_hap2.bed --out "$OUT/score_R2_flagger_hifi_conservative.json" > /dev/null
python3 "$HERE/score_flagged_cases.py" qc-hg01943 --bed fl_hifi_full_hap1.bed,fl_hifi_full_hap2.bed --out "$OUT/score_R2_flagger_hifi_full.json" > /dev/null
python3 "$HERE/score_flagged_cases.py" qc-hg01943 --bed nucflag_pat.bed,nucflag_mat.bed --out "$OUT/score_R2_nucflag_v0.3.3_hifi.json" > /dev/null
for h in hap1 hap2; do curl -sS -o fl_ont_cons_$h.bed "$FL/v1.2.0_ont/HG01943.R941_minimap2_2.28.corrected.downsample_1.0.hmm_flagger_v1.2.0_prediction.conservative.canonical.no_Hap.$h.bed"
  curl -sS -o fl_ont_full_$h.bed "$FL/v1.2.0_ont/HG01943.R941_minimap2_2.28.corrected.downsample_1.0.hmm_flagger_v1.2.0_prediction.canonical.no_Hap.$h.bed"; done
python3 "$HERE/score_flagged_cases.py" qc-hg01943 --bed fl_ont_cons_hap1.bed,fl_ont_cons_hap2.bed --out "$OUT/score_R2_flagger_ont_conservative.json" > /dev/null
python3 "$HERE/score_flagged_cases.py" qc-hg01943 --bed fl_ont_full_hap1.bed,fl_ont_full_hap2.bed --out "$OUT/score_R2_flagger_ont_full.json" > /dev/null
echo "done: $OUT"
