#!/usr/bin/env bash
# L2/L3 pilot: realistic personalized references for NA19240 chr20.
#  (a) pers1KG : 1kGP high-coverage phased panel (SNV/indel/sequence-resolved SV; symbolic
#                SV alleles dropped), haplotype 1 and 2 -> diploid consensus (bcftools consensus -H 1/-H 2)
#  (b) pers2p  : two-pass from the same Illumina reads: bcftools mpileup/call on the GRCh38
#                minimap2 alignment, hard filter, then an UNPHASED diploid pair:
#                pers2p_A = all ALT alleles (het + hom), pers2p_R = hom-ALT only.
# usage: 04_two_pass_personalize.sh MODE ...
#   04_two_pass_personalize.sh panel  GRCh38.fa PANEL_chr20.vcf.gz OUTDIR
#   04_two_pass_personalize.sh twopass GRCh38.fa ilmn.GRCh38.sorted.bam OUTDIR [threads]
set -euo pipefail
mode=$1; REF=$2; IN=$3; O=$4; T=${5:-4}
mkdir -p $O
if [ "$mode" = panel ]; then
  bcftools view -s NA19240 -c 1 $IN -Ou | bcftools view -e 'ALT~"<"' -Oz -o $O/NA19240.panel.chr20.seq.vcf.gz
  bcftools index -t $O/NA19240.panel.chr20.seq.vcf.gz
  for h in 1 2; do
    bcftools consensus -f $REF -H $h -c $O/pers1KG_$h.chain $O/NA19240.panel.chr20.seq.vcf.gz > $O/pers1KG_$h.fa 2> $O/pers1KG_$h.log
    samtools faidx $O/pers1KG_$h.fa
  done
else
  # chr20 split into T chunks, called in parallel (one thread each)
  L=$(samtools view -H $IN | awk -F'\t' '$2=="SN:chr20"{sub("LN:","",$3); print $3}')
  step=$(( (L + T - 1) / T )); parts=""
  for i in $(seq 0 $((T-1))); do s=$((i*step+1)); e=$(( (i+1)*step )); [ $e -gt $L ] && e=$L
    ( bcftools mpileup -Ou -f $REF -r chr20:$s-$e -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP --max-depth 1000 $IN 2>/dev/null | \
      bcftools call -mv -Oz -o $O/twopass.part$i.vcf.gz 2>/dev/null ) &
    parts="$parts $O/twopass.part$i.vcf.gz"
  done; wait
  bcftools concat $parts -Ou | bcftools norm -f $REF -m -any -Ou | \
    bcftools view -i 'QUAL>=30 && FORMAT/DP>=8 && FORMAT/DP<=90' -Oz -o $O/NA19240.twopass.chr20.vcf.gz
  bcftools index -t $O/NA19240.twopass.chr20.vcf.gz; rm -f $parts
  for h in A R; do
    bcftools consensus -f $REF -H $h -c $O/pers2p_$h.chain $O/NA19240.twopass.chr20.vcf.gz > $O/pers2p_$h.fa 2> $O/pers2p_$h.log
    samtools faidx $O/pers2p_$h.fa
  done
fi
