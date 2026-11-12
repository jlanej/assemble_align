#!/usr/bin/env bash
# L2/L3 pilot: call-set evaluation, long-read residual-error profiles, metric tables, figures.
# usage: 08_summarize.sh WORKDIR SCRIPTDIR OUTDIR
set -euo pipefail
S=$1; SD=$2; OUT=$3; M=$S/map
mkdir -p $OUT
# 1) personalization call sets vs assembly truth (chr20 confident region minus centromere model +-1 Mbp)
python $SD/l2l3_calleval.py --truth-prefix $S/truth/truth.chr20 --vcf $S/pers/NA19240.panel.chr20.seq.vcf.gz \
   --name 1kGP_phased_panel --exclude-bed $S/strata/GRCh38.CEN.bed.gz --out $S/eval/calleval.1KG.tsv
python $SD/l2l3_calleval.py --truth-prefix $S/truth/truth.chr20 --vcf $S/pers/NA19240.twopass.chr20.vcf.gz \
   --name twopass_bcftools_ilmn --exclude-bed $S/strata/GRCh38.CEN.bed.gz --out $S/eval/calleval.twopass.tsv
python - $S/eval/calleval.1KG.tsv $S/eval/calleval.twopass.tsv $OUT/L2L3_callsets_vs_truth.tsv <<'EOF'
import sys, pandas as pd
pd.concat([pd.read_csv(f, sep='\t') for f in sys.argv[1:-1]]).to_csv(sys.argv[-1], sep='\t', index=False)
EOF
cp $S/truth/truth.chr20.summary.json $OUT/L2L3_truth_summary.json
# 2) long-read residual errors vs own assembly (oracle best-of-two alignments on the chr20 sequences)
for tech in hifi ont; do
  python $SD/l2l3_lrerr.py --pair $M/$tech.oracle.c1.sorted.bam $S/asm/hap1.fa "NA19240#1#CM099588.1" \
                           --pair $M/$tech.oracle.c2.sorted.bam $S/asm/hap2.fa "NA19240#2#CM099626.1" \
                           --minmapq 20 --out $S/eval/lrerr.$tech
  cp $S/eval/lrerr.$tech.context.tsv $OUT/L2L3_lr_error_context_$tech.tsv
  cp $S/eval/lrerr.$tech.qv_calibration.tsv $OUT/L2L3_lr_qv_calibration_$tech.tsv
done
# 3) NYGC production alignment (bwa-mem 0.7.15, GRCh38 + alt/decoy/HLA) of the same Illumina reads, as context
[ -s $M/ilmn.NYGC.npz ] || python $SD/l2l3_reads.py --bam $S/ilmn/chr20.pairs.bam --out $M/ilmn.NYGC.npz
# 4) 'clean' het-SNV annotation (no indel-representation conflicts with the personalization call sets)
python $SD/l2l3_sites.py --truth $S/truth/truth.chr20.het_snv.tsv.gz --strata $S/strata --outdir $S/sites --annot-only \
   --indel-vcf $S/pers/NA19240.panel.chr20.seq.vcf.gz --indel-vcf $S/pers/NA19240.twopass.chr20.vcf.gz \
   --chain pers1KG_1=$S/pers/pers1KG_1.chain --chain pers1KG_2=$S/pers/pers1KG_2.chain \
   --chain pers2p_A=$S/pers/pers2p_A.chain --chain pers2p_R=$S/pers/pers2p_R.chain
# 5) manifest + metric tables
H1=NA19240#1#CM099588.1; H2=NA19240#2#CM099626.1
{
  printf 'tech\tref\tkind\tnpz\ttarget\tpileups\n'
  for t in ilmn hifi ont; do
    printf '%s\tGRCh38\thaploid\t%s\tchr20\t%s\n' $t $M/$t.GRCh38.npz $M/$t.GRCh38.ab.tsv.gz
    printf '%s\tCHM13\thaploid\t%s\tchr20\t%s\n' $t $M/$t.CHM13.npz $M/$t.CHM13.ab.tsv.gz
    printf '%s\toracle\tdiploid\t%s\th1_%s,h2_%s\t%s,%s\n' $t $M/$t.oracle.npz $H1 $H2 $M/$t.oracle.ab1.tsv.gz $M/$t.oracle.ab2.tsv.gz
    for p in pers1KG pers2p; do
      printf '%s\t%s\tdiploid\t%s\th1_chr20,h2_chr20\t%s,%s\n' $t $p $M/$t.$p.npz $M/$t.$p.ab1.tsv.gz $M/$t.$p.ab2.tsv.gz
    done
  done
  printf 'ilmn\tNYGC_bwamem_GRCh38alt\texternal\t%s\tchr20\t\n' $M/ilmn.NYGC.npz
} > $S/eval/manifest.tsv
python $SD/l2l3_summarize.py --manifest $S/eval/manifest.tsv --strata $S/strata --sites $S/sites --outdir $OUT
python $SD/l2l3_recall_strata.py --truth $S/truth/truth.chr20.het_snv.tsv.gz --annot $S/sites/site_annot.tsv \
   --vcf 1kGP_phased_panel=$S/pers/NA19240.panel.chr20.seq.vcf.gz --vcf twopass_bcftools_ilmn=$S/pers/NA19240.twopass.chr20.vcf.gz \
   --out $OUT/L2L3_callset_recall_by_stratum.tsv
python $SD/l2l3_ci.py --manifest $S/eval/manifest.tsv --strata $S/strata --sites $S/sites --outdir $OUT --B 1000
# 6) figures
python $SD/l2l3_figures.py --results $OUT --outdir $OUT
