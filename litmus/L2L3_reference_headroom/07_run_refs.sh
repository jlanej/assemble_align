#!/usr/bin/env bash
# L2/L3 pilot: map + post-process every reference sequentially (one minimap2 at a time, -t 4).
# usage: 07_run_refs.sh WORKDIR SCRIPTDIR [threads] [steps...]
# WORKDIR layout: ref/{GRCh38,CHM13}.fa  asm/hap{1,2}.fa  pers/pers{1KG_1,1KG_2,2p_A,2p_R}.fa
#                 ilmn/R{1,2}.fq.gz  lr/{hifi,ont}.cand.fq.gz  sites/*.sites.tsv  map/ (output)
# steps (default all): GRCh38 CHM13 oracle pers1KG pers2p
set -euo pipefail
S=$1; SD=$2; T=${3:-4}; shift 3 || true
steps=${@:-GRCh38 CHM13 oracle pers1KG pers2p}
M=$S/map; mkdir -p $M
hap(){ # label fasta
  [ -s $M/ont.$1.bam ] || bash $SD/05_map_reference.sh $1 $2 $M $S/ilmn $S/lr $T 0
  for tech in ilmn hifi ont; do bash $SD/06_post.sh hap $tech $1 $M $2 $S/sites/$1.sites.tsv $T; done
}
dip(){ # label h1label h2label h1fasta h2fasta h1sites h2sites
  [ -s $M/ont.$2.bam ] || bash $SD/05_map_reference.sh $2 $4 $M $S/ilmn $S/lr $T 0
  [ -s $M/ont.$3.bam ] || bash $SD/05_map_reference.sh $3 $5 $M $S/ilmn $S/lr $T 0
  for tech in ilmn hifi ont; do bash $SD/06_post.sh dip $tech $1 $M $2 $3 $4 $5 $6 $7 $T; done
}
for st in $steps; do
  t0=$(date +%s)
  case $st in
    GRCh38)  hap GRCh38 $S/ref/GRCh38.fa ;;
    CHM13)   hap CHM13 $S/ref/CHM13.fa ;;
    oracle)  dip oracle hap1 hap2 $S/asm/hap1.fa $S/asm/hap2.fa $S/sites/hap1.sites.tsv $S/sites/hap2.sites.tsv ;;
    pers1KG) dip pers1KG pers1KG_1 pers1KG_2 $S/pers/pers1KG_1.fa $S/pers/pers1KG_2.fa $S/sites/pers1KG_1.sites.tsv $S/sites/pers1KG_2.sites.tsv ;;
    pers2p)  dip pers2p pers2p_A pers2p_R $S/pers/pers2p_A.fa $S/pers/pers2p_R.fa $S/sites/pers2p_A.sites.tsv $S/sites/pers2p_R.sites.tsv ;;
  esac
  echo "step $st done in $(( $(date +%s) - t0 )) s"
done
