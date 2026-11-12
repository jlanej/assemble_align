#!/usr/bin/env bash
# L2/L3 pilot: fetch public inputs for NA19240 into WORKDIR (all public, no credentials).
# usage: 00_fetch.sh WORKDIR [lr_bases_per_file=3100000000]
# MD5s verified in the 2026-11-11 run are listed next to each file (see also ../PROVENANCE.tsv).
set -euo pipefail
S=$1; LRB=${2:-3100000000}
mkdir -p $S/ref $S/asm $S/vcf $S/ann $S/ilmn $S/reads
HPRC=https://s3-us-west-2.amazonaws.com/human-pangenomics
G1K=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp
G1KS3=https://1000genomes.s3.amazonaws.com
dl(){ curl -sS -L --retry 5 --retry-delay 10 -o "$1" "$2"; }
# references
dl $S/ref/GRCh38_noalt.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz   # md5 a08035b6a6e31780e96a34008ff21bd6
dl $S/ref/chm13v2.0_maskedY_rCRS.fa.gz $HPRC/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz  # md5 b17cd2bdd8eb5544e195d6a8882e417b
dl $S/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa $G1KS3/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa  # CRAM decoding; chr20 M5 b18e6c531b0bd70e949a7fc20859cb01
dl $S/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai $G1KS3/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
gzip -dc $S/ref/GRCh38_noalt.fna.gz > $S/ref/GRCh38.fa && samtools faidx $S/ref/GRCh38.fa
gzip -dc $S/ref/chm13v2.0_maskedY_rCRS.fa.gz > $S/ref/CHM13.fa && samtools faidx $S/ref/CHM13.fa
# HPRC release 2 NA19240 trio-phased hifiasm 0.19.9 assembly (hap1 = paternal, hap2 = maternal)
dl $S/asm/hap1.fa.gz $HPRC/working/HPRC_PLUS/NA19240/assemblies/release2/NA19240_pat_hprc_r2_v1.0.1.fa.gz   # md5 2cd406c57e30c7b5a911b276bb255b1e
dl $S/asm/hap2.fa.gz $HPRC/working/HPRC_PLUS/NA19240/assemblies/release2/NA19240_mat_hprc_r2_v1.0.1.fa.gz   # md5 b9ed3f11aeb037a0d83f23fb93f8d1fc
for h in hap1 hap2; do gzip -dc $S/asm/$h.fa.gz > $S/asm/$h.fa && samtools faidx $S/asm/$h.fa; done
# 1kGP high-coverage phased panel (SNV/INDEL/SV), chr20
dl $S/vcf/panel.chr20.vcf.gz $G1K/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
dl $S/vcf/panel.chr20.vcf.gz.tbi $G1K/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
# UCSC hg38 annotation tables
for t in genomicSuperDups simpleRepeat centromeres gap; do dl $S/ann/hg38.$t.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/$t.txt.gz; done
# NYGC 30x CRAM (ERR3989410; md5 e6d4b9e1e19d1f4d4f5de6a0d16d9c10), downloaded whole in 10 parallel ranges
U=$G1KS3/1000G_2504_high_coverage/additional_698_related/data/ERR3989410/NA19240.final.cram
SIZE=$(curl -sI $U | awk 'tolower($1)=="content-length:"{print $2+0}'); N=10; CH=$(( (SIZE + N - 1) / N ))
for i in $(seq 0 $((N-1))); do s=$((i*CH)); e=$(( (i+1)*CH - 1 )); [ $e -ge $SIZE ] && e=$((SIZE-1)); curl -sS --retry 8 -r $s-$e -o $S/ilmn/part.$i $U & done; wait
cat $(for i in $(seq 0 $((N-1))); do echo $S/ilmn/part.$i; done) > $S/ilmn/NA19240.final.cram && rm -f $S/ilmn/part.*
dl $S/ilmn/NA19240.final.cram.crai $U.crai
# long reads: first LRB bases of each HPRC release-2 file (DeepConsensus v1.2 q20 HiFi; R10.4.1 dorado sup4.1.0 ONT)
SD=$(cd "$(dirname "$0")" && pwd)
bash $SD/01_stream_long_reads.sh $S/reads $LRB
