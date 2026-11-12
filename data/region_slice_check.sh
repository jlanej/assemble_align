#!/usr/bin/env bash
# Laptop smoke test: region-sliced HTTPS access to the main 1kGP-linked data types (no full downloads).
# Usage: region_slice_check.sh [region=chr20:40000000-40100000]
# Requires samtools/bcftools >= 1.10 with libcurl, curl. Prints mean depth / record counts and elapsed seconds.
set -euo pipefail
REG=${1:-chr20:40000000-40100000}
d=$(mktemp -d "${TMPDIR:-/tmp}/slice.XXXXXX"); trap 'rm -rf "$d"' EXIT; cd "$d"
REF=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
t() { local t0=$(date +%s); "$@"; echo "  elapsed $(( $(date +%s)-t0 )) s"; }
echo "[1] NYGC 30x CRAM NA19240 (remote reference)";      t bash -c "samtools depth -a -r $REG --reference $REF https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram | awk '{s+=\$3;n++} END{printf \"mean_depth=%.2f over %d bp\n\", s/n, n}'"
echo "[2] 1KGP PacBio HiFi HG01103 (S3 BAM)";             t bash -c "samtools depth -a -r $REG https://1000g-ont.s3.amazonaws.com/1KGP_PacBio_WGS/Alignment/GRCh38/HG01103-W-0.GRCh38_no_alt.longphase.haplotagged.bam | awk '{s+=\$3;n++} END{printf \"mean_depth=%.2f over %d bp\n\", s/n, n}'"
echo "[3] 1kG-ONT R10 HG01258 (S3 BAM; large .bai)";      t bash -c "samtools depth -a -r $REG https://1000g-ont.s3.amazonaws.com/PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/HG01258-ONT-hg38-R10-LSK114-5mCG_5hmCG.phased.bam | awk '{s+=\$3;n++} END{printf \"mean_depth=%.2f over %d bp\n\", s/n, n}'"
echo "[4] GIAB HG002 v5.0q smvar VCF records";            t bash -c "bcftools view -H -r $REG https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz | wc -l"
echo "[5] HPRC R2 NA19240 maternal assembly faidx";
A=https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA19240/assemblies/release2/NA19240_mat_hprc_r2_v1.0.1.fa.gz
curl -sfLO "$A.fai"; curl -sfLO "$A.gzi"; c=$(head -1 "$(basename $A).fai" | cut -f1)
t samtools faidx --fai-idx "$(basename $A).fai" --gzi-idx "$(basename $A).gzi" "$A" "${c}:1000001-1000060"
