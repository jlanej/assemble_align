#!/usr/bin/env bash
# Reproduce the region-slice byte fractions for the NA19240 NYGC 30x CRAM (1kGP related-samples collection).
# Usage: run_region_slice_bytes.sh [OUT_JSON] [CRAM_URL] [MD5]
# Only the .crai index, the first 1 MB of the CRAM (header) and two small region tables are downloaded.
set -euo pipefail
OUT=${1:-region_slice_bytes_NA19240.json}
CRAM_URL=${2:-https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram}
MD5=${3:-e6d4b9e1e19d1f4d4f5de6a0d16d9c10}
d=$(mktemp -d "$PWD/scratch.XXXXXX"); trap 'rm -rf "$d"' EXIT
python3 "$(dirname "$0")/region_slice_bytes.py" --cram-url "$CRAM_URL" --md5 "$MD5" \
  --ucsc-sd GRCh38_SD_merged=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
  --bed MRJD_19_regions=https://1000g-ont.s3.amazonaws.com/Illumina_TruPath_WGS/HG00097/MRJD_regions_38.bed.gz \
  --pad 0 --pad 1000 --scratch "$d" --out "$OUT"
