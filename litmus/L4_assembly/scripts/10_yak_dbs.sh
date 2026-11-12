#!/usr/bin/env bash
# Build yak k-mer databases (k=31) from chromosome-sliced Illumina reads of child and parents.
# Usage: 10_yak_dbs.sh <OUT_DIR> <THREADS> <child.fq.gz> <father.fq.gz> <mother.fq.gz>
set -euo pipefail
out=$1; t=$2; mkdir -p "$out"
yak count -k31 -b32 -t "$t" -o "$out/child.yak" "$3"
yak count -k31 -b32 -t "$t" -o "$out/pat.yak" "$4"
yak count -k31 -b32 -t "$t" -o "$out/mat.yak" "$5"
ls -la "$out"
