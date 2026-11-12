#!/usr/bin/env bash
# Fetch sequences from a remote bgzipped, faidx-indexed assembly (HPRC release 2 *.fa.gz + .fai + .gzi)
# Usage: 02_fetch_truth_seq.sh <FA_GZ_URL> <OUT.fa> <SEQNAME[:start-end]> [more regions...]
set -euo pipefail
url=$1; out=$2; shift 2
tmp=$(mktemp -d "${TMPDIR:-$PWD}/scratch.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
cd "$tmp"   # htslib caches the remote .fai/.gzi in the working directory
samtools faidx "$url" "$@" > seq.fa
cd - >/dev/null; mv "$tmp/seq.fa" "$out"; samtools faidx "$out"; cut -f1,2 "$out.fai"
