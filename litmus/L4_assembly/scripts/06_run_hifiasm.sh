#!/usr/bin/env bash
# Run hifiasm (tested with 0.25.0-r726) in one of three modes and export haplotype contigs as FASTA.
# -f0 disables hifiasm's initial 2^37-bit (16 GB) Bloom filter, as recommended by the hifiasm docs for small
# genomes; without it a single-chromosome run peaked at 18.7 GB RSS on this machine (budget 12 GB).
# Usage: 06_run_hifiasm.sh <hifi|hifi_ul|ont> <OUT_DIR> <THREADS> <READS.fq.gz> [ONT_UL.fq.gz] [extra hifiasm args...]
# Outputs in OUT_DIR: hap1.fa hap2.fa primary.fa hifiasm.log cmd.txt time.txt (hifiasm's own real/CPU/peak-RSS line)
set -euo pipefail
mode=$1; out=$2; t=$3; reads=$4; ul=${5:-}; shift $(( $# < 5 ? $# : 5 )); extra=("$@")
export OMP_NUM_THREADS=$t
mkdir -p "$out"; out=$(cd "$out" && pwd)
tmp=$(mktemp -d "$(cd "${TMPDIR:-$PWD}" && pwd)/scratch.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
base=(hifiasm -o "$tmp/asm" -t "$t" -f0)
case "$mode" in
  hifi)    cmd=("${base[@]}" ${extra[@]+"${extra[@]}"} "$reads") ;;
  hifi_ul) cmd=("${base[@]}" --ul "$ul" ${extra[@]+"${extra[@]}"} "$reads") ;;
  ont)     cmd=("${base[@]}" --ont ${extra[@]+"${extra[@]}"} "$reads") ;;
  *) echo "unknown mode $mode" >&2; exit 1 ;;
esac
echo "${cmd[*]}" > "$out/cmd.txt"; hifiasm --version >> "$out/cmd.txt"
rc=0; "${cmd[@]}" 2> "$out/hifiasm.log" || rc=$?
grep -E "Real time" "$out/hifiasm.log" | tail -1 > "$out/time.txt" || true
ls "$tmp" > "$out/outputs.txt"
gfa2fa() { awk '/^S/{print ">"$2; print $3}' "$1"; }
h1=$(ls "$tmp"/asm*.hap1.p_ctg.gfa 2>/dev/null | head -1 || true); h2=$(ls "$tmp"/asm*.hap2.p_ctg.gfa 2>/dev/null | head -1 || true)
p=$(ls "$tmp"/asm*.p_ctg.gfa 2>/dev/null | grep -v hap | head -1 || true)
[[ -n "$h1" ]] && gfa2fa "$h1" > "$out/hap1.fa"; [[ -n "$h2" ]] && gfa2fa "$h2" > "$out/hap2.fa"; [[ -n "$p" ]] && gfa2fa "$p" > "$out/primary.fa"
for f in "$tmp"/asm*.noseq.gfa; do [[ -e "$f" ]] && gzip -c "$f" > "$out/$(basename "$f").gz"; done
grep -E "peak_hom|peak_het" "$out/hifiasm.log" | head -2 >> "$out/cmd.txt" || true
echo "exit=$rc" >> "$out/cmd.txt"
[[ $rc -eq 0 && -s "$out/hap1.fa" ]] || { echo "hifiasm failed rc=$rc ($out)" >&2; tail -5 "$out/hifiasm.log" >&2; exit 1; }
