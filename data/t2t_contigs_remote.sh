#!/usr/bin/env bash
# Count telomere-to-telomere (T2T) contigs in a bgzip-compressed assembly hosted over HTTPS
# using htslib random access only (downloads .fai + .gzi + the requested end windows; no full download).
# Usage: t2t_contigs_remote.sh <assembly.fa.gz URL> <out.tsv> [min_contig_bp=1000000] [end_window_bp=2000] [min_telo_frac=0.5]
# A contig end is "telomeric" if >= min_telo_frac of the end window is covered by (TTAGGG)n or (CCCTAA)n
# hexamers (vertebrate telomere motif; strand-agnostic). T2T contig = both ends telomeric.
# Requires samtools >= 1.10 built with libcurl; curl; awk.
set -euo pipefail
URL=${1:?assembly .fa.gz URL (bgzip, with .fai and .gzi alongside)}; OUT=${2:?out.tsv}
MINLEN=${3:-1000000}; WIN=${4:-2000}; MINF=${5:-0.5}
d=$(mktemp -d "${TMPDIR:-/tmp}/t2t.XXXXXX"); trap 'rm -rf "$d"' EXIT
curl -sfL "$URL.fai" -o "$d/asm.fai"; curl -sfL "$URL.gzi" -o "$d/asm.gzi"
awk -v m="$MINLEN" -v w="$WIN" '$2>=m{e=(w<$2?w:$2); s=$2-w+1; if(s<1)s=1; print $1":1-"e; print $1":"s"-"$2}' "$d/asm.fai" > "$d/regions.txt"
samtools faidx --fai-idx "$d/asm.fai" --gzi-idx "$d/asm.gzi" -r "$d/regions.txt" "$URL" > "$d/ends.fa"
awk -v W="$WIN" -v M="$MINF" '
  function frac(s,   a,b,t1,t2){ s=toupper(s); t1=s; a=gsub(/TTAGGG/,"",t1); t2=s; b=gsub(/CCCTAA/,"",t2); return ((a>b?a:b)*6)/length(s) }
  /^>/{ if(seq!=""){ F[name]=frac(seq) } name=substr($0,2); seq=""; order[++n]=name; next }
  { seq=seq $0 }
  END{ if(seq!="") F[name]=frac(seq);
       print "contig\tlength\tstart_telo_frac\tend_telo_frac\tt2t";
       for(i=1;i<=n;i+=2){ r1=order[i]; r2=order[i+1]; c=r1; sub(/:[0-9]+-[0-9]+$/,"",c);
         split(r2,p,/[:-]/); L=p[length(p)];
         t=(F[r1]>=M && F[r2]>=M)?"yes":"no";
         printf "%s\t%s\t%.3f\t%.3f\t%s\n", c, L, F[r1], F[r2], t } }' "$d/ends.fa" > "$OUT"
awk -F'\t' 'NR>1{n++; if($5=="yes") t++} END{printf "contigs_ge_minlen=%d\tT2T_contigs=%d\n", n, t}' "$OUT"
