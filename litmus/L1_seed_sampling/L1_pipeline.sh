#!/usr/bin/env bash
# L1_pipeline.sh -- end-to-end reproduction of litmus pilot L1 (k-mer sampling schemes for minimap2).
# usage: L1_pipeline.sh WORK_DIR [--wholegenome]
#   WORK_DIR receives data/, results/, fq/ and the deliverables (L1_results.tsv, L1_density_sensitivity.png, ...).
# env:   python3 with numpy>=2, numba, pandas, matplotlib; samtools, minimap2 2.31, pbsim3 (pbsim), k8, cc.
#        Single-threaded by design (NUMBA/OMP threads = 1); peak RSS of the Python steps was 3-6 GB.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$1"; W=$(cd "$1" && pwd); WG=${2:-}
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export L1_GM_DIR="$W/data/gm" L1_SCRATCH="$W"
cd "$W"
[ -s data/chm13_chr20.fa ] || bash "$HERE/L1_fetch_data.sh" data 2 ${WG:+--keep-genome}
# 1. C harness around minimap2's own sketch.c (stock and patched) for bit-identity checks
scr=$(mktemp -d "$W/scratch.XXXXXX"); trap 'rm -rf "$scr"' EXIT
for v in stock modoc; do
  mkdir -p "$scr/$v" harness_$v
  tar xzf data/minimap2-2.31.tar.gz -C "$scr/$v" --strip-components 1 --exclude='*/.gitmodules' 2>/dev/null || true
  [ $v = modoc ] && ( cd "$scr/$v" && patch -p1 -s < "$HERE/minimap2_sampling.patch" )
  for f in sketch.c kalloc.c kalloc.h kvec.h mmpriv.h minimap.h bseq.h kseq.h; do cp "$scr/$v/$f" harness_$v/; done
  cp "$HERE/L1_sketch_harness.c" harness_$v/sketch_harness.c
  ( cd harness_$v && cc -O2 -I. -o sketch_harness sketch_harness.c sketch.c kalloc.c )
done
# 2. validation (port == sketch.c; python mod-OC == patched sketch.c; published densities; strand symmetry)
python3 "$HERE/L1_validate.py" data/chm13_chr20.fa L1_validation.tsv harness_stock/sketch_harness harness_modoc/sketch_harness
# 3. experiments (CHM13 chr20 + NA19240 reads/tiles; GRCh38 index stats; block counts; curves; hash fix)
for part in main grch38 blocks curves hashw; do python3 "$HERE/L1_run.py" data results $part; done
# 4. statistics, tables, figure
python3 "$HERE/L1_stats.py" results L1_results.tsv
python3 "$HERE/L1_figure.py" results L1_results.tsv L1_density_sensitivity.png
# 5. stretch: patched minimap2 vs stock on the chr20 index (mapeval), optional whole genome
python3 "$HERE/L1_make_fastq.py" data fq
bash "$HERE/L1_patch_eval.sh" data/chm13_chr20.fa fq "$HERE/minimap2_sampling.patch" results/patch 1
python3 "$HERE/L1_tables.py" L1_results.tsv results L1_tables.md
if [ "$WG" = "--wholegenome" ]; then
  B=$scr/bin; mkdir -p $B
  for v in stock modoc; do
    mkdir -p "$B/$v"; tar xzf data/minimap2-2.31.tar.gz -C "$B/$v" --strip-components 1 --exclude='*/.gitmodules' 2>/dev/null || true
    [ $v = modoc ] && ( cd "$B/$v" && patch -p1 -s < "$HERE/minimap2_sampling.patch" )
    if [ "$(uname -m)" = arm64 ] || [ "$(uname -m)" = aarch64 ]; then ( cd "$B/$v" && make -s arm_neon=1 aarch64=1 >/dev/null 2>&1 ); else ( cd "$B/$v" && make -s >/dev/null 2>&1 ); fi
  done
  bash "$HERE/L1_wholegenome.sh" data/chm13v2.0.fa.gz "$B/stock/minimap2" "$B/modoc/minimap2" fq/hifi.fq "$B/stock/misc/paftools.js" results/wholegenome 1
fi
echo "L1 pipeline finished in $W"
