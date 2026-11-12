#!/usr/bin/env bash
# L5 MAPQ-calibration pilot - end-to-end driver.
# usage: l5_run_all.sh OUTDIR [THREADS=2]
# Requires (conda/bioconda): minimap2=2.31 bwa=0.7.19 samtools pysam scikit-learn>=1.4 numpy pandas scipy matplotlib
# All large inputs/intermediates (~70 GB measured: 51 GB indexes, 10 GB inputs, 7 GB SAM) live in a mktemp
# scratch dir that is removed on exit; only small result tables and the figure are copied to OUTDIR.
# Peak RSS: minimap2 single-part human index builds 13.0-15.7 GB; mapping 9.6-14.3 GB; everything else <6 GB.
set -euo pipefail
OUT=$(mkdir -p "$1" && cd "$1" && pwd); T=${2:-2}
S=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T VECLIB_MAXIMUM_THREADS=$T NUMBA_NUM_THREADS=$T
W=$(mktemp -d "${TMPDIR:-/tmp}/l5.XXXXXX"); trap 'rm -rf "$W"' EXIT
cd "$W"; mkdir -p data idx sim aln feat feat_strict recal recal_strict diag models
# 1. inputs (md5-checked where the provider publishes md5s)
bash "$S/l5_fetch_data.sh" data
curl -sSfL -o data/chm13v2.0_SD.full.bed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.full.bed
samtools faidx data/chm13v2.0_maskedY_rCRS.fa.gz
gzip -dc data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > data/grch38.fna && samtools faidx data/grch38.fna
# 2. bwa indexes: NCBI prebuilt (GRCh38), built locally (CHM13, ~50 min, 1 thread)
tar -xzf data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz -C idx
bwa index -a bwtsw -p idx/chm13v2.0_maskedY_rCRS data/chm13v2.0_maskedY_rCRS.fa.gz 2> idx/bwa_index_chm13.log
for e in amb ann bwt pac sa; do
  ln -sf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.$e idx/grch38.$e; ln -sf chm13v2.0_maskedY_rCRS.$e idx/chm13.$e
done
# 3. simulation with lifted truth
C="--hap data/NA19240_pat_hprc_r2_v1.0.1.fa.gz --hap data/NA19240_mat_hprc_r2_v1.0.1.fa.gz --paf grch38=data/NA19240.vs_grch38.paf.gz --paf chm13=data/NA19240.vs_chm13.paf.gz --sd grch38=data/hg38.genomicSuperDups.txt.gz --sd chm13=data/chm13v2.0_SD.full.bed"
python3 "$S/l5_simulate.py" --tech hifi $C --n 40000 --pool 300000 --frac-e 0.4 --seed 202 --out sim/hifi
python3 "$S/l5_simulate.py" --tech ont $C --n 24000 --pool 300000 --frac-e 0.4 --seed 303 --out sim/ont
python3 "$S/l5_simulate.py" --tech illumina $C --n 600000 --seed 101 --out sim/ill
# 4. minimap2 index builds + mapping + feature extraction (one memory-heavy step at a time)
bash "$S/l5_pipeline.sh" . "$T" "chm13 grch38"
# 5. recalibration (default truth) + strict colinear truth for long reads
python3 "$S/l5_recal.py" --feat feat --outdir recal --save-models models
for run in grch38.mm2_hifi grch38.mm2_lrhq grch38.mm2_ont chm13.mm2_hifi chm13.mm2_lrhq chm13.mm2_ont; do
  ref=${run%%.*}; case ${run#*.} in mm2_hifi) t=hifi;; *) t=ont;; esac
  case $ref in grch38) sd=data/hg38.genomicSuperDups.txt.gz;; chm13) sd=data/chm13v2.0_SD.full.bed;; esac
  python3 "$S/l5_extract.py" --bam aln/$run.sam --truth sim/$t.truth.tsv.gz --ref $ref --sd $sd --run $run --min-cov-frac 0.9 --out feat_strict/$run.feat.tsv.gz
done
python3 "$S/l5_recal.py" --feat feat_strict --outdir recal_strict
# 6. diagnosis of wrong MAPQ-60 placements
for run in grch38.bwa_sr grch38.mm2_sr grch38.mm2_hifi grch38.mm2_lrhq grch38.mm2_ont chm13.bwa_sr chm13.mm2_sr chm13.mm2_hifi chm13.mm2_lrhq chm13.mm2_ont; do
  ref=${run%%.*}
  case ${run#*.} in bwa_sr|mm2_sr) r=sim/ill_R1.fq.gz,sim/ill_R2.fq.gz; p=sr;; mm2_hifi) r=sim/hifi.fa.gz; p=map-hifi;; mm2_lrhq) r=sim/ont.fa.gz; p=lr:hq;; mm2_ont) r=sim/ont.fa.gz; p=map-ont;; esac
  case $ref in grch38) fa=data/grch38.fna;; chm13) fa=data/chm13v2.0_maskedY_rCRS.fa.gz;; esac
  python3 "$S/l5_q60_diag.py" --feat feat/$run.feat.tsv.gz --reads $r --ref $fa --preset $p --out diag/$run.q60diag.tsv
done
python3 - <<'PY'
import glob, numpy as np, pandas as pd
rows, allr = [], []
for f in sorted(glob.glob('diag/*.q60diag.tsv')):
    d = pd.read_csv(f, sep='\t'); run = f.split('/')[-1].replace('.q60diag.tsv', ''); d['run'] = run; allr.append(d)
    gap = (d.AS_mapped - d.AS_truth)[d.cls == 'mapped_better']
    rows.append(dict(run=run, n_wrong_mapq60=len(d), **{k: int((d.cls == k).sum()) for k in ('mapped_better', 'tie', 'truth_better', 'truth_unaligned', 'mapped_unaligned')},
                     frac_mapped_better=round(float((d.cls == 'mapped_better').mean()), 3), frac_in_sd=round(float(d.in_sd.mean()), 3),
                     frac_same_chrom=round(float(d.same_chrom.mean()), 3), median_AS_gap_mapped_better=float(gap.median()) if len(gap) else np.nan,
                     median_truth_sd_identity=round(float(d.t_sdid[d.in_sd == 1].median()), 4)))
pd.DataFrame(rows).to_csv('recal/L5_q60_error_diagnosis_summary.tsv', sep='\t', index=False)
pd.concat(allr).to_csv('recal/L5_q60_error_diagnosis_reads.tsv.gz', sep='\t', index=False, compression='gzip')
PY
# 7. tables, figure, small outputs
python3 "$S/l5_tables.py" recal recal_strict > "$OUT/L5_tables.md"
python3 "$S/l5_figure.py" recal "$OUT/L5_mapq_calibration.png"
cp recal/L5_calibration.tsv recal/L5_recalibrator_metrics.tsv recal/L5_q60_error_diagnosis_summary.tsv recal/L5_q60_error_diagnosis_reads.tsv.gz recal/L5_run_meta.json "$OUT/"
cp recal_strict/L5_recalibrator_metrics.tsv "$OUT/L5_strict_truth_recalibrator_metrics.tsv"
cp sim/*.log aln/rss.log data/md5_report.tsv "$OUT/" 2>/dev/null || true
echo "L5 done -> $OUT"
