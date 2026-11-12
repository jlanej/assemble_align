#!/usr/bin/env python3
"""l5_extract.py - per-read alignment features + correctness labels for simulated reads.

Record selection and correctness follow `paftools.js mapeval` (minimap2 v2.31) mode 0:
  * secondary records (0x100) are ignored for selection; among the remaining mapped records the one
    with most matching bases (M/X/= length - (NM - gap bases)) is evaluated; MAPQ capped at 60;
  * correct  <=>  same chromosome AND same strand AND overlap/union > 0.1 with the lifted truth.
    For reads whose origin is covered by two assembly-alignment records of near-equal coverage
    (runner-up >= 0.9 x best; 'ambiguous truth') a match to the runner-up locus also counts.
Features use only SAM fields/tags and reference annotation (SD overlap of the mapped interval).

usage: l5_extract.py --bam ALN.bam --truth READS.truth.tsv.gz --ref grch38|chm13 --sd SD_FILE
                     --run NAME --out FEATURES.tsv.gz [--min-cov-frac 0.9]
"""
import argparse, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # PYTHONSAFEPATH-safe
import numpy as np
import pandas as pd
import pysam
from l5_common import SDIndex, read_sd_table, fold_of_chrom
from l5_features import new_state, update, rows_from, add_mate_and_sd, mate_of


def log(*a):
    sys.stderr.write(time.strftime('[%H:%M:%S] ') + ' '.join(map(str, a)) + '\n'); sys.stderr.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True); ap.add_argument('--truth', required=True)
    ap.add_argument('--ref', required=True); ap.add_argument('--sd', required=True)
    ap.add_argument('--run', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--min-cov-frac', type=float, default=None,
                    help='truth valid if lifted aligned bases >= frac x read length (default 0.9 short, 0.5 long)')
    a = ap.parse_args()
    sd = SDIndex(read_sd_table(a.sd, a.ref))
    states = {}   # (qname, mate) -> per-read state (l5_features)
    n_rec = 0
    with pysam.AlignmentFile(a.bam, 'rb' if a.bam.endswith('.bam') else 'r', check_sq=False) as fh:
        for r in fh:
            n_rec += 1
            key = (r.query_name, mate_of(r))
            if key not in states:
                states[key] = new_state()
            update(states[key], r)
    log(a.run, 'records', n_rec, 'reads', len(states))
    df = add_mate_and_sd(rows_from(states), sd)
    # truth
    t = pd.read_csv(a.truth, sep='\t')
    R = a.ref
    keep = ['qname', 'mate', 'hap', 'len', 'weight', 'stratumE', R + '_chrom', R + '_start', R + '_end', R + '_strand',
            R + '_cov', R + '_cov2', R + '_chrom2', R + '_start2', R + '_end2', R + '_strand2', R + '_sdfrac', R + '_sdid']
    t = t[[c for c in keep if c in t.columns]].rename(columns=lambda c: c.replace(R + '_', 't_'))
    df = df.merge(t, on=['qname', 'mate'], how='right', validate='one_to_one')
    df['mapped'] = df['mapped'].fillna(0).astype(int)
    frac = a.min_cov_frac if a.min_cov_frac is not None else (0.9 if (df.mate > 0).any() else 0.5)
    df['truth_ok'] = (df.t_cov >= frac * df.len).astype(int)
    df['ambig'] = ((df.t_cov2 >= 0.9 * df.t_cov) & (df.t_cov > 0)).astype(int)

    def ovl(c1, s1, e1, st1, c2, s2, e2, st2):
        same = (c1 == c2) & (st1 == st2)
        o = np.minimum(e1, e2) - np.maximum(s1, s2)
        u = np.maximum(e1, e2) - np.minimum(s1, s2)
        return same & (o > 0) & (o / np.maximum(u, 1) > 0.1)
    c1 = ovl(df.chrom.values, df.pos.values, df.end.values, df.strand.values,
             df.t_chrom.values, df.t_start.values, df.t_end.values, df.t_strand.values)
    c2 = ovl(df.chrom.values, df.pos.values, df.end.values, df.strand.values,
             df.t_chrom2.values, df.t_start2.values, df.t_end2.values, df.t_strand2.values) & (df.ambig.values == 1)
    df['correct'] = ((c1 | c2) & (df.mapped.values == 1)).astype(int)
    df['correct_strict'] = (c1 & (df.mapped.values == 1)).astype(int)
    df['fold'] = [fold_of_chrom(c) if isinstance(c, str) and c else -1 for c in df.t_chrom.values]
    df['in_sd'] = (df.t_sdfrac >= 0.5).astype(int)
    df['run'] = a.run
    df.to_csv(a.out, sep='\t', index=False, compression='gzip')
    ev = df[(df.truth_ok == 1) & (df.mapped == 1)]
    log(a.run, 'reads', len(df), 'mapped', int(df.mapped.sum()), 'evaluable', len(ev),
        'wrong', int((ev.correct == 0).sum()), 'wrong@MAPQ60', int(((ev.correct == 0) & (ev.mapq == 60)).sum()))


if __name__ == '__main__':
    main()
