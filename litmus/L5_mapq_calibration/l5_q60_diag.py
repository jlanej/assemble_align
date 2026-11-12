#!/usr/bin/env python3
"""l5_q60_diag.py - for wrongly placed reads with MAPQ >= Q (default 60), realign the read to
(i) its lifted truth locus and (ii) its reported locus (each +/- FLANK bp) with the same minimap2 preset
and compare alignment scores (AS). Classes:
  truth_better   AS(truth) >  AS(mapped) + tol  -> aligner placed the read sub-optimally (search/heuristic error)
  tie            |AS difference| <= tol          -> the two loci are indistinguishable for this read
  mapped_better  AS(mapped) > AS(truth) + tol    -> read genuinely resembles the other paralog
                                                    (haplotype/paralog divergence; truth-definition limit)
usage: l5_q60_diag.py --feat F.feat.tsv.gz --reads READS.fa.gz|R1.fq.gz[,R2.fq.gz] --ref REF.fa(.gz)
                      --preset map-hifi --out OUT.tsv [--minq 60 --flank 2000 --tol-frac 0.001]
"""
import argparse, gzip, os, subprocess, sys, tempfile
import numpy as np
import pandas as pd
import pysam


def read_seqs(paths, names):
    out = {}
    for mate, p in enumerate(paths.split(','), 1 if ',' in paths else 0):
        with pysam.FastxFile(p) as fh:
            for r in fh:
                if (r.name, mate) in names:
                    out[(r.name, mate)] = r.sequence
    return out


def best_as(paf_lines):
    best = {}
    for l in paf_lines:
        f = l.rstrip('\n').split('\t')
        a = [int(x[5:]) for x in f[12:] if x.startswith('AS:i:')]
        if a:
            best[f[5]] = max(best.get(f[5], -10 ** 9), a[0])
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--feat', required=True); ap.add_argument('--reads', required=True)
    ap.add_argument('--ref', required=True); ap.add_argument('--preset', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--minq', type=int, default=60)
    ap.add_argument('--flank', type=int, default=2000); ap.add_argument('--tol-frac', type=float, default=0.001)
    a = ap.parse_args()
    d = pd.read_csv(a.feat, sep='\t', low_memory=False)
    w = d[(d.truth_ok == 1) & (d.mapped == 1) & (d.correct == 0) & (d.mapq >= a.minq)].copy()
    seqs = read_seqs(a.reads, set(zip(w.qname, w.mate)))
    fa = pysam.FastaFile(a.ref)
    rows = []
    with tempfile.TemporaryDirectory(dir='.') as td:
        for r in w.itertuples():
            s = seqs.get((r.qname, r.mate))
            if s is None:
                continue
            wins = {'truth': (r.t_chrom, int(r.t_start), int(r.t_end)), 'mapped': (r.chrom, int(r.pos), int(r.end))}
            tf = os.path.join(td, 't.fa'); qf = os.path.join(td, 'q.fa')
            with open(tf, 'w') as fo:
                for k, (c, s0, e0) in wins.items():
                    L = fa.get_reference_length(c)
                    fo.write('>%s\n%s\n' % (k, fa.fetch(c, max(0, s0 - a.flank), min(L, e0 + a.flank))))
            with open(qf, 'w') as fo:
                fo.write('>q\n%s\n' % s)
            res = subprocess.run(['minimap2', '-t', '1', '-c', '-x', a.preset, '--secondary=yes', '-N', '10', '-p', '0.1', tf, qf],
                                 capture_output=True, text=True)
            b = best_as(res.stdout.splitlines())
            at, am = b.get('truth', np.nan), b.get('mapped', np.nan)
            tol = max(2, a.tol_frac * len(s))
            cls = ('truth_unaligned' if np.isnan(at) else 'mapped_unaligned' if np.isnan(am) else
                   'truth_better' if at > am + tol else 'mapped_better' if am > at + tol else 'tie')
            rows.append(dict(qname=r.qname, mate=r.mate, mapq=r.mapq, in_sd=r.in_sd, t_sdid=r.t_sdid, ambig=r.ambig,
                             same_chrom=int(r.chrom == r.t_chrom), dist=(abs(r.pos - r.t_start) if r.chrom == r.t_chrom else -1),
                             AS_truth=at, AS_mapped=am, AS_reported=r.AS, cls=cls))
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep='\t', index=False)
    if len(out):
        print(out.cls.value_counts().to_dict(), 'n', len(out))


if __name__ == '__main__':
    main()
