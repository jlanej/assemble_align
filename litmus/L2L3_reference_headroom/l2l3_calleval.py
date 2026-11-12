#!/usr/bin/env python3
"""L2/L3 pilot: compare a call set used to personalize the reference against the
assembly-based truth (chr20 SNVs inside the dipcall-like confident region).

  l2l3_calleval.py --truth-prefix truth.chr20 --vcf calls.vcf.gz --name pers2p --out eval.tsv

Reports het/hom-ALT SNV recall (allele match; genotype match), precision, and for phased
call sets the agreement of phase with the trio-phased assembly (hap1=paternal) and the
switch rate between consecutive heterozygous SNVs.
"""
import argparse, gzip, sys
import numpy as np
import pandas as pd
import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth-prefix', required=True)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--name', required=True)
    ap.add_argument('--chrom', default='chr20')
    ap.add_argument('--out', required=True)
    ap.add_argument('--exclude-bed', help='BED(.gz) of regions removed from the confident region (e.g. centromere models)')
    a = ap.parse_args()
    het = pd.read_csv(a.truth_prefix + '.het_snv.tsv.gz', sep='\t')
    hom = pd.read_csv(a.truth_prefix + '.homalt_snv.tsv.gz', sep='\t')
    conf = pd.read_csv(a.truth_prefix + '.confident.bed.gz', sep='\t', header=None, names=['c', 's', 'e'])
    if a.exclude_bed:
        ex = pd.read_csv(a.exclude_bed, sep='\t', header=None, names=['c', 's', 'e'])
        ex = ex[ex.c == a.chrom]
        keep = []
        for r in conf.itertuples():
            segs = [(r.s, r.e)]
            for x in ex.itertuples():
                segs = [y for (s0, e0) in segs for y in ((s0, min(e0, x.s)), (max(s0, x.e), e0)) if y[1] > y[0]]
            keep.extend(segs)
        conf = pd.DataFrame(keep, columns=['s', 'e'])
    st, en = conf.s.values, conf.e.values

    def inconf(p):
        i = np.searchsorted(st, p, side='right') - 1
        return (i >= 0) & (p < en[np.maximum(i, 0)])

    truth = {}
    for r in het.itertuples():
        truth[r.pos0] = (r.alt, 'het', r.alt_hap)
    for r in hom.itertuples():
        truth[r.pos0] = (r.alt, 'hom', 0)
    calls = {}
    v = pysam.VariantFile(a.vcf)
    smp = list(v.header.samples)[0]
    phased_any = False
    for rec in v.fetch(a.chrom):
        if len(rec.ref) != 1 or rec.alts is None:
            continue
        s = rec.samples[smp]
        gt = s['GT']
        if gt is None or None in gt:
            continue
        for ai, alt in enumerate(rec.alts, start=1):
            if len(alt) != 1 or alt.startswith('<'):
                continue
            n = sum(1 for g in gt if g == ai)
            if n == 0:
                continue
            ph = s.phased
            phased_any |= ph
            first = 1 if (ph and gt[0] == ai) else (2 if (ph and gt[1] == ai) else 0)
            calls[rec.pos - 1] = (alt.upper(), 'hom' if n == 2 else 'het', first if n == 1 else 0)
    cpos = np.array(sorted(calls))
    cin = cpos[inconf(cpos)] if len(cpos) else cpos
    tpos = np.array(sorted(truth))
    tin = tpos[inconf(tpos)]
    rows = {}
    for cls in ('het', 'hom'):
        tp_allele = tp_gt = n = 0
        for p in tin:
            alt, c, _ = truth[p]
            if c != cls:
                continue
            n += 1
            if p in calls and calls[p][0] == alt:
                tp_allele += 1
                if calls[p][1] == cls:
                    tp_gt += 1
        rows['truth_%s_snv' % cls] = n
        rows['recall_allele_%s' % cls] = tp_allele / max(n, 1)
        rows['recall_genotype_%s' % cls] = tp_gt / max(n, 1)
    fp = sum(1 for p in cin if not (p in truth and truth[p][0] == calls[p][0]))
    rows['calls_snv_in_confident'] = int(len(cin))
    rows['confident_bp'] = int((en - st).sum())
    rows['precision_allele'] = 1 - fp / max(len(cin), 1)
    if phased_any:
        same = []
        for p in tin:
            alt, c, ah = truth[p]
            if c == 'het' and p in calls and calls[p][1] == 'het' and calls[p][2] > 0 and calls[p][0] == alt:
                same.append(calls[p][2] == ah)
        same = np.array(same)
        rows['phased_het_compared'] = int(len(same))
        rows['phase_agree_frac'] = float(same.mean()) if len(same) else np.nan
        rows['switch_rate_adjacent'] = float(np.mean(same[1:] != same[:-1])) if len(same) > 1 else np.nan
    df = pd.DataFrame([dict(callset=a.name, **rows)])
    df.to_csv(a.out, sep='\t', index=False, float_format='%.6g')
    print(df.T.to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()
