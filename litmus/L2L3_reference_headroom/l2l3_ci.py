#!/usr/bin/env python3
"""L2/L3 pilot: block-bootstrap confidence intervals for the litmus quantities.

  l2l3_ci.py --manifest manifest.tsv --strata strata/ --sites sites/ --outdir results/ [--B 1000]

Blocks are 1 Mbp windows of GRCh38 chr20 (reads: GRCh38 primary placement; reads unmapped on
GRCh38 form one extra block; sites: their GRCh38 position). The same block resampling weights
are applied to every reference within a technology, so differences and 'fraction of the
oracle gap recovered' are paired.
Outputs: L2L3_readmetric_ci.tsv, L2L3_allele_balance_ci.tsv
  gap_recovered = (GRCh38 - X) / (GRCh38 - oracle)   (1 = realistic method reaches the oracle)
"""
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2l3_summarize import load, stratum_labels
from l2l3_strata import Intervals

BLK = 1_000_000


def boot(num, den, W):
    """num/den: per-block sums (refs x blocks); W: (B x blocks) multinomial weights"""
    return (W @ num.T) / np.maximum(W @ den.T, 1e-12)   # B x refs


def q(x):
    return np.nanpercentile(x, 2.5, axis=0), np.nanpercentile(x, 97.5, axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--sites', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--B', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=20261111)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    man = pd.read_csv(a.manifest, sep='\t')
    man = man[man.kind != 'external']
    sd = Intervals(os.path.join(a.strata, 'GRCh38.SD.bed.gz'))
    tr = Intervals(os.path.join(a.strata, 'GRCh38.TR.bed.gz'))
    cen = Intervals(os.path.join(a.strata, 'GRCh38.CEN.bed.gz'))
    out_r, out_a = [], []
    for tech, g in man.groupby('tech', sort=False):
        lr = tech != 'ilmn'
        refs = list(g.ref)
        data = {x.ref: load(x.npz) + (x.target.split(','),) for x in g.itertuples()}
        r0, n0, _ = data['GRCh38']
        keep = np.ones(len(r0), bool)
        if lr:
            keep[:] = False
            for k in ('GRCh38', 'CHM13', 'oracle'):
                r, names, tg = data[k]
                idx = [i for i, nme in enumerate(names) if nme in set(tg)]
                keep |= (r['mapped'] == 1) & np.isin(r['tid'], idx)
        lab = stratum_labels(r0, n0, sd, tr, lr, cen)
        # block id from GRCh38 placement: chr20 1-Mbp windows; other contigs one block each; unmapped one block
        tid20 = n0.index('chr20')
        blk = np.where(r0['mapped'] == 1,
                       np.where(r0['tid'] == tid20, r0['pos'] // BLK, 1000 + r0['tid']), 999)
        ublk, binv = np.unique(blk, return_inverse=True)
        nb = len(ublk)
        W = rng.multinomial(nb, np.full(nb, 1.0 / nb), size=a.B).astype(float)
        for s in ['all', 'all_noCEN'] + sorted(set(lab[keep])):
            m = keep if s == 'all' else (keep & (lab != 'CEN') if s == 'all_noCEN' else keep & (lab == s))
            if m.sum() == 0:
                continue
            den = np.bincount(binv[m], minlength=nb).astype(float)
            mets = {}
            for k in refs:
                r = data[k][0]
                bad = (r['mapped'][m] == 0) | (r['mapq'][m] < 20)
                mets.setdefault('frac_mapq_lt20', []).append(np.bincount(binv[m], weights=bad, minlength=nb))
                cl = np.maximum(r['clip_l'][m], r['clip_r'][m]) > (1000 if lr else 20)
                mets.setdefault('frac_clipped', []).append(np.bincount(binv[m], weights=cl & (r['mapped'][m] == 1), minlength=nb))
            for met, nums in mets.items():
                num = np.vstack(nums)
                est = num.sum(1) / den.sum()
                bs = boot(num, np.tile(den, (len(refs), 1)), W)
                lo, hi = q(bs)
                gi, oi = refs.index('GRCh38'), refs.index('oracle')
                for j, k in enumerate(refs):
                    row = dict(tech=tech, stratum=s, metric=met, ref=k, n_reads=int(m.sum()), n_blocks=int((den > 0).sum()),
                               value=est[j], ci_lo=lo[j], ci_hi=hi[j])
                    if k not in ('GRCh38', 'oracle'):
                        gap = bs[:, gi] - bs[:, oi]
                        rec = np.where(np.abs(gap) > 1e-12, (bs[:, gi] - bs[:, j]) / gap, np.nan)
                        row['gap_recovered'] = (est[gi] - est[j]) / (est[gi] - est[oi]) if est[gi] != est[oi] else np.nan
                        row['gap_recovered_ci_lo'], row['gap_recovered_ci_hi'] = q(rec[:, None])
                        row['gap_recovered_ci_lo'] = float(row['gap_recovered_ci_lo'][0])
                        row['gap_recovered_ci_hi'] = float(row['gap_recovered_ci_hi'][0])
                    out_r.append(row)
        # allele balance at truth het SNVs (MAPQ>=20)
        ann = pd.read_csv(os.path.join(a.sites, 'site_annot.tsv'), sep='\t').set_index('site_id')
        pos = ann.index.str.split(':').str[1].astype(int).values
        sblk = pos // BLK
        usb, sinv = np.unique(sblk, return_inverse=True)
        Ws = rng.multinomial(len(usb), np.full(len(usb), 1.0 / len(usb)), size=a.B).astype(float)
        cnt = {}
        for x in g.itertuples():
            tabs = [pd.read_csv(p, sep='\t').set_index('site_id') for p in x.pileups.split(',')]
            c = tabs[0][['ref_q20', 'alt_q20']].copy()
            for t in tabs[1:]:
                c = c.add(t[['ref_q20', 'alt_q20']], fill_value=0)
            cnt[x.ref] = c.reindex(ann.index).fillna(0)
        clean = (ann['clean'].values == 1) if 'clean' in ann.columns else np.ones(len(ann), bool)
        for siteset, sm, s in [(ss, smask, st_) for ss, smask in (('all', np.ones(len(ann), bool)), ('clean', clean))
                               for st_ in ['ALL', 'SD', 'TR', 'other', 'CEN']]:
            m = ((ann.stratum != 'CEN').values if s == 'ALL' else (ann.stratum == s).values) & sm
            if m.sum() == 0:
                continue
            num = np.vstack([np.bincount(sinv[m], weights=cnt[k].ref_q20.values[m], minlength=len(usb)) for k in refs])
            den = np.vstack([np.bincount(sinv[m], weights=(cnt[k].ref_q20.values + cnt[k].alt_q20.values)[m], minlength=len(usb)) for k in refs])
            est = num.sum(1) / np.maximum(den.sum(1), 1)
            bs = (Ws @ num.T) / np.maximum(Ws @ den.T, 1e-12)
            lo, hi = q(bs)
            gi, oi = refs.index('GRCh38'), refs.index('oracle')
            for j, k in enumerate(refs):
                dif = bs[:, j] - bs[:, gi]
                exc = np.abs(bs[:, gi] - 0.5) - np.abs(bs[:, j] - 0.5)
                row = dict(tech=tech, siteset=siteset, stratum=s, ref=k, n_sites=int(m.sum()), n_blocks=int(len(np.unique(sinv[m]))),
                           informative_reads=int(den[j].sum()), ref_frac=est[j], ci_lo=lo[j], ci_hi=hi[j],
                           diff_vs_grch38=est[j] - est[gi], diff_ci_lo=np.percentile(dif, 2.5), diff_ci_hi=np.percentile(dif, 97.5),
                           bias_reduction_vs_grch38=abs(est[gi] - 0.5) - abs(est[j] - 0.5),
                           bias_reduction_ci_lo=np.percentile(exc, 2.5), bias_reduction_ci_hi=np.percentile(exc, 97.5))
                out_a.append(row)
    pd.DataFrame(out_r).to_csv(os.path.join(a.outdir, 'L2L3_readmetric_ci.tsv'), sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(out_a).to_csv(os.path.join(a.outdir, 'L2L3_allele_balance_ci.tsv'), sep='\t', index=False, float_format='%.6g')


if __name__ == '__main__':
    main()
