#!/usr/bin/env python3
"""L2/L3 pilot: summarize per-read metrics and allele balance across references.

  l2l3_summarize.py --manifest manifest.tsv --strata strata/ --sites sites/ --outdir results/

manifest.tsv columns: tech ref kind npz target pileups
  tech   : ilmn | hifi | ont
  ref    : GRCh38 | CHM13 | oracle | pers1KG | pers2p ...  (GRCh38, CHM13, oracle required)
  kind   : haploid | diploid | external (independent alignment of the same reads, any order)
  npz    : output of l2l3_reads.py
  target : comma-separated contig names that represent chr20 in this reference
  pileups: comma-separated l2l3_pileup.py outputs (one per haplotype for diploid refs)
Read sets: ilmn = all templates of the NYGC chr20 slice; hifi/ont = prefilter candidates
whose primary alignment lies on chr20 in GRCh38, CHM13 or the oracle (union).
Strata (per read, from its GRCh38 primary alignment): CEN (>=50% of span in a centromere
model +-1 Mbp), SD (>=50% of aligned span in UCSC segmental duplications), TR (any overlap with TRF/homopolymer>=7 intervals; short
reads only), other, unplaced (unmapped on GRCh38). Allele-balance 'ALL' excludes CEN sites.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2l3_strata import Intervals


def load(npz):
    z = np.load(npz, allow_pickle=False)
    return z['reads'], [str(x) for x in z['refnames']]


def stratum_labels(r, names, sd, tr, long_reads, cen=None):
    lab = np.full(len(r), 'unplaced', dtype=object)
    mapped = r['mapped'] == 1
    tid = r['tid']
    for t in np.unique(tid[mapped]):
        c = names[t]
        m = mapped & (tid == t)
        s, e = r['pos'][m], r['end'][m]
        span = np.maximum(e - s, 1)
        osd = sd.overlap(c, s, e)
        lab_m = np.where(osd >= 0.5 * span, 'SD', 'other').astype(object)
        if not long_reads:
            otr = tr.overlap(c, s, e)
            lab_m = np.where((lab_m == 'other') & (otr > 0), 'TR', lab_m)
        else:
            lab_m = np.where(lab_m == 'other', 'nonSD', lab_m)
        if cen is not None:
            ocen = cen.overlap(c, s, e)
            lab_m = np.where(ocen >= 0.5 * span, 'CEN', lab_m)
        lab[m] = lab_m
    return lab


def metrics(r, names, target, long_reads):
    n = len(r)
    if n == 0:
        return {}
    mapped = r['mapped'] == 1
    mq = r['mapq']
    clipthr = 1000 if long_reads else 20
    clip = np.maximum(r['clip_l'], r['clip_r'])
    cols = r['cols'].astype(np.int64)
    nm = r['nm'].astype(np.int64)
    tset = set(target)
    on_t = np.array([names[t] in tset if t >= 0 else False for t in r['tid']]) if n < 2_000_000 else \
        np.isin(r['tid'], [i for i, x in enumerate(names) if x in tset])
    er = np.where(cols > 0, nm / np.maximum(cols, 1), np.nan)
    d = dict(n_reads=n,
             frac_unmapped=float((~mapped).mean()),
             frac_mapq0=float((mapped & (mq == 0)).mean()),
             frac_mapq_lt20=float((~mapped | (mq < 20)).mean()),
             frac_mapq_ge20=float((mapped & (mq >= 20)).mean()),
             frac_clipped=float((mapped & (clip > clipthr)).sum() / max(mapped.sum(), 1)),
             frac_with_supp=float((mapped & (r['nsupp'] > 0)).sum() / max(mapped.sum(), 1)),
             edit_rate_pooled=float(nm[mapped].sum() / max(cols[mapped].sum(), 1)),
             edit_rate_mean=float(np.nanmean(er[mapped])) if mapped.any() else np.nan,
             edit_rate_median=float(np.nanmedian(er[mapped])) if mapped.any() else np.nan,
             frac_on_chr20=float(on_t.mean()))
    if 'big_gap_bases' in r.dtype.names:
        bg = r['big_gap_bases'].astype(np.int64)
        d['edit_rate_small_pooled'] = float((nm[mapped] - bg[mapped]).sum() / max((cols[mapped] - bg[mapped]).sum(), 1))
        d['frac_with_indel_ge50'] = float((mapped & (bg > 0)).sum() / max(mapped.sum(), 1))
    if long_reads:
        de = r['de'][mapped]
        d['gc_divergence_mean'] = float(np.nanmean(de)) if mapped.any() else np.nan
        d['gc_divergence_median'] = float(np.nanmedian(de)) if mapped.any() else np.nan
        q = np.maximum(r['qlen'], 1)
        d['aligned_frac_mean'] = float(np.mean(np.minimum((r['qaln'] + r['supp_qaln']) / q, 1.0)))
        d['frac_lt90pct_aligned'] = float(np.mean(((r['qaln'] + r['supp_qaln']) / q) < 0.9))
    else:
        d['frac_proper_pair'] = float((r['proper'] == 1).mean())
        d['frac_perfect'] = float((mapped & (nm == 0) & (clip == 0)).mean())
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--sites', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    man = pd.read_csv(a.manifest, sep='\t')
    sd = Intervals(os.path.join(a.strata, 'GRCh38.SD.bed.gz'))
    tr = Intervals(os.path.join(a.strata, 'GRCh38.TR.bed.gz'))
    cen = Intervals(os.path.join(a.strata, 'GRCh38.CEN.bed.gz'))
    rows, srows, trows = [], [], []
    for tech, g in man.groupby('tech', sort=False):
        lr = tech != 'ilmn'
        data = {}
        external = []
        for x in g.itertuples():
            r, names = load(x.npz)
            if x.kind == 'external':
                external.append((x.ref, r, names, x.target.split(',')))
                continue
            data[x.ref] = (r, names, x.target.split(','))
        for k, r, names, tg in external:
            # independent alignment of the same templates (e.g. the original NYGC bwa-mem CRAM);
            # read order differs, so strata come from its own GRCh38 placement
            labx = stratum_labels(r, names, sd, tr, lr, cen)
            rows.append(dict(tech=tech, ref=k, readset='nygc_chr20_slice(own placement)', **metrics(r, names, tg, lr)))
            for s_ in np.unique(labx):
                srows.append(dict(tech=tech, ref=k, stratum=s_, **metrics(r[labx == s_], names, tg, lr)))
        base_r, base_n, _ = data['GRCh38']
        for k, (r, _, _) in data.items():
            if len(r) != len(base_r) or not np.array_equal(r['crc'], base_r['crc']):
                raise SystemExit('read order mismatch %s %s' % (tech, k))
        keep = np.ones(len(base_r), bool)
        if lr:
            keep = np.zeros(len(base_r), bool)
            for k in ('GRCh38', 'CHM13', 'oracle'):
                r, names, tg = data[k]
                idx = [i for i, nme in enumerate(names) if nme in set(tg)]
                keep |= (r['mapped'] == 1) & np.isin(r['tid'], idx)
        lab = stratum_labels(base_r, base_n, sd, tr, lr, cen)
        mq_g = base_r['mapq'] * (base_r['mapped'] == 1)
        ok_g = (base_r['mapped'] == 1) & (mq_g >= 20)
        for k, (r, names, tg) in data.items():
            d = metrics(r[keep], names, tg, lr)
            rows.append(dict(tech=tech, ref=k, readset='chr20_union' if lr else 'nygc_chr20_slice', **d))
            for s in np.unique(lab[keep]):
                m = keep & (lab == s)
                d = metrics(r[m], names, tg, lr)
                srows.append(dict(tech=tech, ref=k, stratum=s, **d))
            ok = (r['mapped'] == 1) & (r['mapq'] >= 20)
            trows.append(dict(tech=tech, ref=k,
                              grch38_lt20=int((keep & ~ok_g).sum()),
                              rescued_to_ge20=int((keep & ~ok_g & ok).sum()),
                              grch38_ge20=int((keep & ok_g).sum()),
                              lost_to_lt20=int((keep & ok_g & ~ok).sum())))
    pd.DataFrame(rows).to_csv(os.path.join(a.outdir, 'L2L3_metrics.tsv'), sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(srows).to_csv(os.path.join(a.outdir, 'L2L3_metrics_by_stratum.tsv'), sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(trows).to_csv(os.path.join(a.outdir, 'L2L3_mapq_transitions.tsv'), sep='\t', index=False)

    # allele balance
    ann = pd.read_csv(os.path.join(a.sites, 'site_annot.tsv'), sep='\t').set_index('site_id')
    arows = []
    for x in man.itertuples():
        if not isinstance(x.pileups, str) or not x.pileups:
            continue
        tabs = [pd.read_csv(p, sep='\t').set_index('site_id') for p in x.pileups.split(',')]
        num = [c for c in tabs[0].columns if c != 'xref_class']
        c = tabs[0][num].copy()
        for t in tabs[1:]:
            c = c.add(t[num], fill_value=0)
        xcls = tabs[0]['xref_class'] if len(tabs) == 1 else None
        c = c.reindex(ann.index).fillna(0)
        c['stratum'] = ann['stratum']
        c['near_indel'] = ann['near_indel']
        clean = ann['clean'].values == 1 if 'clean' in ann.columns else np.ones(len(ann), bool)
        for siteset, sm, s in [(ss, smask, st_) for ss, smask in (('all', np.ones(len(c), bool)), ('clean', clean))
                               for st_ in ['ALL', 'SD', 'TR', 'other', 'CEN']]:
            m = ((c['stratum'] != 'CEN').values if s == 'ALL' else (c['stratum'] == s).values) & sm
            cc = c[m]
            for q in ('q0', 'q20'):
                ra, aa = cc['ref_' + q].sum(), cc['alt_' + q].sum()
                dep = cc['ref_' + q] + cc['alt_' + q]
                hi = dep >= (10 if x.tech == 'ilmn' else 2)
                persite = (cc['ref_' + q][hi] / dep[hi]).mean() if hi.any() else np.nan
                d = dict(tech=x.tech, ref=x.ref, siteset=siteset, stratum=s, mapq_min=q[1:], n_sites=int(m.sum()),
                         n_sites_covered=int((dep > 0).sum()), ref_reads=int(ra), alt_reads=int(aa),
                         ref_frac_pooled=ra / max(ra + aa, 1), ref_frac_mean_sites=persite,
                         mean_depth_ref_alt=float(dep.mean()),
                         frac_sites_alt0_depth10=float(((cc['alt_' + q] == 0) & (dep >= 10)).sum() / max((dep >= 10).sum(), 1)))
                if x.ref in ('GRCh38', 'CHM13'):
                    xc = xcls.reindex(cc.index)
                    match = cc['ref_' + q].where(xc == 'R', 0).sum() + cc['alt_' + q].where(xc == 'A', 0).sum()
                    nonm = cc['alt_' + q].where(xc == 'R', 0).sum() + cc['ref_' + q].where(xc == 'A', 0).sum()
                    d['refgenome_allele_frac'] = match / max(match + nonm, 1)
                    d['frac_sites_refgenome_carries_alt'] = float((xc == 'A').mean())
                d['other_or_del_frac'] = float((cc['oth_' + q].sum() + cc['del_' + q].sum()) / max(ra + aa + cc['oth_' + q].sum() + cc['del_' + q].sum(), 1))
                arows.append(d)
    pd.DataFrame(arows).to_csv(os.path.join(a.outdir, 'L2L3_allele_balance.tsv'), sep='\t', index=False, float_format='%.6g')


if __name__ == '__main__':
    main()
