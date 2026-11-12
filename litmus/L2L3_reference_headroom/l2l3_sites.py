#!/usr/bin/env python3
"""L2/L3 pilot: per-reference site tables for allele-balance pileups.

  l2l3_sites.py --truth truth.chr20.het_snv.tsv.gz --strata strata/ \
     --chain pers1KG_1=pers1KG_1.chain --chain pers1KG_2=pers1KG_2.chain ... --outdir sites/

Writes sites/<REF>.sites.tsv (site_id contig pos0 ref_x alt_x) for GRCh38, CHM13, hap1,
hap2 and every chain-lifted personalized haplotype, plus sites/site_annot.tsv
(site_id stratum near_indel alt_hap) where stratum is CEN (centromere model +-1 Mbp) /
SD / TR / other by GRCh38 position (in that order of precedence).
"""
import argparse, gzip, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2l3_strata import Intervals

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def read_chain(path, chrom):
    """return list of (t_start, q_start, size) ungapped blocks for chrom (+ strand only)"""
    blocks = []
    cur = None
    with open(path) as fh:
        for l in fh:
            f = l.split()
            if not f:
                continue
            if f[0] == 'chain':
                cur = (f[2] == chrom and f[7] == chrom)
                t, q = int(f[5]), int(f[10])
                if f[4] != '+' or f[9] != '+':
                    raise SystemExit('unexpected strand in chain')
                continue
            if not cur:
                continue
            size = int(f[0])
            blocks.append((t, q, size))
            t += size
            q += size
            if len(f) == 3:
                t += int(f[1])
                q += int(f[2])
    return blocks


def lift(blocks, pos):
    ts = np.array([b[0] for b in blocks])
    qs = np.array([b[1] for b in blocks])
    sz = np.array([b[2] for b in blocks])
    i = np.searchsorted(ts, pos, side='right') - 1
    ok = (i >= 0)
    ii = np.where(ok, i, 0)
    inside = ok & (pos < ts[ii] + sz[ii])
    return np.where(inside, qs[ii] + (pos - ts[ii]), -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True)
    ap.add_argument('--strata', required=True)
    ap.add_argument('--chain', action='append', default=[])
    ap.add_argument('--chrom', default='chr20')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--annot-only', action='store_true')
    ap.add_argument('--indel-vcf', action='append', default=[],
                    help='call sets used for personalization; truth SNVs within --indel-pad bp of any of their non-SNV records are not "clean"')
    ap.add_argument('--indel-pad', type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    t = pd.read_csv(a.truth, sep='\t')
    t['site_id'] = t.chrom + ':' + t.pos0.astype(str)
    sd = Intervals(os.path.join(a.strata, 'GRCh38.SD.bed.gz'))
    tr = Intervals(os.path.join(a.strata, 'GRCh38.TR.bed.gz'))
    p = t.pos0.values
    insd = sd.overlap(a.chrom, p, p + 1) > 0
    intr = tr.overlap(a.chrom, p, p + 1) > 0
    t['stratum'] = np.where(insd, 'SD', np.where(intr, 'TR', 'other'))
    cenf = os.path.join(a.strata, 'GRCh38.CEN.bed.gz')
    if os.path.exists(cenf):
        incen = Intervals(cenf).overlap(a.chrom, p, p + 1) > 0
        t.loc[incen, 'stratum'] = 'CEN'
    # 'clean' sites: no truth indel within 5 bp (near_indel), no non-SNV record of any personalization
    # call set within indel-pad bp, and liftable to every personalized haplotype (chains given)
    bad = np.zeros(len(t), bool)
    import pysam
    for vf in a.indel_vcf:
        iv = []
        for rec in pysam.VariantFile(vf).fetch(a.chrom):
            if len(rec.ref) != 1 or any(len(x) != 1 for x in (rec.alts or ())):
                iv.append((rec.pos - 1 - a.indel_pad, rec.pos - 1 + len(rec.ref) + a.indel_pad))
        if iv:
            iv = np.array(sorted(iv))
            st = iv[:, 0]
            en = np.maximum.accumulate(iv[:, 1])
            i = np.searchsorted(st, p, side='right') - 1
            bad |= (i >= 0) & (p < en[np.maximum(i, 0)])
    for c in a.chain:
        name, path = c.split('=')
        bad |= lift(read_chain(path, a.chrom), p) < 0
    t['clean'] = ((t.near_indel == 0) & ~bad).astype(int)
    if a.annot_only:
        t[['site_id', 'stratum', 'near_indel', 'alt_hap', 'clean']].to_csv(os.path.join(a.outdir, 'site_annot.tsv'), sep='\t', index=False)
        return

    def comp(b):
        return COMP.get(b, 'N')

    def write(name, contig, pos, orient, mask=None):
        df = pd.DataFrame(dict(site_id=t.site_id, contig=contig, pos0=pos,
                               ref_x=[r if o >= 0 else comp(r) for r, o in zip(t.ref, orient)],
                               alt_x=[x if o >= 0 else comp(x) for x, o in zip(t.alt, orient)]))
        if mask is not None:
            df = df[mask]
        df.to_csv(os.path.join(a.outdir, name + '.sites.tsv'), sep='\t', index=False)
        print(name, len(df), file=sys.stderr)

    one = np.ones(len(t), int)
    write('GRCh38', t.chrom, t.pos0, one)
    write('CHM13', t.chm13_chrom, t.chm13_pos0, t.chm13_orient.values, mask=(t.chm13_pos0 >= 0).values)
    write('hap1', t.h1_contig, t.h1_pos0, t.h1_orient.values)
    write('hap2', t.h2_contig, t.h2_pos0, t.h2_orient.values)
    for c in a.chain:
        name, path = c.split('=')
        q = lift(read_chain(path, a.chrom), p)
        write(name, t.chrom, q, one, mask=(q >= 0))
    # CHM13 allele (from the projection: CHM13 base checked later in pileup ref column); annotate
    t[['site_id', 'stratum', 'near_indel', 'alt_hap', 'clean']].to_csv(os.path.join(a.outdir, 'site_annot.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
