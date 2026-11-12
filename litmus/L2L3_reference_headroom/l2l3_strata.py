#!/usr/bin/env python3
"""L2/L3 pilot: GRCh38 stratification intervals and interval-overlap helpers.

  l2l3_strata.py --fasta GRCh38.fa --sd genomicSuperDups.txt.gz --srep simpleRepeat.txt.gz \
                 --outdir strata/
Writes strata/GRCh38.SD.bed.gz (UCSC genomicSuperDups, merged) and
strata/GRCh38.TR.bed.gz (UCSC simpleRepeat (TRF) union homopolymers >=7 bp, each padded
by 5 bp, merged) and strata/GRCh38.CEN.bed.gz (UCSC centromere models, min..max per
chromosome, padded by 1 Mbp), restricted to chr1-22,X,Y. Labels are assigned elsewhere with
Intervals.overlap() (prefix-sum overlap of disjoint sorted intervals).
"""
import argparse, gzip, os
import numpy as np


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


class Intervals:
    """disjoint sorted intervals per chromosome with O(log n) overlap queries"""

    def __init__(self, bed_gz):
        d = {}
        op = gzip.open if bed_gz.endswith('.gz') else open
        with op(bed_gz, 'rt') as fh:
            for l in fh:
                c, s, e = l.split('\t')[:3]
                d.setdefault(c, []).append((int(s), int(e)))
        self.iv = {}
        for c, v in d.items():
            v = np.array(merge(v), dtype=np.int64)
            st, en = v[:, 0], v[:, 1]
            cum = np.concatenate([[0], np.cumsum(en - st)])
            self.iv[c] = (st, en, cum)

    def _F(self, c, x):
        st, en, cum = self.iv[c]
        i = np.searchsorted(st, x, side='right') - 1
        ok = i >= 0
        ii = np.where(ok, i, 0)
        v = cum[ii] + np.minimum(x, en[ii]) - st[ii]
        return np.where(ok, v, 0)

    def overlap(self, c, s, e):
        s = np.asarray(s, np.int64)
        e = np.asarray(e, np.int64)
        if c not in self.iv:
            return np.zeros(len(s), np.int64)
        return self._F(c, e) - self._F(c, s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--sd', required=True)
    ap.add_argument('--srep', required=True)
    ap.add_argument('--cen', help='UCSC centromeres.txt.gz (hg38); span min..max per chromosome, padded to flanking gaps by --cenpad')
    ap.add_argument('--cenpad', type=int, default=1000000)
    ap.add_argument('--hp', type=int, default=7)
    ap.add_argument('--pad', type=int, default=5)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    import pysam
    os.makedirs(a.outdir, exist_ok=True)
    chroms = ['chr%s' % x for x in list(range(1, 23)) + ['X', 'Y']]
    fa = pysam.FastaFile(a.fasta)
    sd = {}
    with gzip.open(a.sd, 'rt') as fh:
        for l in fh:
            f = l.split('\t')
            if f[1] in chroms:
                sd.setdefault(f[1], []).append((int(f[2]), int(f[3])))
    with gzip.open(os.path.join(a.outdir, 'GRCh38.SD.bed.gz'), 'wt') as fo:
        for c in chroms:
            for s, e in merge(sd.get(c, [])):
                fo.write('%s\t%d\t%d\n' % (c, s, e))
    tr = {}
    with gzip.open(a.srep, 'rt') as fh:
        for l in fh:
            f = l.split('\t')
            if f[1] in chroms:
                tr.setdefault(f[1], []).append((max(0, int(f[2]) - a.pad), int(f[3]) + a.pad))
    nhp = 0
    for c in chroms:
        s = np.frombuffer(fa.fetch(c).upper().encode(), dtype=np.uint8)
        n = len(s)
        ch = np.empty(n, bool)
        ch[0] = True
        ch[1:] = s[1:] != s[:-1]
        st = np.flatnonzero(ch)
        ln = np.diff(np.append(st, n))
        keep = (ln >= a.hp) & (s[st] != ord('N'))
        for b, l in zip(st[keep], ln[keep]):
            tr.setdefault(c, []).append((max(0, int(b) - a.pad), int(b + l) + a.pad))
        nhp += int(keep.sum())
    with gzip.open(os.path.join(a.outdir, 'GRCh38.TR.bed.gz'), 'wt') as fo:
        for c in chroms:
            for s, e in merge(tr.get(c, [])):
                fo.write('%s\t%d\t%d\n' % (c, s, e))
    print('homopolymers>=%d: %d' % (a.hp, nhp))
    if a.cen:
        cen = {}
        with gzip.open(a.cen, 'rt') as fh:
            for l in fh:
                f = l.split('\t')
                if f[1] in chroms:
                    s0, e0 = int(f[2]), int(f[3])
                    c0 = cen.get(f[1], (s0, e0))
                    cen[f[1]] = (min(c0[0], s0), max(c0[1], e0))
        with gzip.open(os.path.join(a.outdir, 'GRCh38.CEN.bed.gz'), 'wt') as fo:
            for c in chroms:
                if c in cen:
                    fo.write('%s\t%d\t%d\n' % (c, max(0, cen[c][0] - a.cenpad), cen[c][1] + a.cenpad))


if __name__ == '__main__':
    main()
