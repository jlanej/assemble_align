#!/usr/bin/env python3
"""l5_common.py - shared helpers for the L5 MAPQ-calibration pilot.
SDIndex: merged segmental-duplication intervals per chromosome with the max pairwise identity
(fracMatch) of the SDs merged into each interval; vectorised overlap fraction / identity queries."""
import gzip
import numpy as np
import pandas as pd


def read_sd_table(path, ref):
    """GRCh38: UCSC hg38 genomicSuperDups.txt.gz (cols chrom=1,start=2,end=3,fracMatch=26).
    CHM13: T2T chm13v2.0_SD.bed (SEDEF; BED with fracMatch column, detected by header/name)."""
    if ref == 'grch38':
        df = pd.read_csv(path, sep='\t', header=None, usecols=[1, 2, 3, 26],
                         names=['chrom', 'start', 'end', 'ident'], compression='gzip')
    else:
        op = gzip.open if str(path).endswith('.gz') else open
        with op(path, 'rt') as fh:
            first = fh.readline()
        hdr = first.lstrip('#').rstrip('\n').split('\t')
        if 'fracMatch' not in hdr:
            raise ValueError('need the SEDEF table with a fracMatch column (chm13v2.0_SD.full.bed)')
        # header has duplicated names (start1/end1 twice) -> read positionally
        fm = hdr.index('fracMatch')
        cols = [0, 1, 2, fm]
        both = 'chr2' in hdr
        if both:
            c2 = hdr.index('chr2'); cols += [c2, c2 + 1, c2 + 2]
        raw = pd.read_csv(path, sep='\t', skiprows=1, header=None, usecols=cols)
        a = raw[[0, 1, 2, fm]].copy(); a.columns = ['chrom', 'start', 'end', 'ident']
        if both:
            b = raw[[c2, c2 + 1, c2 + 2, fm]].copy(); b.columns = ['chrom', 'start', 'end', 'ident']
            a = pd.concat([a, b], ignore_index=True)
        df = a
    df['ident'] = pd.to_numeric(df['ident'], errors='coerce').fillna(0.0)
    return df


class SDIndex:
    def __init__(self, df):
        self.idx = {}
        for c, g in df.groupby('chrom'):
            g = g.sort_values('start')
            s, e, idn = g.start.values.astype(np.int64), g.end.values.astype(np.int64), g.ident.values
            ms, me, mid = [], [], []
            cs, ce, ci = s[0], e[0], idn[0]
            for a, b, i in zip(s[1:], e[1:], idn[1:]):
                if a <= ce:
                    ce = max(ce, b); ci = max(ci, i)
                else:
                    ms.append(cs); me.append(ce); mid.append(ci); cs, ce, ci = a, b, i
            ms.append(cs); me.append(ce); mid.append(ci)
            ms, me = np.array(ms), np.array(me)
            self.idx[c] = (ms, me, np.r_[0, np.cumsum(me - ms)], np.array(mid))
        self.total_bp = sum(int(v[2][-1]) for v in self.idx.values())

    def _covered_before(self, c, x):
        ms, me, cum, _ = self.idx[c]
        i = np.searchsorted(ms, x, side='right') - 1
        out = np.zeros(len(x), np.int64)
        ok = i >= 0
        out[ok] = cum[i[ok]] + np.minimum(x[ok], me[i[ok]]) - ms[i[ok]]
        return out

    def query(self, chrom, a, b):
        """chrom: array of str; a,b: int arrays. returns (overlap_fraction, max_identity)."""
        chrom = np.asarray(chrom, dtype=object)
        a = np.asarray(a, np.int64); b = np.asarray(b, np.int64)
        frac = np.zeros(len(a)); ident = np.zeros(len(a))
        for c in pd.unique(chrom):
            if c not in self.idx:
                continue
            m = np.flatnonzero(chrom == c)
            aa, bb = a[m], np.maximum(b[m], a[m])
            ov = self._covered_before(c, bb) - self._covered_before(c, aa)
            frac[m] = ov / np.maximum(1, bb - aa)
            ms, me, _, mid = self.idx[c]
            i1 = np.searchsorted(me, aa, side='right')
            i2 = np.searchsorted(ms, bb, side='left')
            has = (i2 > i1) & (ov > 0)
            if has.any():
                # max identity over overlapped merged intervals (ranges are short; loop is cheap)
                vals = np.zeros(has.sum())
                for k, (x1, x2) in enumerate(zip(i1[has], i2[has])):
                    vals[k] = mid[x1:x2].max()
                ident[m[has]] = vals
        return frac, ident


def fold_of_chrom(chrom):
    """Chromosome hold-out fold: odd autosomes + chrX -> 1, even autosomes, chrY, chrM, unplaced -> 0.
    Unlocalised '_random' contigs inherit their parent chromosome."""
    c = str(chrom)
    base = c.split('_')[0].replace('chr', '')
    if base.isdigit():
        return int(base) % 2
    return 1 if base == 'X' else 0
