#!/usr/bin/env python3
"""l5_simulate.py - simulate reads with known origin from the NA19240 HPRC release-2 haplotypes and
lift every read's origin to GRCh38 and CHM13 through the impg/wfmash assembly-to-reference PAFs.

Custom vectorised simulator (pbsim3 3.0.x was installable, but its template mode did not deliver
the requested per-read accuracy: QSHMM-ONT-HQ at --accuracy-mean 0.99 gave mean 3.3% / max 15%
error on a 30-read test). Error models:
  illumina : 2x150 paired-end, fragment ~ N(450,60) clipped [200,1000]; uniform substitutions 0.2%
             (equivalent to `wgsim -e 0.002 -r 0 -R 0 -X 0 -d 450 -s 60`).
  hifi     : length ~ N(17000,2500) clipped [12000,24000]; per-read Phred Q ~ N(30,3) clipped [22,40];
             error events: 30% substitution, 50% homopolymer run-length +/-1 (60% shortening),
             10% random insertion, 10% random deletion.
  ont      : length ~ lognormal(median 30 kb, sigma 0.3) clipped [20000,50000]; per-read Q ~ N(20,2.5)
             clipped [12,28]; events: 35% substitution, 35% homopolymer +/-1 (60% shortening),
             10% random insertion, 20% random deletion. Homopolymer events hit runs >=3 bp with
             probability proportional to (run length - 2).
Long reads use SD-stratified sampling: a uniform candidate pool is lifted first; candidates whose truth
interval is >=20% SD on GRCh38 or CHM13 ('E' stratum) are kept with probability p_E, the rest with
p_N; every read carries weight = 1/p so that weighted statistics equal uniform-sampling statistics.

usage: l5_simulate.py --tech {illumina,hifi,ont} --hap HAP1.fa.gz --hap HAP2.fa.gz
         --paf grch38=NA19240.vs_grch38.paf.gz --paf chm13=NA19240.vs_chm13.paf.gz
         --sd grch38=hg38.genomicSuperDups.txt.gz --sd chm13=chm13v2.0_SD.bed
         --n N [--pool M --frac-e 0.4] --seed S --out PREFIX
"""
import argparse, gzip, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # PYTHONSAFEPATH-safe
import numpy as np
import pandas as pd
import pysam
from l5_liftover import load_paf, lift_intervals, strip_pansn
from l5_common import SDIndex, read_sd_table

ENC = np.full(256, 4, np.uint8)
for i, ch in enumerate('ACGT'):
    ENC[ord(ch)] = i; ENC[ord(ch.lower())] = i
DEC = np.frombuffer(b'ACGTN', np.uint8)

MODELS = {
    'hifi': dict(qmu=30.0, qsd=3.0, qlo=22, qhi=40, comp=[0.30, 0.50, 0.10, 0.10], hp_del=0.6),
    'ont': dict(qmu=20.0, qsd=2.5, qlo=12, qhi=28, comp=[0.35, 0.35, 0.10, 0.20], hp_del=0.6),
}


def log(*a):
    sys.stderr.write(time.strftime('[%H:%M:%S] ') + ' '.join(map(str, a)) + '\n'); sys.stderr.flush()


def read_fai(path):
    fai = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1], names=['name', 'len'])
    return list(fai.name), fai.len.values.astype(np.int64)


def sample_positions(lens, L, rng):
    """contig ~ length, start ~ U[0, len-L]; resample reads that do not fit."""
    n = len(L)
    ctg = np.zeros(n, np.int64); st = np.zeros(n, np.int64)
    todo = np.arange(n)
    p = lens / lens.sum()
    while len(todo):
        c = rng.choice(len(lens), size=len(todo), p=p)
        room = lens[c] - L[todo]
        ok = room >= 0
        ctg[todo[ok]] = c[ok]
        st[todo[ok]] = (rng.random(ok.sum()) * (room[ok] + 1)).astype(np.int64)
        todo = todo[~ok]
    return ctg, st


def revcomp(x):
    y = x[::-1].copy(); m = y < 4; y[m] = 3 - y[m]; return y


def long_errors(x, q, m, rng):
    L = len(x)
    n_ev = rng.poisson(10 ** (-q / 10) * L)
    n_sub, n_hp, n_ins, n_del = rng.multinomial(n_ev, m['comp'])
    x = x.copy()
    if n_sub:
        pos = rng.integers(0, L, n_sub); ok = x[pos] < 4
        x[pos[ok]] = (x[pos[ok]] + rng.integers(1, 4, ok.sum())) % 4
    cnt = np.ones(L, np.int16)
    if n_del:
        cnt[rng.integers(0, L, n_del)] = 0
    if n_hp:
        brk = np.flatnonzero(np.diff(x) != 0) + 1
        st = np.r_[0, brk]; en = np.r_[brk, L]; rl = en - st
        ok = np.flatnonzero((rl >= 3) & (x[st] < 4))
        if len(ok):
            w = (rl[ok] - 2).astype(float); w /= w.sum()
            pick = rng.choice(ok, size=n_hp, p=w)
            isdel = rng.random(n_hp) < m['hp_del']
            np.add.at(cnt, en[pick] - 1, np.where(isdel, -1, 1).astype(np.int16))
            np.maximum(cnt, 0, out=cnt)
    ins = rng.integers(0, L, n_ins) if n_ins else None
    if n_ins:
        np.add.at(cnt, ins, 1)
    out = np.repeat(x, cnt)
    if n_ins:
        last = np.cumsum(cnt.astype(np.int64)) - 1
        tg = last[ins]
        out[tg] = rng.integers(0, 4, len(tg))
    return out, n_ev


def lift_all(df, pafs, sds):
    """df has contig (query name), a, b, strand. adds <ref>_* truth columns."""
    for ref, alns in pafs.items():
        cols = {k: np.empty(len(df), dtype=object if k in ('tname', 'tname2') else np.int64)
                for k in ('tname', 't0', 't1', 'tstrand', 'cov', 'cov2', 'nrec', 'tname2', 't0_2', 't1_2', 'tstrand2')}
        for c, idx in df.groupby('contig').indices.items():
            r = lift_intervals(alns, c, df.a.values[idx], df.b.values[idx], df.strand.values[idx])
            for k in cols:
                cols[k][idx] = r[k]
        tn = np.array([strip_pansn(x) if x else '' for x in cols['tname']], dtype=object)
        tn2 = np.array([strip_pansn(x) if x else '' for x in cols['tname2']], dtype=object)
        df[ref + '_chrom'] = tn; df[ref + '_start'] = cols['t0']; df[ref + '_end'] = cols['t1']
        df[ref + '_strand'] = cols['tstrand']; df[ref + '_cov'] = cols['cov']; df[ref + '_cov2'] = cols['cov2']
        df[ref + '_nrec'] = cols['nrec']
        df[ref + '_chrom2'] = tn2; df[ref + '_start2'] = cols['t0_2']; df[ref + '_end2'] = cols['t1_2']
        df[ref + '_strand2'] = cols['tstrand2']
        has = df[ref + '_cov'].values > 0
        fr = np.zeros(len(df)); idn = np.zeros(len(df))
        if has.any():
            f_, i_ = sds[ref].query(tn[has], cols['t0'][has], cols['t1'][has])
            fr[has], idn[has] = f_, i_
        df[ref + '_sdfrac'] = np.round(fr, 4); df[ref + '_sdid'] = np.round(idn, 4)
        log(ref, 'lifted', int(has.sum()), 'of', len(df), 'reads; SD>=0.5:', int((fr >= 0.5).sum()))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tech', required=True, choices=['illumina', 'hifi', 'ont'])
    ap.add_argument('--hap', action='append', required=True)
    ap.add_argument('--paf', action='append', required=True)
    ap.add_argument('--sd', action='append', required=True)
    ap.add_argument('--n', type=int, required=True, help='pairs (illumina) or reads (long)')
    ap.add_argument('--pool', type=int, default=0, help='candidate pool size for long reads')
    ap.add_argument('--frac-e', type=float, default=0.4, help='target fraction of SD-stratum reads')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    pafs = {k: load_paf(v, qprefix=b'NA19240#') for k, v in (x.split('=', 1) for x in a.paf)}
    sds = {k: SDIndex(read_sd_table(v, k)) for k, v in (x.split('=', 1) for x in a.sd)}
    log('loaded PAFs', {k: sum(len(v) for v in d.values()) for k, d in pafs.items()},
        'SD bp', {k: v.total_bp for k, v in sds.items()})
    names, lens, hapid = [], [], []
    for h, fa in enumerate(a.hap, 1):
        n_, l_ = read_fai(fa + '.fai')
        names += n_; lens.append(l_); hapid += [h] * len(n_)
    lens = np.concatenate(lens); names = np.array(names, dtype=object); hapid = np.array(hapid)
    fa_of_hap = {h: fa for h, fa in enumerate(a.hap, 1)}

    if a.tech == 'illumina':
        n = a.n
        F = np.clip(np.rint(rng.normal(450, 60, n)), 200, 1000).astype(np.int64)
        ctg, x = sample_positions(lens, F, rng)
        flip = rng.random(n) < 0.5
        RL = 150
        # R1/R2 intervals
        left = np.c_[x, x + RL]; right = np.c_[x + F - RL, x + F]
        r1 = np.where(flip[:, None], right, left); r2 = np.where(flip[:, None], left, right)
        s1 = np.where(flip, -1, 1); s2 = -s1
        qn = np.array(['I%07d' % i for i in range(n)], dtype=object)
        df = pd.DataFrame(dict(qname=np.r_[qn, qn], mate=np.r_[np.ones(n, int), np.full(n, 2)],
                               hap=np.r_[hapid[ctg], hapid[ctg]], contig=np.r_[names[ctg], names[ctg]],
                               a=np.r_[r1[:, 0], r2[:, 0]], b=np.r_[r1[:, 1], r2[:, 1]],
                               strand=np.r_[s1, s2], frag=np.r_[F, F]))
        df['len'] = RL; df['weight'] = 1.0; df['stratumE'] = False
        df = lift_all(df, pafs, sds)
        # sequences
        seqs = np.zeros((2 * n, RL), np.uint8)
        for c, idx in df.groupby('contig').indices.items():
            h = int(df.hap.values[idx[0]])
            with pysam.FastaFile(fa_of_hap[h]) as fa:
                s = np.frombuffer(fa.fetch(c).encode(), np.uint8)
            codes = ENC[s]
            st = df.a.values[idx]
            blk = codes[st[:, None] + np.arange(RL)[None, :]]
            neg = df.strand.values[idx] < 0
            if neg.any():
                rb = blk[neg][:, ::-1].copy(); m = rb < 4; rb[m] = 3 - rb[m]; blk[neg] = rb
            seqs[idx] = blk
        nerr = rng.random(seqs.shape) < 0.002
        ok = nerr & (seqs < 4)
        seqs[ok] = (seqs[ok] + rng.integers(1, 4, ok.sum())) % 4
        df['n_err'] = nerr.sum(1)
        df['nfrac'] = (seqs == 4).mean(1)
        qual = b'?' * RL
        for mate in (1, 2):
            with gzip.open('%s_R%d.fq.gz' % (a.out, mate), 'wb', compresslevel=1) as fo:
                sel = np.flatnonzero(df.mate.values == mate)
                buf = []
                for i in sel:
                    buf.append(b'@%s\n%s\n+\n%s\n' % (df.qname.values[i].encode(), DEC[seqs[i]].tobytes(), qual))
                    if len(buf) >= 100000:
                        fo.write(b''.join(buf)); buf = []
                fo.write(b''.join(buf))
    else:
        m = MODELS[a.tech]
        npool = a.pool or 5 * a.n
        if a.tech == 'hifi':
            L = np.clip(np.rint(rng.normal(17000, 2500, npool)), 12000, 24000).astype(np.int64)
        else:
            L = np.clip(np.rint(np.exp(rng.normal(np.log(30000), 0.3, npool))), 20000, 50000).astype(np.int64)
        ctg, x = sample_positions(lens, L, rng)
        strand = np.where(rng.random(npool) < 0.5, 1, -1)
        df = pd.DataFrame(dict(contig=names[ctg], hap=hapid[ctg], a=x, b=x + L, strand=strand, len=L))
        df = lift_all(df, pafs, sds)
        E = ((df.grch38_sdfrac.values >= 0.2) | (df.chm13_sdfrac.values >= 0.2))
        nE, nN = int(E.sum()), int((~E).sum())
        tE = min(nE, int(round(a.frac_e * a.n))); tN = min(nN, a.n - tE)
        pE, pN = tE / max(1, nE), tN / max(1, nN)
        keep = np.where(E, rng.random(npool) < pE, rng.random(npool) < pN)
        log('pool', npool, 'E', nE, 'nonE', nN, 'pE %.4f pN %.4f kept %d' % (pE, pN, keep.sum()))
        df['stratumE'] = E; df['weight'] = np.where(E, 1 / pE, 1 / pN)
        df = df[keep].reset_index(drop=True)
        df['qname'] = ['%s%07d' % ('H' if a.tech == 'hifi' else 'O', i) for i in range(len(df))]
        df['mate'] = 0
        df['q'] = np.clip(rng.normal(m['qmu'], m['qsd'], len(df)), m['qlo'], m['qhi']).round(2)
        df['n_err'] = 0; df['nfrac'] = 0.0; df['out_len'] = 0
        with gzip.open(a.out + '.fa.gz', 'wb', compresslevel=1) as fo:
            for c, idx in df.groupby('contig').indices.items():
                h = int(df.hap.values[idx[0]])
                with pysam.FastaFile(fa_of_hap[h]) as fa:
                    s = ENC[np.frombuffer(fa.fetch(c).encode(), np.uint8)]
                buf = []
                for i in idx:
                    t = s[df.a.values[i]:df.b.values[i]]
                    if df.strand.values[i] < 0:
                        t = revcomp(t)
                    r, ne = long_errors(t, df.q.values[i], m, rng)
                    df.at[i, 'n_err'] = ne; df.at[i, 'nfrac'] = float((t == 4).mean()); df.at[i, 'out_len'] = len(r)
                    buf.append(b'>%s\n%s\n' % (df.qname.values[i].encode(), DEC[r].tobytes()))
                fo.write(b''.join(buf))
    cols = ['qname', 'mate', 'hap', 'contig', 'a', 'b', 'strand', 'len', 'weight', 'stratumE', 'n_err', 'nfrac']
    cols += [c for c in ('frag', 'q', 'out_len') if c in df.columns]
    cols += [c for c in df.columns if c.startswith(('grch38_', 'chm13_'))]
    df[cols].to_csv(a.out + '.truth.tsv.gz', sep='\t', index=False, compression='gzip')
    log('wrote', len(df), 'truth rows ->', a.out + '.truth.tsv.gz')


if __name__ == '__main__':
    main()
