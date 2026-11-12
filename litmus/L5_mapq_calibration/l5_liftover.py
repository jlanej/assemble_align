#!/usr/bin/env python3
"""l5_liftover.py - base-level liftover of assembly (query) intervals to reference (target)
coordinates through PAF alignments that carry a cg:Z CIGAR (=/X/M/I/D ops), e.g. the HPRC
release-2 impg/wfmash 'all-vs-1' PAFs (query = NA19240#1|2#contig, target = GRCh38/CHM13).

Truth definition used by the L5 pilot (documented in L5_mapq_calibration.md):
  * every alignment record covering a read's origin interval [a,b) on the haplotype is evaluated;
  * aligned bases = query bases of [a,b) inside =/X/M runs of that record;
  * the record with the most aligned bases is the truth; its first/last aligned query bases are
    projected to the target to give the truth interval; truth strand = read strand x record strand;
  * 'second' = aligned bases of the next-best record (flags paralogous / overlapping assembly
    alignments, i.e. ambiguous truth).
Everything is vectorised with numpy; no per-base Python loops.
"""
import gzip
import numpy as np

_OP = {ord('M'): 0, ord('='): 0, ord('X'): 0, ord('I'): 1, ord('D'): 2, ord('N'): 2}


def parse_cigar(cg: bytes):
    """Vectorised CIGAR parser -> (lengths int64, opclass int8: 0 match,1 ins(query),2 del(target))."""
    b = np.frombuffer(cg, dtype=np.uint8)
    isd = (b >= 48) & (b <= 57)
    op_pos = np.flatnonzero(~isd)
    dig_pos = np.flatnonzero(isd)
    idx = np.searchsorted(op_pos, dig_pos)            # number index of each digit
    expo = op_pos[idx] - 1 - dig_pos                   # decimal exponent of each digit
    vals = np.bincount(idx, weights=(b[dig_pos] - 48) * np.power(10.0, expo), minlength=len(op_pos))
    lut = np.full(256, -1, dtype=np.int8)
    for k, v in _OP.items():
        lut[k] = v
    ops = lut[b[op_pos]]
    if (ops < 0).any():
        raise ValueError('unsupported CIGAR op in cg:Z')
    return np.rint(vals).astype(np.int64), ops


def match_segments(L, ops):
    """Merge consecutive match ops; return (q_offset, t_offset, length) of each gap-free run,
    offsets measured from the alignment start along the CIGAR walk."""
    qc = np.where(ops != 2, L, 0)
    tc = np.where(ops != 1, L, 0)
    qoff = np.cumsum(qc) - qc
    toff = np.cumsum(tc) - tc
    ism = ops == 0
    prev = np.r_[False, ism[:-1]]
    newg = ism & ~prev
    gid = np.cumsum(newg) - 1
    mi = np.flatnonzero(ism)
    g = gid[mi]
    first = mi[np.r_[True, g[1:] != g[:-1]]] if len(mi) else mi
    seglen = np.bincount(g, weights=L[mi]).astype(np.int64) if len(mi) else np.zeros(0, np.int64)
    return qoff[first], toff[first], seglen


class Alignment:
    __slots__ = ('qname', 'qlen', 'qs', 'qe', 'strand', 'tname', 'tlen', 'ts', 'te', 'q0', 'qend',
                 'tq0', 'cum', 'ident', 'nmatch_cols')

    def __init__(self, f, cg):
        self.qname, self.qlen, self.qs, self.qe = f[0], int(f[1]), int(f[2]), int(f[3])
        self.strand = 1 if f[4] == '+' else -1
        self.tname, self.tlen, self.ts, self.te = f[5], int(f[6]), int(f[7]), int(f[8])
        L, ops = parse_cigar(cg)
        qo, to, sl = match_segments(L, ops)
        if self.strand == 1:
            q0 = self.qs + qo
            tq0 = self.ts + to                     # target pos of query pos q0
        else:                                      # CIGAR walks target forward, query backward
            qend = self.qe - qo
            q0 = qend - sl
            tq0 = self.ts + to + sl - 1            # target pos of query pos q0 (decreasing map)
            q0, tq0, sl = q0[::-1].copy(), tq0[::-1].copy(), sl[::-1].copy()
        self.q0, self.qend, self.tq0 = q0, q0 + sl, tq0
        self.cum = np.r_[0, np.cumsum(sl)]
        self.ident = None

    def t_of_q(self, i, q):
        """target coordinate of query coordinate q lying in match segment i (vectorised)."""
        return self.tq0[i] + self.strand * (q - self.q0[i])

    def lift(self, a, b):
        """For query intervals [a,b) (arrays) return aligned bases, target start, target end."""
        i1 = np.searchsorted(self.qend, a, side='right')      # first seg with qend > a
        i2 = np.searchsorted(self.q0, b, side='left') - 1     # last seg with q0 < b
        ok = i1 <= i2
        cov = np.zeros(len(a), np.int64)
        t0 = np.full(len(a), -1, np.int64)
        t1 = np.full(len(a), -1, np.int64)
        if not ok.any():
            return cov, t0, t1
        j1, j2, aa, bb = i1[ok], i2[ok], a[ok], b[ok]
        c = self.cum[j2 + 1] - self.cum[j1]
        c -= np.maximum(0, aa - self.q0[j1])
        c -= np.maximum(0, self.qend[j2] - bb)
        qa = np.maximum(aa, self.q0[j1])
        qb = np.minimum(bb, self.qend[j2]) - 1
        ta, tb = self.t_of_q(j1, qa), self.t_of_q(j2, qb)
        cov[ok] = c
        t0[ok] = np.minimum(ta, tb)
        t1[ok] = np.maximum(ta, tb) + 1
        return cov, t0, t1


def load_paf(path, qprefix=None, min_len=0):
    """Load PAF records (with cg:Z) grouped by query name. Returns dict qname -> [Alignment]."""
    out = {}
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rb') as fh:
        for line in fh:
            if qprefix and not line.startswith(qprefix):
                continue
            f = line.rstrip(b'\n').split(b'\t')
            cg = None
            ident = None
            for x in f[12:]:
                if x.startswith(b'cg:Z:'):
                    cg = x[5:]
                elif x.startswith(b'gi:f:'):
                    ident = float(x[5:])
            if cg is None:
                raise ValueError('PAF record without cg:Z tag')
            fs = [y.decode() for y in f[:12]]
            if int(fs[3]) - int(fs[2]) < min_len:
                continue
            aln = Alignment(fs, cg)
            aln.ident = ident
            out.setdefault(aln.qname, []).append(aln)
    for k in out:
        out[k].sort(key=lambda r: r.qs)
    return out


def lift_intervals(alns_by_q, qname, a, b, read_strand):
    """Lift reads on query contig `qname` (arrays a,b,read_strand=+1/-1).
    Returns dict of arrays: best record (tname, t0, t1, tstrand, cov), runner-up record
    (tname2, t0_2, t1_2, tstrand2, cov2) and the number of overlapping records (nrec)."""
    n = len(a)
    res = dict(tname=np.full(n, '', dtype=object), t0=np.full(n, -1, np.int64),
               t1=np.full(n, -1, np.int64), tstrand=np.zeros(n, np.int8),
               cov=np.zeros(n, np.int64), cov2=np.zeros(n, np.int64), nrec=np.zeros(n, np.int32),
               ident=np.full(n, np.nan),
               tname2=np.full(n, '', dtype=object), t0_2=np.full(n, -1, np.int64),
               t1_2=np.full(n, -1, np.int64), tstrand2=np.zeros(n, np.int8))
    alns = alns_by_q.get(qname, [])
    if not alns or n == 0:
        return res
    order = np.argsort(a, kind='stable')
    sa = a[order]
    maxlen = int((b - a).max())
    for aln in alns:
        lo = np.searchsorted(sa, aln.qs - maxlen, side='left')
        hi = np.searchsorted(sa, aln.qe, side='left')
        if hi <= lo:
            continue
        cand = order[lo:hi]
        cand = cand[b[cand] > aln.qs]
        if len(cand) == 0:
            continue
        cov, t0, t1 = aln.lift(a[cand], b[cand])
        has = cov > 0
        cand, cov, t0, t1 = cand[has], cov[has], t0[has], t1[has]
        if len(cand) == 0:
            continue
        res['nrec'][cand] += 1
        better = cov > res['cov'][cand]
        cb = cand[better]
        # demote current best to second
        for k1, k2 in (('tname', 'tname2'), ('t0', 't0_2'), ('t1', 't1_2'), ('tstrand', 'tstrand2')):
            res[k2][cb] = res[k1][cb]
        res['cov2'][cb] = res['cov'][cb]
        res['cov'][cb] = cov[better]
        res['t0'][cb] = t0[better]
        res['t1'][cb] = t1[better]
        res['tname'][cb] = aln.tname
        res['tstrand'][cb] = (read_strand[cb] * aln.strand).astype(np.int8)
        res['ident'][cb] = aln.ident if aln.ident is not None else np.nan
        worse = ~better
        sec = worse & (cov > res['cov2'][cand])
        cs_ = cand[sec]
        res['cov2'][cs_] = cov[sec]
        res['t0_2'][cs_] = t0[sec]
        res['t1_2'][cs_] = t1[sec]
        res['tname2'][cs_] = aln.tname
        res['tstrand2'][cs_] = (read_strand[cs_] * aln.strand).astype(np.int8)
    return res


def strip_pansn(name):
    """GRCh38#0#chr1 -> chr1 ; CHM13#0#chr1 -> chr1"""
    return name.split('#')[-1] if '#' in name else name


if __name__ == '__main__':   # tiny self-test on synthetic CIGARs
    f = ['q', '100', '10', '60', '+', 't', '1000', '500', '552', '0', '0', '60']
    aln = Alignment(f, b'10=2I8=1X9=4D20=')   # query 50 bases (10+2+8+1+9+20), target 52
    cov, t0, t1 = aln.lift(np.array([10, 20, 30]), np.array([20, 32, 60]))
    print('fwd', cov.tolist(), t0.tolist(), t1.tolist())   # expect [10,10,30] [500,510,522]...
    f[4] = '-'
    aln = Alignment(f, b'10=2I8=1X9=4D20=')
    cov, t0, t1 = aln.lift(np.array([40, 50]), np.array([50, 60]))
    print('rev', cov.tolist(), t0.tolist(), t1.tolist())   # q[50,60) <-> first 10 target bases
