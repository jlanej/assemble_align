"""L1_eval.py -- index construction (minimap2-style occurrence filter) and per-read seed evaluation.

Index: canonical 2-bit codes of sampled reference k-mers, sorted; occurrence = number of index
entries with the same code.  mid_occ follows minimap2 v2.31 mm_idx_cal_max_occ()+mm_mapopt_update():
   mid_occ = (1-f)-quantile (ks_ksmall at index (1-f)*n) of per-distinct-minimizer occurrence + 1,
   clamped to [min_mid_occ, max_mid_occ]; f = 2e-4.  Presets: map-ont default [10, 1e6];
   map-hifi/lr:hq and all asm* presets (asm5/asm10/asm20) [50, 500]; lr:hqae [50, 500]; sr fixed mid_occ = 1000.
   (minimap2's occ_dist rescue of a few high-occurrence seeds per 500 bp is NOT modelled.)

Truth for a read seed (k-mer ending at read base P): source coordinates of its two end bases e0,e1 are
lifted to the reference (identity for CHM13-derived reads; PAF/cg liftover for NA19240 haplotypes);
expected reference end = max(lift(e0), lift(e1)); expected relative strand = read orientation XOR
alignment strand.  A seed is TRUE if any of its index hits has that relative strand and lies within
`tol` bp (20; k+10 if one end is unliftable) of the expected end.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import numba as nb
import L1_schemes as S

PRESET_KW = {'sr': (21, 11), 'map-ont': (15, 10), 'map-hifi': (19, 19), 'asm20': (19, 10),
             'asm5': (19, 19), 'lr:hqae': (25, 51)}
PRESET_OCC = {'sr': ('fixed', 1000, 1000), 'map-ont': ('frac', 10, 1000000), 'asm20': ('frac', 50, 500),
              'map-hifi': ('frac', 50, 500), 'asm5': ('frac', 50, 500), 'lr:hqae': ('frac', 50, 500)}
MID_OCC_FRAC = 2e-4


def parse_cigar(cg):
    import re
    return [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cg)]


def load_lift(paf, qlen):
    """Per-base liftover of query (haplotype) coordinates onto target (reference) from a PAF with cg:Z.
    Overlaps resolved by (mapq, aligned length) priority.  Returns pos(int32,-1), strand(int8)."""
    recs = []
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if cg is None:
            continue
        recs.append((int(f[11]), int(f[10]), int(f[2]), int(f[3]), f[4], int(f[7]), cg))
    recs.sort(key=lambda r: (-r[0], -r[1]))
    pos = np.full(qlen, -1, np.int32)
    st = np.zeros(qlen, np.int8)
    for mapq, alen, qs, qe, strand, ts, cg in recs:
        ops = parse_cigar(cg)
        t = ts
        if strand == '+':
            q = qs
            for n, op in ops:
                if op in 'M=X':
                    sl = slice(q, q + n)
                    unset = pos[sl] < 0
                    seg = np.arange(t, t + n, dtype=np.int32)
                    pos[sl][unset] = seg[unset]; st[sl][unset] = 0
                    q += n; t += n
                elif op == 'I':
                    q += n
                elif op in 'DN':
                    t += n
        else:
            q = qe            # walk query backwards
            for n, op in ops:
                if op in 'M=X':
                    qq = np.arange(q - 1, q - 1 - n, -1)
                    seg = np.arange(t, t + n, dtype=np.int32)
                    unset = pos[qq] < 0
                    pos[qq[unset]] = seg[unset]; st[qq[unset]] = 1
                    q -= n; t += n
                elif op == 'I':
                    q -= n
                elif op in 'DN':
                    t += n
    return pos, st


def mid_occ_from_counts(counts, preset):
    mode, lo, hi = PRESET_OCC[preset]
    if mode == 'fixed':
        return lo
    n = len(counts)
    kth = int((1.0 - MID_OCC_FRAC) * n)
    thr = int(np.partition(counts, kth)[kth]) + 1
    thr = max(thr, lo)
    if hi > lo and thr > hi:
        thr = hi
    return thr


class Index:
    def __init__(self, codes, scheme, k, w, preset, s=None, greedy=None, chunk=3_000_000):
        pos, st, cn = S.sample2_chunked(codes, scheme, k, w, s=s, greedy=greedy, chunk=chunk)
        ps = np.sort(pos)
        gaps = np.diff(ps)
        self.max_gap = int(gaps.max()) if len(gaps) else 0
        self.mean_gap = float(gaps.mean()) if len(gaps) else 0.0
        del ps, gaps
        o = np.argsort(cn, kind='stable')
        self.key = cn[o]; del cn
        self.pos = pos[o].astype(np.int32); del pos
        self.st = st[o]; del st, o
        self.k, self.w, self.scheme, self.preset = k, w, scheme, preset
        brk = np.flatnonzero(self.key[1:] != self.key[:-1]) + 1          # keys are sorted: run lengths = counts
        cnt = np.diff(np.concatenate(([0], brk, [len(self.key)])))
        del brk
        self.n_entries = len(self.key); self.n_distinct = len(cnt)
        self.mid_occ = mid_occ_from_counts(cnt, preset)
        nvalid = int(np.count_nonzero(np.convolve((codes < 4).astype(np.int8), np.ones(k, np.int8), 'valid') == k))
        self.density = self.n_entries / nvalid
        hf = cnt > self.mid_occ
        self.frac_entries_hf = float(cnt[hf].sum() / self.n_entries)
        self.frac_distinct_hf = float(hf.mean())
        self.occ_sq = float((cnt.astype(np.float64) ** 2).sum() / self.n_entries)   # mean occ per entry

    def summary(self):
        return dict(scheme=self.scheme, preset=self.preset, k=self.k, w=self.w, entries=self.n_entries,
                    distinct=self.n_distinct, density=self.density, dfactor=self.density * (self.w + 1),
                    mid_occ=self.mid_occ, frac_entries_hf=self.frac_entries_hf,
                    frac_distinct_hf=self.frac_distinct_hf, mean_occ_per_entry=self.occ_sq,
                    ref_max_gap=self.max_gap, ref_mean_gap=self.mean_gap)


@nb.njit(cache=True)
def _eval_core(key, rpos, rst, mid_occ, seed_key, seed_P, seed_z, seed_read, offs, lens, orient,
               src, lift_pos, lift_st, identity, k, tol, n_reads):
    n_seeds = np.zeros(n_reads, np.int64); n_hf = np.zeros(n_reads, np.int64)
    n_true = np.zeros(n_reads, np.int64); n_true_any = np.zeros(n_reads, np.int64)
    anchors = np.zeros(n_reads, np.int64); max_gap = np.zeros(n_reads, np.int64)
    liftable = np.zeros(n_reads, np.int64)
    last_cov = np.full(n_reads, -1, np.int64)      # last covered read base (qpos) of true seeds
    nk = key.shape[0]
    for i in range(seed_key.shape[0]):
        r = seed_read[i]
        P = seed_P[i]
        q = P - offs[r]
        n_seeds[r] += 1
        c = seed_key[i]
        lo = np.searchsorted(key, c, 'left'); hi = np.searchsorted(key, c, 'right')
        occ = hi - lo
        # truth
        e1 = src[P]; e0 = src[P - k + 1]
        exp = -1; so = 0; one = False
        if identity:
            if e1 >= 0 and e0 >= 0:
                exp = max(e0, e1)
            elif e1 >= 0 or e0 >= 0:
                exp = max(e0, e1); one = True
        else:
            l1 = lift_pos[e1] if e1 >= 0 else -1
            l0 = lift_pos[e0] if e0 >= 0 else -1
            if l1 >= 0 and l0 >= 0:
                exp = max(l0, l1); so = lift_st[e1]
            elif l1 >= 0:
                exp = l1; so = lift_st[e1]; one = True
            elif l0 >= 0:
                exp = l0; so = lift_st[e0]; one = True
        if exp >= 0:
            liftable[r] += 1
        rel = orient[r] ^ so
        t = tol + (k if one else 0)
        is_true = False
        if exp >= 0 and occ > 0:
            for j in range(lo, hi):
                if (seed_z[i] ^ rst[j]) == rel and abs(rpos[j] - exp) <= t:
                    is_true = True
                    break
        if occ > mid_occ:
            n_hf[r] += 1
            if is_true:
                n_true_any[r] += 1
        else:
            anchors[r] += occ
            if is_true:
                n_true[r] += 1; n_true_any[r] += 1
                # seedless gap bookkeeping (seeds arrive in increasing P within a read)
                st_ = q - k + 1
                g = st_ - (last_cov[r] + 1)
                if g > max_gap[r]:
                    max_gap[r] = g
                if q > last_cov[r]:
                    last_cov[r] = q
    for r in range(n_reads):
        g = lens[r] - (last_cov[r] + 1)
        if g > max_gap[r]:
            max_gap[r] = g
    return n_seeds, n_hf, n_true, n_true_any, anchors, max_gap, liftable


def eval_reads(index, rs, scheme, k, w, lift=None, tol=20, s=None, greedy=None, chunk=3_000_000):
    P, z, key = S.sample2_chunked(rs['codes'], scheme, k, w, s=s, greedy=greedy, chunk=chunk)
    o = np.argsort(P, kind='stable'); P, z, key = P[o], z[o], key[o]
    rid = np.searchsorted(rs['offs'], P, 'right') - 1
    ok = (P - rs['offs'][rid]) < rs['lens'][rid]
    P, z, key, rid = P[ok], z[ok], key[ok], rid[ok]
    if lift is None:
        lp = np.zeros(1, np.int32); ls = np.zeros(1, np.int8); ident = True
    else:
        lp, ls = lift; ident = False
    out = _eval_core(index.key, index.pos, index.st, index.mid_occ, key, P, z, rid,
                     rs['offs'], rs['lens'], rs['meta']['orient'].astype(np.int64), rs['src'], lp, ls,
                     ident, k, tol, len(rs['offs']))
    names = ['n_seeds', 'n_hf', 'n_true', 'n_true_any', 'anchors', 'max_gap', 'liftable']
    return dict(zip(names, out))
