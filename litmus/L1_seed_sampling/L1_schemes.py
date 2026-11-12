"""L1_schemes.py -- k-mer sampling schemes on 2-bit DNA, strand-symmetric (canonical) by default.

All schemes emit positions as minimap2 does: the coordinate of the LAST base of the sampled k-mer,
plus the strand bit z (0 if the forward 2-bit code < reverse-complement code), plus the canonical
2-bit k-mer code (the index key; minimap2 keys on hash64(canonical), an invertible map of it).

Schemes
  mm2        exact port of minimap2 v2.31 mm_sketch() (hash64 of canonical k-mer, robust winnowing:
             every copy of the window minimum is emitted; palindromic k-mers skipped; N resets).
  <order>    same window logic as mm2 but with the k-mer order replaced by a scheme-specific rank:
             miniception (closed syncmers first, r=max(4,k-w)), ddecyc (double decycling, Pellow et al. 2023),
             oc (open-closed minimizer, Groot Koerkamp et al. 2025), greedy (GreedyMini binary order,
             Golan et al. 2025, with k'-prefix extension for k>16). Random hash64 tie-break inside classes.
  mod / modoc  mod-sampling (Groot Koerkamp & Pibiri 2024): smallest t-mer (random, or open-closed order)
             at window offset x selects the k-mer at x mod w; t = r + ((k-r) mod w), r = 4.
  osync/csync  open / closed syncmers (Edgar 2021); context-free; offset v=floor((k-s)/2) (Shaw & Yu 2022).

Strand: canonical (strand-symmetric) versions are used for indexing. Syncmer-type tests use canonical
s-mer hashes and the inclusive rule "the minimum value occurs at position p", with p and its mirror
k-s-p both accepted, so the decision is identical on both strands. Orders without a symmetric
construction (ddecyc, greedy) use rank_c(x) = min(rank(x), rank(rc(x))). Mod-sampling with canonical
t-mer ranks is strand-consistent because k-t is a multiple of w (see L1_sampling_schemes.md).
"""
import math
import numpy as np
import numba as nb

U64 = np.uint64
U64MAX = np.uint64(0xFFFFFFFFFFFFFFFF)

NT4 = np.full(256, 4, np.uint8)
for _c, _v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    NT4[_c] = _v
NT4[ord('U')] = 3; NT4[ord('u')] = 3          # as in minimap2's seq_nt4_table


def encode(seq):
    if isinstance(seq, str):
        seq = seq.encode()
    return NT4[np.frombuffer(seq, np.uint8)]


def read_fasta(path):
    """Return list of (name, bytes) for a (possibly multi-record) uncompressed FASTA."""
    out, name, buf = [], None, []
    with open(path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                if name is not None:
                    out.append((name, b''.join(buf)))
                name, buf = line[1:].split()[0].decode(), []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, b''.join(buf)))
    return out


@nb.njit(inline='always')
def hash64(key, mask):
    key = (~key + (key << U64(21))) & mask
    key = key ^ (key >> U64(24))
    key = ((key + (key << U64(3))) + (key << U64(8))) & mask
    key = key ^ (key >> U64(14))
    key = ((key + (key << U64(2))) + (key << U64(4))) & mask
    key = key ^ (key >> U64(28))
    key = (key + (key << U64(31))) & mask
    return key


# ----------------------------------------------------------------------------------------------
# Per-position k-mer codes.  Position i = last base of the k-mer.  valid[i]=1 iff the k bases are ACGT.
# ----------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def kmer_codes(codes, k):
    n = codes.shape[0]
    mask = (U64(1) << U64(2 * k)) - U64(1)
    shift1 = U64(2 * (k - 1))
    fwd = np.zeros(n, np.uint64)
    rev = np.zeros(n, np.uint64)
    valid = np.zeros(n, np.uint8)
    f = U64(0); r = U64(0); l = 0
    for i in range(n):
        c = codes[i]
        if c < 4:
            f = ((f << U64(2)) | U64(c)) & mask
            r = (r >> U64(2)) | (U64(3 - c) << shift1)
            l += 1
            if l >= k:
                fwd[i] = f; rev[i] = r; valid[i] = 1
        else:
            l = 0
    return fwd, rev, valid


@nb.njit(cache=True)
def canon_hash(fwd, rev, valid, k):
    """minimap2 hash64 of the canonical k-mer (UINT64_MAX where invalid)."""
    n = fwd.shape[0]
    mask = (U64(1) << U64(2 * k)) - U64(1)
    out = np.full(n, U64MAX, np.uint64)
    for i in range(n):
        if valid[i]:
            c = fwd[i] if fwd[i] < rev[i] else rev[i]
            out[i] = hash64(c, mask)
    return out


# ----------------------------------------------------------------------------------------------
# Exact port of minimap2 v2.31 mm_sketch() (non-HPC) with optional replacement of the k-mer order.
# If use_rank: info.x = rank[i] instead of hash64(kmer[z])<<8 | span.  rank must be < UINT64_MAX and
# injective on canonical k-mers (so that equal x <=> identical canonical k-mer, as in minimap2).
# Returns x array, y array (pos<<1 | strand), as minimap2 would push them.
# ----------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def mm_sketch_generic(codes, w, k, rank, use_rank):
    n = codes.shape[0]
    shift1 = U64(2 * (k - 1))
    mask = (U64(1) << U64(2 * k)) - U64(1)
    km0 = U64(0); km1 = U64(0)
    bufx = np.full(w, U64MAX, np.uint64)
    bufy = np.full(w, U64MAX, np.uint64)
    cap = n // max(1, w // 2) + 1024
    outx = np.empty(cap, np.uint64); outy = np.empty(cap, np.uint64)
    no = 0
    minx = U64MAX; miny = U64MAX
    l = 0; buf_pos = 0; min_pos = 0
    kspan = k
    for i in range(n):
        c = codes[i]
        infox = U64MAX; infoy = U64MAX
        if c < 4:
            km0 = ((km0 << U64(2)) | U64(c)) & mask
            km1 = (km1 >> U64(2)) | (U64(3 - c) << shift1)
            if km0 == km1:
                continue          # symmetric k-mer: skipped entirely (no buffer slot), as in minimap2
            z = 0 if km0 < km1 else 1
            l += 1
            if l >= k:
                if use_rank:
                    infox = rank[i]
                else:
                    kk = km0 if z == 0 else km1
                    infox = (hash64(kk, mask) << U64(8)) | U64(kspan)
                infoy = (U64(i) << U64(1)) | U64(z)
        else:
            l = 0
        bufx[buf_pos] = infox; bufy[buf_pos] = infoy
        if l == w + k - 1 and minx != U64MAX:
            for j in range(buf_pos + 1, w):
                if minx == bufx[j] and bufy[j] != miny:
                    outx[no] = bufx[j]; outy[no] = bufy[j]; no += 1
            for j in range(0, buf_pos):
                if minx == bufx[j] and bufy[j] != miny:
                    outx[no] = bufx[j]; outy[no] = bufy[j]; no += 1
        if infox <= minx:
            if l >= w + k and minx != U64MAX:
                outx[no] = minx; outy[no] = miny; no += 1
            minx = infox; miny = infoy; min_pos = buf_pos
        elif buf_pos == min_pos:
            if l >= w + k - 1 and minx != U64MAX:
                outx[no] = minx; outy[no] = miny; no += 1
            minx = U64MAX
            for j in range(buf_pos + 1, w):
                if minx >= bufx[j]:
                    minx = bufx[j]; miny = bufy[j]; min_pos = j
            for j in range(0, buf_pos + 1):
                if minx >= bufx[j]:
                    minx = bufx[j]; miny = bufy[j]; min_pos = j
            if l >= w + k - 1 and minx != U64MAX:
                for j in range(buf_pos + 1, w):
                    if minx == bufx[j] and miny != bufy[j]:
                        outx[no] = bufx[j]; outy[no] = bufy[j]; no += 1
                for j in range(0, buf_pos + 1):
                    if minx == bufx[j] and miny != bufy[j]:
                        outx[no] = bufx[j]; outy[no] = bufy[j]; no += 1
        buf_pos += 1
        if buf_pos == w:
            buf_pos = 0
        if no > cap - 2 * w - 4:      # grow
            nx = np.empty(cap * 2, np.uint64); ny = np.empty(cap * 2, np.uint64)
            nx[:no] = outx[:no]; ny[:no] = outy[:no]
            outx = nx; outy = ny; cap *= 2
    if minx != U64MAX:
        outx[no] = minx; outy[no] = miny; no += 1
    return outx[:no], outy[:no]


# ----------------------------------------------------------------------------------------------
# s-mer based classes (syncmer tests) with canonical s-mer hashes and inclusive tie rule.
# For k-mer ending at i, its s-mers end at i-k+s .. i (relative r = 0..k-s).
#   closed: min occurs at r=0 or r=k-s ; open: min occurs at r=v or r=k-s-v, v=floor((k-s)/2)
# Sliding-window minimum (monotone deque) over the s-mer hash array.
# If canonical=False: forward s-mer hashes, leftmost-min position (as in the reference Rust code).
# ----------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def smer_classes(codes, k, s, canonical, seed):
    n = codes.shape[0]
    fs, rs, vs = kmer_codes(codes, s)
    smask = (U64(1) << U64(2 * s)) - U64(1)
    hs = np.full(n, U64MAX, np.uint64)
    for j in range(n):
        if vs[j]:
            c = fs[j]
            if canonical and rs[j] < c:
                c = rs[j]
            hs[j] = hash64(c ^ (U64(seed) & smask), smask)
    W = k - s + 1
    v = (k - s) // 2
    is_open = np.zeros(n, np.uint8); is_closed = np.zeros(n, np.uint8)
    dq = np.empty(n, np.int64); head = 0; tail = 0
    for j in range(n):
        # push s-mer ending at j
        if canonical:
            while tail > head and hs[dq[tail - 1]] > hs[j]:
                tail -= 1
        else:
            while tail > head and hs[dq[tail - 1]] > hs[j]:   # keep equal -> front is leftmost min
                tail -= 1
        dq[tail] = j; tail += 1
        first = j - W + 1          # s-mer index (end pos) of r=0 for the k-mer ending at i=j
        while dq[head] < first:
            head += 1
        if first >= 0 and hs[j] != U64MAX:
            # k-mer ending at j must be valid: all its s-mers valid <=> codes ok; check ends
            m = hs[dq[head]]
            if m == U64MAX:
                continue
            # quick validity: the k-mer is valid iff no N within; s-mers all valid and W of them
            ok = True
            # (validity of the whole k-mer is checked by caller via kmer valid array)
            if canonical:
                if hs[first] == m or hs[j] == m:
                    is_closed[j] = 1
                if hs[first + v] == m or hs[first + (k - s - v)] == m:
                    is_open[j] = 1
            else:
                p = dq[head] - first      # leftmost minimum position
                if p == 0 or p == k - s:
                    is_closed[j] = 1
                if p == v:
                    is_open[j] = 1
    return is_open, is_closed


@nb.njit(cache=True)
def class_rank(chash, valid, cls):
    """rank = class<<58 | hash (hash < 2^56 for k<=28); invalid -> UINT64_MAX."""
    n = chash.shape[0]
    out = np.full(n, U64MAX, np.uint64)
    for i in range(n):
        if valid[i]:
            out[i] = (U64(cls[i]) << U64(58)) | chash[i]
    return out


@nb.njit(cache=True)
def decycling_class(km, k):
    """Double-decycling class of one forward 2-bit k-mer (A0 C1 G2 T3), as in RagnarGrootKoerkamp/minimizers:
    x = sum_j c_j * exp(2*pi*i*j/k);  arg in [pi-2pi/k, pi) -> 0 ; arg in [-2pi/k, 0) -> 1 ; else 2."""
    re = 0.0; im = 0.0
    for j in range(k):
        c = (km >> U64(2 * (k - 1 - j))) & U64(3)
        a = 2.0 * math.pi * j / k
        re += float(c) * math.cos(a); im += float(c) * math.sin(a)
    arg = math.atan2(im, re)
    if math.pi - 2.0 * math.pi / k <= arg and arg < math.pi:
        return 0
    if -2.0 * math.pi / k <= arg and arg < 0.0:
        return 1
    return 2


@nb.njit(cache=True)
def decycling_rank(fwd, rev, valid, chash, k, canonical):
    n = fwd.shape[0]
    out = np.full(n, U64MAX, np.uint64)
    for i in range(n):
        if valid[i]:
            c = decycling_class(fwd[i], k)
            if canonical:
                c2 = decycling_class(rev[i], k)
                if c2 < c:
                    c = c2
            out[i] = (U64(c) << U64(58)) | chash[i]
    return out


@nb.njit(cache=True)
def _greedy_rank1(km, k, kp, order):
    """GreedyMini DNA rank of one forward k-mer: binary order on the high bits (base>>1) of the first
    kp bases (k'-prefix extension for k>kp), then the full high-bit pattern (lex), then the low bits (lex)."""
    hi = U64(0); lo = U64(0)
    for j in range(k):
        c = (km >> U64(2 * (k - 1 - j))) & U64(3)
        hi = (hi << U64(1)) | (c >> U64(1))
        lo = (lo << U64(1)) | (c & U64(1))
    pre = hi >> U64(k - kp)
    return (U64(order[pre]) << U64(2 * k)) | (hi << U64(k)) | lo


@nb.njit(cache=True)
def greedy_rank(fwd, rev, valid, k, kp, order, canonical):
    n = fwd.shape[0]
    out = np.full(n, U64MAX, np.uint64)
    for i in range(n):
        if valid[i]:
            r = _greedy_rank1(fwd[i], k, kp, order)
            if canonical:
                r2 = _greedy_rank1(rev[i], k, kp, order)
                if r2 < r:
                    r = r2
            out[i] = r
    return out


def load_greedy_order(path):
    """GreedyMini 'ranks' text file: python-style list of 2^k' ints (value 2^k' = unranked)."""
    txt = open(path).read().strip()
    arr = np.array([int(x) for x in txt.strip('[]').split(',')], dtype=np.uint64)
    kp = int(round(math.log2(len(arr))))
    assert len(arr) == 1 << kp
    return arr, kp


# ----------------------------------------------------------------------------------------------
# Mod-sampling.  Window = l = w+k-1 bases (w k-mers, w+k-t t-mers).  x = offset of the smallest t-mer
# (by trank); the k-mer at offset (x mod w) is sampled.  Ties (identical canonical t-mers) are broken
# by the canonical hash of the k-mer each candidate would select (strand-symmetric), then leftmost.
# Windows containing an invalid base are skipped.  Returns sorted unique k-mer end positions.
# ----------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def mod_sample(codes, w, k, t, trank, khash, kvalid):
    n = codes.shape[0]
    L = w + k - 1
    nt = w + k - t
    out = np.empty(n // max(1, w // 2) + 1024, np.int64); no = 0
    dq = np.empty(n, np.int64); head = 0; tail = 0
    run = 0                         # current N-free run length
    last = -1
    for e in range(n):
        if codes[e] >= 4:
            run = 0; head = 0; tail = 0
            continue
        run += 1
        if run >= t:
            # t-mer ending at e enters
            while tail > head and trank[dq[tail - 1]] > trank[e]:
                tail -= 1
            dq[tail] = e; tail += 1
        if run >= L:
            first_t = e - nt + 1                # end position of the t-mer at offset x=0
            while dq[head] < first_t:
                head += 1
            m = trank[dq[head]]
            wstart = e - L + 1
            best = -1; bh = U64MAX
            h = head
            while h < tail and trank[dq[h]] == m:
                x = dq[h] - (t - 1) - wstart      # offset of this t-mer in the window
                q = x % w
                kend = wstart + q + k - 1
                if kvalid[kend] and (best < 0 or khash[kend] < bh):
                    best = kend; bh = khash[kend]
                h += 1
            if best >= 0 and best != last:
                out[no] = best; no += 1; last = best
                if no >= out.shape[0] - 2:
                    nout = np.empty(out.shape[0] * 2, np.int64); nout[:no] = out[:no]; out = nout
    res = np.unique(out[:no])
    return res


def tmer_rank(codes, t, kind, s=4, canonical=True):
    fw, rv, va = kmer_codes(codes, t)
    th = canon_hash(fw, rv, va, t) if canonical else fwd_hash(fw, va, t)
    if kind == 'random':
        return th, va
    if kind == 'oc':
        so, sc = smer_classes(codes, t, s, canonical, 0x5bd1e995)
        cls = np.where(so == 1, 0, np.where(sc == 1, 1, 2)).astype(np.uint8)
        return class_rank(th, va, cls), va
    raise ValueError(kind)


@nb.njit(cache=True)
def fwd_hash(fwd, valid, k):
    mask = (U64(1) << U64(2 * k)) - U64(1)
    out = np.full(fwd.shape[0], U64MAX, np.uint64)
    for i in range(fwd.shape[0]):
        if valid[i]:
            out[i] = hash64(fwd[i], mask)
    return out


# ----------------------------------------------------------------------------------------------
# Driver: return (pos[int64], strand[uint8], canon[uint64]) of sampled k-mers for a scheme.
# ----------------------------------------------------------------------------------------------
def mod_t(k, w, r=4):
    return k if k < r else r + ((k - r) % w)


def sample(codes, scheme, k, w, s=None, greedy=None, canonical=True):
    fwd, rev, valid = kmer_codes(codes, k)
    can = np.minimum(fwd, rev)
    strand = (fwd > rev).astype(np.uint8)
    if scheme == 'mm2':
        x, y = mm_sketch_generic(codes, w, k, np.empty(1, np.uint64), False)
        pos = (y >> U64(1)).astype(np.int64)
        return pos, (y & U64(1)).astype(np.uint8), can[pos]
    chash = canon_hash(fwd, rev, valid, k) if canonical else fwd_hash(fwd, valid, k)
    if scheme in ('random', 'miniception', 'oc', 'ddecyc', 'greedy'):
        if scheme == 'random':
            rank = chash
        elif scheme in ('miniception', 'oc'):
            r = max(4, k - w) if scheme == 'miniception' else (s if s is not None else 4)
            so, sc = smer_classes(codes, k, r, canonical, 0x5bd1e995)
            if scheme == 'miniception':
                cls = np.where(sc == 1, 0, 1).astype(np.uint8)
            else:
                cls = np.where(so == 1, 0, np.where(sc == 1, 1, 2)).astype(np.uint8)
            rank = class_rank(chash, valid, cls)
        elif scheme == 'ddecyc':
            rank = decycling_rank(fwd, rev, valid, chash, k, canonical)
        else:
            order, kp = greedy
            kp_use = min(kp, k)
            rank = greedy_rank(fwd, rev, valid, k, kp_use, order, canonical)
        if not canonical:
            raise NotImplementedError('forward-only order schemes: use sample_forward()')
        x, y = mm_sketch_generic(codes, w, k, rank, True)
        pos = (y >> U64(1)).astype(np.int64)
        return pos, (y & U64(1)).astype(np.uint8), can[pos]
    if scheme in ('mod', 'modoc'):
        t = mod_t(k, w)
        trank, tval = tmer_rank(codes, t, 'random' if scheme == 'mod' else 'oc', 4, canonical)
        pos = mod_sample(codes, w, k, t, trank, chash, valid)
        return pos, strand[pos], can[pos]
    if scheme in ('osync', 'csync'):
        so, sc = smer_classes(codes, k, s, canonical, 0x5bd1e995)
        sel = (so if scheme == 'osync' else sc).astype(bool) & valid.astype(bool)
        pos = np.nonzero(sel)[0].astype(np.int64)
        return pos, strand[pos], can[pos]
    raise ValueError(scheme)


@nb.njit(cache=True)
def window_min_select_fwd(rank, valid, w, k, n):
    """Forward (non-canonical) leftmost-min window selection, one k-mer per window (reference semantics)."""
    out = np.empty(n, np.int64); no = 0
    dq = np.empty(n, np.int64); head = 0; tail = 0
    run = 0; last = -1
    for i in range(n):
        if not valid[i]:
            run = 0; head = 0; tail = 0
            continue
        run += 1
        while tail > head and rank[dq[tail - 1]] > rank[i]:
            tail -= 1
        dq[tail] = i; tail += 1
        if run >= w:
            while dq[head] < i - w + 1:
                head += 1
            p = dq[head]
            if p != last:
                out[no] = p; no += 1; last = p
    return np.unique(out[:no])


def sample_forward(codes, scheme, k, w, s=None, greedy=None):
    """Forward-strand-only version of each scheme (paper semantics: leftmost tie-break, one k-mer per window)."""
    fwd, rev, valid = kmer_codes(codes, k)
    fh = fwd_hash(fwd, valid, k)
    n = codes.shape[0]
    if scheme in ('random', 'mm2'):
        rank = fh
    elif scheme in ('miniception', 'oc'):
        r = max(4, k - w) if scheme == 'miniception' else (s if s is not None else 4)
        so, sc = smer_classes(codes, k, r, False, 0x5bd1e995)
        cls = np.where(sc == 1, 0, 1) if scheme == 'miniception' else np.where(so == 1, 0, np.where(sc == 1, 1, 2))
        rank = class_rank(fh, valid, cls.astype(np.uint8))
    elif scheme == 'ddecyc':
        rank = decycling_rank(fwd, rev, valid, fh, k, False)
    elif scheme == 'greedy':
        order, kp = greedy
        rank = greedy_rank(fwd, rev, valid, k, min(kp, k), order, False)
    elif scheme in ('mod', 'modoc'):
        t = mod_t(k, w)
        trank, tval = tmer_rank(codes, t, 'random' if scheme == 'mod' else 'oc', 4, False)
        return mod_sample_fwd(codes, w, k, t, trank, valid)
    elif scheme in ('osync', 'csync'):
        so, sc = smer_classes(codes, k, s, False, 0x5bd1e995)
        sel = (so if scheme == 'osync' else sc).astype(bool) & valid.astype(bool)
        return np.nonzero(sel)[0].astype(np.int64)
    else:
        raise ValueError(scheme)
    return window_min_select_fwd(rank, valid, w, k, n)


@nb.njit(cache=True)
def mod_sample_fwd(codes, w, k, t, trank, kvalid):
    n = codes.shape[0]
    L = w + k - 1
    nt = w + k - t
    out = np.empty(n, np.int64); no = 0
    dq = np.empty(n, np.int64); head = 0; tail = 0
    run = 0; last = -1
    for e in range(n):
        if codes[e] >= 4:
            run = 0; head = 0; tail = 0
            continue
        run += 1
        if run >= t:
            while tail > head and trank[dq[tail - 1]] > trank[e]:
                tail -= 1
            dq[tail] = e; tail += 1
        if run >= L:
            first_t = e - nt + 1
            while dq[head] < first_t:
                head += 1
            wstart = e - L + 1
            x = dq[head] - (t - 1) - wstart
            kend = wstart + (x % w) + k - 1
            if kend != last:
                out[no] = kend; no += 1; last = kend
    return np.unique(out[:no])


def sample_chunked(codes, scheme, k, w, s=None, greedy=None, chunk=6_000_000, forward=False):
    """Memory-bounded sampling: process overlapping chunks, keep positions owned by each chunk.
    Padding of 4*(w+k) bases on the left makes every kept window fully determined (verified
    bit-identical to the unchunked result in L1_validate.py)."""
    n = codes.shape[0]
    pad = 4 * (w + k) + 64
    P, ST, CN = [], [], []
    a = 0
    while a < n:
        b = min(n, a + chunk)
        lo = max(0, a - pad)
        hi = min(n, b + pad)
        sub = codes[lo:hi]
        if forward:
            pos = sample_forward(sub, scheme, k, w, s=s, greedy=greedy)
            st = np.zeros(len(pos), np.uint8)
            fwd, rev, valid = kmer_codes(sub, k)
            cn = fwd[pos]
        else:
            pos, st, cn = sample(sub, scheme, k, w, s=s, greedy=greedy)
        pos = pos + lo
        keep = (pos >= a) & (pos < b)
        P.append(pos[keep]); ST.append(st[keep]); CN.append(cn[keep])
        a = b
    return np.concatenate(P), np.concatenate(ST), np.concatenate(CN)


def revcomp_codes(codes):
    rc = codes[::-1].copy()
    m = rc < 4
    rc[m] = 3 - rc[m]
    return rc


# ----------------------------------------------------------------------------------------------
# Strong 64-bit mixer (splitmix64 finalizer) as an alternative random order.  Output >> 8 (56 bits)
# so that class prefixes still fit; ties at 56 bits are negligible and never cause false index hits
# (index keys are canonical 2-bit codes, not ranks).
# ----------------------------------------------------------------------------------------------
@nb.njit(inline='always')
def mix64(x):
    x = x ^ U64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> U64(30))) * U64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> U64(27))) * U64(0x94D049BB133111EB)
    x = x ^ (x >> U64(31))
    return x >> U64(8)


@nb.njit(cache=True)
def canon_mix(fwd, rev, valid):
    out = np.full(fwd.shape[0], U64MAX, np.uint64)
    for i in range(fwd.shape[0]):
        if valid[i]:
            c = fwd[i] if fwd[i] < rev[i] else rev[i]
            out[i] = mix64(c)
    return out


@nb.njit(cache=True)
def fwd_mix(fwd, valid):
    out = np.full(fwd.shape[0], U64MAX, np.uint64)
    for i in range(fwd.shape[0]):
        if valid[i]:
            out[i] = mix64(fwd[i])
    return out


@nb.njit(cache=True)
def smer_classes_mix(codes, k, s, canonical):
    """As smer_classes but with the strong mixer for the s-mer order."""
    n = codes.shape[0]
    fs, rs, vs = kmer_codes(codes, s)
    hs = np.full(n, U64MAX, np.uint64)
    for j in range(n):
        if vs[j]:
            c = fs[j]
            if canonical and rs[j] < c:
                c = rs[j]
            hs[j] = mix64(c ^ U64(0x5bd1e995))
    W = k - s + 1
    v = (k - s) // 2
    is_open = np.zeros(n, np.uint8); is_closed = np.zeros(n, np.uint8)
    dq = np.empty(n, np.int64); head = 0; tail = 0
    for j in range(n):
        while tail > head and hs[dq[tail - 1]] > hs[j]:
            tail -= 1
        dq[tail] = j; tail += 1
        first = j - W + 1
        while dq[head] < first:
            head += 1
        if first >= 0:
            m = hs[dq[head]]
            if m == U64MAX:
                continue
            if canonical:
                if hs[first] == m or hs[j] == m:
                    is_closed[j] = 1
                if hs[first + v] == m or hs[first + (k - s - v)] == m:
                    is_open[j] = 1
            else:
                p = dq[head] - first
                if p == 0 or p == k - s:
                    is_closed[j] = 1
                if p == v:
                    is_open[j] = 1
    return is_open, is_closed


# ----------------------------------------------------------------------------------------------
# Consolidated canonical driver used by the experiments.  All non-baseline schemes use the strong
# mixer (mix64) as their random order; 'mm2' is the bit-exact minimap2 v2.31 baseline.
#   s for oc/osync/csync and the inner s of modoc are chosen with (k - s) even (resp. t - s even) so
#   that the open-syncmer offset is the unique mirror-symmetric centre (canonical without loss).
# ----------------------------------------------------------------------------------------------
def parity_s(k, s_min=4):
    s = s_min
    while (k - s) % 2 != 0:
        s += 1
    return s


def sample2(codes, scheme, k, w, s=None, greedy=None):
    fwd, rev, valid = kmer_codes(codes, k)
    can = np.minimum(fwd, rev)
    strand = (fwd > rev).astype(np.uint8)
    if scheme == 'mm2':
        x, y = mm_sketch_generic(codes, w, k, np.empty(1, np.uint64), False)
        pos = (y >> U64(1)).astype(np.int64)
        return pos, (y & U64(1)).astype(np.uint8), can[pos]
    if scheme == 'mm2_1pick':      # hash64 order, but no emission of tied copies (classic winnowing)
        ch = canon_hash(fwd, rev, valid, k)
        pos = window_min_select_fwd(ch, valid, w, k, codes.shape[0])
        return pos, strand[pos], can[pos]
    if scheme == 'hashw':          # minimap2 window logic + width-adapted invertible 2k-bit hash (L1-B)
        x, y = mm_sketch_generic(codes, w, k, canon_hash_w(fwd, rev, valid, k), True)
        pos = (y >> U64(1)).astype(np.int64)
        return pos, (y & U64(1)).astype(np.uint8), can[pos]
    ch = canon_mix(fwd, rev, valid)
    if scheme in ('rand_mix', 'miniception', 'oc', 'ddecyc', 'greedy'):
        if scheme == 'rand_mix':
            rank = ch
        elif scheme == 'miniception':
            r = max(4, k - w)
            so, sc = smer_classes_mix(codes, k, r, True)
            rank = class_rank(ch, valid, np.where(sc == 1, 0, 1).astype(np.uint8))
        elif scheme == 'oc':
            ss = s if s is not None else parity_s(k, 5)
            so, sc = smer_classes_mix(codes, k, ss, True)
            rank = class_rank(ch, valid, np.where(so == 1, 0, np.where(sc == 1, 1, 2)).astype(np.uint8))
        elif scheme == 'ddecyc':
            rank = decycling_rank(fwd, rev, valid, ch, k, True)
        else:
            order, kp = greedy
            rank = greedy_rank(fwd, rev, valid, k, min(kp, k), order, True)
        x, y = mm_sketch_generic(codes, w, k, rank, True)
        pos = (y >> U64(1)).astype(np.int64)
        return pos, (y & U64(1)).astype(np.uint8), can[pos]
    if scheme in ('mod', 'modoc'):
        t = mod_t(k, w)
        tf, tr, tv = kmer_codes(codes, t)
        th = canon_mix(tf, tr, tv)
        if scheme == 'modoc' and t > 4:
            st = s if s is not None else parity_s(t, 4 if t >= 6 else 3)
            so, sc = smer_classes_mix(codes, t, st, True)
            th = class_rank(th, tv, np.where(so == 1, 0, np.where(sc == 1, 1, 2)).astype(np.uint8))
        pos = mod_sample(codes, w, k, t, th, ch, valid)
        return pos, strand[pos], can[pos]
    if scheme in ('osync', 'csync'):
        so, sc = smer_classes_mix(codes, k, s, True)
        sel = (so if scheme == 'osync' else sc).astype(bool) & valid.astype(bool)
        pos = np.nonzero(sel)[0].astype(np.int64)
        return pos, strand[pos], can[pos]
    raise ValueError(scheme)


def sample2_chunked(codes, scheme, k, w, s=None, greedy=None, chunk=6_000_000):
    n = codes.shape[0]
    pad = 4 * (w + k) + 64
    P, ST, CN = [], [], []
    a = 0
    while a < n:
        b = min(n, a + chunk)
        lo = max(0, a - pad); hi = min(n, b + pad)
        pos, st, cn = sample2(codes[lo:hi], scheme, k, w, s=s, greedy=greedy)
        pos = pos + lo
        keep = (pos >= a) & (pos < b)
        P.append(pos[keep]); ST.append(st[keep]); CN.append(cn[keep])
        a = b
    return np.concatenate(P), np.concatenate(ST), np.concatenate(CN)


# ----------------------------------------------------------------------------------------------
# Width-adapted invertible mixer: a bijection on b = 2k bits (xorshift-right by b/2 and odd
# multiplications modulo 2^b), proposed drop-in replacement for minimap2's hash64(key, mask), whose
# <<31 step and >>24/>>28 steps are (partly) no-ops when b <= 31 (k <= 15).
# ----------------------------------------------------------------------------------------------
@nb.njit(inline='always')
def hash_w(x, mask, b):
    h = U64(b // 2)
    x = (x ^ (x >> h)) & mask
    x = (x * U64(0xBF58476D1CE4E5B9)) & mask
    x = x ^ (x >> h)
    x = (x * U64(0x94D049BB133111EB)) & mask
    x = x ^ (x >> h)
    return x


@nb.njit(cache=True)
def canon_hash_w(fwd, rev, valid, k):
    mask = (U64(1) << U64(2 * k)) - U64(1)
    out = np.full(fwd.shape[0], U64MAX, np.uint64)
    for i in range(fwd.shape[0]):
        if valid[i]:
            c = fwd[i] if fwd[i] < rev[i] else rev[i]
            out[i] = hash_w(c, mask, 2 * k)
    return out


@nb.njit(cache=True)
def hash_w_is_bijection(k):
    mask = (U64(1) << U64(2 * k)) - U64(1)
    n = 1 << (2 * k)
    seen = np.zeros(n, np.uint8)
    for x in range(n):
        y = hash_w(U64(x), mask, 2 * k)
        if seen[y]:
            return False
        seen[y] = 1
    return True
