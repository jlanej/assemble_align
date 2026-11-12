"""L1_simulate.py -- read simulator with exact per-base truth (read base -> source coordinate).

Error models (per-base rates; HP = homopolymer-length change of a run of length L>=2 with
probability c*(L-1)^1.5, capped at 0.5, c calibrated on the source so the HP-error rate per base
equals the target):
  illumina : 2x150, insert N(400,60); substitutions rising quadratically 0.1% -> 0.5% along the
             read (mean 0.23%); indels 1e-5/base each.
  hifi     : ~Q30 total 1.0e-3/base split as in Wenger et al. 2019 (CCS discordances 92.0% homopolymer
             indels, 4.6% other indels, 3.4% mismatches): HP 9.2e-4 (65% contractions) + subs 3.4e-5
             + ins 2.3e-5 + del 2.3e-5; per-read rate multiplier lognormal(sigma 0.5, mean 1);
             length lognormal mean 16 kb, sd 4 kb (N50 ~18 kb, cf. NA19240 HPRC HiFi N50 17.8 kb).
  ont      : ~Q20 total 1.0e-2/base = HP 4.0e-3 (70% contractions) + subs 3.5e-3 + ins 1.0e-3
             + del 1.5e-3; multiplier lognormal(sigma 0.4); length lognormal mean 30 kb, sd 25 kb.
             (Illustrative split; homopolymer indels as the dominant residual R10.4 error is
             qualitative, Sereika et al. 2022.)
  ontpb    : pbsim3 3.0.5 errhmm ERRHMM-ONT-HQ.model, --accuracy-mean 0.99, length mean 30 kb sd 25 kb
             (trained FIC-HMM error model, Ono et al. 2022); truth parsed from pbsim3's MAF.  NOTE: the realised
             error rate is ~2.95% (sub 0.84%, ins 0.88%, del 1.24%; measured on chr20 reads), i.e. harder than
             Q20 -- used as a stress model, not as an R10.4.1 Q20 model.
  asm      : error-free tiles of a haplotype (real divergence only), fixed length.
Read orientation is random (50/50).  src[j] = source coordinate of read base j (-1 for inserted base).
"""
import math
import numpy as np
import numba as nb


@nb.njit(cache=True)
def _hp_rate_per_base(codes, c, gamma):
    n = codes.shape[0]
    tot = 0.0
    i = 0
    while i < n:
        j = i
        while j + 1 < n and codes[j + 1] == codes[i]:
            j += 1
        L = j - i + 1
        if L >= 2 and codes[i] < 4:
            p = c * (L - 1) ** gamma
            if p > 0.5:
                p = 0.5
            tot += p
        i = j + 1
    return tot / n


def calibrate_hp(codes, target, gamma=1.5):
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _hp_rate_per_base(codes, mid, gamma) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@nb.njit(cache=True)
def _apply_errors(frag, fsrc, p_sub_a, p_sub_b, p_ins, p_del, c_hp, gamma, p_contr, seed, out, osrc):
    """frag: oriented fragment codes; fsrc: source coordinate of each fragment base.
    Substitution probability for emitted base j of an expected-length-m read:
       p_sub_a + p_sub_b*(j/m)^2 (Illumina); for long reads p_sub_b=0.
    Returns number of bases written."""
    np.random.seed(seed)
    n = frag.shape[0]
    m = n
    no = 0
    i = 0
    while i < n:
        b = frag[i]
        j = i
        while j + 1 < n and frag[j + 1] == b:
            j += 1
        L = j - i + 1
        emit = L
        if L >= 2 and c_hp > 0.0:
            p = c_hp * (L - 1) ** gamma
            if p > 0.5:
                p = 0.5
            if np.random.random() < p:
                if np.random.random() < p_contr:
                    emit = L - 1
                else:
                    emit = L + 1
        for e in range(emit):
            q = i + (e if e < L else L - 1)          # source base backing this emitted base
            if np.random.random() < p_del:
                continue
            x = frag[q]
            ps = p_sub_a + p_sub_b * (no / m) ** 2
            if np.random.random() < ps:
                x = (x + 1 + np.random.randint(0, 3)) % 4
            out[no] = x
            osrc[no] = fsrc[q] if e < L else -1
            no += 1
            if np.random.random() < p_ins:
                out[no] = np.random.randint(0, 4)
                osrc[no] = -1
                no += 1
        i = j + 1
    return no


MODELS = {
    'illumina': dict(p_sub_a=0.001, p_sub_b=0.004, p_ins=1e-5, p_del=1e-5, hp=0.0, p_contr=0.5, sigma=0.0),
    'hifi': dict(p_sub_a=3.4e-5, p_sub_b=0.0, p_ins=2.3e-5, p_del=2.3e-5, hp=9.2e-4, p_contr=0.65, sigma=0.5),
    'ont': dict(p_sub_a=3.5e-3, p_sub_b=0.0, p_ins=1.0e-3, p_del=1.5e-3, hp=4.0e-3, p_contr=0.70, sigma=0.4),
}


def _lengths(kind, n, rng):
    if kind == 'hifi':
        mu, sd = 16000.0, 4000.0
        lo, hi = 2000, 40000
    elif kind == 'ont':
        mu, sd = 30000.0, 25000.0
        lo, hi = 2000, 300000
    else:
        raise ValueError(kind)
    s2 = math.log(1 + (sd / mu) ** 2)
    L = rng.lognormal(math.log(mu) - s2 / 2, math.sqrt(s2), n)
    return np.clip(L, lo, hi).astype(np.int64)


def revcomp(c):
    r = c[::-1].copy()
    m = r < 4
    r[m] = 3 - r[m]
    return r


def simulate(src_codes, kind, n, seed, c_hp=None, src_id=0):
    """Returns dict with concatenated read codes (reads separated by one N), read offsets/lengths,
    src coordinates (int32, -1 = inserted / separator), orientation, source interval, per-read error mult."""
    rng = np.random.default_rng(seed)
    N = src_codes.shape[0]
    reads, srcs, meta = [], [], []
    if kind == 'illumina':
        mdl = MODELS['illumina']
        npairs = n
        for p in range(npairs):
            for _try in range(100):
                ins = int(max(200, rng.normal(400, 60)))
                a = int(rng.integers(0, N - ins))
                if (src_codes[a:a + ins] < 4).all():
                    break
            frag_o = int(rng.integers(0, 2))
            # mate 1 forward from a (or reverse from end), mate 2 the opposite strand: FR library
            for mate in (0, 1):
                if (mate == 0) == (frag_o == 0):
                    fr = src_codes[a:a + 150].copy(); fs = np.arange(a, a + 150, dtype=np.int64); o = 0
                else:
                    fr = revcomp(src_codes[a + ins - 150:a + ins]); fs = np.arange(a + ins - 1, a + ins - 151, -1, dtype=np.int64); o = 1
                out = np.empty(400, np.uint8); osrc = np.empty(400, np.int64)
                m = _apply_errors(fr, fs, mdl['p_sub_a'], mdl['p_sub_b'], mdl['p_ins'], mdl['p_del'], 0.0, 1.5, 0.5,
                                  int(rng.integers(0, 2**31 - 1)), out, osrc)
                reads.append(out[:m]); srcs.append(osrc[:m])
                meta.append((p, mate, o, a, a + ins, 1.0, src_id))
    else:
        mdl = MODELS[kind]
        Ls = _lengths(kind, n, rng)
        mult = rng.lognormal(-mdl['sigma'] ** 2 / 2, mdl['sigma'], n) if mdl['sigma'] > 0 else np.ones(n)
        for r in range(n):
            L = int(min(Ls[r], N - 1))
            for _try in range(100):
                a = int(rng.integers(0, N - L))
                if (src_codes[a:a + L] < 4).all():
                    break
            o = int(rng.integers(0, 2))
            if o == 0:
                fr = src_codes[a:a + L].copy(); fs = np.arange(a, a + L, dtype=np.int64)
            else:
                fr = revcomp(src_codes[a:a + L]); fs = np.arange(a + L - 1, a - 1, -1, dtype=np.int64)
            f = float(mult[r])
            out = np.empty(int(L * 1.2) + 100, np.uint8); osrc = np.empty(int(L * 1.2) + 100, np.int64)
            m = _apply_errors(fr, fs, mdl['p_sub_a'] * f, 0.0, mdl['p_ins'] * f, mdl['p_del'] * f, c_hp * f, 1.5,
                              mdl['p_contr'], int(rng.integers(0, 2**31 - 1)), out, osrc)
            reads.append(out[:m]); srcs.append(osrc[:m])
            meta.append((r, 0, o, a, a + L, f, src_id))
    return pack(reads, srcs, meta)


def tiles(src_codes, tile, src_id=0, step=None):
    """Error-free haplotype tiles (assembly-to-reference scenario)."""
    step = step or tile
    N = src_codes.shape[0]
    reads, srcs, meta = [], [], []
    for r, a in enumerate(range(0, N - tile + 1, step)):
        reads.append(src_codes[a:a + tile].copy()); srcs.append(np.arange(a, a + tile, dtype=np.int64))
        meta.append((r, 0, 0, a, a + tile, 1.0, src_id))
    return pack(reads, srcs, meta)


def pack(reads, srcs, meta):
    lens = np.array([len(x) for x in reads], np.int64)
    offs = np.zeros(len(reads), np.int64)
    offs[1:] = np.cumsum(lens + 1)[:-1]
    tot = int(lens.sum() + len(reads))
    codes = np.full(tot, 4, np.uint8)
    src = np.full(tot, -1, np.int32)
    for i, (x, s) in enumerate(zip(reads, srcs)):
        codes[offs[i]:offs[i] + lens[i]] = x
        src[offs[i]:offs[i] + lens[i]] = s
    meta = np.array(meta, dtype=[('rid', 'i8'), ('mate', 'i8'), ('orient', 'i8'), ('st', 'i8'), ('en', 'i8'),
                                 ('mult', 'f8'), ('src_id', 'i8')])
    return dict(codes=codes, src=src, offs=offs, lens=lens, meta=meta)


def error_summary(rs, src_codes):
    """Empirical error composition: compare each read to its source via src mapping."""
    codes, src = rs['codes'], rs['src']
    m = codes < 4
    ins = int(((src == -1) & m).sum())
    mapped = (src >= 0) & m
    sub = int((codes[mapped] != np.where(rs_orient_mask(rs)[mapped], 3 - src_codes[src[mapped]], src_codes[src[mapped]])).sum())
    tot_bases = int(m.sum())
    # deletions: source bases not represented = span length - mapped count
    span = int((rs['meta']['en'] - rs['meta']['st']).sum())
    return dict(bases=tot_bases, sub=sub / tot_bases, ins=ins / tot_bases)


def rs_orient_mask(rs):
    om = np.zeros(rs['codes'].shape[0], bool)
    for i in range(len(rs['offs'])):
        if rs['meta']['orient'][i] == 1:
            om[rs['offs'][i]:rs['offs'][i] + rs['lens'][i]] = True
    return om


def simulate_pbsim3(fasta, src_codes, depth, seed, src_id=0, model=None, acc=0.99, mean=30000, sd=25000):
    """Run pbsim3 (errhmm) on `fasta` in a self-cleaning temp dir and import reads with per-base truth."""
    import gzip, os, subprocess, tempfile, shutil
    model = model or os.path.join(os.environ.get('CONDA_PREFIX', ''), 'data', 'ERRHMM-ONT-HQ.model')
    d = tempfile.mkdtemp(prefix='pbsim3_', dir=os.environ.get('L1_SCRATCH', '.'))
    try:
        subprocess.run(['pbsim', '--strategy', 'wgs', '--method', 'errhmm', '--errhmm', model, '--depth', str(depth),
                        '--genome', os.path.abspath(fasta), '--length-mean', str(mean), '--length-sd', str(sd),
                        '--accuracy-mean', str(acc), '--seed', str(seed), '--prefix', 'p'], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        enc = np.full(256, 4, np.uint8)
        for c, v in zip(b'ACGTacgt', [0, 1, 2, 3, 0, 1, 2, 3]):
            enc[c] = v
        reads, srcs, meta = [], [], []
        with gzip.open(os.path.join(d, 'p_0001.maf.gz'), 'rt') as f:
            rows = []
            for line in f:
                if line.startswith('s '):
                    rows.append(line.split())
                    if len(rows) == 2:
                        ref, rd = rows; rows = []
                        rstart = int(ref[2]); rstr = rd[4]
                        ra = np.frombuffer(ref[6].encode(), np.uint8); qa = np.frombuffer(rd[6].encode(), np.uint8)
                        rgap = ra == ord('-'); qgap = qa == ord('-')
                        rpos = rstart + np.cumsum(~rgap) - 1
                        col_src = np.where(rgap, -1, rpos)[~qgap]
                        seq = enc[qa[~qgap]]
                        if rstr == '-':      # MAF row is in reference orientation: read = revcomp(row)
                            seq = revcomp(seq); col_src = col_src[::-1].copy(); o = 1
                        else:
                            o = 0
                        reads.append(seq); srcs.append(col_src.astype(np.int64))
                        meta.append((len(meta), 0, o, rstart, rstart + int(ref[3]), 1.0, src_id))
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return pack(reads, srcs, meta)
